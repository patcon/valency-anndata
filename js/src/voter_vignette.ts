import type { RenderProps } from "@anywidget/types";
import { axisBottom, axisLeft } from "d3-axis";
import { scaleLinear, scaleUtc, type ScaleTime } from "d3-scale";
import { select } from "d3-selection";
import { utcDay, utcMinute, utcMonth, utcYear } from "d3-time";
import { utcFormat } from "d3-time-format";
import { zoom as d3zoom, zoomIdentity, zoomTransform, type ZoomTransform } from "d3-zoom";
import "./voter_vignette.css";

interface Vote {
  t: number; // epoch ms
  vote: -1 | 0 | 1;
  statement_id: string;
  content: string | null;
  moderation_state: number | null;
}

interface Statement {
  t: number; // epoch ms
  statement_id: string;
  content: string | null;
  moderation_state: number | null;
}

interface UserData {
  user_id: string;
  votes: Vote[];
  statements: Statement[];
}

interface Model {
  all_users: string[];
  voters: string[];
  commenters: string[];
  user_id: string;
  user_data: UserData;
}

const VOTE_LABEL: Record<number, string> = { 1: "Agree", 0: "Pass", [-1]: "Disagree" };
const VOTE_CLASS: Record<number, string> = { 1: "agree", 0: "pass", [-1]: "disagree" };
const MOD_CLASS: Record<number, string> = { 1: "accepted", 0: "unmoderated", [-1]: "rejected" };
const MOD_LABEL: Record<number, string> = { 1: "accepted", 0: "unmoderated", [-1]: "moderated out" };

const HEIGHT = 240;
const MARGIN = { top: 12, right: 16, bottom: 28, left: 72 };
const MIN_WINDOW_MS = 60 * 1000; // max zoom: ~1 minute across the plot
const FOCUS_HALF_WINDOW_MS = 30 * 60 * 1000; // ±30 min around a clicked statement
const HOVER_RADIUS_PX = 8;
const MAX_TOOLTIP_VOTES = 8;

const fmtFull = utcFormat("%Y-%m-%d %H:%M:%S UTC");

// Multi-scale 24h UTC tick labels (d3's default mixes "01:40" with "02 AM").
const fmtSecond = utcFormat(":%S");
const fmtMinute = utcFormat("%H:%M");
const fmtDay = utcFormat("%b %d");
const fmtMonth = utcFormat("%B");
const fmtYear = utcFormat("%Y");
function fmtTick(date: Date): string {
  if (utcMinute(date) < date) return fmtSecond(date);
  if (utcDay(date) < date) return fmtMinute(date);
  if (utcMonth(date) < date) return fmtDay(date);
  if (utcYear(date) < date) return fmtMonth(date);
  return fmtYear(date);
}

function el<K extends keyof HTMLElementTagNameMap>(
  tag: K,
  className?: string,
  text?: string,
): HTMLElementTagNameMap[K] {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
}

/** Same adaptive rules as v1: minutes under an hour, hours under a day, else days. */
function formatDuration(ms: number): string {
  const minutes = ms / 60000;
  if (minutes < 60) return `${minutes.toFixed(1)} minutes`;
  const hours = minutes / 60;
  if (hours < 24) return `${hours.toFixed(1)} hours`;
  return `${Math.floor(hours / 24)} days`;
}

function render({ model, el: root }: RenderProps<Model>) {
  root.classList.add("vv-root");

  // ── Header: user picker + random buttons ──────────────────────────────
  const header = el("div", "vv-header");
  const label = el("label", "vv-label", "User ID:");
  const input = el("input", "vv-input");
  const listId = `vv-users-${Math.random().toString(36).slice(2)}`;
  input.setAttribute("list", listId);
  input.placeholder = "type or pick a user id";
  const datalist = el("datalist");
  datalist.id = listId;
  const randomVoterBtn = el("button", "vv-btn", "Random voter");
  const randomCommenterBtn = el("button", "vv-btn", "Random commenter");
  const resetBtn = el("button", "vv-btn", "Reset zoom");
  label.append(input);
  header.append(label, datalist, randomVoterBtn, randomCommenterBtn, resetBtn);

  const title = el("div", "vv-title");

  // ── Chart ──────────────────────────────────────────────────────────────
  const chart = el("div", "vv-chart");
  const tooltip = el("div", "vv-tooltip");
  tooltip.hidden = true;
  chart.append(tooltip);

  const legend = el("div", "vv-legend");
  const list = el("div", "vv-statements");
  const hint = el(
    "div",
    "vv-hint",
    "Scroll to zoom · drag to pan · double-click to reset · hover a vote to see its statement",
  );

  root.append(header, title, chart, legend, hint, list);

  const clipId = `vv-clip-${Math.random().toString(36).slice(2)}`;
  const svg = select(chart).append("svg").attr("class", "vv-svg").attr("height", HEIGHT);
  svg.append("defs").append("clipPath").attr("id", clipId).append("rect");
  const xAxisG = svg.append("g").attr("class", "vv-axis vv-axis-x");
  const yAxisG = svg.append("g").attr("class", "vv-axis vv-axis-y");
  const plot = svg.append("g").attr("clip-path", `url(#${clipId})`);
  const statementsG = plot.append("g");
  const votesG = plot.append("g");

  const y = scaleLinear().domain([-1.5, 1.5]);
  let x0: ScaleTime<number, number> = scaleUtc();
  let x: ScaleTime<number, number> = x0;
  let width = 0;
  let data: UserData = { user_id: "", votes: [], statements: [] };
  let hovered = new Set<Vote>();

  const zoom = d3zoom<SVGSVGElement, unknown>().on("zoom", (event) => {
    x = (event.transform as ZoomTransform).rescaleX(x0);
    draw();
  });
  svg.call(zoom).on("dblclick.zoom", null).on("dblclick", () => resetZoom());

  function plotWidth() {
    return Math.max(0, width - MARGIN.left - MARGIN.right);
  }

  function resetZoom() {
    svg.call(zoom.transform, zoomIdentity);
  }

  function updateExtent() {
    const w = plotWidth();
    const [d0, d1] = x0.domain().map((d) => d.getTime());
    const maxK = Math.max(1, (d1 - d0) / MIN_WINDOW_MS);
    zoom
      .scaleExtent([1, maxK])
      .extent([
        [MARGIN.left, 0],
        [MARGIN.left + w, HEIGHT],
      ])
      .translateExtent([
        [MARGIN.left, 0],
        [MARGIN.left + w, HEIGHT],
      ]);
  }

  /** Recompute the base (unzoomed) domain from the current user's data. */
  function setDomain() {
    const times = [...data.votes.map((v) => v.t), ...data.statements.map((s) => s.t)];
    let lo: number;
    let hi: number;
    if (times.length === 0) {
      lo = Date.now() - FOCUS_HALF_WINDOW_MS;
      hi = Date.now() + FOCUS_HALF_WINDOW_MS;
    } else {
      lo = Math.min(...times);
      hi = Math.max(...times);
      const pad = Math.max((hi - lo) * 0.03, 60 * 1000);
      lo -= pad;
      hi += pad;
    }
    x0 = scaleUtc()
      .domain([new Date(lo), new Date(hi)])
      .range([MARGIN.left, MARGIN.left + plotWidth()]);
  }

  function draw() {
    const w = plotWidth();
    y.range([HEIGHT - MARGIN.bottom, MARGIN.top]);
    svg.attr("width", width);
    svg
      .select(`#${clipId} rect`)
      .attr("x", MARGIN.left)
      .attr("y", 0)
      .attr("width", w)
      .attr("height", HEIGHT);

    xAxisG
      .attr("transform", `translate(0,${HEIGHT - MARGIN.bottom})`)
      .call(
        axisBottom(x)
          .ticks(Math.max(2, Math.floor(w / 110)))
          .tickFormat((d) => fmtTick(d as Date)),
      );
    yAxisG
      .attr("transform", `translate(${MARGIN.left},0)`)
      .call(
        axisLeft(y)
          .tickValues([-1, 0, 1])
          .tickFormat((d) => VOTE_LABEL[d as number]),
      );

    statementsG
      .selectAll<SVGLineElement, Statement>("line")
      .data(data.statements)
      .join("line")
      .attr("class", (s) => `vv-statement-line mod-${MOD_CLASS[s.moderation_state ?? 0]}`)
      .attr("x1", (s) => x(s.t))
      .attr("x2", (s) => x(s.t))
      .attr("y1", MARGIN.top)
      .attr("y2", HEIGHT - MARGIN.bottom);

    votesG
      .selectAll<SVGCircleElement, Vote>("circle")
      .data(data.votes)
      .join("circle")
      .attr("class", (v) => `vv-vote vote-${VOTE_CLASS[v.vote]}`)
      .classed("is-hovered", (v) => hovered.has(v))
      .attr("cx", (v) => x(v.t))
      .attr("cy", (v) => y(v.vote))
      .attr("r", 5);
  }

  // ── Hover: nearest votes within a small radius, else a statement line ──
  function showTooltip(html: HTMLElement[], px: number, py: number) {
    tooltip.replaceChildren(...html);
    tooltip.hidden = false;
    const cw = chart.clientWidth;
    const tw = tooltip.offsetWidth;
    const left = px + 14 + tw > cw ? Math.max(0, px - 14 - tw) : px + 14;
    tooltip.style.left = `${left}px`;
    tooltip.style.top = `${Math.max(0, py - 10)}px`;
  }

  function hideTooltip() {
    tooltip.hidden = true;
    if (hovered.size) {
      hovered = new Set();
      draw();
    }
  }

  function voteRow(v: Vote): HTMLElement {
    const row = el("div", "vv-tip-row");
    const head = el("div", "vv-tip-head");
    head.append(
      el("span", `vv-chip vote-${VOTE_CLASS[v.vote]}`, VOTE_LABEL[v.vote]),
      el("span", "vv-tip-meta", ` #${v.statement_id} · ${fmtFull(new Date(v.t))}`),
    );
    row.append(head, el("div", "vv-tip-text", v.content ?? "(statement text unavailable)"));
    if (v.moderation_state === -1) {
      row.append(el("div", "vv-tip-note", "This statement was moderated out."));
    }
    return row;
  }

  svg.on("pointermove", (event: PointerEvent) => {
    if (event.buttons) return hideTooltip(); // dragging
    const [mx, my] = [event.offsetX, event.offsetY];
    if (mx < MARGIN.left || mx > MARGIN.left + plotWidth()) return hideTooltip();

    const near = data.votes
      .map((v) => ({ v, d: Math.hypot(x(v.t) - mx, y(v.vote) - my) }))
      .filter(({ d }) => d <= HOVER_RADIUS_PX)
      .sort((a, b) => a.v.t - b.v.t);

    if (near.length) {
      hovered = new Set(near.map(({ v }) => v));
      draw();
      const rows = near.slice(0, MAX_TOOLTIP_VOTES).map(({ v }) => voteRow(v));
      if (near.length > MAX_TOOLTIP_VOTES) {
        rows.push(el("div", "vv-tip-note", `+${near.length - MAX_TOOLTIP_VOTES} more — zoom in to separate`));
      }
      return showTooltip(rows, mx, my);
    }

    const stmt = data.statements.find((s) => Math.abs(x(s.t) - mx) <= 4);
    if (stmt) {
      if (hovered.size) {
        hovered = new Set();
        draw();
      }
      const row = el("div", "vv-tip-row");
      const head = el("div", "vv-tip-head");
      head.append(
        el("span", `vv-chip mod-${MOD_CLASS[stmt.moderation_state ?? 0]}`, "Authored"),
        el("span", "vv-tip-meta", ` #${stmt.statement_id} · ${fmtFull(new Date(stmt.t))}`),
      );
      row.append(head, el("div", "vv-tip-text", stmt.content ?? ""));
      return showTooltip([row], mx, my);
    }

    hideTooltip();
  });
  svg.on("pointerleave", hideTooltip);

  // ── Text sections ──────────────────────────────────────────────────────
  function focusOn(t: number) {
    const [d0, d1] = x0.domain().map((d) => d.getTime());
    const lo = Math.max(d0, t - FOCUS_HALF_WINDOW_MS);
    const hi = Math.min(d1, t + FOCUS_HALF_WINDOW_MS);
    const k = (x0(d1) - x0(d0)) / Math.max(1, x0(hi) - x0(lo));
    const transform = zoomIdentity
      .translate(MARGIN.left, 0)
      .scale(k)
      .translate(-x0(lo), 0);
    svg.call(zoom.transform, transform);
  }

  function renderText() {
    const { votes, statements, user_id } = data;
    if (votes.length) {
      const first = votes[0].t;
      const last = votes[votes.length - 1].t;
      title.textContent =
        `User ${user_id} activity | ${fmtFull(new Date(first))} → ${fmtFull(new Date(last))}` +
        ` (${formatDuration(last - first)})`;
    } else {
      title.textContent = `User ${user_id} activity | No votes`;
    }

    legend.replaceChildren(
      el("span", "vv-legend-item vv-legend-vote", `Votes (${votes.length})`),
      el("span", "vv-legend-item vv-legend-statement", `Statements (${statements.length})`),
    );

    list.replaceChildren();
    if (statements.length) {
      list.append(el("div", "vv-statements-title", `Statements by ${user_id} in submission order:`));
      const ul = el("ul");
      for (const s of statements) {
        const li = el("li", `vv-statement mod-${MOD_CLASS[s.moderation_state ?? 0]}`);
        li.title = `Zoom to ±30 min around this statement (${MOD_LABEL[s.moderation_state ?? 0]})`;
        li.append(el("span", "vv-statement-time", fmtFull(new Date(s.t))), el("span", "", ` ${s.content ?? ""}`));
        li.addEventListener("click", () => focusOn(s.t));
        ul.append(li);
      }
      list.append(ul);
    }
  }

  // ── Model sync ─────────────────────────────────────────────────────────
  function setUser(id: string) {
    model.set("user_id", id);
    model.save_changes();
  }

  function pick(ids: string[]) {
    if (ids.length) setUser(ids[Math.floor(Math.random() * ids.length)]);
  }

  function syncUsers() {
    datalist.replaceChildren(
      ...model.get("all_users").map((id) => {
        const opt = document.createElement("option");
        opt.value = id;
        return opt;
      }),
    );
  }

  function syncData() {
    data = model.get("user_data") ?? { user_id: "", votes: [], statements: [] };
    input.value = model.get("user_id");
    hovered = new Set();
    tooltip.hidden = true;
    setDomain();
    updateExtent();
    renderText();
    x = x0;
    resetZoom(); // fires the zoom handler, which draws
  }

  input.addEventListener("change", () => {
    const id = input.value.trim();
    if (model.get("all_users").includes(id)) setUser(id);
    else input.value = model.get("user_id");
  });
  randomVoterBtn.addEventListener("click", () => pick(model.get("voters")));
  randomCommenterBtn.addEventListener("click", () => pick(model.get("commenters")));
  resetBtn.addEventListener("click", resetZoom);

  const onUsers = () => syncUsers();
  const onData = () => syncData();
  model.on("change:all_users", onUsers);
  model.on("change:user_data", onData);

  // Re-layout on resize, keeping the current zoom transform.
  const resize = new ResizeObserver(() => {
    const w = chart.clientWidth;
    if (w === width || w === 0) return;
    width = w;
    x0.range([MARGIN.left, MARGIN.left + plotWidth()]);
    updateExtent();
    const node = svg.node();
    x = (node ? zoomTransform(node) : zoomIdentity).rescaleX(x0);
    draw();
  });
  resize.observe(chart);

  width = chart.clientWidth || 800;
  syncUsers();
  syncData();

  return () => {
    resize.disconnect();
    model.off("change:all_users", onUsers);
    model.off("change:user_data", onData);
  };
}

export default { render };
