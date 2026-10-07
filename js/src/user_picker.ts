/**
 * Searchable participant picker. Shows every participant (virtualized, so tens
 * of thousands of rows stay fast) with vote/statement counts, All/Voters/
 * Commenters filters, sorting, keyboard navigation and prev/next stepping.
 */

export type UserFilter = "all" | "voters" | "commenters";
type UserSort = "id" | "votes" | "statements";

interface Row {
  id: string;
  votes: number;
  statements: number;
}

const ROW_HEIGHT = 26;
const LIST_HEIGHT = 260;
const OVERSCAN = 6;

const fmtInt = new Intl.NumberFormat();
const plural = (n: number, word: string) => `${fmtInt.format(n)} ${word}${n === 1 ? "" : "s"}`;

function el<K extends keyof HTMLElementTagNameMap>(tag: K, className?: string, text?: string) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
}

export class UserPicker {
  readonly root = el("div", "vv-picker");

  private rows: Row[] = [];
  private byId = new Map<string, Row>();
  private shown: Row[] = []; // after filter + search + sort
  private current = "";
  private active = -1; // keyboard-highlighted index into `shown`
  private filter: UserFilter = "all";
  private sort: UserSort = "id";
  private query = "";

  private input = el("input", "vv-input vv-picker-input");
  private panel = el("div", "vv-picker-panel");
  private tabs = el("div", "vv-picker-tabs");
  private status = el("div", "vv-picker-status");
  private list = el("div", "vv-picker-list");
  private spacer = el("div", "vv-picker-spacer");
  private sortSelect = el("select", "vv-picker-sort");

  constructor(private onSelect: (id: string) => void) {
    const label = el("label", "vv-label", "User ID:");
    const field = el("div", "vv-picker-field");
    const toggle = el("button", "vv-btn vv-picker-toggle", "▾");
    toggle.title = "Browse all participants";
    toggle.type = "button";
    this.input.placeholder = "search id…";
    this.input.autocomplete = "off";
    this.input.spellcheck = false;
    field.append(this.input, toggle);
    label.append(field);

    const prev = el("button", "vv-btn", "◀");
    const next = el("button", "vv-btn", "▶");
    prev.title = "Previous participant in the list";
    next.title = "Next participant in the list";
    prev.addEventListener("click", () => this.step(-1));
    next.addEventListener("click", () => this.step(1));

    for (const [value, text] of [
      ["id", "Sort: ID"],
      ["votes", "Sort: most votes"],
      ["statements", "Sort: most statements"],
    ] as const) {
      const opt = el("option", undefined, text);
      opt.value = value;
      this.sortSelect.append(opt);
    }
    this.sortSelect.addEventListener("change", () => {
      this.sort = this.sortSelect.value as UserSort;
      this.refresh(true);
    });

    const controls = el("div", "vv-picker-controls");
    controls.append(this.tabs, this.sortSelect);
    this.list.style.height = `${LIST_HEIGHT}px`;
    this.list.append(this.spacer);
    this.panel.append(controls, this.status, this.list);
    this.panel.hidden = true;

    this.root.append(label, prev, next, this.panel);

    this.input.addEventListener("focus", () => {
      this.input.select();
      this.open();
    });
    this.input.addEventListener("input", () => {
      this.query = this.input.value.trim();
      this.open();
      this.refresh(true);
    });
    this.input.addEventListener("keydown", (e) => this.onKeydown(e));
    toggle.addEventListener("click", () => (this.panel.hidden ? this.open(true) : this.close()));
    this.list.addEventListener("scroll", () => this.paint());
    this.list.addEventListener("pointerdown", (e) => e.preventDefault()); // keep input focus
    this.list.addEventListener("click", (e) => {
      const id = (e.target as HTMLElement).closest<HTMLElement>(".vv-picker-row")?.dataset.id;
      if (id) this.choose(id);
    });
    document.addEventListener("pointerdown", this.onOutside);
  }

  destroy() {
    document.removeEventListener("pointerdown", this.onOutside);
  }

  setUsers(ids: string[], votes: number[], statements: number[]) {
    this.rows = ids.map((id, i) => ({ id, votes: votes[i] ?? 0, statements: statements[i] ?? 0 }));
    this.byId = new Map(this.rows.map((r) => [r.id, r]));
    this.renderTabs();
    this.refresh(false);
  }

  setCurrent(id: string) {
    this.current = id;
    if (document.activeElement !== this.input) this.input.value = id;
    this.paint();
  }

  /** Ids for the random buttons. */
  ids(kind: UserFilter): string[] {
    return this.rows.filter((r) => this.matchesFilter(r, kind)).map((r) => r.id);
  }

  private matchesFilter(r: Row, kind: UserFilter) {
    return kind === "all" || (kind === "voters" ? r.votes > 0 : r.statements > 0);
  }

  private renderTabs() {
    const count = (kind: UserFilter) => this.rows.filter((r) => this.matchesFilter(r, kind)).length;
    this.tabs.replaceChildren(
      ...(["all", "voters", "commenters"] as const).map((kind) => {
        const name = { all: "All", voters: "Voters", commenters: "Commenters" }[kind];
        const tab = el("button", "vv-picker-tab", `${name} (${fmtInt.format(count(kind))})`);
        tab.type = "button";
        tab.classList.toggle("is-active", kind === this.filter);
        tab.addEventListener("pointerdown", (e) => e.preventDefault()); // keep input focus
        tab.addEventListener("click", () => {
          this.filter = kind;
          this.renderTabs();
          this.refresh(true);
        });
        return tab;
      }),
    );
  }

  private refresh(resetScroll: boolean) {
    const q = this.query;
    let shown = this.rows.filter((r) => this.matchesFilter(r, this.filter) && (!q || r.id.includes(q)));
    if (q) {
      // Exact match first, then ids starting with the query.
      shown = shown
        .map((r, i) => ({ r, i, rank: r.id === q ? 0 : r.id.startsWith(q) ? 1 : 2 }))
        .sort((a, b) => a.rank - b.rank || a.i - b.i)
        .map(({ r }) => r);
    }
    if (this.sort !== "id") {
      const key = this.sort;
      shown = [...shown].sort((a, b) => b[key] - a[key]); // stable: ties keep id order
    }
    this.shown = shown;
    this.active = q && shown.length ? 0 : shown.findIndex((r) => r.id === this.current);
    this.status.textContent =
      shown.length === this.rows.length
        ? `${fmtInt.format(shown.length)} participants`
        : `${fmtInt.format(shown.length)} of ${fmtInt.format(this.rows.length)} participants`;
    this.spacer.style.height = `${shown.length * ROW_HEIGHT}px`;
    if (resetScroll) this.list.scrollTop = 0;
    this.scrollToActive();
    this.paint();
  }

  /** Render only the rows in (and just around) the visible window. */
  private paint() {
    if (this.panel.hidden) return;
    const first = Math.max(0, Math.floor(this.list.scrollTop / ROW_HEIGHT) - OVERSCAN);
    const last = Math.min(this.shown.length, first + Math.ceil(LIST_HEIGHT / ROW_HEIGHT) + OVERSCAN * 2);
    const nodes: HTMLElement[] = [];
    for (let i = first; i < last; i++) {
      const r = this.shown[i];
      const row = el("div", "vv-picker-row");
      row.dataset.id = r.id;
      row.style.top = `${i * ROW_HEIGHT}px`;
      row.classList.toggle("is-current", r.id === this.current);
      row.classList.toggle("is-active", i === this.active);
      row.append(
        el("span", "vv-picker-id", r.id),
        el("span", "vv-picker-meta", `${plural(r.votes, "vote")} · ${plural(r.statements, "statement")}`),
      );
      nodes.push(row);
    }
    this.spacer.replaceChildren(...nodes);
  }

  private scrollToActive() {
    if (this.active < 0) return;
    const top = this.active * ROW_HEIGHT;
    if (top < this.list.scrollTop) this.list.scrollTop = top;
    else if (top + ROW_HEIGHT > this.list.scrollTop + LIST_HEIGHT) {
      this.list.scrollTop = top + ROW_HEIGHT - LIST_HEIGHT;
    }
  }

  private open(focusInput = false) {
    if (focusInput) this.input.focus();
    if (!this.panel.hidden) return;
    this.panel.hidden = false;
    this.refresh(false);
    // Centre the current participant when opening.
    if (this.active >= 0) this.list.scrollTop = this.active * ROW_HEIGHT - LIST_HEIGHT / 2;
    this.paint();
  }

  private close() {
    this.panel.hidden = true;
    this.query = "";
    this.input.value = this.current;
  }

  private choose(id: string) {
    this.close();
    this.input.blur();
    if (id !== this.current) this.onSelect(id);
  }

  /** Move to the previous/next participant in the current filter + sort order. */
  private step(delta: number) {
    if (!this.shown.length) return;
    const i = this.shown.findIndex((r) => r.id === this.current);
    const next = i < 0 ? 0 : (i + delta + this.shown.length) % this.shown.length;
    this.onSelect(this.shown[next].id);
  }

  private onKeydown(e: KeyboardEvent) {
    if (e.key === "ArrowDown" || e.key === "ArrowUp") {
      e.preventDefault();
      this.open();
      if (!this.shown.length) return;
      const delta = e.key === "ArrowDown" ? 1 : -1;
      this.active = Math.min(this.shown.length - 1, Math.max(0, this.active + delta));
      this.scrollToActive();
      this.paint();
    } else if (e.key === "Enter") {
      e.preventDefault();
      const typed = this.input.value.trim();
      const target = this.shown[this.active]?.id ?? (this.byId.has(typed) ? typed : undefined);
      if (target) this.choose(target);
    } else if (e.key === "Escape") {
      this.close();
      this.input.blur();
    }
  }

  private onOutside = (e: PointerEvent) => {
    if (!this.panel.hidden && !this.root.contains(e.target as Node)) this.close();
  };
}
