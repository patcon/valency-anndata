# ------------------------------------------------------------
# SVG primitives
# ------------------------------------------------------------
def draw_grid_block(
    dwg,
    *,
    x,
    y,
    width,
    height,
    rows,
    cols,
    label,
    stroke="#333",
    grid_stroke="#ccc",
):
    group = dwg.g()

    # Outer rectangle
    group.add(
        dwg.rect(
            insert=(x, y),
            size=(width, height),
            fill="white",
            # fill_opacity=0.9,
            stroke=stroke,
            stroke_width=2,
        )
    )

    # Horizontal grid
    if rows > 1:
        row_h = height / rows
        for i in range(1, rows):
            group.add(
                dwg.line(
                    start=(x, y + i * row_h),
                    end=(x + width, y + i * row_h),
                    stroke=grid_stroke,
                )
            )

    # Vertical grid
    if cols > 1:
        col_w = width / cols
        for j in range(1, cols):
            group.add(
                dwg.line(
                    start=(x + j * col_w, y),
                    end=(x + j * col_w, y + height),
                    stroke=grid_stroke,
                )
            )

    # Label (centered)
    lines = label.split("\n")
    line_height = 14
    start_y = y + height / 2 - (len(lines) - 1) * line_height / 2

    for i, line in enumerate(lines):
        group.add(
            dwg.text(
                line,
                insert=(x + width / 2, start_y + i * line_height),
                text_anchor="middle",
                font_size=12,
                font_family="sans-serif",
            )
        )

    dwg.add(group)