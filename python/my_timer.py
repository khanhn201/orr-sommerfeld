import time
from collections import defaultdict

_time   = defaultdict(float)   # accumulated seconds per name
_count  = defaultdict(int)     # number of completed intervals per name
_start  = {}                   # active starts per name (single slot = simplest)

# ---- api ----
def tic(name="default"):
    _start[name] = time.perf_counter()
    return _start[name]

def toc(name="default", t0=None, message=None, out=False):
    if t0 is None:
        t0 = _start.pop(name, None)
        if t0 is None:
            print(f"tic('{name}') was not called!")
            return None
    dt = time.perf_counter() - t0
    _time[name]  += dt
    _count[name] += 1
    if out:
        print(f"{message or name}: {dt:.6f} s")
    return dt

def reset_timers():
    _time.clear(); _count.clear(); _start.clear()

# ---- printing helpers (sum-at-print-time only) ----
def _children_of(prefix):
    """Immediate children under 'prefix' ('' = top level)."""
    plen = len(prefix)
    kids = set()
    for k in _time.keys():
        if prefix:
            if not (k == prefix or k.startswith(prefix + ":")):
                continue
            tail = k[plen+1:] if k != prefix else ""
        else:
            tail = k
        if ":" in tail:
            kids.add((tail.split(":")[0]))
        elif tail != "":
            kids.add(tail)
    return sorted(kids)

def _subkey(prefix, child):
    return f"{prefix}:{child}" if prefix else child

def _subtree_sum(prefix):
    tot_t = 0.0; tot_c = 0
    if prefix:
        pfx = prefix + ":"
        for k in _time.keys():
            if k == prefix or k.startswith(pfx):
                tot_t += _time[k]; tot_c += _count[k]
    return tot_t, tot_c


# --- helpers you already have: _children_of, _subkey, _subtree_sum ---

def _print_node(prefix, parent_total, level, lines,
                name_w, time_w, pct_w, count_w,
                name_indent=2, col_indent=2, base_col_levels=2,
                count_start_col=0):
    """
    Layout per row:
      [name_indent*level spaces][label padded to (name_w - name_indent*level)]
      ' ' +
      [col_indent*(base_col_levels + level) spaces][time right][ ' s  ' ][% right]
      '  ' +
      [Count left-aligned at fixed start column]
    """
    tot_t, tot_c = _subtree_sum(prefix)
    label = prefix.split(":")[-1]

    # visible name indent
    name_pad = " " * (name_indent * level)
    # column pad shares the same baseline as header: base_col_levels
    col_pad  = " " * (col_indent * (base_col_levels + level))

    # keep the time/% columns aligned by fixing the label field width
    avail = max(1, name_w - len(name_pad))
    label_field = f"{label:<{avail}}"

    # % of parent
    pct_str = ""
    if parent_total and parent_total > 0:
        pct = 100.0 * tot_t / parent_total
        pct_str = f"{pct:>{pct_w-1}.1f}%"

    left = (
        f"{name_pad}{label_field} "
        f"{col_pad}{tot_t:>{time_w}.6f} s  "
        f"{pct_str:>{pct_w}}  "   # two spaces before Count
    )

    # fixed Count start column
    pad_spaces = max(0, count_start_col - len(left))
    lines.append(left + (" " * pad_spaces) + f"{tot_c:<{count_w}}")

    # children
    kids = _children_of(prefix)
    kids_sorted = sorted(kids, key=lambda ch: _subtree_sum(_subkey(prefix, ch))[0], reverse=True)
    for ch in kids_sorted:
        _print_node(_subkey(prefix, ch), tot_t, level + 1, lines,
                    name_w, time_w, pct_w, count_w,
                    name_indent, col_indent, base_col_levels, count_start_col)
# ---- assumes you already have: _children_of, _subkey, _subtree_sum, and the tic/toc state ----

def _row_left(prefix, parent_total, level,
              name_w, time_w, pct_w,
              name_indent=2, col_indent=2, base_col_levels=2):
    """Build the left chunk (name + time + %), no Count; return string and subtree count."""
    tot_t, tot_c = _subtree_sum(prefix)
    label = prefix.split(":")[-1]

    # visible indents
    name_pad = " " * (name_indent * level)
    col_pad  = " " * (col_indent * (base_col_levels + level))

    # label width so time column stays consistent within this row’s baseline
    avail = max(1, name_w - len(name_pad))
    label_field = f"{label:<{avail}}"

    # % of parent
    if parent_total and parent_total > 0:
        pct = 100.0 * tot_t / parent_total
        pct_str = f"{pct:>{pct_w-1}.1f}%"
    else:
        pct_str = ""

    left = (
        f"{name_pad}{label_field} "
        f"{col_pad}{tot_t:>{time_w}.6f} s  "
        f"{pct_str:>{pct_w}}  "  # two spaces before Count (we'll pad further later)
    )
    return left, tot_c, tot_t

def _collect_rows(prefix, parent_total, level, rows,
                  name_w, time_w, pct_w,
                  name_indent, col_indent, base_col_levels):
    left, cnt, _ = _row_left(prefix, parent_total, level,
                             name_w, time_w, pct_w,
                             name_indent, col_indent, base_col_levels)
    rows.append((left, cnt))
    # children
    kids = _children_of(prefix)
    kids_sorted = sorted(kids, key=lambda ch: _subtree_sum(_subkey(prefix, ch))[0], reverse=True)
    for ch in kids_sorted:
        full = _subkey(prefix, ch)
        _collect_rows(full,
                      _subtree_sum(prefix)[0],   # parent_total = this node's subtree total
                      level + 1, rows,
                      name_w, time_w, pct_w,
                      name_indent, col_indent, base_col_levels)

def print_all_timers(name_width=30, time_width=12, pct_width=8, count_width=7,
                     name_indent=2, col_indent=2, header_col_levels=2,
                     count_align="right"):
    """
    Two-indent scheme; header and rows share the same fixed Count start column.
    """
    if not _time:
        print("(no timers)"); return

    # --- Build header LEFT chunk (no Count yet) with header's column baseline ---
    col_pad_header = " " * (col_indent * header_col_levels)
    header_left = (
        f"{'Timer':<{name_width}} "
        f"{col_pad_header}{'Total Time':>{time_width}}  "
        f"{'     % parent':>{pct_width}}  "
    )

    # --- Collect all rows' left chunks first (no Count yet) ---
    tops = sorted({k.split(':',1)[0] for k in _time.keys()},
                  key=lambda t: _subtree_sum(t)[0], reverse=True)

    rows = []
    for top in tops:
        _collect_rows(top, parent_total=None, level=0, rows=rows,
                      name_w=name_width, time_w=time_width, pct_w=pct_width,
                      name_indent=name_indent, col_indent=col_indent,
                      base_col_levels=header_col_levels)

    # --- Decide fixed Count start column from max of header_left and row lefts ---
    max_left = max([len(header_left)] + [len(L) for (L, _) in rows])
    count_start_col = max_left
    total_line_len  = count_start_col + count_width

    # --- Build header with padding so "Count" aligns with data column ---
    pad_header = " " * max(0, count_start_col - len(header_left))
    if count_align == "right":
        count_head = f"{'Count':>{count_width}}"
    else:  # left
        count_head = f"{'Count':<{count_width}}"
    header = header_left + pad_header + count_head

    # --- Render table ---
    lines = [header, "-" * total_line_len]
    for left, cnt in rows:
        pad = " " * max(0, count_start_col - len(left))
        cnt_field = f"{cnt:>{count_width}d}" if count_align == "right" else f"{cnt:<{count_width}d}"
        lines.append(left + pad + cnt_field)

    # Footer (keep Count blank, but aligned)
    grand_total = sum(_time.values())
    footer_left = (
        f"{'ALL':<{name_width}} "
        f"{col_pad_header}{grand_total:>{time_width}.6f} s  "
        f"{'':>{pct_width}}  "
    )
    pad_footer = " " * max(0, count_start_col - len(footer_left))
    lines.append("-" * total_line_len)
    lines.append(footer_left + pad_footer + " " * count_width)
    print("\n")
    print("\n".join(lines))
