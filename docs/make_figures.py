#!/usr/bin/env python3
"""Draw the figures of the documentation's "How gPlace works" page (docs/src/_static/figures/*.svg).

    python docs/make_figures.py        (or: make figures, in docs/)

The SVG files are committed: building the documentation does not run this script. Animations are
CSS inside each SVG, so they also play in an <img> tag; every figure's resting state is the complete
picture, which is what readers who ask for reduced motion see. The toy tree's tiles come from
gPlace's tiling code (tiled.tile_layout); the search example is illustrative but follows pplacer's
rules (strike box, strikes, keep at most / keep factor).
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from gplace.engine.tiled import tile_layout  # noqa: E402

OUT = os.path.join(HERE, 'src', '_static', 'figures')

INK, MUTED, FAINT, GRID, PAPER = '#141414', '#6b6b6b', '#d6d3cc', '#ece9e3', '#f7f6f2'
TREE = '#141414'                                            # the reference tree: the logo's branches
ACCENT = '#8f5e12'                                          # gPlace's changes and checks: the logo's amber, darkened for text
BOX = '#e3efdf'                                             # the strike box: the logo's green ring, lighter
ORANGE, ORANGE_T, ORANGE_L = '#d4972f', '#8f5e12', '#fdf3e4'  # queries and placements: the logo's placement path
RED = '#b5452f'
VEC_D, VEC_P, VEC_S = '#5b83b5', '#5f9a5a', '#c26a61'      # distal, proximal, supernode: the rings' blue, green, pink, deeper
TILES = ['#5b83b5', '#c26a61', '#5f9a5a', '#d4972f', '#8e7bbe', '#141414']   # the rings' hues, deeper, then amber and black
FONT = '"Helvetica Neue", Helvetica, Arial, "Liberation Sans", sans-serif'
MONO = 'Menlo, Consolas, "DejaVu Sans Mono", "Liberation Mono", monospace'


def esc(s):
    return str(s).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def num(x):
    return f'{x:.2f}'.rstrip('0').rstrip('.')


def attrs(**kw):
    out = []
    for k, v in kw.items():
        if v is None:
            continue
        out.append(f'{k.rstrip("_").replace("_", "-")}="{num(v) if isinstance(v, float) else v}"')
    return ' '.join(out)


def T(x, y, s, size=13, anchor=None, weight=None, fill=INK, cls=None, raw=False, **kw):
    cls = ' '.join(filter(None, [cls, kw.pop('class_', None)])) or None
    a = attrs(x=float(x), y=float(y), font_size=size, fill=fill, text_anchor=anchor,
              font_weight=weight, class_=cls, **kw)
    return f'<text {a}>{s if raw else esc(s)}</text>'


def span(s, fill=None, weight=None, cls=None, **kw):
    return f'<tspan {attrs(fill=fill, font_weight=weight, class_=cls, **kw)}>{esc(s)}</tspan>'


def R(x, y, w, h, fill='none', **kw):
    return f'<rect {attrs(x=float(x), y=float(y), width=float(w), height=float(h), fill=fill, **kw)}/>'


def C(cx, cy, r, fill='none', **kw):
    return f'<circle {attrs(cx=float(cx), cy=float(cy), r=float(r), fill=fill, **kw)}/>'


def L(x1, y1, x2, y2, stroke=INK, **kw):
    return f'<line {attrs(x1=float(x1), y1=float(y1), x2=float(x2), y2=float(y2), stroke=stroke, **kw)}/>'


def P(d, fill='none', **kw):
    return f'<path {attrs(d=d, fill=fill, **kw)}/>'


def G(items, **kw):
    return f'<g {attrs(**kw)}>' + ''.join(items) + '</g>'


def check(x, y, r=8, color=ACCENT):
    """A check mark in a filled circle centred on (x, y)."""
    return (C(x, y, r, fill=color) +
            P(f'M{num(x - r * .45)} {num(y)} l{num(r * .32)} {num(r * .32)} l{num(r * .58)} {num(-r * .62)}',
              stroke='#ffffff', stroke_width=float(r * .26), stroke_linecap='round', stroke_linejoin='round'))


class Figure:
    def __init__(self, name, width, height, title, desc, period=None):
        self.name, self.w, self.h = name, width, height
        self.title, self.desc, self.period = title, desc, period
        self.body, self.css, self.defs, self.n = [], [], [], 0

    def add(self, *items):
        self.body.extend(items)

    def _cls(self):
        self.n += 1
        return f'a{self.n}'

    def keyframes(self, frames, timing='ease-in-out', extra=''):
        """A looping animation over the figure's period. frames: [(seconds, 'declarations')], each
        value reached at its time; returns the class to put on the element."""
        p = self.period
        frames = sorted(frames, key=lambda f: f[0])
        if frames[0][0] > 0:
            frames.insert(0, (0, frames[0][1]))
        if frames[-1][0] < p:
            frames.append((p, frames[-1][1]))
        c = self._cls()
        steps = ''.join(f'{100 * min(max(t, 0), p) / p:.3f}%{{{d}}}' for t, d in frames)
        self.css.append(f'.{c}{{{extra}animation:k{c} {p}s {timing} infinite}}@keyframes k{c}{{{steps}}}')
        return c

    def show(self, on, off=None, fade=0.3):
        """Visible from `on` to `off` seconds of each loop (default: until the loop's end)."""
        off = self.period - 0.6 if off is None else off
        if on <= 0:
            return self.keyframes([(0, 'opacity:1'), (off, 'opacity:1'), (off + fade, 'opacity:0')])
        return self.keyframes([(0, 'opacity:0'), (on, 'opacity:0'), (on + fade, 'opacity:1'),
                               (off, 'opacity:1'), (off + fade, 'opacity:0')])

    def grow_loop(self, on, axis='X', origin='left center', dur=0.6):
        """Grows from zero at `on` seconds of each loop; fades out at the loop's end."""
        off = self.period - 0.6
        s = f'scale{axis}'
        return self.keyframes([(0, f'opacity:1;transform:{s}(0)'), (on, f'opacity:1;transform:{s}(0)'),
                               (on + dur, f'opacity:1;transform:{s}(1)'), (off, f'opacity:1;transform:{s}(1)'),
                               (off + 0.3, f'opacity:0;transform:{s}(1)')],
                              extra=f'transform-box:fill-box;transform-origin:{origin};')

    def grow(self, delay, axis='X', origin='left center', dur=0.9):
        """Grows from zero once, `delay` seconds after the figure loads."""
        c = self._cls()
        self.css.append(f'.{c}{{transform-box:fill-box;transform-origin:{origin};'
                        f'animation:grow{axis} {dur}s cubic-bezier(.2,.7,.3,1) {delay}s both}}')
        return c

    def appear(self, delay, dur=0.5):
        c = self._cls()
        self.css.append(f'.{c}{{animation:fadein {dur}s ease-out {delay}s both}}')
        return c

    def draw(self, delay, dur=1.6):
        """A path (pathLength="1", stroke-dasharray="1") drawn once."""
        c = self._cls()
        self.css.append(f'.{c}{{animation:draw {dur}s ease-in-out {delay}s both}}')
        return c

    def save(self):
        css = (f'text{{font-family:{FONT}}}.mono{{font-family:{MONO}}}.static-only{{display:none}}'
               '@keyframes growX{from{transform:scaleX(0)}}@keyframes growY{from{transform:scaleY(0)}}'
               '@keyframes fadein{from{opacity:0}}@keyframes draw{from{stroke-dashoffset:1}}'
               + ''.join(self.css) +
               '@media (prefers-reduced-motion:reduce){*{animation:none!important}'
               '.transient{display:none}.static-only{display:inline}}')
        defs = f'<defs>{"".join(self.defs)}</defs>\n' if self.defs else ''
        svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {self.w} {self.h}" width="{self.w}" '
               f'height="{self.h}" role="img" aria-labelledby="title desc">\n'
               f'<title id="title">{esc(self.title)}</title>\n<desc id="desc">{esc(self.desc)}</desc>\n'
               f'<style>{css}</style>\n{defs}' + '\n'.join(self.body) + '\n</svg>\n')
        with open(os.path.join(OUT, self.name + '.svg'), 'w') as fh:
            fh.write(svg)
        return self.name


# Trees ------------------------------------------------------------------------------------------

def parse(nested):
    """children lists and root of a tree written as nested tuples; node ids are in post-order,
    which is also pplacer's edge numbering (the edge above node v is edge v)."""
    children = []

    def walk(t):
        ids = [walk(c) for c in t] if isinstance(t, tuple) else []
        children.append(ids)
        return len(children) - 1
    root = walk(nested)
    return children, root


def cladogram(children):
    """Leaf slot (0, 1, ...), height above the tips and parent of every node."""
    n = len(children)
    pos, height, parent = [0.0] * n, [0] * n, [None] * n
    leaf = 0
    for v in range(n):                                   # post-order
        ch = children[v]
        if ch:
            pos[v] = (pos[ch[0]] + pos[ch[-1]]) / 2
            height[v] = 1 + max(height[c] for c in ch)
            for c in ch:
                parent[c] = v
        else:
            pos[v], leaf = leaf, leaf + 1
    return pos, height, parent


# Figures ----------------------------------------------------------------------------------------

def fig_overview():
    f = Figure('overview', 900, 300, 'gPlace places query genomes on a reference tree',
               'Three aligned query genomes move onto edges of a reference tree. Each placement is '
               'written to the jplace file, with every number identical to pplacer\'s.', period=12)
    children, root = parse(((('a', 'b'), 'c'), (('d', ('e', 'f')), ('g', 'h'))))
    pos, height, parent = cladogram(children)
    top = max(height)
    X = lambda v: 610 - height[v] * 72                   # noqa: E731
    Y = lambda v: 84 + pos[v] * 25.5                     # noqa: E731

    f.add(T(30, 34, 'Query alignment', 15, weight='bold'), T(30, 52, 'aligned genomes to place', 12, fill=MUTED),
          T(300, 34, 'Reference tree', 15, weight='bold'), T(300, 52, 'from the reference package', 12, fill=MUTED),
          T(660, 34, 'jplace file', 15, weight='bold'), T(660, 52, 'pplacer\'s output format', 12, fill=MUTED))
    # the tree
    for v in range(len(children)):
        if v != root:
            p = parent[v]
            f.add(P(f'M{num(X(p))} {num(Y(p))} V{num(Y(v))} H{num(X(v))}', stroke=TREE, stroke_width=3,
                    stroke_linecap='round', stroke_linejoin='round'))
        if not children[v]:
            f.add(C(X(v), Y(v), 3, fill=TREE))
    f.add(L(X(root) - 14, Y(root), X(root), Y(root), stroke=TREE, stroke_width=3, stroke_linecap='round'))
    del top

    # queries: (name, target edge, fraction along the branch, side of the pendant branch)
    targets = [('genome A', 3, .3, -1), ('genome B', 9, .5, 1), ('genome C', 12, .5, 1)]
    seed = 7
    palette = ['#cadbef', '#f8ddba', '#cfe6c8', '#f4ccc6', '#dcd2ea']   # residues: the logo's five rings
    rows = []
    for i, (name, edge, frac, side) in enumerate(targets):
        sy = 104 + 66 * i
        f.add(T(30, sy - 14, name, 11.5, fill=MUTED, cls='mono'))
        for k in range(22):
            seed = (seed * 1103515245 + 12345) % 2 ** 31
            colour = palette[seed % 5] if seed % 9 else '#f1efea'
            f.add(R(30 + k * 7, sy - 7, 6, 14, fill=colour, rx=1))
        p = parent[edge]
        lx, ly = X(p) + frac * (X(edge) - X(p)), Y(edge)
        ex, ey = lx, ly + side * 17
        t0 = 0.8 + 2.3 * i
        t1 = t0 + 1.3
        pend = f.show(t1 - 0.05)
        f.add(L(lx, ly, ex, ey, stroke=ORANGE, stroke_width=2.5, stroke_dasharray='3 3', class_=pend))
        dx, dy = 196 - ex, sy - ey
        move = f.keyframes([(0, f'opacity:0;transform:translate({num(dx)}px,{num(dy)}px)'),
                            (0.3, f'opacity:1;transform:translate({num(dx)}px,{num(dy)}px)'),
                            (t0, f'opacity:1;transform:translate({num(dx)}px,{num(dy)}px)'),
                            (t1, 'opacity:1;transform:translate(0px,0px)'),
                            (f.period - 0.6, 'opacity:1;transform:translate(0px,0px)'),
                            (f.period - 0.3, 'opacity:0;transform:translate(0px,0px)')])
        f.add(C(ex, ey, 5.5, fill=ORANGE, class_=move))
        rows.append((name, edge, t1 + 0.25))

    # the jplace file
    f.add(R(660, 72, 220, 196, fill=PAPER, stroke=GRID, rx=8),
          T(676, 98, '"placements": [', 12, cls='mono', fill=MUTED))
    for i, (name, edge, t) in enumerate(rows):
        y = 126 + 28 * i
        f.add(G([C(682, y - 4, 4.5, fill=ORANGE),
                 T(694, y, f'{name}: edge {edge}', 12, cls='mono', raw=True)], class_=f.show(t)))
    f.add(G([check(684, 222), T(700, 218, 'identical to pplacer', 13, weight='bold', fill=ACCENT),
             T(700, 236, 'every number, digit for digit', 11.5, fill=MUTED)], class_=f.show(8.2)))
    return f.save()


def fig_pipeline():
    f = Figure('pipeline', 900, 300, 'What pplacer computes, step by step',
               'Six steps: pre-mask, model and cache once per run; prefilter, baseball search and '
               'refinement for each query. gPlace changes only how the cache and the prefilter store '
               'their vectors.', period=6)
    cards = [
        ('Pre-mask', 'keep the informative columns', None),
        ('Model', 'WAG + CAT rates; P(t) for every rate', None),
        ('Cache', 'distal and proximal vectors per node', 'gPlace: internal nodes only'),
        ('Prefilter', 'score every edge, rank the edges', 'gPlace: supernodes on the fly'),
        ('Baseball', 'optimise ranked edges until 6 strikes', None),
        ('Keep and refine', 're-optimise the best ≤ 7 at 1e-5', None),
    ]
    W, H = 276, 82
    xs = [20, 312, 604]
    ys = [56, 196]
    f.add(T(20, 40, 'Once per run', 13, weight='bold', fill=MUTED),
          T(20, 180, '↻  For each query', 13, weight='bold', fill=ORANGE_T))
    for i, (title, what, note) in enumerate(cards):
        x, y = xs[i % 3], ys[i // 3]
        f.add(R(x, y, W, H, fill='#ffffff', stroke=FAINT, stroke_width=1.2, rx=8))
        if i >= 3:
            hl = f.show(1.5 * (i - 3), 1.5 * (i - 2) - 0.25, fade=0.25)
            f.add(R(x, y, W, H, fill=ORANGE_L, stroke=ORANGE, stroke_width=2, rx=8, class_=f'{hl} transient'))
        f.add(C(x + 22, y + 24, 12, fill=INK), T(x + 22, y + 28.5, str(i + 1), 13, anchor='middle', weight='bold',
                                                   fill='#ffffff'),
              T(x + 42, y + 29, title, 15, weight='bold'),
              T(x + 16, y + 52, what, 12.5, fill=MUTED),
              T(x + 16, y + 71, note or 'gPlace: unchanged', 12, fill=ACCENT if note else '#9b9288',
                weight='bold' if note else None))
        if i % 3 < 2:
            ax = x + W + 3
            f.add(P(f'M{ax + 3} {y + H / 2 - 6} l7 6 l-7 6', stroke=FAINT, stroke_width=2, stroke_linecap='round',
                    stroke_linejoin='round'))
    # from step 3 down to step 4
    x3 = xs[2] + W / 2
    x4 = xs[0] + W / 2
    f.add(P(f'M{x3} {ys[0] + H + 2} V{ys[0] + H + 22} H{x4} V{ys[1] - 4}', stroke=FAINT, stroke_width=2,
            stroke_linejoin='round'),
          P(f'M{x4 - 6} {ys[1] - 11} l6 7 l6 -7', stroke=FAINT, stroke_width=2, stroke_linecap='round',
            stroke_linejoin='round'))
    return f.save()


def fig_storage():
    f = Figure('storage', 900, 350, 'Likelihood vectors kept by pplacer and by gPlace',
               'The same small tree twice. pplacer keeps distal, proximal and supernode vectors for every '
               'node; gPlace keeps distal and proximal vectors for internal nodes only and rebuilds leaf '
               'vectors when they are needed.', period=9)
    children, root = parse(((('A', 'B'), 'C'), (('D', 'E'), 'F')))
    pos, height, parent = cladogram(children)
    sq, gap = 11, 2.5
    for panel, px in (('pplacer', 20), ('gPlace', 470)):
        X = lambda v: px + 45 + pos[v] * 62                       # noqa: E731
        Y = lambda v: 236 - height[v] * 46                        # noqa: E731
        f.add(T(px, 32, panel, 16, weight='bold'),
              T(px, 51, '3 vectors for every node' if panel == 'pplacer' else
                '2 vectors for each internal node; leaves rebuilt when needed', 12.5, fill=MUTED))
        for v in range(len(children)):
            if v != root:
                p = parent[v]
                f.add(P(f'M{num(X(p))} {num(Y(p))} H{num(X(v))} V{num(Y(v))}', stroke=TREE, stroke_width=2.5,
                        stroke_linejoin='round', stroke_linecap='round'))
        f.add(L(X(root), Y(root) - 12, X(root), Y(root), stroke=TREE, stroke_width=2.5, stroke_linecap='round'))
        leaf_i = 0
        for v in range(len(children)):
            leaf = not children[v]
            if leaf:
                bx, by = X(v) - (3 * sq + 2 * gap) / 2, Y(v) + 8
            else:
                bx, by = X(v) + 7, Y(v) - sq - 6
            kinds = [VEC_D, VEC_P, VEC_S] if panel == 'pplacer' else [VEC_D, VEC_P]
            for j, colour in enumerate(kinds):
                x = bx + j * (sq + gap)
                if panel == 'gPlace' and leaf:
                    f.add(R(x + .75, by + .75, sq - 1.5, sq - 1.5, stroke=colour, stroke_width=1.5,
                            stroke_dasharray='2.2 1.6', rx=1.5))
                    t = 0.8 + 1.25 * leaf_i
                    f.add(R(x, by, sq, sq, fill=colour, rx=1.5, class_=f'{f.show(t, t + 0.9, fade=0.2)} transient'))
                else:
                    f.add(R(x, by, sq, sq, fill=colour, rx=1.5))
            f.add(C(X(v), Y(v), 3.5 if not leaf else 3, fill=TREE))
            leaf_i += leaf
        n_nodes, n_internal = len(children), sum(1 for c in children if c)
        if panel == 'pplacer':
            f.add(T(px, 284, f'{3 * n_nodes} vectors kept ({n_nodes} nodes × 3)', 13, weight='bold'),
                  T(px, 304, '3 × nodes × sites × 20 × 8 bytes', 12, fill=MUTED, cls='mono'))
        else:
            f.add(T(px, 284, f'{2 * n_internal} vectors kept ({n_internal} internal nodes × 2)', 13, weight='bold',
                    fill=ACCENT),
                  T(px, 304, '2 × internal nodes × sites × 20 × 8 bytes', 12, fill=MUTED, cls='mono'))
    f.add(L(450, 20, 450, 310, stroke=GRID, stroke_width=1.5))
    legend = [(VEC_D, 'distal: the tree below the node'), (VEC_P, 'proximal: the rest of the tree'),
              (VEC_S, 'supernode: prefilter (gPlace: 8 edges at a time, then discarded)')]
    x = 20
    for colour, label in legend:
        f.add(R(x, 330, sq, sq, fill=colour, rx=1.5), T(x + 17, 340, label, 12, fill=MUTED))
        x += 17 + len(label) * 6.1 + 26
    return f.save()


def fig_baseball():
    f = Figure('baseball', 900, 334, 'pplacer\'s baseball search for one query',
               'Edges are optimised in prefilter order. An edge more than 3 log-likelihood units below the '
               'best so far is a strike; after 6 strikes the search stops. The best 7 edges within a factor '
               '0.01 of the best are kept and get like-weight ratios.', period=22)
    vals = [-1.5, -3.5, 0, -4.5, -1.7, -7.5, -5.5, -1.3, -11.5, -4.0, -9.5]
    n_slots, x0, x1 = 16, 70, 598
    w = (x1 - x0) / n_slots
    xr = lambda i: x0 + (i - .5) * w                      # noqa: E731
    yv = lambda v: 72 + -v * 18                           # noqa: E731
    ts = [1.0 + 1.0 * i for i in range(len(vals))]
    f.add(T(20, 30, 'Baseball search', 15, weight='bold'),
          T(158, 30, 'one query, pplacer\'s default settings (illustration)', 12.5, fill=MUTED))
    for v in (0, -3, -6, -9, -12):
        f.add(L(x0, yv(v), x1, yv(v), stroke=GRID), T(x0 - 8, yv(v) + 4, f'{v}'.replace('-', '−'), 11,
                                                         anchor='end', fill=MUTED))
    f.add(T(18, 180, 'log-likelihood (best = 0)', 11.5, anchor='middle', fill=MUTED, transform='rotate(-90 18 180)'))
    for i in range(1, n_slots + 1):
        f.add(T(xr(i), 306, str(i), 11, anchor='middle', fill=MUTED))
    f.add(T((x0 + x1) / 2, 324, 'edges in prefilter rank order', 11.5, anchor='middle', fill=MUTED))
    # cursor, best so far and strike box
    for i, t in enumerate(ts, 1):
        f.add(R(xr(i) - w / 2 + 1, 62, w - 2, 230, fill=ORANGE_L, class_=f'{f.show(t, t + 0.95, fade=0.12)} transient'))
    best, strikes, events = -math.inf, 0, []
    for i, (v, t) in enumerate(zip(vals, ts), 1):
        if v > best:
            events.append(('best', i, v, t))
            best = v
        elif v < best - 3:
            strikes += 1
            events.append(('strike', i, v, t, strikes))
        else:
            events.append(('in', i, v, t))
    bests = [(e[3], e[2]) for e in events if e[0] == 'best']
    for k, (t, b) in enumerate(bests):
        last = k == len(bests) - 1
        off = None if last else bests[k + 1][0] + 0.35
        cls = f.show(t + 0.35, off, fade=0.2) + ('' if last else ' transient')
        f.add(G([R(x0, yv(b), x1 - x0, 3 * 18, fill=BOX, opacity=.85),
                 L(x0, yv(b), x1, yv(b), stroke=INK, stroke_width=1.5),
                 T(x0 + 6, yv(b) + 3 * 18 - 5, 'strike box (3 units)', 11, fill=MUTED)], class_=cls))
    # the strike counter
    f.add(T(390, 52, 'strikes', 12, fill=MUTED))
    for k in range(6):
        f.add(C(442 + 16 * k, 48, 5.5, stroke=RED, stroke_width=1.5, fill='#ffffff'))
    for e in events:
        kind, i, v, t = e[:4]
        cls = f.show(t + 0.35, fade=0.2)
        if kind == 'strike':
            f.add(G([C(xr(i), yv(v), 6, fill=RED), T(xr(i) + 11, yv(v) - 8, str(e[4]), 10.5, weight='bold', fill=RED)],
                    class_=cls),
                  C(442 + 16 * (e[4] - 1), 48, 5.5, fill=RED, class_=f.show(t + 0.45, fade=0.2)))
        else:
            f.add(C(xr(i), yv(v), 6, fill=TREE, class_=cls))
            if kind == 'best' and v == 0:
                f.add(T(xr(i), yv(v) - 15, 'best', 10.5, anchor='middle', weight='bold', fill=TREE, class_=cls))
    t_stop = ts[-1] + 1.2
    f.add(G([R(xr(len(vals) + 1) - w / 2 + 1, 62, w * (n_slots - len(vals)) - 2, 230, fill=PAPER, fill_opacity=.8),
             T((xr(len(vals) + 1) + xr(n_slots)) / 2, 172, 'not optimised', 11.5, anchor='middle', fill=MUTED)],
            class_=f.show(t_stop)),
          T(540, 52, '6th strike: stop', 12, weight='bold', fill=RED, class_=f.show(t_stop)))
    # kept edges and their like-weight ratios
    keep = sorted(range(len(vals)), key=lambda i: -vals[i])
    keep = [i for i in keep[:7] if vals[i] >= math.log(0.01) + max(vals)]
    t_keep = t_stop + 1.1
    f.add(G([P(f'M{x1 + 6} {yv(0)} h5 V{yv(math.log(0.01))} h-5', stroke=ORANGE, stroke_width=2),
             T(x1 + 15, yv(-2.3) + 4, 'kept', 11.5, weight='bold', fill=ORANGE_T)], class_=f.show(t_keep)))
    for j, i in enumerate(keep):
        f.add(C(xr(i + 1), yv(vals[i]), 9.5, stroke=ORANGE, stroke_width=2.5, class_=f.show(t_keep + 0.08 * j)))
    lw = [math.exp(vals[i]) for i in keep]
    lw = [x / sum(lw) for x in lw]
    f.add(T(660, 82, 'Like-weight ratios', 14, weight='bold', class_=f.show(t_keep + 0.6)),
          T(660, 100, 'the 7 kept edges, after refinement', 11.5, fill=MUTED, class_=f.show(t_keep + 0.6)))
    for j, (i, r) in enumerate(zip(keep, lw)):
        y = 122 + 24 * j
        t = t_keep + 1.0 + 0.2 * j
        f.add(T(660, y + 12, f'rank {i + 1}', 11.5, fill=MUTED, class_=f.show(t)),
              R(714, y, max(r * 140, 1.5), 16, fill=ORANGE, rx=2, class_=f.grow_loop(t)),
              T(720 + max(r * 140, 1.5), y + 12.5, f'{r:.3f}', 11.5, weight='bold', class_=f.show(t + 0.4)))
    return f.save()


def fig_tiling(tile_edges=12):
    A = (('a', 'b'), ('c', ('d', 'e')))
    B = (('f', 'g'), (('h', 'i'), 'j'))
    Cc = ((('k', 'l'), 'm'), (('n', 'o'), ('p', ('q', 'r'))))
    D = (('s', 't'), (('u', 'v'), ('w', 'x')))
    children, root = parse(((A, B), (Cc, D)))
    pos, height, parent = cladogram(children)
    tiles = tile_layout(children, list(range(len(children))), root, tile_edges)
    roots = [t[0] for t in tiles]
    tile_of = [None] * len(children)
    for v in range(len(children) - 1, -1, -1):           # parents before children
        tile_of[v] = roots.index(v) if v in roots else tile_of[parent[v]]
    k = len(tiles)
    f = Figure('tiling', 900, 410, 'Tiled placement on a small tree',
               f'A 24-leaf tree cut into {k} tiles of at most {tile_edges} edges. A bottom-up pass keeps the '
               'distal vector of each tile root, a top-down pass its proximal vector; placement then holds '
               'one tile in memory at a time.', period=21)
    X = lambda v: 60 + pos[v] * 33.9                      # noqa: E731
    Y = lambda v: 352 - height[v] * 41                    # noqa: E731
    edge = lambda v: f'M{num(X(parent[v]))} {num(Y(parent[v]))} H{num(X(v))} V{num(Y(v))}'  # noqa: E731
    t1 = [0.8 + 1.1 * j for j in range(k)]
    t2 = [t1[-1] + 1.6 + 1.0 * j for j in range(k)]       # top-down: reverse order of the tiles
    t3 = [t2[-1] + 1.6 + 1.2 * j for j in range(k)]
    phases = [
        (t1[0] - 0.4, t2[0] - 0.6, '1 · Bottom-up pass', 'build the tiles in turn, keep the distal vector D of each tile root'),
        (t2[0] - 0.4, t3[0] - 0.6, '2 · Top-down pass', 'keep the proximal vector P of each tile root'),
        (t3[0] - 0.4, None, '3 · Place', 'one tile in memory at a time; its boundary vectors are all it needs'),
    ]
    for on, off, title, what in phases:
        f.add(T(20, 34, span(title, weight='bold') + '   ' + span(what, fill=MUTED), 14, raw=True,
                cls=f'{f.show(on, off, fade=0.25)} transient'))
    f.add(T(20, 34, span(f'{k} tiles of at most {tile_edges} edges', weight='bold') + '   ' +
            span('D, P: the boundary vectors kept for each tile root', fill=MUTED), 14, raw=True, cls='static-only'))
    for v in range(len(children)):
        if v != root:
            f.add(P(edge(v), stroke=FAINT, stroke_width=3, stroke_linejoin='round', stroke_linecap='round'))
    for j, (r, _, _) in enumerate(tiles):
        members = [v for v in range(len(children)) if tile_of[v] == j and v != root]
        colour = TILES[j % len(TILES)]
        path = ' '.join(edge(v) for v in members)
        rev = k - 1 - j
        f.add(P(path, stroke=colour, stroke_width=11, stroke_opacity=.2, stroke_linejoin='round',
                stroke_linecap='round', class_=f'{f.show(t2[rev], t2[rev] + 0.9, fade=0.2)} transient'),
              P(path, stroke=colour, stroke_width=11, stroke_opacity=.25, stroke_linejoin='round',
                stroke_linecap='round', class_=f'{f.show(t3[j], t3[j] + 1.1, fade=0.2)} transient'),
              P(path, stroke=colour, stroke_width=3.2, stroke_linejoin='round', stroke_linecap='round',
                class_=f.show(t1[j])))
    for v in range(len(children)):
        if not children[v]:
            f.add(C(X(v), Y(v), 2.6, fill=TILES[tile_of[v] % len(TILES)], class_=f.show(t1[tile_of[v]])))
    for j, (r, _, _) in enumerate(tiles):
        if r == root:
            continue
        x, y = X(r), Y(r)
        rev = k - 1 - j
        f.add(C(x, y, 5.5, fill='#ffffff', stroke=INK, stroke_width=2),
              G([R(x + 7, y - 36, 15, 14, fill=VEC_D, rx=2), T(x + 14.5, y - 25.5, 'D', 10.5, anchor='middle',
                                                                weight='bold', fill='#ffffff')],
                class_=f.show(t1[j] + 0.6)),
              G([R(x + 7, y - 20, 15, 14, fill=VEC_P, rx=2), T(x + 14.5, y - 9.5, 'P', 10.5, anchor='middle',
                                                                weight='bold', fill='#ffffff')],
                class_=f.show(t2[rev] + 0.4)))
    x = 20
    for j, (r, size, _) in enumerate(tiles):
        f.add(R(x, 384, 12, 12, fill=TILES[j % len(TILES)], rx=2),
              T(x + 17, 394, f'tile {j + 1}: {size - (1 if r == root else 0)} edges', 12, fill=MUTED))
        x += 130
    print(f'  tiling: {k} tiles, sizes {[t[1] for t in tiles]}')
    return f.save()


def fig_engines():
    f = Figure('engines', 900, 250, 'The three engine variants compute the same bits',
               'avx2 computes 4 results per instruction, neon 2 and portable 1. All three perform the same '
               'operations in the same order for every result, so their results are identical.', period=8)
    rows = [('avx2', 'x86-64 with AVX2 + FMA', 4), ('neon', 'ARM64', 2), ('portable', 'any 64-bit CPU', 1)]
    tick, start = 0.28, 0.8
    f.add(T(20, 30, '16 results computed by each engine variant', 15, weight='bold'),
          T(20, 49, 'one box per result; boxes filled together come from one instruction', 12.5, fill=MUTED))
    for r, (name, cpu, lanes) in enumerate(rows):
        y = 78 + 58 * r
        f.add(T(20, y + 16, name, 15, weight='bold', cls='mono'), T(20, y + 34, cpu, 12, fill=MUTED))
        x = 210
        for i in range(16):
            if i and i % lanes == 0:
                x += 9
            t = start + tick * (i // lanes)
            f.add(R(x, y, 22, 22, fill='#ffffff', stroke=FAINT, stroke_width=1.2, rx=3),
                  R(x, y, 22, 22, fill=ORANGE, rx=3, class_=f.show(t, fade=0.12)))
            x += 25
        done = start + tick * (16 // lanes) + 0.2
        f.add(G([check(772, y + 11), T(788, y + 16, 'identical bits', 12.5, weight='bold', fill=ACCENT)],
                class_=f.show(done)))
    return f.save()


# Exponent pull ----------------------------------------------------------------------------------

EXP_FILL, EXP_LINE = ORANGE_L, ORANGE                      # what a pull changes: the exponents
MANT_FILL, MANT_LINE = '#e8eff8', VEC_D                    # what it keeps: the mantissa bits


def pow2(e, size=9.5):
    """2 to the power e (raw SVG text; the exponent is raised, so it must end the text)."""
    return '2' + f'<tspan dy="-6" font-size="{size}">{str(e).replace("-", "−")}</tspan>'


def sci(x):
    """x as 3.1e−9."""
    m, e = f'{x:.1e}'.split('e')
    return f'{m}e{int(e)}'.replace('-', '−')


def pull_site():
    """Four values of one site and pplacer's pull of that site (Glv.perhaps_pull_exponent): for each
    value, (value, binary exponent, the 52 mantissa bits in hex) before and after; and the power of
    two the site is multiplied by."""
    import struct
    vals = [3.1e-9, 1.2e-10, 7.4e-12, 2.0e-9]
    _, twoexp = math.frexp(max(vals))
    rows = []
    for v in vals:
        w = v * 2.0 ** -twoexp
        row = []
        for x in (v, w):
            b = struct.unpack('>Q', struct.pack('>d', x))[0]
            row.append((x, ((b >> 52) & 0x7ff) - 1023, f'{b & ((1 << 52) - 1):013x}'))
        assert row[0][2] == row[1][2]                      # the pull keeps every mantissa bit
        rows.append(row)
    return rows, -twoexp


def fig_pull():
    rows, k = pull_site()
    f = Figure('pull', 900, 352, 'One exponent pull, bit by bit',
               f'Four values of one site before and after a pull. Every value is multiplied by 2 to the power '
               f'{k}: its exponent changes, its 52 mantissa bits are identical, and the site exponent goes from 0 '
               f'to -{k}, so every value times 2 to the power e is unchanged.', period=10)
    f.add(T(20, 30, 'One site of a likelihood vector, before and after a pull', 15, weight='bold'),
          T(20, 49, '4 of its 20 values; each double is drawn as its binary exponent and its 52 mantissa bits',
            12.5, fill=MUTED))
    t_row = [1.4 + 0.45 * i for i in range(len(rows))]
    for side, dx in ((0, 0), (1, 456)):
        f.add(T(20 + dx, 82, 'Before the pull' if side == 0 else 'After the pull', 13.5, weight='bold'),
              T(20 + dx, 102, 'value', 11.5, fill=MUTED),
              T(152 + dx, 102, 'exponent', 11.5, anchor='middle', fill=MUTED),
              T(304 + dx, 102, 'mantissa bits', 11.5, anchor='middle', fill=MUTED))
        for i, row in enumerate(rows):
            x, e, m = row[side]
            y = 112 + 40 * i
            items = [T(20 + dx, y + 20, sci(x) if side == 0 else f'{x:.4f}', 13, cls='mono'),
                     R(108 + dx, y, 88, 30, fill=EXP_FILL, stroke=EXP_LINE, stroke_width=1.2, rx=5),
                     T(152 + dx, y + 20, pow2(e), 13, anchor='middle', raw=True),
                     R(204 + dx, y, 200, 30, fill=MANT_FILL, stroke=MANT_LINE, stroke_width=1.2, rx=5),
                     T(304 + dx, y + 19.5, m, 12.5, anchor='middle', cls='mono')]
            f.add(*(items if side == 0 else [G(items, class_=f.show(t_row[i]))]))
        e_box = [R(108 + dx, 280, 296, 30, fill=EXP_FILL, stroke=EXP_LINE, stroke_width=1.2, rx=5),
                 T(256 + dx, 300, f'site exponent e = {0 if side == 0 else "−" + str(k)}', 13, anchor='middle',
                   weight='bold')]
        f.add(*(e_box if side == 0 else [G(e_box, class_=f.show(t_row[-1] + 0.7))]))
    f.add(G([T(440, 176, '× ' + pow2(k), 13, anchor='middle', weight='bold', raw=True),
             P('M414 187 H462', stroke=INK, stroke_width=1.6),
             P('M456 181 l6 6 l-6 6', stroke=INK, stroke_width=1.6, stroke_linecap='round', stroke_linejoin='round')],
            class_=f.show(0.6)))
    t_end = t_row[-1] + 1.6
    f.add(R(20, 326, 13, 13, fill=EXP_FILL, stroke=EXP_LINE, stroke_width=1.2, rx=3),
          T(40, 337, f'changes: every exponent +{k}, the site exponent e −{k}', 12.5, fill=MUTED),
          G([check(484, 332, 7.5), T(498, 337, 'identical: the 52 mantissa bits, so nothing is rounded', 12.5,
                                      weight='bold', fill=ACCENT)], class_=f.show(t_end)))
    return f.save()


def fig_pull_order():
    f = Figure('pull_order', 900, 334, 'When the final exponent pull happens, in pplacer and in gPlace',
               'pplacer builds every vector, leaves included, then pulls every stored vector, then places. '
               'gPlace builds the internal-node vectors, stores them as built, and pulls a copy when a vector '
               'is used; a leaf vector is rebuilt first, from the inputs pplacer used when it built it. The '
               'vectors used are bit-identical.', period=9)
    xs, W, H = [130, 390, 650], 230, 72
    for x, head in zip(xs, ['1 · Build the cache', '2 · Final pull (−10)', '3 · Place the queries']):
        f.add(T(x, 50, head, 13, weight='bold', fill=MUTED))
    lanes = [('pplacer', 70, [('Build all vectors', 'every node, leaves too', False),
                              ('Pull every vector', 'the cache is overwritten', True),
                              ('Place', 'reads the pulled vectors', False)]),
             ('gPlace', 200, [('Build vectors', 'internal nodes only', False),
                              ('Store them as built', 'no pull here', False),
                              ('Copy, then pull', 'a leaf: rebuilt first', True)])]
    for li, (name, y, cards) in enumerate(lanes):
        f.add(T(20, y + H / 2 + 5, name, 15, weight='bold'))
        for j, (title, what, pull) in enumerate(cards):
            x = xs[j]
            f.add(R(x, y, W, H, fill=EXP_FILL if pull else '#ffffff', stroke=EXP_LINE if pull else FAINT,
                    stroke_width=1.2, rx=8),
                  T(x + 16, y + 30, title, 14.5, weight='bold', fill=ORANGE_T if pull else INK),
                  T(x + 16, y + 53, what, 12.5, fill=MUTED))
            t = 0.6 + 1.0 * (3 * li + j)
            f.add(R(x, y, W, H, stroke=ORANGE, stroke_width=2.5, rx=8,
                    class_=f'{f.show(t, t + 0.95, fade=0.2)} transient'))
            if j < 2:
                ax = x + W + 9
                f.add(P(f'M{ax} {y + H / 2 - 6} l7 6 l-7 6', stroke=FAINT, stroke_width=2, stroke_linecap='round',
                        stroke_linejoin='round'))
    cx1, cx3 = xs[0] + W / 2, xs[2] + W / 2
    f.add(L(cx1, 146, cx1, 196, stroke=MUTED, stroke_width=1.2, stroke_dasharray='4 3'),
          T(cx1 + 9, 175, 'same operations', 12, fill=MUTED),
          L(cx3, 146, cx3, 196, stroke=MUTED, stroke_width=1.2, stroke_dasharray='4 3'),
          G([check(cx3 + 17, 171, 7.5), T(cx3 + 31, 175.5, 'same bits', 12.5, weight='bold', fill=ACCENT)],
            class_=f.show(6.8)))
    cx2 = xs[1] + W / 2
    f.add(P(f'M{cx2} 276 V298 H{cx3} V281', stroke=ACCENT, stroke_width=1.6, stroke_linejoin='round'),
          P(f'M{cx3 - 6} 287 l6 -7 l6 7', stroke=ACCENT, stroke_width=1.6, stroke_linecap='round',
            stroke_linejoin='round'),
          T((cx2 + cx3) / 2, 320, 'a rebuilt leaf vector gets exactly the inputs pplacer used in step 1', 12.5,
            anchor='middle', fill=MUTED))
    return f.save()


# Tiled placement, step by step ------------------------------------------------------------------

TILE_ABC = [('A', TILES[0]), ('B', TILES[2]), ('C', TILES[4])]
GREY_CELL = '#e4e1da'
# prefilter scores of one query on a 30-edge tree in three tiles of 10 edges (edges 1-10: tile A, ...);
# edges 13 and 27 tie
PREF_SCORES = {1: -412, 2: -407, 3: -419, 4: -401, 5: -416, 6: -425, 7: -403, 8: -430, 9: -414, 10: -421,
               11: -413, 12: -418, 13: -409, 14: -422, 15: -404, 16: -427, 17: -431, 18: -410, 19: -424,
               20: -433, 21: -415, 22: -406, 23: -420, 24: -426, 25: -429, 26: -417, 27: -409, 28: -432,
               29: -423, 30: -435}


def best_edges(tiles_seen, k=8):
    """The query's first k edges in pplacer's ranking (score descending, ties by ascending edge
    number) among the tiles seen so far."""
    return sorted((e for e in PREF_SCORES if (e - 1) // 10 in tiles_seen), key=lambda e: (-PREF_SCORES[e], e))[:k]


def fig_tiled_prefilter():
    t = [0.6, 4.0, 7.4]
    f = Figure('tiled_prefilter', 900, 352, 'The prefilter of tiled placement, tile by tile',
               'One query is scored against a 30-edge tree in three tiles, one tile in memory at a time. After '
               'each tile its best edges are merged into a short list; after the last tile the list is the '
               'start of pplacer\'s ranking of all the edges. Edges 13 and 27 have equal scores, so the lower '
               'edge number comes first, as in pplacer.', period=12)
    f.add(T(20, 30, 'Prefilter, tile by tile', 15, weight='bold'),
          T(206, 30, 'one query shown; every query of a batch (up to 20,000) is scored at once', 12.5, fill=MUTED))
    lists = [best_edges(range(k + 1)) for k in range(3)]
    for k, (name, colour) in enumerate(TILE_ABC):
        x = 20 + 295 * k
        on, off = t[k], (t[k + 1] if k < 2 else None)
        f.add(R(x, 66, 270, 164, fill=PAPER, stroke=colour, stroke_width=2, rx=8,
                class_=f'{f.show(on, off - 0.3 if off else None, fade=0.25)} transient'),
              R(x, 66, 270, 164, stroke=FAINT, stroke_width=1.2, rx=8),
              T(x + 14, 88, f'tile {name}', 13.5, weight='bold', fill=colour),
              T(x + 256, 88, 'in memory', 12, anchor='end', weight='bold', fill=ACCENT,
                class_=f'{f.show(on, off - 0.3 if off else None, fade=0.25)} transient'))
        for j in range(10):
            e = 10 * k + j + 1
            h = 10 + (PREF_SCORES[e] + 440) * 2.5
            bx = x + 22 + 24 * j
            f.add(R(bx, 210 - h, 14, h, fill=colour, opacity=.28, rx=2),
                  T(bx + 7, 224, str(e), 10.5, anchor='middle', fill=MUTED))
            member = [p for p in range(3) if e in lists[p]]
            if member:
                last = member[-1]
                cls = f.show(t[member[0]], t[last + 1] - 0.3 if last < 2 else None, fade=0.25)
                f.add(R(bx, 210 - h, 14, h, fill=colour, rx=2, class_=cls + ('' if last == 2 else ' transient')))
    f.add(T(20, 262, 'The query\'s best edges so far', 13, weight='bold'),
          T(880, 262, '8 shown; gPlace keeps 164', 12, anchor='end', fill=MUTED))
    for p, edges in enumerate(lists):
        final = p == 2
        cls = f.show(t[p], None if final else t[p + 1] - 0.3, fade=0.25)
        items = [T(244, 262, 'after the last tile: the start of pplacer\'s ranking' if final else
                   f'after tile {TILE_ABC[p][0]}', 12.5, weight='bold', fill=ACCENT if final else MUTED)]
        for i, e in enumerate(edges):
            x = 20 + 108 * i
            colour = TILE_ABC[(e - 1) // 10][1]
            items += [R(x, 274, 100, 48, fill=colour, rx=6),
                      T(x + 50, 295, f'edge {e}', 13, anchor='middle', weight='bold', fill='#ffffff'),
                      T(x + 50, 313, f'−{-PREF_SCORES[e]}', 12, anchor='middle', fill='#ffffff')]
        if final:
            i13, i27 = edges.index(13), edges.index(27)
            x1, x2 = 20 + 108 * i13 + 10, 20 + 108 * i27 + 90
            items += [P(f'M{x1} 328 v6 H{x2} v-6', stroke=MUTED, stroke_width=1.2),
                      T((x1 + x2) / 2, 348, 'equal scores: lower edge number first', 11.5, anchor='middle',
                        fill=MUTED)]
        f.add(G(items, class_=cls + ('' if final else ' transient')))
    return f.save()


# Each query's tiles along its ranking, and the first-stage log-likelihoods of the edges pplacer's
# search reaches (illustrative values; edges optimised past the stop get no value)
ROUND_QUERIES = [
    ('query 1', 'AABCABCBCABBCAAC', [-1000, -1003.4, -1006, -1005.1, -1008, -1007, -1009]),
    ('query 2', 'CBCABCABACBACABC', [-1200, -1198.5, -1199.4, -1205.2, -1200.9, -1206, -1204.8, -1196,
                                     -1203.3, -1197.9, -1208, -1210]),
    ('query 3', 'BACBACABACCBABCA', [-800, -804.9, -801, -806, -807, -805.5, -809, -808]),
]


def replay(lls, max_strikes=6, strike_box=3.0):
    """pplacer's baseball rule on first-stage log-likelihoods in rank order: a mark per edge and the
    number of edges the search uses."""
    best, strikes, marks = -math.inf, 0, []
    for i, v in enumerate(lls):
        if v > best:
            best, mark = v, 'best'
        elif v < best - strike_box:
            strikes += 1
            mark = f'×{strikes}'
        else:
            mark = 'ok'
        marks.append(mark)
        if mark[0] == '×' and strikes >= max_strikes:
            return marks, i + 1
    raise ValueError('the search does not stop')


def fig_tiled_rounds(round_size=8):
    f = Figure('tiled_rounds', 900, 400, 'The baseball search of tiled placement, in rounds',
               'Three queries\' ranked edges, coloured by tile. Each round optimises the next 8 ranked edges of '
               'every query still searching, tile by tile. pplacer\'s stopping rule is then replayed in rank '
               'order: query 1 stops at rank 7, so its rank 8 is discarded; query 3 stops at rank 8; query 2 '
               'needs a second round and stops at rank 12, so its ranks 13 to 16 are discarded. The kept edges '
               'are refined.', period=22)
    tile_colour = dict(TILE_ABC)
    X = lambda r: 112 + 44 * (r - 1)                       # noqa: E731
    f.add(T(20, 24, 'Baseball in rounds', 15, weight='bold'),
          T(198, 24, f'the next {round_size} ranked edges of every query still searching, run tile by tile; '
            'then pplacer\'s rule is replayed', 12.5, fill=MUTED))
    f.add(P('M112 64 v-6 H458 v6', stroke=MUTED, stroke_width=1.2),
          T(285, 48, 'round 1', 12.5, anchor='middle', weight='bold'),
          P('M464 64 v-6 H810 v6', stroke=MUTED, stroke_width=1.2),
          T(637, 48, 'round 2', 12.5, anchor='middle', weight='bold'))
    for r in range(1, 17):
        f.add(T(X(r) + 19, 80, str(r), 11, anchor='middle', fill=MUTED))
    t_fill = {1: [2.0 + 1.3 * k for k in range(3)], 2: [10.4 + 1.0 * k for k in range(3)]}
    t_queue = {1: 0.6, 2: 9.4}
    t_mark = {1: 6.4, 2: 13.6}
    t_keep = 16.6
    for q, (name, tiles, lls) in enumerate(ROUND_QUERIES):
        y = 92 + 72 * q
        marks, used = replay(lls)
        rounds = -(-used // round_size)
        best = max(lls[:used])
        kept = sorted(range(used), key=lambda i: -lls[i])[:7]
        kept = {i + 1 for i in kept if lls[i] >= best + math.log(0.01)}
        f.add(T(20, y + 24, name, 13, weight='bold'))
        for r in range(1, 17):
            x, tile = X(r), tiles[r - 1]
            colour = tile_colour[tile]
            rnd = -(-r // round_size)
            f.add(R(x, y, 38, 38, fill='#ffffff', stroke=FAINT, stroke_width=1.2, rx=5),
                  T(x + 19, y + 24, tile, 12.5, anchor='middle', fill=FAINT))
            if rnd > rounds:
                continue
            tf = t_fill[rnd]['ABC'.index(tile)]
            tm = t_mark[rnd] + 0.25 * ((r - 1) % round_size)
            extra = r > used
            t_extra = t_mark[rnd] + 0.25 * round_size
            f.add(R(x + .8, y + .8, 36.4, 36.4, stroke=colour, stroke_width=1.6, stroke_dasharray='3 2', rx=5,
                    class_=f'{f.show(t_queue[rnd], tf, fade=0.2)} transient'),
                  G([R(x, y, 38, 38, fill=colour, rx=5),
                     T(x + 19, y + 24, tile, 12.5, anchor='middle', weight='bold', fill='#ffffff')],
                    class_=f.show(tf, t_extra if extra else None, fade=0.25) + (' transient' if extra else '')))
            if extra:
                f.add(G([R(x, y, 38, 38, fill=GREY_CELL, rx=5),
                         T(x + 19, y + 24, tile, 12.5, anchor='middle', fill=MUTED),
                         T(x + 19, y + 55, 'extra', 11.5, anchor='middle', fill=MUTED)], class_=f.show(t_extra)))
                continue
            mark = marks[r - 1]
            colour_m, weight = (INK, 'bold') if mark == 'best' else (RED, 'bold') if mark[0] == '×' else (MUTED, None)
            f.add(T(x + 19, y + 55, mark, 11.5, anchor='middle', weight=weight, fill=colour_m, class_=f.show(tm)))
            if r in kept:
                f.add(R(x - 3, y - 3, 44, 44, stroke=ORANGE, stroke_width=2.5, rx=7,
                        class_=f.show(t_keep + 0.12 * r)))
    f.add(T(20, 344, 'tile in memory', 12, fill=MUTED))
    chip_times = {k: [t_fill[1][k], t_fill[2][k], t_keep + 0.9 * k] for k in range(3)}
    for k, (name, colour) in enumerate(TILE_ABC):
        x = 112 + 96 * k
        f.add(R(x, 326, 86, 28, fill='#ffffff', stroke=colour, stroke_width=1.6, rx=6),
              T(x + 43, 345, f'tile {name}', 12.5, anchor='middle', weight='bold', fill=colour))
        for tc in chip_times[k]:
            f.add(G([R(x, 326, 86, 28, fill=colour, rx=6),
                     T(x + 43, 345, f'tile {name}', 12.5, anchor='middle', weight='bold', fill='#ffffff')],
                    class_=f'{f.show(tc, tc + 0.9, fade=0.2)} transient'))
    f.add(T(412, 345, 'each tile is rebuilt once per round, for every query', 12.5, fill=MUTED))
    f.add(T(20, 386, span('best', weight='bold') + span(': new best so far', fill=MUTED), 12, raw=True),
          T(170, 386, span('×n', weight='bold', fill=RED) + span(': strike n, more than 3 below the best', fill=MUTED),
            12, raw=True),
          T(440, 386, span('extra', weight='bold', fill=MUTED) + span(': optimised, then discarded', fill=MUTED), 12,
            raw=True),
          R(664, 374, 16, 16, stroke=ORANGE, stroke_width=2.5, rx=4),
          T(688, 386, 'kept: refined at 1e-5', 12, fill=MUTED))
    return f.save()


def fig_ranked_edges(n_cand=164, max_pitches=40):
    used = max_pitches + 1
    f = Figure('ranked_edges', 900, 206, 'How many ranked edges each query keeps',
               f'Each query keeps {n_cand} ranked edges. pplacer\'s search uses at most {used} successful '
               f'optimisations; failed ones are skipped, so {n_cand - used} spare edges remain. A query that uses '
               f'all {n_cand} is ranked again with {4 * n_cand} edges, and its first {n_cand} are checked to be '
               'unchanged.')
    W = 860 / n_cand
    f.add(T(20, 30, f'Each query\'s {n_cand} ranked edges (4 × {used})', 15, weight='bold'),
          R(20, 44, used * W, 34, fill=VEC_D, rx=5),
          T(20 + used * W / 2, 66, str(used), 13.5, anchor='middle', weight='bold', fill='#ffffff'),
          R(23 + used * W, 44, 857 - used * W, 34, fill=GRID, rx=5),
          T(23 + used * W + (857 - used * W) / 2, 66, f'{n_cand - used} spare', 13.5, anchor='middle', weight='bold'),
          T(20, 98, f'at most {used} successful optimisations', 12.5, weight='bold'),
          T(20, 115, f'pplacer: --max_pitches {max_pitches}, plus the last one', 12, fill=MUTED),
          T(880, 98, 'room for failed optimisations', 12.5, anchor='end', weight='bold'),
          T(880, 115, 'pplacer skips them: they do not count', 12, anchor='end', fill=MUTED))
    W2 = 860 / (4 * n_cand)
    f.add(T(20, 152, f'A query that uses all {n_cand} without stopping is ranked again with 4 × {n_cand} = '
            f'{4 * n_cand} edges', 13, weight='bold'),
          R(20, 164, n_cand * W2, 30, fill=GRID, rx=5),
          T(20 + n_cand * W2 / 2, 184, f'first {n_cand}: checked unchanged', 12, anchor='middle'),
          R(23 + n_cand * W2, 164, 857 - n_cand * W2, 30, fill=VEC_D, rx=5),
          T(23 + n_cand * W2 + (857 - n_cand * W2) / 2, 184, f'{3 * n_cand} more ranked edges: the search goes on',
            12, anchor='middle', weight='bold', fill='#ffffff'))
    return f.save()


def main():
    os.makedirs(OUT, exist_ok=True)
    for fig in (fig_overview, fig_pipeline, fig_storage, fig_baseball, fig_tiling, fig_engines, fig_pull,
                fig_pull_order, fig_tiled_prefilter, fig_tiled_rounds, fig_ranked_edges):
        print(f'wrote {fig()}.svg')


if __name__ == '__main__':
    main()
