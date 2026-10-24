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


def main():
    os.makedirs(OUT, exist_ok=True)
    for fig in (fig_overview, fig_pipeline, fig_storage, fig_baseball, fig_tiling, fig_engines):
        print(f'wrote {fig()}.svg')


if __name__ == '__main__':
    main()
