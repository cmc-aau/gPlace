###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

"""Newick reading and writing with pplacer v1.1.alpha19 semantics.

Reading follows ``newick_lexer.mll`` / ``newick_parser.mly``: quoted labels use
backslash escapes (``\\'`` and ``\\\\``), branch lengths go through ``float()``
(correctly rounded, as OCaml's ``float_of_string``), and node ids are assigned in
post-order, which makes them the jplace edge numbers.

Writing follows ``Newick_bark.to_newick_string`` / ``Newick_gtree.to_string``:
``label:%g{edge}``, labels quoted only when they contain special characters.
"""

import re

from gplace.exceptions import GPlaceException


class Tree(object):
    """A rooted tree with pplacer's post-order node ids.

    Attributes
    ----------
    children : list[list[int]]
        Child ids of every node.
    bl : list[float | None]
        Branch length above every node (None when absent in the Newick string).
    label : list[str | None]
        Node labels (leaf names, or internal labels such as support/taxonomy).
    edge_label : list[str | None]
        Optional ``[...]`` edge labels.
    root : int
        Id of the root node (always the last id).
    """

    def __init__(self):
        self.children, self.bl, self.label, self.edge_label = [], [], [], []
        self.root = None

    @property
    def n_nodes(self):
        return len(self.children)

    def leaves(self):
        return [v for v in range(self.n_nodes) if not self.children[v]]

    def branch_lengths_for_placement(self):
        """Branch lengths as pplacer uses them: missing lengths are errors except at the
        root, which gets 0.0 (``Newick_gtree.add_zero_root_bl``)."""
        out = list(self.bl)
        for v, b in enumerate(out):
            if b is None:
                if v == self.root:
                    out[v] = 0.0
                else:
                    raise GPlaceException(f'Branch length missing for node {v} ({self.label[v]}).')
        out[self.root] = 0.0
        return out

    def postorder(self):
        order, stack = [], [(self.root, False)]
        while stack:
            v, done = stack.pop()
            if done or not self.children[v]:
                order.append(v)
                continue
            stack.append((v, True))
            for c in reversed(self.children[v]):
                stack.append((c, False))
        return order


_UNQUOTED = re.compile(r"[^ \t\n()\[\]{}':;,]+")
_QUOTED = re.compile(r"'(?:\\['\\]|[^ '])*'")
_ESCAPES = re.compile(r"\\(['\\])")


def _unquote(s):
    return _ESCAPES.sub(r"\1", s[1:-1])


def _tokens(s):
    """Tokenise like newick_lexer.mll: yields (kind, text)."""
    pos, n = 0, len(s)
    while pos < n:
        c = s[pos]
        if c in ' \t\n\r':
            pos += 1
            continue
        if c in '():;,[]{}':
            yield c, c
            pos += 1
            continue
        if c == "'":
            m = _QUOTED.match(s, pos)
            if not m:
                raise GPlaceException(f'Malformed quoted label at position {pos} of the Newick string.')
            yield 'LABEL', _unquote(m.group(0))
            pos = m.end()
            continue
        m = _UNQUOTED.match(s, pos)
        if not m:
            raise GPlaceException(f'Unexpected character {c!r} at position {pos} of the Newick string.')
        yield 'LABEL', m.group(0)
        pos = m.end()


def parse_newick(text):
    """Parse a Newick string into a :class:`Tree` with pplacer's post-order ids."""
    toks = list(_tokens(text.strip()))
    t = Tree()
    i = 0

    def peek():
        return toks[i][0] if i < len(toks) else None

    def new_node(ch):
        t.children.append(ch); t.bl.append(None); t.label.append(None); t.edge_label.append(None)
        return len(t.children) - 1

    def bark(v):
        nonlocal i
        if peek() == 'LABEL':
            t.label[v] = toks[i][1]; i += 1
        if peek() == ':':
            i += 1
            if peek() != 'LABEL':
                raise GPlaceException('Branch length expected after ":" in the Newick string.')
            try:
                t.bl[v] = float(toks[i][1])
            except ValueError:
                raise GPlaceException(f'Invalid branch length {toks[i][1]!r} in the Newick string.')
            i += 1
        if peek() == '{':                      # node numbers in the input are ignored
            i += 3
        if peek() == '[':
            t.edge_label[v] = toks[i + 1][1]; i += 3

    def subtree():
        nonlocal i
        if peek() == '(':
            i += 1
            kids = [subtree()]
            while peek() == ',':
                i += 1
                kids.append(subtree())
            if peek() != ')':
                raise GPlaceException('Unbalanced parentheses in the Newick string.')
            i += 1
            v = new_node(kids)
        else:
            v = new_node([])
        bark(v)
        return v

    # iterative-safe recursion depth: GTDB trees are deep but well within Python limits
    import sys
    old = sys.getrecursionlimit()
    sys.setrecursionlimit(max(old, 100000))
    try:
        t.root = subtree()
    finally:
        sys.setrecursionlimit(old)
    if peek() == ';':
        i += 1
    if i != len(toks):
        raise GPlaceException('Unexpected text after the end of the Newick tree.')
    return t


_NEEDS_QUOTES = re.compile(r"([ \t\n(){}:;,\[\]])|(['\\])")


def quote_label(s):
    """Newick_bark.quote_label: escape ' and \\ with a backslash; quote when anything special."""
    replaced = False

    def rep(m):
        nonlocal replaced
        replaced = True
        return m.group(1) if m.group(1) is not None else '\\' + m.group(2)

    s2 = _NEEDS_QUOTES.sub(rep, s)
    return f"'{s2}'" if replaced else s2


def node_string(t, v, with_node_numbers=True):
    s = quote_label(t.label[v]) if t.label[v] is not None else ''
    if t.bl[v] is not None:
        s += ':%g' % t.bl[v]
    if with_node_numbers:
        s += '{%d}' % v
    if t.edge_label[v] is not None:
        s += '[%s]' % quote_label(t.edge_label[v])
    return s


def to_newick(t, with_node_numbers=True):
    """Newick_gtree.to_string: children joined by ',', then the node's bark, then ';'."""
    out = {}
    for v in t.postorder():
        if t.children[v]:
            out[v] = '(' + ','.join(out.pop(c) for c in t.children[v]) + ')' + node_string(t, v, with_node_numbers)
        else:
            out[v] = node_string(t, v, with_node_numbers)
    return out[t.root] + ';'
