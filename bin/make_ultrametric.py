#!/usr/bin/env python3
"""
Make an ultrametric tree from a rooted input tree.

Algorithm adapted from OrthoFinder's make_ultrametric.py
(https://github.com/davidemms/OrthoFinder/blob/master/tools/make_ultrametric.py)
by David Emms. Reimplemented here as a standalone script using only the Python
standard library — no ete3, numpy, or other external dependencies required.

Adjusts branch lengths so that all root-to-leaf distances equal the average
root-to-leaf distance of the input tree, preserving relative topology.

Usage:
    make_ultrametric.py <tree_file> [root_age] > output.nwk

Arguments:
    tree_file   Rooted newick tree (e.g. from OrthoFinder SpeciesTree)
    root_age    Optional: rescale so root-to-tip distance equals this value
                (default: 1.0)
"""
import sys
import argparse


# ---------------------------------------------------------------------------
# Minimal newick parser — no external dependencies
# ---------------------------------------------------------------------------

class Node:
    __slots__ = ('name', 'dist', 'children', 'parent')

    def __init__(self, name='', dist=0.0):
        self.name = name
        self.dist = dist
        self.children = []
        self.parent = None

    def is_leaf(self):
        return not self.children

    def is_root(self):
        return self.parent is None

    def get_leaves(self):
        if self.is_leaf():
            return [self]
        leaves = []
        for child in self.children:
            leaves.extend(child.get_leaves())
        return leaves

    def dist_to_root(self):
        """Sum of branch lengths from this node up to (but not including) root."""
        d, node = 0.0, self
        while node.parent is not None:
            d += node.dist
            node = node.parent
        return d

    def dist_to_descendant(self, desc):
        """Sum of branch lengths from self down to a descendant node."""
        d, node = 0.0, desc
        while node is not self:
            d += node.dist
            node = node.parent
        return d

    def ave_dist_to_leaves(self):
        leaves = self.get_leaves()
        return sum(self.dist_to_descendant(leaf) for leaf in leaves) / len(leaves)

    def traverse_preorder(self):
        yield self
        for child in self.children:
            yield from child.traverse_preorder()


def _find_outer_colon(s):
    """Return index of the last colon at parenthesis depth 0, or -1."""
    depth = 0
    for i in range(len(s) - 1, -1, -1):
        c = s[i]
        if c == ')':
            depth += 1
        elif c == '(':
            depth -= 1
        elif c == ':' and depth == 0:
            return i
    return -1


def _split_children(inner):
    """Split comma-separated newick children respecting nested parentheses."""
    children, depth, buf = [], 0, []
    for c in inner:
        if c == '(':
            depth += 1
            buf.append(c)
        elif c == ')':
            depth -= 1
            buf.append(c)
        elif c == ',' and depth == 0:
            children.append(''.join(buf))
            buf = []
        else:
            buf.append(c)
    if buf:
        children.append(''.join(buf))
    return children


def _parse(s, parent=None):
    node = Node(parent=parent)
    s = s.strip()

    # Strip branch length from the end
    colon = _find_outer_colon(s)
    if colon != -1:
        try:
            node.dist = float(s[colon + 1:])
        except ValueError:
            node.dist = 0.0
        s = s[:colon]

    if s.startswith('('):
        last_paren = s.rfind(')')
        node.name = s[last_paren + 1:].strip()
        for child_s in _split_children(s[1:last_paren]):
            child = _parse(child_s, parent=node)
            node.children.append(child)
    else:
        node.name = s.strip()

    return node


def parse_newick(text):
    """Parse a newick string and return the root Node."""
    return _parse(text.strip().rstrip(';'))


def write_newick(node):
    """Serialise a Node tree back to newick format."""
    if node.is_leaf():
        return f"{node.name}:{node.dist:.10g}"
    children_str = ','.join(write_newick(c) for c in node.children)
    label = node.name or ''
    if node.is_root():
        return f"({children_str}){label};"
    return f"({children_str}){label}:{node.dist:.10g}"


# ---------------------------------------------------------------------------
# Ultrametricisation — same algorithm as OrthoFinder's make_ultrametric.py
# ---------------------------------------------------------------------------

def make_ultrametric(root, root_age=1.0):
    """
    Adjust branch lengths in-place so every leaf is equidistant from the root.
    Uses the average-distance algorithm from OrthoFinder (Emms & Kelly).
    """
    if len(root.children) != 2:
        sys.exit(
            f"ERROR: input tree must be rooted (root has {len(root.children)} "
            "children; expected 2)"
        )

    d = root.ave_dist_to_leaves()
    if d == 0.0:
        sys.exit("ERROR: all branch lengths are zero in the input tree")

    for node in root.traverse_preorder():
        if node.is_root():
            node.dist = 0.0
            continue
        x = node.dist_to_root() - node.dist   # depth of parent (using already-updated ancestors)
        y = node.dist                           # current branch length
        z = 0.0 if node.is_leaf() else node.ave_dist_to_leaves()
        node.dist = 0.0 if (y + z) == 0.0 else y * (d - x) / (y + z)

    scale = root_age / d
    for node in root.traverse_preorder():
        if not node.is_root():
            node.dist *= scale


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Convert a rooted species tree to an ultrametric tree."
    )
    parser.add_argument("tree_fn", help="Input newick tree file (must be rooted)")
    parser.add_argument(
        "root_age",
        type=float,
        nargs="?",
        default=1.0,
        help="Root-to-tip distance in output tree (default: 1.0)",
    )
    args = parser.parse_args()

    with open(args.tree_fn) as fh:
        newick_text = fh.read()

    root = parse_newick(newick_text)
    make_ultrametric(root, root_age=args.root_age)
    print(write_newick(root))


if __name__ == "__main__":
    main()
