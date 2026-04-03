#!/usr/bin/env python3
"""
Make an ultrametric tree from a rooted input tree.

Algorithm adapted from OrthoFinder's make_ultrametric.py
(https://github.com/davidemms/OrthoFinder/blob/master/tools/make_ultrametric.py)
by David Emms. Reimplemented here as a standalone script using ete3 directly,
without requiring OrthoFinder's internal libraries.

Adjusts branch lengths so that all root-to-leaf distances are equal to the
average root-to-leaf distance, preserving relative topology.

Usage:
    make_ultrametric.py <tree_file> [root_age] > output.nwk

Arguments:
    tree_file   Rooted newick tree (e.g. from OrthoFinder SpeciesTree)
    root_age    Optional: rescale so root-to-tip distance equals this value
                (default: 1.0 — gives relative branch lengths summing to 1)
"""
import sys
import argparse

import numpy as np
from ete3 import Tree


def ave_dist(node):
    """Average distance from node to all its leaf descendants."""
    return float(np.mean([node.get_distance(l) for l in node.get_leaf_names()]))


def make_ultrametric(t, root_age=1.0):
    """
    Adjust branch lengths of tree t in-place so it is ultrametric.
    After this call every leaf is at distance root_age from the root.
    """
    if len(t.get_children()) != 2:
        sys.exit("ERROR: input tree must be rooted (binary root node)")

    # Average root-to-leaf distance in the original tree
    d = ave_dist(t)
    if d == 0:
        sys.exit("ERROR: all branch lengths are zero in the input tree")

    for n in t.traverse('preorder'):
        if n.is_root():
            n.dist = 0.0
            continue
        # distance from root to parent of n
        x = t.get_distance(n) - n.dist
        y = n.dist
        z = 0.0 if n.is_leaf() else ave_dist(n)
        if (y + z) == 0.0:
            n.dist = 0.0
        else:
            n.dist = y * (d - x) / (y + z)

    # Rescale to requested root age
    scale = root_age / d
    for n in t.traverse():
        if not n.is_root():
            n.dist *= scale


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

    t = Tree(args.tree_fn, format=1)
    make_ultrametric(t, root_age=args.root_age)
    # format=5: internal branch lengths + leaf names, no internal node names
    print(t.write(format=5))


if __name__ == "__main__":
    main()
