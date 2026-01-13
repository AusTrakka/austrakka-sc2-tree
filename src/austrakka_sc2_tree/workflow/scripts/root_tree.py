#!/usr/bin/env python3

import argparse

from Bio import Phylo


def root_tree(input, output, root_on):
    # Read the tree from a Newick file
    tree = Phylo.read(input, "newick")

    # Root the tree on the specified outgroup
    tree.root_with_outgroup({"name": root_on})

    # Write the rooted tree back to the same Newick file
    Phylo.write(tree, output, "newick")

def main():
    parser = argparse.ArgumentParser(description="Ladderize a Newick tree.")
    parser.add_argument("input", help="The Newick file to be ladderized.")
    parser.add_argument("output", help="Output file name.")
    parser.add_argument("root_on", help="The outgroup to root the tree on.")

    args = parser.parse_args()

    root_tree(args.input, args.output, args.root_on)

if __name__ == "__main__":
    main()
