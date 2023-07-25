#!/usr/bin/env python3

import argparse

from Bio import Phylo


def ladderize_tree(input, output):
    # Read the tree from a Newick file
    tree = Phylo.read(input, "newick")

    # Ladderize the tree
    tree.ladderize()

    # Write the ladderized tree back to the same Newick file
    Phylo.write(tree, output, "newick")

def main():
    parser = argparse.ArgumentParser(description="Ladderize a Newick tree.")
    parser.add_argument("input", help="The Newick file to be ladderized.")
    parser.add_argument("output", help="Output file name.")

    args = parser.parse_args()

    ladderize_tree(args.input, args.output)

if __name__ == "__main__":
    main()
