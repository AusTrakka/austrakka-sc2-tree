#!/usr/bin/env python3

import argparse
import sys

from phytest import Data, Tree, main


def test_for_long_branches(tree: Tree):
    """
    We expect that the tree is rooted.
    """
    tree.assert_branch_lengths(min=0, max=100)

def test_tree_is_not_rooted(tree: Tree):
    """
    We expect that the tree is not rooted.
    """
    assert not tree.rooted  # noqa: S101

def test_all_leaves_are_unique(tree: Tree):
    """
    We expect that all leaves are unique.
    """
    tree.assert_unique_tips()

def test_number_of_tips(tree: Tree, data: Data):
    """
    We expect that the tree has less tips then the data.
    """
    tree.assert_number_of_tips(max=len(data))

def test_validate_data(data: Data):
    """
    We expect that the tree passes qc.
    """
    data.assert_columns(["Seq_ID", "Coverage", "Lineage", "Lineage_family", "Lineage_full", "Lineage_note", "QC"])
    data.assert_range("Coverage", min=0, max=100)
    data.assert_match("Lineage", r"^[A-Z0-9\.]+$|Unassigned")
    data_nonan = Data(data.dropna())
    data_nonan.assert_match("Lineage_family", r"^[A-Z0-9\.]+$|Recombinant|Unassigned")
    data_nonan.assert_match("Lineage_full", r"^[A-Z0-9\.]+$|Unassigned")
    data.assert_match("Lineage_note", r"^nextclade")
    data.assert_values("QC", ["PASS", "FAIL"])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run phytests on a pipeline results.")
    parser.add_argument("--tree")
    parser.add_argument("--data")
    parser.add_argument("--report")
    parser.add_argument("--cores")

    args = parser.parse_args()
    res = main(
        tree=args.tree,
        data=args.data,
        report=args.report,
        cores=args.cores,
    )
    sys.exit(res)
