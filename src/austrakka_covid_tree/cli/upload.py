import os
from pathlib import Path
from typing import Optional

import typer
from rich import print

app = typer.Typer()


@app.command()
def upload(
    tree: Optional[Path] = typer.Option(None, "--tree", "-t", help="Newick file to upload to AT2", exists=True, dir_okay=False),
    lineage_file: Optional[Path] = typer.Option(
        None, "--lineage", "-l", help="Lineage file to upload to AT2", exists=True, dir_okay=False
    ),
    analysis: Optional[str] = typer.Option("SC2-ANZ-test", "--analysis-id", "-a", help="AT2 analysis ID"),
    pro_forma: Optional[str] = typer.Option("SC2-ANZ-lineage", "--pro-forma", "-p", help="AT2 pro forma ID"),
    last_lineage_file: Optional[Path] = typer.Option(
        None,
        "--last-lineage",
        "-ll",
        help="Previous lineage TSV file previously uploaded to AT2. Any identical rows in this TSV will be filtered out of the lineage file before upload.",
        exists=True,
        dir_okay=False,
    ),
):
    """Upload a tree and lineage to AT2"""
    if lineage_file:
        # get diff of lineage
        if last_lineage_file:
            print("Comparing to previous lineage")
            diff_file = f"lineage-update-{analysis}.tsv"

            # Read the file content only once
            with open(last_lineage_file) as f:
                old_lineage_data = f.readlines()

            with open(lineage_file) as f:
                new_lineage_data = f.readlines()

            # Check headers
            if old_lineage_data[0] != new_lineage_data[0]:
                msg = "Header rows of old and new lineage TSVs do not match"
                raise ValueError(msg)

            # Write the header to diff_file
            with open(diff_file, "w") as f:
                f.write(new_lineage_data[0])

            # Filter and write the rows that are not in old lineage file
            with open(diff_file, "a") as f:
                for line in new_lineage_data[1:]:
                    if line not in old_lineage_data:
                        f.write(line)

            # Count rows in new and diff files
            num_rows_md = len(new_lineage_data) - 1  # Subtract 1 for header
            with open(diff_file) as f:
                num_rows_diff = len(f.readlines()) - 1  # Subtract 1 for header

            print(f"Kept {num_rows_diff} of {num_rows_md} rows")
            lineage_file = diff_file

        # push analysis
        os.system(
            f'austrakka-covid-tree env run at "austrakka metadata append -p {pro_forma} {lineage_file}"'  # noqa: S605
        )

    if tree:
        # push tree
        os.system(f'austrakka-covid-tree env run at "austrakka tree add -a {analysis} {tree}"')  # noqa: S605
