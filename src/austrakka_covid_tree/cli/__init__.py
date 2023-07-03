from pathlib import Path
from typing import Optional

import typer
from snk.cli import CLI

from austrakka_covid_tree.__about__ import __version__

austrakka_covid_tree = CLI(pipeline_dir_path=Path(__file__).parent.parent)

def _print_pipline_version(ctx: typer.Context, value: bool):
    if value:
        typer.echo(__version__)
        raise typer.Exit()

def callback(
        ctx: typer.Context,
        version: Optional[bool] = typer.Option(  # noqa: B008
            None,
            "-v",
            "--version",
            help="Show the pipeline version.",
            is_eager=True,
            callback=_print_pipline_version,
            show_default=False,
        )
    ):
    if ctx.invoked_subcommand is None:
        typer.echo(f"{ctx.get_help()}")

callback.__doc__ = austrakka_covid_tree.logo

austrakka_covid_tree.register_callback(
    callback,
    invoke_without_command=True,
    context_settings={"help_option_names": ["-h", "--help"]},
)
