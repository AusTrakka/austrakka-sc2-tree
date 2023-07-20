from pathlib import Path

from snk.cli import CLI

from austrakka_covid_tree.__about__ import __version__
from austrakka_covid_tree.cli.upload import upload

austrakka_covid_tree = CLI(pipeline_dir_path=Path(__file__).parent.parent)
austrakka_covid_tree.version = __version__
austrakka_covid_tree.register_command(upload)
