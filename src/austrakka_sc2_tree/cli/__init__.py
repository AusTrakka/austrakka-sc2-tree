from pathlib import Path

from snk_cli import CLI

from austrakka_sc2_tree.__about__ import __version__
from austrakka_sc2_tree.cli.upload import upload

austrakka_sc2_tree = CLI(Path(__file__).parent.parent)
austrakka_sc2_tree.workflow.name = "austrakka-sc2-tree"
austrakka_sc2_tree.version = __version__
austrakka_sc2_tree.register_command(upload)
