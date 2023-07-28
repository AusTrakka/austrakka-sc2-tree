import shutil

import pytest
from typer.testing import CliRunner

from austrakka_sc2_tree.cli import austrakka_sc2_tree


@pytest.fixture
def runner():
    return CliRunner()

@pytest.fixture
def app_runner(runner, tmp_path):
    data_dir = tmp_path / "tests" / "data"
    print(f"Copying tests/data to {data_dir}")
    shutil.copytree("tests/data", str(data_dir))
    def run_command(*args):
        with runner.isolated_filesystem(temp_dir=str(tmp_path)):
            return runner.invoke(austrakka_sc2_tree.app, args)
    return run_command
