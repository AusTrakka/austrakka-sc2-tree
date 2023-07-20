import pytest


@pytest.mark.parametrize("args, exit_code, expected_stdout", [
    (["run", "--fasta", "../tests/data/test.fasta"], 0, "(100%) done"),
    (["run", "--data", "../tests/data/test.metadata.csv"], 2, ""),
])
def test_command(app_runner, args, exit_code, expected_stdout):
    result = app_runner(*args)
    assert result.exit_code == exit_code
    assert expected_stdout in result.stdout
