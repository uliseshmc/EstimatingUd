from click.testing import CliRunner

from estimating_ud.cli import main

COMMANDS = [
    ["align-cds"], ["filter-introns"], ["filter-noncds"], ["human-seq"], ["fit"], ["constraint"],
    ["quality"], ["quality", "masked-databank"], ["quality", "masked-concat"],
    ["quality", "filter-spscore"], ["quality", "fit-spscore"],
]


def test_help_of_every_command():
    runner = CliRunner()
    assert runner.invoke(main, ["--help"]).exit_code == 0
    for cmd in COMMANDS:
        result = runner.invoke(main, [*cmd, "--help"])
        assert result.exit_code == 0, (cmd, result.output)


def test_single_dash_aliases_accepted(data_root):
    result = CliRunner().invoke(main, ["filter-noncds", "-chrm", "22", "-submodel", "singlent", "--data-root", str(data_root)])
    # distalIG and intergenicAR exist in the fixture, the other three regions do not
    assert "does not exist" in str(result.exception) or result.exit_code == 0


def test_data_root_from_env(data_root, monkeypatch):
    monkeypatch.setenv("DATA_HUMCHIMPORANGOR114", str(data_root))
    result = CliRunner().invoke(main, ["human-seq", "--chrm", "22"])
    assert result.exit_code != 0  # fixture lacks intronsAR etc., but the env var was picked up
    assert "Missing option" not in result.output


def test_bad_chromosome_rejected():
    result = CliRunner().invoke(main, ["fit", "--chrm", "23", "--submodel", "trinuc", "--data-root", "."])
    assert result.exit_code == 2
    assert "Invalid value" in result.output
