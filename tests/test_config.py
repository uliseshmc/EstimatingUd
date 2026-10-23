import pytest

from estimating_ud import config


def test_env_var_wins(monkeypatch, tmp_path):
    monkeypatch.setenv(config.DATA_ENV, str(tmp_path))
    assert config.data_root() == tmp_path


def test_dotenv_fallback(monkeypatch, tmp_path):
    monkeypatch.delenv(config.DATA_ENV, raising=False)
    dotenv = tmp_path / ".env"
    dotenv.write_text(f"# comment\nexport {config.DATA_ENV}='{tmp_path / 'data'}'\nOTHER=1\n")
    monkeypatch.setattr(config, "dotenv_path", lambda: dotenv)
    assert config.data_root() == tmp_path / "data"


def test_missing_everywhere_raises(monkeypatch, tmp_path):
    monkeypatch.delenv(config.DATA_ENV, raising=False)
    monkeypatch.setattr(config, "dotenv_path", lambda: tmp_path / "absent.env")
    with pytest.raises(RuntimeError, match=config.DATA_ENV):
        config.data_root()


def test_override_beats_both(monkeypatch, tmp_path):
    monkeypatch.setenv(config.DATA_ENV, "/somewhere/else")
    assert config.data_root(tmp_path) == tmp_path
