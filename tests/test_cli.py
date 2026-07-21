import pytest

from trf_pipeline.cli import main


def test_help_smoke():
    with pytest.raises(SystemExit) as exc:
        main(["--help"])
    assert exc.value.code == 0
