"""Entry point for ``uv run python -m tools.track2``."""

from __future__ import annotations

from pathlib import Path

import typer

from tools.track2.drugs_to_avoid import write as write_drugs_to_avoid
from tools.track2.ranking import write as write_ranking

app = typer.Typer(add_completion=False, help="Regenerate the Track 2 generated documents.")

_REPO_ROOT = Path(__file__).resolve().parents[2]


@app.callback()
def main() -> None:
    """Present so that commands stay named subcommands as more are added.

    Typer collapses a single-command app into a bare invocation, which would have
    made `python -m tools.track2 drugs-to-avoid` an error while that was the only
    generator and a working command the day a second one landed. `ranking` is that
    second generator; the shape was pinned before it arrived, so the justfile recipe
    and every document citing the command stayed correct across the change.
    """


@app.command("drugs-to-avoid")
def drugs_to_avoid() -> None:
    """Rewrite docs/track2-drugs-to-avoid.md from a live triage run."""
    typer.echo(f"wrote {write_drugs_to_avoid(_REPO_ROOT).relative_to(_REPO_ROOT).as_posix()}")


@app.command("ranking")
def ranking() -> None:
    """Rewrite docs/track2-ranking.md from a live triage run."""
    typer.echo(f"wrote {write_ranking(_REPO_ROOT).relative_to(_REPO_ROOT).as_posix()}")


if __name__ == "__main__":
    app()
