"""Entry point for ``uv run python -m tools.track2``."""

from __future__ import annotations

from pathlib import Path

import typer

from tools.track2.drugs_to_avoid import write

app = typer.Typer(add_completion=False, help="Regenerate the Track 2 generated documents.")

_REPO_ROOT = Path(__file__).resolve().parents[2]


@app.callback()
def main() -> None:
    """Present so that commands stay named subcommands as more are added.

    Typer collapses a single-command app into a bare invocation, which would make
    `python -m tools.track2 drugs-to-avoid` an error today and a working command
    the day a second generator lands. Pinning the shape now keeps the justfile
    recipe and the docs stable.
    """


@app.command("drugs-to-avoid")
def drugs_to_avoid() -> None:
    """Rewrite docs/track2-drugs-to-avoid.md from a live triage run."""
    typer.echo(f"wrote {write(_REPO_ROOT).relative_to(_REPO_ROOT).as_posix()}")


if __name__ == "__main__":
    app()
