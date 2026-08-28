"""``python -m tools.citations`` -- refresh the ledger, or check without one.

Two commands, and the split between them is the PRIV-05 boundary made operational:
``check`` is offline and is what CI runs; ``refresh`` is the only thing here that
opens a socket, and it is run by a human (or by ``just verify-citations``).
"""

from __future__ import annotations

from pathlib import Path

import typer

from mva.clock import SystemClock
from tools.citations.extract import collect_written_citations, unique_pmids
from tools.citations.ledger import (
    LEDGER_PATH,
    LEDGER_SCHEMA,
    Ledger,
    compare,
    form_counts,
    format_mismatches,
    load_ledger,
    write_ledger,
)

app = typer.Typer(
    add_completion=False,
    help="Verify that every PMID this project writes down resolves to the paper we cite.",
)

_REPO_ROOT = Path(__file__).resolve().parents[2]


def _report(repo_root: Path) -> int:
    citations = collect_written_citations(repo_root)
    ledger = load_ledger(repo_root / LEDGER_PATH)
    mismatches = compare(citations, ledger)
    counts = form_counts(citations)
    typer.echo(
        f"{len(citations)} citation(s), {len(unique_pmids(citations))} unique PMID(s); "
        f"ledger resolved {ledger.resolved_utc}"
    )
    typer.echo("  written detail: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    if mismatches:
        typer.echo(f"\n{len(mismatches)} MISMATCH(ES):\n{format_mismatches(mismatches)}")
        return 1
    typer.echo("  no mismatches")
    return 0


@app.command()
def check() -> None:
    """Compare our prose against the committed ledger. Offline."""
    raise typer.Exit(_report(_REPO_ROOT))


@app.command()
def refresh() -> None:
    """Resolve every cited PMID against PubMed and rewrite the ledger. ONLINE."""
    # Deferred on purpose (PLC0415): importing the resolver at module scope would
    # pull `urllib` into `check`, which is the offline command CI runs.
    from tools.citations.resolve import ESUMMARY_URL, resolve_pmids  # noqa: PLC0415

    citations = collect_written_citations(_REPO_ROOT)
    pmids = unique_pmids(citations)
    typer.echo(f"resolving {len(pmids)} PMID(s) against {ESUMMARY_URL}")
    records = resolve_pmids(pmids)
    write_ledger(
        _REPO_ROOT / LEDGER_PATH,
        Ledger(
            schema=LEDGER_SCHEMA,
            resolved_utc=SystemClock().now().isoformat().replace("+00:00", "Z"),
            source=f"NCBI E-utilities esummary (db=pubmed) via {ESUMMARY_URL}",
            records=records,
        ),
    )
    typer.echo(f"wrote {LEDGER_PATH}")
    raise typer.Exit(_report(_REPO_ROOT))
