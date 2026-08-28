"""Machine verification of the PMIDs this project writes down.

``docs/track2-hypothesis.md`` §0 asserts that "every citation was resolved against
PubMed before being written down." That is an assertion about a past action, and an
assertion is exactly the thing a reader cannot check. This package turns it into a
property that CI re-proves.

The split across the three modules is a PRIV-05 boundary, not an aesthetic one:

* :mod:`tools.citations.extract` reads our own files and never opens a socket.
* :mod:`tools.citations.ledger` compares what we wrote against a committed record
  of what PubMed returned, and never opens a socket.
* :mod:`tools.citations.resolve` is the ONLY module here that touches the network,
  imports ``urllib`` at module scope so the import lint can see it, and is imported
  by nothing on the always-on path.

The always-on offline test therefore runs with zero network. The online test is
opt-in (``-m network`` and/or ``MVA_VERIFY_CITATIONS=1``) and refreshes the ledger.
"""

from __future__ import annotations
