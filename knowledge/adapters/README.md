# The annotation adapter boundary

Everything the pipeline knows about a variant beyond its own VCF record arrives
through an **adapter**: a small object satisfying one of three Protocols in
[`src/mva/annotation/base.py`](../../src/mva/annotation/base.py). The
orchestration in `src/mva/annotation/service.py` never names a tool, a file
format or a URL. That is the whole point of the boundary — the pipeline's
scientific logic is written once, against typed models, and the question of
*where the annotation came from* is answered in exactly one swappable place.

| Protocol | Method | Returns | Real-world equivalent |
|---|---|---|---|
| `ConsequenceAdapter` | `annotate(variant_ids)` | `Mapping[str, tuple[ConsequenceAnnotation, ...]]` | VEP, SnpEff, Nirvana, SpliceAI |
| `FrequencyAdapter` | `frequencies(variant_ids)` | `Mapping[str, tuple[PopulationFrequency, ...]]` | gnomAD, TOPMed, local cohort |
| `ClinicalAdapter` | `assertions(variant_ids)` | `Mapping[str, tuple[ClinicalAssertion, ...]]` | ClinVar, ClinGen, an internal curation database |

All three also expose `name` and `version`, which are mandatory: they are
stamped onto every `EvidenceItem` the stage emits and into the run manifest, so
any claim in a report can be traced back to the tool and release that produced
it.

## What ships in this repository

A **synthetic substitute**, graded as such in `docs/maturity-ledger.md`:

| Adapter | `name` | `version` | Backing file |
|---|---|---|---|
| `LocalConsequenceAdapter` | `local-tsv-consequence` | `synthetic-v0.0` | `knowledge/public/consequences.tsv` |
| `LocalFrequencyAdapter` | `local-tsv-frequency` | `synthetic-v0.0` | `knowledge/public/frequencies.tsv` |
| `NullClinicalAdapter` | `null-clinical` | `unavailable-v0.0` | — (nothing on record) |

The genes are fictional, the allele frequencies are invented, and there is no
clinical source at all. This is not hidden behind a config flag: the adapter
names and versions say it, every `EvidenceItem.limitations` says it, and the run
emits a `GP-20` warning naming each mocked adapter. A synthetic substitute is
never described as biologically valid.

`NullClinicalAdapter` exists rather than leaving the slot empty because an
adapter that returns nothing *and says who it is* keeps the gap visible in the
coverage table. Silence and absence-of-a-source look identical downstream
otherwise, and neither is evidence of benignity.

## The three rules an adapter must obey

**1. Omit what you do not know (GP-14).** The return type is a `Mapping`, and a
variant the source has never heard of must be **missing from it** — not present
with an empty tuple, and absolutely not present with `allele_frequency=0.0`.
`knowledge/public/frequencies.tsv` demonstrates the distinction the type is
protecting: `chr15:40200000:C:T` is listed with `allele_frequency=0.0` over
`allele_number=152312`, which is a real observation of zero carriers in a large
cohort. `chr11:5000000:A:GT` is simply absent, which means nothing is known. The
first is strong evidence of rarity; the second is no evidence at all. A default
of zero would silently convert the second into the first, and would do it most
often for sites that reference cohorts cover worst.

**2. Return typed models (GP-02).** Parse at the boundary, once, inside the
adapter. No dicts, no `Any`, no deferred parsing. `PopulationFrequency` cannot
be constructed without `source`, `version` and `population` (GP-18) — take those
from the source data, never from the adapter's own identity.

**3. Preserve every transcript.** `annotate` returns *all* transcript
annotations the source holds for a variant. Collapsing to the canonical or
MANE-Select transcript is a data-loss bug: a variant can be benign on
MANE-Select and splice-disrupting on the tissue-relevant isoform. The service
attaches whatever it is given, in full; ordering is presentation, never
selection.

## Slotting in a real adapter

Nothing in `service.py` changes. The work is: implement the Protocol, declare
your maturity, wire it at the composition root, and document the grade.

```python
# src/mva/annotation/vep_local.py  (illustrative)
class VepConsequenceAdapter:
    """Consequences from a locally-installed VEP + offline cache."""

    def __init__(self, vep_binary: Path, cache_dir: Path, *, version: str) -> None: ...

    @property
    def name(self) -> str:
        return "ensembl-vep"

    @property
    def version(self) -> str:
        return self._version          # e.g. "112.0/GRCh38/cache-112"

    @property
    def synthetic(self) -> bool:
        return False                  # opt out of the mock disclosure, deliberately

    def annotate(
        self, variant_ids: Sequence[str]
    ) -> Mapping[str, tuple[ConsequenceAnnotation, ...]]:
        ...                           # subprocess -> parse -> typed models
```

Then bind it:

```python
adapters = AdapterSet(
    consequence=VepConsequenceAdapter(...),
    frequency=GnomadTabixAdapter(...),
    clinical=ClinvarVcfAdapter(...),
)
result = annotate_variants(records, adapters=adapters, clock=clock)
```

Notes on each field of that swap:

* **`synthetic`** is an optional property, and `is_synthetic()` **fails closed**:
  an adapter that does not declare it is treated as a mock and gets the "not
  biologically valid" disclosure on every evidence item it produces. Declaring
  `synthetic = False` is how a real adapter opts out, and it should be as
  deliberate as it looks (GP-20).
* **`version`** must identify the artifact that produced the answer, not the
  Python package: VEP's cache release, the gnomAD release, the ClinVar weekly
  VCF date. `EvidenceItem` refuses a `DATABASE_ASSERTION` whose citation carries
  no version, because "AF = 0.0001" is unreproducible without knowing whether it
  came from gnomAD v2.1.1 or v4.1.
* **SpliceAI** is a `ConsequenceAdapter` concern, not a fourth Protocol: its
  score lands in `ConsequenceAnnotation.splice_ai_delta_max`, and other
  pathogenicity scores (CADD, REVEL) go in `pathogenicity_scores` under their
  own names. They remain predictions, so the evidence stays
  `COMPUTATIONAL_PREDICTION` / `IN_SILICO_PREDICTION` however impressive the
  number (GP-12).
* **Add a row to `docs/maturity-ledger.md`** grading the adapter `real`,
  `synthetic-substitute` or `stub`, and register the resource in the knowledge
  manifest (below). Reports surface the grade of everything they depend on.

## The knowledge manifest and integrity

`knowledge/manifests/knowledge.yaml` pins every local table by sha256. It is a
generated file:

```bash
uv run python -c "from pathlib import Path; \
  from mva.annotation import render_manifest_yaml; \
  Path('knowledge/manifests/knowledge.yaml').write_text( \
    render_manifest_yaml(Path('knowledge')), encoding='utf-8')"
```

`load_default_adapters(knowledge_root, manifest_path)` parses the manifest,
verifies **every** declared table's hash, and only then opens anything as data.
A mismatch raises `AdapterUnavailableError` naming the file — never echoing its
contents (PRIV-09). Regenerate the manifest after any table change and review
the diff; hand-editing a hash turns an integrity check into decoration.

Each entry carries `name`, `path` (relative to the knowledge root), `version`,
`sha256`, `retrieved`, `description` and `synthetic`. A real resource is added
with `synthetic: false`, its true release string and a real retrieval date.

## Why no adapter here may touch the network (PRIV-05)

`src/mva/annotation/` is structurally forbidden from importing `requests`,
`httpx`, `urllib`, `http`, `aiohttp`, `ftplib` or `smtplib`. This is enforced by
`tests/unit/test_architecture.py::test_no_network_clients_in_sensitive_stages`,
not by convention, and the same rule covers `ingestion`, `phenotype` and
`prioritization`.

The reason is specific rather than general. A variant coordinate is not
metadata: a handful of rare coordinates identifies an individual, and their
parents. Posting a proband's coordinates to a public annotation endpoint —
Ensembl REST, a variant-normalisation API, a "just this once" curl during
debugging — discloses patient genetic data to a third party, permanently and
irreversibly, with no consent, no data-processing agreement and no way to recall
it. A crashing request that logs its URL leaks the same data into wherever logs
go. There is no version of this that is safe under a hackathon deadline, so the
capability is removed rather than governed.

Two further consequences: a live API cannot give byte-identical repeat runs
(GP-30), and it makes every result contingent on someone else's uptime.

**If a remote source is genuinely needed**, it belongs in a separate,
public-only acquisition step that lives **outside** `src/mva/annotation/` — a
tool that downloads the gnomAD sites VCF, the ClinVar release or a VEP cache,
hashes it, records it in the knowledge manifest, and exits. That step:

* runs before, and separately from, any patient data being loaded;
* sends only *public reference identifiers* — a dataset name, a release tag, a
  URL — and **never** a proband coordinate, genotype, sample ID or pedigree
  detail, in a URL, a query parameter, a request body, a header or a log line;
* writes into `knowledge/`, which is public and committed, so what was fetched
  is reviewable;
* is what the `network_profile` config setting governs. A non-synthetic case may
  not run with `network_profile: online` at all (`CaseConfig` rejects it).

Annotation then runs against the downloaded, hash-pinned artifact, locally and
offline, exactly as it does today against these synthetic TSVs — which is why
swapping the tables for real ones changes the adapter and the manifest, and
nothing else.

## Table schemas

The header row of each TSV *is* the schema; comment lines start with `#`. Empty
cells parse to `None`, never to `0` or `""`.

**`knowledge/public/consequences.tsv`** — one row per (variant, transcript):
`variant_id`, `gene_symbol`, `gene_id`, `transcript_id`, `transcript_biotype`,
`is_canonical`, `is_mane_select`, `consequence_terms` (comma-separated, most
severe first), `impact` (`high|moderate|low|modifier`), `hgvs_c`, `hgvs_p`,
`exon`, `protein_position`, `amino_acids`, `splice_ai_delta_max`, `cadd_phred`,
`revel`. Multiple rows for one `variant_id` are multiple transcripts and are all
retained — see `chr15:40200000:C:T`.

**`knowledge/public/frequencies.tsv`** — one row per (variant, source, version,
population): `variant_id`, `source`, `version`, `population`,
`allele_frequency`, `allele_count`, `allele_number`, `homozygote_count`,
`filter_status`. A variant with no row has **no frequency data**.

`variant_id` in both is the canonical, build-qualified join key
`{build}:{contig}:{pos}:{ref}:{alt}`, e.g. `GRCh38:chr15:40200000:C:T`.
