# Repo follow-up backlog

Written 2026-08-28, after Track 1 attempt 1 banked 100.0 / F-max 1.000.

Everything here is deferred work, recorded now so it survives the gap between
"we noticed this" and "we had time for it". Submissions close **2026-10-24**,
so none of it is urgent and all of it is affordable. Nothing in this file blocks
a submission; if something did, it would be in `docs/next-actions.md` instead.

Ordered by what a judge would notice first, not by effort.

---

## 0. Why this list is ordered the way it is

Two things set the order:

1. **The automated metric is saturated.** Attempt 1 banked 100.0 rank points and
   F-max 1.000 on its first try — the mathematical ceiling for both. A score that
   cannot move is a format check, not a differentiator. See
   `docs/submission-ledger.md` for the full argument.
2. **Therefore every remaining distinction is qualitative.** What a reviewer
   reads, clones, and can re-run is the whole contest.

That reframes this backlog: these are not engineering niceties, they are the
only axis left.

### What the repo already has

| Property | This repo |
|---|---|
| Snakemake DAG over the whole pipeline | 5 rule files, containers + envs pinned |
| Architecture decision records | 27 |
| Test functions | 1,341 across 52 files (~1,520 cases) |
| Privacy audit with self-proving detectors | 13/13 rules re-proved each run |
| Differential oracle against the real scorer | 2,000 randomised trials |

Making the whole pipeline a reproducible DAG is the most unusual of these and
the hardest to retrofit, which is why it is worth saying out loud.

### What is missing

| Gap | Why it matters |
|---|---|
| **Blind ablation** — which evidence sources were load-bearing | It is the only answer to "did the pipeline do the work, or did ClinVar?" |
| Discard ledger — what was considered and rejected | GP-13 is enforced inside the pipeline and invisible from outside it |
| Claims ledger — every assertion traced to its source | Turns a stylistic virtue into a checkable one |
| Negative fixtures — inputs that *must* fail | The suite proves good input works and infers the rest |
| A "start here" path for a reviewer with 10 minutes | 28 ADRs and 275 files is not a front door |

The ablation is the important one. It is also the thing the real-callset run
currently in flight exists to produce.

### Track 2: where the argument actually stands

**Direction-of-effect reasoning is not by itself a differentiator.** The
inversion argument — that aneuploidy-selective oncology agents lose their
therapeutic window when the patient's own tissue *is* the aneuploid compartment —
is reachable by anyone who works the problem carefully, so leading with "we
noticed the trap" says less than it looks like.

**What does hold up:**

- The nicotinamide-versus-niacin call, and the fact that it **flipped on its own
  evidence, 2026-08-29**. The §5 table used to show niacin scoring **higher**
  (0.410 vs 0.365) and losing on a gate about whether the cancer-susceptibility
  question had been *asked*. Both `worsens_cin` cells turned out to be wrong:
  niacin's malignancy counts by randomised arm are public (HPS2-THRIVE and
  AIM-HIGH) and it has been through a human micronucleus assay, while
  nicotinamide's cell had been filled from a skin-cancer trial that measures no
  chromosome number. Correcting the two cells — no code, no weight, no gate
  touched — inverted the ranking: nicotinic acid accepted at 0.500, nicotinamide
  rejected at 0.275. The molecule is not the differentiator; the fact that we can
  show the rule, the condition that flipped it, and the run that flipped on its
  own, is.
- **§5.2** — 7 of 12 rejections are `oncogenic_risk` where the risk was never
  measured rather than found. "The bottleneck in repurposing for
  chromosomal-instability disorders is not finding targets; it is that the one
  safety question that matters has not been asked."
- The blindness constraint (ADR 0008 / GP-44): the agent that wrote the
  hypothesis could not read patient data, which forced a genotype-conditional
  framing.

**Practices worth having that we do not yet** — each is a backlog item below:
CI-verified citations (§2.5), a standalone drugs-to-avoid artifact (§2.6), a
stated probability of benefit, a withdrawn-claims table (§1.4), and a positive
control on the science pipeline (§2.7).

---

## 1. Reviewer-facing — highest value per hour

### 1.1 A "start here" path
A reviewer arrives with ten minutes and no context. Today the README asks them
to orient themselves across 28 ADRs and 275 files. Needed: one path that runs in
minutes on synthetic data, produces a visible artifact, and links the three docs
worth reading in order. Everything else stays where it is; this is a front door,
not a reorganisation.

### 1.2 Publish the ablation
Blocked on the real-callset run. Once it lands: re-run with each evidence source
withheld in turn, and report which candidates survive. This closes the gap §0
names as the important one, and it is the strongest possible answer to "did the
pipeline do the work, or did ClinVar?"

### 1.3 Discard ledger
GP-13 says a deleted candidate is invisible and unfalsifiable. We enforce that
inside the pipeline but do not surface it. A reviewer cannot currently see what
was considered and dropped, or why. Render the discards as a first-class
artifact, ranked, with the rule that killed each one.

### 1.4 Claims ledger
Every factual assertion in the Track 2 report already carries a PMID. Nothing
mechanically checks that. A generated table — claim, source, resolution status —
turns a stylistic virtue into a verifiable one.

---

## 2. Correctness machinery we designed but did not build

### 2.1 Harvest-the-real-secrets differential scan
**This is the one that matters.** The privacy audit scans for *patterns*. It did
not catch the IRB number that leaked, because nothing told it that specific
string was sensitive. The fix: take the actual secrets from the case volume,
scan the repo for them directly, and fail on any hit. It is a different question
from "does this look like an MRN" and it is the question that was actually wrong.

Constraint: the scanner must never write a found secret into any output. It
reports *that* a match occurred, its file and offset, never its content
(PRIV-09).

### 2.2 Two-sided self-test
`check_detector_liveness` proves all 13 rules still fire (the false-negative
half). Nothing proves they *don't* fire on clean input — the false-positive half.
A detector battery that matches everything passes the liveness check perfectly.

### 2.3 Negative fixtures
Inputs that must be rejected: a truncated VCF, a wrong-build contig set, a
mismatched reference allele, a callset with no PASS variants. Today the suite
proves the pipeline works on good input and infers the rest.

### 2.4 Digest lockfile with a policy block
Resource digests are pinned in several places. One lockfile, one policy
statement about what may change without review.

### 2.5 Citations verified in CI, not asserted in prose
`docs/track2-hypothesis.md` §0 states that "every citation was resolved against
PubMed before being written down" — 47 unique PMIDs. That is an assertion, and an
assertion does not catch a hallucinated citation. A CI check *enforces* it, which
is the whole difference. Existence is not enough: a PMID that resolves to a different
paper than the one cited is the failure mode worth catching, so first author,
year and journal must match too. Offline-by-default against a committed ledger,
with an opt-in online refresh, because PRIV-05 forbids network on the sensitive
path. **In flight on `feat/track2-controls`.**

### 2.5b The citation parser is fragile in a way that fails quietly
Two separate incidents, both found by luck rather than design:

1. **Asterisk parity.** A citation inside an italic- or bold-wrapped span
   (`*"… *Mol Cell* 2005 …"*`) does not parse, because the wrapper's asterisks
   break the `*Journal*` match. The citation looks perfect to a human.
2. **Author position.** `(HPS2-THRIVE, Landray et al., *N Engl J Med* 2014)`
   parses the first author as "THRIVE" — the parser takes the last name before
   the journal, so any parenthetical preceding the surname captures it.

Both degrade a citation to `UNANNOTATED`, which means **existence-only checking**
— the identifier resolves, the content is never compared. `MINIMUM_FORM_COUNTS`
and the `unannotated == 0` assertion are what surfaced these, and they are the
only reason the class of failure is visible at all. Worth either teaching the
parser these two shapes or documenting the house citation form so writers do not
have to rediscover it. **The house form that works:
`Surname et al., *Journal* YYYY, PMID nnnnnnn`, on one line, not nested inside
emphasis markers.**

### 2.6 A standalone drugs-to-avoid artifact
The 12 rejected agents live inline in a §5 table. A separate "do not use" page is
more useful to a clinician and more legible to a judge, and must be *generated*
from the pipeline's rejection output so it cannot drift. The distinction that
matters on it is GP-14: **OPPOSES** (measured wrong direction) versus
`oncogenic_risk` (never measured). **In flight on `feat/track2-controls`.**

### 2.7 A positive control for the science pipeline
We built a positive control for the privacy scanner (`detector_canaries()`, 13/13
re-proved each run) and none for the science. Nothing currently proves the
intervention pipeline can still accept a true positive — a pipeline that rejects
everything scores identically to a careful one on the present tests. This is our
own idiom, unapplied where it matters most. **In flight on
`feat/track2-controls`.**

---

## 3. Performance and ergonomics

### 3.1 Parallel annotation
In flight on branch `perf/parallel-annotate`. The annotate stage pins one core
of ten on a 4.96 M-variant callset. Hard requirement: bit-identical output
(GP-30 extended — parallelism must not change the answer any more than batch
size may).

### 3.2 Progress output on long stages
A multi-hour stage prints nothing. The only way to confirm it was alive was to
sample its stack with `sample(1)`. Being handled alongside 3.1.

### 3.3a Measured: a genome-wide run cannot complete, and never could

A real whole-genome run was started 2026-08-28 16:39 and stopped 7 h 30 m later,
by hand, because it was about to exhaust the case volume. It was not stopped for
impatience. Measured at the point it was killed:

| | |
|---|---|
| Variants to annotate | 4,962,060 |
| Annotated in 7 h 30 m | **~1,911,224 (38.5%)** — estimated at 6,103 bytes/record |
| `annotated.json` | 11.2 GB, projecting **~30 GB** complete |
| Evidence ledger | 24 GB, growing *faster* than the artifact |
| Case volume | 60 GB total, 17 GB free at kill time |

**Annotation alone projects to ~19.5 hours and roughly 95 GB of intermediates on
a 60 GB volume.** The run could not have finished. Left alone it would have
filled the encrypted volume that holds the patient callset and the submitted
files, some hours after midnight, and produced nothing.

Two things follow, and the second is the important one:

1. The reproduction claim in `docs/track1-method.md` — that the submitted ranking
   came from an ad-hoc script rather than the library — **cannot be closed by
   waiting**. It needs the design change in §3.5.
2. Nothing in the repository knew this. There is no size estimate, no disk
   pre-flight, no progress output, and no projection from a first batch. A stage
   that cannot possibly complete should say so in its first minute, not its
   twentieth hour. **That check is worth more than the speedup in §3.1**, which
   would have turned an impossible 19.5-hour run into an impossible 3-hour one.

For contrast, the ad-hoc script completed the same callset genome-wide in
**14 minutes**, because it filtered before annotating. That is the trade §3.5
exists to decide.

### 3.3 The evidence ledger is the real bottleneck at genome scale
Measured on the 4.96 M-variant run: the SQLite evidence ledger reached **11 GB**
and was still growing at ~57 MiB/min, while the process sat at **8.5% CPU**. The
annotate stage stopped being CPU-bound some time ago; it is now blocked writing
evidence to a sparsebundle-backed volume. Parallelising annotation (§3.1) will
not help past this point — the next run to profile is the ledger writer, not the
adapters. Options: batch the inserts, write the ledger outside the encrypted
image (only if it holds no patient-derived content — check before assuming), or
make the per-variant evidence rows optional for a whole-genome pass.

Disk headroom was fine (38 GB free, ~11 h at that rate), but a genome-scale run
needing ~18 GB of scratch is worth stating in the runbook rather than
discovering.

### 3.4 Resume from artifacts
There is no `--start-from`. Any change to a late stage costs a full re-run
including a 31-minute ingest. `normalised.json` and `annotated.json` already
exist on disk; the stages just cannot be entered from them.

### 3.5 The prefilter question — a design decision, not a bug
The ad-hoc script finished genome-wide in 14 minutes by filtering *before*
annotating. The library annotates all 4.96 M variants first, then filters, and
takes hours. That ordering is deliberate: GP-13 makes pre-annotation deletion
unauditable. But the gap is ~50×, far larger than parallelism buys, so the
trade deserves an explicit ADR rather than an implicit inheritance. Options
include annotating everything but *streaming* the discards to the ledger, so
the audit trail survives a cheaper order.

---

## 4. Documentation

- **`docs/next-actions.md` is stale.** It still opens with "the pipeline is not
  yet wired to it". Rewrite or retire.
- **Runbook:** add the `⌘⇧G` note — `nobrowse` APFS mounts are invisible to file
  pickers, which is why submission files had to be staged elsewhere. This cost
  real time and will cost it again.
- **Reconcile the submitted report.** The file that reached the organizers
  (`neonbhai_track1_report.as-submitted.md`, preserved on the case volume) is
  ten lines shorter than the working copy and is missing the confidentiality
  banner. Decide which is canonical.

---

## 5. Deliberately not doing

- **Spending Track 1 attempts 2–6 on the score.** It is at the ceiling. See
  `docs/submission-ledger.md`.
- **Tuning anything against the known answer.** The real-callset run is a
  held-out test with known ground truth. Its value is entirely destroyed by
  adjusting weights toward it.
