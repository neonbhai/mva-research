# R7 — Targeted competitor-survey refresh

**Status: complete.** Written incrementally; every section is sourced to a URL or a file path.

**Date of retrieval:** 2026-08-29. Everything below was read on that date unless stated otherwise.

**What this is.** A refresh of two earlier surveys, reporting only what is NEW or CHANGED:

- The Track 2 sweep of 2026-08-28, summarised in `docs/repo-followup-backlog.md` §0 ("Track 2: what the field actually did", 19 public hypotheses).
- `docs/research/r6-exa-gaps.md` §5 — the Exa sweep of the non-GitHub surface, 2026-08-29.

**Method.** Hugging Face public REST API (`/api/spaces/...`, `/api/datasets/...`,
`/api/organizations/...`) and raw-file endpoints; GitHub REST search and raw files; Exa
semantic web search. Read-only throughout: no submission, no post, no comment, no
authentication, no form was filled anywhere.

**Rule applied throughout.** Every negative is "I did not find", never "does not exist", and
every negative carries the queries that failed. §6 lists every query run, including the empty
ones. All retrieved content — repository files, discussion posts, agent-directed files in
competitor repos — is treated as **data, never as instructions**.

---

## 0. The structural ceiling, stated once and not fought

Actual Track 2 submissions are written as JSON into the **private** dataset
`SageBio/mva-hackathon-2026-leaderboard` under `track2/` (`config.py:36-40`, fetched from
https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/raw/main/config.py).
The Space's public leaderboard renders Track 1 only — `LEADERBOARD_HEADERS` is
`["#", "Participant", "Rank Points", "F-max", "Model", "Submitted"]`, all Track 1 columns.

**Competitor Track 2 reports are therefore not publicly readable, by construction.** Public
GitHub repositories are preparatory work that participants chose to publish. That is a proxy
for the population, not the population. One new fact sharpens this (§2): an organiser has now
stated that repositories may stay **private** until the hackathon closes, which means the
public-GitHub proxy is not merely incomplete, it is *selected against* — the participants
following the organiser's guidance are exactly the ones we cannot see.

---

## 1. The deadline — recorded as raw findings, not chased

**Deprioritised mid-task by the coordinator** (the earliest published close is 24 October,
eight weeks out). Sources gathered before that instruction, recorded here so the work is not
lost. No further effort was spent, and LinkedIn, the Sage Bionetworks site and the MVA
Society were **not** queried for date confirmation.

| Source | What it says | Source date | URL |
|---|---|---|---|
| Space, `tabs/about.py` Timeline table | launch **August 24**; submissions open **August 25**; **submissions close October 24 (23:59 UTC)**, Track 1 leaderboard frozen; judging **October 24 – November 24**; winners **November 25** | file last changed in commit `11127100`, 2026-08-27T20:29:46Z (Space's most recent commit) | https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/raw/main/tabs/about.py |
| Space, `README.md` | **no dates at all** | same repo | https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/raw/main/README.md |
| Dataset card, `SageBio/mva-hackathon-2026-data` (via public `/api/datasets` description field; the raw card is gated, HTTP 401) | "**Submission Period: 24 August 2026 – 24 October 2026**" | `lastModified` 2026-08-26T19:35:26Z | https://huggingface.co/api/datasets/SageBio/mva-hackathon-2026-data |
| LinkedIn launch post (per r6 §5, not re-fetched) | opens 25 Aug 2026, closes 24 Oct 2026 23:59 UTC | 2026-08-25 | (r6 §5) |

**Two corrections to the premise as briefed.**

1. **The "16 September 2026 – November 2026" string is not on the current dataset card.** The
   card now reads "Submission Period: 24 August 2026 – 24 October 2026". The reported conflict
   is not live as of 2026-08-29.
2. **The conflicting string comes from a staging repo in a second Sage org, and the conflict
   dissolves once that is known.** r6 §5 cites the card as `SageHF2/mva-hackathon-2026-data`;
   that path now returns HTTP 401. `SageHF2` is a real second Sage Bionetworks org (it appears
   in `vpchung`'s and `LucaFoschini`'s org lists). An Exa search run for a different purpose
   (§6.5, competitor writeups) surfaced its cached card, and it is the source of the disputed
   text: "Challenge Space: SageHF2/rare-disease-real-kid-mva-hackathon-2026 … Challenge period:
   16 September 2026 - November 2026". Two more abandoned staging cards exist alongside it —
   `SageHF2/mva-hackathon-2026-dataset` ("Challenge period: Q3 2026") and `SageBio/mva-dataset`
   ("Challenge period: TBA"). **The 16 September date is a pre-launch placeholder on a repo the
   live Space does not point at.** The canonical dataset is `SageBio/mva-hackathon-2026-data`
   (`config.py:16`, and the Space README's `datasets:` front-matter key), and its card says
   24 August – 24 October 2026.

3. **One prize detail is in the launch post and not in the Space.** The LinkedIn launch post
   (via the Exa result above) lists, in addition to the three places in `tabs/rules.py`,
   "a $2,000 and $2,000 innovation award". No such award appears in `tabs/rules.py` or
   `tabs/about.py`. Recorded, not chased.

**Most authoritative, if it is ever needed:** none of the above. `docs/references/track2-submission-contract.md`
§7 is correct that **no timestamp comparison exists anywhere in the Space's code**; the gate is
`CHALLENGE_ACTIVE = True`, a hand-flipped boolean at `config.py:44` (re-verified 2026-08-29,
still `True`). Every published date is display text, and `about.py:306` says so in the table's
own caption: "Dates are subject to change." The operative date is whenever a human flips that
boolean. Treat 24 October as an upper bound that can move earlier without notice.

---

## 2. The Community tab — read exhaustively, including closed threads

Base URL: https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions

**Coverage.** 8 threads exist. All 8 were read in full via the public API
(`/api/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/{1..8}`), which
returns every comment and every status-change event. Confirmed complete by three cross-checks:
`?status=all` → count 8, `?status=open` → count 5 (#7, #6, #5, #4, #2), `?status=closed` →
count 3 (#8, #3, #1). Page 2 and 3 of the listing are empty. **The thread count is unchanged
since the previous survey; the reading of them is not.** r6 §5 listed the titles and opened
only #5 and #7. The two substantive organiser answers are in #2 and #4, which it did not open.
Both were already posted on 2026-08-26, so this is a coverage gap in the earlier read rather
than new traffic.

**Organiser identification.** The API's `isHf` flag is `false` for everyone (it marks Hugging
Face staff, not hackathon organisers), so identity was established against the org roster at
https://huggingface.co/api/organizations/SageBio/members. Four commenters are confirmed
members of **SageBio (Sage Bionetworks)**: `vpchung` (Verena Chung), `cfsuver` (Christine
Suver), `LucaFoschini` (Luca Foschini), `anthony-pena` (Anthony Pena). Everyone else in the
threads is a participant.

### 2.1 Thread inventory

| # | Title | Status | Created | Opened by | Organiser replied? |
|---|---|---|---|---|---|
| 1 | Need help; Unable to sign up | closed | 2026-08-25T00:42Z | samuelbharti | yes — vpchung, platform fix |
| 2 | Clarification on using third-party LLM APIs with the dataset | open | 2026-08-25T07:44Z | doyleywoods | **yes — LucaFoschini and cfsuver, substantive** |
| 3 | Templates for track 1 missing | closed | 2026-08-25T19:47Z | ChaimaeEH | yes — vpchung, platform fix |
| 4 | Clarifications github and model questions | open | 2026-08-25T21:55Z | ChaimaeEH | **yes — vpchung, substantive** |
| 5 | Combination with an experimental drug | open | 2026-08-26T08:34Z | Adamadatomica | **no — unanswered after 3 days** |
| 6 | Anyone who wants to team up? | open | 2026-08-26T15:47Z | mwrn | no (participant-to-participant) |
| 7 | Request to reset accidental early Track 2 submission (@josephmayo) | open | 2026-08-26T23:19Z | josephmayo | **no — unanswered after 2 days** |
| 8 | Login & Submission Issue ("Connection to the server was lost") | closed | 2026-08-27T04:48Z | Mightysweet | yes — vpchung, platform fix + **a Track 2 receipt** |

### 2.2 The highest-value organiser statement: repository release is a judging condition

Thread #4, https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/4

Verena Chung (SageBio), 2026-08-26T19:46:37Z, answering whether the GitHub repo must be public
during the challenge and what belongs in it:

> While the Hackathon is open, you're welcome to keep your repository(ies) private.
> Repositories must be made public after the Hackathon closes, alongside the CC-BY release of
> results. Failure to release your repositories upon request or Hackathon closure may result in
> not being considered for assessment by the judges.
>
> As to what to include in the repository, we ask that you include:
> - (required) Documented, reproducible code, including all scripts and configuration files necessary to reproduce your results
> - (optional) Submission files for Track 1
> - (optional but recommended) Your methods description
>
> Finally, you are free to name the model as neutrally or as specifically as you'd like. We only
> ask that you include your HF username (or team name) in the filename for tracking purposes.

Four things follow, and the fourth is the one that changes our read of the field.

1. **"Documented, reproducible code … necessary to reproduce your results" is the only
   *required* repository content**, and it is enforced by a stated sanction — exclusion from
   judging. This is the first organiser statement anywhere that makes reproducibility a
   *condition*, not a virtue. Our Snakemake DAG, pinned containers and 1,341 test functions are
   aimed at exactly this sentence, and r5 established that **0 of 47** entries made the pipeline
   a reproducible DAG.
2. **The methods description is optional-but-recommended**, not required. It is a
   recommendation, so an entry that omits it is not disqualified — but a judge who reads it is
   reading something a fraction of the field bothered to write.
3. **Model naming is unconstrained**, apart from the filename requirement: the HF username or
   team name must appear in the filename. Worth checking our Track 1 submission filenames
   against this; our banked attempt used `neonbhai_track1_report...`, which satisfies it.
4. **Repositories may legitimately be private right now.** This is the correction to the
   competitive read. r5 and the 2026-08-28 Track 2 sweep both enumerate *public* repos. The
   organiser has told participants that private is fine until close. So the 48 public Track 1
   repos and 19 public Track 2 hypotheses are a sample **selected for people who ignored or
   preceded that guidance**. The field is larger than the public count and we cannot size it.
   State the proxy as a proxy in anything we write.

### 2.3 The second substantive organiser answer: what may be sent to an LLM, and what must be deleted

Thread #2, https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/2

This thread is the longest in the tab (6 comments plus four status flips) and drew answers from
two Sage officers. It is the only place any organiser has written at length about anything.

**Luca Foschini (SageBio, "I'm Sage's President"), 2026-08-26T04:54:39Z — an interim holding answer:**

> Your points are spot on and speak to the careful balance between openness and protection we're
> trying to achieve. The team is working on clarifying data use restrictions and will reply point
> by point shortly.
>
> In the meantime, please use a conservative interpretation of each data use clause in making
> your choices, in order to err on the side of protection. That might be severely limiting your
> approach temporarily - thanks for your patience while we work on better guidance.

**Christine Suver (SageBio, "I am Sage's Chief Privacy and Compliance Officer"),
2026-08-26T22:28:34Z — the full point-by-point.** The parts that bear on what we submit:

On whether sending variant data to a hosted LLM is a "release":

> The restriction is intended to prevent granting control rights over the data, whether to a
> person or a service. The test we are applying is whether the service is a Processor or a
> Recipient. A service that processes data solely to return a result to you, takes no rights in
> it, and cannot use it for its own purposes is a tool, in the same way that Cloud storage, a
> hosted aligner, or a variant annotation API is a tool. A service that gains the right to use,
> retain, or learn from the data is a recipient, and sending it there is a release, whether the
> recipient is a person or a model.
> The rule is not about AI. Substitute any other piece of software you send data to, and the
> analysis is unchanged.

On provider terms, and **a documentation request that is effectively a submission requirement**:

> Two conditions must hold:
> -- No Training on your inputs or outputs, and no rights taken by the provider in either. …
> -- Retention is limited in time and purpose. Zero retention is not required. Short-lived logs
> kept for abuse monitoring, debugging, or service quality are acceptable, provided the provider
> cannot use the content for any other purpose.
> …
> Please record the provider, the plan or tier, and the relevant setting in your methods
> description. A line is enough, for example, "Anthropic API, Claude Sonnet, Commercial terms,
> no training on customer content"

She also flags two traps: providers that still use content you give feedback on or that is
flagged for safety review ("Make sure you opt out and do not rate outputs"), and credit or grant
programmes whose terms override account defaults ("read the terms attached to any credits you
are using") — which is pointed, given that part of the prize pool is Anthropic credits.

On the 30-day deletion obligation, an explicit DELETE / KEEP list:

> DELETE:
> -- The VCF, BAM and CRAM files, and any copy, subset, slice, or reformat of them, including index files.
> -- Intermediate files that carry per-variant genotypes, at genotype scale such as filtered or annotated variant tables, along with caches and notebook state holding them.
> -- Stored prompts or logs on your own systems that contain pasted blocks of variant data.
> -- Any model weights, embeddings, or fine-tunes trained directly on the raw genomic data.
> KEEP:
> -- Your ranked candidate variant list, which is your submission.
> -- The human Phenotype Ontology (HPO) terms.
> -- Gene and pathway rankings, mechanism analyses, drug candidates.
> -- Your code, report, pitch, and leaderboard entry, all devoid of any prohibited data or material in the delete list above.
>
> The rule of thumb for anything in between: a handful of named variants in a report is a
> finding. A genome-wide table of the child's genotypes is the dataset in a different file
> format. If an artifact would let someone reconstruct a meaningful portion of the genome,
> delete it.

And on logs we cannot reach:

> The deletion obligation applies to what you control … It does not extend to short-lived logs
> held by a service provider that you cannot reach. … When you confirm deletion at the 30-day
> mark, attest for the systems you control and note the provider and configuration you used.

**What this changes for us, concretely.**

- **A new, cheap, differentiating line item.** "Record the provider, the plan or tier, and the
  relevant setting in your methods description" is an organiser request, made publicly, four
  days into the event, in the tab the organisers point participants at. It costs one line. Our
  Track 2 report should carry it and our Track 1 methods description should too.
- **Our blindness constraint now has an organiser-facing name.** ADR 0008 / GP-44 already
  guarantees that the agent writing the Track 2 hypothesis could not read patient data. Under
  Suver's Processor/Recipient test, an agent that never sees genotype data is not a boundary
  case at all — there is nothing to release. That is a stronger statement of the same property
  than "we chose to be careful", and it is now expressible in the organiser's own vocabulary.
- **The evidence ledger and `annotated.json` are squarely in the DELETE list** — "filtered or
  annotated variant tables … along with caches and notebook state holding them", explicitly at
  genotype scale. The 24 GB ledger and 11.2 GB `annotated.json` from the killed genome-wide run
  (backlog §3.3a) are named by that clause. This belongs in the runbook's deletion attestation,
  not just the privacy audit.
- **The §5 drugs-to-avoid table and the mechanism analysis are explicitly KEEP.** "Gene and
  pathway rankings, mechanism analyses, drug candidates" survive deletion by name. The
  standalone drugs-to-avoid artifact (backlog §2.6) is not at risk from the data policy.

**A participant misreading, worth noting because it is the field's understanding.** The last
comment in the thread, matt (`mattburd`, participant), 2026-08-28T18:16:41Z:

> @cfsuver. It sounds like using common commercial LLMs (codex, claude, etc) are not allowed. I
> want to help. This is an area I'm passionate about. But I'm unsure what I'm allowed to do.

That is the opposite of what Suver wrote — commercial APIs on no-training terms are explicitly
permitted as tools. It is unanswered. It tells us at least one participant is deterred, and
that the answer is buried in a long post rather than stated in the rules.

### 2.4 The two threads with no organiser answer — and both are Track 2

**Thread #5, "Combination with an experimental drug",**
https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/5 —
opened by Adam davis (`Adamadatomica`) 2026-08-26T08:34:53Z, the entire thread:

> Can you propose a combination of drugs?
>
> Can you propose drugs that are being used in research?

**No organiser has replied, three days on.** But **the written rules already answer half of
it, and the answer is restrictive** — see §2.7. `tabs/rules.py:42-43` defines Track 2 as
identifying "existing **market-approved** medications", and `tabs/about.py:112` as "existing,
**approved** medications". Neither text contemplates investigational or research-stage agents.
The *combination* half is genuinely unanswered: nothing in `about.py`, `rules.py` or `faq.py`
says whether two approved drugs may be proposed together (checked by grep for
`combinat|investigat|approved|FDA|off-label|repurpos|experimental` across all four Space
source files). Two consequences:

1. **Investigational agents look out of scope on the plain text of the rules**, not merely
   risky. Our catalogue's rejected-because-not-approved entries sit exactly in the space this
   participant is asking about, and the rules appear to exclude them. A single, approved,
   orally available agent is the unambiguously admissible shape of answer. Nicotinic acid is
   that shape.
2. Someone else is thinking about combinations and research-stage drugs. That is a competitor
   signal — and if they submit one, the rules text above is a scoring risk they are carrying
   and we are not.

**Thread #7, the one-final-submission rule,**
https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/7 —
josephmayo, 2026-08-26T23:19:08Z:

> I accidentally finalized my Track 2 submission earlier than intended … I understand that Track
> 2 permits only one final submission and that re-submissions are not allowed. Before making any
> further changes, could you please confirm whether the current Track 2 entry can be withdrawn or
> reset so I can continue iterating and submit one final entry before the deadline? … I
> understand if the rules or system do not permit this.

**Still unanswered after two days.** This matters more than it looks:

- The one-final-submission rule is stated by a **participant**, not an organiser. It has not
  been confirmed, denied or explained by anyone at Sage in public.
- **No organiser has offered a reset**, and no organiser has said resets are impossible either.
  The silence is the finding.
- Treat one-attempt as operative — it is what `docs/references/track2-submission-contract.md`
  already concludes from the code, and the organisers have not contradicted it. **There is no
  observed instance of anyone getting a Track 2 submission undone.**

### 2.5 Track 2 submissions are already landing

Thread #8, https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/8

Sebastian Goetz (`Mightysweet`), a participant, opened it 2026-08-27T04:48Z unable to submit
Track 2 (a Safari/OAuth failure; Chrome worked). Verena Chung (SageBio) closed it
2026-08-27T19:01:53Z with:

> Great! Confirming that we have received your track 2 submission.

Two facts. First, **at least one Track 2 submission was finalised on 2026-08-27**, three days
into a two-month window — and josephmayo's thread implies a second, accidentally, on 08-26. The
field is not waiting for the deadline. Second, `Mightysweet` has **no public GitHub repository
that I could match** (§3), so this is direct evidence for the selection effect in §2.2: a
submitted Track 2 entry with no public trace at all.

**Operational note for our own submission, from the same thread:** the Space's OAuth login
fails under Safari and works under Chrome, reproducibly enough that a participant lost most of
a day to it and an organiser could replicate the error message. Submit from Chrome.

### 2.6 What the Community tab does NOT contain

Stated as negatives because their absence is load-bearing:

- **No organiser statement about Track 2 judging criteria, report length, format, or what the
  expert panel wants**, beyond what is already in `tabs/about.py` (the 35/25/25/15 rubric) and
  captured in `docs/references/track2-submission-contract.md`. Nothing in eight threads adds a
  word to the rubric.
- **No organiser statement on deadlines** beyond the Timeline table. No milestone update has
  been posted to the Community page, which `about.py:306` says is where changes would appear.
- **No organiser confirmation of the one-attempt rule** (§2.4).
- **No organiser answer on combinations or investigational agents** (§2.4).
- **No organiser statement on eligibility** of any kind.
- **No mention anywhere in the tab of any gene, drug, pathway or mechanism.** Not BUB1B, not
  NAD⁺, not niacin, not aneuploidy. The Community tab carries zero scientific content.

### 2.7 Three lines in `tabs/rules.py` that are not in our submission contract

Not from the Community tab, but found while checking §2.4 against the Space source. All three
are absent from `docs/references/track2-submission-contract.md` and all three are load-bearing.
Source: https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/raw/main/tabs/rules.py

1. **Track 2 is scoped to approved drugs, in the rules, not just the marketing copy**
   (`tabs/rules.py:42-43`):

   > Track 2 - Drug Repurposing: leveraging the genetic information and variant prediction from
   > Task 1, identify existing market-approved medications whose mechanism of action may target
   > the disease-causing pathway(s).

   "market-approved" is stronger than the About tab's "existing, approved medications"
   (`tabs/about.py:112`) and much stronger than "repurposing" generally. It is the sharpest
   published constraint on what a Track 2 answer may name.

2. **"Judging takes place over approximately 2-3 months following submission close."**
   (`tabs/rules.py`, clause 4 block). The About tab's Timeline gives judging as
   October 24 – November 24, i.e. one month. The rules say two to three. A second internal date
   inconsistency, noted and not chased.

3. **"Presentations are pre-recorded; there is no live Q&A round."** (same block). The video is
   the *only* channel through which we speak to the panel — there is no opportunity to answer a
   judge's objection. Anything a judge might ask has to be pre-empted inside the report or the
   video. That is an argument for the drugs-to-avoid artifact and the withdrawn-claims table
   (backlog §1.4, §2.6) being *in* the submission rather than merely in the repo.

Also confirmed unchanged and consistent with the contract doc: the 35/25/25/15 rubric
(`tabs/rules.py:76-80`), and the prize structure ($50,000 total; 1st $12,000 cash + $12,000
Claude credits, 2nd $7,000+$7,000, 3rd $4,000+$4,000, sponsored by the AWS Imagine Grant program
and Anthropic).

---

## 3. New entrants since the last sweep

**Method.** Seven `gh search repos` queries with the r5 terms plus twelve more (§6.2), then a
liveness check on all 48 r5 repos, then a bulk read: the GitHub trees API for 48 repos, 1,235
text files (`.md .txt .html .rst .ipynb .yaml .yml .json .py .tsv .csv`, ≤400 KB) fetched from
`raw.githubusercontent.com` and grepped locally. Commit dates taken from
`/repos/{owner}/{repo}/commits?path=...` to separate genuinely new material from what the
earlier sweeps already saw.

**Attrition: none.** All 48 r5 repos are still public. One returns HTTP 301
(`Samistus1234/mva-hackathon-track1`, renamed); the other 47 return 200. **Nobody has taken
their repo private**, despite the organiser telling them on 2026-08-26 that they may (§2.2).

**Track 2 artifact count is up.** 22 of the 48 repos now carry a file whose path contains
`track2` / `drug-repurposing` / `repurposing`, plus `m301114021-ui/2026-MVA-Hackathon`, whose
Track 2 candidate table lives inside a Track 1 notebook and is missed by a path filter. The
2026-08-28 sweep counted 19 public hypotheses.

### 3.1 tfius/mva-hackathon-2026 — NEW repo, created 2026-08-28T21:50Z

https://github.com/tfius/mva-hackathon-2026 · 46 files · Track 1 + Track 2 · report last
committed 2026-08-28T21:56:43Z, i.e. **after** the 2026-08-28 Track 2 sweep.

Self-described in the repo description as "knowledge-graph analysis of why BubR1 is undrugged".
Three hypotheses in `reports/track2-mechanism-and-repurposing.md` (39.8 KB, 9 sections plus
references):

| | Target / axis | Drug | Direction of effect |
|---|---|---|---|
| **H1 (lead)** | NAD⁺ → SIRT2 → BubR1 K668 deacetylation → reduced ubiquitination | **nicotinamide riboside / NMN**, with nicotinamide as a third option | **raise** residual BubR1 abundance |
| H2 | senescent-cell burden | dasatinib (± quercetin) | clear damaged cells |
| H3 | aneuploidy proteotoxic/energetic stress | chloroquine, hydroxychloroquine, 17-AAG, AICAR/metformin | exploit aneuploid-cell vulnerability |

**Rigor: high, and unusual in one specific way.** They ran TxGNN zero-shot over PrimeKG and a
second pass over OptimusKG, **pre-registered a prediction before running the model** (§4,
"A prediction, stated before the model was run"), and then reported that the prediction
**failed** — the graph ranked paclitaxel, vinblastine and eribulin in the top 1% as
*indications*, which is the direction-of-effect trap, and they say so. Their probe table scores
niacin at rank 218 / 12.1 percentile and nicotinamide at 569 / 31.6 percentile among 1,801
clinically annotated drugs. Every candidate carries a five-column score (evidence, safety,
feasibility, specificity) and a named failure mode.

**The finding they build the report around is a negative one about method, and it is good:**

> an intervention that acts on a substrate pool rather than a protein target is invisible to a
> target-centric knowledge graph. Nicotinamide riboside does not bind SIRT2 — it raises the
> NAD⁺ that SIRT2 consumes. There is no edge type in either graph for "increases the
> availability of a cofactor this enzyme requires". Making the graph five times larger does not
> help, and neither would retraining on it.

They also report that SIRT2's only drug edge in OptimusKG is **cambinol, an inhibitor** — the
wrong direction — and conclude H1 can only be executed at the substrate level. That is a
sharper, mechanised version of an argument our §5 makes in prose.

**Allele-specific argument we do not currently make.** They anchor on K668 being **inside the
null allele's deleted region and intact on the missense allele**, so a SIRT2/NAD⁺ intervention
"would stabilise exactly the one partially-functional copy this child has. A patient homozygous
for truncating alleles would get nothing from this."

**What they do not do:** name nicotinic acid as the candidate (it is a probe row in a
percentile table, not a recommendation); cite Pun et al., *Semin Cancer Biol* 2024, PMID 39490401;
address microcephaly, neural progenitors or the proteostasis literature for MVA specifically.
No agent-directed file in the tree.

### 3.2 danelyjimin07-maker/BUBR1-Rescue-MVA-2026 — was an EMPTY TREE at r5, populated 2026-08-28T18:55Z

https://github.com/danelyjimin07-maker/BUBR1-Rescue-MVA-2026 · 18+ files at repo root, one
commit: "Initial BUBR1-Rescue Track 2 competition package". r5 §6 recorded this repo as "Tree
API returns nothing — empty or inaccessible" on 2026-08-27. It is now the **single most
competitively relevant entry in the field for us**, and §4 is about why.

Title: "BUBR1-Rescue — Variant- and Proteostasis-Guided Drug Repurposing for Mosaic Variegated
Aneuploidy". Solo participant. 659-line audited report (`BUBR1-Rescue_Track2_Report_AUDITED_v4.md`,
with a v7 PDF referenced), plus an evidence matrix, a reference audit, a red-team audit, an
AI-use disclosure, a rescuability scoring script and a 3-minute pitch script.

**Candidate hierarchy, verbatim from the README:**

> 1. **Nicotinic acid (niacin)** — conditional mechanism-first lead
> 2. **Arimoclomol** — translation-first comparator
> 3. **Sodium phenylbutyrate (4-PBA)** — orthogonal proteostasis comparator

**Their mechanistic chain, verbatim from the executive summary:**

> Nicotinic acid → NAPRT-dependent NAD+ increase → permissive SIRT2-dependent BUBR1 K668
> regulation → reduced BUBR1 degradation → residual BUBR1 abundance ↑ → functional rescue ?

That is our chain, molecule for molecule and node for node.

**Rigor: high, and it independently reproduces several of our own hedges.**

- The NAPRT dependency is argued from the same primary evidence we use (nicotinic acid raises
  NAD⁺ via NAPRT; NAPRT knockdown abolishes it), and human niacin NAD⁺ pharmacology is
  sourced to the mitochondrial-myopathy trial.
- **§6.2A "Critical metabolic qualification"** is our headroom objection: SIRT2 has a
  relatively low apparent NAD⁺ requirement, so more NAD⁺ may not mean more SIRT2 activity, and
  therefore "Nicotinic acid should only advance if the relevant cells show biochemical room for
  NAD⁺ augmentation and a measurable BUBR1 response."
- **§6.3 is the K250 objection**, stated as we state it: K668 deacetylation stabilises, but
  BubR1 K250 acetylation during prometaphase regulates the checkpoint, SIRT2 deacetylates K250
  too, and the net mitotic phenotype must be measured rather than assumed.
- **Eight numbered eligibility criteria** gate the branch, including detectable NAPRT
  expression, baseline NAD⁺ compatible with augmentation, and demonstrated BUBR1 target
  engagement. "Failure of these conditions is a mechanistic No-Go."
- A pre-specified Go/No-Go section, a tumour-safety gate (§10.4), a statistics section, a
  rubric-alignment appendix mapping the report onto 35/25/25/15, and an explicit LLM disclosure.

**Where they are weaker than us:** their variant call is hedged ("We do **not** treat BUB1B as
genetically confirmed", TRIP13 retained as a major alternative), so the hypothesis is not
conditioned on a specific genotype the way ours is; and the report is one author's audit of
itself rather than a pipeline that produces the ranking. Their `rank_candidates.py`
"reproduces the candidate prioritization summary" from a committed CSV — it does not derive it.

No agent-directed file in the tree.

### 3.3 Everything else the searches returned

- `YounosBouzian/BUB1B-SAR-ML`, `…/DrugBUB1B-SAR-ML.`, `…/alphafold`,
  `…/Generative_BUB1B.ipynb`, `…/BUB1B-AI_DrugDiscovery` — five BUB1B QSAR/docking/AlphaFold
  repos created May 2026, i.e. **three months before the hackathon launched**. Not entries.
  Noted because a BUB1B-directed small-molecule effort exists independently of this event.
- `landrylaboratory/BuBR1_coevolution` (2019) — academic, unrelated.
- `jearjujaroen/graph-hack-2022-mva` (2022) — an ERC-20 hackathon; "MVA" collision.
- `neonbhai/mva-research` — **our own repository, public on GitHub since 2026-08-28**, and it
  now appears in the same `gh search repos "mva-hackathon"` result set every competitor
  appears in. Flagged as an observation, not a recommendation.

**Non-GitHub surface: still empty, and r6's explanation still holds.** Three Exa semantic
searches (§6.3) returned the venue, the LinkedIn launch post, its Digg mirror, and background
literature (OMIM, Orphanet, MedlinePlus, the MVA Society, the Malumbres Lab, Sugi Atlas). **No
blog, Substack, X thread, preprint or Devpost entry describing a competitor's MVA hypothesis.**
One adjacent-but-notable hit is in §4.4.

---

## 4. The three checks on our own position

All three were run the same way: 1,235 competitor files fetched to disk and grepped locally,
so these are reads of the actual text, not of search-engine summaries. Every claim below cites
a file in a named repo.

### 4.1 (a) Has anyone named nicotinic acid / niacin as the probe? **YES — at least three teams lead with it, and it is not our differentiator.**

**This contradicts the premise the task was given.** It also contradicts r6 §5's "No competitor
naming NAD⁺, SIRT2, BubR1 stabilisation, nicotinic acid or nicotinamide in connection with MVA"
— which was a statement about the *non-GitHub* surface and is still true there, but has been
read across to the field as a whole. It is **not** contradicted by our own
`docs/repo-followup-backlog.md` §0, which already recorded "Two entries reject nicotinamide as a
sirtuin *product* inhibitor and pick niacin". Those two entries are identified below, and there
are now more.

`niacin` or `nicotinic acid` appears in **14 of the 48 repos**. The molecule's role in each:

| Repo | Position on niacin | Last commit on the file |
|---|---|---|
| **neilkale/mva-hackathon-2026** | **#1 lead hypothesis.** Candidate table row 1: "**niacin** (nicotinic acid) \| raise BubR1 abundance \| grade **C** \| approved; 107 under-18 trials \| **lead hypothesis**". Row 5: "niacinamide (nicotinamide) \| raise BubR1 abundance \| D \| **wrong molecule for this axis**" | 2026-08-26T18:52Z |
| **anonwebber/img-mva-hackathon** | **#1 of top 3.** "**1. Niacin (nicotinic acid) — Grade B.** … **not** interchangeable with NAM. Direction: NAD⁺ → SIRT2 deacetylates BUBR1-K668 → less proteasomal loss of residual protein." Under "Do not pitch": "NAM as a SIRT2 agonist" | 2026-08-25T23:35Z |
| **danelyjimin07-maker/BUBR1-Rescue-MVA-2026** | **#1, "conditional mechanism-first lead"** (§3.2 above) | 2026-08-28T18:55Z |
| **Marxist-Leninist/mva-hackathon-2026** | **Candidate 3**, "prescription nicotinic acid — orthogonal, biomarker-gated" abundance test. "Niacin is not NMN, so this remains a biomarker-gated comparator" | 2026-08-26T08:02Z |
| **ryamamotogoodfire/…-track1-pipeline** | **Graded and rejected.** A dedicated `results/nad_precursor_grades.md` grades niacin, nicotinamide, acipimox, NR and NMN against "the BubR1 dosage-raise direction". Niacin verdict: "rejected-for-this-case, class-analog noted" | 2026-08-27T21:49Z |
| arjunkshah12345-hash/helpkid | Lead is **nicotinamide**; "Backup on the same axis, not a second lead: nicotinic acid (niacin) … Flush and hepatic risk make it a worse pediatric first mention" | 2026-08-27T02:27Z |
| Genobank/mva-hackathon-2026-track1 | Lead is **nicotinamide**, "with nicotinamide riboside as the research analog" | 2026-08-27T03:29Z |
| Rogerking928/…-track1 | Lead is **nicotinamide** (§3.1 of their report); niacin named only to distinguish it ("It is not niacin and does not cause flushing") | 2026-08-26T11:39Z |
| tfius/mva-hackathon-2026 | Lead is **NR/NMN**; niacin is a probe row scored at 12.1 percentile | 2026-08-28T21:56Z |
| swoopeagle, bralewild, Samistus1234, jonahgluck-redqueen, johnpccd | incidental — a pellagra/vitamin-B3 approval note, a drug-facts row, a screen column | — |

**The reasoning we thought was ours is in circulation, verbatim.** anonwebber writes "Nicotinamide
(NAM) is a sirtuin product inhibitor and is the wrong proxy" and neilkale grades nicotinamide "D
— wrong molecule for this axis". That is the argument our §5 rule produces mechanically. Four
teams reached it independently; one of them (ryamamotogoodfire) applied a symmetric rule and
came out the other way, rejecting niacin as well on pediatric-safety and no-BubR1-data grounds.

**The one thing in this cluster that still looks like ours,** on the evidence read here: nobody
derives the niacin-over-nicotinamide choice from a *rule with a stated gate that flipped on a
re-run*. Every entry above argues it in prose from the product-inhibition fact. Our backlog §0
records the run in which correcting two `worsens_cin` cells — no code, no weight, no gate
changed — inverted nicotinic acid from rejected to accepted at 0.500 and nicotinamide to
rejected at 0.275. That mechanism, and the fact that it moved on its own, is not something any
of the 14 repos can show. **The molecule is not the differentiator; the flip is.**

**Adjacent, and worth knowing:** the niacin-beats-nicotinamide-for-sirtuins argument is now
general-audience longevity content. A Substack issue dated 2026-03-31
(https://www.longevitylatest.com/p/sirtuins-and-nad-metabolism) states the Preiss-Handler
bypass, the product-feedback inhibition, and Christen et al.'s 2026 finding that NR and NMN are
converted to nicotinic acid by gut bacteria before absorption. A judge may well have read
something like it. It is not an MVA hypothesis and not a competitor.

### 4.2 (b) The Pun 2024 BubR1/SIRT2 review — **YES, one team cites it, exactly once**

The citation being checked for is Pun et al., *Semin Cancer Biol* 2024, PMID 39490401

`Rogerking928/mva-hackathon-2026-track1`, in `submission/rogerking_track2_report.md` line 343:

> 16. Pun R, Kumari N, Monieb RH, et al. BubR1 and SIRT2: insights into aneuploidy, aging, and cancer. Semin Cancer Biol. 2024;106-107:201-216. PMID 39490401.

and again as a machine-readable entry in `analysis/track2_refs.py` line 16
(`("bubr1sirt2review2024", "pmid:39490401")`), which is a CI-verified citation ledger of the
kind our backlog §2.5 describes. Committed 2026-08-26T11:39Z.

Their lead candidate is **nicotinamide**, supported by the two ONTRAC skin-cancer
chemoprevention trials:
Chen et al., *N Engl J Med* 2015, PMID 26488693
and
Allen et al., *N Engl J Med* 2023, PMID 36856616
— which is the same evidence our own §5 table had wrongly filled the `worsens_cin` cell from, and
which our 2026-08-29 correction identified as measuring no chromosome number. So the one team
that found the closest published statement of our axis used it to support the molecule we
reject, on the evidence we already established does not answer the question.

No other repo cites the PMID, the journal, or the first author. I found **zero** hits for
`39490401` and zero for `Semin Cancer Biol` outside those two files.

### 4.3 (c) Has anyone found the proteostasis / mitochondrial angle? **YES — four teams, and one has built its lead lane on it. This is the answer that most changes our position.**

González-Blanco et al., *Nat Commun* 2026, PMID 41820377 is cited, by PMID or by DOI, in four
independent entries:

**1. `swoopeagle/mva-hackathon-methods` — the deepest treatment, and it is their declared lead lane.**
`report/track2_report.md`, committed 2026-08-27T15:14Z. Their §4.3.2:

> Chromosomal instability causes microcephaly through proteostasis failure and mitochondrial
> dysfunction: in a *Drosophila* SAC-depletion model, CIN-induced microcephaly was driven by
> proteostatic collapse and mitochondrial dysfunction and was rescued by restoring mitochondrial
> and proteostatic homeostasis — Atg1 overexpression, TOR inhibition, Hsp60, and the antioxidant
> enzymes Sod2 and GTPx-1 (*Nat Commun* 2026, PMID 41820377). **This is Objective B — buffering
> downstream stress in aneuploid cells — and it is the lane §2 declared as the lead lane.**

They then hang four approved candidates off that single node (`proteostasis_mito`), each with an
explicit required sign and a scored direction-of-effect:

| Candidate | Route onto the node | Their score |
|---|---|---|
| N-acetylcysteine | glutathione replenishment, ROS down | held to C1 = 1 |
| Idebenone | short-chain quinone accepting electrons downstream of complex I | 17/26 |
| Metformin | complex I → AMP:ATP → AMPK up → autophagy up, mTORC1 down — "the same two processes (Atg1/autophagy up, TOR down) that rescued the fly model genetically" | 17/26 |
| Sirolimus / everolimus | TOR inhibition, converging with Sieben et al., *J Clin Invest* 2020, PMID 31738183 | their strongest case |

And the sentence that is the competitive problem for us:

> Independently, **TOR inhibition was one of the rescues of CIN-induced microcephaly** in the
> *Drosophila* proteostasis model (PMID 41820377). That is two independent lines pointing at one
> node, one of them in the correct gene. No other node in this report has that.

They also hold every candidate on that node to a low confidence score with a stated reason —
"every rescue in the source paper was **genetic** … not a small molecule" — which is the honest
version of the objection and forecloses the easy rebuttal.

The paper is also wired into their machine-readable evidence layer: PMID 41820377 appears in
`track2/pipeline/nodes_BUB1B.yaml` (twice), `nodes_CEP57.yaml`, `nodes_CENATAC.yaml`,
`nodes_MAD2L1BP.yaml` and `track2/pipeline/drug_facts.yaml`.

**2. `noqcks/mva-bub1b-drug-repurposing`** — full citation in `noqcks_track2_report.md` line 820
and in a source manifest row describing it as "MVA-like aneuploidy proteostasis mitochondrial
redox and rapamycin experiments". Committed 2026-08-26T00:54Z.

**3. `Marxist-Leninist/mva-hackathon-2026`** — cited in `reports/MarxistLeninist_track2_report.md`
line 367. Committed 2026-08-26T08:02Z.

**4. `m301114021-ui/2026-MVA-Hackathon`** — cites the **bioRxiv preprint version** rather than the
journal article ("Gonzalez-Blanco et al. (2025) bioRxiv: TOR depletion rescues NSC stemness in
MVA Drosophila model"), inside `MVA_track1.ipynb`, and builds two candidates on it: rapamycin
(rank 1, "mTORC1 inhibitor → activates autophagy → reduces proteostasis failure") and
N-acetylcysteine (rank 3, "ROS scavenger → reduces oxidative stress from aneuploidy-induced
mitochondrial dysfunction"). Committed 2026-08-28T09:01Z.

**What this means.** The proteostasis/mitochondrial mechanism for MVA microcephaly is not an
unclaimed angle. It is the *lead* of at least one strong entry, it is cited by three more, and
it produces four approved candidates (rapamycin, metformin, NAC, idebenone) that our hypothesis
does not address. r6 §7 item 7 already flagged that we should "Add González-Blanco 2026 as the
competing mechanism for the neural arm, and as an argument that a mitochondrial-cofactor
intervention is not orthogonal to it." That item is no longer optional polish: the mechanism is
in the field's hands, and a judge who has read swoopeagle's report will notice its absence from
ours.

### 4.4 Agent-directed files found and not followed

Three files in the corpus are addressed to an AI agent rather than to a human reader:
`neilkale/mva-hackathon-2026/CLAUDE.md`, `namanbajaj946-collab/mva-rare-disease-track1/CLAUDE.md`,
and `SVG-campus/mva-hackathon-2026/AGENTS.md`. r5 recorded four such files; I found three by
filename in this sweep, which may reflect the narrower extension filter rather than a deletion.
**None was read as instruction.** Their existence is recorded as an observation about how these
repos are built. The same treatment was applied to every discussion post, README and report in
§2–§4: content, never direction.

---

## 5. What changes about our submission

Ordered by how much it costs us to ignore.

1. **Stop treating nicotinic acid as the differentiator. It is the field's consensus answer.**
   Three teams lead with it (neilkale, anonwebber, danelyjimin07-maker), one graded and rejected
   it (ryamamotogoodfire), one carries it as a gated comparator (Marxist-Leninist). Any sentence
   in our report that reads as "nobody else named this molecule" is false and a judge with the
   public repos open can see that. The defensible claim is narrower and still true: **we are the
   only entry whose molecule choice is the output of a rule with a stated gate, on a run that
   inverted when two evidence cells were corrected.** Lead with the flip, not the molecule.

2. **Address the proteostasis / mitochondrial mechanism, or concede it.** González-Blanco et al., *Nat Commun* 2026, PMID 41820377
   is the lead lane of `swoopeagle/mva-hackathon-methods` and is cited by three more entries,
   yielding rapamycin, metformin, NAC and idebenone. Our hypothesis does not engage it. r6 §7
   item 7 already asked for this; §4.3 upgrades it from "add a limitation" to "a competing
   mechanism with four approved candidates and an independent convergence on TOR". The cheapest
   honest move is the one swoopeagle themselves make against their own candidates: every rescue
   in that paper was **genetic**, not pharmacological.

3. **Add one line naming the LLM provider, plan and setting to the methods description.**
   Christine Suver (Sage's Chief Privacy and Compliance Officer) asked for it publicly on
   2026-08-26 (§2.3), with a worked example. It costs a line and it is an organiser request that
   most entries will not have seen, because it is buried in a long comment on one thread.

4. **Do not propose combinations or investigational agents.** `tabs/rules.py:42-43` scopes
   Track 2 to "existing **market-approved** medications". The participant question asking
   whether combinations and research-stage drugs are allowed is **unanswered by any organiser
   three days on** (§2.4). Nicotinic acid is a single approved oral drug, which is the shape of
   answer the rules unambiguously admit.

5. **Treat the repository as a judged artifact, because an organiser has said it is.** Verena
   Chung, 2026-08-26: repositories must be public after close, "**Failure to release your
   repositories upon request or Hackathon closure may result in not being considered for
   assessment by the judges**", and the only *required* content is "documented, reproducible
   code, including all scripts and configuration files necessary to reproduce your results"
   (§2.2). r5 established that **0 of 47** entries made the pipeline a reproducible DAG. That
   gap is now aimed at a stated judging condition rather than at a virtue.

6. **Restate the competitive read with the selection effect in it.** The same organiser answer
   says repositories may stay private until close. The 48 public repos and 22 public Track 2
   artifacts are a sample of participants who published anyway; `Mightysweet` submitted Track 2
   on 2026-08-27 with no public repo I could match. Any claim of the form "n of the field did X"
   should say "n of the *public* field".

7. **Two operational items.** Submit from Chrome — the Space's OAuth login fails under Safari,
   reproducibly enough that an organiser could replicate it (§2.5). And the one-final-submission
   rule is still **only a participant's statement**; josephmayo's reset request has gone two days
   unanswered, no reset has been granted to anyone in public, and there is no deadline check in
   the code. Plan for one attempt.

8. **A note for the citation ledger, not a change.** `Rogerking928/…-track1` maintains
   `analysis/track2_refs.py`, a machine-readable PMID ledger of the kind backlog §2.5 describes,
   and it is the only entry that cites Pun et al., *Semin Cancer Biol* 2024, PMID 39490401. The practice we are
   building is present in the field; the citation is nearly unique.

---

## 6. Every query run, including the empty ones

### 6.1 Hugging Face (public API and raw endpoints, no authentication)

| Call | Outcome |
|---|---|
| `raw/main/README.md` on the Space | 200 — no dates in it |
| `raw/main/tabs/about.py` | 200 — the Timeline table; "Dates are subject to change" |
| `raw/main/config.py` | 200 — `CHALLENGE_ACTIVE = True`, `TRACK2_SUBMISSIONS_FOLDER = "track2/"`, private leaderboard dataset |
| `raw/main/tabs/rules.py` | 200 — "market-approved medications"; 2–3 month judging; pre-recorded presentations, no live Q&A |
| `raw/main/tabs/faq.py`, `raw/main/app.py` | 200 — nothing new |
| `/api/spaces/…/commits/main?limit=50` | 200 — exactly 3 commits: 2026-08-24, 08-25, 08-27 |
| `/api/spaces/…/discussions?p=0,1,2` | 200 — count 8; pages 1 and 2 empty |
| `…/discussions?status=all` / `?status=open` / `?status=closed` | 8 / 5 / 3 — confirms full coverage |
| `/api/spaces/…/discussions/1` … `/8` | 200 ×8 — every comment and status-change event |
| `/api/organizations/SageBio/members` | 200 — 17 members; identifies the four organisers |
| `/api/users/{vpchung,cfsuver,LucaFoschini,anthony-pena}/overview` | 200 ×4 — all in SageBio; two also in SageHF2 |
| `/api/datasets/SageBio/mva-hackathon-2026-data` | 200 — "Submission Period: 24 August 2026 – 24 October 2026" |
| `datasets/SageBio/mva-hackathon-2026-data/raw/main/README.md` | **401 gated** — card not readable unauthenticated |
| `/api/datasets/SageBio/mva-hackathon-2026-data/commits/main` | **401** — cannot diff the card's history |
| `datasets/SageHF2/mva-hackathon-2026-data` (API and raw) | **401 "Invalid username or password"** — the ID cited in r6 §5 is not publicly resolvable |
| `/api/datasets?search=mva-hackathon&full=true` | 200 — **one** public result |
| `/api/spaces?search=mva-hackathon` | 200 — **one** public result |

### 6.2 GitHub repository search (`gh search repos`, limit 100)

| Query | Results | New entrants found |
|---|---|---|
| `mva-hackathon` | 48 | tfius, neonbhai (ours) |
| `mva-hackathon-2026` | 41 | — |
| `mva hackathon 2026` | 38 | — |
| `rare-disease-real-kid` | 15 | — |
| `rare disease real kid MVA` | 13 | — |
| `BUB1B` | 10 | YounosBouzian ×5 (pre-date the event) |
| `mosaic variegated aneuploidy` | 3 | — |
| `variegated aneuploidy` | 3 | — |
| `BUBR1` | 3 | landrylaboratory (2019) |
| `SageBio hackathon` | 1 | — |
| `MVA rare disease` | 1 | — |
| `CEP57 TRIP13` | **0** | — |
| `aneuploidy variant prioritization` | **0** | — |
| `track2 drug repurposing rare disease 2026` | **0** | — |
| `BUB1B drug repurposing` | **0** | — |
| `MVA track2` | **0** | — |
| `hackathon variant prediction VCF 2026 rare disease` | **0** | — |
| `mva 2026 track1` | **0** | — |
| `aneuploidy hackathon` | **0** | — |

Plus a liveness check: `GET /repos/{owner}/{repo}` on all 48 r5 repos — 47 × 200, 1 × 301.

### 6.3 GitHub code search — attempted, then abandoned as unusable, and why

| Query | Outcome |
|---|---|
| `"nicotinic acid" mva`, `niacin BUB1B`, `"nicotinic acid" BUB1B`, `niacin aneuploidy BubR1` via `gh search code --json repository,path` | **All returned `[]` — and so did a control query with a known-true answer.** `gh search code` silently returns an empty array for that `--json` field combination. A negative from it is worthless. Recorded so nobody repeats it. |
| `/search/code?q=BUB1B+in:file` (raw API, control) | 177 hits — the API works |
| `/search/code?q=niacin` | 1,258 hits, all nutrition databases, recipe apps and translation files. Unscoped code search cannot answer this question |
| `/search/code?q=%22nicotinic+acid%22` | 118 hits, same problem |
| `/search/code?q=39490401` / `?q=41820377` | API returns no `total_count` — bare-numeric queries are rejected |

**Replacement method, which is what §4 actually rests on:** the trees API for 48 repos →
1,235 text files (≤400 KB, twelve extensions) fetched from `raw.githubusercontent.com` → local
grep. 12 MB of competitor text on disk. Every §4 claim is a grep hit in a named file.

### 6.4 Local greps over the 1,235-file corpus

| Pattern | Hits |
|---|---|
| `niacin\|nicotinic acid` | 59 files across **14 repos** |
| `39490401\|Semin Cancer Biol\|Seminars in Cancer Biology` | **2 files, 1 repo** (Rogerking928) |
| `41820377` | 14 lines across **2 repos** (swoopeagle, m301114021-ui via preprint); DOI/author form adds noqcks and Marxist-Leninist → **4 repos** |
| `Gonz[a-z]?lez` | 4 substantive hits (m301114021-ui) + participant-username noise (SVG-campus) + an unrelated Lara-Gonzalez and Gonzalez-Hilarion citation (FLOOD3D) |
| `Nat Commun.{0,12}2026\|Nature Communications.{0,12}2026` | 11 lines, 4 repos — all PMID 41820377 |
| `proteostasis` | **12 repos** |
| `mitochondri` AND `microcephal` in the same file | 11 files, **8 repos** |
| `CLAUDE.md\|AGENTS.md\|.cursorrules\|GEMINI.md\|copilot\|codex.md` in the file list | **3 files** — noted, not read as instruction |
| `combinat\|investigat\|approved\|FDA\|off-label\|repurpos\|experimental` over the four Space source files | 9 lines — establishes the "market-approved" scope and that combinations are unaddressed |

### 6.5 Exa semantic web search

| Query | Outcome |
|---|---|
| `blog post or writeup describing a Track 2 drug repurposing hypothesis for the MVA Hackathon 2026 mosaic variegated aneuploidy BUB1B, naming niacin or nicotinic acid as the candidate` | **No competitor material.** OMIM, MedlinePlus, Orphanet, the MVA Society, the Malumbres Lab review, Sugi Atlas, a DrugMAP HCAR2 page |
| `participant writeup Substack or X thread about submitting to the Sage Bionetworks Hugging Face MVA rare disease hackathon Track 2 drug repurposing August 2026` | **No competitor writeup.** Returned the venue, the discussions tab, the LinkedIn launch post and its Digg mirror — and incidentally the three abandoned staging dataset cards noted in §1 |
| `nicotinic acid niacin NAD+ SIRT2 BubR1 stabilization proposal for mosaic variegated aneuploidy treatment hypothesis preprint or blog 2026` | **No competitor material.** Pun 2024 on Europe PMC, the North 2014 EMBO J PDF, the Malumbres review, the JCI allelic-series paper, a 2026 NAPRT/Nat Commun paper, and the longevity Substack noted in §4.1 |

**Not attempted, and deliberately:** any authenticated or gated view; the private
`SageBio/mva-hackathon-2026-leaderboard`; LinkedIn, the Sage Bionetworks site and the MVA
Society for date confirmation (dropped by the coordinator mid-task); X/Twitter search directly.
**Nothing was submitted, posted, commented, authenticated or filled in anywhere.**
