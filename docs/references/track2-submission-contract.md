# Track 2 submission contract (verified)

Vendored into the repo deliberately, same reasoning as
`track1-submission-contract.md`: the Space is mutable and Track 2 gives us
**one attempt with no resubmission**, so a rushed re-derivation at submit time
is exactly how the attempt gets wasted. Treat this file as the contract of
record and re-verify before the one shot.

**Source:** raw files read from
`https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/raw/main/`
— `tabs/submit_track2.py`, `app.py`, `config.py`, `utils.py`, `tabs/rules.py`,
`tabs/faq.py`, `tabs/about.py`, `tabs/submit_track1.py`, `README.md`,
`requirements.txt`, and `static/templates/methods_description_form.xlsx`
(unzipped and read). Every file's byte size was checked against
`https://huggingface.co/api/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/tree/main?recursive=true`.
Verified **2026-08-28** against Space commit
`1112710080520b3a2848d11b6ce1327dafbe79cf` (`lastModified 2026-08-27T20:29:46Z`,
runtime stage `RUNNING`, hardware `cpu-upgrade`).

All line numbers below are for `tabs/submit_track2.py` unless another file is
named. Read-only verification: nothing was submitted, POSTed, or authenticated.

**Correction to a standing assumption:** the Space is **not gated**. The HF API
reports `"private": false, "gated": false` for the Space; only
`SageBio/mva-hackathon-2026-data` is gated. You need to *sign in* to the Space
(HF OAuth, `README.md:10` `hf_oauth: true`), not to accept terms on it.

---

## 1. The tab exists, and it is called "Submit - Track 2"

```python
183:    with gr.Tab("Submit - Track 2") as tab:
```

It is rendered only when `CHALLENGE_ACTIVE` is true:

```python
# app.py
106:        if CHALLENGE_ACTIVE:
107:            quota_t1 = submit_track1.render()
108:            status_t2 = submit_track2.render()
```

`config.py:44` currently has `CHALLENGE_ACTIVE = True`. There is no separate
Track 2 endpoint — it is a Gradio tab in the same app, submitting through the
same session. The button is labelled **"Submit"** (line 223), not "Submit &
Score" like Track 1 (`tabs/submit_track1.py:270`), and there is no score
returned.

---

## 2. Exact fields, in declared order

From `render()`, lines 197–222. Fields 1–4 are stacked in the left column
(line 198); field 5 is alone in the right column (line 218).

| # | Widget | Label as shown | Required | Line |
|---|---|---|---|---|
| 1 | `gr.Textbox` | `Team / Display Name (optional)` | no | 199–204 |
| 2 | `gr.Textbox` | `GitHub repo URL *` | **yes** | 205–208 |
| 3 | `gr.Textbox` | `Pitch video URL *` | **yes** | 209–212 |
| 4 | `gr.Textbox` (3 lines) | `Notes for judges (optional)` | no | 213–217 |
| 5 | `gr.File` | `Report file (PDF or Markdown) *` | **yes** | 219–222 |

The `*` in the labels is the author's own required-marker; the enforcement is
in `_handle_submit`, quoted in §3.

The callback wiring passes them in a *different* order from the visual layout —
harmless, but it is why the function signature reads report-before-notes:

```python
227:            inputs=[team_name_box, github_box, video_box, report_file, notes_box],
```

Field 1 default, when left blank (line 131):

```python
131:    display_name = (display_name_input or "").strip() or hf_username_to_display_slug(hf_username)
```

`hf_username_to_display_slug` (`utils.py:53-56`) lowercases and hyphenates. If
you are submitting as a team, type the team name — the blank default is your
own username, slugified.

**There is no other field.** No abstract box, no proband field, no track-1
linkage, no declaration checkbox, no email field, no affiliation. `gr.Checkbox`
does not appear anywhere in the Space source.

---

## 3. Validation, in the order it actually runs

`_handle_submit`, lines 127–152. Each check returns a warning string and stops;
nothing is written.

```python
127:    hf_username = get_hf_username(request, oauth_profile)
128:    if not hf_username:
129:        return "⚠️ Please sign in with your Hugging Face account to submit."
...
133:    count = _track2_submissions_by_user(hf_username)
134:    if count >= 1:
135:        return "⚠️ You have already submitted a Track 2 entry. ..."
136:
137:    github_url = (github_url or "").strip()
138:    if not github_url.startswith("https://github.com/"):
139:        return "⚠️ GitHub URL must start with `https://github.com/`."
140:    video_url = (video_url or "").strip()
141:    if not video_url:
142:        return "⚠️ A pitch video URL is required. ..."
143:    if report_file is None:
144:        return "⚠️ Please upload your report (PDF or Markdown)."
145:
146:    report_path = report_file if isinstance(report_file, str) else report_file.name
147:    if Path(report_path).suffix.lower() not in {".pdf", ".md"}:
148:        return "⚠️ Report must be a PDF (.pdf) or Markdown (.md) file."
```

### GitHub URL
A single **case-sensitive literal prefix test**, line 138. That is the whole
validation. Consequences, all read directly off `str.startswith`:

- `http://github.com/...` — **rejected**.
- `https://www.github.com/...` — **rejected**.
- `https://GitHub.com/...` — **rejected** (`startswith` is case-sensitive).
- `git@github.com:org/repo.git` — **rejected**.
- `https://github.com/` with nothing after it — **accepted**.
- A URL to a **private** repo, a deleted repo, a 404, or a gist — **accepted**.
  Nothing fetches the URL. "public" appears only in prose (line 43).

### Video URL
```python
141:    if not video_url:
```
**Non-empty after `.strip()` is the entire test.** No scheme check, no host
allow-list, no regex, no length check, no reachability check. "YouTube or
Vimeo" appears only inside the *error message* (line 142) and the instructions
(lines 44–45) — neither is enforced. The string `banana` would be accepted and
stored. A typo, an unlisted/private link, a Drive link, or a link that expires
all pass silently and there is no second attempt to fix them.

### Report file
A **file upload** (`gr.File`, line 219), not a text paste. Accepted extensions
are exactly `{".pdf", ".md"}` compared after `.lower()` (line 147), so `.PDF`
and `.MD` are fine and `.markdown`, `.txt`, `.docx`, `.html`, `.zip` are not.
`file_types=[".pdf", ".md"]` (line 221) is the picker filter; line 147 is the
real gate.

**No size limit is set anywhere in the Space source.** `app.launch(...)`
(`app.py:117-123`) passes no `max_file_size`, and `_handle_submit` never
inspects the size. Whatever limit exists is HF/Gradio infrastructure, which is
not readable from here — see §9.

---

## 4. What is stored, and where

```python
154:    entry = {
155:        "hf_username": hf_username,
156:        "display_name": display_name,
157:        "github_url": github_url,
158:        "video_url": (video_url or "").strip(),
159:        "report_filename": report_name,
160:        "notes": (notes or "").strip(),
161:        "submitted_at": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
162:    }
```

Written to the private dataset `SageBio/mva-hackathon-2026-leaderboard`
(`config.py:31`) under `track2/` (`config.py:38`), in a per-user subfolder
(lines 55–57), as two separate commits:

```python
62:        upload_file(                      # 1. the metadata JSON
64:            path_in_repo=f"{subfolder}{filename}",     # sub_<UTC-timestamp>.json
70:        upload_file(                      # 2. the report itself
72:            path_in_repo=f"{subfolder}{entry['report_filename']}",
```

### The report filename is stored verbatim — this differs from Track 1

Track 1 renames the report before upload:

```python
# utils.py
77:        ext = Path(entry["report_filename"]).suffix
78:        entry["report_filename"] = f"report_{timestamp}{ext}"
```

Track 2 does **not**. Line 72 uses `entry['report_filename']` directly, and
that value is the raw uploaded basename (lines 150–152). `_sanitize` is applied
to the *username* (line 55) and never to the filename. So:

- the filename you upload is the filename the judging panel sees;
- it is also a path component in a git repo, unsanitised.

Which is why the tab asks you to embed your identity in it:

```
36: When saving your report, please include your username (or team name) in the filename. ...
39: * `jane-doe_track2_report.pdf` (`.md` also accepted)
```

Use a plain `[A-Za-z0-9._-]` filename. No spaces, no `#`, no unicode, no
leading `-`.

---

## 5. The one-attempt limit IS enforced in code

Guard at lines 133–135, counter at lines 85–105:

```python
93:            return sum(
94:                1 for f in api.list_repo_files(
95:                    LEADERBOARD_DATASET, repo_type="dataset", token=token
96:                )
97:                if f.startswith(f"{TRACK2_SUBMISSIONS_FOLDER}{safe_user}/sub_") and f.endswith(".json")
98:            )
99:        except Exception:
100:            return 0
```

Four things follow, and they matter:

1. **A second submit returns the string at line 135 and writes nothing.** It is
   a refusal, not a silent overwrite and not an error page.
2. **The limit is per HF username**, keyed on `_sanitize(hf_username)`
   (`utils.py:32-34`). The code cannot see teams. "One per team" is a *rule*,
   not a mechanism — lines 22–23 and `tabs/faq.py:23-25` say duplicate
   submissions from other team members "will not be reviewed" / "will be
   ignored". Do not treat a teammate's account as a backup attempt.
3. **A rejected submit does not consume the attempt.** `_save_track2_submission`
   is only reached at line 163, after every check in §3 has passed. A bad
   GitHub prefix or a `.docx` report costs you nothing.
4. **The guard fails open** (lines 99–100). If the HF listing call raises, the
   count reports `0`. It is a best-effort read, not a transactional lock.

### The status card does not update after you submit

```python
225:        submit_btn.click(
226:            fn=_handle_submit,
227:            inputs=[...],
228:            outputs=result_md,          # <- result_md ONLY
229:        )
230:        tab.select(fn=_submission_status, inputs=[], outputs=status_md)
```

`_submission_status` is wired to tab-select (line 230) and app-load
(`app.py:110`), never to the submit button. Track 1 refreshes its quota card on
submit (`tabs/submit_track1.py:275` outputs `[result_md, quota_md]`); Track 2
does not. **After a successful Track 2 submit the card still reads
`0 / 1 submissions used`.** That is a stale widget, not a failed submission.
The only success signal is the confirmation block:

```python
168: ## Track 2 submission received ✓
...
176: - GitHub repository: [{github_url}]({github_url})
177: - Video: {video_line}
178: - Report file: `{report_name}`
```

That block echoes back exactly what was stored. It is the receipt — read it and
screenshot it.

---

## 6. No `proband_id`, no Track 1 linkage

The Track 2 `entry` dict (lines 154–162) has no `proband_id`, no reference to a
Track 1 submission, and no score fields. Compare Track 1, which does carry one:

```python
# tabs/submit_track1.py
178:        "proband_id": proband_id,
```

`tabs/rules.py:42-43` describes Track 2 as "leveraging the genetic information
and variant prediction from Task 1", but nothing joins the two submissions
mechanically. The linkage is entirely something your *report* has to make.

---

## 7. No deadline check exists in the code

There is **no timestamp comparison anywhere in the Space**. `datetime.now(timezone.utc)`
appears exactly twice in Track 2 — line 54 (filename stamp) and line 161
(`submitted_at`) — and both are `strftime` calls, never compared against
anything. Grepping the whole source for a cutoff constant finds nothing.

The only gate is the hand-flipped boolean:

```python
# config.py
44: CHALLENGE_ACTIVE = True
```

When it is set to `False`, `app.py:106-110` simply does not render either submit
tab, and `app.py:97-101` swaps the login CTA for "_Registration is now closed._".

The close date is **display-only HTML** in the About tab:

```html
<!-- tabs/about.py -->
326:          <td>October 24 (23:59 UTC)</td>
327:          <td>Submissions close; Track 1 leaderboard frozen</td>
```

…printed under an explicit caveat two lines above the table:

```html
306:      Dates are subject to change. Any updates to key milestones will be reflected in this table
307:      posted to the <a href="{DISCUSSIONS_URL}" ...>Community page</a>.
```

That table does not state a year (2026 is inferred from context) and the tab
closes when a human edits `config.py`, which could be before or after 23:59 on
the 24th. **Do not plan to submit near the deadline.** There is no code that
will politely accept a submission at 23:58.

---

## 8. The undocumented requirement we did not know about

The tab offers one template, and it is not the Track 1 CSV:

```python
189:                gr.DownloadButton(
190:                    label="📄 Methods description template (Excel)",
191:                    value="static/templates/methods_description_form.xlsx",
```

That workbook has two sheets — `Track 1 methods` and **`Track 2 methods`**. The
Track 2 sheet (cells A1–A16, read by unzipping the xlsx) says, verbatim in A2:

> Please submit your answers with your Track 2 report. All questions are
> voluntary; however, submissions without a methods description may be judged
> less favorably on the **Scientific Rigor and Innovation** criteria, since the
> panel can't assess what it can't see. Track 2 accepts one final submission per
> team.

Scientific Rigor and Innovation are **60% of the Track 2 rubric**
(`tabs/rules.py:76-80`). "Voluntary" with a stated scoring penalty is not
optional in any sense that matters.

The Track 2 questions, in sheet order (A7–A16):

1. Team name
2. Describe your approach in detail (how you went from variant/mechanism to
   candidate medication(s)).
3. Was candidate identification automated (drug-target database query,
   LLM-assisted search) or primarily manual literature review and expert
   curation?
4. Describe any manual literature review or expert curation in detail.
5. Did you rely solely on publicly available drug/target databases and
   literature (**preferred**), or were proprietary sources also used?
6. Describe the public data sources used (e.g. DrugBank, ChEMBL,
   ClinicalTrials.gov, OMIM, PubMed).
7. Describe any proprietary data sources used.
8. How did you characterize the variant's mechanism (LoF vs GoF, pathway
   disrupted, downstream biological consequence)?
9. Estimate of time or effort spent on this analysis.
10. **Method abstract, up to 500 words** — include strengths *and limitations*
    as well as methodological detail.

That 500-word abstract is the **only word limit anywhere in the Space source**.

There is no upload slot for the workbook. The instruction is "submit your
answers *with* your Track 2 report" — so these answers have to live **inside the
report file itself**, since the report is the only artifact that gets uploaded.

### The report's stated content contract

```
29: Prepare a written report (PDF or Markdown) proposing repositioned drug candidates supported by your
30: analysis. Include a characterization of the variant's mechanism (loss-of-function / gain-of-function,
31: pathway disrupted, downstream biological consequence) as the basis for your repurposing rationale.
```

Not enforced by code, but it is the tab's own statement of what the report must
contain, and it maps one-to-one onto the Scientific Rigor criterion described in
`tabs/about.py:150-153` ("Is the variant mechanism characterization sound, and is
the proposed drug repurposing candidate well-supported by that mechanism?").

The video is stated as "a 3-minute pitch video walking through your reasoning
and proposed candidates" uploaded "to YouTube or Vimeo" (lines 44–45). Nothing
checks the duration or the host.

---

## 9. Where `docs/submission-runbook.md` §5 is wrong or incomplete

§5 is currently lines 244–263. Nothing in it is factually *wrong* about the
three artifacts — the gap is that it stops well short of the mechanics.

**Confirmed correct:**
- Report is PDF or Markdown (line 147). ✓
- Public GitHub URL required (line 138 — prefix only; "public" is unchecked). ✓
- 3-minute video URL required (line 141 — presence only). ✓
- Exactly one attempt, enforced (lines 133–135). ✓
- Rubric 35/25/25/15 (`tabs/rules.py:76-80`, `tabs/about.py:141-153`). ✓
- "Expect an analogous 'Submit - Track 2' tab on the same Space, same HF
  sign-in" — correct; the tab is literally `"Submit - Track 2"` (line 183),
  same OAuth. ✓

**Missing — add:**
- **Two optional fields exist that §5 never mentions:** *Team / Display Name*
  and *Notes for judges*. The display name is what appears in the public
  announcement (line 202–203); leaving it blank publishes your slugified HF
  username instead of the team name.
- **The methods description template and its 500-word abstract** (§8). Omitting
  it is stated to cost you on 60% of the rubric. This is the single largest
  omission in §5.
- **The report must characterize the variant mechanism** (lines 29–31) — an
  explicit content requirement, not a stylistic suggestion.
- **The report filename is stored verbatim** and should embed the team/HF
  username (lines 36–39, 72). Track 1 renames; Track 2 does not.
- **The video URL is not validated at all.** §5 lists "video URL plays" as a
  human check — it is now the *only* check that will ever happen. Same for
  "GitHub URL is actually public".
- **The GitHub prefix test is exact and case-sensitive.** `http://`,
  `www.github.com`, `GitHub.com` and SSH URLs are all rejected.
- **A rejected submit does NOT burn the attempt** (§5 line 251–255 implies one
  click is all you get). Validation failures are free; only a fully valid
  payload consumes the attempt.
- **The `0 / 1` status card does not refresh after submitting** (lines 225–230).
  §5 gives no success criterion; without this note an operator may read the
  stale card as a failure.
- **The write is not atomic** — the attempt is consumed by the metadata upload
  before the report upload runs (§10, failure mode 1).
- **The one-per-team limit is enforced per HF username, not per team**, and the
  code fails open on an HF API error (lines 99–100).
- **There is no deadline check in code** (§7). §6 of the runbook presents
  "2026-10-24 23:59 UTC" as a hard limit; it is display text in `about.py`
  under a "dates are subject to change" caveat, and the real gate is a boolean
  someone flips by hand. Treat the date as an upper bound that may move earlier,
  not a countdown to run down.

**Outside §5, one correction:** the runbook's header (lines 15–26) says the
Space and the dataset are "both gated" and instructs accepting terms on both.
The HF API reports `"gated": false` for the Space. Only the dataset is gated.
The Space needs a sign-in, not a terms click.

---

## 10. What would have wasted the attempt

Ranked by how expensive and how silent each one is.

**1. The two-phase write, with no error handling.** The metadata JSON is
committed first (lines 62–69), the report second (lines 70–77), and
`_save_track2_submission` is called at line 163 with no `try`/`except` around
it anywhere in `_handle_submit`. If the *report* upload fails — oversize,
transient HF error, a filename HF rejects as a path — the `sub_*.json` is
already committed, so `_track2_submissions_by_user` now returns 1 and you are
locked out, with a submission that has no report attached and no confirmation
block rendered. Mitigation: a modest-sized PDF and a plain ASCII filename.

**2. A video URL that does not work.** Any non-empty string is accepted
(line 141). An unlisted-but-actually-private YouTube link, a typo, or a link
that expires is indistinguishable from a good one at submit time, and there is
no edit path.

**3. A GitHub repo that is not public yet.** Line 138 checks a prefix, nothing
more. If the repo is private, or the branch isn't pushed, or the URL 404s, the
submission still succeeds and the judges see nothing.

**4. Omitting the methods description.** Costs rubric points on 60% of the
weighting (§8) without producing any error at all.

**5. Leaving Team / Display Name blank when submitting for a team.** Your
slugified HF username is published instead (line 131).

**6. Treating a teammate's account as a second attempt.** Explicitly not
reviewed (lines 22–23, `tabs/faq.py:23-25`).

**7. Trusting the stale `0 / 1` card and hunting for a bug.** Harmless in
itself — the resubmit is refused at line 135 — but if the count check happens to
fail open (lines 99–100) a second entry gets written and it becomes ambiguous
which one the panel reviews.

---

## 11. Pre-flight checklist — run immediately before the one shot

Do the artifact checks **before** touching the tab. Everything below is either
enforced in code at the line cited or unrecoverable after the click.

Artifacts:

- [ ] Report is **`.pdf` or `.md`** (line 147). Not `.docx`, `.txt`, `.html`,
      `.markdown`.
- [ ] Report **filename is plain ASCII `[A-Za-z0-9._-]`, no spaces**, and
      contains the team or HF username, e.g. `<team>_track2_report.pdf`
      (lines 36–39). It is stored verbatim as a repo path (line 72).
- [ ] Report **opens and renders cleanly** from the exact file you are about to
      upload — not from a source doc you exported from.
- [ ] Report **characterizes the variant mechanism** (LoF/GoF, pathway
      disrupted, downstream consequence) as the basis for the repurposing
      rationale (lines 29–31).
- [ ] Report **contains the methods description answers** — all ten Track 2
      questions in §8, including a **≤500-word method abstract with strengths
      and limitations**. There is no separate upload slot; if it isn't in the
      report, the panel never sees it.
- [ ] Report file size is sane (a few MB). No limit is set in code, and an
      upload failure here burns the attempt (§10.1).
- [ ] **GitHub repo is public** — open the URL in a private/logged-out window
      and confirm the code is actually there. Nothing in the code checks this.
- [ ] GitHub URL starts with the **exact lowercase literal
      `https://github.com/`** (line 138). No `www.`, no `http://`, no SSH form.
- [ ] **Video plays from a logged-out browser**, is the right video, and is
      within 3 minutes (lines 44–45). Nothing checks any of this. An
      "unlisted" YouTube link is fine; a "private" one is not.
- [ ] Video URL pasted is the **share/watch URL**, not a studio/editor URL.

At the tab:

- [ ] Signed into the Space with the **designated submitting account** — the one
      attempt is spent against that HF username (lines 133–135, `utils.py:32-34`).
- [ ] Status card reads **`0 / 1 submissions used`** *before* you fill anything
      in. If it reads "Track 2 submission received", stop — the attempt is gone.
- [ ] **Team / Display Name** filled in with the exact team name. Blank
      publishes your slugified username (line 131).
- [ ] **Notes for judges** — optional, free text, stored in the JSON (line 160).
      Use it or leave it, but decide deliberately.
- [ ] All four inputs re-read once, by eye, against the artifacts above.

After clicking Submit:

- [ ] You get the **`## Track 2 submission received ✓`** block (lines 167–179).
      That block, not the status card, is the confirmation.
- [ ] The block's echoed **GitHub URL, video URL, and report filename** match
      what you intended (lines 176–178). Screenshot it.
- [ ] **Ignore the `0 / 1` card** — it is not refreshed by the submit button
      (lines 225–230) and will still say zero. Re-select the tab to refresh it
      if you want independent confirmation; it should then read "Track 2
      submission received".
- [ ] **Do not click Submit again.** Record the outcome in
      `docs/submission-ledger.md`.

---

## 12. Still not verifiable from here

- **Upload size limits.** None are set in the Space source
  (`app.py:117-123` passes no `max_file_size`), so whatever cap exists belongs
  to Gradio 6.24.0 / HF Spaces infrastructure and to the leaderboard dataset's
  own `.gitattributes` — that dataset is private and unreadable. Treated as a
  risk in §10.1 rather than a known number.
- **Whether `gr.File` preserves the original basename.** Lines 150–152 take
  `Path(...).name` of whatever Gradio hands back. Gradio's own temp-file naming
  is not part of the Space source and was not verified; the filename guidance in
  §11 is written to be safe either way.
- **The judging panel's actual reading of the report.** Rubric weights are
  stated (`tabs/rules.py:76-80`); nothing about how a report is scored is in
  code, because Track 2 has no scorer.
- **Whether `CHALLENGE_ACTIVE` will be flipped before or after the announced
  close.** It is a manual edit to `config.py:44`.
