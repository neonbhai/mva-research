# Reviewer brief — privacy and security

You are a privacy engineer reviewing `mva-research` adversarially. Read
`src/mva/privacy/`, `src/mva/config.py`, `.gitignore`, `docs/privacy-model.md`.
Do NOT modify files. Do NOT create any file containing realistic patient data.

## What to attack
1. **The scanner's own output (GP-41).** This is the highest-value target: the
   scanner is run by an agent and its output enters model context. Find ANY path
   by which matched content, a filename fragment carrying a sample ID, or a
   decode error message carrying bytes reaches stdout, a log, or a returned
   object. `UnicodeDecodeError` embeds the offending bytes — is every read
   wrapped?
2. **`.gitignore` correctness.** Verify with `git check-ignore -v --no-index`.
   Remember exit 0 means "a rule matched", including a **negation** — a check that
   only tests the exit code is wrong. Are any negations too broad? Can a genomic
   file reach a negated path?
3. **Tracked files.** `.gitignore` does not protect already-tracked files. Is
   there a check for that, and does it work?
4. **Workspace containment.** Try to defeat it: symlinks, hardlinks, `..`,
   relative paths, a symlinked intermediate directory, case-insensitive HFS+
   collisions.
5. **The log filter (GP-42).** Is it attached to every handler, not just loggers?
   Records propagating from child loggers skip ancestor filters — is that hole
   closed? Are `exc_info`, `exc_text` and `stack_info` scrubbed?
6. **Network denial.** Are the stated limits honest? Verify the audit hook gates
   correctly and cannot be defeated by `from socket import socket`.
7. **The export gate (GP-43).** Can a SENSITIVE artifact reach a public path? Is
   the post-render re-scan actually run, or is classification trusted alone?
8. **Regex false negatives.** For each rule, construct patient data it would MISS.
   Under-detection matters more here than over-detection.
9. **Error and warning messages** across the whole codebase — do any interpolate a
   record, a genotype, a sample ID or a workspace path?

## Output
Severity-ranked findings with file:line, the concrete disclosure scenario, and the
fix. Mark anything that could disclose real patient data as CRITICAL regardless of
likelihood. Every valid finding must map to a test.
