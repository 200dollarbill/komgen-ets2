# Design: Noncoding / Start-codon / Stop-codon states for the splicing HMM

## Context

`unit1.pas` currently models exon/intron splicing with an 8-state deterministic
finite-state automaton (`THMMState`): `hsExon`, `hsIntronStart1..3`,
`hsIntron`, `hsIntronEnd1..3`. It treats the *entire* input sequence as coding
(exon/intron only) starting from position 1, with no concept of untranslated
regions or where translation actually begins/ends.

The professor's feedback: a real intron-snipper needs to account for the fact
that genes have noncoding flanking regions, and that the coding region is
bounded by a start codon (ATG) and a stop codon (TAA/TAG/TGA) — not just
"whatever is left after removing introns."

This design adds 7 new states (Noncoding ×1, Start-codon ×3, Stop-codon ×3) to
the existing 8, for **15 states total**, wrapping the existing exon/intron
machine inside a new outer Noncoding → Start → Coding → Stop → Noncoding
cycle that can repeat for multiple ORFs in one sequence.

## Decisions (confirmed with user)

1. **Literal motif matching**, not true reading-frame tracking across spliced
   exons. Start/Stop codons are detected as literal 3-character substrings
   while accumulating exon bases — same simplification level the existing
   code already uses for the GTA donor / CAG acceptor sites. A stop codon
   split across an intron boundary is *not* detected; this matches the
   existing program's documented teaching-level simplification.
2. **Spliced mRNA includes Start and Stop codon bases**, not just exon
   bases — biologically those bases ARE coding sequence. Introns and
   noncoding/UTR regions remain excluded.
3. **Multiple ORFs per sequence are supported** — after a Stop codon is
   confirmed, the machine returns to `hsNoncoding` and can begin a new
   Start→Coding→Stop cycle if another ATG appears later in the sequence.

## State machine

```pascal
THMMState = (
  hsNoncoding,
  hsStartC1, hsStartC2, hsStartC3,
  hsExon, hsIntronStart1, hsIntronStart2, hsIntronStart3,
  hsIntron, hsIntronEnd1, hsIntronEnd2, hsIntronEnd3,
  hsStopC1, hsStopC2, hsStopC3
);
```

### Transitions

- **`hsNoncoding`**: any base stays; `A` begins a start-codon attempt →
  `hsStartC1`.
- **`hsStartC1`** (saw `A`): `T` → `hsStartC2`; else abort to `hsNoncoding`
  (rechecking current base for a new `A` trigger, same abort-and-recheck
  pattern as the existing donor/acceptor logic).
- **`hsStartC2`** (saw `AT`): `G` → `hsStartC3`; else abort to `hsNoncoding`
  (recheck for `A`).
- **`hsStartC3`** (saw `ATG`, confirmed): commit any pending noncoding
  segment (`skNoncoding`), commit `ATG` as a `skStartCodon` segment, set
  `ExonStart` to the next position, transition to `hsExon`.
- **`hsExon`** (now doubles as "currently in a coding region, accumulating
  exon bases"): `G` begins a donor attempt → `hsIntronStart1` (unchanged);
  **new:** `T` begins a stop-codon attempt → `hsStopC1`. Both checks coexist
  since they key off different first letters — no ambiguity.
- **`hsIntronStart1..3` / `hsIntron` / `hsIntronEnd1..3`**: unchanged from
  today. Stop-codon scanning does not apply inside introns.
- **`hsStopC1`** (saw `T`): `A` or `G` → `hsStopC2`; else abort to `hsExon`
  (recheck current base for `G`/`T` triggers).
- **`hsStopC2`** (saw `TA` or `TG`): if buffer is `TA`, accept `A` or `G`
  (completing `TAA`/`TAG`); if buffer is `TG`, accept only `A` (completing
  `TGA`); else abort to `hsExon` (recheck current base).
- **`hsStopC3`** (confirmed `TAA`/`TAG`/`TGA`): commit pending exon segment
  (`skExon`), commit the 3 bases as a `skStopCodon` segment, transition to
  `hsNoncoding` with a new noncoding region starting at the next position.

### End-of-sequence flush

Extends the existing trailing `case State of` block: whatever partial
state/segment remains when the sequence ends is flushed as the best-fit
segment kind — `hsNoncoding`/mid-start-attempt → flush pending noncoding;
`hsExon`/intron sub-states → flush pending exon/intron (as today);
mid-stop-attempt → flush pending exon (stop not confirmed, so those bases
revert to exon, same "abort" semantics used elsewhere).

## Data structures

```pascal
TSegmentKind = (skNoncoding, skStartCodon, skExon, skIntron, skStopCodon);
```

New display tag constants `TAG_NONCODING`, `TAG_START`, `TAG_STOP` with new
`CLR_*` background colors, added alongside the existing `TAG_HEADER`,
`TAG_EXON`, `TAG_INTRON`. `lbOutputDrawItem`'s tag→color `case` extends to
cover the two new tags.

`DisplayResults` gains `NONCODING #n`, `START CODON #n`, `STOP CODON #n`
output lines in the same `[pos x..y] len=z` format as existing EXON/INTRON
lines, in sequence order alongside them.

## Spliced mRNA

`SplicedExon` accumulation in `DisplayResults` appends `skStartCodon`,
`skExon`, and `skStopCodon` segments (in order); `skIntron` and
`skNoncoding` segments are excluded, same exclusion logic as today just
extended to the two new non-coding-adjacent kinds.

## Stats (`ComputeProbStats`)

Add `P_START` and `P_STOP` prior-cost constants (same role as `P_DONOR` /
`P_ACCEPTOR`), each charged once per confirmed Start/Stop codon segment, and
fold into `PriorCost`/`TotalLL`. Noncoding bases contribute 0 to
log-likelihood (no emission model assigned to noncoding region in this
design). `NullLL` (whole-sequence-as-exon baseline) is unchanged.

## Out of scope

- True reading-frame tracking across intron-spliced exon boundaries.
- Reverse-strand / reverse-complement scanning.
- Updating `hmm_diagram.svg` (can be regenerated separately if wanted).
- Updating `EXPLANATION.md` narrative doc (can be done as a follow-up).
