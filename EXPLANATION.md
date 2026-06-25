# Program 2 (B) — Exon/Intron Splicing using a Hidden Markov Model

This program takes a raw nucleotide sequence and splits it into **exons**
(coding regions, kept) and **introns** (non-coding regions, spliced out),
mimicking the splicing step of mRNA processing. It then displays the
"spliced mRNA" (exons concatenated).

## 1. The HMM structure

### States (`THMMState`)

```
hsExon            — inside an exon
hsIntronStart1..3  — matching the 3-base donor splice site "GTA"
hsIntron          — inside an intron, scanning for the acceptor site
hsIntronEnd1..3    — matching the 3-base acceptor splice site "CAG"
```

This is a **finite-state automaton shaped like an HMM with deterministic
("Mealy-style") emissions** — instead of probabilistic emission/transition
matrices like Program 1, the states here directly correspond to recognizing
fixed biological motifs:

- **Donor site** (exon→intron boundary): canonical splice donor consensus
  `5' GT...` — here simplified/fixed to the literal 3-mer **`GTA`**.
  - `hsExon` --(`G`)--> `hsIntronStart1` --(`T`)--> `hsIntronStart2` --(`A`)--> `hsIntronStart3`
- **Acceptor site** (intron→exon boundary): canonical splice acceptor
  consensus `...AG 3'` — here simplified/fixed to the literal 3-mer **`CAG`**.
  - `hsIntron` --(`C`)--> `hsIntronEnd1` --(`A`)--> `hsIntronEnd2` --(`G`)--> `hsIntronEnd3`

Each state transition is **deterministic given the next nucleotide** — this
is effectively a Viterbi path through an HMM whose emission probabilities are
0 or 1 (a base either matches the expected motif character or it doesn't), so
the "decoding" collapses to simple pattern matching with backtracking-like
resets (see below). This is a common simplification of real splice-site HMMs
(which in practice use weight matrices/profile HMMs over a window of bases
with non-zero probabilities for all 4 bases) into a deterministic recognizer
for a teaching example.

### "Transition probabilities" (implicit)

There's no explicit probability matrix in this program (unlike Program 1) —
transitions are encoded directly in the `case State of` logic of `RunHMM`:

- From `hsExon`: any base keeps you in `hsExon`, **except** `G`, which begins
  a donor-site attempt (`hsIntronStart1`).
- From `hsIntronStart1`/`hsIntronStart2`: if the next expected base of `GTA`
  matches, advance; otherwise **abort back to `hsExon`** (but re-check
  whether the current base is itself a new `G`, to avoid missing an
  immediately-following donor site — see `IntronStart` reuse in the `else`
  branches).
- From `hsIntronStart3` (i.e., `GTA` fully matched): the donor site is
  confirmed, so the preceding bases since `ExonStart` are committed as an
  exon segment, and the state moves to `hsIntron` to scan for the acceptor.
- From `hsIntron`: any base stays in `hsIntron` except `C`, which begins an
  acceptor-site attempt (`hsIntronEnd1`).
- From `hsIntronEnd1`/`hsIntronEnd2`: same match-or-abort-to-`hsIntron` logic
  as the donor side, for `CAG`.
- From `hsIntronEnd3` (i.e., `CAG` fully matched): the intron (from
  `IntronStart` through current position) is committed as an intron segment,
  and the state returns to `hsExon`, with a new exon beginning right after.

### "Emission probabilities" (implicit)

Equivalent to a deterministic 0/1 emission model:
- `hsExon`/`hsIntron` "emit" any of A/C/G/T with equal likelihood (no
  preference) **except** for the specific trigger base (`G` or `C`
  respectively) that has special handling.
- `hsIntronStart1..3` and `hsIntronEnd1..3` only "accept" exactly one
  specific base each (`T`,`A` for donor; `A`,`G` for acceptor) — any other
  base is a mismatch causing reset.

### Initial state and termination

- The scan always starts in `hsExon` at position 1 (`ExonStart := 1`).
- At the end of the sequence, whatever partial state the automaton is left in
  is flushed as a final exon or intron segment (handled in the final `case
  State of` block after the main loop) — e.g. if the sequence ends mid-exon,
  or even mid-attempt at a donor site (`hsIntronStart1/2`), the partial bases
  are still counted as exon.

## 2. Code walkthrough

### Types
- `TSegmentKind` — `skExon` or `skIntron`.
- `TSegment` — one identified segment: start/end index (1-based, inclusive),
  kind, and the actual substring.

### Fields (`TForm1`)
- `FSegments` / `FSegCount` — dynamic array of all segments found by `RunHMM`,
  in sequence order.
- `FOutputTags` / `FOutputTagCnt` — parallel array to `lbOutput.Items`,
  recording which "color tag" (`TAG_HEADER`/`TAG_EXON`/`TAG_INTRON`) each
  output line should be drawn with (used by the owner-draw handlers for
  color-coded display).

### `FormCreate`
Pre-populates the input list box with a hardcoded example sequence
(comma-separated nucleotides) so the program has a working demo out of the
box.

### `CleanSequence`
Strips everything except `A/C/G/T` (uppercased) from the raw input — this is
what allows the comma-separated demo string to be parsed into a clean
sequence.

### `AddSegment`
Appends a `TSegment` to `FSegments` (growing the array in chunks of 32),
skipping zero/negative-length segments.

### `RunHMM(Seq)`
The core splicing automaton described above. Walks the sequence once,
character by character, maintaining `State`, `ExonStart` (start of the
current pending exon), `IntronStart` (start of the current pending intron,
also doubles as where the donor `GTA` match began), and emits a segment via
`AddSegment` every time a donor site (exon end) or acceptor site (intron end)
is fully confirmed. The trailing `case State of` after the loop flushes any
final partial segment.

### `DisplayResults(Seq)`
Builds the human-readable output:
- Header lines (sequence length, donor/acceptor motif legend).
- For every segment in order: an `EXON #n` or `INTRON #n` line with
  position/length, followed by the actual subsequence — each tagged for
  color-coded rendering.
- A summary line with total exon/intron counts.
- Separately, builds `lbSpliced` (the final "mRNA"): exons concatenated in
  order (introns removed), with its length.
- Uses `AddOutputLine` to keep `lbOutput.Items` and `FOutputTags` in sync.

### Owner-draw handlers
- `lbOutputDrawItem` — colors each line of the output list box according to
  its tag: orange/tan for exons (`CLR_EXON_BG`), light blue for introns
  (`CLR_INTRON_BG`), gray for headers (`CLR_HDR_BG`), with normal selection
  highlighting preserved.
- `lbSplicedDrawItem` — simpler version for the "spliced mRNA" box; hardcodes
  index 2 (the actual sequence line, after the two header lines) to be drawn
  with the exon color.

### GUI event handlers
- `btnRunClick` — concatenates all lines from `lbInput`, cleans the
  sequence, runs `RunHMM`, then `DisplayResults`. Shows an error message if
  no valid bases remain after cleaning.
- `btnClearClick` — clears all three list boxes and resets segment/tag
  counters.
