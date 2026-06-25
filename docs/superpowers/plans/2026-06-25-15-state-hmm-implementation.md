# 15-State HMM Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend the splicing FSM in this Lazarus project from 8 states to the 15-state design (Noncoding / Start-codon / Coding-core / Stop-codon) from `docs/superpowers/specs/2026-06-25-start-stop-codon-states-design.md`, with automated tests for the core engine.

**Architecture:** Extract the FSM core (states, segment types, `RunHMM`, `ComputeProbStats`) out of `unit1.pas` into a new dependency-free unit `HMMCore.pas` (no LCL imports), so it can be exercised by a plain console test program compiled with `fpc` directly (no GUI toolchain needed for tests). `unit1.pas` keeps only GUI wiring (display, colors, button handlers) and calls into `HMMCore`.

**Tech Stack:** Free Pascal 3.2.2 (`fpc`), Lazarus 4.2 LCL (`lazbuild` at `/home/pingwalk/lazarus-4.2-0/lazarus/lazbuild`), no external test framework — a console program with manual assertions serves as the test runner.

## Global Constraints

- Literal motif matching only (no true reading-frame tracking across spliced exons) — per the approved design.
- Spliced mRNA = StartCodon + Exon segments + StopCodon segments, in order (Introns and Noncoding excluded).
- Multiple ORFs per sequence must be supported (cycle repeats after a Stop codon).
- The existing 8-state exon/intron transition logic must be preserved byte-for-byte in behavior (same donor/acceptor matching), just relocated.

---

### Task 1: Core FSM engine (`HMMCore.pas`) with console tests

**Files:**
- Create: `HMMCore.pas`
- Create: `test_hmmcore.pas` (standalone console test program, plain `fpc`, no LCL)

**Interfaces:**
- Produces: `THMMState` (15 values), `TSegmentKind = (skNoncoding, skStartCodon, skExon, skIntron, skStopCodon)`, `TSegment` record (`StartIdx, EndIdx: Integer; Kind: TSegmentKind; Sequence: string`), `TSegmentArray = array of TSegment`.
- Produces: `procedure RunHMM(const Seq: string; var Segments: TSegmentArray; out SegCount: Integer);`
- Produces: `procedure ComputeProbStats(const Seq: string; const Segments: TSegmentArray; SegCount: Integer; out ExonLL, IntronLL, PriorCost, TotalLL, NullLL: Double);`
- Produces: constants `NEG_INF, P_DONOR, P_ACCEPTOR, P_START, P_STOP`, arrays `EmitExon, EmitIntron`.

- [ ] **Step 1: Write the test program (will fail to compile — `HMMCore` doesn't exist yet)**

Create `test_hmmcore.pas`:

```pascal
program test_hmmcore;

{$mode objfpc}{$H+}

uses
  SysUtils, HMMCore;

var
  FailCount: Integer = 0;

procedure Check(Cond: Boolean; const Msg: string);
begin
  if Cond then
    WriteLn('PASS: ', Msg)
  else
  begin
    WriteLn('FAIL: ', Msg);
    Inc(FailCount);
  end;
end;

procedure CheckSeg(const Segs: TSegmentArray; Idx, SegCount: Integer;
  ExpKind: TSegmentKind; const ExpSeq: string; const Label_: string);
begin
  if Idx >= SegCount then
  begin
    WriteLn('FAIL: ', Label_, ' (segment index ', Idx, ' out of range, only ', SegCount, ' segments)');
    Inc(FailCount);
    Exit;
  end;
  Check(Segs[Idx].Kind = ExpKind, Label_ + ' kind');
  Check(Segs[Idx].Sequence = ExpSeq, Label_ + ' sequence (got "' + Segs[Idx].Sequence + '", want "' + ExpSeq + '")');
end;

procedure Test_PureNoncoding;
var
  Segs: TSegmentArray;
  Cnt: Integer;
begin
  WriteLn('-- Test_PureNoncoding --');
  RunHMM('CCCCCCGGGGTTTT', Segs, Cnt);
  Check(Cnt = 1, 'segment count = 1');
  CheckSeg(Segs, 0, Cnt, skNoncoding, 'CCCCCCGGGGTTTT', 'whole-seq noncoding');
end;

procedure Test_FullGene;
var
  Segs: TSegmentArray;
  Cnt: Integer;
begin
  WriteLn('-- Test_FullGene --');
  // CC | ATG | AAA | GTACCCCCAG (intron, GTA..CAG) | CC | TAA | GG
  RunHMM('CCATGAAAGTACCCCCAGCCTAAGG', Segs, Cnt);
  Check(Cnt = 6, 'segment count = 6 (got ' + IntToStr(Cnt) + ')');
  CheckSeg(Segs, 0, Cnt, skNoncoding,  'CC',            'leading UTR');
  CheckSeg(Segs, 1, Cnt, skStartCodon, 'ATG',           'start codon');
  CheckSeg(Segs, 2, Cnt, skExon,       'AAA',           'exon 1');
  CheckSeg(Segs, 3, Cnt, skIntron,     'GTACCCCCAG',    'intron');
  CheckSeg(Segs, 4, Cnt, skExon,       'CC',             'exon 2');
  CheckSeg(Segs, 5, Cnt, skStopCodon,  'TAA',           'stop codon');
end;

procedure Test_ZeroLengthExonBetweenStartAndStop;
var
  Segs: TSegmentArray;
  Cnt: Integer;
begin
  WriteLn('-- Test_ZeroLengthExonBetweenStartAndStop --');
  RunHMM('ATGTAA', Segs, Cnt);
  Check(Cnt = 2, 'segment count = 2 (no zero-length exon emitted)');
  CheckSeg(Segs, 0, Cnt, skStartCodon, 'ATG', 'start codon');
  CheckSeg(Segs, 1, Cnt, skStopCodon,  'TAA', 'stop codon');
end;

procedure Test_SequenceEndsExactlyAtStartCodon;
var
  Segs: TSegmentArray;
  Cnt: Integer;
begin
  WriteLn('-- Test_SequenceEndsExactlyAtStartCodon --');
  RunHMM('CCATG', Segs, Cnt);
  Check(Cnt = 2, 'segment count = 2 (got ' + IntToStr(Cnt) + ')');
  CheckSeg(Segs, 0, Cnt, skNoncoding,  'CC',  'leading UTR');
  CheckSeg(Segs, 1, Cnt, skStartCodon, 'ATG', 'start codon flushed at EOS');
end;

procedure Test_SequenceEndsExactlyAtStopCodon;
var
  Segs: TSegmentArray;
  Cnt: Integer;
begin
  WriteLn('-- Test_SequenceEndsExactlyAtStopCodon --');
  RunHMM('ATGAAATAG', Segs, Cnt);
  Check(Cnt = 3, 'segment count = 3 (got ' + IntToStr(Cnt) + ')');
  CheckSeg(Segs, 0, Cnt, skStartCodon, 'ATG', 'start codon');
  CheckSeg(Segs, 1, Cnt, skExon,       'AAA', 'exon');
  CheckSeg(Segs, 2, Cnt, skStopCodon,  'TAG', 'stop codon flushed at EOS');
end;

procedure Test_MultipleOrfsRepeatCycle;
var
  Segs: TSegmentArray;
  Cnt: Integer;
begin
  WriteLn('-- Test_MultipleOrfsRepeatCycle --');
  // ORF1: ATG AAA TAA, then noncoding CC, then ORF2: ATG CCC TGA
  RunHMM('ATGAAATAACCATGCCCTGA', Segs, Cnt);
  Check(Cnt = 7, 'segment count = 7 (got ' + IntToStr(Cnt) + ')');
  CheckSeg(Segs, 0, Cnt, skStartCodon, 'ATG', 'ORF1 start');
  CheckSeg(Segs, 1, Cnt, skExon,       'AAA', 'ORF1 exon');
  CheckSeg(Segs, 2, Cnt, skStopCodon,  'TAA', 'ORF1 stop');
  CheckSeg(Segs, 3, Cnt, skNoncoding,  'CC',  'inter-ORF UTR');
  CheckSeg(Segs, 4, Cnt, skStartCodon, 'ATG', 'ORF2 start');
  CheckSeg(Segs, 5, Cnt, skExon,       'CCC', 'ORF2 exon');
  CheckSeg(Segs, 6, Cnt, skStopCodon,  'TGA', 'ORF2 stop');
end;

procedure Test_FalseStartAbortsBackToNoncoding;
var
  Segs: TSegmentArray;
  Cnt: Integer;
begin
  WriteLn('-- Test_FalseStartAbortsBackToNoncoding --');
  // "AT" then 'C' (not G) aborts the start attempt; whole thing stays noncoding
  RunHMM('ATCATG', Segs, Cnt);
  Check(Cnt = 2, 'segment count = 2 (got ' + IntToStr(Cnt) + ')');
  CheckSeg(Segs, 0, Cnt, skNoncoding,  'ATC', 'failed start attempt stays noncoding');
  CheckSeg(Segs, 1, Cnt, skStartCodon, 'ATG', 'real start codon follows');
end;

procedure Test_FalseStopAbortsBackToExon;
var
  Segs: TSegmentArray;
  Cnt: Integer;
begin
  WriteLn('-- Test_FalseStopAbortsBackToExon --');
  // After ATG: "TCC" looks like a stop attempt (T) but fails (C is not A/G) -> stays exon
  // then "TGA" really stops.
  RunHMM('ATGTCCTGA', Segs, Cnt);
  Check(Cnt = 3, 'segment count = 3 (got ' + IntToStr(Cnt) + ')');
  CheckSeg(Segs, 0, Cnt, skStartCodon, 'ATG', 'start codon');
  CheckSeg(Segs, 1, Cnt, skExon,       'TCC', 'failed stop attempt stays exon');
  CheckSeg(Segs, 2, Cnt, skStopCodon,  'TGA', 'real stop codon');
end;

procedure Test_ProbStatsIncludeStartStopPriors;
var
  Segs: TSegmentArray;
  Cnt: Integer;
  ExonLL, IntronLL, PriorCost, TotalLL, NullLL: Double;
begin
  WriteLn('-- Test_ProbStatsIncludeStartStopPriors --');
  RunHMM('ATGAAATAG', Segs, Cnt);
  ComputeProbStats('ATGAAATAG', Segs, Cnt, ExonLL, IntronLL, PriorCost, TotalLL, NullLL);
  Check(Abs(PriorCost - (Ln(P_START) + Ln(P_STOP))) < 1e-9,
    'PriorCost = ln(P_START) + ln(P_STOP) (got ' + FloatToStr(PriorCost) + ')');
end;

begin
  Test_PureNoncoding;
  Test_FullGene;
  Test_ZeroLengthExonBetweenStartAndStop;
  Test_SequenceEndsExactlyAtStartCodon;
  Test_SequenceEndsExactlyAtStopCodon;
  Test_MultipleOrfsRepeatCycle;
  Test_FalseStartAbortsBackToNoncoding;
  Test_FalseStopAbortsBackToExon;
  Test_ProbStatsIncludeStartStopPriors;

  WriteLn('');
  if FailCount = 0 then
    WriteLn('ALL TESTS PASSED')
  else
  begin
    WriteLn(FailCount, ' TEST(S) FAILED');
    Halt(1);
  end;
end.
```

- [ ] **Step 2: Run it to confirm it fails (HMMCore unit missing)**

Run: `cd "/mnt/Data/yep/Kuliah/Tugas/Semester 6 bosen/Komgen Finals/Komgen ETS/Program 2 (B)" && fpc test_hmmcore.pas`
Expected: compile error, `Fatal: Cannot find HMMCore used by test_hmmcore`

- [ ] **Step 3: Write `HMMCore.pas`**

Create `HMMCore.pas`:

```pascal
unit HMMCore;

{$mode objfpc}{$H+}

interface

uses
  SysUtils;

type
  THMMState = (
    hsNoncoding,
    hsStartC1, hsStartC2, hsStartC3,
    hsExon, hsIntronStart1, hsIntronStart2, hsIntronStart3,
    hsIntron, hsIntronEnd1, hsIntronEnd2, hsIntronEnd3,
    hsStopC1, hsStopC2, hsStopC3
  );

  TSegmentKind = (skNoncoding, skStartCodon, skExon, skIntron, skStopCodon);

  TSegment = record
    StartIdx : Integer;
    EndIdx   : Integer;
    Kind     : TSegmentKind;
    Sequence : string;
  end;

  TSegmentArray = array of TSegment;

const
  NEG_INF = -1e300;
  P_DONOR    = 0.02;
  P_ACCEPTOR = 0.05;
  P_START    = 1.0 / 64.0;  // chance of 3 random bases matching ATG
  P_STOP     = 3.0 / 64.0;  // chance of 3 random bases matching TAA/TAG/TGA

  EmitExon   : array[0..3] of Double = (0.25, 0.25, 0.25, 0.25);
  EmitIntron : array[0..3] of Double = (0.32, 0.32, 0.18, 0.18);

function SafeLog(const X: Double): Double;
function NtIndex(C: Char): Integer;
procedure AddSeg(var Segments: TSegmentArray; var SegCount: Integer;
  AStart, AEnd: Integer; AKind: TSegmentKind; const ASeq: string);
procedure RunHMM(const Seq: string; var Segments: TSegmentArray; out SegCount: Integer);
procedure ComputeProbStats(const Seq: string; const Segments: TSegmentArray;
  SegCount: Integer; out ExonLL, IntronLL, PriorCost, TotalLL, NullLL: Double);

implementation

function SafeLog(const X: Double): Double;
begin
  if X <= 0.0 then Result := NEG_INF
  else Result := Ln(X);
end;

function NtIndex(C: Char): Integer;
begin
  case UpCase(C) of
    'A': Result := 0;
    'T': Result := 1;
    'C': Result := 2;
    'G': Result := 3;
  else
    Result := -1;
  end;
end;

procedure AddSeg(var Segments: TSegmentArray; var SegCount: Integer;
  AStart, AEnd: Integer; AKind: TSegmentKind; const ASeq: string);
begin
  if AEnd < AStart then Exit;
  if SegCount >= Length(Segments) then
    SetLength(Segments, SegCount + 32);
  Segments[SegCount].StartIdx := AStart;
  Segments[SegCount].EndIdx   := AEnd;
  Segments[SegCount].Kind     := AKind;
  Segments[SegCount].Sequence := ASeq;
  Inc(SegCount);
end;

procedure RunHMM(const Seq: string; var Segments: TSegmentArray; out SegCount: Integer);
var
  State          : THMMState;
  i              : Integer;
  nt             : Char;
  NoncodingStart : Integer;
  StartPos       : Integer;
  ExonStart      : Integer;
  IntronStart    : Integer;
  StopStart      : Integer;
  DonorBuf       : string;
  AccBuf         : string;
  StartBuf       : string;
  StopBuf        : string;
begin
  SegCount := 0;
  SetLength(Segments, 64);
  if Length(Seq) = 0 then Exit;

  State          := hsNoncoding;
  NoncodingStart := 1;
  StartPos       := 0;
  ExonStart      := 0;
  IntronStart    := 0;
  StopStart      := 0;
  DonorBuf       := '';
  AccBuf         := '';
  StartBuf       := '';
  StopBuf        := '';

  i := 1;
  while i <= Length(Seq) do
  begin
    nt := Seq[i];
    case State of

      hsNoncoding:
      begin
        if nt = 'A' then
        begin
          StartBuf := 'A';
          StartPos := i;
          State    := hsStartC1;
        end;
      end;

      hsStartC1:
      begin
        if nt = 'T' then
        begin
          StartBuf := StartBuf + 'T';
          State    := hsStartC2;
        end
        else
        begin
          State    := hsNoncoding;
          StartBuf := '';
          if nt = 'A' then
          begin
            StartBuf := 'A';
            StartPos := i;
            State    := hsStartC1;
          end;
        end;
      end;

      hsStartC2:
      begin
        if nt = 'G' then
        begin
          StartBuf := StartBuf + 'G';
          State    := hsStartC3;
        end
        else
        begin
          State    := hsNoncoding;
          StartBuf := '';
          if nt = 'A' then
          begin
            StartBuf := 'A';
            StartPos := i;
            State    := hsStartC1;
          end;
        end;
      end;

      hsStartC3:
      begin
        if StartPos - 1 >= NoncodingStart then
          AddSeg(Segments, SegCount, NoncodingStart, StartPos - 1, skNoncoding,
                 Copy(Seq, NoncodingStart, StartPos - NoncodingStart));
        AddSeg(Segments, SegCount, StartPos, StartPos + 2, skStartCodon,
               Copy(Seq, StartPos, 3));
        ExonStart := StartPos + 3;
        StartBuf  := '';
        State     := hsExon;
        Continue;
      end;

      hsExon:
      begin
        if nt = 'G' then
        begin
          DonorBuf    := 'G';
          IntronStart := i;
          State       := hsIntronStart1;
        end
        else if nt = 'T' then
        begin
          StopBuf   := 'T';
          StopStart := i;
          State     := hsStopC1;
        end;
      end;

      hsIntronStart1:
      begin
        if nt = 'T' then
        begin
          DonorBuf := DonorBuf + 'T';
          State    := hsIntronStart2;
        end
        else
        begin
          State    := hsExon;
          DonorBuf := '';
          if nt = 'G' then
          begin
            DonorBuf    := 'G';
            IntronStart := i;
            State       := hsIntronStart1;
          end
          else if nt = 'T' then
          begin
            StopBuf   := 'T';
            StopStart := i;
            State     := hsStopC1;
          end;
        end;
      end;

      hsIntronStart2:
      begin
        if nt = 'A' then
        begin
          DonorBuf := DonorBuf + 'A';
          State    := hsIntronStart3;
        end
        else
        begin
          State    := hsExon;
          DonorBuf := '';
          if nt = 'G' then
          begin
            DonorBuf    := 'G';
            IntronStart := i;
            State       := hsIntronStart1;
          end
          else if nt = 'T' then
          begin
            StopBuf   := 'T';
            StopStart := i;
            State     := hsStopC1;
          end;
        end;
      end;

      hsIntronStart3:
      begin
        if IntronStart - 1 >= ExonStart then
          AddSeg(Segments, SegCount, ExonStart, IntronStart - 1, skExon,
                 Copy(Seq, ExonStart, IntronStart - ExonStart));
        AccBuf := '';
        State  := hsIntron;
        Continue;
      end;

      hsIntron:
      begin
        if nt = 'C' then
        begin
          AccBuf := 'C';
          State  := hsIntronEnd1;
        end;
      end;

      hsIntronEnd1:
      begin
        if nt = 'A' then
        begin
          AccBuf := AccBuf + 'A';
          State  := hsIntronEnd2;
        end
        else
        begin
          AccBuf := '';
          State  := hsIntron;
          if nt = 'C' then
          begin
            AccBuf := 'C';
            State  := hsIntronEnd1;
          end;
        end;
      end;

      hsIntronEnd2:
      begin
        if nt = 'G' then
        begin
          AccBuf := AccBuf + 'G';
          State  := hsIntronEnd3;
        end
        else
        begin
          AccBuf := '';
          State  := hsIntron;
          if nt = 'C' then
          begin
            AccBuf := 'C';
            State  := hsIntronEnd1;
          end;
        end;
      end;

      hsIntronEnd3:
      begin
        AddSeg(Segments, SegCount, IntronStart, i, skIntron,
               Copy(Seq, IntronStart, i - IntronStart + 1));
        ExonStart := i + 1;
        AccBuf    := '';
        State     := hsExon;
      end;

      hsStopC1:
      begin
        if (nt = 'A') or (nt = 'G') then
        begin
          StopBuf := StopBuf + nt;
          State   := hsStopC2;
        end
        else
        begin
          State   := hsExon;
          StopBuf := '';
          if nt = 'G' then
          begin
            DonorBuf    := 'G';
            IntronStart := i;
            State       := hsIntronStart1;
          end
          else if nt = 'T' then
          begin
            StopBuf   := 'T';
            StopStart := i;
            State     := hsStopC1;
          end;
        end;
      end;

      hsStopC2:
      begin
        if ((StopBuf = 'TA') and ((nt = 'A') or (nt = 'G'))) or
           ((StopBuf = 'TG') and (nt = 'A')) then
        begin
          StopBuf := StopBuf + nt;
          State   := hsStopC3;
        end
        else
        begin
          State   := hsExon;
          StopBuf := '';
          if nt = 'G' then
          begin
            DonorBuf    := 'G';
            IntronStart := i;
            State       := hsIntronStart1;
          end
          else if nt = 'T' then
          begin
            StopBuf   := 'T';
            StopStart := i;
            State     := hsStopC1;
          end;
        end;
      end;

      hsStopC3:
      begin
        if StopStart - 1 >= ExonStart then
          AddSeg(Segments, SegCount, ExonStart, StopStart - 1, skExon,
                 Copy(Seq, ExonStart, StopStart - ExonStart));
        AddSeg(Segments, SegCount, StopStart, i, skStopCodon,
               Copy(Seq, StopStart, i - StopStart + 1));
        NoncodingStart := i + 1;
        StopBuf        := '';
        State          := hsNoncoding;
      end;

    end;
    Inc(i);
  end;

  case State of
    hsNoncoding, hsStartC1, hsStartC2:
      if NoncodingStart <= Length(Seq) then
        AddSeg(Segments, SegCount, NoncodingStart, Length(Seq), skNoncoding,
               Copy(Seq, NoncodingStart, Length(Seq) - NoncodingStart + 1));
    hsStartC3:
      begin
        if StartPos - 1 >= NoncodingStart then
          AddSeg(Segments, SegCount, NoncodingStart, StartPos - 1, skNoncoding,
                 Copy(Seq, NoncodingStart, StartPos - NoncodingStart));
        AddSeg(Segments, SegCount, StartPos, StartPos + 2, skStartCodon,
               Copy(Seq, StartPos, 3));
      end;
    hsExon, hsIntronStart1, hsIntronStart2, hsStopC1, hsStopC2:
      if ExonStart <= Length(Seq) then
        AddSeg(Segments, SegCount, ExonStart, Length(Seq), skExon,
               Copy(Seq, ExonStart, Length(Seq) - ExonStart + 1));
    hsStopC3:
      begin
        if StopStart - 1 >= ExonStart then
          AddSeg(Segments, SegCount, ExonStart, StopStart - 1, skExon,
                 Copy(Seq, ExonStart, StopStart - ExonStart));
        AddSeg(Segments, SegCount, StopStart, Length(Seq), skStopCodon,
               Copy(Seq, StopStart, Length(Seq) - StopStart + 1));
      end;
    hsIntron, hsIntronEnd1, hsIntronEnd2:
      if IntronStart <= Length(Seq) then
        AddSeg(Segments, SegCount, IntronStart, Length(Seq), skIntron,
               Copy(Seq, IntronStart, Length(Seq) - IntronStart + 1));
    else ;
  end;
end;

procedure ComputeProbStats(const Seq: string; const Segments: TSegmentArray;
  SegCount: Integer; out ExonLL, IntronLL, PriorCost, TotalLL, NullLL: Double);
var
  i, j, NtIdx : Integer;
  seg         : TSegment;
begin
  ExonLL    := 0.0;
  IntronLL  := 0.0;
  PriorCost := 0.0;

  for i := 0 to SegCount - 1 do
  begin
    seg := Segments[i];
    for j := 1 to Length(seg.Sequence) do
    begin
      NtIdx := NtIndex(seg.Sequence[j]);
      if NtIdx < 0 then Continue;
      if (seg.Kind = skExon) or (seg.Kind = skStartCodon) or (seg.Kind = skStopCodon) then
        ExonLL := ExonLL + SafeLog(EmitExon[NtIdx])
      else if seg.Kind = skIntron then
        IntronLL := IntronLL + SafeLog(EmitIntron[NtIdx]);
    end;
    if seg.Kind = skIntron then
      PriorCost := PriorCost + SafeLog(P_DONOR) + SafeLog(P_ACCEPTOR);
    if seg.Kind = skStartCodon then
      PriorCost := PriorCost + SafeLog(P_START);
    if seg.Kind = skStopCodon then
      PriorCost := PriorCost + SafeLog(P_STOP);
  end;

  TotalLL := ExonLL + IntronLL + PriorCost;
  NullLL := 0.0;
  for i := 1 to Length(Seq) do
  begin
    NtIdx := NtIndex(Seq[i]);
    if NtIdx >= 0 then
      NullLL := NullLL + SafeLog(EmitExon[NtIdx]);
  end;
end;

end.
```

- [ ] **Step 4: Run the tests, fix until all pass**

Run: `cd "/mnt/Data/yep/Kuliah/Tugas/Semester 6 bosen/Komgen Finals/Komgen ETS/Program 2 (B)" && fpc test_hmmcore.pas && ./test_hmmcore`
Expected: `ALL TESTS PASSED`, exit code 0.

- [ ] **Step 5: Commit**

```bash
cd "/mnt/Data/yep/Kuliah/Tugas/Semester 6 bosen/Komgen Finals/Komgen ETS/Program 2 (B)"
git add HMMCore.pas test_hmmcore.pas
git commit -m "Add 15-state HMM core engine with console tests"
```

---

### Task 2: Wire `unit1.pas` (GUI) to `HMMCore`

**Files:**
- Modify: `unit1.pas`

**Interfaces:**
- Consumes: everything produced by Task 1 (`THMMState`, `TSegmentKind`, `TSegment`, `TSegmentArray`, `RunHMM`, `ComputeProbStats`, constants).

- [ ] **Step 1: Remove duplicated types/constants/functions from `unit1.pas`, add `HMMCore` to `uses`**

In the `interface` section, replace:
```pascal
uses
  Classes, SysUtils, Forms, Controls, Graphics, Dialogs, StdCtrls, ExtCtrls,
  ComCtrls, LCLType;

type
  THMMState = (hsExon, hsIntronStart1, hsIntronStart2, hsIntronStart3,
               hsIntron, hsIntronEnd1, hsIntronEnd2, hsIntronEnd3);

  TSegmentKind = (skExon, skIntron);

  TSegment = record
    StartIdx : Integer;
    EndIdx   : Integer;
    Kind     : TSegmentKind;
    Sequence : string;
  end;
```
with:
```pascal
uses
  Classes, SysUtils, Forms, Controls, Graphics, Dialogs, StdCtrls, ExtCtrls,
  ComCtrls, LCLType, HMMCore;
```
(the `TSegment`/`TSegmentKind`/`THMMState` declarations are now provided by `HMMCore` — delete them here).

- [ ] **Step 2: Remove `RunHMM`, `AddSegment`, `ComputeProbStats` methods from `TForm1`'s method list and constants that moved to `HMMCore`**

In the `TForm1` private section, replace:
```pascal
  private
    FSegments    : array of TSegment;
    FSegCount    : Integer;
    FOutputTags  : array of Integer;
    FOutputTagCnt: Integer;
    function  CleanSequence(const Raw: string): string;
    procedure RunHMM(const Seq: string);
    procedure DisplayResults(const Seq: string);
    procedure AddSegment(AStart, AEnd: Integer; AKind: TSegmentKind;
                         const ASeq: string);
    procedure AddOutputLine(const S: string; ATag: Integer);
    procedure ComputeProbStats(const Seq: string; out ExonLL, IntronLL,
                               PriorCost, TotalLL, NullLL: Double);
```
with:
```pascal
  private
    FSegments    : TSegmentArray;
    FSegCount    : Integer;
    FOutputTags  : array of Integer;
    FOutputTagCnt: Integer;
    function  CleanSequence(const Raw: string): string;
    procedure DisplayResults(const Seq: string);
    procedure AddOutputLine(const S: string; ATag: Integer);
```

In the `const` section, replace:
```pascal
const
  TAG_HEADER = 0;
  TAG_EXON   = 1;
  TAG_INTRON = 2;

  CLR_EXON_BG   = $00FFE8B0;
  CLR_INTRON_BG = $00C8E8FF;
  CLR_HDR_BG    = $00E8E8E8;

  NEG_INF = -1e300;
  P_DONOR    = 0.02;
  P_ACCEPTOR = 0.05;

  EmitExon   : array[0..3] of Double = (0.25, 0.25, 0.25, 0.25);
  EmitIntron : array[0..3] of Double = (0.32, 0.32, 0.18, 0.18);
```
with:
```pascal
const
  TAG_HEADER    = 0;
  TAG_EXON      = 1;
  TAG_INTRON    = 2;
  TAG_NONCODING = 3;
  TAG_START     = 4;
  TAG_STOP      = 5;

  CLR_EXON_BG      = $00FFE8B0;
  CLR_INTRON_BG    = $00C8E8FF;
  CLR_HDR_BG       = $00E8E8E8;
  CLR_NONCODING_BG = $00D8D8D8;
  CLR_START_BG     = $00B0FFC8;
  CLR_STOP_BG      = $00B0B0FF;
```

In the `implementation` section, delete the standalone `function SafeLog` and `function NtIndex` (now from `HMMCore`), and delete the entire `procedure TForm1.AddSegment`, `procedure TForm1.RunHMM`, and `procedure TForm1.ComputeProbStats` bodies.

- [ ] **Step 3: Update `lbOutputDrawItem` to color the 3 new tags**

Replace:
```pascal
  case ATag of
    TAG_EXON  : BG := CLR_EXON_BG;
    TAG_INTRON: BG := CLR_INTRON_BG;
  else          BG := CLR_HDR_BG;
  end;
```
with:
```pascal
  case ATag of
    TAG_EXON     : BG := CLR_EXON_BG;
    TAG_INTRON   : BG := CLR_INTRON_BG;
    TAG_NONCODING: BG := CLR_NONCODING_BG;
    TAG_START    : BG := CLR_START_BG;
    TAG_STOP     : BG := CLR_STOP_BG;
  else             BG := CLR_HDR_BG;
  end;
```

- [ ] **Step 4: Update `btnRunClick` to call `HMMCore.RunHMM`**

Replace:
```pascal
  RunHMM(Seq);
  DisplayResults(Seq);
```
with:
```pascal
  HMMCore.RunHMM(Seq, FSegments, FSegCount);
  DisplayResults(Seq);
```

- [ ] **Step 5: Update `DisplayResults` for the 3 new segment kinds and Start/Stop-inclusive spliced mRNA**

Replace the whole `DisplayResults` procedure body with:
```pascal
procedure TForm1.DisplayResults(const Seq: string);
var
  i           : Integer;
  seg         : TSegment;
  SplicedExon : string;
  NoncodingNum, StartNum, ExonNum, IntronNum, StopNum : Integer;
  ExonLL, IntronLL, PriorCost, TotalLL, NullLL : Double;
begin
  lbOutput.Items.BeginUpdate;
  lbOutput.Items.Clear;
  FOutputTagCnt := 0;
  SetLength(FOutputTags, 128);

  AddOutputLine('HMM mRNA SPLICING RESULTS', TAG_HEADER);
  AddOutputLine('Sequence length : ' + IntToStr(Length(Seq)) + ' nt', TAG_HEADER);
  AddOutputLine('Start codon     : 5'' ATG', TAG_HEADER);
  AddOutputLine('Donor site      : 5'' GTA', TAG_HEADER);
  AddOutputLine('Acceptor site   : CAG 3''', TAG_HEADER);
  AddOutputLine('Stop codon      : TAA / TAG / TGA 3''', TAG_HEADER);
  AddOutputLine('', TAG_HEADER);

  NoncodingNum := 0;
  StartNum     := 0;
  ExonNum      := 0;
  IntronNum    := 0;
  StopNum      := 0;
  SplicedExon  := '';

  for i := 0 to FSegCount - 1 do
  begin
    seg := FSegments[i];
    case seg.Kind of
      skNoncoding:
      begin
        Inc(NoncodingNum);
        AddOutputLine(Format('NONCODING #%d   [pos %d..%d]   len=%d',
          [NoncodingNum, seg.StartIdx, seg.EndIdx, seg.EndIdx - seg.StartIdx + 1]),
          TAG_NONCODING);
        AddOutputLine('  ' + seg.Sequence, TAG_NONCODING);
      end;
      skStartCodon:
      begin
        Inc(StartNum);
        AddOutputLine(Format('START CODON #%d   [pos %d..%d]   len=%d',
          [StartNum, seg.StartIdx, seg.EndIdx, seg.EndIdx - seg.StartIdx + 1]),
          TAG_START);
        AddOutputLine('  ' + seg.Sequence, TAG_START);
        SplicedExon := SplicedExon + seg.Sequence;
      end;
      skExon:
      begin
        Inc(ExonNum);
        AddOutputLine(Format('EXON   #%d   [pos %d..%d]   len=%d',
          [ExonNum, seg.StartIdx, seg.EndIdx, seg.EndIdx - seg.StartIdx + 1]),
          TAG_EXON);
        AddOutputLine('  ' + seg.Sequence, TAG_EXON);
        SplicedExon := SplicedExon + seg.Sequence;
      end;
      skIntron:
      begin
        Inc(IntronNum);
        AddOutputLine(Format('INTRON #%d   [pos %d..%d]   len=%d spliced',
          [IntronNum, seg.StartIdx, seg.EndIdx, seg.EndIdx - seg.StartIdx + 1]),
          TAG_INTRON);
        AddOutputLine('  ' + seg.Sequence, TAG_INTRON);
      end;
      skStopCodon:
      begin
        Inc(StopNum);
        AddOutputLine(Format('STOP CODON #%d   [pos %d..%d]   len=%d',
          [StopNum, seg.StartIdx, seg.EndIdx, seg.EndIdx - seg.StartIdx + 1]),
          TAG_STOP);
        AddOutputLine('  ' + seg.Sequence, TAG_STOP);
        SplicedExon := SplicedExon + seg.Sequence;
      end;
    end;
  end;

  AddOutputLine('', TAG_HEADER);
  AddOutputLine(Format('Summary: %d noncoding, %d start codon(s), %d exon(s), %d intron(s), %d stop codon(s).',
                        [NoncodingNum, StartNum, ExonNum, IntronNum, StopNum]), TAG_HEADER);
  lbOutput.Items.EndUpdate;

  lbSpliced.Items.Clear;
  lbSpliced.Items.Add('SPLICED mRNA');
  lbSpliced.Items.Add('Length: ' + IntToStr(Length(SplicedExon)) + ' nt');
  lbSpliced.Items.Add(SplicedExon);

  ComputeProbStats(Seq, FSegments, FSegCount, ExonLL, IntronLL, PriorCost, TotalLL, NullLL);

  lbSpliced.Items.Add('');
  lbSpliced.Items.Add('HMM / Viterbi Path Probability Statistics');
  lbSpliced.Items.Add(Format('Exon emission log-likelihood     : %.4f', [ExonLL]));
  lbSpliced.Items.Add(Format('Intron emission log-likelihood   : %.4f', [IntronLL]));
  lbSpliced.Items.Add(Format('Path probability  exp(total LL)   : %.6e', [Exp(TotalLL)]));
end;
```

- [ ] **Step 6: Compile via lazbuild**

Run: `cd "/mnt/Data/yep/Kuliah/Tugas/Semester 6 bosen/Komgen Finals/Komgen ETS/Program 2 (B)" && /home/pingwalk/lazarus-4.2-0/lazarus/lazbuild project2.lpi`
Expected: `(1008) NNN lines compiled` with no `Fatal:`/`Error:` lines.

- [ ] **Step 7: Commit**

```bash
cd "/mnt/Data/yep/Kuliah/Tugas/Semester 6 bosen/Komgen Finals/Komgen ETS/Program 2 (B)"
git add unit1.pas
git commit -m "Wire GUI to 15-state HMMCore engine (Noncoding/Start/Stop codon display)"
```

## Out of scope (carried over from the design doc)

- Updating the static legend panel (`pnlLegend` labels in `unit1.lfm`) to show Noncoding/Start/Stop swatches — the output list already color-codes these via `lbOutputDrawItem`; the legend panel still only mentions Exon/Intron. Cosmetic, can be a follow-up.
- True reading-frame tracking across intron-spliced exon boundaries (per design decision).
- Reverse-strand scanning.
