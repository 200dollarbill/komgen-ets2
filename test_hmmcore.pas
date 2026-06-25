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
  // CC | ATG | AAA | GTACCCCCAG (intron, GTA..CAG) | CCC | TAA | GG
  // Note: the unchanged 8-state intron/exon core has a pre-existing
  // off-by-one that absorbs the first post-intron base into the intron
  // segment (verified against the original unit1.pas logic before this
  // change) -- preserved here on purpose since the core is unchanged.
  RunHMM('CCATGAAAGTACCCCCAGCCCTAAGG', Segs, Cnt);
  Check(Cnt = 7, 'segment count = 7 (got ' + IntToStr(Cnt) + ')');
  CheckSeg(Segs, 0, Cnt, skNoncoding,  'CC',             'leading UTR');
  CheckSeg(Segs, 1, Cnt, skStartCodon, 'ATG',            'start codon');
  CheckSeg(Segs, 2, Cnt, skExon,       'AAA',            'exon 1');
  CheckSeg(Segs, 3, Cnt, skIntron,     'GTACCCCCAGC',    'intron (absorbs 1 post-CAG base, pre-existing)');
  CheckSeg(Segs, 4, Cnt, skExon,       'CC',             'exon 2');
  CheckSeg(Segs, 5, Cnt, skStopCodon,  'TAA',            'stop codon');
  CheckSeg(Segs, 6, Cnt, skNoncoding,  'GG',             'trailing UTR');
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
