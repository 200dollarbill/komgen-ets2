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
  P_START    = 1.0 / 64.0;  // chance of 3 ATG
  P_STOP     = 3.0 / 64.0;  // chance of 3 TAA/TAG/TGA

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
        AddSeg(Segments, SegCount, StopStart, StopStart + 2, skStopCodon,
               Copy(Seq, StopStart, 3));
        NoncodingStart := StopStart + 3;
        StopBuf        := '';
        State          := hsNoncoding;
        Continue;
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
        AddSeg(Segments, SegCount, StopStart, StopStart + 2, skStopCodon,
               Copy(Seq, StopStart, 3));
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
