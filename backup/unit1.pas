unit Unit1;

{$mode objfpc}{$H+}

interface

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

  { TForm1 }

  TForm1 = class(TForm)
    btnClear: TButton;
    btnRun: TButton;
    lbInput: TListBox;
    lblInput: TLabel;
    lblOutput: TLabel;
    lblSpliced: TLabel;
    lbOutput: TListBox;
    lbSpliced: TListBox;
    pnlLegend: TPanel;
    pnlTop       : TPanel;
    lblTitle     : TLabel;
    lblLegend    : TLabel;
    lblExonLeg   : TLabel;
    lblIntronLeg : TLabel;
    procedure FormCreate(Sender: TObject);
    procedure btnRunClick(Sender: TObject);
    procedure btnClearClick(Sender: TObject);
    procedure lbOutputDrawItem(Control: TWinControl; Index: Integer;
                               ARect: TRect; State: TOwnerDrawState);
    procedure lbSplicedDrawItem(Control: TWinControl; Index: Integer;
                                ARect: TRect; State: TOwnerDrawState);
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
  public
  end;

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

var
  Form1: TForm1;

implementation

{$R *.lfm}

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

procedure TForm1.FormCreate(Sender: TObject);
const
  DEFAULT_SEQ =
    'C,C,G,G,C,A,C,T,G,T,T,C,A,T,G,G,G,C,A,A,T,G,C,A,A,G,G,T,A,C,G,G,T,G,A,G,C,' +
    'A,G,G,T,A,A,G,T,G,A,T,T,A,A,T,G,C,A,T,T,T,C,T,C,G,C,C,A,G,T,G,G,C,T,A,G,A,C,' +
    'G,A,T,G,C,A,T,A,G,G,A,G,A,T,C,A,T,T,G,A,C,G,A,T,G,C,A,T,A,G,A,C,C,G,G,A,A,G,' +
    'C';
begin
  lbInput.Items.Add(DEFAULT_SEQ);
end;

procedure TForm1.AddOutputLine(const S: string; ATag: Integer);
begin
  lbOutput.Items.Add(S);
  if FOutputTagCnt >= Length(FOutputTags) then
    SetLength(FOutputTags, FOutputTagCnt + 64);
  FOutputTags[FOutputTagCnt] := ATag;
  Inc(FOutputTagCnt);
end;

procedure TForm1.lbOutputDrawItem(Control: TWinControl; Index: Integer;
                                   ARect: TRect; State: TOwnerDrawState);
var
  LB   : TListBox;
  ATag : Integer;
  BG   : TColor;
begin
  LB := Control as TListBox;
  if (Index < 0) or (Index >= FOutputTagCnt) then ATag := TAG_HEADER
  else ATag := FOutputTags[Index];

  case ATag of
    TAG_EXON  : BG := CLR_EXON_BG;
    TAG_INTRON: BG := CLR_INTRON_BG;
  else          BG := CLR_HDR_BG;
  end;

  LB.Canvas.Brush.Color := BG;
  LB.Canvas.FillRect(ARect);
  LB.Canvas.Font.Color := clBlack;
  if odSelected in State then
  begin
    LB.Canvas.Brush.Color := clHighlight;
    LB.Canvas.Font.Color  := clHighlightText;
    LB.Canvas.FillRect(ARect);
  end;
  LB.Canvas.TextOut(ARect.Left + 4, ARect.Top + 2, LB.Items[Index]);
end;

procedure TForm1.lbSplicedDrawItem(Control: TWinControl; Index: Integer;
                                    ARect: TRect; State: TOwnerDrawState);
var
  LB : TListBox;
  BG : TColor;
begin
  LB := Control as TListBox;
  if Index = 2 then BG := CLR_EXON_BG
  else BG := CLR_HDR_BG;

  LB.Canvas.Brush.Color := BG;
  LB.Canvas.FillRect(ARect);
  LB.Canvas.Font.Color := clBlack;
  if odSelected in State then
  begin
    LB.Canvas.Brush.Color := clHighlight;
    LB.Canvas.Font.Color  := clHighlightText;
    LB.Canvas.FillRect(ARect);
  end;
  LB.Canvas.TextOut(ARect.Left + 4, ARect.Top + 2, LB.Items[Index]);
end;

function TForm1.CleanSequence(const Raw: string): string;
var
  i : Integer;
  c : Char;
begin
  Result := '';
  for i := 1 to Length(Raw) do
  begin
    c := UpCase(Raw[i]);
    if c in ['A','C','G','T'] then
      Result := Result + c;
  end;
end;

procedure TForm1.AddSegment(AStart, AEnd: Integer; AKind: TSegmentKind;
                             const ASeq: string);
begin
  if AEnd < AStart then Exit;
  if FSegCount >= Length(FSegments) then
    SetLength(FSegments, FSegCount + 32);
  FSegments[FSegCount].StartIdx := AStart;
  FSegments[FSegCount].EndIdx   := AEnd;
  FSegments[FSegCount].Kind     := AKind;
  FSegments[FSegCount].Sequence := ASeq;
  Inc(FSegCount);
end;

procedure TForm1.RunHMM(const Seq: string);
var
  State       : THMMState;
  i           : Integer;
  nt          : Char;
  ExonStart   : Integer;
  IntronStart : Integer;
  DonorBuf    : string;
  AccBuf      : string;
begin
  FSegCount := 0;
  SetLength(FSegments, 64);
  if Length(Seq) = 0 then Exit;

  State       := hsExon;
  ExonStart   := 1;
  IntronStart := 0;
  DonorBuf    := '';
  AccBuf      := '';

  i := 1;
  while i <= Length(Seq) do
  begin
    nt := Seq[i];
    case State of

      hsExon:
      begin
        if nt = 'G' then
        begin
          DonorBuf    := 'G';
          IntronStart := i;
          State       := hsIntronStart1;
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
          end;
        end;
      end;

      hsIntronStart3:
      begin
        if IntronStart - 1 >= ExonStart then
          AddSegment(ExonStart, IntronStart - 1, skExon,
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
        AddSegment(IntronStart, i, skIntron,
                   Copy(Seq, IntronStart, i - IntronStart + 1));
        ExonStart := i + 1;
        AccBuf    := '';
        State     := hsExon;
      end;

    end;
    Inc(i);
  end;

  case State of
    hsExon, hsIntronStart1, hsIntronStart2:
      if ExonStart <= Length(Seq) then
        AddSegment(ExonStart, Length(Seq), skExon,
                   Copy(Seq, ExonStart, Length(Seq) - ExonStart + 1));
    hsIntron, hsIntronEnd1, hsIntronEnd2:
      if IntronStart <= Length(Seq) then
        AddSegment(IntronStart, Length(Seq), skIntron,
                   Copy(Seq, IntronStart, Length(Seq) - IntronStart + 1));
    else ;
  end;
end;

procedure TForm1.ComputeProbStats(const Seq: string; out ExonLL, IntronLL,
  PriorCost, TotalLL, NullLL: Double);
var
  i, j, NtIdx : Integer;
  seg         : TSegment;
begin
  ExonLL    := 0.0;
  IntronLL  := 0.0;
  PriorCost := 0.0;

  for i := 0 to FSegCount - 1 do
  begin
    seg := FSegments[i];
    for j := 1 to Length(seg.Sequence) do
    begin
      NtIdx := NtIndex(seg.Sequence[j]);
      if NtIdx < 0 then Continue;
      if seg.Kind = skExon then
        ExonLL := ExonLL + SafeLog(EmitExon[NtIdx])
      else
        IntronLL := IntronLL + SafeLog(EmitIntron[NtIdx]);
    end;
    if seg.Kind = skIntron then
      PriorCost := PriorCost + SafeLog(P_DONOR) + SafeLog(P_ACCEPTOR);
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

procedure TForm1.DisplayResults(const Seq: string);
var
  i           : Integer;
  seg         : TSegment;
  SplicedExon : string;
  ExonNum, IntronNum : Integer;
  ExonLL, IntronLL, PriorCost, TotalLL, NullLL : Double;
begin
  lbOutput.Items.BeginUpdate;
  lbOutput.Items.Clear;
  FOutputTagCnt := 0;
  SetLength(FOutputTags, 128);

  AddOutputLine('HMM mRNA SPLICING RESULTS', TAG_HEADER);
  AddOutputLine('Sequence length : ' + IntToStr(Length(Seq)) + ' nt', TAG_HEADER);
  AddOutputLine('Donor site      : 5'' GTA', TAG_HEADER);
  AddOutputLine('Acceptor site   : CAG 3''', TAG_HEADER);
  AddOutputLine('', TAG_HEADER);

  ExonNum     := 0;
  IntronNum   := 0;
  SplicedExon := '';

  for i := 0 to FSegCount - 1 do
  begin
    seg := FSegments[i];
    case seg.Kind of
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
        AddOutputLine(Format('INTRON #%d   [pos %d..%d]   len=%d   << spliced out >>',
          [IntronNum, seg.StartIdx, seg.EndIdx, seg.EndIdx - seg.StartIdx + 1]),
          TAG_INTRON);
        AddOutputLine('  ' + seg.Sequence, TAG_INTRON);
      end;
    end;
  end;

  AddOutputLine('', TAG_HEADER);
  AddOutputLine(Format('Summary: %d exon(s),  %d intron(s) found.',
                        [ExonNum, IntronNum]), TAG_HEADER);
  lbOutput.Items.EndUpdate;

  lbSpliced.Items.Clear;
  lbSpliced.Items.Add('SPLICED mRNA  (exons concatenated)');
  lbSpliced.Items.Add('Length: ' + IntToStr(Length(SplicedExon)) + ' nt');
  lbSpliced.Items.Add(SplicedExon);

  ComputeProbStats(Seq, ExonLL, IntronLL, PriorCost, TotalLL, NullLL);

  lbSpliced.Items.Add('');
  lbSpliced.Items.Add('HMM / Viterbi Path Probability Statistics');
  lbSpliced.Items.Add(Format('Exon emission log-likelihood     : %.4f', [ExonLL]));
  lbSpliced.Items.Add(Format('Intron emission log-likelihood   : %.4f', [IntronLL]));
  //lbSpliced.Items.Add(Format('Splice-site prior cost (%d intron(s)) : %.4f',
  //                            [IntronNum, PriorCost]));
  //lbSpliced.Items.Add(Format('Total path log-likelihood         : %.4f', [TotalLL]));
  //lbSpliced.Items.Add(Format('Null model (no splicing) log-lik  : %.4f', [NullLL]));
  //lbSpliced.Items.Add(Format('Log-odds score (path vs. null)    : %.4f', [TotalLL - NullLL]));
  lbSpliced.Items.Add(Format('Path probability  exp(total LL)   : %.6e', [Exp(TotalLL)]));
end;

procedure TForm1.btnRunClick(Sender: TObject);
var
  Raw, Seq : string;
  i        : Integer;
begin
  Raw := '';
  for i := 0 to lbInput.Items.Count - 1 do
    Raw := Raw + lbInput.Items[i];
  Seq := CleanSequence(Raw);

  if Length(Seq) = 0 then
  begin
    ShowMessage('Please enter a valid nucleotide sequence (A/C/G/T).');
    Exit;
  end;

  RunHMM(Seq);
  DisplayResults(Seq);
end;

procedure TForm1.btnClearClick(Sender: TObject);
begin
  lbInput.Items.Clear;
  lbOutput.Items.Clear;
  lbSpliced.Items.Clear;
  FSegCount     := 0;
  FOutputTagCnt := 0;
end;

end.
