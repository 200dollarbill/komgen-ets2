unit Unit1;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Forms, Controls, Graphics, Dialogs, StdCtrls, ExtCtrls,
  ComCtrls, LCLType, HMMCore;

type
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
    FSegments    : TSegmentArray;
    FSegCount    : Integer;
    FOutputTags  : array of Integer;
    FOutputTagCnt: Integer;
    function  CleanSequence(const Raw: string): string;
    procedure DisplayResults(const Seq: string);
    procedure AddOutputLine(const S: string; ATag: Integer);
  public
  end;

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

var
  Form1: TForm1;

implementation

{$R *.lfm}

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
    TAG_EXON     : BG := CLR_EXON_BG;
    TAG_INTRON   : BG := CLR_INTRON_BG;
    TAG_NONCODING: BG := CLR_NONCODING_BG;
    TAG_START    : BG := CLR_START_BG;
    TAG_STOP     : BG := CLR_STOP_BG;
  else             BG := CLR_HDR_BG;
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

  HMMCore.RunHMM(Seq, FSegments, FSegCount);
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
