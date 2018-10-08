// "PALcolour": software PAL (maybe in future also NTSC) composite-video decode-to-colour image
//   Performs 2D subcarrier filtering to process stand-alone fields of video signal

// Copyright (C) 2018  William Andrew Steer


//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
//   Contact the author at palcolour@techmind.org



// Long history, goes back to integer math version to process off-air composite captures from
// the London (UK) Crystal Palace TV transmitters in the early 2000's.

// Then dormant for several years, then re-vamped to do "Colour Recovery" (restoring PAL colour from 1970's
// black-and-white film-recordings of PAL-coded shows, with all the distortions of film and CRTs, and no
// chroma-burst reference). This was quite 'experimental' and again not tidied to the point where other
// people could use the code - but served its purpose as a proof-of-principle and attracting more interest
// to that project.

// Then dormant for almost another decade, now using for Laser Disc and Domesday
// disc colour decoding (with proper burst-based PAL decoding).

// October 2018, quickly modified file-open routines to take raw ".tbc" (timebase corrected) files from
// the open-source ld-decode LaserDisc / Domesday project.
// The .tbc files are header-less files of 16bit greylevel data and an assumed frame size of 1052x610
// with the odd and even fields interleaved by previous processing steps (which is not how raw composite video is)
// The PAL tbc files are sampled at 4Fsc (4x PAL subcarrier frequency of 4.43361875 MHz, so 17.734475 MHz.)
// The .tbc files may have multiple frames following each other - you can only could how many by the filesize.
// At present (2018/10/08) the program attempts to open about 10 frames from the start of the file, but I hope
// to make it give the user some options...
// The .tbc files from ld-decode do not include the horizontal sync pulse, so the line-length is less than
// the 1135 pixels you'd expect from 64us (PAL line duration) at 17.73MHz.
// While the .tbc files are internally 16-bit, this version of the code takes just the top 8-bits to
// copy into a regular 8-bit image buffer and processes from there. It would not be difficult to make the
// PALcolour processing work from the original 16-bit digitised values if this was desired.

// Otherwise the program will open regular 8-bit greyscale Windows .bmp bitmaps.

// The program needs some guidance on where to look for the chroma burst (colour signal phase reference),
// back porch (black level reference), etc.
// Historically when I assumed the 64us lines were digitised, this was easy, given the PAL spec. With parts missing,
// it's not so straightforward. This version has some of these numbers hardcoded assuming 17.73MHz samples with the
// sync missing. I hope to tidy it to make it clearer, and add UI options to guide it in odd-cases again soon.

// Note that fundamentally PALdecode can cope with arbitrary sample rates, they do not need to be multiples of Fsc, though
// at sample rates much below 12-13MHz, aliassing will increasingly tend to cause image quality problems.

// Also added quick function to show "oscilloscope" trace
// of video line - something I'd been meaning to do for years.

// Still needs more documentation - and refactoring to tidy up!
// there's still a lot of "magic numbers" in the code, which really ought to be banished, and replaced
// with #defines or constants, and calculated values.
// also there's a number of points in the code you can make aesthetic tweaks (on filter bandwidth etc)
// which really ought to be parameterised and the tradeoffs explained.

// Note that this code has never done NTSC "properly" (there were a few nasty hacks) - though in principle
// it wouldn't be difficult to do an NTSC version.

// The "proper" PAL decoding is presently all 8-bit (8-bit CVBS, 8-bit RGB). The Colour Recovery version
// did work from professional/broadcast files so it worked with higher bit-depths. Not difficult to change.

// All the interesting/clever stuff happens in the function PALdecode1Click()
// everything else is just housekeeping, file opening/closing, UI, and utils.

// Bear in mind that much of this code was originally written in the year 2000, yet comments have
// been added more recently. With the passage of time, I've forgotten exactly how and why some of it works!

// For anyone wanting a better insight into what's going on, I recommend
//  - getting hold of the "Rec.470" PAL TV standard https://www.itu.int/rec/R-REC-BT.470/en
//  - being familiar with the principles of Fourier analysis and Fourier integrals
//  - and the RF principle of "product detection".

// It's also very informative to construct a zoneplate image (overlaying a 'chroma burst' in the
// correct place down the left hand edge, then processing that. This will visually reveal show
// what the filters are doing in 2D frequency-space.

// Its also worth bearing in mind that fundamentally the PAL (and NTSC) systems add a chroma subcarrier
// patterning to the image, but do not at production filter out genuine luma patterning ("pin striped shirts")
// at the same frequency. Thus in the general case, there is no "perfect" filter, and fine luminance patterning
// will cause "cross colour" effects (false/moving colour patches) and abrupt colour transitions will cause
// luminance patterning in the decoded images. Modern comb filtering techniques are vastly superior to the
// early analog "delay-line" decoders, but image-artifacts are still liable to pop up, and various aesthetic
// compromises are wrapped up within any such colour-decoding algorithms. As things stand, parameters have been
// tuned to what *I* consider to be the least objectionable artifacts for my purposes and aesthetics, across
// a handful of test images. YMMV.
// The algorithm embodied in this code is a non-adaptive 2D filter.

// However when PAL is coded, the vertical "bandwidth" is not limited, i.e. there is no requirement for adjacent
// scan lines to be colour-correlated. 2D filters gain their performance improvement by making assumptions that
// colour on adjacent lines often is correlated. When it isn't, especially at strong colour transitions, or
// with overlayed graphics or tickers, then artifacts such as "hanging dots" appear.

// Do a web search for "cross colour" "cross luma" "dot crawl" "hanging dots" for examples of the sort of
// artifacts to watch out for in PAL/NTSC decoding. These issues tend to be most noticeable on cartoons or
// computer generated graphics and graphic overlays. Natural images tend to be less problematic, though the
// coloured balls on a green snooker table make plenty of dot-crawl dots.

// More sophisticated adaptive filters can be developed, but have more complex failure-modes too.
// This algorithm/filter aims to be pretty much as good as you can get with a non-adaptive filter.
// That said, it may still want some tweaking depending on the quality (and timing stability) of the source material.
// Again, the sophisticated "tighter" filters may tend to only work well with high quality source material,
// and may be worse than a more basic filter with poor material.


// The present full project will compile in the long-obsolete Borland C++ Builder version 6 for Windows.
// However it'd be relatively straghtforward to port/fork the core function PALdecode1Click() to some other
// environment, then build up the bitmap-handling and GUI (or command-line UI) around it.


// Below are some historical notes in reverse chromological order

// Mods 12 October 2005

// FLOATING POINT version
// INCOMPLETE MOD needs tidying/checking/simplifying
//  Does not compensate for HF rolloff by reference to burst

//  Trying to remove the empirical brightness and saturation hacks, and
//  get back to a mathematically-accurate baseline, at least.
//  From Rec.470, levels for PAL-I are
//    sync      -43       0              [  0]
//    blanking    0      30              [ 62]
//    pk-white  100     100              [208]
//    pk s/c    133     123*             [255]

//    *  '123' is not quoted, I calculated it
//  in [], derived digitised nominal PAL-I CVBS with no spare head/footroom

//  https://tech.ebu.ch/docs/tech/tech3280.pdf  also has advice/standards on
//  levels for digitised composite video. It says:
//
//  For digitised composite (PAL), in 8-bit:
//    sync-tip should be at code 01h
//    blanking level should be at code 40h  (64 decimal)
//    peak-white should be at code D3h (211 decimal)

//  Both these references are pretty closely aligned, and will give
//  sufficient headroom for the chroma on, for example, the yellow bar of the
//  chroma bar, to go above peak-white, as it needs to.

//  Thus in my decodes to regular full-range RGB computer-land images, the
//   nominal 'contrast' gain needs to be x1.75 (256/(211-64))
//  If you were going to Rec.601 digital images, which have black at 16 and
//   white at 240, then you'd need slightly less gain, and to add the set-up.
//  Saturation should look after itself (if maths is all unity-gain)
//  ...unless there's HF-rolloff - in which case could scale relative to burst
//        but care required since actual peak-white level not known
//        could base it on burst amplitude compared to blanking level?

// Modified 23 November 2002
//  Increased line-buffer size to enable operation with 27MHz sampled signals

// Now with image in a scroll-box, a StatusBar, and a ProgressBar indicator

// 2D filtering, with separate U,V, and Y filters.
// Works extremely well!
//  FAST (3.3s) for 16MHz img with INTEGER MATHS... even on a PENTIUM II !!!


//---------------------------------------------------------------------------
#include <vcl\vcl.h>
#pragma hdrstop

#define MAX_WIDTH 1800

//Manual includes
#include <math.h>
#include <stdlib.h> // for atoi
//eoMI

#include "PALcolour_30.h"
#include "colctrls.h"
#include "unit3.h"

// palette structure for rebuilding 256 colour greyscale .bmps
typedef struct {
  TLogPalette lpal;
  TPaletteEntry entry[256];
} LogPal;



//---------------------------------------------------------------------------
#pragma resource "*.dfm"
TForm1 *Form1;
//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
	: TForm(Owner)
{
 image=NULL;
 pColImage=NULL;
 uyvybuf=NULL;

 //int digRate;
 //int Fsc;
 //char colsys;
 //long sat;
 //int brightness;

 digRate=17734476; // really it should pull this from the form at the start
 Fsc=4433619; colsys='P';
 sat=100; // colour saturation control: 0-100, 64typ.
 brightness=100; // 70ish for a perfect-level capture
 filtHarmonics=FALSE;

 // We will manually fill any excess space!
 ScrollBox1->Brush->Style=bsClear;

 ProgressBar1->Parent=StatusBar1;
 StatusBar1->Panels->Items[1]->Style=psOwnerDraw;
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Exit1Click(TObject *Sender)
{
 if (uyvybuf!=NULL) GlobalFree(uyvybuf);

 Close();
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Open1Click(TObject *Sender)
{
 OpenDialog1->FilterIndex=1; // set to default to ".bmp" files (using .tbc will cause a mess!)
 if (OpenDialog1->Execute())
 {
  if (image!=NULL) delete image;
  image = new Graphics::TBitmap();
  image->LoadFromFile(OpenDialog1->FileName.c_str());

  if (pColImage!=NULL) delete pColImage;
  pColImage = new Graphics::TBitmap();
  pColImage->Width=image->Width;
  pColImage->Height=image->Height;
  pColImage->PixelFormat=pf24bit;
  pColImage->HandleType=bmDIB;

  for (int y = 0; y < image->Height; y++)
  {
   BYTE *ptr = (unsigned char*)pColImage->ScanLine[y];
   BYTE *buffer = (unsigned char*)image->ScanLine[y];

   for (int x = 0; x < image->Width; x++)
   {
    ptr[x*3] = ptr[x*3+1] = ptr[x*3+2] = buffer[x];
   }
  }

//  HorzScrollBar->Range=pColImage->Width;

//  Repaint();

  // just to keep the colourisation happy - April 2008
  if (uyvybuf!=NULL) GlobalFree(uyvybuf);
//  uyvybuf=GlobalAlloc(GMEM_FIXED, 1920L*2*1080);
    uyvybuf=GlobalAlloc(GMEM_FIXED, image->Width*2L*image->Height);


 MenuSaveC->Enabled=FALSE;
Image1->Picture->Graphic=pColImage;
ScrollBox1->HorzScrollBar->Range=pColImage->Width;
ScrollBox1->VertScrollBar->Range=pColImage->Height;
ScrollBox1Resize(this);
 }
}
//---------------------------------------------------------------------------
void __fastcall TForm1::FormDestroy(TObject *Sender)
{
 if (image!=NULL) delete image;
// if (pColImage!=NULL) delete pColImage;
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Settings1Click(TObject *Sender)
{
 if (Form2->ShowModal() == mrOk)  // which it always will!
 {
  digRate=atoi(Form2->ComboBox1->Text.c_str());
  Fsc=atoi(Form2->ComboBox3->Text.c_str());
  if (Form2->ComboBox2->ItemIndex==0)
   colsys='P';
  else
   colsys='N';

  sat=Form2->SatScroll->Position;
  brightness=Form2->BrightScroll->Position;
  filtHarmonics=false; // Form2->CheckBox2->Checked;
 }
}
//---------------------------------------------------------------------------
void __fastcall TForm1::MenuSaveCClick(TObject *Sender)
{
 if (SaveDialog1->Execute())
 {
  pColImage->SaveToFile(SaveDialog1->FileName.c_str());
 }
}
//---------------------------------------------------------------------------

void __fastcall TForm1::ScrollBox1Resize(TObject *Sender)
{
 int bottomfill; int rightfill;

 if (pColImage==NULL)
 {
  bottomfill=ScrollBox1->ClientHeight;
  rightfill=ScrollBox1->ClientWidth;
 }
 else
 {
  bottomfill=ScrollBox1->ClientHeight-pColImage->Height;
  if (bottomfill<0) bottomfill=0;
  rightfill=ScrollBox1->ClientWidth-pColImage->Width;
  if (rightfill<0) rightfill=0;
 }
 Shape1->Height=bottomfill;
 Shape2->Width=rightfill;
}
//---------------------------------------------------------------------------

void __fastcall TForm1::StatusBar1DrawPanel(TStatusBar *StatusBar,
      TStatusPanel *Panel, const TRect &Rect)
{
 ProgressBar1->BoundsRect=Rect;
}
//---------------------------------------------------------------------------

void __fastcall TForm1::About1Click(TObject *Sender)
{
 Application->MessageBox("Disclaimers!\nPALcolour decodes the colour in PAL-coded composite video images.\nThe code is messy and non-optimised! Expl: filter-widths now (kind-of) scale with samplerate, so bandwidth should be constant irrespective of sample-rate!\n\nI am well-aware that the code is far from optimal efficiency-wise.\n\n(C)W.A. Steer / www.techmind.org", "About PALcolour - October 2018", IDOK);
}
//---------------------------------------------------------------------------


void __fastcall TForm1::SaveUYVY1Click(TObject *Sender)
{
 SaveUYVYframe(2); // note frames go from 0 to 856
}
//---------------------------------------------------------------------------

void TForm1::SaveUYVYframe(int framenum)
{
  BOOL ioresult;
  HANDLE theFile;

  char filename[60]="D:\\coloured2.uyvy";

  // Open the file
  theFile=CreateFile(filename, GENERIC_WRITE, FILE_SHARE_READ, NULL,
                     OPEN_ALWAYS, FILE_FLAG_SEQUENTIAL_SCAN, NULL);

  if (theFile!=INVALID_HANDLE_VALUE)
  {
   unsigned long bytesWritten;

   DWORD frameno=framenum;  // note frames go from 0 to 856
   const DWORD framesz = 1920L*2*1080;
   long highword=0;   // silly hack to get from 2GB to 4GB. See also SetFilePointerEx
   SetFilePointer(theFile, framesz*frameno, &highword, FILE_BEGIN);

   ioresult=WriteFile(theFile, uyvybuf, 1920L*2*1080, &bytesWritten, NULL);
  }
  else
  {
   Application->MessageBox("Failed to open file to write", "Information", IDOK);
  }

  CloseHandle(theFile);
}



//---------------------------------------------------------------------------

/*
void __fastcall TForm1::Batchprocess1Click(TObject *Sender)
{
 for (int frame=300; frame<857; frame++)
 {
  char msg[50];
  wsprintf(msg,"Opening fr:%i",frame);
  StatusBar1->Panels->Items[0]->Text=msg;

  OpenUYVYframe(frame);
  StatusBar1->Panels->Items[0]->Text="";

  Colourise1Click(NULL);

  Repaint();
  wsprintf(msg,"Saving fr:%i",frame);
  StatusBar1->Panels->Items[0]->Text=msg;
  SaveUYVYframe(frame);
  StatusBar1->Panels->Items[0]->Text="";
 }
}
*/

//---------------------------------------------------------------------------


void __fastcall TForm1::Opentbc16bit1052x1Click(TObject *Sender)
{
 OpenDialog1->FilterIndex=2; // set to default to ".tbc" files (using bmp will cause a mess!)
 if (OpenDialog1->Execute())
 {
  if (image!=NULL) delete image;
  int framestoopen=50; // 15
  image = new Graphics::TBitmap();
  image->Width = 1052;
  image->Height = 610*framestoopen;
  image->PixelFormat=pf8bit;
  image->HandleType=bmDIB;

  // Palette malarky
  LogPal mypal;
  mypal.lpal.palVersion = 0x300;
  mypal.lpal.palNumEntries = 256;

  for (int pe=0; pe<mypal.lpal.palNumEntries; pe++)
  {
   mypal.entry[pe].peRed = pe;
   mypal.entry[pe].peGreen = pe;
   mypal.entry[pe].peBlue = pe;
   mypal.entry[pe].peFlags = NULL;
  }

  image->Palette = CreatePalette((const tagLOGPALETTE *)&mypal.lpal);
  // end palette malarky


  int frameno=0;


  BOOL ioresult;
  HANDLE theFile;

  // Open the file
  theFile=CreateFile(OpenDialog1->FileName.c_str(), GENERIC_READ, FILE_SHARE_READ, NULL,
                     OPEN_EXISTING, FILE_FLAG_SEQUENTIAL_SCAN, NULL);

  if (uyvybuf!=NULL) GlobalFree(uyvybuf);
  uyvybuf=GlobalAlloc(GMEM_FIXED, 1052L*610*2*framestoopen);

  if (theFile!=INVALID_HANDLE_VALUE)
  {
   unsigned long bytesRead;

//   DWORD frameno=0; // framenum; // 751;//750;  // note frames go from 0 to 856
   const DWORD framesz = 1052L*610*2;
   long highword=0;   // silly hack to get from 2GB to 4GB. See also SetFilePointerEx
   SetFilePointer(theFile, framesz*frameno, &highword, FILE_BEGIN);

   ioresult=ReadFile(theFile, uyvybuf, 1052L*610*2*framestoopen, &bytesRead, NULL);
  }
  else
  {
   Application->MessageBox("Failed to open file", "Information", IDOK);
  }

  CloseHandle(theFile);



  for (int y = 0; y < 610*framestoopen; y++)
  {
   BYTE *ptr = (unsigned char*)(uyvybuf)+1052L*2*y;
   BYTE *buffer = (unsigned char*)image->ScanLine[y];

   for (int x = 0; x < image->Width; x++)
   {
    buffer[x]= ptr[x*2+1]; // be lazy and just take the high-byte and work in 8-bit for now
   }
  }

  if (pColImage!=NULL) delete pColImage;
  pColImage = new Graphics::TBitmap();
  pColImage->Width=image->Width;
  pColImage->Height=image->Height;
  pColImage->PixelFormat=pf24bit;
  pColImage->HandleType=bmDIB;

  // We copy the composite (black and white) source into the coloured-image (display) buffer purely
  // so we can display it on the screen easily
  for (int y = 0; y < image->Height; y++)
  {
   BYTE *ptr = (unsigned char*)pColImage->ScanLine[y];
//   BYTE *ptruyvy = (unsigned char*)(uyvybuf) +1052L*2*y;
   BYTE *buffer = (unsigned char*)image->ScanLine[y];

   for (int x = 0; x < image->Width; x++)
   {
    ptr[x*3] = ptr[x*3+1] = ptr[x*3+2] = buffer[x];
   }
  }

  // just to keep the colourisation happy - April 2008
  if (uyvybuf!=NULL) GlobalFree(uyvybuf);
//  uyvybuf=GlobalAlloc(GMEM_FIXED, 1920L*2*1080);
    uyvybuf=GlobalAlloc(GMEM_FIXED, image->Width*2L*image->Height);


 MenuSaveC->Enabled=TRUE; //FALSE;
Image1->Picture->Graphic=pColImage;
ScrollBox1->HorzScrollBar->Range=pColImage->Width;
ScrollBox1->VertScrollBar->Range=pColImage->Height;
ScrollBox1Resize(this);
 }

}
//---------------------------------------------------------------------------

void __fastcall TForm1::PALdecode1Click(TObject *Sender)
{
 // Fsc=4433619;    //   4.43MHz

 float bright=1.75 * brightness/100.0;
 // NB 1.75 is nominal scaling factor for full-range digitised composite (with sync at code 0 or 1,
 // blanking at code 64 (40h), and peak white at code 211 (d3h) to give 0-255 RGB.


 if (image!=NULL)
 {
  int w=image->Width;

  // Step 1: create sine/cosine lookups
  float sine[MAX_WIDTH], cosine[MAX_WIDTH];    // formerly short int
  short int refAmpl=128;    // original scaling for integer math

  float rad;
  for (int i=0; i<w; i++)
  {
   rad=(2*M_PI*i*Fsc/digRate);
   sine[i]=refAmpl*sin(rad); cosine[i]=refAmpl*cos(rad);
  }

  // Next create filter-profiles for colour filtering.
  //  One can argue over merits of different filters, but I stick with simple raised cosine
  //  unless there's compelling reason to do otherwise.
  // PAL-I colour bandwidth should be around 1.1 or 1.2 MHz
  //  acc to Rec.470, +1066 or -1300kHz span of colour sidebands!

  // width of filter-window should therefore scale with samplerate

  // Create filter-profile lookup
  // chromaBandwidthHz values between 1.1MHz and 1.3MHz can be tried. Some specific values in that range may work best at minimising residual
  // dot pattern at given sample rates due to the discrete nature of the filters. It'd be good to find ways to optimise this more rigourously
  float chromaBandwidthHz=1100000.0 /0.93; // the 0.93 is a bit empirical for the 4Fsc sampled LaserDisc scans
  float ca=0.5*digRate/chromaBandwidthHz, ya=0.5*digRate/chromaBandwidthHz; // where does the 0.5* come from?
  // note in principle you could have different bandwidths for extracting the luma and chroma, according to aesthetic tradeoffs. Not really very justifyable though.

  const int a=18; // 'a' is the array-size, corresponding to at least half the filter-width, and should be at least Fsampling(max supported by build)/colourfilterBandwidth(min supported by build)
                  //  'a' must be greater than or equal to the bigger of 'ca' and 'ya' above
  float cfilt[4][a+1]; float yfilt[4][a+1];
  float cdiv=0; float ydiv=0;

  // Note that we choose to make the y-filter *much* less selective in the vertical direction:
  // - this is to prevent castellation on horizontal colour boundaries.

  // may wish to broaden vertical bandwidth *slightly* so as to better pass
  // one- or two-line colour bars - underlines/graphics etc.

  // Note also that if Y-bandwidth was made the same as C,
  // and that 'lines' of the masks were equivalent, then
  // significant time-savings could be made.

  for (short int f=0; f<=a; f++)
  {
   float  fc=f; if (fc>ca) fc=ca;
   float  ff=sqrt(f*f+2*2); if ( ff>ca)  ff=ca;  // 2 -- 4 -- 6 sequence
   float fff=sqrt(f*f+4*4); if (fff>ca) fff=ca;  // because only one FIELD!
   float ffff=sqrt(f*f+6*6); if (ffff>ca) ffff=ca;

   int d;
   if (f==0) d=2; else d=1; // divider because we're only making half a filter-kernel and the zero-th point is counted twice later.

   cfilt[0][f]=256*(1+cos(M_PI*fc/ca))/d;
   cfilt[1][f]=256*(1+cos(M_PI*ff/ca))/d;
   cfilt[2][f]=256*(1+cos(M_PI*fff/ca))/d;
   cfilt[3][f]=256*(1+cos(M_PI*ffff/ca))/d;

   cdiv+=cfilt[0][f]+2*cfilt[1][f]+2*cfilt[2][f]+2*cfilt[3][f];

   float  fy=f; if (fy>ya) fy=ya;
   float fffy=sqrt(f*f+4*4); if (fffy>ya) fffy=ya;

   yfilt[0][f]=256*(1+cos(M_PI*fy/ya))/d;
   yfilt[1][f]=0;
   yfilt[2][f]=0.2*256*(1+cos(M_PI*fffy/ya))/d;  // 0.2 makes much less sensitive to adjacent lines and reduces castellations and residual dot patterning
   yfilt[3][f]=0;

   ydiv+=yfilt[0][f]+2*yfilt[2][f];
  }
  cdiv*=2; ydiv*=2;


  // Step 2:
  BYTE *buffer, *b1, *b2, *b3, *b4, *b5, *b6, Y[MAX_WIDTH];
// were all short ints
  float pu[MAX_WIDTH], qu[MAX_WIDTH], pv[MAX_WIDTH], qv[MAX_WIDTH], py[MAX_WIDTH], qy[MAX_WIDTH];
  float m[MAX_WIDTH], n[MAX_WIDTH];
  float m1[MAX_WIDTH], n1[MAX_WIDTH], m2[MAX_WIDTH], n2[MAX_WIDTH];
  float m3[MAX_WIDTH], n3[MAX_WIDTH], m4[MAX_WIDTH], n4[MAX_WIDTH];
  float m5[MAX_WIDTH], n5[MAX_WIDTH], m6[MAX_WIDTH], n6[MAX_WIDTH];


  // lines below try to set X-positions of start of burst, end of burst, start of back-porch (blanking reference) etc
  // they SHOULD be derived from sample-rate, Rec.470 timings, and where X=0 corresponds to in the composite signal...
  // in practice these tend to get empirically set for given source images
  int Pstart=w*5.2/64;       // was 5.2
  //int Ystart=w*10.0/64; // blanking officially ends 10.5us after sync-start
  int Ystart=96; // for tbc images

  // Btstart/end:  8-20 (14MHz), 85-115 (15MHz), 90-120 (16MHz), 110-140 (20MHz)
  //int Btstart=Pstart; int Btend=(Pstart+Ystart)/2;

  //int Btstart=w*5.6/64; int Btend=w*(5.6+2.25)/64;
  int Btstart=w*1.2/64; int Btend=w*(1.2+2.25)/64; // burst start and burst end

  //  int Bkstart=Pstart, Bkend=Ystart; // rem 50-60 (14MHz), 120-150 (15MHz), 130-160 (16MHz), 160-185 (20MHz)
  //int Bkstart=w*(5.6+2.25 +0.4)/64; int Bkend=w*(10.5 -0.4)/64;
  int Bkstart=64; int Bkend=96; // empirical


  short int H=image->Height;

  int interleaved=1; // use 0 if raw CVBS (field-sequential), use 1 if signal has line the field-lines interleaved to make a single image
  if (Form2->Interleaved->Checked) interleaved=1; else interleaved=0;

  int Vsw; // this will represent the PAL Vswitch state later on...

  StatusBar1->Panels->Items[0]->Text="Colourising...";

  for (short int fi=0; fi<(interleaved+1); fi++)
  for (short int l=3*(interleaved+1)+fi; l<(H-3*(interleaved+1)); l+=(interleaved+1))
  {
   // show progress
   if (l%8) ProgressBar1->Position=(100*(l-3))/(H-6);  // not interlaced-mode aware

   buffer = (BYTE*)image->ScanLine[l];
   b1 = (BYTE*)image->ScanLine[l-1*(interleaved+1)];
   b2 = (BYTE*)image->ScanLine[l+1*(interleaved+1)];
   b3 = (BYTE*)image->ScanLine[l-2*(interleaved+1)];
   b4 = (BYTE*)image->ScanLine[l+2*(interleaved+1)];
   b5 = (BYTE*)image->ScanLine[l-3*(interleaved+1)];
   b6 = (BYTE*)image->ScanLine[l+3*(interleaved+1)];

   for (short int i=0; i<w; i++)
   {
    m[i]=buffer[i]*sine[i]; n[i]=buffer[i]*cosine[i];

    m1[i]=b1[i]*sine[i];  n1[i]=b1[i]*cosine[i];
    m2[i]=b2[i]*sine[i];  n2[i]=b2[i]*cosine[i];
    m3[i]=b3[i]*sine[i];  n3[i]=b3[i]*cosine[i];
    m4[i]=b4[i]*sine[i];  n4[i]=b4[i]*cosine[i];
    m5[i]=b5[i]*sine[i];  n5[i]=b5[i]*cosine[i];
    m6[i]=b6[i]*sine[i];  n6[i]=b6[i]*cosine[i];
   }


   // Find absolute burst phase

   // Does "some maths" on the burst to determine the Vswitch phase
   int bp=0, bq=0, bpo=0, bqo=0;
   for (short int i=Btstart; i<Btend; i++) { bp+=(m[i]-(m3[i]+m4[i])/2)/2; bq+=(n[i]-(n3[i]+n4[i])/2)/2; bpo+=(m2[i]-m1[i])/2; bqo+=(n2[i]-n1[i])/2; }

   bp/=(Btend-Btstart);  bq/=(Btend-Btstart);  // normalises those sums
   bpo/=(Btend-Btstart); bqo/=(Btend-Btstart); // normalises those sums

   // Generate V-switch phase
   if (((bp-bpo)*(bp-bpo)+(bq-bqo)*(bq-bqo))<(bp*bp+bq*bq)*2) Vsw=1; else Vsw=-1;
   if (colsys=='N') Vsw=1; // NTSC fixup!

   // NB bp and bq will be of the order of 1000.
   bp=(bp-bqo)/2; bq=(bq+bpo)/2; // ave two lines to get -U phase out
/*
   // Rotate burst phase according to V-switch
   int tbp=(bp*0.707-bq*Vsw*0.707);
   bq=bq*0.707+bp*Vsw*0.707;
   bp=tbp;
*/

   //int norm=sqrt(bp*bp+bq*bq)*refAmpl*16; // 16 empirical scaling factor
   float norm=sqrt(bp*bp+bq*bq); // TRIAL - 7 Oct 2005

   if (norm<130000) norm=130000;   // kill colour if burst too weak!

   // p & q should be sine/cosine components' amplitudes
   // NB: Multiline averaging/filtering assumes perfect
   //     inter-line phase registration...

   int PU,QU, PV,QV, PY,QY;
   for (short int i=Ystart; i<w-a; i++)
   {
    PU=QU=0; PV=QV=0; PY=QY=0;

    // Carry out 2D filtering. P and Q are the two arbitrary SINE & COS
    // phases components. U filters for U, V for V, and Y for Y
    // U and V are the same for lines n, n+/-2, but differ in sine for
    // n+/-1, n+/-3 owing to the forward/backward axis slant  /   \
    // For Y, only use lines n, n+/-2: the others cancel!!!
    //  *have tried* using lines +/-1 & 3 --- can be made to work, but
    //  introduces *phase-sensitivity* to the filter -> leaks too much
    //  subcarrier if *any* phase-shifts!

    register short int l,r;
    for (short int b=0; b<=a; b++)
    {
     l=i-b; r=i+b;

     PU+=(m[r]+m[l])*cfilt[0][b]+(+n1[r]+n1[l]-n2[l]-n2[r])*cfilt[1][b]-(m3[l]+m3[r]+m4[l]+m4[r])*cfilt[2][b]+(-n5[r]-n5[l]+n6[l]+n6[r])*cfilt[3][b];
     QU+=(n[r]+n[l])*cfilt[0][b]+(-m1[r]-m1[l]+m2[l]+m2[r])*cfilt[1][b]-(n3[l]+n3[r]+n4[l]+n4[r])*cfilt[2][b]+(+m5[r]+m5[l]-m6[l]-m6[r])*cfilt[3][b];
     PV+=(m[r]+m[l])*cfilt[0][b]+(-n1[r]-n1[l]+n2[l]+n2[r])*cfilt[1][b]-(m3[l]+m3[r]+m4[l]+m4[r])*cfilt[2][b]+(+n5[r]+n5[l]-n6[l]-n6[r])*cfilt[3][b];
     QV+=(n[r]+n[l])*cfilt[0][b]+(+m1[r]+m1[l]-m2[l]-m2[r])*cfilt[1][b]-(n3[l]+n3[r]+n4[l]+n4[r])*cfilt[2][b]+(-m5[r]-m5[l]+m6[l]+m6[r])*cfilt[3][b];

     PY+=(m[r]+m[l])*yfilt[0][b]-(m3[l]+m3[r]+m4[l]+m4[r])*yfilt[2][b];
     QY+=(n[r]+n[l])*yfilt[0][b]-(n3[l]+n3[r]+n4[l]+n4[r])*yfilt[2][b];
    }
    pu[i]=PU/cdiv; qu[i]=QU/cdiv;
    pv[i]=PV/cdiv; qv[i]=QV/cdiv;
    if (colsys=='N') { pv[i]=PU/cdiv; qv[i]=QU/cdiv; } // NTSC fixup
    py[i]=PY/ydiv; qy[i]=QY/ydiv;
   }

   // Obtain the black level from the "back porch"
   // Bkstart and Bkend define the zone of the back-porch, used for blacklevel reference
   int blacklevel=0;
   for (short int i=Bkstart; i<Bkend; i++) blacklevel+=buffer[i]+b1[i]+b2[i]+b3[i]+b4[i];
    if (colsys=='N') blacklevel*=(47.5/40); // NTSC fixup!
   blacklevel/=(Bkend-Bkstart)*5;

   int normalise=refAmpl*refAmpl/2;     // refAmpl is the integer sinewave amplitude

   // Generate the luminance (Y), by filtering out Fsc
   for (short int i=Ystart; i<w; i++)
   {
    short int tmp=buffer[i]-(py[i]*sine[i]+qy[i]*cosine[i])/normalise -blacklevel;
    if (tmp<0) tmp=0; if (tmp>255) tmp=255;
    Y[i]=tmp;
   }

   BYTE *ptr = (BYTE*)pColImage->ScanLine[l];

   // 'sat' is a user saturation control from scrollbar, nom. 100%
   float colscale=(sat/100.0)/(norm/2.0);  // 'norm' normalises bp and bq to 1
                                           // the '2' is arbitrary until I can work out a physical reason for it!
   for (short int i=Ystart; i<w-a; i++)
   {
    short int R, G, B;
    float U, V;

    U=-((pu[i]*bp+qu[i]*bq)) *colscale;
    V=-(Vsw*(qv[i]*bp-pv[i]*bq)) *colscale;

    // These magic numbers below come from the PAL matrices (I ought to have a reference for these. Tancock and/or Rec.470, I expect)
    R=bright*(Y[i]+1.14*V);            if (R<0) R=0; if (R>255) R=255;
    G=bright*(Y[i]-0.581*V-0.394*U);   if (G<0) G=0; if (G>255) G=255;  // coefficients corrected 10 Sept 2004
    B=bright*(Y[i]+2.03*U);            if (B<0) B=0; if (B>255) B=255;

    int pp=i*3;
    ptr[pp]=(BYTE)B; ptr[pp+1]=(BYTE)G; ptr[pp+2]=(BYTE)R;
   }
   ptr[Bkstart*3+1]=255; // show where black-level taken from
   ptr[Bkend*3+1]=255;  //
   ptr[Btstart*3+2]=255; // show where colour burst taken from
   ptr[Btend*3+2]=255;  //
  }
 }

 Image1->Picture->Graphic=pColImage;   // (pColImage has changed!)
 StatusBar1->Panels->Items[0]->Text=""; // restore std caption
 ProgressBar1->Position=0;
 MenuSaveC->Enabled=TRUE;
 ScrollBox1->Repaint();
}


void __fastcall TForm1::Image1MouseMove(TObject *Sender, TShiftState Shift,
      int X, int Y)
{
 char tmpstring[64];
 BYTE *ptr = (BYTE*)image->ScanLine[Y]; // in the original image
 BYTE *ptrcol = (BYTE*)pColImage->ScanLine[Y]; // in the original image

 wsprintf(tmpstring,"(%4d,%4d) C:%3d, R:%3d G:%3d B:%3d",X,Y, int(ptr[X]), int(ptrcol[X*3+2]), int(ptrcol[X*3+1]), int(ptrcol[X*3]));
 StatusBar1->Panels->Items[0]->Text=tmpstring;
}


//---------------------------------------------------------------------------

void __fastcall TForm1::Image1MouseDown(TObject *Sender,
      TMouseButton Button, TShiftState Shift, int X, int Y)
{
 BYTE *ptr = (BYTE*)image->ScanLine[Y]; // in the original image
 BYTE *ptrcol = (BYTE*)pColImage->ScanLine[Y]; // in the coloured image

 // erase the canvas - draw black rectangle
 Form3->Image1->Canvas->Brush->Color=0x000000;
 Form3->Image1->Canvas->FillRect(Rect(0,0,1052,256));

 // select a grey pen and draw markers at 64 (nominal blanking (black) level) and at 211 (nominal peak white level)
 Form3->Image1->Canvas->Pen->Color=0x808080;
 Form3->Image1->Canvas->MoveTo(0,255-64);
 Form3->Image1->Canvas->LineTo(1052,255-64);
 Form3->Image1->Canvas->MoveTo(0,255-211);
 Form3->Image1->Canvas->LineTo(1052,255-211);

 // select a white pen
 Form3->Image1->Canvas->Pen->Color=0xffffff;

 // and draw the 'oscilloscope' view of the composite video scan-line
 for (int x=0; x<((image->Width)-1); x++)
 {
  Form3->Image1->Canvas->MoveTo(x,255-int(ptr[x]));
  Form3->Image1->Canvas->LineTo(x+1,255-int(ptr[x+1]));
 }

 /*
 Image1->Canvas->Pen->Color=0x0000ff;

 for (int x=0; x<((image->Width)-1); x++)
 {
  Image1->Canvas->MoveTo(x,512-int(ptrcol[x*3+2]));
  Image1->Canvas->LineTo(x+1,512-int(ptrcol[(x+1)*3+2]));
 }

 Image1->Canvas->Pen->Color=0x00ff00;

 for (int x=0; x<((image->Width)-1); x++)
 {
  Image1->Canvas->MoveTo(x,512-int(ptrcol[x*3+1]));
  Image1->Canvas->LineTo(x+1,512-int(ptrcol[(x+1)*3+1]));
 }

 Image1->Canvas->Pen->Color=0xff0000;

 for (int x=0; x<((image->Width)-1); x++)
 {
  Image1->Canvas->MoveTo(x,512-int(ptrcol[x*3]));
  Image1->Canvas->LineTo(x+1,512-int(ptrcol[(x+1)*3]));
 }
 */
}
//---------------------------------------------------------------------------

