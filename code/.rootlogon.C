{
//cout << endl << "Welcome to the Alida rootlogon.C" << endl;
//

TStyle *style= new TStyle("style","Alida style");

// use plain black on white colors
Int_t icol=0; // WHITE
//style->SetFrameBorderMode(icol);
style->SetFrameFillColor(icol);
style->SetCanvasBorderMode(icol);
style->SetCanvasColor(icol);
style->SetPadBorderMode(icol);
//style->SetPadBorderMode(1);
style->SetPadColor(icol);
style->SetStatColor(icol);
style->SetPaperSize(20,26);
style->SetPadTopMargin(0.05);
style->SetPadRightMargin(0.10);
style->SetPadBottomMargin(0.16);
style->SetPadLeftMargin(0.16);
style->SetTitleXOffset(1.4);
style->SetTitleYOffset(1.4);
style->SetLegendBorderSize(0);

Int_t font=42; // Helvetica
Double_t tsize=0.04;
style->SetTextFont(font);
style->SetTextSize(tsize);
style->SetLabelFont(font,"x");
style->SetTitleFont(font,"x");
style->SetLabelFont(font,"y");
style->SetTitleFont(font,"y");
style->SetLabelFont(font,"z");
style->SetTitleFont(font,"z");
style->SetLabelSize(tsize,"x");
style->SetTitleSize(tsize,"x");
style->SetLabelSize(tsize,"y");
style->SetTitleSize(tsize,"y");
style->SetLabelSize(tsize,"z");
style->SetTitleSize(tsize,"z");
style->SetMarkerStyle(20);
style->SetMarkerSize(1.2);
style->SetHistLineWidth(2.);
style->SetLineStyleString(2,"[12 12]"); // postscript dashes
style->SetEndErrorSize(0.);
style->SetOptTitle(0);
style->SetOptStat(0);
style->SetOptFit(0);
style->SetPalette(75);
//style->SetBorderSize(0)
//style->SetPadTickX(1);
//style->SetPadTickY(1);
gROOT->SetStyle("style");
gROOT->ForceStyle();
}
