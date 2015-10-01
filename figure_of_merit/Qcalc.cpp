#include <iostream>
#include <cmath>
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TAttMarker.h"

using namespace std;


int main()
{
TChain *tree = new TChain("albero");
tree->Add("../Samples_BDT/BDT_NOCUT_DiJetPt_BJetPlusX_Run2012B-13Jul.root");
tree->Add("../Samples_BDT/BDT_NOCUT_DiJetPt_BJetPlusX_Run2012C-24Aug.root");
tree->Add("../Samples_BDT/BDT_NOCUT_DiJetPt_BJetPlusX_Run2012C-Promp.root");
tree->Add("../Samples_BDT/BDT_NOCUT_DiJetPt_BJetPlusX_Run2012D-Promp.root");

float data_BDT;
tree->SetBranchAddress("BDT", &data_BDT);
unsigned int nEvents = tree->GetEntries();

TH1D *BDT_B = new TH1D("bkg_bdt", "bkg_bdt", 100, 0, 1);
TH1D *BDT_S = new TH1D("mc_bdt", "mc_bdt", 100, 0 , 1);

for(unsigned int i=0; i<nEvents; ++i)
{
	tree->GetEntry(i);
	BDT_B->Fill((double)data_BDT);
}


TFile *mc_input = new TFile("../Samples_BDT/BDT_NOCUT_hh4b-nores-8Tev300k_step2_L1y1.root", "READ");
TTree* tree2 = (TTree*)mc_input->Get("albero");
unsigned int nnEvents = tree2->GetEntries();
float mc_BDT;
tree2->SetBranchAddress("BDT", &mc_BDT);
for(unsigned int i=0; i<nnEvents; ++i)
{
	tree2->GetEntry(i);
	BDT_S->Fill((double)mc_BDT);
}

TH1D *merit_BDT = new TH1D("Q_BDT", "BDT Q", 101, 0, 1);
TH1D *merit_SL = new TH1D("SL_BDT", "BDT SL", 101, 0, 1);
//	double integral = BDT_B->Integral();
	//bkg_bdt->Scale(2.83222748815166);
	BDT_S->Scale(2.01461E-4);		//branchin ratio updated
	double scale = 2.01461E-4;
	double scale_err = 0.623183E-4; 

	double cut = 0.0;
	double QBDT, s, b, SS2, SB2, QBDT_err, Lsb, Lb, SL, SL_err;

	while(cut <= 1)
	{
		s = BDT_S->Integral(BDT_S->FindBin(cut), BDT_S->FindBin(1.0));
		b = BDT_B->Integral(BDT_B->FindBin(cut), BDT_B->FindBin(1.0));
		SS2 = s*s*scale_err*scale_err + scale*scale*s;
		SB2 = b;

		if(s == 0 || b == 0) QBDT_err = 0.0;
		else{QBDT_err = 10*sqrt(4*pow(-1/(2.*sqrt(b)) + 1/(2.*sqrt(b + s)),2)*SB2 + SS2/(b + s));}
		
		if(s != 0 && b != 0) QBDT = 2*(sqrt(s+b)-sqrt(b));	
		else QBDT = 0;
		merit_BDT->SetBinContent(merit_BDT->FindBin(cut), QBDT);
		merit_BDT->SetBinError(merit_BDT->FindBin(cut), QBDT_err);


		if(s != 0 && b != 0){ 
			cout << s << "\t" << b << "\t" << s+b << "\t" << SL << endl;
		//	if((-s+(b+s)*log((b+s)/b))<0) cout << "PORCADDIO" <<endl;
			SL = sqrt(2.0)*sqrt(-s + (b+s)*log(1.0+(s/b)));
			SL_err = 10*sqrt((SS2*pow(log((b + s)/b),2))/(2.*(-s + (b + s)*log((b + s)/b))) + (SB2*pow(b*(1/b - (b + s)/pow(b,2)) + log((b + s)/b),2))/(2.*(-s + (b + s)*log((b + s)/b))));
			}
		else {SL = 0; SL_err = 0;}

		merit_SL->SetBinContent(merit_SL->FindBin(cut), SL);
		merit_SL->SetBinError(merit_SL->FindBin(cut), SL_err);
		cut+=0.01;

	}
			//merit_SL->Print("all");
	gStyle->SetOptStat(0);

	TCanvas *c = new TCanvas("c","c", 1000,400);
	c->Divide(2,1);
	merit_SL->SetLineWidth(1);
	merit_SL->SetMarkerStyle(8);
	merit_SL->SetMarkerColor(1);
	merit_SL->SetMarkerSize(0.7);
	merit_SL->SetLineColor(1);
	merit_SL->SetMinimum(0.0);
	merit_SL->SetTitle("Likelihood FOM");
	merit_SL->SetXTitle("BDT response cut");
	merit_SL->SetYTitle("Q = #sqrt{2 ln(L_{s+b}/L_{b})}");
	merit_SL->GetYaxis()->SetTitleOffset(1.4);


	merit_BDT->SetLineWidth(1);
	merit_BDT->SetMarkerStyle(8);
	merit_BDT->SetMarkerColor(1);
	merit_BDT->SetMarkerSize(0.7);
	merit_BDT->SetLineColor(1);
	merit_BDT->SetMinimum(0.0);
	merit_BDT->SetTitle("Bityukov-Krasnikov FOM");
	merit_BDT->SetXTitle("BDT response cut");
	merit_BDT->SetYTitle("Q = 2(#sqrt{s+b} - #sqrt{b})");
	merit_BDT->GetYaxis()->SetTitleOffset(1.4);
	
	c->cd(1);
	TLatex * mylatex = new TLatex(0.04,0.005, "10x error bars");
	mylatex->SetTextSize(0.04);
	merit_BDT->Draw("E");
	mylatex->Draw();

	c->cd(2);
	TLatex * mylatex1 = new TLatex(0.04,0.005, "10x error bars");
	mylatex1->SetTextSize(0.04);
	merit_SL->Draw("E");
	mylatex1->Draw();
	gStyle->SetTitleFontSize(0.06);
	c->Print("Q.pdf");

TFile * Qfile = new TFile("qfile.root", "RECREATE");
BDT_S->Write();
BDT_B->Write();
merit_BDT->Write();
merit_SL->Write();
Qfile->Close();

return 0;
}
