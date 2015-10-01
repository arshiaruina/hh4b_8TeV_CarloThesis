#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChainIndex.h>
#include <TSystem.h>
#include <TChain.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <TLorentzVector.h>
#include <fstream>

#include "myheader.h"
#include "HelperFunctions.h"

using namespace std;

void plottini(string inputfilename, int TriggerFlagSel, bool check_comb, float in1, float in2, string filename)
{
const double pi=3.14159265358979;
float bTagCSV_tightCut=0.898;
float bTagCSV_mediumCut=0.679;
float bTagCSV_looseCut=0.244;
float bTagCSV_noCut=0.;


TChain *tree=new TChain("tree");		
 tree->Add(inputfilename.c_str()); 

struct EventInfo
{
  int run;
  int lumi;
  int event;
  int json;
};


int nhJets, naJets, vType;
bool triggerFlags[500];
float hJetE[100], hJetpT[100], hJeteta[100], hJetphi[100], hJetCSV[100], hJet_genpT[100], hJet_JECUnc[100], hJet_tracks[100];
float aJetE[100], aJetpT[100], aJeteta[100], aJetphi[100], aJetCSV[100], aJet_genpT[100], aJet_JECUnc[100], aJet_tracks[100];
METInfo metObj;
float sigmaJERUnc = 0.;
float sigmaJECUnc = 0.;
EventInfo EVENT;
tree->SetBranchAddress("EVENT", &EVENT);
tree->SetBranchAddress("nhJets", &(nhJets));
tree->SetBranchAddress("hJet_e", &(hJetE)); 
tree->SetBranchAddress("hJet_pt", &(hJetpT));
tree->SetBranchAddress("hJet_eta", &(hJeteta));
tree->SetBranchAddress("hJet_phi", &(hJetphi));
tree->SetBranchAddress("hJet_csv", &(hJetCSV));
tree->SetBranchAddress("hJet_JECUnc", &(hJet_JECUnc));
tree->SetBranchAddress("hJet_genPt", &(hJet_genpT));
tree->SetBranchAddress("MET", &(metObj));
tree->SetBranchAddress("naJets", &(naJets));
tree->SetBranchAddress("aJet_e", &(aJetE));
tree->SetBranchAddress("aJet_pt", &(aJetpT));
tree->SetBranchAddress("aJet_eta", &(aJeteta));
tree->SetBranchAddress("aJet_phi", &(aJetphi));
tree->SetBranchAddress("aJet_csv", &(aJetCSV));
tree->SetBranchAddress("aJet_JECUnc", &(aJet_JECUnc));
tree->SetBranchAddress("aJet_genPt", &(aJet_genpT));
tree->SetBranchAddress("Vtype", &(vType));
tree->SetBranchAddress("triggerFlags", &(triggerFlags));
tree->SetBranchAddress("aJet_nconstituents", &(aJet_tracks));
tree->SetBranchAddress("hJet_nconstituents", &(hJet_tracks));

int nEvents=tree->GetEntries();
	float pTcut = 20.;
	int outcounter = 0; 
	int FourJetEvents = 0;
	int FourJetEventsTrigger = 0;
	int aftercut = 0;
	int zero = 0;
	int uno = 0;
	const int trf = TriggerFlagSel;

//_____O_U_T_P_U_T__H_I_S_T_O_S_______________________________________________________
//string out_branch_name = "branch_3csv_"+filename; //change name here for NOCUT
TH1F * pt2count = new TH1F("pt2", "pt2", 200, 0, 600); 
TH1F * pt3count = new TH1F("pt3", "pt3", 200, 0, 600); 
TH1F * eta2count = new TH1F("eta2", "eta2", 100, 0, 4.5); 
TH1F * eta3count = new TH1F("eta3", "eta3", 100, 0, 4.5); 
TH1F * ntracks2count = new TH1F("ntracks2", "ntracks2", 50, 0, 50); 
TH1F * ntracks3count = new TH1F("ntracks3", "ntracks3", 50, 0, 50); 

//________________________________________________________________________________
for (int i=0; i<nEvents; ++i)
{	
	//if (trf==0){ if(i>=207883 && i<=208307) continue;} //Pixel misaligment
	if (trf != 0){ if ((207883<=EVENT.run && EVENT.run<=208307)) continue;}
	float progress = 100.0*((float)i)/((float)nEvents);
	if(!(i%100)) cout << setprecision(3) << progress << "% \r" ; 
	
  	tree->GetEvent(i);
  	if(nhJets+naJets < 4) {continue;}
	FourJetEvents++;
  	vector<JET> EVJets;
	float numj = 0;
//_____T_R_I_G_G_E_R___________________________________________________________

	if (i==5) cout << "confirm trigger flag " << trf << endl;
	if (triggerFlags[trf]==false) continue;
//_____________________________________________________________________________
	FourJetEventsTrigger++;
	double b_temp;

  	for (int j=0; j<nhJets; ++j)
    {

    	if (trf==0) {hJetpT[j] = smear_pt_resErr(hJetpT[j], hJet_genpT[j], hJeteta[j], sigmaJERUnc)+sigmaJECUnc*hJet_JECUnc[j];}
    	if (hJetpT[j] < pTcut) continue;	//preliminary cuts
    	EVJets.push_back(JET());
    	EVJets[numj].pt = hJetpT[j];
    	EVJets[numj].eta = hJeteta[j];
    	EVJets[numj].phi = hJetphi[j];
    	EVJets[numj].e = hJetE[j];
    	EVJets[numj].nTr = hJet_tracks[j];
    	b_temp = hJetCSV[j];
		if (b_temp < 0 ) {b_temp=0;}
		EVJets[numj].btag = b_temp;
		numj++;
    }

    for (int j=0; j<naJets; ++j)
    {
    	if (trf==0) {aJetpT[j] = smear_pt_resErr(aJetpT[j], aJet_genpT[j], aJeteta[j], sigmaJERUnc)+sigmaJECUnc*aJet_JECUnc[j];}
    	if (aJetpT[j] < pTcut) continue; //preliminary cuts
    	EVJets.push_back(JET());
    	EVJets[numj].pt = aJetpT[j];
		EVJets[numj].eta = aJeteta[j];
		EVJets[numj].phi = aJetphi[j];
		EVJets[numj].e = aJetE[j];
		EVJets[numj].nTr = aJet_tracks[j];
		b_temp = aJetCSV[j];
		if (b_temp < 0 ) {b_temp=0;}
		EVJets[numj].btag = b_temp;
		numj++;
	}
	



//_S_O_R_T_I_N_G_____________________________________________________________________________________________________________
	if(EVJets.size()<3) continue;
	std::sort(EVJets.begin(), EVJets.end(), ordinamento_csv);
		if(EVJets[0].btag>bTagCSV_mediumCut && EVJets[1].btag>bTagCSV_mediumCut) 
			{
			pt2count->Fill(EVJets[2].pt);
			eta2count->Fill(fabs(EVJets[2].eta));
			ntracks2count->Fill(EVJets[2].nTr);
				if(EVJets[2].btag>bTagCSV_mediumCut) 
				{
				pt3count->Fill(EVJets[2].pt);
				eta3count->Fill(fabs(EVJets[2].eta));
				ntracks3count->Fill(EVJets[2].nTr);
				}
			}

} //close event loop
TH1F * pt_histo = new TH1F("pt_histo", "pt_histo", 100, 0, 600);
TH1F * eta_histo = new TH1F("eta_histo", "eta_histo", 70,0, 3.5);
TH1F * ntr_histo = new TH1F("ntr_histo", "ntr_histo", 50, 0, 50);

for(int i=1; i<200;i++)
	{	float due = pt2count->GetBinContent(i); float tre = pt3count->GetBinContent(i);
		if(due!=0) {pt_histo->SetBinContent(i, tre/due);}
		else pt_histo->SetBinContent(i,0);
	}
for(int i=1; i<100;i++)
	{	float due = eta2count->GetBinContent(i); float tre = eta3count->GetBinContent(i);
		if(due!=0) {eta_histo->SetBinContent(i, tre/due);}
		else eta_histo->SetBinContent(i,0);
	}
for(int i=1; i<50;i++)
	{	float due = ntracks2count->GetBinContent(i); float tre = ntracks3count->GetBinContent(i);
		if(due!=0) {ntr_histo->SetBinContent(i, tre/due);}
		else ntr_histo->SetBinContent(i,0);
	}


TCanvas *ratios = new TCanvas("ratios","ratios", 1024, 300);
TCanvas *check = new TCanvas("check","check",1024,600);
ratios->Divide(3,1,0.008,0.008);
gStyle->SetOptStat(0);
pt_histo->SetTitle("(ev 3rd jet CSV>M / ev 2CSV>M) on pt; pT (GeV); N3/N2" );
eta_histo->SetTitle("(ev 3rd jet CSV>M / ev 2CSV>M) on eta ; #eta; N3/N2" );
ntr_histo->SetTitle("(ev 3rd jet CSV>M / ev 2CSV>M) on n tracks ; constituent tracks; N3/N2" );
pt_histo->GetYaxis()->SetLabelOffset(1);
eta_histo->GetYaxis()->SetLabelOffset(1);
ntr_histo->GetYaxis()->SetLabelOffset(1);
pt_histo->SetFillColor(kBlue-7);
eta_histo->SetFillColor(kBlue-7);
ntr_histo->SetFillColor(kBlue-7);

ratios->cd(1); pt_histo->Draw(); ratios->cd(2); eta_histo->Draw(); ratios->cd(3); ntr_histo->Draw();
string outname1 = "ratios"+filename+".pdf";
ratios->SaveAs(outname1.c_str());
check->Divide(3,3,0.008,0.008);
gStyle->SetOptStat(0);
check->cd(1); pt2count->Draw(); check->cd(2); eta2count->Draw(); check->cd(3); ntracks2count->Draw(); 
check->cd(4); pt3count->Draw(); check->cd(5); eta3count->Draw(); check->cd(6); ntracks3count->Draw(); 
string outname2 = "check"+filename+".pdf";
check->SaveAs(outname2.c_str());
return;
}


