#include <iostream>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TAxis.h"
#include <vector>
#include <cmath>
#include <string>
#include <fstream>

using namespace std;
//float Validation_cut[11] = {0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52};
//float lower_cut[14] = {0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51};


float Validation_cut[1] = {0.45};
float Q_cut = 0.58;
float lower_cut = 0.39;
//float Validation_cut = 0.50;
float CSVcut = 0.679;

struct abcd
{
float nTr;
float eta;
float pt;
float bv;
float qcsv;
float mass;
int region; // A=1, A_val=2, B=3, C = 4 , C_val = 5, D = 6; 
abcd()
	{
		nTr=0;
		eta=0;
		pt=0;
		bv=0;
		qcsv=0;
		region = 0; // 0 = unassigned
		mass = 0;
	}
};
vector <double> Dcount(vector<abcd> &, int, int, int, int, float *);



int main()
{

float N[10];
vector <abcd> evU;
vector <abcd> evL;
int massbins = 23;
int bdtbins = 20;
TH2F *mass_scatt_ALL = new TH2F("BDT_vs_mass", "BDT vs mass", bdtbins, 0.25, 0.75, massbins, 40, 500); 
TH2F *mass_scatt_3CSV = new TH2F("BDT_vs_mass_3CSV", "BDT vs mass (3btag)", bdtbins, 0.25, 0.75, massbins, 40, 500); 
TH2F *mass_scatt_4CSV = new TH2F("BDT_vs_mass_4CSV", "BDT vs mass (4btag)", bdtbins, 0.25, 0.75, massbins, 40, 500); 
TH2F *mass_scatt_2CSV = new TH2F("BDT_vs_mass_2CSV", "BDT vs mass (2btag)", bdtbins, 0.25, 0.75, massbins, 40, 500);
TH2F *mass_scatt_TEST = new TH2F("BDT_vs_mass_TEST", "2-3 btag syst. error / 4 btag stat. unc. ", bdtbins, 40.25, 0.75, massbins, 40, 500); 
TH2F *MC_mass_scatt_4CSV = new TH2F("MC_BDT_vs_mass_4CSV", "MC BDT vs mass (4btag)", bdtbins, 0.25, 0.75, massbins, 40, 500); 
TH2F *MC_mass_scatt_3CSV = new TH2F("MC_BDT_vs_mass_3CSV", "MC BDT vs mass (3btag)", bdtbins, 0.25, 0.75, massbins, 40, 500);
TH2F *MC_mass_scatt_ALL = new TH2F("MC_BDT_vs_mass_TEST", "MC BDT vs mass (ALL)", bdtbins, 0.25, 0.75, massbins, 40, 500); 
gStyle->SetOptStat(0);
mass_scatt_ALL->GetXaxis()->SetTitle("BDT response"); mass_scatt_ALL->GetYaxis()->SetTitle("h mass");
mass_scatt_3CSV->GetXaxis()->SetTitle("BDT response"); mass_scatt_3CSV->GetYaxis()->SetTitle("h mass");
mass_scatt_4CSV->GetXaxis()->SetTitle("BDT response"); mass_scatt_4CSV->GetYaxis()->SetTitle("h mass");
mass_scatt_2CSV->GetXaxis()->SetTitle("BDT response"); mass_scatt_2CSV->GetYaxis()->SetTitle("h mass");
mass_scatt_TEST->GetXaxis()->SetTitle("BDT response"); mass_scatt_TEST->GetYaxis()->SetTitle("h mass");
MC_mass_scatt_4CSV->GetXaxis()->SetTitle("BDT response"); MC_mass_scatt_4CSV->GetYaxis()->SetTitle("h mass");
MC_mass_scatt_3CSV->GetXaxis()->SetTitle("BDT response"); MC_mass_scatt_3CSV->GetYaxis()->SetTitle("h mass");
MC_mass_scatt_ALL->GetXaxis()->SetTitle("BDT response"); MC_mass_scatt_ALL->GetYaxis()->SetTitle("h mass");



//-------------------Reading from files---------
TChain *upperTree = new TChain("albero");
upperTree->Add("../Samples_BDT/BDT_CR_NOCUT_DiJetPt_BJetPlusX_Run2012B-13.root");
upperTree->Add("../Samples_BDT/BDT_CR_NOCUT_DiJetPt_BJetPlusX_Run2012C-24.root");
upperTree->Add("../Samples_BDT/BDT_CR_NOCUT_DiJetPt_BJetPlusX_Run2012C-Pr.root");
upperTree->Add("../Samples_BDT/BDT_CR_NOCUT_DiJetPt_BJetPlusX_Run2012D-Pr.root");
	float u_tr, u_eta, u_pt, u_BDT, u_csv, u_m1, u_m2;
//	upperTree->SetBranchAddress("QPt_3", &u_pt);
//	upperTree->SetBranchAddress("QEta_3", &u_eta);
//	upperTree->SetBranchAddress("ThirdJetTracks", &u_tr);
	upperTree->SetBranchAddress("QCSV_4", &u_csv);
	upperTree->SetBranchAddress("BDT", &u_BDT);
	upperTree->SetBranchAddress("DJ_1_mass", &u_m1);
	upperTree->SetBranchAddress("DJ_2_mass", & u_m2);
	for(int i=0; i<upperTree->GetEntries(); i++)
	{
		upperTree->GetEntry(i);
		evU.push_back(abcd());
//		evU[i].nTr = u_tr;
//		evU[i].eta = u_eta;
//		evU[i].pt = u_pt;
		evU[i].bv = u_BDT;
		evU[i].qcsv = u_csv;
		evU[i].mass = (u_m1+u_m2)/2.0;
		mass_scatt_ALL->Fill(evU[i].bv, evU[i].mass);
		mass_scatt_2CSV->Fill(evU[i].bv, evU[i].mass);
	}

TChain *lowerTree = new TChain("albero");
lowerTree->Add("../Samples_BDT/BDT_NOCUT_DiJetPt_BJetPlusX_Run2012B-13Jul.root");
lowerTree->Add("../Samples_BDT/BDT_NOCUT_DiJetPt_BJetPlusX_Run2012C-24Aug.root");
lowerTree->Add("../Samples_BDT/BDT_NOCUT_DiJetPt_BJetPlusX_Run2012C-Promp.root");
lowerTree->Add("../Samples_BDT/BDT_NOCUT_DiJetPt_BJetPlusX_Run2012D-Promp.root");
	
	float d_tr, d_eta, d_pt, d_BDT, d_csv, d_m1, d_m2;
//	lowerTree->SetBranchAddress("QPt_3", &d_pt);
//	lowerTree->SetBranchAddress("QEta_3", &d_eta);
//	lowerTree->SetBranchAddress("ThirdJetTracks", &d_tr);
	lowerTree->SetBranchAddress("QCSV_4", &d_csv);
	lowerTree->SetBranchAddress("BDT", &d_BDT);
	lowerTree->SetBranchAddress("DJ_1_mass", & d_m1);
	lowerTree->SetBranchAddress("DJ_2_mass", & d_m2);

	for(int i=0; i<lowerTree->GetEntries(); i++)
	{
		lowerTree->GetEntry(i);
		evL.push_back(abcd());
//		evL[i].nTr = d_tr;
//		evL[i].eta = fabs(d_eta);
//		evL[i].pt = d_pt;
		evL[i].bv = d_BDT;
		evL[i].qcsv = d_csv;
		evL[i].mass = (d_m1+d_m2)/2.0;
		mass_scatt_ALL->Fill(evL[i].bv, evL[i].mass);
		if(evL[i].qcsv < CSVcut) mass_scatt_3CSV->Fill(evL[i].bv, evL[i].mass);
		if(evL[i].qcsv >= CSVcut) mass_scatt_4CSV->Fill(evL[i].bv, evL[i].mass);
	}

TH2F *diff = new TH2F("twothree_btag_Residuals", "2-3 btag Residuals", bdtbins, 0.25, 0.75, massbins, 40, 500);
int TOTbins = mass_scatt_2CSV->GetXaxis()->GetNbins()*mass_scatt_2CSV->GetYaxis()->GetNbins();

float delta, unc;
float N2 = mass_scatt_2CSV->GetEntries();
float N3 = mass_scatt_3CSV->GetEntries();
float bin2, bin3, bin4;
float test;
for(int i=1; i<=TOTbins; i++)
{
	bin2 = mass_scatt_2CSV->GetBinContent(i);
	bin3 = mass_scatt_3CSV->GetBinContent(i);
	bin4 = mass_scatt_4CSV->GetBinContent(i);

	if (bin2 == 0 || bin3 == 0) {delta = 0.0; unc = 1.1;}
	else { 
		delta = bin2 - (bin3*(N2/N3));
		unc = sqrt(bin2 + (N2/N3)*(N2/N3)*bin3 + (bin3/N3)*(bin3/N3)*N2 + ((bin3*N2)/(N3*N3))*((bin3*N2)/(N3*N3))*N3 );
		if (bin4 != 0) {test = ((((N3/N2)*bin2)-bin3)/bin3)*sqrt(bin4);}
		if (bin4 == 0) {test = 0.0;}
		}
		diff->SetBinContent(i,delta/unc);
		mass_scatt_TEST->SetBinContent(i,test);
}


TFile * MC = new TFile("../Samples_BDT/BDT_NOCUT_hh4b-nores-8Tev300k_step2_L1y1.root", "READ");
TTree * leggo = (TTree*)MC->Get("albero");

vector <abcd> evMC;

	float mc_BDT, mc_csv, mc_m1, mc_m2;
//	lowerTree->SetBranchAddress("QPt_3", &d_pt);
//	lowerTree->SetBranchAddress("QEta_3", &d_eta);
//	lowerTree->SetBranchAddress("ThirdJetTracks", &d_tr);
	leggo->SetBranchAddress("QCSV_4", &mc_csv);
	leggo->SetBranchAddress("BDT", &mc_BDT);
	leggo->SetBranchAddress("DJ_1_mass", & mc_m1);
	leggo->SetBranchAddress("DJ_2_mass", & mc_m2);

	for(int i=0; i<leggo->GetEntries(); i++)
	{
		leggo->GetEntry(i);
		evMC.push_back(abcd());
		evMC[i].bv = mc_BDT;
		evMC[i].qcsv = mc_csv;
		evMC[i].mass = (mc_m1+mc_m2)/2.0;
		MC_mass_scatt_ALL->Fill(evMC[i].bv, evMC[i].mass);
		if(evMC[i].qcsv < CSVcut) MC_mass_scatt_3CSV->Fill(evMC[i].bv, evMC[i].mass);
		if(evMC[i].qcsv >= CSVcut) MC_mass_scatt_4CSV->Fill(evMC[i].bv, evMC[i].mass);
	}
MC->Close();


TFile *outhisto = new TFile("scatters.root", "RECREATE");
mass_scatt_ALL->Write();
mass_scatt_2CSV->Write();
mass_scatt_3CSV->Write();
mass_scatt_4CSV->Write();
mass_scatt_TEST->Write();
MC_mass_scatt_3CSV->Write();
MC_mass_scatt_4CSV->Write();
MC_mass_scatt_ALL->Write();
diff->Write();
outhisto->Close();

return 0;
}

/*


// Region definition
// |	1	|	2	|	3	|	2btags
// |		|		|		|
// |--------|-------|-------|
// |	4	|	5	|	6	|   >= 3btags
// |		|		|		|
// 	       Vcut    Qcut


for(int nn=0; nn<10; nn++) N[nn]=0.0;

for (int l=0; l<1; l++)
{
	if (lower_cut >= Validation_cut[l]) {continue;}
	for(int i=0; i<evU.size(); ++i){
		if(evU[i].bv > lower_cut && evU[i].bv < Validation_cut[l]) {evU[i].region = 1; N[1]++;}
		if(evU[i].bv < Q_cut && evU[i].bv > Validation_cut[l] ) {evU[i].region = 2; N[2]++;}
		if(evU[i].bv >= Q_cut) {evU[i].region = 3; N[3]++;}
	}
	for(int i=0; i<evL.size(); ++i){
		//if(evL[i].qcsv < CSVcut){
			if(evL[i].bv > lower_cut && evL[i].bv < Validation_cut[l]) {evL[i].region = 4; N[4]++;}
			if(evL[i].bv < Q_cut && evL[i].bv > Validation_cut[l] ) {evL[i].region = 5; N[5]++;}
			if(evL[i].bv >= Q_cut) {evL[i].region = 6; N[6]++; }
		//}
		//if(evL[i].qcsv >=CSVcut){
		//	if(evL[i].bv > lower_cut && evL[i].bv < Validation_cut) {evL[i].region = 7; N[7]++;}
		// 	if(evL[i].bv < Q_cut && evL[i].bv > Validation_cut[m] ) {evL[i].region = 8; N[8]++;}
		//	if(evL[i].bv >= Q_cut) {evL[i].region = 9; N[9]++; }
		//}
	}
vector <abcd> ev;
ev.reserve(evU.size()+evL.size());
ev.insert(ev.end(), evU.begin(), evU.end()); 
ev.insert(ev.end(), evL.begin(), evL.end()); 

vector <double> app_results(3);
vector <double> check_results(3);
vector <double> temp_results(3);
//-------------------
int A,B,C,D;
cout << "\n\n\n=================================================================" << endl;
cout << "Lower: " << lower_cut << "\t Validation: "<< Validation_cut[l] << endl;

	A=1; B=2; C=4; D=5;
	check_results = Dcount(ev, A,B,C,D,N);

	A=1; B=3; C=4; D=6;
	app_results = Dcount(ev, A,B,C,D,N);

double c = check_results[1];
double d = check_results[2];
cout << "c " << c << " d " << d << endl;
double outcome = app_results[0];
double syst_err = app_results[0]*sqrt(fabs(c*c - d*d));
double stat_err =  app_results[0]*app_results[1]; 



cout << "\n\n===========================================" << endl;
cout << "Events in D region:\t" << N[6] << endl;
cout << "Background expected:\t" << outcome << " +/- " << stat_err << " (stat.) +/- " << syst_err << " (syst.)" << endl;

	for(int i=1; i<7; i++)
	{
		cout << i << " " << N[i] << endl;
	}
cout << "_____________END_OF_RUN__________ " << endl;


string number = std::to_string(Validation_cut[l]);
string othernumber = std::to_string(lower_cut);
string canvasname =  "masses_L" + othernumber.substr(2,2) +  "_Vcut" + number.substr (2,2) +".pdf";

TCanvas *canv = new TCanvas("canv","canv",800,800);
TH1F * Amass = new TH1F("h_mass_A", "h_mass_A", 100, 0, 500);
TH1F * Bmass = new TH1F("h_mass_B", "h_mass_B", 100, 0, 500);
TH1F * Cmass = new TH1F("h_mass_C", "h_mass_C", 100, 0, 500);
TH1F * Dmass = new TH1F("h_mass_D", "h_mass_D", 100, 0, 500);

for(int i=0; i<ev.size(); i++)
{
	if(ev[i].region == 1) Amass->Fill(ev[i].mass);
	if(ev[i].region == 3) Bmass->Fill(ev[i].mass);
	if(ev[i].region == 4) Cmass->Fill(ev[i].mass);
	if(ev[i].region == 6) Dmass->Fill(ev[i].mass);
}

canv->Divide(2,2);
canv->cd(1);
Amass->Draw();
canv->cd(2);
Bmass->Draw();
canv->cd(3);
Cmass->Draw();
canv->cd(4);
Dmass->Draw();
canv->Print(canvasname.c_str());
delete canv; delete Amass; delete Bmass; delete Cmass; delete Dmass;
ev.clear();
} // end for on val cut
return 0;
}




vector<double> Dcount(vector<abcd> & ev, int A_reg, int B_reg, int C_reg, int D_reg, float *N)
{
// TH3: 25 bins for |eta|, 20 for tracks no, 125 for pt
	int etaBins = 20; int nTrBins = 14; int ptBins = 150;
	TH3F * bin_A = new TH3F("A_binned", "A_binned", etaBins, 0, 3.5, nTrBins, 0, 80, ptBins, 20, 620);
	TH3F * bin_B_val = new TH3F("B_val", "B_val", etaBins, 0, 3.5, nTrBins, 0, 80, ptBins, 20, 620);
	TH3F * bin_B = new TH3F("B_binned", "B_binned", etaBins, 0, 3.5, nTrBins, 0, 80, ptBins, 20, 620);
	TH3F * bin_C = new TH3F("C_binned", "C_binned", etaBins, 0, 3.5, nTrBins, 0, 80, ptBins, 20, 620);
	TH3F * bin_D_val = new TH3F("D_val", "D_val", etaBins, 0, 3.5, nTrBins, 0, 80, ptBins, 20, 620);
	TH3F * bin_D = new TH3F("D_binned", "D_binned", etaBins, 0, 3.5, nTrBins, 0, 80, ptBins, 20, 620);

	for(int i=0; i<ev.size(); i++)
	{
		if(ev[i].region==A_reg) bin_A->Fill(fabs(ev[i].eta), ev[i].nTr, ev[i].pt);
		if(ev[i].region==B_reg) bin_B->Fill(fabs(ev[i].eta), ev[i].nTr, ev[i].pt);
		if(ev[i].region==C_reg) bin_C->Fill(fabs(ev[i].eta), ev[i].nTr, ev[i].pt);
		if(ev[i].region==D_reg) bin_D->Fill(fabs(ev[i].eta), ev[i].nTr, ev[i].pt);
	}

double ND_exp = 0.0;
for(int i=1; i<=etaBins; i++) { // eta bin
	for(int j=1; j<=nTrBins; j++) {
		for(int k=1; k<=ptBins; k++){
			if(bin_A->GetBinContent(i,j,k) != 0){
				ND_exp += ((bin_C->GetBinContent(i,j,k)/bin_A->GetBinContent(i,j,k)))*bin_B->GetBinContent(i,j,k);}
		}
	}
}

double Ndx = N[B_reg] + N[D_reg];
double Nsx = N[A_reg]+N[C_reg];
double rat = Ndx/(Nsx+Ndx);
double perc = 100.0*fabs(N[D_reg]-ND_exp)/(N[D_reg]+ND_exp);
//float stat = ND_exp*sqrt( (1./(N[A_reg]+N[B_reg])) + (1./N[B_reg]) + (1./(N[C_reg]+N[D_reg])));

double stat = sqrt( (1./N[A_reg]) + (1./N[B_reg]) + (1./N[C_reg])  );
cout << "_______________________________________" << endl;
cout << " A = " << A_reg << "; B = " << B_reg << "; C = " << C_reg << "; D = " << D_reg << endl;
cout << "BD/(AC+BD) events ratio " << rat << endl;
cout << "ND \t NDexp \t stat \t differenza \%" << endl;
cout << N[D_reg] <<"\t"<< ND_exp << "\t" << stat << "\t" << perc << endl;
cout << "_______________________________________" << endl;

vector<double> result(3);
result[0] = ND_exp;						//ND expected
result[1] = fabs(N[D_reg]-ND_exp)/ND_exp; //errore relativo a NDexp
result[2] = stat;						//stat relative to NDexp
delete bin_A; delete bin_B_val; delete bin_B; delete bin_C; delete bin_D_val; delete bin_D;
return result;
}


*/
