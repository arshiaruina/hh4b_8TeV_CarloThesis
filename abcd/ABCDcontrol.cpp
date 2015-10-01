
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
#include <algorithm>
#include <iterator>
#include <functional>
#include <cmath>
#include <fstream>
#include <omp.h>

using namespace std;
float Validation_cut[16] = {0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53};
//float lower_cut[14] = {0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51};


// float Validation_cut = 0.52;
float Q_cut = 0.55;
// float lower_cut = 0.50;
float CSVcut = 0.679;

struct abcd
{
float nTr;
float eta;
float pt;
float bv;
float qcsv;
int region; // A=1, A_val=2, B=3, C = 4 , C_val = 5, D = 6; 
abcd()
	{
		nTr=0;
		eta=0;
		pt=0;
		bv=0;
		qcsv=0;
		region = 0; // 0 = unassigned
	}
};
vector <double> Dcount(vector<abcd> &, int, int, int, int, double *);


int main()
{
double N[10];
vector <abcd> evU;
vector <abcd> evL;

//-------------------Reading from files---------
TChain *upperTree = new TChain("albero");
upperTree->Add("../Samples_BDT/BDT_CR_NOCUT_DiJetPt_BJetPlusX_Run2012B-13.root");
upperTree->Add("../Samples_BDT/BDT_CR_NOCUT_DiJetPt_BJetPlusX_Run2012C-24.root");
upperTree->Add("../Samples_BDT/BDT_CR_NOCUT_DiJetPt_BJetPlusX_Run2012C-Pr.root");
upperTree->Add("../Samples_BDT/BDT_CR_NOCUT_DiJetPt_BJetPlusX_Run2012D-Pr.root");

	float u_tr, u_eta, u_pt, u_BDT, u_csv;
	upperTree->SetBranchAddress("QPt_3", &u_pt);
	upperTree->SetBranchAddress("QEta_3", &u_eta);
	upperTree->SetBranchAddress("ThirdJetTracks", &u_tr);
	upperTree->SetBranchAddress("QCSV_4", &u_csv);
	upperTree->SetBranchAddress("BDT", &u_BDT);
	for(int i=0; i<upperTree->GetEntries(); i++)
	{
		upperTree->GetEntry(i);
		evU.push_back(abcd());
		evU[i].nTr = u_tr;
		evU[i].eta = u_eta;
		evU[i].pt = u_pt;
		evU[i].bv = u_BDT;
		evU[i].qcsv = u_csv;
	}

TChain *lowerTree = new TChain("albero");
lowerTree->Add("../Samples_BDT/BDT_NOCUT_DiJetPt_BJetPlusX_Run2012B-13Jul.root");
lowerTree->Add("../Samples_BDT/BDT_NOCUT_DiJetPt_BJetPlusX_Run2012C-24Aug.root");
lowerTree->Add("../Samples_BDT/BDT_NOCUT_DiJetPt_BJetPlusX_Run2012C-Promp.root");
lowerTree->Add("../Samples_BDT/BDT_NOCUT_DiJetPt_BJetPlusX_Run2012D-Promp.root");
	
	float d_tr, d_eta, d_pt, d_BDT, d_csv;
	lowerTree->SetBranchAddress("QPt_3", &d_pt);
	lowerTree->SetBranchAddress("QEta_3", &d_eta);
	lowerTree->SetBranchAddress("ThirdJetTracks", &d_tr);
	lowerTree->SetBranchAddress("QCSV_4", &d_csv);
	lowerTree->SetBranchAddress("BDT", &d_BDT);

	for(int i=0; i<lowerTree->GetEntries(); i++)
	{
		lowerTree->GetEntry(i);
		evL.push_back(abcd());
		evL[i].nTr = d_tr;
		evL[i].eta = fabs(d_eta);
		evL[i].pt = d_pt;
		evL[i].bv = d_BDT;
		evL[i].qcsv = d_csv;
	}


// Region definition
// |	1	|	2	|	3	|	    2btags
// |		|		|		|
// |-------	|-------|-------|
// |	4	|	5	|	6	|		3btags
// |		|		|		|
// |--------|-------|-------|
// |	7	|	8	|	9	|		4btags
// |		|		|		|
// 	      Vcut 		Qcut
float rightbinsV[17] = {0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54};
//float rightbinsL[15] = {0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51,0.52};



TFile * storage = new TFile("ValidationCut.root", "RECREATE");
TH1D * H1245 = new TH1D("H1245", "H1245 (NDexp-ND)/NDexp",14, rightbinsV); //, 14, rightbinsL );
TH1D * H1346 = new TH1D("H1346", "H1346 (NDexp-ND)/NDexp",14, rightbinsV); //, 14, rightbinsL );
TH1D * H4578 = new TH1D("H4578", "H4578 (NDexp-ND)/NDexp",14, rightbinsV); // 14, rightbinsL );
TH1D * H4679 = new TH1D("H4679", "H4679 syst. error",14, rightbinsV); // 14, rightbinsL );
H1245->SetTitle("1245"); H1346->GetXaxis()->SetTitle("Validation cut"); //H1346->GetYaxis()->SetTitle("Lower Cut");
H1346->SetTitle("1346"); H1346->GetXaxis()->SetTitle("Validation cut"); //H1346->GetYaxis()->SetTitle("Lower Cut");
H4578->SetTitle("4578"); H4578->GetXaxis()->SetTitle("Validation cut"); //H4578->GetYaxis()->SetTitle("Lower Cut");
H4679->SetTitle("4679"); H4679->GetXaxis()->SetTitle("Validation cut"); //H4679->GetYaxis()->SetTitle("Lower Cut");

for(int m = 0; m < 16; m++) {
	//for(int l = 0; l < 14; l++){
	//	if(lower_cut[l] >= Validation_cut[m]) {continue;}
		
		for(int nn=0; nn<10; nn++) N[nn]=0.0;
	for(int i=0; i<evU.size(); ++i){
		if(evU[i].bv <= Validation_cut[m]) {evU[i].region = 1; N[1]++;}
		if(evU[i].bv < Q_cut && evU[i].bv > Validation_cut[m] ) {evU[i].region = 2; N[2]++;}
		if(evU[i].bv >= Q_cut) {evU[i].region = 3; N[3]++;}
	}
	for(int i=0; i<evL.size(); ++i){
		if(evL[i].qcsv < CSVcut){
			if(evL[i].bv <= Validation_cut[m]) {evL[i].region = 4; N[4]++;}
			if(evL[i].bv < Q_cut && evL[i].bv > Validation_cut[m] ) {evL[i].region = 5; N[5]++;}
			if(evL[i].bv >= Q_cut) {evL[i].region = 6; N[6]++; }
		}
		if(evL[i].qcsv >=CSVcut){
			if(evL[i].bv <= Validation_cut[m]) {evL[i].region = 7; N[7]++;}
			if(evL[i].bv < Q_cut && evL[i].bv > Validation_cut[m] ) {evL[i].region = 8; N[8]++;}
			if(evL[i].bv >= Q_cut) {evL[i].region = 9; N[9]++; }
		}
	}
vector <abcd> ev;
ev.reserve(evU.size()+evL.size());
ev.insert(ev.end(), evU.begin(), evU.end()); 
ev.insert(ev.end(), evL.begin(), evL.end()); 

vector <double> app_results(4);
vector <double> check_results(4);
vector <double> temp_results(4);
//-------------------
int A,B,C,D;
cout << "\n\n\n=================================================================" << endl;
//cout << "Lower: " << lower_cut[l] << "\t Validation: "<< Validation_cut[m] << endl;
cout << "Validation: "<< Validation_cut[m] << endl;
	// A=1; B=2; C=4; D=5;
	// temp_results = Dcount(ev, A,B,C,D,N);

	// A=1; B=2; C=7; D=8;
	// temp_results = Dcount(ev, A,B,C,D,N);

	// A=2; B=3; C=5; D=6;
	// temp_results = Dcount(ev, A,B,C,D,N);

int binx = H4679->GetXaxis()->FindBin(Validation_cut[m]);
// int biny_last = H4679->GetYaxis()->FindBin(lower_cut[l]);


	A=1; B=2; C=4; D=5;
	temp_results = Dcount(ev, A,B,C,D,N);
	H1245->SetBinContent(binx, temp_results[2]);

	A=4; B=5; C=7; D=8;
	temp_results = Dcount(ev, A,B,C,D,N);
	H4578->SetBinContent(binx, temp_results[2]);

	A=1; B=3; C=4; D=6;
	check_results = Dcount(ev, A,B,C,D,N);	
	H1346->SetBinContent(binx, check_results[2]);

	A=4; B=6; C=7; D=9;
	app_results = Dcount(ev, A,B,C,D,N);

	double c = check_results[2];
	double d = check_results[3];
	double syst_err = 0.0;
	if ((c*c-d*d) > 0.0) {syst_err = app_results[0]*sqrt(c*c - d*d);}
	double stat_err =  app_results[1]; 
	double combined = sqrt(stat_err*stat_err+syst_err*syst_err);
	H4679->SetBinContent(binx, combined);


//0	ND expected											    APP: ND expected
//1 VAL: è l'errore statistico su ND_exp [relativo]			APP: errore stat su NDexp [relativo]
//2 VAL: è l'errore sistematico c [relativo]				APP: -----
//3 VAL: è d, l'errore sull'errore syst c [relativo] 		APP: -----


cout << "\n\n===========================================" << endl;
cout << "Events in D region:\t" << N[9] << endl;
cout << "Background expected:\t" << app_results[0] << " +/- " << stat_err << " (stat.) +/- " << syst_err << " (syst.)" << endl;

	for(int i=1; i<10; i++)
	{
		cout << i << " " << N[i] << endl;
	}
cout << "_____________END_OF_RUN__________ " << endl;

ev.clear();
}//}
H1245->Write();
H1346->Write();
H4578->Write();
H4679->Write();
storage->Close();
return 0;
}




vector<double> Dcount(vector<abcd> & ev, int A_reg, int B_reg, int C_reg, int D_reg, double * N)
{

 float bineta[10] = {0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.5, 3.2};
 float bintr[10] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90};	
 float binpt[11] = {20, 30, 40, 50, 60, 70, 80, 90, 100, 180, 620};

// TH3: 25 bins for |eta|, 20 for tracks no, 125 for pt
	TH3F * bin_A = new TH3F("A_binned", "A_binned", 9, bineta, 9, bintr, 10, binpt);
	TH3F * bin_B = new TH3F("B_binned", "B_binned", 9, bineta, 9, bintr, 10, binpt);
	TH3F * bin_C = new TH3F("C_binned", "C_binned", 9, bineta, 9, bintr, 10, binpt);
	TH3F * bin_D = new TH3F("D_binned", "D_binned", 9, bineta, 9, bintr, 10, binpt);

	for(int i=0; i<ev.size(); i++)
	{
		if(ev[i].region==A_reg) bin_A->Fill(fabs(ev[i].eta), ev[i].nTr, ev[i].pt);
		if(ev[i].region==B_reg) bin_B->Fill(fabs(ev[i].eta), ev[i].nTr, ev[i].pt);
		if(ev[i].region==C_reg) bin_C->Fill(fabs(ev[i].eta), ev[i].nTr, ev[i].pt);
		if(ev[i].region==D_reg) bin_D->Fill(fabs(ev[i].eta), ev[i].nTr, ev[i].pt);
	}

//bin_A->Print("all");
float alt_a = 0.0;
float culo = 0;
double ND_exp = 0.0;
int fx, bx, fy, by, fz, bz;
double statW = 0.;
for(int i=1; i<=9; i++) { // eta bin
	for(int j=1; j<=9; j++) {
		for(int k=1; k<=10; k++){

			if(bin_C->GetBinContent(i,j,k) == 0){continue;}
			if(bin_B->GetBinContent(i,j,k) == 0){continue;}
			if(bin_A->GetBinContent(i,j,k) != 0)
			{
				ND_exp += (bin_C->GetBinContent(i,j,k)/bin_A->GetBinContent(i,j,k))*bin_B->GetBinContent(i,j,k);
				statW += ((bin_C->GetBinContent(i,j,k)/bin_A->GetBinContent(i,j,k))*bin_B->GetBinContent(i,j,k))*
						((bin_C->GetBinContent(i,j,k)/bin_A->GetBinContent(i,j,k))*bin_B->GetBinContent(i,j,k))*
						((1./bin_B->GetBinContent(i,j,k)) + (1./bin_C->GetBinContent(i,j,k)) + (1./bin_A->GetBinContent(i,j,k)) );
			}
			if(bin_A->GetBinContent(i,j,k) == 0){

			bx=i-1; fx=i+1; by=j-1; fy=j+1; bz=k-1; fz=k+1; 
			if(bx==0) {bx=1;} if(fx==10){fx=9;} if(by==0) {by=1;} if(fy==10){fy=9;} if(bz==0) {bz=1;} if(fz==11){fz=10;}
			alt_a = (1./26.0)*(bin_A->GetBinContent(bx,j,k)+bin_A->GetBinContent(fx,j,k)+ //primi vicini
							bin_A->GetBinContent(i,by,k)+bin_A->GetBinContent(i,fy,k)+
							bin_A->GetBinContent(i,j,bz)+bin_A->GetBinContent(i,j,fz)+
							bin_A->GetBinContent(fx,fy,k)+bin_A->GetBinContent(bx,by,k)+ // piano z =0
							bin_A->GetBinContent(fx,by,k)+bin_A->GetBinContent(bx,fy,k)+ //
							bin_A->GetBinContent(fx,j,fz)+bin_A->GetBinContent(bx,j,fz)+
							bin_A->GetBinContent(fx,j,bz)+bin_A->GetBinContent(bx,j,bz)+
							bin_A->GetBinContent(i,fy,fz)+bin_A->GetBinContent(i,by,fz)+ //
							bin_A->GetBinContent(i,fy,bz)+bin_A->GetBinContent(i,by,bz)+  //
							bin_A->GetBinContent(fx,fy,fz)+bin_A->GetBinContent(bx,by,bz)+
							bin_A->GetBinContent(fx,by,fz)+bin_A->GetBinContent(fx,by,bz)+
							bin_A->GetBinContent(bx,fy,fz)+bin_A->GetBinContent(bx,by,fz)+
							bin_A->GetBinContent(bx,j,fz)+bin_A->GetBinContent(fx,j,bz)+
							bin_A->GetBinContent(fx,j,fz)+bin_A->GetBinContent(bx,j,bz));
			ND_exp += (bin_C->GetBinContent(i,j,k)/alt_a)*bin_B->GetBinContent(i,j,k);

			statW += ((bin_C->GetBinContent(i,j,k)/alt_a)*bin_B->GetBinContent(i,j,k))*
					((bin_C->GetBinContent(i,j,k)/alt_a)*bin_B->GetBinContent(i,j,k))*
					((1./bin_B->GetBinContent(i,j,k)) + (1./bin_C->GetBinContent(i,j,k)) + (1./alt_a) );

			if(alt_a==0) culo++;
				}
		}
	}
}

double Ndx = N[B_reg] + N[D_reg];
double Nsx = N[A_reg]+N[C_reg];
double rat = Ndx/(Nsx+Ndx);
double perc = 100.0*fabs(N[D_reg]-ND_exp)/(ND_exp);
//float stat = ND_exp*sqrt( (1./(N[A_reg]+N[B_reg])) + (1./N[B_reg]) + (1./(N[C_reg]+N[D_reg])));

double stat = sqrt( (1./N[A_reg]) + (1./N[B_reg]) + (1./N[C_reg])  );
cout << "_______________________________________" << endl;
cout << " A = " << A_reg << "; B = " << B_reg << "; C = " << C_reg << "; D = " << D_reg << endl;
cout << "BD/(AC+BD) events ratio " << rat << endl;
cout << "ND \t NDexp \t stat [rel] \t differenza \%" << endl;
cout << "number of bin_A = 0 " << culo << endl;
cout << N[D_reg] <<"\t"<< ND_exp << "\t" << sqrt(statW) << "\t" << perc << endl;
//cout << "VALUE CHECK - VALUE CHECK - VALUE CHECK - VALUE CHECK" << endl;
//for(int i=1; i<10; i++){ cout << i << " " << N[i] << endl;}

cout << "_______________________________________" << endl;

vector<double> result(4);
result[0] = ND_exp;							//ND expected											APP: ND expected
result[1] = sqrt(statW);								//VAL: è l'errore statistico su ND_exp [relativo]		APP: errore stat su ND [relativo]
result[2] = fabs(N[D_reg]-ND_exp)/ND_exp; 		//VAL: è l'errore sistematico c [relativo]				APP: -----
result[3] = sqrt((N[D_reg]/(ND_exp*ND_exp))+(stat*stat)); 	//VAL: è d, l'errore sull'errore syst c [relativo] 		APP: -----

delete bin_A; delete bin_B; delete bin_C; delete bin_D;
return result;
}
