#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChainIndex.h>
#include <TSystem.h>
#include <TChain.h>
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


void hhev_maker(string inputfilename, int TriggerFlagSel, bool check_comb, float in1, float in2, string filename)
{
const double pi=3.14159265358979;
// CSV cut value, I use them on CMVA as there are not customary working points. The medium cut should be 0.7.
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
tree->SetBranchAddress("hJet_cmva", &(hJetCSV));
tree->SetBranchAddress("hJet_JECUnc", &(hJet_JECUnc));
tree->SetBranchAddress("hJet_genPt", &(hJet_genpT));
tree->SetBranchAddress("MET", &(metObj));
tree->SetBranchAddress("naJets", &(naJets));
tree->SetBranchAddress("aJet_e", &(aJetE));
tree->SetBranchAddress("aJet_pt", &(aJetpT));
tree->SetBranchAddress("aJet_eta", &(aJeteta));
tree->SetBranchAddress("aJet_phi", &(aJetphi));
tree->SetBranchAddress("aJet_cmva", &(aJetCSV));
tree->SetBranchAddress("aJet_JECUnc", &(aJet_JECUnc));
tree->SetBranchAddress("aJet_genPt", &(aJet_genpT));
tree->SetBranchAddress("Vtype", &(vType));
tree->SetBranchAddress("triggerFlags", &(triggerFlags));
tree->SetBranchAddress("aJet_nconstituents", &(aJet_tracks));
tree->SetBranchAddress("hJet_nconstituents", &(hJet_tracks));

int nEvents=tree->GetEntries();
	float pTcut = 20.;
	int incounter = 0;
	int outcounter = 0; 
	int FourJetEvents = 0;
	int FourJetEventsTrigger = 0;
	int aftercut = 0;
	int zero = 0;
	int uno = 0;
	const int trf = TriggerFlagSel;

//_____O_U_T_P_U_T__T_R_E_E_______________________________________________________
//string out_branch_name = "branch_3csv_"+filename; //change name here for NOCUT
string out_branch_name = "NOCUT_"+filename; 
TFile out_branches(out_branch_name.c_str(),"recreate");
TTree dodo("albero","a Tree");
HHevent hhev;
// "All" variables
dodo.Branch("APt_min",&hhev.APt_min,"APt_min/F"); dodo.Branch("APt_mean",&hhev.APt_mean,"APt_mean/F"); dodo.Branch("APt_max",&hhev.APt_max,"APt_max/F");
dodo.Branch("AEta_min",&hhev.AEta_min,"AEta_min/F"); dodo.Branch("AEta_mean",&hhev.AEta_mean,"AEta_mean/F"); dodo.Branch("AEta_max",&hhev.AEta_max,"AEta_max/F");
dodo.Branch("ACSV_min",&hhev.ACSV_min,"ACSV_min/F"); dodo.Branch("ACSV_mean",&hhev.ACSV_mean,"ACSV_mean/F"); dodo.Branch("ACSV_max",&hhev.ACSV_max,"ACSV_max/F");
dodo.Branch("Acent",&hhev.Acent,"Acent/F");
// "4" In_variables
dodo.Branch("QPt_1",&hhev.QPt_1,"QPt_1/F"); dodo.Branch("QPt_2",&hhev.QPt_2,"QPt_2/F"); dodo.Branch("QPt_3",&hhev.QPt_3,"QPt_3/F"); dodo.Branch("QPt_4",&hhev.QPt_4,"QPt_4/F");
dodo.Branch("QEta_1",&hhev.QEta_1,"QEta_1/F"); dodo.Branch("QEta_2",&hhev.QEta_2,"QEta_2/F"); dodo.Branch("QEta_3",&hhev.QEta_3,"QEta_3/F"); dodo.Branch("QEta_4",&hhev.QEta_4,"QEta_4/F");
dodo.Branch("QCSV_1",&hhev.QCSV_1,"QCSV_1/F"); dodo.Branch("QCSV_2",&hhev.QCSV_2,"QCSV_2/F"); dodo.Branch("QCSV_3",&hhev.QCSV_3,"QCSV_3/F"); dodo.Branch("QCSV_4",&hhev.QCSV_4,"QCSV_4/F");
dodo.Branch("Qcent",&hhev.Qcent,"Qcent/F");

// "dijet" In_variables -> pick just dijet_1? pick mean?
dodo.Branch("DJ_1_Phi_aperture",&hhev.DJ_1_Phi_aperture,"DJ_1_Phi_aperture/F");
dodo.Branch("DJ_1_Eta_aperture",&hhev.DJ_1_Eta_aperture,"DJ_1_Eta_aperture/F");
dodo.Branch("DJ_1_R_aperture",&hhev.DJ_1_R_aperture,"DJ_1_R_aperture/F");
dodo.Branch("DJ_1_mass",&hhev.DJ_1_mass,"DJ_1_mass/F");
dodo.Branch("DJ_1_pt",&hhev.DJ_1_pt,"DJ_1_pt/F");
dodo.Branch("DJ_2_Phi_aperture",&hhev.DJ_2_Phi_aperture,"DJ_2_Phi_aperture/F");
dodo.Branch("DJ_2_Eta_aperture",&hhev.DJ_2_Eta_aperture,"DJ_2_Eta_aperture/F");
dodo.Branch("DJ_2_R_aperture",&hhev.DJ_2_R_aperture,"DJ_2_R_aperture/F");
dodo.Branch("DJ_2_mass",&hhev.DJ_2_mass,"DJ_2_mass/F");
dodo.Branch("DJ_2_pt",&hhev.DJ_2_pt,"DJ_2_pt/F");
//"2 dijets" variables
dodo.Branch("TDJ_pt",&hhev.TDJ_pt,"TDJ_pt/F"); //vector sum
dodo.Branch("TDJ_deltaPt",&hhev.TDJ_deltaPt,"TDJ_deltaPt/F");
dodo.Branch("TDJ_deltaPhi",&hhev.TDJ_deltaPhi,"TDJ_deltaPhi/F");
dodo.Branch("TDJ_deltaEta",&hhev.TDJ_deltaEta,"TDJ_deltaEta/F");
dodo.Branch("TDJ_deltaR",&hhev.TDJ_deltaR,"TDJ_deltaR/F");
// already one dimensional variables 
dodo.Branch("met",&hhev.met,"met/F");
dodo.Branch("min_3csv",&hhev.min_3csv,"min_3csv/F");
dodo.Branch("avg_3csv",&hhev.avg_3csv,"avg_3csv/F");
dodo.Branch("costhetast", &hhev.costhetast, "costhetastst/F" );
dodo.Branch("costhetaCS", &hhev.costhetaCS, "costhetastCS/F" );
dodo.Branch("tau_1", &hhev.tau_1, "tau_1/F" );
dodo.Branch("tau_2", &hhev.tau_2, "tau_2/F" );
// number of jets
dodo.Branch("HHM",&hhev.HHM,"HHM/F");
dodo.Branch("ThirdJetTracks",&hhev.ThirdJetTracks,"ThirdJetTracks/F");
dodo.Branch("FourthJetTracks",&hhev.FourthJetTracks,"FourthJetTracks/F");
dodo.Branch("JetsN",&hhev.JetsN,"JetsN/F");


//________________________________________________________________________________
for (int i=0; i<nEvents; ++i)
{	
	//skip events affected by pixel misaligment
	if (trf != 0){ if ((207883<=EVENT.run && EVENT.run<=208307)) continue;}
	float progress = 100.0*((float)i)/((float)nEvents);
	if(!(i%100)) cout << setprecision(3) << progress << "% \r" ; 
	
  	tree->GetEvent(i);
  	if(nhJets+naJets < 4) {continue;}
	FourJetEvents++;
  	vector<JET> EVJets;
	float numj = 0;
//_____T_R_I_G_G_E_R___________________________________________________________
	if (triggerFlags[trf]==true) uno++;
	if (triggerFlags[trf]==false) zero++;

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
	

	if(numj < 4) {continue;}


//_S_O_R_T_I_N_G______________________________________________________________________________________________________________
// Sort by csv, loop on list, pick first three that meet the csv cut
	std::sort(EVJets.begin(), EVJets.end(), ordinamento_csv);
	TLorentzVector JETvec[4];
	int jcount = 0;
	float btag_backup[4];
		for(int l=0; l<numj; l++) {
			if (jcount>=3) {break;} 
			if(EVJets[l].btag > bTagCSV_mediumCut) {
				JETvec[jcount]=fillTLorentzVector(EVJets[l].pt, EVJets[l].eta, EVJets[l].phi, EVJets[l].e);
				btag_backup[jcount]=EVJets[l].btag;
				EVJets[l].flag=true;			//already token flag!
				jcount++;
				if(jcount==3) {hhev.ThirdJetTracks = EVJets[l].nTr;}	
				}
		}
	if(jcount<3) continue;
	aftercut++;

//__M_A_T_C_H_I_N_G_____________________________________________________________________________________________________________
// Loop on all the remaining jets, find and save the best matching in pairs for each.

vector <double> list_leastsqrtS;
vector <int> list_min_k;
vector <int> list_jj;
for(int jj=0; jj<EVJets.size(); jj++)
{
	if(EVJets[jj].flag == true ) { continue; }
	JETvec[3] = fillTLorentzVector(EVJets[jj].pt, EVJets[jj].eta, EVJets[jj].phi, EVJets[jj].e);
	double s[4][4];
	for(int l=0; l<4; l++)
	{
		for(int j=0; j<4; j++)
		{
			s[l][j]=JETvec[l]*JETvec[j];
		} //even in the massless case s=2(p1*p2)
	}
// 3 possible couples - uso sqrt(s)
	vector <double> s_diff;
	s_diff.push_back(fabs(sqrt(s[0][1])-sqrt(s[2][3])));
	s_diff.push_back(fabs(sqrt(s[0][2])-sqrt(s[1][3])));
	s_diff.push_back(fabs(sqrt(s[0][3])-sqrt(s[1][2])));
	int min_k = distance(s_diff.begin(), min_element(s_diff.begin(), s_diff.end()));		//min_k individua miglior coppia
    list_leastsqrtS.push_back(s_diff[min_k]);												//ricorda minor ∆s con questo jj
    list_min_k.push_back(min_k);															//ricorda il min_k (miglior coppia) con jj
    list_jj.push_back(jj);																	//ricorda chi è jj

}
	int bestmatch = distance(list_leastsqrtS.begin(), min_element(list_leastsqrtS.begin(), list_leastsqrtS.end())); //indice di minor ∆S
	int bestjj = list_jj[bestmatch];											//..da cui ricavo chi era il jj analizzato
	int best_min_k = list_min_k[bestmatch];										// ...e qual'era il suo accoppiamento migliore.
	
	JETvec[3] = fillTLorentzVector(EVJets[bestjj].pt, EVJets[bestjj].eta, EVJets[bestjj].phi, EVJets[bestjj].e);
	btag_backup[3] = EVJets[bestjj].btag;
	hhev.FourthJetTracks = EVJets[bestjj].nTr;

	TLorentzVector diJet_1;
   	TLorentzVector diJet_2;

    int sel[4];
    if (best_min_k==0) {sel[0]=0; sel[1]=1; sel[2]=2; sel[3]=3;} 
	if (best_min_k==1) {sel[0]=0; sel[1]=2; sel[2]=1; sel[3]=3;} 
	if (best_min_k==2) {sel[0]=0; sel[1]=3; sel[2]=1; sel[3]=2;}
	if (best_min_k>3 || best_min_k<0) cout << "non va bene" << endl;

	diJet_1 = JETvec[sel[0]]+JETvec[sel[1]]; diJet_2 = JETvec[sel[2]]+JETvec[sel[3]];
// sort diJets by decreasing pt. DiJet_1 is the one with highest pt
	TLorentzVector temp;
	if (diJet_2.Pt()>diJet_1.Pt()) {temp=diJet_1; diJet_1 = diJet_2; diJet_2=temp;}
// jets and dijets are now defined. It's just a matter of calculare and print the kinematical quantities.

float allE4centrality;

	for(int l=0; l<numj; l++)
	{
		allE4centrality+=EVJets[l].e;
	}
//______________________________________________________________________________________________________________________
// I remove the 4 hh jets from all the others, from now on "all" = "all the others"
//______________________________________________________________________________________________________________________

//removal of hh EVJets 
EVJets.erase(EVJets.begin()+bestjj);

for(int jj=EVJets.size(); jj>=0; jj--)
{
	if(EVJets[jj].flag == true ) { EVJets.erase(EVJets.begin()+jj); }
}

	float allPTsum = 0; float allEsum = 0; float allEtasum = 0; float bb = 0;

	for(int l=0; l<EVJets.size(); l++)
	{
		allEtasum += EVJets[l].eta;
		allPTsum += EVJets[l].pt;
		allEsum += EVJets[l].e;
		bb += EVJets[l].btag;
	} 
//______________________________________________________________________________________________________________________	   

	
	float qPT_alg_sum = JETvec[0].Pt()+JETvec[1].Pt()+JETvec[2].Pt()+JETvec[3].Pt();
   // float inv_mass = sqrt(m1*m1 + m2*m2 + 2*(diJet_1*diJet_2));


//__M_A_S_S___C_U_T_S_________________________________________________________________
	if (true) {

	//all jets variables
	hhev.APt_max = std::max_element(EVJets.begin(), EVJets.end(), [](const JET& a, const JET& b){return (a.pt<b.pt);}) -> pt;	
	hhev.APt_min = std::min_element(EVJets.begin(), EVJets.end(), [](const JET& a, const JET& b){return (a.pt<b.pt);}) -> pt;
	hhev.APt_mean = allPTsum/(float)numj;

	hhev.AEta_max = std::max_element(EVJets.begin(), EVJets.end(), [](const JET& a, const JET& b){return (a.eta<b.eta);}) -> eta;	
	hhev.AEta_min = std::min_element(EVJets.begin(), EVJets.end(), [](const JET& a, const JET& b){return (a.eta<b.eta);}) -> eta;
	hhev.AEta_mean = allEtasum/(float)numj;

	hhev.ACSV_max = std::max_element(EVJets.begin(), EVJets.end(), [](const JET& a, const JET& b){return (a.btag<b.btag);}) -> btag;	
	hhev.ACSV_min = std::min_element(EVJets.begin(), EVJets.end(), [](const JET& a, const JET& b){return (a.btag<b.btag);}) -> btag;
	hhev.ACSV_mean = bb/(float)numj;
	hhev.Acent = allPTsum/allEsum;

	
// 4 jets variables, single jets sorted by csv
	hhev.QPt_1 = JETvec[0].Pt();	hhev.QPt_2 = JETvec[1].Pt();	hhev.QPt_3 = JETvec[2].Pt();	hhev.QPt_4 = JETvec[3].Pt();
	hhev.QEta_1 = JETvec[0].Eta();	hhev.QEta_2 = JETvec[1].Eta();	hhev.QEta_3 = JETvec[2].Eta();	hhev.QEta_4 = JETvec[3].Eta();
	hhev.QCSV_1 = btag_backup[0];	hhev.QCSV_2 = btag_backup[1];	hhev.QCSV_3 = btag_backup[2];	hhev.QCSV_4 = btag_backup[3];	
	hhev.avg_3csv = (btag_backup[0]+btag_backup[1]+btag_backup[2])/3.0;
	hhev.min_3csv = (float)*std::min_element(btag_backup, btag_backup+2);
// centrality four jets
	float pttot = JETvec[0].Pt()+JETvec[1].Pt()+JETvec[2].Pt()+JETvec[3].Pt();
	float Etot = JETvec[0].E()+JETvec[1].E()+JETvec[2].E()+JETvec[3].E();
	hhev.Qcent = pttot/Etot;

// dijet variables	
// mass
	hhev.DJ_1_mass=diJet_1.M();
	hhev.DJ_2_mass=diJet_2.M();
//pt
	hhev.DJ_1_pt = diJet_1.Pt();
	hhev.DJ_2_pt = diJet_2.Pt();
//phi aperture
	hhev.DJ_1_Phi_aperture=(float)(pi-fabs(fabs(JETvec[sel[0]].Phi()-JETvec[sel[1]].Phi())-pi));
	hhev.DJ_2_Phi_aperture=(float)(pi-fabs(fabs(JETvec[sel[2]].Phi()-JETvec[sel[3]].Phi())-pi));
//eta aperture
	hhev.DJ_1_Eta_aperture=fabs(JETvec[sel[0]].Eta()-JETvec[sel[1]].Eta());
	hhev.DJ_2_Eta_aperture=fabs(JETvec[sel[2]].Eta()-JETvec[sel[3]].Eta());
//twist (tau)
	hhev.tau_1 = atan(hhev.DJ_1_Phi_aperture/hhev.DJ_1_Eta_aperture);
	hhev.tau_2 = atan(hhev.DJ_2_Phi_aperture/hhev.DJ_2_Eta_aperture);
//R aperture
	hhev.DJ_1_R_aperture=JETvec[sel[0]].DeltaR(JETvec[sel[1]]);
	hhev.DJ_2_R_aperture=JETvec[sel[2]].DeltaR(JETvec[sel[3]]);
//2dj pt vect sum
	hhev.TDJ_pt = (diJet_1+diJet_2).Pt();
// hh inv mass
	hhev.HHM = (diJet_1+diJet_2).M();
//delta phi
	hhev.TDJ_deltaPhi = pi-fabs(fabs(diJet_1.Phi()-diJet_2.Phi())-pi);
//delta eta
	hhev.TDJ_deltaEta = fabs(diJet_1.Eta()-diJet_2.Eta());
//delta R
	hhev.TDJ_deltaR = diJet_1.DeltaR(diJet_2);
//delta pt
	hhev.TDJ_deltaPt = fabs(diJet_1.Pt()-diJet_2.Pt());
// single variables 
	hhev.JetsN = numj; 
	hhev.met = metObj.et;


	//angles in other frames
	hhev.costhetast = fabs(func_Dcosthetast(diJet_1, diJet_2));
	hhev.costhetaCS = fabs(func_costhetaCS(diJet_1, diJet_2));



		incounter++;
		dodo.Fill();
	}


	else{ outcounter++; }


} //close event loop
dodo.Write();
out_branches.Close();
//begin comment here for NOCUT

ofstream asciifile;
string asciifilename = "stats_"+filename+".txt";
//asciifile.open (asciifilename.c_str(), ios::out | ios::app);
asciifile.open (asciifilename.c_str());
asciifile << "THREE CSV" << endl;
asciifile << "total events:\t" << nEvents << endl;
asciifile << "with 4 jets:\t" << FourJetEvents << endl;
asciifile << "...after trigger\t" << FourJetEventsTrigger << endl;
asciifile << "...after cuts\t" << aftercut << endl;
asciifile << "inside m cut:\t" << incounter << endl;
asciifile << "outside m cut:\t" << outcounter << endl;
asciifile << "zero " << zero << endl;
asciifile << "uno " << uno << endl;
asciifile << "somma " << zero+uno << endl;
asciifile.close();


return;
}






