
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TChain.h>
#include <TChainIndex.h>
#include <TLorentzVector.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include "myheader.h"

// #include <omp.h>

using namespace std;

//void fourpt(std::string, int, bool, float, float, std::string);
void hhev_maker(std::string, int, bool, float, float, std::string);
void plottini(std::string, int, bool, float, float, std::string);
void hhev_control_region(std::string, int, bool, float, float, std::string);

TLorentzVector fillTLorentzVector(double pT, double eta, double phi, double E)
{
  TLorentzVector jet;					
  jet.SetPtEtaPhiE(pT, eta, phi, E);	
  return jet;							
}

bool ordinamento_pt(const JET & a, const JET & b) 
{return a.pt > b.pt;}

bool ordinamento_csv(const JET & a, const JET & b) 
{return a.btag > b.btag;}


int main()
{
//	omp_set_num_threads(2);
	cout << "Carlo hhbbbb Analyzer 2.0" << endl;
	cout << "Enter 1 for MC, 2 for data ";
	int switcher;
	int TriggerFlagSel;
	cin >> switcher;
	switch (switcher)
	{
	case(1): TriggerFlagSel = 0; break;
	case(2): TriggerFlagSel = 54; break;
	default: cout << "error" << endl;
	}
	cout << "selected " << switcher << endl;

  	bool checkcomb = false;
  	float m_in1 = 100.0; 
  	float m_in2 = 150.0;

  	string dir;
  	string filename;
  	string inputfilename;
//______T_E_S_T__S_E_C_T_I_O_N_________________
 /*	dir = "";
	filename="OfficialStep2_RadionToHH_4b_M-650_TuneZ2star_8TeV_FULLSIM.root";
	inputfilename=dir+filename;
	fourpt(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename);
  	threecsv(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); */
 //_____________________________________________

dir = "/lustre/cmswork/tosi/ana/hh2bbbb/ntuple/slc5/";
for (int fn=0; fn<1; fn++)
{
	if(fn==0) {filename = "DiJetPt_BJetPlusX_Run2012B-13Jul2012-v1_EDMNtuple_V42_ProcV1.root"; 
				cout << "processing DATA file 1 of 2" << endl; inputfilename=dir+filename;
  	//hhev_maker(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "signal selection done" << endl;
  	plottini(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "ratios tests done" << endl;}
  //	hhev_control_region(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "control region done" << endl;}
	if(fn==1) {filename = "DiJetPt_BJetPlusX_Run2012C-PromptReco-v2_EDMNtuple_V42_ProcV1.root"; 
	//filename = "DiJetPt_BJetPlusX_Run2012D-PromptReco-v1_EDMNtuple_V42_ProcV1.root"
	cout << "processing DATA file 2 of 2" << endl; inputfilename=dir+filename;
  	hhev_maker(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "signal selection done" << endl;
  	plottini(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "ratios tests done" << endl;
  	hhev_control_region(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "control region done" << endl;}
}

/*
  	dir = "/lustre/cmswork/dallosso/hh2bbbb/non-resonant/event_generation/CMSSW_5_3_3_patch3/src/hh4b_8Tev_step2/";
	for (int fn=0; fn<8; fn++)
	{

	if(fn==0) {filename = "hh4b-nores-8Tev300k_step2_Lm10y23.root"; cout << "processing MC file 1 of 9" << endl; inputfilename=dir+filename;
  	hhev_maker(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "signal selection done" << endl;
  	plottini(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "ratios tests done" << endl;
	hhev_control_region(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename);}
	
	if(fn==1) {filename = "hh4b-nores-8Tev300k_step2_L2y16.root"; cout << "processing MC file 2 of 9" << endl; inputfilename=dir+filename;
  	hhev_maker(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "signal selection done" << endl;
  	plottini(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "ratios tests done" << endl;
hhev_control_region(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename);}
	
	if(fn==2) {filename = "hh4b-nores-8Tev300k_step2_L2y1.root"; cout << "processing MC file 3 of 9" << endl; inputfilename=dir+filename;
  	hhev_maker(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "signal selection done" << endl;
  	plottini(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "ratios tests done" << endl;
hhev_control_region(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename);}
	
	if(fn==3) {filename = "hh4b-nores-8Tev300k_step2_L2y075.root"; cout << "processing MC file 4 of 9" << endl; inputfilename=dir+filename;
  	hhev_maker(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "signal selection done" << endl;
  	plottini(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "ratios tests done" << endl;
hhev_control_region(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename);}
	
	if(fn==4) {filename = "hh4b-nores-8Tev300k_step2_L2y05.root"; cout << "processing MC file 5 of 9" << endl; inputfilename=dir+filename;
  	hhev_maker(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "signal selection done" << endl;
  	plottini(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "ratios tests done" << endl;
hhev_control_region(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename);}
	
	if(fn==5) {filename = "hh4b-nores-8Tev300k_step2_L20y1.root"; cout << "processing MC file 6 of 9" << endl; inputfilename=dir+filename;
  	hhev_maker(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "signal selection done" << endl;
  	plottini(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "ratios tests done" << endl;
hhev_control_region(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename);}
	
	if(fn==6) {filename = "hh4b-nores-8Tev300k_step2_L1y16.root"; cout << "processing MC file 7 of 9" << endl; inputfilename=dir+filename;
  	hhev_maker(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "signal selection done" << endl;
  	plottini(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "ratios tests done" << endl;
hhev_control_region(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename);}


	if(fn==40) {filename = "hh4b-nores-8Tev300k_step2_L1y1.root"; cout << "processing MC file 8 of 9" << endl; inputfilename=dir+filename;
	hhev_control_region(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename);
  	hhev_maker(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "signal selection done" << endl;
  	plottini(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "ratios tests done" << endl;}
	
	if(fn==7) {filename = "hh4b-nores-8Tev300k_step2_L1y05.root"; cout << "processing MC file 9 of 9" << endl; inputfilename=dir+filename;
  	hhev_maker(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "signal selection done" << endl;
	hhev_control_region(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename);
  	plottini(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename); cout << "ratios tests done" << endl;}

  } 
*/
return 0;
}


	//tree->Add("OfficialStep2_RadionToHH_4b_M-650_TuneZ2star_8TeV_FULLSIM.root"); 
  	//tree->Add("/lustre/cmswork/dallosso/hh2bbbb/non-resonant/event_generation/CMSSW_5_3_3_patch3/src/hh4b_8Tev_step2/hh4b-nores-8Tev300k_step2_L1y1.root"); 
	//tree->Add("/lustre/cmswork/tosi/ana/hh2bbbb/ntuple/slc5/DiJetPt_BJetPlusX_Run2012B-13Jul2012-v1_EDMNtuple_V42_ProcV1.root"); 
  	//std::cout<<"Opened input file "<<inputfilename<<std::endl;


	
/*
// MC file list
dir = "/lustre/cmswork/dallosso/hh2bbbb/non-resonant/event_generation/CMSSW_5_3_3_patch3/src/hh4b_8Tev_step2/";
#pragma omp parallel for
for (int fn=0; fn<9; fn++)
{
	if(fn==0) {filename = "hh4b-nores-8Tev300k_step2_Lm10y23.root"; cout << "processing MC file 1 of 9" << endl;}
	if(fn==1) {filename = "hh4b-nores-8Tev300k_step2_L2y16.root"; cout << "processing MC file 2 of 9" << endl;}
	if(fn==2) {filename = "hh4b-nores-8Tev300k_step2_L2y1.root"; cout << "processing MC file 3 of 9" << endl;}
	if(fn==3) {filename = "hh4b-nores-8Tev300k_step2_L2y075.root"; cout << "processing MC file 4 of 9" << endl;}
	if(fn==4) {filename = "hh4b-nores-8Tev300k_step2_L2y05.root"; cout << "processing MC file 5 of 9" << endl;}
	if(fn==5) {filename = "hh4b-nores-8Tev300k_step2_L20y1.root"; cout << "processing MC file 6 of 9" << endl;}
	if(fn==6) {filename = "hh4b-nores-8Tev300k_step2_L1y16.root"; cout << "processing MC file 7 of 9" << endl;}
	if(fn==7) {filename = "hh4b-nores-8Tev300k_step2_L1y1.root"; cout << "processing MC file 8 of 9" << endl;}
	if(fn==8) {filename = "hh4b-nores-8Tev300k_step2_L1y05.root"; cout << "processing MC file 9 of 9" << endl;}
}
// DATA file list
dir = "/lustre/cmswork/tosi/ana/hh2bbbb/ntuple/slc5/";

for (int fn=0; fn<2; fn++)
{
	if(fn==0) {filename = "DiJetPt_BJetPlusX_Run2012B-13Jul2012-v1_EDMNtuple_V42_ProcV1.root"; 
				cout << "processind DATA file 1 of 2" << endl; inputfilename=dir+filename;
  				fourpt(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename);
  				threecsv(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename);} 
	if(fn==1) {filename = "DiJetPt_BJetPlusX_Run2012C-24Aug2012-v2_EDMNtuple_V42_ProcV1.root"; 
				cout << "processind DATA file 2 of 2" << endl;
  				fourpt(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename);
  				threecsv(inputfilename, TriggerFlagSel, checkcomb, m_in1, m_in2, filename);} 
}
*/













