#include <iostream>
#include <cstring>
#include <algorithm>
#include <TLorentzVector.h>


struct JET
	{
		float pt;
		float eta;
		float phi;
		float e;
		float btag;
		float nTr;
		bool flag;

		JET()
		{
			pt = 0;
			eta = 0;
			phi = 0;
			e = 0;
			btag = 0;
			nTr = 0;
			flag = false;
		}
	};

struct METInfo
	{
		float et;
		float sumet;
		float sig;
		float phi;
	};

struct HHevent
	{
		// "All" variables
		float APt_min; float APt_mean; float APt_max;
		float AEta_min; float AEta_mean; float AEta_max;
		float ACSV_min; float ACSV_mean; float ACSV_max; 
		float Acent;
		float min_3csv; float avg_3csv;
		// "4" In_variables
		float QPt_1;   float QPt_2;  float QPt_3;  float QPt_4;
		float QEta_1;  float QEta_2; float QEta_3; float QEta_4;
		float QCSV_1;  float QCSV_2; float QCSV_3; float QCSV_4;
		float Qcent;
		// "dijet" In_variables 
		float DJ_1_Phi_aperture;
		float DJ_1_Eta_aperture; 
		float DJ_1_R_aperture;
		float DJ_1_mass;
		float DJ_1_pt;
		float tau_1;
		float DJ_2_Phi_aperture;
		float DJ_2_Eta_aperture; 
		float DJ_2_R_aperture;
		float DJ_2_mass;
		float DJ_2_pt;
		float tau_2;

		//"2 dijets" variables
		float TDJ_pt;	//vector sum
		float TDJ_deltaPt; 
		float TDJ_deltaPhi;
		float TDJ_deltaEta;
		float TDJ_deltaR;
		float HHM;

		// already one dimensional variables 
		float met;
 		float JetsN;
 		float ThirdJetTracks;
		float FourthJetTracks;
 		float costhetast;
 		float costhetaCS;

	};

bool ordinamento_pt(const JET & a, const JET & b);
bool ordinamento_csv(const JET & a, const JET & b);

TLorentzVector fillTLorentzVector(double pT, double eta, double phi, double E);




