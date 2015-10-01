#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <string>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h" 
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "Riostream.h"
#include "TVirtualFitter.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TColor.h"
#include "TMarker.h"

using namespace std;

struct ranking
{
	string label;
	float like;
	double kolgo;
	double ander;
	int KS_index;
	int AD_index;

	ranking()
	{
		label = "";
		like = 0;
		kolgo = 0;
		ander=0;
	}
};

bool ordinamento_like(const ranking & a, const ranking & b) 
{return a.like < b.like;}
bool ordinamento_kolgo(const ranking & a, const ranking & b) 
{return a.kolgo < b.kolgo;}
bool ordinamento_ander(const ranking & a, const ranking & b) 
{return a.ander < b.ander;}
double Factorial(int n) {return(n<=1. ? 1. : n*Factorial(n-1.));}
double Probability(double n, double mu) {return( n<25 ? pow(mu, n)*exp(-mu)/(double)Factorial(n) : (1./sqrt(2*3.1415926*mu))*exp(-0.5*pow((double)n-mu,2)/mu) );}
double MyAndersonDarlingTest (const TH1 *h1, const TH1 *h2);

vector<string> setTitle(int);

inline void plotCMStext(double fracerr = 0.1) {
  TLatex * latex = new TLatex();
  latex->SetNDC();                  //Set Normalized Device Coordinates
  latex->SetTextSize(0.03);         
  latex->DrawLatex(0.4,0.85,Form("Fract. error = %.4f",fracerr)); 
  delete latex;                                                    
}

string histo; string branch;
TH1D * S_templ_in;
TH1D * B_templ_in;
TH1D * ToySample_in;
int Nbins;
double minXvalue;
double maxXvalue;
double minXvalue_cfr;
double maxXvalue_cfr;

void Likelihood(int& npar, double* grad, double& fval, double* xval, int flag) {
  //  cout << "In likelihood" << endl;
  double fs=(1+sin(xval[0]))/2.;  // i cazzi nel culo     
  double flike=0.;
  double mu, n, lik;    
  double Ntot = ToySample_in->Integral();
  // cout << "ToySample integral = " << Ntot << " fs= " << fs << endl;
  for (int i=1; i<=Nbins; i++) {
    if (S_templ_in->GetBinContent(i)+B_templ_in->GetBinContent(i)<=0.) continue;
    n=ToySample_in->GetBinContent(i);
    double ns = S_templ_in->GetBinContent(i);
    double nb = B_templ_in->GetBinContent(i);
    mu=(fs*ns+(1.-fs)*nb)*Ntot; 
    lik=log(Probability(n, mu));
    if((lik!=lik)) { 
    cout << "Problems in likelihood" << endl;
   cout << "i=" << i << " n=" << n << " ns="<<ns << " nb=" << nb << " fs=" << fs << " mu=" << mu << endl;
   } else flike+=lik; //
  }
  fval=-flike; 
}

int main()
{
vector <ranking> rank;
vector <float> out_error;
vector <float> out_mKS;
vector <float> out_mAD;
bool override;
for(int nh=1; nh<47; nh++){
rank.push_back(ranking());
override = false;
if(nh==1) {histo = "Pt_all_jets"; branch = "APt_min" ;override=true; Nbins=200; minXvalue = 0.; maxXvalue=200; }
if(nh==2) {histo = "Pt_all_jets"; branch = "APt_mean";override=true; Nbins=200; minXvalue = 0.; maxXvalue=200; }
if(nh==3) {histo = "Pt_all_jets"; branch = "APt_max" ;}
if(nh==4) {histo = "Eta_all_jets"; branch = "AEta_min" ;}
if(nh==5) {histo = "Eta_all_jets"; branch = "AEta_mean" ;}
if(nh==6) {histo = "Eta_all_jets"; branch = "AEta_max" ;}
if(nh==7) {histo = "CSV_all_jets"; branch = "ACSV_min";}// override=true; Nbins=50; minXvalue = 0.; maxXvalue=0.5; }
if(nh==8) {histo = "CSV_all_jets"; branch = "ACSV_mean";}
if(nh==9) {histo = "CSV_all_jets"; branch = "ACSV_max"; override=true; Nbins=100; minXvalue = 0.0; maxXvalue=1.0; }
if(nh==10) {histo = "InCentrality"; branch = "Acent";}
if(nh==11) {histo = "InCentrality"; branch = "Qcent";}
if(nh==12) {histo = "In_Pt_4_jets"; branch = "QPt_1";}
if(nh==13) {histo = "In_Pt_4_jets"; branch = "QPt_2";}
if(nh==14) {histo = "In_Pt_4_jets"; branch = "QPt_3";}
if(nh==15) {histo = "In_Pt_4_jets"; branch = "QPt_4";}
if(nh==16) {histo = "In_Eta_4_jets"; branch = "QEta_1";}
if(nh==17) {histo = "In_Eta_4_jets"; branch = "QEta_2";}
if(nh==18) {histo = "In_Eta_4_jets"; branch = "QEta_3";}
if(nh==19) {histo = "In_Eta_4_jets"; branch = "QEta_4";}
if(nh==20) {histo = "In_CSV_4_jets"; branch = "QCSV_1"; override=true; Nbins=35; minXvalue = 0.65; maxXvalue=1.0; }
if(nh==21) {histo = "In_CSV_4_jets"; branch = "QCSV_2"; override=true; Nbins=35; minXvalue = 0.65; maxXvalue=1.0; }
if(nh==22) {histo = "In_CSV_4_jets"; branch = "QCSV_3"; override=true; Nbins=35; minXvalue = 0.65; maxXvalue=1.0; }
if(nh==23) {histo = "In_CSV_4_jets"; branch = "QCSV_4";}
if(nh==24) {histo = "In_dj_mass";  	branch = "DJ_1_mass";}
if(nh==25) {histo = "dj_Pt"; 		branch = "DJ_1_pt";}
if(nh==26) {histo = "In_dj_Phi_aperture";branch = "DJ_1_Phi_aperture";}
if(nh==27) {histo = "In_dj_Eta_aperture"; branch = "DJ_1_Eta_aperture";}
if(nh==28) {histo = "In_dj_R_aperture"; branch = "DJ_1_R_aperture";}
if(nh==29) {histo = "In_dj_mass";  	branch = "DJ_2_mass";}
if(nh==30) {histo = "dj_Pt"; 		branch = "DJ_2_pt";}
if(nh==31) {histo = "In_dj_Phi_aperture";branch = "DJ_2_Phi_aperture";}
if(nh==32) {histo = "In_dj_Eta_aperture"; branch = "DJ_2_Eta_aperture";}
if(nh==33) {histo = "In_dj_R_aperture"; branch = "DJ_2_R_aperture";}
if(nh==34) {histo = "2dj_pt_sum";  	branch = "TDJ_pt";}
if(nh==35) {histo = "In_2dj_deltaPhi"; branch = "TDJ_deltaPhi";}
if(nh==36) {histo = "In_2dj_deltaEta"; branch = "TDJ_deltaEta";}
if(nh==37) {histo = "In_2dj_deltaR"; branch = "TDJ_deltaR";}
if(nh==38) {histo = "In_hh_inv_mass";  branch = "HHM";}
if(nh==39) {histo ="h_met"; 		branch = "met";}
if(nh==40) {histo ="In_min_csv"; 	branch = "min_3csv";override=true; Nbins=40; minXvalue = 0.6 ; maxXvalue=1.0; }
if(nh==41) {histo ="In_mean_csv"; 	branch = "avg_3csv";override=true; Nbins=35; minXvalue = 0.65 ; maxXvalue=1.0; }
if(nh==42) {histo ="hh_scatter_angle_CM"; branch = "costhetast";}
if(nh==43) {histo ="hh_scatter_angle_CS"; branch = "costhetaCS";}
if(nh==44) {histo ="tau"; 			branch = "tau_1";}
if(nh==45) {histo ="tau";			branch = "tau_2";}
if(nh==46) {histo ="Jets_no"; 		branch = "JetsN"; override=true; Nbins=10; minXvalue = 4 ; maxXvalue=14;}


TFile *HISTO = new TFile("/Users/Carlo/Desktop/hh2bbbb/BACKUP_hhbbbb/hh_events_closed29May/histo_3csv_hh4b-nores-8Tev300k_step2_L1y1.root","");
TFile *MC = new TFile("/Users/Carlo/Desktop/Final_Software/NOCUT_hh4b-nores-8Tev300k_step2_L1y1.root", "");
TFile *DATA = new TFile("/Users/Carlo/Desktop/Final_Software/NOCUT_DiJetPt_BJetPlusX_Run2012B-13Jul2012-v1_EDMNtuple_V42_ProcV1.root", "");

  if (override == false)
  {
  TH1D *utility = (TH1D*)HISTO->Get(histo.c_str());
  Nbins = utility->GetXaxis()->GetNbins();
  minXvalue = utility->GetXaxis()->GetXmin();
  maxXvalue = utility->GetXaxis()->GetXmax();
  HISTO->Close();
 }
  TTree *ts = (TTree*)MC->Get("albero");
  TTree *tb = (TTree*)DATA->Get("albero");
  float temp_s; float temp_b; float s_m1; float s_m2; float b_m1; float b_m2;
  TBranch *s_branch = ts->GetBranch(branch.c_str());
  TBranch *b_branch = tb->GetBranch(branch.c_str());
  TBranch *s_m1_branch = ts->GetBranch("DJ_1_mass");
  TBranch *s_m2_branch = ts->GetBranch("DJ_2_mass");
  TBranch *b_m1_branch = tb->GetBranch("DJ_1_mass");
  TBranch *b_m2_branch = tb->GetBranch("DJ_2_mass");
  s_branch->SetAddress(&temp_s);
  b_branch->SetAddress(&temp_b);
  s_m1_branch->SetAddress(&s_m1);
  b_m1_branch->SetAddress(&b_m1);
  s_m2_branch->SetAddress(&s_m2);
  b_m2_branch->SetAddress(&b_m2);
  int s_entries = ts->GetEntries();
  int b_entries = tb->GetEntries();
  string s_label = "S_" + histo; 
  string b_label = "B_" + histo;
  TH1D * S_templ = new TH1D(s_label.c_str(), s_label.c_str(), Nbins, minXvalue, maxXvalue);
  TH1D * B_templ = new TH1D(b_label.c_str(), b_label.c_str(), Nbins, minXvalue, maxXvalue);

  for (long i=0;i<s_entries;i++) {
  s_branch->GetEntry(i);
  s_m1_branch->GetEntry(i);
  s_m2_branch->GetEntry(i);
  if(s_m1 <= 150.0 && s_m1>= 100.0 && s_m2 <= 150.0 && s_m2>= 100.0)
  {S_templ->Fill(temp_s); }
}

  for (long i=0;i<b_entries;i++) {
  b_branch->GetEntry(i);
  b_m1_branch->GetEntry(i);
  b_m2_branch->GetEntry(i);
  if(b_m1 <= 150.0 && b_m1>= 100.0 && b_m2 <= 150.0 && b_m2>= 100.0)
  {B_templ->Fill(temp_b); }
}
 // 2-samples Anderson-Darling
 double AD_A = MyAndersonDarlingTest(S_templ, B_templ);

	S_templ_in = (TH1D*)S_templ->Clone();
	B_templ_in = (TH1D*)B_templ->Clone();
	
	gStyle->SetOptStat(0);
 	gStyle->SetOptFit(0);
	TH1D * ToySample =  new TH1D ("ToySample", "Measured distribution", Nbins, minXvalue, maxXvalue);
	TH1D * Fit_result = new TH1D ("Fit_result", "Fit result", Nbins, minXvalue, maxXvalue);

double Ntest1_old = S_templ->GetEntries();
double Ntest2_old = B_templ->GetEntries();
double ffs = 0.20;
double supprB = ((1./ffs)-1)*(Ntest1_old/Ntest2_old); //3. / 0.33333 25 % signal 75 % bkg
double supprS = ((1./(1-ffs))-1)*(Ntest2_old/Ntest1_old); //
double Ntest1, Ntest2;
if (Ntest1_old<Ntest2_old) {Ntest2 = supprB*Ntest2_old; supprS = 1; Ntest1 = Ntest1_old;}
if (Ntest1_old>Ntest2_old) {Ntest1 = supprS*Ntest1_old; supprB = 1; Ntest2 = Ntest2_old;}


	for (int i=1; i<=Nbins; i++) {
    	double s = supprS*S_templ->GetBinContent(i);
    	double b = supprB*B_templ->GetBinContent(i);
    	ToySample->SetBinContent(i,gRandom->Poisson(s+b));
  	}
ToySample_in = (TH1D*)ToySample->Clone();

TH1D * cdf_S = (TH1D*)S_templ_in->GetCumulative();
TH1D * cdf_B = (TH1D*)B_templ_in->GetCumulative();
cdf_S->Scale(1./S_templ_in->GetEntries());
cdf_B->Scale(1./B_templ_in->GetEntries());
double n1 = S_templ_in->GetEntries(); double n2 = B_templ_in->GetEntries();

S_templ_in->Scale(1./Ntest1_old);
B_templ_in->Scale(1./Ntest2_old);

 
// Likelihood maximization
// ----------------------- 
	double Fsin = (double)Ntest1/(double)(Ntest1+Ntest2);
// Minuit routine
  	TMinuit rmin(1);
  	rmin.SetFCN(Likelihood);
  	double arglis[4];
  	arglis[0]=1;
  	int error_flag=0;
  	rmin.mnexcm("SET ERR",arglis,1,error_flag);
  	rmin.mninit(5,6,7);
  	double f_null=0.0;
  	error_flag = 0;
		TString parnamEj[1]={"fs"};
		double Start[1]={asin(2.*Fsin-1.)};
		double Step[1]={0.00001};
		double Min[1]={-3.1415926};
		double Max[1]= {3.1415926}; 
	rmin.mnparm(0, parnamEj[0],  Start[0], Step[0], Min[0], Max[0], error_flag);
	rmin.mnexcm( "MIGRAD",0,0,error_flag);
// Read results
	double a[1], err[1], pmin, pmax; int ivar;
	rmin.mnpout (0, parnamEj[0], a[0], err[0], pmin, pmax, ivar);
	long double Fraction1, Error1;  
  	Fraction1=(1.+sin(a[0]))/2.;
  	Error1=fabs(cos(a[0])*err[0])/2.;
  for (int i=1; i<=Nbins; i++) {
    double cont = S_templ_in->GetBinContent(i)*Fraction1*ToySample->Integral()+B_templ_in->GetBinContent(i)*(1.-Fraction1)*ToySample->Integral();
    Fit_result->SetBinContent(i,cont);
    Fit_result->SetBinError(i,cont*Error1/Fraction1);
    double error = Error1*sqrt(pow(S_templ_in->GetBinContent(i)*ToySample->Integral(),2)+pow(B_templ_in->GetBinContent(i)*ToySample->Integral(),2));
    Fit_result->SetBinError(i,error);
}

// 2-samples Kolmogorov Smirnoff test
    vector <double> cdf_diff;
    for(int j=1; j<=Nbins; j++) {
    	cdf_diff.push_back( fabs(cdf_S->GetBinContent(j)-cdf_B->GetBinContent(j)) );
    	}
   double factor = 1./sqrt((n1*n2)/(n1+n2));
   factor = 1.;
   cout << "factor " << factor << endl;
   double maxD = factor*(*max_element(cdf_diff.begin(), cdf_diff.end()));
   cout << "KS distance: " << maxD << endl;
   cout << "AD test " << AD_A << endl;
//--------------------------------------
// Tests correlation plot
//--------------------------------------

  out_error.push_back(Error1/Fraction1);
  out_mKS.push_back(-maxD);
  out_mAD.push_back(-AD_A);


//--------------------------------------
// New output
//--------------------------------------
rank[nh-1].label = branch;
rank[nh-1].like = Error1/Fraction1;
rank[nh-1].kolgo = -maxD;
rank[nh-1].ander = -AD_A;

if(nh==25)
{
TCanvas * pony = new TCanvas ("pony","pony", 800, 800);
cdf_S->Draw();
cdf_B->Draw("SAME");
pony->SaveAs("Cumulative_example.pdf");
}
  TCanvas * Res = new TCanvas ("Res","Res", 800, 800);
  Res->cd();
  ToySample_in->SetMinimum(0);
  ToySample_in->SetMarkerStyle(5);
  ToySample_in->SetMarkerSize(1.);
  ToySample_in->SetMarkerColor(kMagenta);
  ToySample_in->SetLineWidth(2);
  ToySample_in->SetTitle(setTitle(nh)[0].c_str()); 
  ToySample_in->SetXTitle(setTitle(nh)[1].c_str());
  //ToySample_in->SetYTitle(setTitle(nh)[2].c_str());
  ToySample_in->Draw("PE");
  S_templ_in->Scale(ToySample->Integral()*Fraction1);
  B_templ_in->Scale(ToySample->Integral()*(1-Fraction1));
  S_templ_in->SetLineColor(kRed);
  S_templ_in->SetLineWidth(2);
  S_templ_in->Draw("SAME");
  B_templ_in->SetLineColor(kBlue);
  B_templ_in->SetLineWidth(2);
  B_templ_in->Draw("SAME");
  Fit_result->SetLineWidth(2);
  Fit_result->SetLineColor(kGreen);
  Fit_result->Draw("SAMEHISTO");  
plotCMStext(Error1/Fraction1);
  TLegend *legend=new TLegend(0.7,0.75,0.9,0.9);
  //legend1->SetTextFont(50);
  legend->SetTextSize(0.03);
  legend->AddEntry(S_templ_in,"Signal" ,"f");
  legend->AddEntry(B_templ_in,"Background" ,"f");
  legend->AddEntry(ToySample_in,"ToySample"  ,"pe");
  legend->AddEntry(Fit_result,"Fit"  ,"f");
  legend->Draw();
  string pic_name_pdf = branch + ".pdf";
  Res->SaveAs(pic_name_pdf.c_str());
}

//--------------------------------------
// Output
//--------------------------------------

TGraph *L_KS = new TGraph(46, &out_error[0], &out_mKS[0]);
L_KS->SetTitle("K-S vs Likelihood meth. ; L ; -KS");
L_KS->SetMarkerSize(1);
L_KS->SetMarkerStyle(20);

TGraph *L_AD = new TGraph(46, &out_error[0], &out_mAD[0]);
L_AD->SetTitle("A-D vs Likelihood meth.; L; -AD");
L_AD->SetMarkerSize(1);
L_AD->SetMarkerStyle(20);

TGraph *KS_AD = new TGraph(46, &out_mKS[0], &out_mAD[0]);
KS_AD->SetTitle("A-D vs K-S; -KS; -AD");
KS_AD->SetMarkerSize(1);
KS_AD->SetMarkerStyle(20);

TCanvas *cor = new TCanvas("ciao", "ciao", 1200,350);
cor->Divide(3,1,0.01,0.01);
cor->cd(1);
L_KS->Draw("AP");
//double *nx = L_KS->GetX(); double *ny = L_KS->GetY();
//for(int nn = 0; nn<38; nn++){
//	TMarker *m = new TMarker(nx[nn], ny[nn],20);
//	m->SetMarkerSize(1);m->SetMarkerColor(kBlue+nn); m->Draw();}
cor->cd(2);
L_AD->Draw("AP");
cor->cd(3);
KS_AD->Draw("AP");
cor->SaveAs("tests.pdf");




std::sort(rank.begin(), rank.end(), ordinamento_ander);
for(int i=0; i<46; i++) {rank[i].AD_index=i+1;}
std::sort(rank.begin(), rank.end(), ordinamento_kolgo);
for(int i=0; i<46; i++) {rank[i].KS_index=i+1;}
std::sort(rank.begin(), rank.end(), ordinamento_like);

string outname = "Results.txt";
ofstream outtxt;
outtxt.open(outname.c_str(), ios::out | ios::app);
for(int i=0; i<46; i++)
{
if (i==0) outtxt << "sorted by likelihood difference \n rank \t like. m. \t KS \t KS_rank \t AD \t AD_rank \t label" << endl;
float score = (i+1+rank[i].KS_index+rank[i].AD_index)/3.0;
outtxt << i+1 << "\t" << rank[i].like << "\t" <<  rank[i].kolgo << "\t" << rank[i].KS_index << "\t" << rank[i].ander 
<< "\t" << rank[i].AD_index << "\t" << rank[i].label << " \t \t "<< score << endl;
}
outtxt.close();


return 0;
}



double MyAndersonDarlingTest (const TH1 *h1, const TH1 *h2) {

  /*
  histograms must have the same binning,
  but not necessarily the same number of events
  */
  int L = h1->GetXaxis()->gsbins();
  double n1 = 0, n2 = 0, N;

  n1 = h1->GetEntries();
  n2 = h2->GetEntries();
  N = n1 + n2;
  if (n1*n2==0) {
    cout << "AD test: h1 or h2 integral is zero" << endl;
    return 0;
  }
  // Error("MyAndersonDarlingTest","Histogram2 %s integral is zero\n",h2->GetName());

  double f1=0, f2=0, l;
  double M1=0, M2=0, B=0, A1=0, A2=0, A=0;

  for (int j=1; j<=L && B<N; j++) { // NNBB as B=N and l=0 we'd get into infinities
    f1 = h1->GetBinContent(j);
    f2 = h2->GetBinContent(j);
    l  = f1 + f2;
    
    M1 += f1/2.0;
    M2 += f2/2.0;
    B  +=  l/2.0;

    if (B>0) { // if B=0 also l is zero and we divide by zero!      
      A1 += l * TMath::Power(N*M1 - n1*B, 2) / (B*(N-B) - N*l/4.0);
      A2 += l * TMath::Power(N*M2 - n2*B, 2) / (B*(N-B) - N*l/4.0);
      M1 += f1/2.0;
      M2 += f2/2.0;
      B  +=  l/2.0;
    }
  }
  A = (N-1)/TMath::Power(N,2) * (A1/n1 + A2/n2);
  return A;
}

vector<string> setTitle(int nh)
{
vector<string> title(3); title[2] = "counts" ;
if(nh==1)  	{title[1] = "pt (GeV)"	;title[0] = "Min pt among discarded jets"; }
if(nh==2)  { title[1] = "pt (GeV)"	;title[0] = "Average pt among discarded jets";}
if(nh==3)  { title[1] = "pt (GeV)"	;title[0] = "Max pt among discarded jets" ;}
if(nh==4)  { title[1] = "#eta"	;title[0] =	"Min #eta among discarded jets";}
if(nh==5)  { title[1] = "#eta"	;title[0] =	"Average #eta among discarded jets";}
if(nh==6)  { title[1] = "#eta"	;title[0] =	"Max #eta among discarded jets" ;}
if(nh==7)  { title[1] = "CMVA"	;title[0] =	"Min CMVA among discarded jets";}
if(nh==8)  { title[1] = "CMVA"	;title[0] =	"Average CMVA among discarded jets";}
if(nh==9)  { title[1] = "CMVA"	;title[0] =	"Max CMVA among discarded jets" ;}
if(nh==10) { title[1] = "centrality"	;title[0] = 	"Centrality (discarded jets)";}
if(nh==11) { title[1] = "centrality"	;title[0] = 	"Centrality (selected jets)";}
if(nh==12) { title[1] = "pt (GeV)"	;title[0] = 	"Jet 1 pt";}
if(nh==13) { title[1] = "pt (GeV)"	;title[0] = 	"Jet 2 pt";}
if(nh==14) { title[1] = "pt (GeV)"	;title[0] = 	"Jet 3 pt";}
if(nh==15) { title[1] = "pt (GeV)"	;title[0] = 	"Jet 4 pt";}
if(nh==16) { title[1] = "#eta"	;title[0] = "Jet 1 #eta";}
if(nh==17) { title[1] = "#eta"	;title[0] = "Jet 2 #eta";}
if(nh==18) { title[1] = "#eta"	;title[0] = "Jet 3 #eta";}
if(nh==19) { title[1] = "#eta"	;title[0] = "Jet 4 #eta";}
if(nh==20) { title[1] = "CMVA"	;title[0] = "Jet 1 CMVA";}
if(nh==21) { title[1] = "CMVA"	;title[0] = "Jet 2 CMVA";}
if(nh==22) { title[1] = "CMVA"	;title[0] = "Jet 3 CMVA";}
if(nh==23) { title[1] = "CMVA"	;title[0] = "Jet 4 CMVA";}
if(nh==24) { title[1] = "mass (GeV)"	;title[0] = "Dijet 1 mass";}
if(nh==25) { title[1] = "pt (GeV)"	;title[0] = "Dijet 1 pt";}
if(nh==26) { title[1] = "#Delta#phi (radians)"	;title[0] =	"Dijet 1 #Delta#phi";}
if(nh==27) { title[1] = "#Delta#eta"	;title[0] = 	"Dijet 1 #Delta#eta";}
if(nh==28) { title[1] = "#Delta R"	;title[0] = 	"Dijet 1 #DeltaR";}
if(nh==29) { title[1] = "mass (GeV)"	;title[0] = 	"Dijet 2 mass";}
if(nh==30) { title[1] = "pt (GeV)"	;title[0] = 	"Dijet 2 pt";}
if(nh==31) { title[1] = "#Delta#phi"	;title[0] = 	"Dijet 2 #Delta#phi";}
if(nh==32) { title[1] = "#Delta#eta"	;title[0] = 	"Dijet 2 #Delta#eta";}
if(nh==33) { title[1] = "#DeltaR"	;title[0] = 	"Dijet 2 #DeltaR";}
if(nh==34) {	title[1] = "pt (GeV)"	;title[0] = 	"Dijets pt vector sum";}
if(nh==35) {	title[1] = "#Delta#phi (radians)"	;title[0] = 	"#Delta#phi between dijets";}
if(nh==36) { title[1] = "#Delta#eta"	;title[0] = 	"#Delta#eta between dijets";}
if(nh==37) { title[1] = "#Delta R"	;title[0] = 	"#DeltaR between dijets";}
if(nh==38) { title[1] = "invariant mass (GeV)"	;title[0] = 	"hh invariant mass";}
if(nh==39) { title[1] = "MET (GeV)"	;title[0] = 	"MET";}
if(nh==40) { title[1] = "CMVA"	;title[0] = 	"minimum CMVA of first 3 selected jets";}
if(nh==41) { title[1] = "CMVA"	;title[0] = 	"average CMVA of first 3 selected jets";}
if(nh==42) { title[1] = "|cos(#theta *)|"	;title[0] = 	"|cos(#theta *)|";}
if(nh==43) { title[1] = "|cos(#theta CS)|"	;title[0] = 	"|cos(#theta CS)|";}
if(nh==44) { title[1] = "#tau"	;title[0] = 	"#tau 1";}
if(nh==45) { title[1] = "#tau"	;title[0] = 	"#tau 2";}
if(nh==46) { title[1] = "Jets Number"	;title[0] = 	"Jets Number";}
return title;
}




