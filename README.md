This repository collects all the source code written for my master thesis entitled "Search for Anomalous Production of Higgs Boson Pairs with the CMS Detector".

Each folder contains the code for a single task, in order:

1. events selection \\
2. variables ranking
3. TMVA deployment 
4. figure of merit calculation
5. ABCD counting experiment
6. shape fit
7. combine datacards and commands

a Mathematica notebook containing the 8 TeV pp->hh cross section parametrization as a function of the anomalous couplings is included.


---------------------------------------------------------
1. EVENTS SELECTION "event_selection/"
---------------------------------------------------------
* description: program to select wanted events and write new n-tuples containing only useful info.
* input: reco-level ntuples
* ouput:- >=3 btags events ntuple (hhev_maker.cc)
	- 2 btags, 1 btag<loosecut events ntuple (hhev_control_region.cc)		
	- plots of 3° jet pt, eta, no. of tracks ratio between the two categories above, check for the ABCD method (plottini.cc)


---------------------------------------------------------
2. VARIABLES RANKING "var_ranking/Minuit_carlo.cpp"
---------------------------------------------------------
* description: ranking of the kinematical variables according to Kolmogorov-Smirnov, Anderson-Darling and Signal-Fraction-Fit. 
* input: reco-level ntuples, one ntuple for background, one for signal, and an optional TFile containing histograms of the variables just to set the titles and ranges for the output histograms.
* output: ascii file containing output values and ranking from each of the three methods, a variable label and a ranking average.

 
---------------------------------------------------------
3. TMVA "tmva/TMVAtrain_app_tfile.cpp"
---------------------------------------------------------
* description: training, testing and application of TMVA methods "Boosted Decision Tree" and "Projective Likelihood".
*input: 1 signal, 1 bkg ntuples, ntuples on which apply the trained tool
*output: - TMVA generated files: weights folder, TMVA_testing.root containing correlations etc.
	 - "TMVAallTraining.root" containig BDT and Likelihood response distribution of the testing phase and of all the applications.
	 - copy of each ntuple on which the BDT is applied that includes a new branch called "BDT" containing the BDT response.
 
---------------------------------------------------------
4. FIGURE OF MERIT CALCULATION "figure_of_merit/Qcalc.cpp"
---------------------------------------------------------
* description: determination of the optimal cut value on the BDT response
* input: bkg and signal ntuples
* output: pdf and .root file containing the plots of the figure of merit 2(sqrt(s+b) - sqrt(b)) and the likelihood-ratio-based one.


---------------------------------------------------------
5) ABCD COUNTING EXPERIMENT "abcd/ABCDcontrol.cpp"
---------------------------------------------------------
* description: 9-region ABCD method parametrized on 3° jet variables (see thesis). It tests 16 different validation cuts. 
* input: data ntuples (2, 3 and 4 btags events)
* output: - on screen: events and expected background in signal region, other related info
	  - ValidationCut.root containing the source if systematic error Delta as a function of the validation cut.

---------------------------------------------------------
6) SHAPE FIT "fit/2-3CSVGOF.cpp" and "fit/histo_to_Combine.cpp"
---------------------------------------------------------
2-3CSVGOF.cpp
* description: creation of three bidimensional shapes, signal template, background template and observed data respectively.
* input: selected events ntuples
* output: "scatters.root" containing the three TH2D objects plus a TH2D for the 2-3btags residual.

histo_to_Combine.cpp
* description: conversion of the 2D histograms into 1D ones with no empty bins, creation of the +1sigma and -1sigma uncertainty shapes.
* input: "scatters.root"
* output: "toCombine.root"
---------------------------------------------------------
7) COMBINE
---------------------------------------------------------
* description: datacard for the counting experiment and the shape fit. Commands in "ToRunCombine.txt"


