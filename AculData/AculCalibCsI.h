#pragma once

//#include "TObject.h"
//#include "TROOT.h"
#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TTree.h>
#include <TCutG.h>
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TH1I.h"
#include <TGraphErrors.h>
//#include "TGraph"
#include "TArrayD.h"
#include "TF1.h"

#define NOCALFILES 5

using std::cout;
using std::endl;

class AculCalibCsI {

private:

	TString detName;
	TString partName;
	
	TClonesArray fr;		//TFile
//	TFile *fData;
	TClonesArray tr;
	TClonesArray colFiles;
	TClonesArray colTrees;
	Int_t nofiles;
	TString fileNames[100];
	TString cutNames[100];

	Int_t energyPoints;
	TArrayD fE;

	TFile *fCuts;
	TString cutsFileName;
	Int_t nCuts;				//number of cuts
	TClonesArray cutsCol;

	TH1I *hfull[NOCALFILES][16];	//!
	TH1I *hcut[NOCALFILES][16];		//!

	Int_t peakMin[NOCALFILES][16];
	Int_t peakMax[NOCALFILES][16];
	Double_t mean[NOCALFILES][16];
	Double_t meanRMS[NOCALFILES][16];

	TGraphErrors *gCal[16];
	TFile *fGraphs;
	TArrayD fA;
	TArrayD fB;

	TString fParFileName;

public:
//	AculCalibCsI() : a(0), b(0), c(0), p(0){};
	AculCalibCsI();
	AculCalibCsI(const char* parfile);
	virtual ~AculCalibCsI();
	// Define the class for the cint dictionary
	ClassDef (AculCalibCsI,1);


	void OpenTrees();
	void LoadCuts();

	void SetParFile(const char* parfile);

	void PrintParameters(const char* option = "");
	void PrintPeakRanges();

	void DrawVariable(const char* variable, Int_t tree, TCanvas *canvas, Int_t lowRange = 0, Int_t upRange = 4096);
	void DrawBeam(TCanvas *canvas, Int_t files, const char* variable);
//	void DrawdEE(const char* variable, Int_t tree, TCanvas *canvas);
	void DrawVariableCut(const char* variable, Int_t tree, TCanvas *canvas, const char* cut1, const char* cut2 = "", Int_t lowRange = 0);
	void GetPeakMean(const char* variable, Int_t tree, Int_t energy, TCanvas *canvas, const char* beamcut, const Int_t nbins = 4096, Int_t lowRange = 0);

	void Calibrate(TCanvas *canvas, Bool_t savefile = 0, const char* filename = "", const char* option = "READ");
	void FillGraph(TGraphErrors *g, Int_t npoints, Double_t *energies, Int_t graphNumber, const char* option = "");
	void WriteClbParameters(const char* filename);
	void SaveClbGraphs(const char* filename, const char* option = "READ");
	void ReadClbParameters(const char* filename);
	void DrawClbGraphs(const char* filename, const char* graphname, TCanvas *canvas);
	void DrawEnergyDeposite(const char* variable, TCanvas *canvas, Int_t tree, const char* option = "");

	void PrintTrees();
	void PrintFiles();
	void PrintCuts();

	Double_t GetA(Int_t i);
	Double_t GetB(Int_t i);

private:

	void SetPars();
};
