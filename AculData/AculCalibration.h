#pragma once

#include <TObject.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1I.h>
#include <TPolyMarker.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TGraph.h>
#include <TObjArray.h>
#include <TRandom3.h>
#include <THStack.h>
#include <TString.h>
#include <TSpectrum.h>
#include <TROOT.h>
#include <TChain.h>

#include <iostream>
#include <fstream>
//#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

#include "../TELoss/TELoss.h"

//#define DEFAULTNOPEAKS 20
#define ADDRESSNUMBER 32

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;
using std::stringstream;
using std::ostringstream;
using std::ifstream;
using std::ofstream;
using std::vector;

struct AculTwoAnglePeakSet {
	TString measurementName;
	TString detectorName;
	Double_t angleDeg;
	vector<Double_t> alphaEnergies;
	vector<vector<Double_t> > stripPeaks;
	Int_t firstStrip;

	AculTwoAnglePeakSet() : angleDeg(0.), firstStrip(0) {}
};

struct AculTwoAngleStripResult {
	Int_t strip;
	Bool_t ok;
	Bool_t fallbackUsed;
	TString status;
	TString message;
	Double_t angleDeg;
	Double_t deadLayer;
	Double_t deadLayerObjectiveRms;
	Double_t temporaryA;
	Double_t temporaryB;
	Double_t finalA;
	Double_t finalB;
	Double_t finalRms;
	Double_t finalMaxResidual;
	vector<Double_t> normalChannels;
	vector<Double_t> angledChannels;
	vector<Double_t> normalEnergies;
	vector<Double_t> angledEnergies;
	vector<Double_t> finalResiduals;

	AculTwoAngleStripResult();
};

class AculCalibration : public TObject
{


public:

	//smysl jako verejne globalni promenne maji:
	//fNOSpectra	- pocet zkalibrovanych spekter
	//????			- pocet spravne zkalibrovanych spekter
	//????			- pocet nespravne zkalibrovanych spekter
	//fEnergy[4]	- tabulka s energiemi piku, nacita se zvenci

private:

//todo: following variables should be deleted
//	TObjArray	*fHRawList;			//list of raw histograms, list is set to owner
	TObjArray	fHRawList;			//list of raw histograms, list is set to owner
	TObjArray	fHAnalyzedList;	//list of fitted and analyzed histograms, list is set to owner
	TObjArray	fHEnergyList;		//list of calibrated histograms, list is set to owner
	THStack		fHEnergyStack;		//some stack
	THStack		*fCurrentHStack;
	//dodelat current histograms
	TObjArray	fCurrentHistList;

	//parameters to be read from file
	Int_t		kRaNOPEAKS;
	TArrayD		fEnergy;		//energy after passing through deadlayer
	TArrayD		fEnergyInput;		//incidental energy, set from .par 
	Bool_t		fQAlphaInput;		//true when fEnergyInput contains Q_alpha values
	Double_t	fFirstAlphaParentMass;	//mass number of the first parent in an alpha decay chain
	Double_t	fLowerChannel;
	Double_t	fUpperChannel;
	Double_t	fLowerPeakRelativeHight;	//pouziva se, private
	Double_t	fUpperPeakRelativeHight;	//pouziva se, private, nastavit nenulovou prednastavenou hodnotu
	Double_t	fPeakPositionTolerance;		//pouziva se, private
	Width_t		fFitFuncLineWidth;			//private
	Double_t	fFitMinSigma;				//pouziva se, private
	Double_t	fFitPeakThreshold;			//pouziva se, private, prozkoumat, k cemu vlastne slouzi ve fci ShowPeaks, popremyslet o vhodnem prednastaveni v konstruktoru
	Double_t	fDeadLayer;					//in mcm; it will be used for recalculation of energies outcomming from the alpha source
	Int_t    fuppersubaddress;
	Int_t    flowersubaddress;
	TString fblock;
	const char* block; //!
	TString ftreename;
	const char* treename; //!

	//these variables are the main for the class
	TArrayD fA;		//calibration parameter, f(x) = fA*x + fB
	TArrayD fB;		//calibration parameter, f(x) = fA*x + fB

	TArrayD fPeak;	//v teto promenne je ulozena momentalni hodnota piku v kanalech, zejmena v jedne fci, mozno udelat ji jako lokalni, bude navratovou hodnotou fce PeaksFitting, predelat delku pole

	TELoss fAlphaSi;

	TString fParFileName;
	TString fInputRootFile;
	TString fWorkDirectory;

	//it is very doubtful that we need this class variable
	TFile		*fCalInformation;

	

public:
	AculCalibration();
	//default constructor
	AculCalibration(const char* calfile);		//
	virtual ~AculCalibration();
	ClassDef(AculCalibration, 1);

	void Init();

	Double_t GetA(Int_t i); //to obtain calib parameter A
	Double_t GetB(Int_t i); //to obtain calib parameter B

	void SetParFileName(const char* parfile = "");
	// Function which loads text file containing parameters for calibration
	//
	//	-inputparfile: file containing information on calibration source
	//	-In file with the data must be preserved systematic
	//	-There cannot be whitelines (probably)
	//
	//Example of "parforcal.par" file
	//................................................................
	//223Ra
	//
	//4		nopeaks		//number of peaks
	//4.415		E1		//in MeV
	//5.153		E2		//in MeV
	//5.683		E3		//in MeV
	//7.419		E4		//in MeV
	//100		lowerchannel	//in channels
	//4096		upperchannel	//in channels
	//0.1		lowerpeakhight		//in relative units
	//0.1		upperpeakhight		//in relative units
	//0.1		peakpositiontolerance	//in relative units
	//2		fitfunctionlinewidth	//integer 1 - 10
	//5		minfitsigma		//minimal sigma of the peaks to be fitted
	//0.4		fithightthreshold	//
	//................................................................
	void SetInputRootFile(const char* rootfile = "");
	void SetWorkDirectory(const char* dir = "");

	Bool_t	SetCalibrationParameters(const char* calparfile);
	//Loads the file with calibration parameters
	//
	//  calparfile: file containing calibration parameters in format: crate number, address, subaddress, fAOld, fBOld, fC, fD
	//
	//	allowed comment characters: *, #, %, //



	void	PrintInputParameters();
	//I hope this is selfunderstanding function

	void	PrintCalibrationParameters();
	//Print calibration parameters fA, fB which are currently saved in object AculCalibration
 

	Bool_t	CalculateCalibParameters(Int_t lowerchannel = 0, Int_t upperchannel = 4095);
	//calculate calibration parameters for given block in given file

	//function is not completely ready to use
	//
	//function for calculation of calibrate parameters for DAQ system based on "Go4"

	//  lowerchannel: minimal channel from which the spectrum will be analysed
	//  upperchannel: maximal channel up to which the spectrum will be analysed


	Bool_t	CalculateCalibParameters5points(const char* inputfile, const char* pedestalfile, const char* block,
			const char* treename, Int_t lowerchannel = 0,
			Int_t upperchannel = 4095, Int_t nEBins = 1000, Int_t lowersubaddress = 0,
			Int_t uppersubaddress = ADDRESSNUMBER-1); //calculate calibration parameters for 5 energy points (4 peaks plus pedestal)
	//function is copy of CalculateCalibParameters. File (.par) with pedestal positions in channels should be provided. Energy of the additional point is set as 0 MeV by default


	Bool_t	Mycalc(	Int_t lowerchannel = 0, Int_t upperchannel = 4095); //calculate calibration parameters for 5 energy points (4 peaks plus pedestal)
	//function is copy of CalculateCalibParameters. File (.par) with pedestal positions in channels should be provided. Energy of the additional point is set as 0 MeV by default
	

	void CalibrateRawSpectra();

	UShort_t SomeFunction(UShort_t x, int len);

	template<typename T> T SomeFunction(T x, int len)
	{
		/*T max = x[0];
		for(int i = 1; i < len; i++)
			if(max < x[i])
				max = x[i];*/

//		return max;
		return x;
	}



	void CalibrateRawSpectra(const char* inputfile, const char* block, const char* treename, Int_t lowerchannel = 0, Int_t upperchannel = 4095, Int_t nEBins = 1000, Int_t lowersubaddress = 0, Int_t uppersubaddress = ADDRESSNUMBER-1);

	Bool_t ShowFullCalibratedSpectra(Int_t lowerchannel = 0, Int_t upperchannel = 4095);
	//Outputs 16 or 32 histograms in *full_calibrated_spectra.root containing full calibrated spectra 
	//from 0 to 10 MeV. Calibration parameters are in calibparfile. 
	//Calibparfile should be generated by function CalculateCalibParameters

	void FindEnergyPeaks(TCanvas *c1, const char* ifile, const char* outfile); 
	// Outputs calibrated energies for each channel in txt file 
 

	void FindAverageEnergies(const char* ifile, const char* outfile);
	// Outputs average values of calibrated energies for the whole detector in txt file


	Int_t	PeaksFitting(TH1* hSpectrum, Option_t* option = "", Double_t sigmamin = 2);
	//return values:
	//	0: OK
	//	1: something wrong
	//	2: something other wrong

	//possible options: "V", "Q", "", ?"gp"?

	Int_t	SearchPeaks(const TH1 *hin, Double_t sigma = 2, Option_t *option = "", const Int_t searchedpeaks = 100);

	void	FillRawSpectraFile(const char* rawdatafile, const char* block, const char* treename, TCanvas* rawCanvas = NULL, Option_t *option = "", Int_t xaxismin = 0, Int_t xaxismax = 4096);

	void	ShowRawSpectra(const char* filename, const Int_t block, TCanvas* rawCanvas = NULL, Int_t xaxismin = 0, Int_t xaxismax = 4096, /*TObjArray* histList = NULL,*/ const Int_t subaddress = 16);
	//this function is written for old format of raw data
	//do not use this function!!!!!!!!!!!
	void	ShowSpectra(const char* filename, TCanvas* rawCanvas = NULL, Option_t *option = "", Int_t xaxismin = 0, Int_t xaxismax = 4096, /*TObjArray* histList = NULL,*/ const Int_t subaddress = ADDRESSNUMBER);
	void	ShowAnalyzedSpectra(const char* filename, TCanvas* fittedRawCanvas = NULL, Int_t xaxismin = 0, Int_t xaxismax = 4096, Int_t subaddress = 16);
	//This function displays analyzed spectrum from a file, divides the canvas into a sufficient number of pads and displays
	//spectrums of each  block subadress on the suitable pads or displays one  selected spectrum .
	//Selects the peaks in the histogram and displays on the histogram how the spectrum were fitted.
	//
	//  filename: file .root containing analysed spectra
	//  fittedRawCanvas: canvas on which one you will see the spectrum
	//  xaxismin: Minimum channel, which will be displayed
	//  xaxismax: Maximum channel, which will be displayed
	//  subaddress:

	void	ShowEnergySpectra(const char *filename, TCanvas* energyCanvas = NULL, const Int_t subaddress = 16, Option_t* option = "", Double_t xaxismin = 0., Double_t xaxismax = 10.); //option: "sum", "c", "+",
	//Displays the spectrum  of the selected subbaddress block in MeV
	//
	//  filename: file .root containing calibrated spectra in MeV
	//  energyCanvas: : canvas on which one you will see the spectrum
	//  subaddress: block subaddress  which will be drawn
	//  option: sum ,+ ,c
	//  xaxismin: Minimum channel, which will be displayed
	//  xaxismax: Maximum channel, which will be displayed

	Bool_t FindPedestals(Int_t lowerchannel, Int_t upperchannel);
//	(const char* inputfile, Int_t lowerchannel, Int_t upperchannel);
	void FindEnergyPedestals();
//	void SetDeadlayer(const Double_t dl);

	AculTwoAngleStripResult AnalyzeTwoAngleStrip(
			const vector<Double_t>& normalChannels,
			const vector<Double_t>& angledChannels,
			const vector<Double_t>& alphaEnergies,
			Double_t angleDeg,
			Int_t stripIndex = -1,
			Double_t maxDeadLayer = 20.,
			Bool_t allowFallback = kFALSE);

	vector<AculTwoAngleStripResult> AnalyzeTwoAnglePreparedPeaks(
			const AculTwoAnglePeakSet& normalMeasurement,
			const AculTwoAnglePeakSet& angledMeasurement,
			Double_t maxDeadLayer = 20.,
			Bool_t allowFallback = kFALSE);

	vector<AculTwoAngleStripResult> AnalyzeTwoAnglePeakFiles(
			const char* normalPeakFile,
			const char* angledPeakFile,
			const char* detectorName,
			Double_t angleDeg,
			const vector<Double_t>& alphaEnergies,
			const char* textReport = "",
			const char* csvReport = "",
			Double_t maxDeadLayer = 20.,
			Bool_t allowFallback = kFALSE);

	Bool_t WriteTwoAngleDiagnostics(
			const char* textFile,
			const char* csvFile,
			const vector<AculTwoAngleStripResult>& results,
			const AculTwoAnglePeakSet& normalMeasurement,
			const AculTwoAnglePeakSet& angledMeasurement);

	Bool_t WriteTwoAngleDiagnosticGraphs(
			const char* rootFile,
			const vector<AculTwoAngleStripResult>& results);

	//dodelat funkce TTree* Get...(...) pro raw, anal i E spectra

	void	ClearHistograms(Option_t* option = "");
	void	DeleteStacks(Option_t* option = "");

	void	Reset();
	//reset calibration parameters fA, fB to zero

private:

	void DivideCanvas(TCanvas *c1, Int_t inputs);

	void SetInputParameters();

	void SetELosses();
	void SetCalEnergies();
	Double_t GetAlphaEnergyFromQAlpha(Double_t qAlpha, Int_t peakIndex) const;
	Double_t GetDetectedAlphaEnergy(Double_t alphaEnergy, Double_t deadLayer);
	Bool_t LoadTwoAnglePeakFile(const char* peakFile, AculTwoAnglePeakSet& measurement) const;


};
