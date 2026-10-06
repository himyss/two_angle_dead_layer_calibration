//////////////////////////////////////////////////////////
//									//
//	AculCalibration					//
//									//
//
//Some description of this very useful class,
//its properties and a short example how to
//use it. This text may be partly used in
//PhD thesis.
//
//Description of the detector itself.
//
//														//
//////////////////////////////////////////////////////////

#include "AculCalibration.h"

#include <algorithm>
#include <cmath>
#include <limits>

//ClassImp(AculCalibration);

namespace {

struct AculLineFit {
	Bool_t ok;
	Double_t a;
	Double_t b;
	Double_t rms;
	Double_t maxResidual;
	vector<Double_t> residuals;

	AculLineFit() :
		ok(kFALSE),
		a(0.),
		b(0.),
		rms(std::numeric_limits<Double_t>::infinity()),
		maxResidual(std::numeric_limits<Double_t>::infinity())
	{}
};

Bool_t IsGoodNumber(Double_t value)
{
	return std::isfinite(value);
}

AculLineFit FitLine(const vector<Double_t>& channels, const vector<Double_t>& energies)
{
	AculLineFit fit;
	const size_t n = channels.size();
	if (n != energies.size() || n < 2) {
		return fit;
	}

	Double_t sumX = 0.;
	Double_t sumY = 0.;
	Double_t sumXX = 0.;
	Double_t sumXY = 0.;
	for (size_t i = 0; i < n; i++) {
		if (!IsGoodNumber(channels[i]) || !IsGoodNumber(energies[i])) {
			return fit;
		}
		sumX += channels[i];
		sumY += energies[i];
		sumXX += channels[i]*channels[i];
		sumXY += channels[i]*energies[i];
	}

	const Double_t denom = static_cast<Double_t>(n)*sumXX - sumX*sumX;
	if (TMath::Abs(denom) < 1.e-12) {
		return fit;
	}

	fit.a = (static_cast<Double_t>(n)*sumXY - sumX*sumY) / denom;
	fit.b = (sumY - fit.a*sumX) / static_cast<Double_t>(n);
	fit.residuals.resize(n);

	Double_t sumResidual2 = 0.;
	fit.maxResidual = 0.;
	for (size_t i = 0; i < n; i++) {
		const Double_t predicted = fit.a*channels[i] + fit.b;
		fit.residuals[i] = energies[i] - predicted;
		sumResidual2 += fit.residuals[i]*fit.residuals[i];
		fit.maxResidual = TMath::Max(fit.maxResidual, TMath::Abs(fit.residuals[i]));
	}
	fit.rms = TMath::Sqrt(sumResidual2 / static_cast<Double_t>(n));
	fit.ok = kTRUE;
	return fit;
}

TString FormatVector(const vector<Double_t>& values)
{
	TString text;
	for (size_t i = 0; i < values.size(); i++) {
		if (i) {
			text += ";";
		}
		text += TString::Format("%.10g", values[i]);
	}
	return text;
}

} // namespace

AculCalibration::AculCalibration() : fEnergy(0), fEnergyInput(0), fA(0), fB(0), fPeak(0)
{
	//default constructor

	fCurrentHStack = NULL;
	fCurrentHistList.IsOwner();
//	todo: change size of fA and fB in some other place
	fA.Set(32);
	fB.Set(32);
//	fEnergy.Set(4);
//	fEnergyInput.Set(4);

	kRaNOPEAKS = 0;
	fLowerPeakRelativeHight = 0.;
	fUpperPeakRelativeHight = 0.;
	fPeakPositionTolerance = 0.;
	fFitFuncLineWidth = 1;
	fFitMinSigma = 0.;
	fFitPeakThreshold = 0.;
	fuppersubaddress = 0.;
	flowersubaddress = 0.;
//	fblock = 'a';

	fDeadLayer = 0.;
	fQAlphaInput = kFALSE;
	fFirstAlphaParentMass = 226.;

	/*for(Int_t i = 0; i < DEFAULTNOPEAKS; i++) {
//		fEnergy[i] = 0.;
//		fEnergyInput[i] = 0.;
//		fPeak[i] = 0.;
	}*/
	fCalInformation = 0;

	Reset();

}

AculCalibration::AculCalibration(const char* calfile) : fEnergy(0), fEnergyInput(0), fPeak(0)
{
	//constructor which fills fAOld, fBOld, fC, fD from file parfile

	fCurrentHStack = NULL;
	fCurrentHistList.IsOwner();

	kRaNOPEAKS = 0;
	fLowerPeakRelativeHight = 0.;
	fUpperPeakRelativeHight = 0.;
	fPeakPositionTolerance = 0.;
	fFitFuncLineWidth = 1;
	fFitMinSigma = 0.;
	fFitPeakThreshold = 0.;
	fuppersubaddress = 0.;
	flowersubaddress = 0.;
//	fblock = 'a';	
	fDeadLayer = 0.;
	fQAlphaInput = kFALSE;
	fFirstAlphaParentMass = 226.;

	/*for(Int_t i = 0; i < DEFAULTNOPEAKS; i++) {
//		fEnergy[i] = 0.;
//		fEnergyInput[i] = 0.;
//		fPeak[i] = 0.;
	}*/

	fCalInformation = 0;

	SetCalibrationParameters(calfile);

}

void AculCalibration::Reset()
{
	for (Int_t j = 0; j < fA.GetSize(); j++) {
			fA[j] = 0;
			fB[j] = 0;
	}

	return;
}

AculCalibration::~AculCalibration()
{

	DeleteStacks();
//	delete fCalInformation;
//	fCalInformation->Close();

}

void AculCalibration::Init() {

	SetELosses();
	SetInputParameters();
	SetCalEnergies();
}

Int_t AculCalibration::SearchPeaks(const TH1 *hin, Double_t sigma, Option_t *option, const Int_t searchedpeaks)
{
	//Function searching peaks in inputed TH1 spectrum and selects the peaks in the histogram.
	//
	//  hin:
	//  sigma:
	//  option:
	//  threshold:
	//  searchedpeaks:

	TSpectrum sc;	//by default for 100 peaks
	Int_t nopeaks = sc.Search(hin, sigma, "goff", fFitPeakThreshold);

	TString opt = option;
	opt.ToLower();

	const Double_t tStep = 0.05;

	while ( nopeaks > searchedpeaks && fFitPeakThreshold <= 1) {
		fFitPeakThreshold = fFitPeakThreshold + tStep;
		nopeaks = sc.Search(hin, sigma, "goff", fFitPeakThreshold);
	}

	if (!nopeaks) {
		return 0;
	}

	if (opt.Contains("goff")) {
		return nopeaks;
	}

	TPolyMarker *pm = (TPolyMarker*)hin->GetListOfFunctions()->FindObject("TPolyMarker");
	if (pm) {
		hin->GetListOfFunctions()->Remove(pm);
		delete pm;
	}
	pm = new TPolyMarker(nopeaks, sc.GetPositionX(), sc.GetPositionY());
	hin->GetListOfFunctions()->Add(pm);
	pm->SetMarkerStyle(23);
	pm->SetMarkerColor(kRed);
	pm->SetMarkerSize(1.3);

	return nopeaks;
}

Int_t AculCalibration::PeaksFitting(TH1* hSpectrum, Option_t* option, Double_t sigmamin)
{

	if (!hSpectrum) {
    cout<< "raw spec was not found " << endl;
    return 1;
  }
	Int_t dimension = hSpectrum->GetDimension();
	if (dimension > 1) {
		Error("PeaksFitting", "Only implemented for 1-d histograms");
		return 1;
	}

	TString	opt = option;
	opt.ToLower();

	if (!kRaNOPEAKS) {
		Error("PeaksFitting", "kRaNOPEAKS is set to zero; calibration spectrum must be set");
		return 1;
	}

	const Int_t peaksNumber =	SearchPeaks(hSpectrum, sigmamin, "", kRaNOPEAKS);

	if (peaksNumber != kRaNOPEAKS) {
		Info("PeaksFitting", "In histogram %s was found %d peaks", hSpectrum->GetName(), peaksNumber);
		return 1;
	}
	//should be optional output
	Info("PeaksFitting", "Number of peaks in %s: %d", hSpectrum->GetName(), peaksNumber);

	//working array for peaks, there are founded in accidental order
	Double_t peak[peaksNumber];
	Double_t *peakPosition;
	Double_t *peakHight;

	TList *functions = hSpectrum->GetListOfFunctions();
	TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");

	peakPosition = pm->GetX();
	peakHight = pm->GetY();

	for (Int_t i = 0; i < peaksNumber; i++) {

		Double_t fitMin = 0;
		Double_t fitMax = 0;
		Double_t fitStep = hSpectrum->GetXaxis()->GetBinWidth(0);
//		cout << fitStep << endl;
//		cout << fLowerPeakRelativeHight << "\t" << fUpperPeakRelativeHight << endl;

		//fitting range:
		//shift a range of fit and search for raw boarder of peak determined by fUpperPeakRelativeHight
		//maximum
		Int_t j = 0;
		Double_t currentHight = peakHight[i];
		while ( currentHight > (peakHight[i]*fUpperPeakRelativeHight) ) {
			j++;
			fitMax = static_cast<Double_t>(peakPosition[i]) + j*fitStep;
			currentHight = hSpectrum->GetBinContent(hSpectrum->GetXaxis()->FindBin(fitMax));
		}

		//minimum
		j = 0;
		currentHight = peakHight[i];
		while ( currentHight > (peakHight[i]*fLowerPeakRelativeHight) ) {
			j++;
			fitMin = static_cast<Double_t>(peakPosition[i]) - j*fitStep;
			currentHight = hSpectrum->GetBinContent(hSpectrum->GetXaxis()->FindBin(fitMin));
		}

		//fitting
		if (opt.Contains("gp")) {
			Info("PeaksFitting", "Option containing gp");
			char fncname[20];
			sprintf(fncname, "gaus_aux_%d", i);
			TF1 *gausAux = new TF1(fncname, "gaus", fitMin - 10, fitMax + 10);		//pomocny gaus
			hSpectrum->Fit(fncname, "0 Q", "", fitMin - 15, fitMax + 15);				//prvotni fitovani

			sprintf(fncname, "auto_gp_%d", i);
			TF1 *fitAuto = new TF1(fncname, "gaus(0) + pol0(3)", fitMin - 15, fitMax + 15);		//fce pro automaticke fitovani
			fitAuto->SetParameter(0, gausAux->GetParameter(0));		//nastavovani parametru fitovaci fce
			fitAuto->SetParameter(1, gausAux->GetParameter(1));
			fitAuto->SetParameter(2, gausAux->GetParameter(2));

			hSpectrum->Fit(fncname, "0 R Q +", "", fitMin - 15, fitMax + 15); //dodelat zapis vsech fci
			hSpectrum->GetFunction(fncname)->ResetBit(TF1::kNotDraw);
			peak[i] = fitAuto->GetParameter(1);			//zapis asi pozice v kanalech do pomocneho pole
			if (opt.Contains("V")) {
				Info("PeaksFitting", "Peak position is\t %4.2f \tresolution is \t %2.1f %%", fitAuto->GetParameter(1), 235*(fitAuto->GetParameter(2))/(fitAuto->GetParameter(1)));
			}
		}
		else {
			char fncname[20];
			sprintf(fncname, "auto_g%d", i);
			TF1 *fitAuto = new TF1(fncname, "gaus", fitMin, fitMax);		//fce pro automaticke fitovani
//			cout << fitMin << "\t" << fitMax << endl;
//			fitAuto->SetParameter(2, fitMax-fitMin);
			fitAuto->SetLineWidth(fFitFuncLineWidth);
			hSpectrum->Fit(fncname, "+ 0 R Q", ""/*, fitMin - 1, fitMax + 1*/);
//			hSpectrum->GetFunction(fncname)->ResetBit(TF1::kNotDraw);
			hSpectrum->GetFunction(fncname)->InvertBit(TF1::kNotDraw);
			peak[i] = fitAuto->GetParameter(1);			//zapis asi pozice v kanalech do pomocneho pole
			if (opt.Contains("v")) {
				Info("PeaksFitting", "Peak position is\t%4.2f\tresolution is \t%2.1f %%", fitAuto->GetParameter(1), 235*(fitAuto->GetParameter(2))/(fitAuto->GetParameter(1)));
			}
		}//else
		//end of fitting
	}//for over all analyzed peaks

	//peaks sorting
	Int_t j[peaksNumber];
	TMath::Sort(peaksNumber, peak, j, kFALSE);
	fPeak.Set(peaksNumber);
	for (Int_t i = 0; i < peaksNumber; i++) {
		fPeak[i] = peak[j[i]];
		//printf("\tPeak peak\t%f\n", fPeak[i]);
	}

	if (!opt.Contains("q") || opt.Contains("v")) {
		Info("PeaksFitting", "Control output:");
		for (Int_t i = 0; i < peaksNumber; i++) {
			printf("\tPeak position is\t%f\n", fPeak[i]);
		}
	}

	//	provest kontrolu pomerne polohy piku,
	//	jestli jsou spatne, provest urcita opatreni,
	//	napr. zapis daneho histogramu do souboru,
	//	zapis do souboru s chybama, vypis na obrazovku, ...
	for (Int_t i = 0; i < peaksNumber; i++) {
			if ( !( (((1-fPeakPositionTolerance)*(fEnergy[0]/fEnergy[i])) < (fPeak[0]/fPeak[i])) && (((1+fPeakPositionTolerance)*(fEnergy[0]/fEnergy[i])) > (fPeak[0]/fPeak[i])) ) ) {
			//printf("\tPeaksFitt fEnergy\t%f\n", fEnergy[i]);		
			if (fCalInformation /* && opt.Contains("writebad")*/) {
				fCalInformation->cd();
				hSpectrum->Write();
			}
			//return 2;*/
		}
	}//for

	return 0;
}


//Bool_t AculCalibration::EnergyPositions(const char* inputfile, const char* block,
//			const Int_t address, const char* treename, Int_t lowerchannel,
//			Int_t upperchannel, Int_t lowersubaddress, Int_t uppersubaddress)
//{
//	TString iFile = inputfile;
//	TFile fr(iFile.Data());
//
//
//
//	return 1;
//}

void AculCalibration::SetParFileName(const char *parfile) {
	fParFileName = parfile;
}
void AculCalibration::SetInputRootFile(const char *rootfile) {
	fInputRootFile = rootfile;
}
void AculCalibration::SetWorkDirectory(const char* dir ){
	fWorkDirectory = dir;
}
void AculCalibration::SetELosses() {

	Info("AculCalibration::SetELosses", "Combination of aplha particle with silicon material only.");
	fAlphaSi.SetEL(1, 2.321); // density in g/cm3
	fAlphaSi.AddEL(14., 28.086, 1);  //Z, mass
//	mSi.SetZP(1., 1.);		//protons
	fAlphaSi.SetZP(2., 4.);		//alphas, Z, A
	fAlphaSi.SetEtab(100000, 200.);	// ?, MeV calculate ranges
	fAlphaSi.SetDeltaEtab(300);
}

void AculCalibration::SetCalEnergies() {

	TArrayD alphaEnergy(kRaNOPEAKS);
	for(Int_t i = 0; i < kRaNOPEAKS; i++) {
		alphaEnergy[i] = fQAlphaInput ? GetAlphaEnergyFromQAlpha(fEnergyInput[i], i) : fEnergyInput[i];
	}

	if (fDeadLayer<=0.) {
		Warning("AculCalibration::SetCalEnergies", "Dead layer was set equal or less than 0.");
		for(Int_t i = 0; i < kRaNOPEAKS; i++) {
			fEnergy[i] = alphaEnergy[i];
		}
		Info("AculCalibration::SetCalEnergies", "Energies used for calibration were set without dead layer correction.");
		return;
	}

	for(Int_t i = 0; i < kRaNOPEAKS; i++) {
		fEnergy[i] = fAlphaSi.GetE(alphaEnergy[i], fDeadLayer);
	}
	Info("AculCalibration::SetCalEnergies", "Energies used for calibration considering %f mcm dead layer were set.", fDeadLayer);

	return;
}

Double_t AculCalibration::GetAlphaEnergyFromQAlpha(Double_t qAlpha, Int_t peakIndex) const {
	Double_t parentMass = fFirstAlphaParentMass - 4. * peakIndex;
	if (parentMass <= 4.) {
		Warning("AculCalibration::GetAlphaEnergyFromQAlpha", "Invalid parent mass %.3f for peak %d; Q_alpha value is used directly.", parentMass, peakIndex);
		return qAlpha;
	}
	return qAlpha * (parentMass - 4.) / parentMass;
}

AculTwoAngleStripResult::AculTwoAngleStripResult() :
	strip(-1),
	ok(kFALSE),
	fallbackUsed(kFALSE),
	status("not_run"),
	message(""),
	angleDeg(0.),
	deadLayer(0.),
	deadLayerObjectiveRms(std::numeric_limits<Double_t>::infinity()),
	temporaryA(0.),
	temporaryB(0.),
	finalA(0.),
	finalB(0.),
	finalRms(std::numeric_limits<Double_t>::infinity()),
	finalMaxResidual(std::numeric_limits<Double_t>::infinity())
{}

Double_t AculCalibration::GetDetectedAlphaEnergy(Double_t alphaEnergy, Double_t deadLayer)
{
	if (deadLayer <= 0.) {
		return alphaEnergy;
	}
	return fAlphaSi.GetE(alphaEnergy, deadLayer);
}

AculTwoAngleStripResult AculCalibration::AnalyzeTwoAngleStrip(
		const vector<Double_t>& normalChannels,
		const vector<Double_t>& angledChannels,
		const vector<Double_t>& alphaEnergies,
		Double_t angleDeg,
		Int_t stripIndex,
		Double_t maxDeadLayer,
		Bool_t allowFallback)
{
	AculTwoAngleStripResult result;
	result.strip = stripIndex;
	result.angleDeg = angleDeg;
	result.normalChannels = normalChannels;
	result.angledChannels = angledChannels;

	const size_t nPeaks = alphaEnergies.size();
	if (nPeaks < 2) {
		result.status = "failed";
		result.message = "at least two alpha peaks are required";
		return result;
	}
	if (normalChannels.size() != nPeaks || angledChannels.size() != nPeaks) {
		result.status = "failed";
		result.message = "normal/angled peak counts do not match alpha energies";
		return result;
	}
	if (angleDeg <= 0. || angleDeg >= 89.) {
		result.status = "failed";
		result.message = "angle must be in (0, 89) degrees";
		return result;
	}
	if (maxDeadLayer <= 0.) {
		result.status = "failed";
		result.message = "maxDeadLayer must be positive";
		return result;
	}

	for (size_t i = 0; i < nPeaks; i++) {
		if (!IsGoodNumber(alphaEnergies[i]) || !IsGoodNumber(normalChannels[i]) || !IsGoodNumber(angledChannels[i])) {
			result.status = "failed";
			result.message = "non-finite input value";
			return result;
		}
	}

	SetELosses();
	const Double_t cosTheta = TMath::Cos(angleDeg*TMath::DegToRad());
	if (cosTheta <= 0.) {
		result.status = "failed";
		result.message = "non-positive cos(angle)";
		return result;
	}

	AculLineFit bestFit;
	Double_t bestDeadLayer = 0.;
	auto evaluate = [&](Double_t deadLayer, AculLineFit* fitOut = 0) -> Double_t {
		vector<Double_t> channels;
		vector<Double_t> energies;
		channels.reserve(2*nPeaks);
		energies.reserve(2*nPeaks);

		for (size_t i = 0; i < nPeaks; i++) {
			const Double_t normalEnergy = GetDetectedAlphaEnergy(alphaEnergies[i], deadLayer);
			const Double_t angledEnergy = GetDetectedAlphaEnergy(alphaEnergies[i], deadLayer/cosTheta);
			if (!IsGoodNumber(normalEnergy) || !IsGoodNumber(angledEnergy) || normalEnergy <= 0. || angledEnergy <= 0.) {
				return std::numeric_limits<Double_t>::infinity();
			}
			channels.push_back(normalChannels[i]);
			energies.push_back(normalEnergy);
			channels.push_back(angledChannels[i]);
			energies.push_back(angledEnergy);
		}

		AculLineFit fit = FitLine(channels, energies);
		if (!fit.ok || !IsGoodNumber(fit.rms)) {
			return std::numeric_limits<Double_t>::infinity();
		}
		if (fitOut) {
			*fitOut = fit;
		}
		return fit.rms;
	};

	const Double_t gr = (TMath::Sqrt(5.) - 1.)/2.;
	Double_t left = 0.;
	Double_t right = maxDeadLayer;
	Double_t x1 = right - gr*(right-left);
	Double_t x2 = left + gr*(right-left);
	Double_t f1 = evaluate(x1);
	Double_t f2 = evaluate(x2);

	for (Int_t iter = 0; iter < 80; iter++) {
		if (f1 > f2) {
			left = x1;
			x1 = x2;
			f1 = f2;
			x2 = left + gr*(right-left);
			f2 = evaluate(x2);
		}
		else {
			right = x2;
			x2 = x1;
			f2 = f1;
			x1 = right - gr*(right-left);
			f1 = evaluate(x1);
		}
		if (TMath::Abs(right-left) < 1.e-5) {
			break;
		}
	}

	bestDeadLayer = 0.5*(left+right);
	result.deadLayerObjectiveRms = evaluate(bestDeadLayer, &bestFit);
	result.deadLayer = bestDeadLayer;
	result.temporaryA = bestFit.a;
	result.temporaryB = bestFit.b;

	if (!bestFit.ok || !IsGoodNumber(result.deadLayerObjectiveRms) || bestDeadLayer <= 1.e-6 || bestDeadLayer >= maxDeadLayer - 1.e-6) {
		result.status = "failed";
		result.message = "dead-layer fit is non-physical or at search boundary";
		if (!allowFallback) {
			return result;
		}
		result.fallbackUsed = kTRUE;
		result.status = "fallback";
		result.message = "using zero dead layer because explicit fallback is enabled";
		bestDeadLayer = 0.;
		result.deadLayer = 0.;
	}

	vector<Double_t> normalEnergies;
	vector<Double_t> angledEnergies;
	normalEnergies.reserve(nPeaks);
	angledEnergies.reserve(nPeaks);
	for (size_t i = 0; i < nPeaks; i++) {
		normalEnergies.push_back(GetDetectedAlphaEnergy(alphaEnergies[i], bestDeadLayer));
		angledEnergies.push_back(GetDetectedAlphaEnergy(alphaEnergies[i], bestDeadLayer/cosTheta));
	}
	result.normalEnergies = normalEnergies;
	result.angledEnergies = angledEnergies;

	AculLineFit finalFit = FitLine(normalChannels, normalEnergies);
	if (!finalFit.ok || !IsGoodNumber(finalFit.rms)) {
		result.ok = kFALSE;
		result.status = "failed";
		result.message = "final normal-incidence calibration fit failed";
		return result;
	}

	result.finalA = finalFit.a;
	result.finalB = finalFit.b;
	result.finalRms = finalFit.rms;
	result.finalMaxResidual = finalFit.maxResidual;
	result.finalResiduals = finalFit.residuals;
	result.ok = !result.fallbackUsed;
	if (result.status != "fallback") {
		result.status = "ok";
		result.message = "two-angle dead-layer fit succeeded";
	}
	return result;
}

vector<AculTwoAngleStripResult> AculCalibration::AnalyzeTwoAnglePreparedPeaks(
		const AculTwoAnglePeakSet& normalMeasurement,
		const AculTwoAnglePeakSet& angledMeasurement,
		Double_t maxDeadLayer,
		Bool_t allowFallback)
{
	vector<AculTwoAngleStripResult> results;
	if (normalMeasurement.detectorName.Length() > 0 &&
			angledMeasurement.detectorName.Length() > 0 &&
			normalMeasurement.detectorName != angledMeasurement.detectorName) {
		AculTwoAngleStripResult failed;
		failed.status = "failed";
		failed.message = "detector names do not match";
		results.push_back(failed);
		return results;
	}

	const vector<Double_t>& alphaEnergies =
		normalMeasurement.alphaEnergies.empty() ? angledMeasurement.alphaEnergies : normalMeasurement.alphaEnergies;
	const size_t nStrips = TMath::Min(normalMeasurement.stripPeaks.size(), angledMeasurement.stripPeaks.size());
	results.reserve(nStrips);

	for (size_t i = 0; i < nStrips; i++) {
		const Int_t strip = normalMeasurement.firstStrip + static_cast<Int_t>(i);
		results.push_back(AnalyzeTwoAngleStrip(
					normalMeasurement.stripPeaks[i],
					angledMeasurement.stripPeaks[i],
					alphaEnergies,
					angledMeasurement.angleDeg,
					strip,
					maxDeadLayer,
					allowFallback));
	}

	if (normalMeasurement.stripPeaks.size() != angledMeasurement.stripPeaks.size()) {
		Warning("AculCalibration::AnalyzeTwoAnglePreparedPeaks",
				"Normal measurement has %zu strips, angled measurement has %zu strips; only %zu common strips were processed.",
				normalMeasurement.stripPeaks.size(), angledMeasurement.stripPeaks.size(), nStrips);
	}
	return results;
}

Bool_t AculCalibration::LoadTwoAnglePeakFile(const char* peakFile, AculTwoAnglePeakSet& measurement) const
{
	if (!peakFile || !strlen(peakFile)) {
		Error("AculCalibration::LoadTwoAnglePeakFile", "Empty peak file name.");
		return kFALSE;
	}

	ifstream input(peakFile);
	if (!input.is_open()) {
		Error("AculCalibration::LoadTwoAnglePeakFile", "Cannot open peak file %s.", peakFile);
		return kFALSE;
	}

	measurement.stripPeaks.clear();
	Int_t expectedStrip = -1;
	Bool_t firstDataLine = kTRUE;
	const Int_t lineLength = 4096;
	char line[lineLength];
	while (input.getline(line, lineLength)) {
		TString text(line);
		text = text.Strip(TString::kBoth);
		if (text.Length() == 0 || text.BeginsWith("#") || text.BeginsWith("//")) {
			continue;
		}

		stringstream parser(line);
		Int_t strip = -1;
		parser >> strip;
		if (!parser.good() && parser.fail()) {
			continue;
		}
		if (firstDataLine) {
			measurement.firstStrip = strip;
			expectedStrip = strip;
			firstDataLine = kFALSE;
		}
		while (expectedStrip < strip) {
			measurement.stripPeaks.push_back(vector<Double_t>());
			expectedStrip++;
		}

		vector<Double_t> peaks;
		Double_t value = 0.;
		while (parser >> value) {
			peaks.push_back(value);
		}
		measurement.stripPeaks.push_back(peaks);
		expectedStrip = strip + 1;
	}

	return !measurement.stripPeaks.empty();
}

vector<AculTwoAngleStripResult> AculCalibration::AnalyzeTwoAnglePeakFiles(
		const char* normalPeakFile,
		const char* angledPeakFile,
		const char* detectorName,
		Double_t angleDeg,
		const vector<Double_t>& alphaEnergies,
		const char* textReport,
		const char* csvReport,
		Double_t maxDeadLayer,
		Bool_t allowFallback)
{
	AculTwoAnglePeakSet normal;
	AculTwoAnglePeakSet angled;
	normal.measurementName = normalPeakFile;
	angled.measurementName = angledPeakFile;
	normal.detectorName = detectorName;
	angled.detectorName = detectorName;
	normal.angleDeg = 0.;
	angled.angleDeg = angleDeg;
	normal.alphaEnergies = alphaEnergies;
	angled.alphaEnergies = alphaEnergies;

	vector<AculTwoAngleStripResult> results;
	if (!LoadTwoAnglePeakFile(normalPeakFile, normal) || !LoadTwoAnglePeakFile(angledPeakFile, angled)) {
		AculTwoAngleStripResult failed;
		failed.status = "failed";
		failed.message = "cannot load prepared peak files";
		results.push_back(failed);
		return results;
	}

	results = AnalyzeTwoAnglePreparedPeaks(normal, angled, maxDeadLayer, allowFallback);
	if ((textReport && strlen(textReport)) || (csvReport && strlen(csvReport))) {
		WriteTwoAngleDiagnostics(textReport, csvReport, results, normal, angled);
	}
	return results;
}

Bool_t AculCalibration::WriteTwoAngleDiagnostics(
		const char* textFile,
		const char* csvFile,
		const vector<AculTwoAngleStripResult>& results,
		const AculTwoAnglePeakSet& normalMeasurement,
		const AculTwoAnglePeakSet& angledMeasurement)
{
	if (textFile && strlen(textFile)) {
		ofstream text(textFile);
		if (!text.is_open()) {
			Error("AculCalibration::WriteTwoAngleDiagnostics", "Cannot open text report %s.", textFile);
			return kFALSE;
		}
		text << "# two-angle dead-layer calibration report\n";
		text << "# detector " << normalMeasurement.detectorName << "\n";
		text << "# normal " << normalMeasurement.measurementName << " angleDeg " << normalMeasurement.angleDeg << "\n";
		text << "# angled " << angledMeasurement.measurementName << " angleDeg " << angledMeasurement.angleDeg << "\n";
		text << "# strip status deadLayer_um objectiveRms tempA tempB finalA finalB finalRms finalMaxResidual normalChannels angledChannels normalEnergies angledEnergies finalResiduals message\n";
		for (size_t i = 0; i < results.size(); i++) {
			const AculTwoAngleStripResult& r = results[i];
			text << r.strip << " "
				<< r.status << " "
				<< setprecision(10) << r.deadLayer << " "
				<< r.deadLayerObjectiveRms << " "
				<< r.temporaryA << " "
				<< r.temporaryB << " "
				<< r.finalA << " "
				<< r.finalB << " "
				<< r.finalRms << " "
				<< r.finalMaxResidual << " "
				<< FormatVector(r.normalChannels) << " "
				<< FormatVector(r.angledChannels) << " "
				<< FormatVector(r.normalEnergies) << " "
				<< FormatVector(r.angledEnergies) << " "
				<< FormatVector(r.finalResiduals) << " "
				<< r.message << "\n";
		}
	}

	if (csvFile && strlen(csvFile)) {
		ofstream csv(csvFile);
		if (!csv.is_open()) {
			Error("AculCalibration::WriteTwoAngleDiagnostics", "Cannot open CSV report %s.", csvFile);
			return kFALSE;
		}
		csv << "strip,status,fallbackUsed,angleDeg,deadLayer_um,objectiveRms,tempA,tempB,finalA,finalB,finalRms,finalMaxResidual,normalChannels,angledChannels,normalEnergies,angledEnergies,finalResiduals,message\n";
		for (size_t i = 0; i < results.size(); i++) {
			const AculTwoAngleStripResult& r = results[i];
			csv << r.strip << ","
				<< r.status << ","
				<< (r.fallbackUsed ? 1 : 0) << ","
				<< r.angleDeg << ","
				<< setprecision(12) << r.deadLayer << ","
				<< r.deadLayerObjectiveRms << ","
				<< r.temporaryA << ","
				<< r.temporaryB << ","
				<< r.finalA << ","
				<< r.finalB << ","
				<< r.finalRms << ","
				<< r.finalMaxResidual << ","
				<< "\"" << FormatVector(r.normalChannels) << "\","
				<< "\"" << FormatVector(r.angledChannels) << "\","
				<< "\"" << FormatVector(r.normalEnergies) << "\","
				<< "\"" << FormatVector(r.angledEnergies) << "\","
				<< "\"" << FormatVector(r.finalResiduals) << "\","
				<< r.message << "\n";
		}
	}
	return kTRUE;
}

Bool_t AculCalibration::WriteTwoAngleDiagnosticGraphs(
		const char* rootFile,
		const vector<AculTwoAngleStripResult>& results)
{
	if (!rootFile || !strlen(rootFile)) {
		Error("AculCalibration::WriteTwoAngleDiagnosticGraphs", "Empty ROOT output file name.");
		return kFALSE;
	}

	TFile output(rootFile, "RECREATE");
	if (output.IsZombie()) {
		Error("AculCalibration::WriteTwoAngleDiagnosticGraphs", "Cannot create ROOT file %s.", rootFile);
		return kFALSE;
	}

	TGraph deadLayerGraph;
	TGraph finalRmsGraph;
	TGraph objectiveGraph;
	deadLayerGraph.SetName("two_angle_dead_layer_um");
	deadLayerGraph.SetTitle("Two-angle fitted dead layer;strip;dead layer, um");
	finalRmsGraph.SetName("two_angle_final_rms");
	finalRmsGraph.SetTitle("Final normal-incidence calibration RMS;strip;RMS, MeV");
	objectiveGraph.SetName("two_angle_objective_rms");
	objectiveGraph.SetTitle("Two-angle objective RMS;strip;RMS, MeV");

	Int_t point = 0;
	for (size_t i = 0; i < results.size(); i++) {
		const AculTwoAngleStripResult& r = results[i];
		if (r.status != "ok" && r.status != "fallback") {
			continue;
		}
		deadLayerGraph.SetPoint(point, r.strip, r.deadLayer);
		finalRmsGraph.SetPoint(point, r.strip, r.finalRms);
		objectiveGraph.SetPoint(point, r.strip, r.deadLayerObjectiveRms);
		point++;
	}

	deadLayerGraph.Write();
	finalRmsGraph.Write();
	objectiveGraph.Write();
	output.Close();
	return kTRUE;
}

void AculCalibration::PrintInputParameters()
{
	//print alpha source parameters

	cout << "AculCalibration::PrintInputParameters:" << endl;
	cout << "\tParFileName: " << fParFileName << endl;
	cout << "\tNumber of peaks: " << kRaNOPEAKS << endl;
	for (Int_t i = 0; i < kRaNOPEAKS; i++) {
		cout << "\t\tfEnergyInput[" << i << "] = " << fEnergyInput[i] << endl;
	}
	cout << "\tqAlphaInput: " << fQAlphaInput << "; firstAlphaParentMass: " << fFirstAlphaParentMass << ";" << endl;
	cout << "\tEnergies used for calibration:" << endl;
	cout << "\t(deadLayer: " << fDeadLayer << " mcm)" << endl;
//	Info("AculCalibration::PrintInputParameters", "Number of peaks: %d", kRaNOPEAKS);
	for (Int_t i = 0; i < kRaNOPEAKS; i++) {
		cout << "\t\tfEnergy[" << i << "] = " << fEnergy[i] << endl;
	}

	cout << "\tlowerChannel: " << fLowerChannel << "; upperChannel: " << fUpperChannel << ";" << endl;
	cout << "\tlowerPeakHight: " << fLowerPeakRelativeHight << "; upperPeakHight: " << fUpperPeakRelativeHight << ";" << endl;
	cout << "\tfitHightThreshold: " << fFitPeakThreshold << "; minFitSigma: " << fFitMinSigma << ";" << endl;
	cout << "\tpeakPositionTolerance: " << fPeakPositionTolerance << ";" << endl;
	cout << "\tfitFunctionLineWidth: " << fFitFuncLineWidth << ";" << endl;
	cout << "\tuppersubaddress: " << fuppersubaddress << ";" << endl;
	cout << "\tlowersubaddress: " << flowersubaddress << ";" << endl;
	cout << "\tblock: " << block << ";" << endl;
	cout << "\ttreename: " << treename << ";" << endl;



	return;

}


Double_t AculCalibration::GetA(Int_t i) 
{
	if (i >= fA.GetSize()) //if i >= number of array element
	{
		return 0.;
	}
	return fA[i];
}

Double_t AculCalibration::GetB(Int_t i) 
{
	if (i >= fB.GetSize()) 
	{
		return 0.;
	}
	return fB[i];
}



Bool_t AculCalibration::SetCalibrationParameters(const char* calparfile)
{

	const Int_t lineLength = 200;
	char line[lineLength];
//	Int_t	crate;
	char	crate[100];
//	Int_t	i, j;
	Int_t	j;
	char	cA[40], cB[40], cC[40], cD[40];

	//open file with calibration parameters
	ifstream calFileR;
	calFileR.open(calparfile);


	if( !calFileR.is_open() ) {
		Error("SetCalibrationParameters", "File %s with calibration data was not opened", calparfile);
		return kFALSE;
	}

	Reset();

	//read calibration parameters from file
	while (!calFileR.eof()) {
		calFileR.getline(line, lineLength);
		if ( line[0] != '*' && line[0] != '#' && line[0] != '%' && (line[0] != '/' && line[1] != '/') ) {		//possible comment characters
			sscanf(line, "%s %d %s %s %s %s", crate, /*&i*,*/ &j, cA, cB, cC, cD); //
			fA[j] = atof(cA);
			fB[j] = atof(cB);
//			fC[0][j] = atof(cC);
//			fD[0][j] = atof(cD);
		}

	}
	calFileR.close();
	return kTRUE;
}

void AculCalibration::PrintCalibrationParameters()
{
	for (Int_t j = 0; j < ADDRESSNUMBER; j++) {
		cout << "C3[" << j << "]" << setw(10) << fA[j] << setw(10) << fB[j] /*<< setw(10) << fC[j] << setw(10) << fD[j]*/ << endl;
	}
}

void AculCalibration::ShowRawSpectra(const char* filename, const Int_t block, TCanvas* rawCanvas, Int_t xaxismin, Int_t xaxismax, /*TObjArray* histList,*/ const Int_t subaddress)
{
	//Displays the spectrum from a file, divides the canvas into a sufficient number of pads and displays spectrums of each  
	//block subaddress on the suitable pads.
	//
        //  filename: input .root file containing spectra to be showed
	//  block: block which will be drawn
	//  rawCanvas: canvas on which one you will see the spectrum
	//  xaxismin: minimal channel, which will be displayed
	//  xaxismax: maximal channel, which will be displayed
	//  subaddress:

	Char_t address[40];
	Char_t histName[40];
	Char_t fillCommand[40];
	Char_t fillCondition[40];


	if (!rawCanvas) {
		//rawCanvas = new TCanvas("RawSpectra", "Raw spectra in channels", 1);
		cout << "You have to assign TCanvas for raw spectra drawing" << endl;
		return;
	}

	rawCanvas->Clear();

//	cout << "hovno" << endl;

	rawCanvas->SetFillColor(10);

//	cout << "hovno" << endl;
	TFile *fr = new TFile(filename);
	if (fr->IsOpen() == 0) {
		cout << endl << "File " << filename << " was not opened and won't be processed" << endl << endl;
		return;
	}
	TH1I *hRead = 0;
	TTree *tr = (TTree*)fr->Get("RAW");
//	cout << "hovno" << endl;


	if (subaddress > 15) {
		rawCanvas->Divide(4, 4);
		rawCanvas->SetFillColor(10);
//		cout << "hovno" << endl;
		for (Int_t i = 0; i < 16; i++) {
			cout << i << endl;
			rawCanvas->cd(i+1);
			hRead = new TH1I("name", "title", 4096, 0, 4095);
			sprintf(address, "C3[%d][%d]", block, i);
			sprintf(histName, "H3[%d][%d]", block, i);
			hRead->SetName(histName);
			sprintf(fillCommand, "%s >> %s", address, hRead->GetName());
			sprintf(fillCondition, "%s > 0", address);
			tr->Draw(fillCommand, fillCondition, "");
			if (hRead) {
				hRead->SetDirectory(0);
//				cout << hRead->GetEntries() << endl;
//				if (fHRawList) {
//					fHRawList->Add(hRead);
					fHRawList.Add(hRead);
//				}
				hRead->SetAxisRange(xaxismin, xaxismax);
			}
		}//for
	}
	else {
		fr->cd();
		hRead = new TH1I("name", "title", 4096, 0, 4095);
		sprintf(address, "C3[%d][%d]", block, subaddress);
		sprintf(histName, "H3[%d][%d]", block, subaddress);
		hRead->SetName(histName);
		sprintf(fillCommand, "%s >> %s", address, hRead->GetName());
		sprintf(fillCondition, "%s > 0", address);
//		cout << fillCommand << setw(20) << fillCondition << endl;
		tr->Draw(fillCommand, fillCondition, "goff");
		if (hRead) {
			hRead->SetDirectory(0);
//			if (fHRawList) {
//				fHRawList->Add(hRead);
//			}
			fHRawList.Add(hRead);
			hRead->Draw();
			hRead->SetAxisRange(xaxismin, xaxismax);
		}
	}//else

	fr->Close();

	rawCanvas->Update();

	return;

}

void AculCalibration::ShowSpectra(const char* filename, TCanvas* rawCanvas, Option_t *option, Int_t xaxismin, Int_t xaxismax, const Int_t subaddress)
{
	//filename: input .root file with saved filled histograms to be showed
	//rawCanvas: canvas on which one you will see the spectrum
	//option: THStack options
	//xaxismin: Minimum channel, which will be displayed
	//xaxismax: Maximum channel, which will be displayed
	//subaddress:
 
	TString opt = option;
	opt.ToLower();

	if (!rawCanvas) {
		Error("ShowRawSpectra", "You have to assign TCanvas for raw spectra drawing");
		return;
	}
	rawCanvas->Clear();

	TFile fr(filename);
	if (fr.IsOpen() == 0) {
		Error("ShowRawSpectra", "File %s was not opened and won't be processed", filename);
		return;
	}

	TList *histList;
	histList = fr.GetListOfKeys();
	Int_t listEntries = histList->GetEntries();
	TH1 *hDraw = 0;
	DeleteStacks();

	if (subaddress >= listEntries) {
		fCurrentHStack = new THStack();
		for (Int_t i = 0; i < listEntries; i++) {	//zkontrolovat hranice
			Info("ShowRawSpectra", "Histogram with spectrum of subaddress %d is loading", i);
			fr.GetObject(histList->At(i)->GetName(), hDraw);
			if (hDraw) {
				hDraw->SetDirectory(0);
				fCurrentHistList.Add(hDraw);
				fCurrentHStack->Add(hDraw);
			}
		}//for
		if ( !fCurrentHStack->GetHists()->IsEmpty() ) {
			Info("ShowRawSpectra", "Histogram stack drawing");
			fCurrentHStack->Draw(opt.Data());
		}
	}//if all subaddresses
	else {
		//zkontrolovat
		fr.GetObject(histList->At(subaddress)->GetName(), hDraw);
		if (hDraw) {
			hDraw->SetAxisRange(xaxismin, xaxismax, "X");
			hDraw->Draw();
			hDraw->SetDirectory(0);
			fCurrentHistList.Add(hDraw);
		}
	}//else

	fr.Close();
	rawCanvas->Update();
	return;
}

void AculCalibration::FillRawSpectraFile(const char* rawdatafile, const char* block, const char* treename, TCanvas* rawCanvas, Option_t *option, Int_t xaxismin, Int_t xaxismax)
{
	//filename: input .root file containing spectra to be showed
	//block:
	//rawCanvas:
	//xaxismin:
	//xaxismax:

	//variables to be became function parameter
	TString opt(option);
	opt.ToLower();

	if (!rawCanvas) {
		Error("ShowRawSpectra", "You have to assign TCanvas for raw spectra drawing");
		return;
	}
	rawCanvas->Clear();

	TFile fr(rawdatafile);
	if (fr.IsOpen() == 0) {
		Error("ShowRawSpectra", "File %s was not opened and won't be processed", rawdatafile);
		return;
	}
	TTree *tr = (TTree*)fr.Get(treename);

	char outputfile[300];
	sprintf(outputfile, "%s[]Raw.root", block);
	TFile fw(outputfile, opt.Data());
	if (fw.IsOpen() == 0) {
		Error("CalculateCalibParameters", "Output file %s was not created.", outputfile);
		return;
	}
	if (fw.IsWritable() == 0) {
		Error("CalculateCalibParameters", "Output file %s is not writable. Set option to \"RECREATE\".", outputfile);
		return;
	}

	char address[40];
	char histName[40];
	char histTitle[40];
	char fillCommand[40];
	char fillCondition[40];

	fw.cd();
	TH1I *hRead = 0;

	for (Int_t i = 0; i < 16; i++) {	//zkontrolovat hranice
		cout << i << endl;	//predelat na info
		hRead = new TH1I("name", "title", 4096, 0, 4095);
		sprintf(address, "%s[%d]", block, i);
		sprintf(histName, "%s[%d]", block, i);
		sprintf(histTitle, "%s : %s", rawdatafile, histName);
		hRead->SetName(histName);
		hRead->SetTitle(histTitle);
		sprintf(fillCommand, "%s >> %s", address, hRead->GetName());
		sprintf(fillCondition, "%s > 0", address);
		tr->Draw(fillCommand, fillCondition, "goff");		//prozkoumat goff
		hRead->Write();
	}//for

	fw.Close();
	fr.cd();
	delete tr;
	fr.Close();

	rawCanvas->Update();

	return;
}

Bool_t AculCalibration::CalculateCalibParameters(/*const char* inputfile, const char* block, const Int_t address,
		const char* treename, */
		Int_t lowerchannel, Int_t upperchannel)
//		, Int_t nEBins, Int_t lowersubaddress,Int_t uppersubaddress)
{
	if (kRaNOPEAKS == 0) {
		Error("CalculateCalibParameters", "Alpha source parameters was not read");
		return 0;
	}

	//muzu nechat
	if ( (fuppersubaddress - flowersubaddress) >= ADDRESSNUMBER ) {
		Error("CalculateCalibParameters", "Possible subaddress values have to be in range 0 - %d", ADDRESSNUMBER - 1);
		return 0;
	}

	//auxiliary variables, particularly for text parameter fields
	TString CalibName;
	CalibName.Form("%s.cal", block);
	TString oFileName =  fWorkDirectory +  CalibName; 

	//creation of the output text file 

	ofstream outcalfile;
	outcalfile.open(oFileName.Data());
	if (!outcalfile.is_open()) {
		Error("CalculateCalibParameters", "Output file %s was not opened", oFileName.Data());
		return 0;
	}//if



	//input file with raw data opening
//	TString iFileName = inputfile;
	TFile *fr = new TFile(fInputRootFile.Data());
	if ( !fr->IsOpen() ) {
		Error("CalculateCalibParameters", "File %s was not opened and won't be processed", fInputRootFile.Data());
		return 0;
	}
	TTree *tr = (TTree*)fr->Get(treename);
	if (!tr) {
		Error("CalculateCalibParameters", "Tree %s was not found in file %s", treename, fInputRootFile.Data());
		return 0;
	}


	//promenne potrebne pro fitovani: presunout nize
	//pohlidat delete
	TF1 	*calFunction = new TF1("calib", "pol1", 0, 1000);	//predelat jako lokalni promennou fce (nebo snad tridy?)
	TGraph 	*calGraph = new TGraph(kRaNOPEAKS, fPeak.GetArray(), fEnergy.GetArray()); //lokalni promenna, dohodit pocet vstupu pomoci parametru

	TString	detectorChannel;
	TString	histName;
	TString histTitle;
	TString fillCommand;
	TString	fillCondition;
	Int_t	fitControl = 0;

	//predelat nazvy histogramu
	//zrusit cyklus, napsat jako fci
	//raw data histogram filling
	TH1I *hRaw = 0;
	TRandom3 ranGen(1);


	for (Int_t i = flowersubaddress; i <= fuppersubaddress; i++) {
		printf("\n\n");
		Info("CalculateCalibParameters", "Calculating calibration parameters for detector channel %s[%d].", block, i);
		//TH1I object preparing
		hRaw = new TH1I("name", "title", 4096, 0, 4095);					//nastavovat hranice histogramu podle parametru fce
		detectorChannel.Form("%s[%d]", block, i);
		histName.Form("Hist%s[%d]", block, i);
		hRaw->SetName(histName.Data());
		fillCommand.Form("%s >> %s", detectorChannel.Data(), hRaw->GetName());
		fillCondition.Form("%s > %d && %s < %d",
				detectorChannel.Data(), lowerchannel, detectorChannel.Data(), upperchannel);
		//filling from the .root raw data file and content arrangement
		tr->Draw(fillCommand.Data(), fillCondition.Data(), "goff");

		//spectrum analysis
		fitControl = PeaksFitting(hRaw, "", fFitMinSigma);
		Info("CalculateCalibParameters", "Value of fitControl is: %d", fitControl);		//ok

		
		//calibration parameters calculation																//ok
		for (Int_t j = 0; j < kRaNOPEAKS; j++) {		//delat podle poctu zkoumanych piku
			calGraph->SetPoint(j, fPeak[j], fEnergy[j]);  //calibration graph filling
			//printf("\tPeak\t%f and energy\t%f\n", fPeak[j], fEnergy[j]);
		}//for
		calGraph->Fit(calFunction, "Q", "goff", 0, 4096);  //omezit hlasitost fitovani, udelat volitelne, dodelat volby rozsahu


			outcalfile
			<< setprecision(4) << calFunction->GetParameter(0)<< "\t"
			<< setprecision(4) << calFunction->GetParameter(1)<< "\t\t" //E=a+bN
			<< fitControl
			<< endl;
	

	}//for; detector channels

	fr->Close();
	outcalfile.close();

	return 1;

}

// 22.12.16
Bool_t AculCalibration::CalculateCalibParameters5points(const char* inputfile, const char* pedestalfile, const char* block, const char* treename, Int_t lowerchannel, Int_t upperchannel, Int_t nEBins, Int_t lowersubaddress, Int_t uppersubaddress) {
	
	if (kRaNOPEAKS == 0) {
		Error("CalculateCalibParameters", "Alpha source parameters was not read");
		return 0;
	}


	if ( (uppersubaddress - lowersubaddress) >= ADDRESSNUMBER ) {
		Error("CalculateCalibParameters", "Possible subaddress values have to be in range 0 - %d", ADDRESSNUMBER - 1);
		return 0;
	}

	TString oFileName;

	//creation of the output text file
	if ( (lowersubaddress == 0) && (uppersubaddress == ADDRESSNUMBER-1) ) { oFileName.Form("%s[].cal", block); }
	else {
		if (lowersubaddress == uppersubaddress) { oFileName.Form("%s[%d].cal", block, lowersubaddress); }
		else { oFileName.Form("%s[%d-%d].cal", block, lowersubaddress, uppersubaddress); }
	}

	ofstream outcalfile;
	outcalfile.open(oFileName.Data());
	if (!outcalfile.is_open()) {
		Error("CalculateCalibParameters", "Output file %s was not opened", oFileName.Data());
		return 0;
	}

	//creation of the output root file
	if ( (lowersubaddress == 0) && (uppersubaddress == ADDRESSNUMBER-1) ) { oFileName.Form("%s[].root", block); }
	else {
		if (lowersubaddress == uppersubaddress) { oFileName.Form("%s[%d].root", block, lowersubaddress); }
		else { oFileName.Form("%s[%d-%d].root", block, lowersubaddress, uppersubaddress); }
	}
	fCalInformation = new TFile(oFileName.Data(), "RECREATE");
	if ( !fCalInformation->IsOpen() ) {
		Error("CalculateCalibParameters", "File %s was not opened and won't be processed", oFileName.Data());
		return 0;
	}

	//input file with raw data opening
	TString iFileName = inputfile;
	TFile *fr = new TFile(iFileName.Data());
	if ( !fr->IsOpen() ) {
		Error("CalculateCalibParameters", "File %s was not opened and won't be processed", iFileName.Data());
		return 0;
	}
	TTree *tr = (TTree*)fr->Get(treename);
	if (!tr) {
		Error("CalculateCalibParameters", "Tree %s was not found in file %s", treename, iFileName.Data());
		return 0;
	}

	TF1 	*calFunction = new TF1("calib", "pol1", 0, 1000);
	TGraph 	*calGraph = new TGraph(kRaNOPEAKS, fPeak.GetArray(), fEnergy.GetArray());
	TString	detectorChannel;
	TString	histName;
	TString histTitle;
	TString fillCommand;
	TString	fillCondition;
	TString line;
	Int_t	fitControl = 0;

	TH1I *hRaw = 0;
	TH1F *hEnergy = 0;

	TRandom3 ranGen(1);

	//outputfile with calibrated spectra
//	if ( (lowersubaddress == 0) && (uppersubaddress == ADDRESSNUMBER-1) ) { oFileName.Form("%s[]E.root", block); }
//	else {
//		if (lowersubaddress == uppersubaddress) { oFileName.Form("%s[%d]E.root", block, lowersubaddress); }
//		else { oFileName.Form("%s[%d-%d]E.root", block, lowersubaddress, uppersubaddress); }

//	}
	TFile *fw = new TFile(oFileName.Data(), "RECREATE");
	if (fw->IsOpen() == 0) {
		Error("CalculateCalibParameters", "File %s was not created and won't be processed\n\n", oFileName.Data());
		return 1;
	}

	std::ifstream inpfile(pedestalfile);

	if( !inpfile.is_open() ) {
		Error("CalculateCalibParameters5points", "File %s was not opened and won't be processed", pedestalfile);
		return 0;

	}



	for (Int_t i = lowersubaddress; i <= uppersubaddress; i++) {
		printf("\n\n");
		Info("CalculateCalibParameters", "Calculating calibration parameters for detector channel %s[%d].", block, i);
		hRaw = new TH1I("name", "title", 4096, 0, 4095);					
		detectorChannel.Form("%s[%d]", block, i);
		histName.Form("Hist%s[%d]", block, i);
		hRaw->SetName(histName.Data());
		fillCommand.Form("%s >> %s", detectorChannel.Data(), hRaw->GetName());
		fillCondition.Form("%s > %d && %s < %d",
				detectorChannel.Data(), lowerchannel, detectorChannel.Data(), upperchannel);
		//filling from the .root raw data file and content arrangement
		tr->Draw(fillCommand.Data(), fillCondition.Data(), "goff");

		//spectrum analysis
		fitControl = PeaksFitting(hRaw, "", fFitMinSigma);
		Info("CalculateCalibParameters", "Value of fitControl is: %d", fitControl);	

		//incorrectly treated spectrum output
//		if (fitControl != 0 && fCalInformation->IsOpen()) {			
//			outcalfile << setw(39) << fitControl << endl;
//			fCalInformation->cd();
//			hRaw->SetLineColor(2);	//red
//			fCalInformation->cd();
//			hRaw->Write();
//			continue;
//		}
		//correctly treated spectrum saving
//		if (fCalInformation->IsOpen()) {
//			fCalInformation->cd();
//			hRaw->SetLineColor(kGreen+3);	//green
//			hRaw->SetFillColor(kGreen+1);
//			hRaw->Write();
//		}
		Double_t pedestal_chan = 0.;
		line.ReadLine(inpfile);
		sscanf(line.Data(), "%lf", &pedestal_chan);
		Info("CalculateCalibParameters5points", "Pedestal position in channels for %s[%d] is %f", block, i, pedestal_chan);
		for (Int_t j = 0; j < kRaNOPEAKS; j++) {		
			calGraph->SetPoint(j, fPeak[j], fEnergy[j]);  //calibration graph filling
			//printf("\tPeak\t%f and energy\t%f\n", fPeak[j], fEnergy[j]);
		}
		calGraph->SetPoint(4, pedestal_chan, 0.); //additional point for the pedestal at 0 MeV
		calGraph->Fit(calFunction, "Q", "goff", 0, 4096);  
	//	outcalfile
	//		<< block << "\t"
			/*<< address << "\t"*/
	/*		<< i << "\t"
			<< setprecision(4) << calFunction->GetParameter(1) << "\t"
			<< setprecision(4) << calFunction->GetParameter(0) << "\t\t"
			<< fitControl
			<< endl;
	*/
			outcalfile
			<< setprecision(4) << calFunction->GetParameter(0)<< "\t"
			<< setprecision(4) << calFunction->GetParameter(1)<< "\t\t" //E=a+bN
			<< fitControl
			<< endl;


		fA[i] = calFunction->GetParameter(1);
		fB[i] = calFunction->GetParameter(0);

		//calibration of raw spectra using obtained parameters
		Info("CalculateCalibParameters", "Energy spectrum from address %s[%d] calculating", block, i);
		histName.Form("%sE", hRaw->GetName());
		histTitle.Form("%s: %s", iFileName.Data(), histName.Data());
		hEnergy = new TH1F(histName.Data(), histTitle.Data(), nEBins, 0., 10.);

		for (Int_t j = lowerchannel; j < upperchannel; j++) {
			Int_t binCont = (Int_t)hRaw->GetBinContent(j);
			//				cout << j << ":\t" << hRaw->GetBinContent(j) << endl;
			for (Int_t k = 0; k < binCont; k++) {
				hEnergy->Fill( fA[i]*( j+ranGen.Uniform(0., 1.0) ) + fB[i] );
			}
		}

		fw->cd();
		hEnergy->Write();


	}//for; detector channels

	fw->Close();
	fr->Close();
	inpfile.close();
	fCalInformation->Close();
	outcalfile.close();

	return 1;



}


Bool_t AculCalibration::Mycalc(Int_t lowerchannel, Int_t upperchannel) {
//	Bool_t AculCalibration::Mycalc(const char* inputfile, const char* pedestalfile, Int_t lowerchannel, Int_t upperchannel) {
	
	if (kRaNOPEAKS == 0) {
		Error("CalculateCalibParameters", "Alpha source parameters was not read");
		return 0;
	}

	//input file with data opening
//	TString iFileName = inputfile;
	TFile *fr = new TFile(fInputRootFile.Data());
	if ( !fr->IsOpen() ) {
		Error("CalculateCalibParameters", "File %s was not opened and won't be processed", fInputRootFile.Data());
		return 0;
	}
	TTree *tr = (TTree*)fr->Get(treename);
	if (!tr) {
		Error("CalculateCalibParameters", "Tree %s was not found in file %s", treename, fInputRootFile.Data());
		return 0;
	}

	TString PedestalName;
	PedestalName.Form("pedestal_%s.par", block);
	TString name =  fWorkDirectory +  PedestalName;  
	

//	std::ifstream inpfile(pedestalfile);
	std::ifstream inpfile(name);

	if( !inpfile.is_open() ) {
		Error("CalculateCalibParameters5points", "File %s was not opened and won't be processed", name.Data());
		return 0;

	}

	TString rawName;
	rawName.Form("%s_Raw.root", block);
	//creation of the raw output root file
  TString rawFileName = fWorkDirectory + rawName;

	fCalInformation = new TFile(rawFileName.Data(), "RECREATE");
	if ( !fCalInformation->IsOpen() ) {
		Error("Mycalc", "File %s was not opened and won't be processed", rawFileName.Data());
		return 0;
	}

	TString EName;
	EName.Form("%s_E.txt", block);
  TString EFileName =  fWorkDirectory +  EName;
	//creating Energies table

	ofstream outEfile;
	outEfile.open(EFileName.Data());
	if (!outEfile.is_open()) {
		Error("CalculateCalibParameters", "Output file %s was not opened", EFileName.Data());
		return 0;
	}


	TString CalibName;
	CalibName.Form("%s.cal", block);
	TString oFileName =  fWorkDirectory +  CalibName; 

	//creation of the output text file

	ofstream outcalfile;
	outcalfile.open(oFileName.Data());
	if (!outcalfile.is_open()) {
		Error("CalculateCalibParameters", "Output file %s was not opened", oFileName.Data());
		return 0;
	}


	TF1 	*calFunction = new TF1("calib", "pol1", 0, 1000);
	TGraph 	*calGraph = new TGraph(kRaNOPEAKS, fPeak.GetArray(), fEnergy.GetArray());
	
	TString	detectorChannel;
	TString	histName;
	TString histTitle;
	TString fillCommand;
	TString	fillCondition;
	TString line;
	Int_t	fitControl = 0;

	TH1I *hRaw = 0;

  outcalfile << calFunction->GetNpar() << endl << fuppersubaddress - flowersubaddress + 1 << endl;

	for (Int_t i = flowersubaddress; i <= fuppersubaddress; i++) {
		printf("\n\n");
		Info("CalculateCalibParameters", "Calculating calibration parameters for detector channel %s[%d].", block, i);
		hRaw = new TH1I("name", "title", 4096, 0, 4095);					
		detectorChannel.Form("%s[%d]", block, i);
		histName.Form("Hist%s[%d]", block, i);
		hRaw->SetName(histName.Data());
		fillCommand.Form("%s >> %s", detectorChannel.Data(), hRaw->GetName());
		fillCondition.Form("%s > %d && %s < %d",
				detectorChannel.Data(), lowerchannel, detectorChannel.Data(), upperchannel);
//------------------------------------------------------------------------------------
 /*   if(i==14) {
      fillCondition.Form("%s > %d && %s < %d",
				detectorChannel.Data(), 120, detectorChannel.Data(), upperchannel);
    }
    if(i==1) {
      fillCondition.Form("%s > %d && %s < %d",
				detectorChannel.Data(), 150, detectorChannel.Data(), upperchannel);
    }    */
//------------------------------------------------------------------------------------
		//filling from the .root raw data file and content arrangement
		tr->Draw(fillCommand.Data(), fillCondition.Data(), "goff");

		//spectrum analysis
		fitControl = PeaksFitting(hRaw, "", fFitMinSigma);
		Info("CalculateCalibParameters", "Value of fitControl is: %d", fitControl);	

		Double_t pedestal_chan = 0.;
		line.ReadLine(inpfile);
		sscanf(line.Data(), "%lf", &pedestal_chan);
		Info("CalculateCalibParameters5points", "Pedestal position in channels for %s[%d] is %f", block, i, pedestal_chan);
		for (Int_t j = 0; j < kRaNOPEAKS; j++) {
			calGraph->SetPoint(j, fPeak[j], fEnergy[j]);  //calibration graph filling
			//printf("\tPeak\t%f and energy\t%f\n", fPeak[j], fEnergy[j]);
      cout << "Peak and energy " << fPeak[j] << " " << fEnergy[j] << endl;
		} 
		//calGraph->SetPoint(4, pedestal_chan, 0.); //additional point for the pedestal at 0 MeV
		calGraph->Fit(calFunction, "Q", "goff", 0, 4096);
  	for (Int_t j = 0; j < kRaNOPEAKS; j++) {
      outEfile << setprecision(4) << fPeak[j]<< "\t";
		}
    outEfile << endl;  
    cout << calFunction->GetParameter(0) << " " <<  calFunction->GetParameter(1) << endl;
		outcalfile
			<< setprecision(4) << calFunction->GetParameter(0)<< "\t"
			<< setprecision(4) << calFunction->GetParameter(1)<< "\t\t" //E=a+bN
			<< fitControl
			<< endl;
    fCalInformation->cd(); // writing hists into rawOutFile
		hRaw->Write();
	}
  outcalfile << endl;
  fCalInformation->Close();
  return 1;
}



//my stupid function, idea stolen from function CalculateCalibParameters, 13.12.2016
Bool_t AculCalibration::ShowFullCalibratedSpectra(//const char* inputfile, 
//	const char* calibparfile, 
	Int_t lowerchannel, Int_t upperchannel) {

	if ( (fuppersubaddress - flowersubaddress) >= ADDRESSNUMBER ) {
		Error("CalculateCalibParameters", "Possible subaddress values have to be in range 0 - %d", ADDRESSNUMBER - 1); 
		return 0;
	}

	//input file with data opening
//	TString iFileName = inputfile;
	TFile *fr = new TFile(fInputRootFile.Data());
	if ( !fr->IsOpen() ) {
		Error("CalculateCalibParameters", "File %s was not opened and won't be processed", fInputRootFile.Data());
		return 0;
	}
	TTree *tr = (TTree*)fr->Get(treename);
	if (!tr) {
		Error("CalculateCalibParameters", "Tree %s was not found in file %s", treename, fInputRootFile.Data());
		return 0;
	}

	TString	histName;
	TString histTitle;
	TString	detectorChannel;
	TString line;
	TString fillCommand;
	TString	fillCondition;
	TH1I *hRaw2 = 0;	
	TH1F *hEnergy = 0;
	TRandom3 ranGen(1);

	TString CalibName;
	CalibName.Form("%s.cal", block);
	TString calibparfile =  fWorkDirectory +  CalibName;  
	std::ifstream infile(calibparfile);

	TString calf;
	calf.Form("%s_full_calibrated_spectra.root", block);
	TString oFileName =  fWorkDirectory +  calf;  
	


	char cA[40], cB[40];
	Int_t nEBins;
	nEBins = upperchannel- lowerchannel;

	
	TFile *fw = new TFile(oFileName.Data(), "RECREATE");
	if (fw->IsOpen() == 0) {
		Error("CalculateCalibParameters", "File %s was not created and won't be processed\n\n", oFileName.Data());
		return 1;
	}

	for(Int_t i = flowersubaddress-2; i <= fuppersubaddress; i++) { //this loop goes from 1 to 16 or 1 to 32 i.e. number of strips with skipping first two lines
    
    if(i<flowersubaddress) {
      line.ReadLine(infile);
      continue;
    }

		Info("ShowFullCalibratedSpectra", "Calculate full spectra for channel %d", i);
		line.ReadLine(infile);
		sscanf(line.Data(), "%s %s", cA, cB);
		fA[i] = atof(cA);
		fB[i] = atof(cB);
		//cout<<fA[i]<<" "<<fB[i]<<"!!!!!!!!!!!!!"<<endl;
		hRaw2 = new TH1I("name", "title", 4096, 0, 4095);
		detectorChannel.Form("%s[%d]", block, i);
		histName.Form("Hist%s[%d]", block, i);
		hRaw2->SetName(histName.Data());
	
		fillCommand.Form("%s >> %s", detectorChannel.Data(), hRaw2->GetName());
		fillCondition.Form("%s > %d && %s < %d", detectorChannel.Data(), lowerchannel, detectorChannel.Data(), upperchannel);

		//filling from the .root raw data file and content arrangement
		tr->Draw(fillCommand.Data(), fillCondition.Data(), "goff");

		histName.Form("%sEfull", hRaw2->GetName());
		histTitle.Form("%s: %s", fInputRootFile.Data(), histName.Data());
		hEnergy = new TH1F(histName.Data(), histTitle.Data(), nEBins, 0., 10.);


		for (Int_t j = lowerchannel; j < upperchannel; j++) {
			Int_t binCont = (Int_t)hRaw2->GetBinContent(j);
			for (Int_t k = 0; k < binCont; k++) {
				//cout<<fA[i]<<" "<<fB[i]<<"!!!!!!!!!!!!!"<<endl;
				hEnergy->Fill(fA[i] + fB[i]*( j+ranGen.Uniform(0., 1.0) ));
			}
		}
		fw->cd();
		hEnergy->Write();
	}
	fw->Close();
	fr->Close();
	infile.close();
	return 1;
}

void AculCalibration::CalibrateRawSpectra() {
	//todo: implement this function
	//function parameters:
	const char* iFileName = "clb01_0001.root";
	const char* treeName = "AnalysisxTree";
	const char* branchName = "LiEvent.SQ22[32]";
	//const Int_t address = 22;
	const Int_t lowerElement = 0;
	const Int_t upperElement = 32;
	const Int_t lowerChannel = 100, upperChannel = 4096;
	const Int_t nEBins = 1000;
	//optional:
	Long64_t nentries = 0;


	//function itself:
	TString iFile = iFileName;
	TFile fr( iFile.Data() );
	if ( !fr .IsOpen() ) {
		Error("CalibrateRawSpectra", "File %s was not opened and won't be processed", iFile.Data());
		return;
	}

	TString tName = treeName;
	TTree *tr = (TTree*)fr.Get(tName.Data());
	if (!tr) {
		Error("CalibrateRawSpectra", "Tree %s was not found in file %s", tName.Data(), iFile.Data());
		return;
	}
	tr->SetMakeClass(1);


	//!!!!!!!!!!!!!!!!!!!!!!!!!!!
	UShort_t variable[32];
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!

	TString bName = branchName;
	tr->SetBranchAddress(bName.Data(), variable);

	if (nentries == 0) nentries = tr->GetEntries();
	Info("CalibrateRawSpectra", "%lld entries from tree %s will be read.", nentries, tr->GetName());

	//make histogram to fill
	const Int_t noElements = upperElement - lowerElement +1;	//number of treated detector elements
	TH1F *hEnergy[noElements];
	for (Int_t i = 0; i < noElements; i++) {
		hEnergy[i] = 0;
	}
	TString	histName;
	TString histTitle;
	for (Int_t i = 0; i < noElements; i++) {
		histName.Form("%sE%d", bName.Data(), i+lowerElement);
		histTitle.Form("%s: %s", iFile.Data(), bName.Data());
		hEnergy[i] = new TH1F(histName.Data(), histTitle.Data(), nEBins, 0., 10.);
	}

	TRandom3 ranGen(1);

	for (Long64_t j = 0; j < nentries; j++) {
		tr->GetEntry(j);

		for (Int_t i = lowerElement; i <= upperElement; i++) {
//			printf("\n\n");
//			Info("CalculateCalibParameters", "Calculating calibration parameters for detector channel %s[%d].", block, i);
			//TH1I object preparing
//			hRaw = new TH1I("name", "title", 4096, 0, 4095);					//nastavovat hranice histogramu podle parametru fce
//			detectorChannel.Form("%s[%d]", block, i);
//			histName.Form("Hist%s[%d]", block, i);
//			hRaw->SetName(histName.Data());
//			fillCommand.Form("%s >> %s", detectorChannel.Data(), hRaw->GetName());
//			fillCondition.Form("%s > %d && %s < %d",
//					detectorChannel.Data(), lowerchannel, detectorChannel.Data(), upperchannel);
//			//filling from the .root raw data file and content arrangement
//			tr->Draw(fillCommand.Data(), fillCondition.Data(), "goff");

			//calibration of raw spectra using obtained parameters
//			Info("CalculateCalibParameters", "Energy spectrum from address %s[%d] calculating", block, i);
//			histName.Form("%sE", hRaw->GetName());
//			histTitle.Form("%s: %s", iFileName.Data(), histName.Data());
			//			detectorChannel.Form("%s[%d]", block, i);
			//			hEnergy->SetName(histName.Data());

			if (variable[i] > lowerChannel && variable[i] < upperChannel) {
				hEnergy[i-lowerElement]->Fill( fA[i]*( variable[i]+ranGen.Uniform(0., 1.0) ) + fB[i] );
			}



		}//for subaddresses

	}//for entries

	fr.Close();

	TString oFileName;
	//outputfile with calibrated spectra
	if ( (lowerElement == 0) && (upperElement == ADDRESSNUMBER-1) ) { oFileName.Form("%sE.root", bName.Data()); }
	else {
		if (lowerElement == upperElement) { oFileName.Form("%s[%d]E.root", bName.Data(), lowerElement); }
		else { oFileName.Form("%s[%d-%d]E.root", bName.Data(), lowerElement, upperElement); }

	}

	TFile fw(oFileName.Data(), "RECREATE");
	if (fw.IsOpen() == 0) {
		Error("CalculateCalibParameters", "File %s was not created and won't be processed\n\n", oFileName.Data());
		return;
	}

	fw.cd();
	for (Int_t i = 0; i < noElements; i++) {
		hEnergy[i]->Write();
	}

	fw.Close();

	return;
}

void AculCalibration::CalibrateRawSpectra(const char* inputfile, const char* block, /*const Int_t address,*/
		const char* treename, Int_t lowerchannel, Int_t upperchannel,
		Int_t nEBins, Int_t lowersubaddress, Int_t uppersubaddress) {


	//input file with raw data opening
	TString iFileName = inputfile;
	TFile *fr = new TFile(iFileName.Data());
	if ( !fr->IsOpen() ) {
		Error("CalculateCalibParameters", "File %s was not opened and won't be processed", iFileName.Data());
		return;
	}
	TTree *tr = (TTree*)fr->Get(treename);
	if (!tr) {
		Error("CalculateCalibParameters", "Tree %s was not found in file %s", treename, iFileName.Data());
		return;
	}


	TH1I *hRaw = 0;
	TH1F *hEnergy = 0;

	TRandom3 ranGen(1);


	TString	detectorChannel;
	TString	histName;
	TString histTitle;
	TString fillCommand;
	TString	fillCondition;
	TString oFileName;

	//outputfile with calibrated spectra
	if ( (lowersubaddress == 0) && (uppersubaddress == ADDRESSNUMBER-1) ) { oFileName.Form("%s[]Ecal.root", block); }
	else {
		if (lowersubaddress == uppersubaddress) { oFileName.Form("%s[%d]Ecal.root", block, lowersubaddress); }
		else { oFileName.Form("%s[%d-%d]Ecal.root", block, lowersubaddress, uppersubaddress); }

	}
	TFile *fw = new TFile(oFileName.Data(), "RECREATE");
	if (fw->IsOpen() == 0) {
		Error("CalculateCalibParameters", "File %s was not created and won't be processed\n\n", oFileName.Data());
		return;
	}

	for (Int_t i = lowersubaddress; i <= uppersubaddress; i++) {
		printf("\n\n");
		Info("CalculateCalibParameters", "Calculating calibration parameters for detector channel %s[%d].", block, i);
		//TH1I object preparing
		hRaw = new TH1I("name", "title", 4096, 0, 4095);					//nastavovat hranice histogramu podle parametru fce
		detectorChannel.Form("%s[%d]", block, i);
		histName.Form("Hist%s[%d]", block, i);
		hRaw->SetName(histName.Data());
		fillCommand.Form("%s >> %s", detectorChannel.Data(), hRaw->GetName());
		fillCondition.Form("%s > %d && %s < %d",
				detectorChannel.Data(), lowerchannel, detectorChannel.Data(), upperchannel);
		//filling from the .root raw data file and content arrangement
		tr->Draw(fillCommand.Data(), fillCondition.Data(), "goff");

		//spectrum analysis
//		fitControl = PeaksFitting(hRaw, "Q", fFitMinSigma);
//		Info("CalculateCalibParameters", "Value of fitControl is: %d", fitControl);		//ok

		//incorrectly treated spectrum output
//		if (fitControl != 0 && fCalInformation->IsOpen()) {															//ok
//			outcalfile << setw(39) << fitControl << endl;
//			fCalInformation->cd();
//			hRaw->SetLineColor(2);	//red
//			fCalInformation->cd();
//			hRaw->Write();
//			continue;
//		}//if

		//correctly treated spectrum saving
//		if (fCalInformation->IsOpen()) {
//			fCalInformation->cd();
//			hRaw->SetLineColor(kGreen+3);	//green
//			hRaw->SetFillColor(kGreen+1);
//			hRaw->Write();
//		}//if

		//calibration parameters calculation																//ok
//		for (Int_t j = 0; j < kRaNOPEAKS; j++) {		//delat podle poctu zkoumanych piku
//			calGraph->SetPoint(j, fPeak[j], fEnergy[j]);  //calibration graph filling
//		}//for
//		calGraph->Fit(calFunction, "Q", "goff", 0, 4096);  //omezit hlasitost fitovani, udelat volitelne, dodelat volby rozsahu
//		outcalfile
//		<< block << "\t"
//		<< address << "\t"
//		<< i << "\t"
//		<< setprecision(4) << calFunction->GetParameter(1) << "\t"
//		<< setprecision(4) << calFunction->GetParameter(0) << "\t\t"
//		<< fitControl
//		<< endl;
//		fAOld[address][i] = calFunction->GetParameter(1);
//		fBOld[address][i] = calFunction->GetParameter(0);


		//calibration of raw spectra using obtained parameters
		Info("CalculateCalibParameters", "Energy spectrum from address %s[%d] calculating", block, i);
		histName.Form("%sE", hRaw->GetName());
		histTitle.Form("%s: %s", iFileName.Data(), histName.Data());
		hEnergy = new TH1F(histName.Data(), histTitle.Data(), nEBins, 0., 10.);
		//			detectorChannel.Form("%s[%d]", block, i);
		//			hEnergy->SetName(histName.Data());

		for (Int_t j = lowerchannel; j < upperchannel; j++) {
			Int_t binCont = (Int_t)hRaw->GetBinContent(j);
			//				cout << j << ":\t" << hRaw->GetBinContent(j) << endl;
			for (Int_t k = 0; k < binCont; k++) {
				hEnergy->Fill( fA[i]*( j+ranGen.Uniform(0., 1.0) ) + fB[i] );
				//				cout << j << ":\t" << fAOld[address][i]*( j+ranGen.Uniform(0., 1.0) ) + fBOld[address][i] << endl;
				//				cout << j << ":\t" << fAOld[address][i] << endl;
			}
		}

		fw->cd();
		hEnergy->Write();

	}//for

}

void AculCalibration::FindEnergyPeaks(TCanvas *c1, const char* ifile, const char* outfile) {

	TString iFile = ifile;
	TFile *fr = new TFile(iFile.Data());
	if ( !fr->IsOpen() ) {
		Error("FindEnergyPeaks", "File %s was not opened and won't be processed.", iFile.Data());
		return;
	}

	TList *histList = fr->GetListOfKeys();
	Info("FindEnergyPeaks", "%d keys found in file %s.", histList->GetEntries(), fr->GetName());

	//creation of output text file with positions of peaks in MeV
	TString workFile = outfile;
	ofstream ofile;
	ofile.open(workFile.Data());
	if (!ofile.is_open()) {
		Error("PeaksFitting", "Output file %s was not opened", workFile.Data());
		return;
	}

	TH1 *hWork = 0;
	c1->Clear();
	c1->Divide(6, 6);

	for (Int_t i = 0; i < histList->GetEntries(); i++) {
		fr->GetObject(histList->At(i)->GetName(), hWork);
		c1->cd(i+1);
		PeaksFitting(hWork);
		hWork->Draw();
		ofile<<i<<"\t";
		for(Int_t j=0; j<kRaNOPEAKS; j++) {
			ofile << fPeak[j] <<"\t";
		}
		ofile<<endl;

	}

	ofile.close();
}

void AculCalibration::FindAverageEnergies(const char* ifile, const char* outfile) {

	TString iFile = ifile;
	TFile *fr = new TFile(iFile.Data());
	if ( !fr->IsOpen() ) {
		Error("FindAverageEnergies", "File %s was not opened and won't be processed.", iFile.Data());
		return;
	}

	TList *histList = fr->GetListOfKeys();
	Info("FindAverageEnergies", "%d keys found in file %s.", histList->GetEntries(), fr->GetName());

	//creation of output text file with average values of peak energies in MeV
	TString workFile = outfile;
	ofstream ofile;
	ofile.open(workFile.Data());
	if (!ofile.is_open()) {
		Error("PeaksFitting", "Output file %s was not opened", workFile.Data());
		return;
	}

	TH1 *hWork = 0;
	Double_t hArray[histList->GetEntries()][kRaNOPEAKS];
	//TString hSumName;
	Double_t hSumE1 = 0.;
	Double_t hAvrE1 = 0.;
	Double_t hSumE2 = 0.;
	Double_t hAvrE2 = 0.;
	Double_t hSumE3 = 0.;
	Double_t hAvrE3 = 0.;
	Double_t hSumE4 = 0.;
	Double_t hAvrE4 = 0.;
//	c1->Clear();
//	c1->Divide(6, 6);

	for (Int_t i = 0; i < histList->GetEntries(); i++) {
		fr->GetObject(histList->At(i)->GetName(), hWork);
		PeaksFitting(hWork);
		for(Int_t j = 0; j < kRaNOPEAKS; j++) {
			hArray[i][j] = fPeak[j];
			if(fPeak[j]==0.){
				Error("FindAverageEnergies", "No peak in channel %i !", histList->GetEntries());			
			}
			//hSumName.Form("hSumE%i",j);
		}

		hSumE1 += hArray[i][0];
		hSumE2 += hArray[i][1];
		hSumE3 += hArray[i][2];
		hSumE4 += hArray[i][3];
//		std::cout<<"i "<<i<<" hSumE1 "<<hSumE1<<std::endl;
	}
	hAvrE1 = hSumE1/histList->GetEntries();
	hAvrE2 = hSumE2/histList->GetEntries();
	hAvrE3 = hSumE3/histList->GetEntries();
	hAvrE4 = hSumE4/histList->GetEntries();
	ofile <<"Average energies are:\t"<<hAvrE1<<"\t"<<hAvrE2<<"\t"<<hAvrE3<<"\t"<<hAvrE4<<std::endl;
	ofile.close();
}


void AculCalibration::ShowAnalyzedSpectra(const char *filename, TCanvas* fittedRawCanvas, Int_t xaxismin, Int_t xaxismax, Int_t subaddress)
{

	if ( subaddress > ADDRESSNUMBER ) {
		Error("ShowAnalyzedSpectra", "Possible subaddress values have to be in range 0 - %d", ADDRESSNUMBER - 1);
		return;
	}

	if (!fittedRawCanvas) {
		Warning("ShowAnalyzedSpectra", "You have to assign TCanvas for fitted raw spectra drawing");
		return;
	}


	TFile *fr = new TFile(filename, "READ");
	if (!fr->IsOpen()) {
		cout << "File " << filename << " was not opened" << endl;
		return;
	}

	TList *histList;
	histList = fr->GetListOfKeys();
	Int_t listEntries = histList->GetEntries();
	TH1I *hDraw = 0;

	fittedRawCanvas->Clear();
//	fittedRawCanvas->SetFillColor(10);

	if ( (listEntries > 1) && (listEntries <= 8) ) {
		fittedRawCanvas->Divide(2, 4);
		fittedRawCanvas->SetFillColor(10);
	}
	if ( (listEntries > 8) && (listEntries <= 16) ) {
		fittedRawCanvas->Divide(4, 4);
		fittedRawCanvas->SetFillColor(10);
	}

	if (subaddress >= listEntries) {
		for (Int_t i = 0; i < listEntries; i++) {
			fittedRawCanvas->cd(i+1);
			fr->GetObject(histList->At(i)->GetName(), hDraw);
			if (hDraw) {
				hDraw->SetAxisRange(xaxismin, xaxismax, "X");
				hDraw->Draw();
				hDraw->SetDirectory(0);
//				if (fHAnalyzedList) {
//					fHAnalyzedList->Add(hDraw);
//				}
				fHAnalyzedList.Add(hDraw);
			}
		}//for
	}
	else {
		fr->GetObject(histList->At(subaddress)->GetName(), hDraw);
		if (hDraw) {
			hDraw->SetAxisRange(xaxismin, xaxismax, "X");
			hDraw->Draw();
			hDraw->SetDirectory(0);
			fHAnalyzedList.Add(hDraw);
		}
	}

	fr->Close();

	fittedRawCanvas->Update();

	return;

}

void AculCalibration::ShowEnergySpectra(const char *filename, TCanvas* energyCanvas, const Int_t subaddress, Option_t* option, Double_t xaxismin, Double_t xaxismax)
{
	if ( subaddress > ADDRESSNUMBER ) {
		Error("ShowEnergySpectra", "Possible subaddress values have to be in range 0 - %d", ADDRESSNUMBER - 1);
		return;
	}

	if (!energyCanvas) {
		Warning("ShowEnergySpectra", "You have to assign TCanvas for fitted raw spectra drawing");
		return;
	}

	TString	opt = option;
	opt.ToLower();

	TFile *fr = new TFile(filename, "READ");
	if (!fr->IsOpen()) {
		cout << "File " << filename << " was not opened" << endl;
		return;
	}

	TList *histList;
	histList = fr->GetListOfKeys();
	Int_t listEntries = histList->GetEntries();
	TH1F *hDraw = 0;

	energyCanvas->Clear();
	energyCanvas->SetFillColor(10);
	if ( (listEntries > 1) && (listEntries <= 8) ) {
		energyCanvas->Divide(2, 4);
		energyCanvas->SetFillColor(10);
	}
	if ( (listEntries > 8) && (listEntries <= 16) ) {
		energyCanvas->Divide(4, 4);
		energyCanvas->SetFillColor(10);
	}
	if ( (listEntries > 16) && (listEntries <= 32) ) {
		energyCanvas->Divide(6, 6);
		energyCanvas->SetFillColor(10);
	}


	if (subaddress >= listEntries) {
		if (opt.Contains("sum")) {
			energyCanvas->cd(0);
			for (Int_t i = 0; i < listEntries; i++) {
				fr->GetObject(histList->At(i)->GetName(), hDraw);
				if (hDraw) {
					hDraw->SetDirectory(0);
					//hDraw->SetAxisRange(xaxismin, xaxismax, "X");
					if (opt.Contains("c")) { hDraw->SetLineColor(i+1); }
					if (opt.Contains("c") && opt.Contains("+")) { hDraw->SetFillColor(i+1); }
//					fHEnergyStack->Add(hDraw);
					fHEnergyStack.Add(hDraw);
				}
			}

			if (opt.Contains("+")) { fHEnergyStack.Draw(); }
			else { fHEnergyStack.Draw("nostack"); }
		}
		else {
			for (Int_t i = 0; i < listEntries; i++) {
				energyCanvas->cd(i+1);
				fr->GetObject(histList->At(i)->GetName(), hDraw);
				if (hDraw) {
					hDraw->SetAxisRange(xaxismin, xaxismax, "X");
					hDraw->Draw();
					hDraw->SetDirectory(0);
//					if (fHEnergyList) {
//						fHEnergyList->Add(hDraw);
//					}
					fHEnergyList.Add(hDraw);
				}
			}//for
		}//else
	}//if
	else {
		fr->GetObject(histList->At(subaddress)->GetName(), hDraw);
		energyCanvas->cd(0);
		if (hDraw) {
			hDraw->SetAxisRange(xaxismin, xaxismax, "X");
			hDraw->Draw();
			hDraw->SetDirectory(0);
			fHEnergyList.Add(hDraw);
		}
	}//else
//energyCanvas->Draw();
	fr->Close();
//	energyCanvas->Draw();
	energyCanvas->Update();

	return;

}

Bool_t AculCalibration::FindPedestals(Int_t lowerchannel, Int_t upperchannel) {



	//input file with data opening
//	TString iFileName = inputfile;
	TFile *fr = new TFile(fInputRootFile.Data());
	if ( !fr->IsOpen() ) {
		Error("FindPedestals", "File %s was not opened and won't be processed", fInputRootFile.Data());
		return 0;
	}

	TTree *tr = (TTree*)fr->Get(treename);
	if (!tr) {
		Error("FindPedestals", "Tree %s was not found in file %s", treename, fInputRootFile.Data());
		return 0;
	}

	


	TString PedestalName;
	PedestalName.Form("pedestal_%s.par", block);
	TString oFileName = fWorkDirectory +  PedestalName;

	TString PedestalCut;
	PedestalCut.Form("PedestalCut_%s.par", block);
	TString oFileNameCut = fWorkDirectory +  PedestalCut;

	//creation of the output text file
//	if ( (lowersubaddress == 0) && (uppersubaddress == ADDRESSNUMBER-1) ) { oFileName.Form("%s[].cal", block); }
//	else {
//		if (lowersubaddress == uppersubaddress) { oFileName.Form("%s[%d].cal", block, lowersubaddress); }
//		else { oFileName.Form("%s[%d-%d].cal", block, lowersubaddress, uppersubaddress); }
//	}
	


	ofstream outcalfile;
	outcalfile.open(oFileName.Data());
	if (!outcalfile.is_open()) {
		Error("FindPedestals", "Output file %s was not opened", oFileName.Data());
		return 0;
	}

	ofstream outcalfilecut;
	outcalfilecut.open(oFileNameCut.Data());
	if (!outcalfilecut.is_open()) {
		Error("FindPedestals", "Output file %s was not opened", oFileNameCut.Data());
		return 0;
	}


	TString	detectorChannel;
	TString	histName;
	TString	funcName;
	TString histTitle;
	TString fillCommand;
	TString	fillCondition;

	TH1I *hRaw = 0;
	TF1 *fun = 0;
	Double_t* para[32];


	for (Int_t i = flowersubaddress; i <= fuppersubaddress; i++) {

		printf("\n\n");
		Info("FindPedestals", "Finding pedestal for detector channel %s[%d].", block, i);
		hRaw = new TH1I(histName.Data(), "title", upperchannel, lowerchannel, upperchannel);	
		fun = new TF1(funcName.Data(),"gaus",lowerchannel,upperchannel);			
		detectorChannel.Form("%s[%d]", block, i);
		histName.Form("Hist%s[%d]", block, i);
		funcName.Form("Func%s[%d]", block, i);
		fillCommand.Form("%s >> %s", detectorChannel.Data(), hRaw->GetName());
		fillCondition.Form("%s > %d && %s < %d",
				detectorChannel.Data(), lowerchannel, detectorChannel.Data(), upperchannel);
		//filling from the .root raw data file and content arrangement
		tr->Draw(fillCommand.Data(), fillCondition.Data(), "goff");

		hRaw->Fit(fun,"R");
		para[i] = new Double_t(3);
		fun->GetParameters(para[i]);
		outcalfile
			<< para[i][1] - 2.35482 * para[i][2]
			<< endl;

		outcalfilecut
			<< para[i][1] + 2.35482 * 6
			<< endl;
			
	}
 	

 	outcalfile.close();
 	outcalfilecut.close();
 	return 1;
}


void AculCalibration::FindEnergyPedestals()
{
	

	TString PedestalName;
	PedestalName.Form("pedestal_%s.par", block);
	TString pedestalfile =  fWorkDirectory +  PedestalName; 

	TString CalibName;
	CalibName.Form("%s.cal", block);
	TString calibparfile =  fWorkDirectory +  CalibName;  	



	TString line1;
	TString line2;

	std::ifstream infile1(calibparfile);
	std::ifstream infile2(pedestalfile);


	
	char ca[32], cb[32], ch[32];	
	Float_t a[32], b[32],fch[32],E[32],Sum, AverageEnergy, err, Err;
	Sum = 0;
	AverageEnergy = 0;
	err = 0;
	Err = 0;
	

	for(Int_t i = flowersubaddress; i < fuppersubaddress; i++) { //this loop goes from 1 to 16 or 1 to 32 i.e. number of strips

	
		line1.ReadLine(infile1);
		sscanf(line1.Data(), "%s %s", ca, cb);
		a[i] = atof(ca);
		b[i] = atof(cb);
		cout<<i<<" A = "<<a[i]<<" "<<"B = "<<b[i]<<endl;
	

		line2.ReadLine(infile2);
		sscanf(line2.Data(), "%s", ch);
		fch[i] = atof(ch);
		cout<<"pedestalchannel = "<<fch[i]<<"     "<<"EnergyPedestal ="<<a[i]+b[i]*fch[i]<<endl;	
		cout<< "_______________________________________________________________________"<<endl;	

		E [i] =a[i]+b[i]*fch[i];
		Sum += E[i];

	}
		AverageEnergy = Sum/(fuppersubaddress-flowersubaddress);


	for(Int_t i = flowersubaddress; i < fuppersubaddress; i++) { 

		err += pow((E[i]-AverageEnergy),2);

	}

		Err = sqrt(err/(fuppersubaddress-flowersubaddress-1));

		cout<< "AverageEnergy = "<< AverageEnergy<< "     "<< "Error = "<< Err <<endl;	
		cout<< "_______________________________________________________________________"<<endl;	

	infile1.close();
	infile2.close();

	

}



void AculCalibration::ClearHistograms(Option_t* option)
{
	//clear THStack and TObjArray members
	//this function will be removed as soon as possible

	TString opt = option;
	opt.ToLower();

//	fHRawList->Add(hRead);
	fHRawList.Clear();
	fHAnalyzedList.Clear();
	fHEnergyList.Clear();
	fHEnergyStack.Clear();

	return;

}

void AculCalibration::DeleteStacks(Option_t* option) {

	if (fCurrentHStack) {
		delete fCurrentHStack;
		fCurrentHStack = NULL;
	}

	fCurrentHistList.Delete();

	return;
}
//void AculCalibration::SetDeadlayer(const Double_t dl){

//	fDeadLayer = static_cast<Double_t>(dl);

//}

void AculCalibration::SetInputParameters() {

//	TString iFile = inputparfile;
	if (fParFileName.Length()==0) {
		Warning("AculCalibration::SetInputsParameters", "File with input parameters was not set.");
		return;
	}

	const Int_t lineLength = 400;
	Char_t	line[lineLength];
	Char_t	parameter[100];
	Char_t	identificator[100];


	ifstream fipr;
	fipr.open(fParFileName.Data());
	if (!fipr.is_open()) {
		Error("AculCalibration::SetInputsParameters", "File with input parameters \"%s\" was not opened.", fParFileName.Data());
		return;
	}

	Info("AculCalibration::SetInputsParameters", "File with input parameters \"%s\" will be processed.", fParFileName.Data());

	while (!fipr.eof()) {

		fipr.getline(line, lineLength);
		if (strlen(line) < 2) {
			continue;
		}

		sscanf(line, "%s %s", parameter, identificator);

		if ( strcmp(identificator, "noPeaks") == 0 ) {
			kRaNOPEAKS = static_cast<Int_t>(atoi(parameter));
			fEnergyInput.Set(kRaNOPEAKS);
			fEnergy.Set(kRaNOPEAKS);
			fPeak.Set(kRaNOPEAKS);
			for (Int_t i = 0; i < kRaNOPEAKS; i++) {
				fipr.getline(line, lineLength);
				sscanf(line, "%s", parameter);
				fEnergyInput[i] = static_cast<Double_t>(atof(parameter));
			}
			continue;
		}//if

		if ( strcmp(identificator, "lowerChannel") == 0 ) {
			sscanf(line, "%s", parameter);
			fLowerChannel = static_cast<Double_t>(atof(parameter));
		}

		if ( strcmp(identificator, "upperChannel") == 0 ) {
			sscanf(line, "%s", parameter);
			fUpperChannel = static_cast<Double_t>(atof(parameter));
		}

		if ( strcmp(identificator, "lowerPeakHight") == 0 ) {
			sscanf(line, "%s", parameter);
			fLowerPeakRelativeHight = static_cast<Double_t>(atof(parameter));
		}

		if ( strcmp(identificator, "upperPeakHight") == 0 ) {
			sscanf(line, "%s", parameter);
			fUpperPeakRelativeHight = static_cast<Double_t>(atof(parameter));
		}

		if ( strcmp(identificator, "peakPositionTolerance") == 0 ) {
			sscanf(line, "%s", parameter);
			fPeakPositionTolerance = static_cast<Double_t>(atof(parameter));
		}

		if ( strcmp(identificator, "fitFunctionLineWidth") == 0 ) {
			sscanf(line, "%s", parameter);
			fFitFuncLineWidth = static_cast<Width_t>(atoi(parameter));
		}

		if ( strcmp(identificator, "minFitSigma") == 0 ) {
			sscanf(line, "%s", parameter);
			fFitMinSigma = static_cast<Double_t>(atof(parameter));
		}

		if ( strcmp(identificator, "fitHightThreshold") == 0 ) {
			sscanf(line, "%s", parameter);
			fFitPeakThreshold = static_cast<Double_t>(atof(parameter));
		}

		if ( strcmp(identificator, "deadLayer") == 0 ) {
			sscanf(line, "%s", parameter);
			fDeadLayer = static_cast<Double_t>(atof(parameter));
		}
		if ( strcmp(identificator, "qAlphaInput") == 0 ) {
			sscanf(line, "%s", parameter);
			fQAlphaInput = static_cast<Bool_t>(atoi(parameter));
		}
		if ( strcmp(identificator, "firstAlphaParentMass") == 0 ) {
			sscanf(line, "%s", parameter);
			fFirstAlphaParentMass = static_cast<Double_t>(atof(parameter));
		}
		if ( strcmp(identificator, "lowersubaddress") == 0 ) {
			sscanf(line, "%s", parameter);
			flowersubaddress = static_cast<Int_t>(atof(parameter));
		}
		if ( strcmp(identificator, "uppersubaddress") == 0 ) {
			sscanf(line, "%s", parameter);
			fuppersubaddress = static_cast<Int_t>(atof(parameter));
		}

		if ( strcmp(identificator, "block") == 0 ) {
			sscanf(line, "%s", parameter);
			fblock = (const char*) parameter;
			block = fblock.Data();

		}

		if ( strcmp(identificator, "treename") == 0 ) {
			sscanf(line, "%s", parameter);
			ftreename = (const char*) parameter;
			treename = ftreename.Data();

		}


	}


	fipr.close();

	return;

}

void AculCalibration::DivideCanvas(TCanvas *c1, Int_t inputs) {}



