#include <vector>

#include "../AculData/AculCalibration.h"
#include "../TELoss/TELoss.h"

AculTwoAnglePeakSet ReadTwoAngleChannelHists(
		const char* directory,
		const char* fileBase,
		const char* detectorName,
		Double_t angleDeg,
		Int_t firstStrip = 0,
		Int_t lastStrip = 15)
{
	AculTwoAnglePeakSet peakSet;
	peakSet.measurementName = fileBase;
	peakSet.detectorName = detectorName;
	peakSet.angleDeg = angleDeg;
	peakSet.firstStrip = firstStrip;

	TString dir(directory);
	if (!dir.EndsWith("/")) {
		dir += "/";
	}

	std::vector<TString> peakNames;
	peakNames.push_back("low");
	peakNames.push_back("middle");
	peakNames.push_back("high");

	std::vector<TH2D*> hists;
	std::vector<TFile*> files;
	for (size_t i = 0; i < peakNames.size(); i++) {
		TString fileName;
		fileName.Form("%schannelsHist_%s_%s_%s.root", dir.Data(), peakNames[i].Data(), fileBase, detectorName);
		TFile* file = TFile::Open(fileName.Data());
		if (!file || file->IsZombie() || file->GetListOfKeys()->GetEntries() == 0) {
			Error("ReadTwoAngleChannelHists", "Cannot read %s", fileName.Data());
			continue;
		}
		TH2D* hist = (TH2D*)file->Get(file->GetListOfKeys()->At(0)->GetName());
		if (!hist) {
			Error("ReadTwoAngleChannelHists", "No TH2D payload in %s", fileName.Data());
			continue;
		}
		files.push_back(file);
		hists.push_back(hist);
	}

	if (hists.size() != peakNames.size()) {
		for (size_t i = 0; i < files.size(); i++) {
			files[i]->Close();
		}
		return peakSet;
	}

	for (Int_t strip = firstStrip; strip <= lastStrip; strip++) {
		std::vector<Double_t> peaks;
		for (size_t i = 0; i < hists.size(); i++) {
			peaks.push_back(hists[i]->GetBinContent(1, strip + 1));
		}
		peakSet.stripPeaks.push_back(peaks);
	}

	for (size_t i = 0; i < files.size(); i++) {
		files[i]->Close();
	}
	return peakSet;
}

void twoAngleDeadLayerFromChannelHists(
		const char* normalDirectory,
		const char* normalFileBase,
		const char* angledDirectory,
		const char* angledFileBase,
		const char* detectorName = "ampDSSD_L_X",
		Double_t angleDeg = 45.,
		Int_t firstStrip = 0,
		Int_t lastStrip = 15,
		const char* textReport = "two_angle_dead_layer_report.txt",
		const char* csvReport = "two_angle_dead_layer_report.csv",
		const char* graphFile = "two_angle_dead_layer_graphs.root")
{
	gSystem->Load("libTELoss.so");
	gSystem->Load("libAculData.so");

	AculTwoAnglePeakSet normal = ReadTwoAngleChannelHists(
			normalDirectory,
			normalFileBase,
			detectorName,
			0.,
			firstStrip,
			lastStrip);
	AculTwoAnglePeakSet angled = ReadTwoAngleChannelHists(
			angledDirectory,
			angledFileBase,
			detectorName,
			angleDeg,
			firstStrip,
			lastStrip);

	normal.alphaEnergies.push_back(4.784);
	normal.alphaEnergies.push_back(6.002);
	normal.alphaEnergies.push_back(7.687);
	angled.alphaEnergies = normal.alphaEnergies;

	AculCalibration calibration;
	std::vector<AculTwoAngleStripResult> results =
		calibration.AnalyzeTwoAnglePreparedPeaks(normal, angled);
	calibration.WriteTwoAngleDiagnostics(textReport, csvReport, results, normal, angled);
	calibration.WriteTwoAngleDiagnosticGraphs(graphFile, results);
}
