#include <fstream>
#include <vector>

#include "../AculData/AculCalibration.h"
#include "../TELoss/TELoss.h"

namespace {

TString FileBaseName(const char* path)
{
	TString name(path);
	if (name.Contains("/")) {
		name.Remove(0, name.Last('/') + 1);
	}
	if (name.Contains(".")) {
		name.Remove(name.Last('.'));
	}
	return name;
}

TString StripText(const char* text)
{
	TString value(text ? text : "");
	return value.Strip(TString::kBoth);
}

TString StationExpression(const char* branchPrefix, const char* detectorName, Int_t strip)
{
	TString prefix = StripText(branchPrefix);
	TString expression;
	if (prefix.Length()) {
		expression.Form("%s.%s[%d]", prefix.Data(), detectorName, strip);
	}
	else {
		expression.Form("%s[%d]", detectorName, strip);
	}
	return expression;
}

AculTwoAnglePeakSet ReadChannelHistsForTwoAngle(
		const char* directory,
		const char* fileBase,
		const char* detectorName,
		Double_t angleDeg,
		Int_t firstStrip,
		Int_t lastStrip)
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

	std::vector<TFile*> files;
	std::vector<TH2D*> hists;
	for (size_t i = 0; i < peakNames.size(); i++) {
		TString fileName;
		fileName.Form("%schannelsHist_%s_%s_%s.root", dir.Data(), peakNames[i].Data(), fileBase, detectorName);
		TFile* file = TFile::Open(fileName.Data());
		if (!file || file->IsZombie() || file->GetListOfKeys()->GetEntries() == 0) {
			Error("ReadChannelHistsForTwoAngle", "Cannot read %s", fileName.Data());
			continue;
		}
		TH2D* hist = (TH2D*)file->Get(file->GetListOfKeys()->At(0)->GetName());
		if (!hist) {
			Error("ReadChannelHistsForTwoAngle", "No TH2D payload in %s", fileName.Data());
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

void DrawPeakLines(const std::vector<Double_t>& peaks, Double_t ymax, Int_t color = kRed + 1)
{
	for (size_t i = 0; i < peaks.size(); i++) {
		TLine* line = new TLine(peaks[i], 0., peaks[i], ymax);
		line->SetLineColor(color);
		line->SetLineWidth(2);
		line->Draw("same");
	}
}

void DrawSpectraWithPeakPositions(
		const char* inputFile,
		const char* treeName,
		const char* branchPrefix,
		const char* detectorName,
		const AculTwoAnglePeakSet& peakSet,
		const std::vector<AculTwoAngleStripResult>& results,
		Bool_t calibrated,
		const char* canvasTitle,
		const char* outputRootFile,
		const char* outputPngFile,
		Int_t bins,
		Double_t xMin,
		Double_t xMax,
		Double_t rawCutMin)
{
	TFile input(inputFile);
	TTree* tree = (TTree*)input.Get(treeName);
	if (!tree) {
		Error("DrawSpectraWithPeakPositions", "Tree %s was not found in %s", treeName, inputFile);
		return;
	}

	TFile output(outputRootFile, "RECREATE");
	TCanvas canvas(canvasTitle, canvasTitle, 1800, 1200);
	canvas.Divide(4, 4);

	for (Int_t strip = peakSet.firstStrip; strip < peakSet.firstStrip + static_cast<Int_t>(peakSet.stripPeaks.size()); strip++) {
		const Int_t localStrip = strip - peakSet.firstStrip;
		canvas.cd(localStrip + 1);
		gPad->SetLogy();

		TString rawExpr = StationExpression(branchPrefix, detectorName, strip);
		TString drawExpr;
		std::vector<Double_t> peakPositions;
		if (calibrated) {
			const AculTwoAngleStripResult& r = results[localStrip];
			drawExpr.Form("(%0.16g)+(%0.16g)*(%s)", r.finalB, r.finalA, rawExpr.Data());
			for (size_t i = 0; i < peakSet.stripPeaks[localStrip].size(); i++) {
				peakPositions.push_back(r.finalB + r.finalA*peakSet.stripPeaks[localStrip][i]);
			}
		}
		else {
			drawExpr = rawExpr;
			peakPositions = peakSet.stripPeaks[localStrip];
		}

		TString histName;
		histName.Form("%s_strip_%02d_%s", detectorName, strip, calibrated ? "calibrated" : "raw");
		TH1D* hist = new TH1D(histName, histName, bins, xMin, xMax);
		hist->GetXaxis()->SetTitle(calibrated ? "Energy, MeV" : "ADC channel");
		hist->GetYaxis()->SetTitle("Counts");

		TString command;
		command.Form("%s>>%s", drawExpr.Data(), histName.Data());
		TString cut;
		cut.Form("%s>%g", rawExpr.Data(), rawCutMin);
		tree->Draw(command.Data(), cut.Data(), "goff");

		hist->SetLineColor(kBlue + 1);
		hist->Draw();
		DrawPeakLines(peakPositions, hist->GetMaximum()*0.9);
		output.cd();
		hist->Write();
	}

	output.cd();
	canvas.Write();
	canvas.SaveAs(outputPngFile);
	output.Close();
	input.Close();
}

void WriteCalibrationCoefficientsLikeOneSpectrum(
		const char* outputDirectory,
		const char* fileBase,
		const char* detectorName,
		const std::vector<AculTwoAngleStripResult>& results)
{
	TString coeffFile;
	TString coeffRevFile;
	coeffFile.Form("%s/calibCoeff_%s_%s.txt", outputDirectory, fileBase, detectorName);
	coeffRevFile.Form("%s/calibCoeff_%s_%s_rev.txt", outputDirectory, fileBase, detectorName);

	std::ofstream coeff(coeffFile.Data());
	std::ofstream coeffRev(coeffRevFile.Data());
	for (size_t i = 0; i < results.size(); i++) {
		coeff << results[i].finalA << " " << results[i].finalB << "\n";
		coeffRev << results[i].finalB << " " << results[i].finalA << "\n";
	}
}

void WriteReportLikeOneSpectrum(
		const char* outputDirectory,
		const char* normalFileBase,
		const char* angledFileBase,
		const char* detectorName,
		const AculTwoAnglePeakSet& normal,
		const AculTwoAnglePeakSet& angled,
		const std::vector<AculTwoAngleStripResult>& results)
{
	TString reportFileName;
	TString coeffFileName;
	TString coeffFileNameRev;
	TString deadFileName;
	reportFileName.Form("%s/report_%s_%s.txt", outputDirectory, normalFileBase, detectorName);
	coeffFileName.Form("calibCoeff_%s_%s.txt", normalFileBase, detectorName);
	coeffFileNameRev.Form("calibCoeff_%s_%s_rev.txt", normalFileBase, detectorName);
	deadFileName.Form("two_angle_dead_layer_graphs.root");

	std::ofstream report(reportFileName.Data());
	report << "Calibration results. " << std::endl << std::endl;
	report << "Two-angle dead-layer algorithm parameters: " << std::endl;
	report << "  normal measurement: " << normalFileBase << std::endl;
	report << "  angled measurement: " << angledFileBase << std::endl;
	report << "  angle: " << angled.angleDeg << " deg" << std::endl;
	report << "  detector: " << detectorName << std::endl;
	report << "  silicon energy-loss model: TELoss" << std::endl;
	report << "  calibration energies [MeV]: ";
	for (size_t i = 0; i < normal.alphaEnergies.size(); i++) {
		report << normal.alphaEnergies[i] << (i + 1 == normal.alphaEnergies.size() ? "" : " ");
	}
	report << std::endl << std::endl;

	report << "ADC-channels counts are stored in files: " << std::endl;
	std::vector<TString> peakNames;
	peakNames.push_back("low");
	peakNames.push_back("middle");
	peakNames.push_back("high");
	for (size_t iPeak = 0; iPeak < peakNames.size() && iPeak < normal.alphaEnergies.size(); iPeak++) {
		TString meanFileName;
		meanFileName.Form("channelsHist_%s_%s_%s.root", peakNames[iPeak].Data(), normalFileBase, detectorName);
		report << " * " << meanFileName << " --- for " << peakNames[iPeak] << " energy peak;" << std::endl;
		for (size_t iStrip = 0; iStrip < normal.stripPeaks.size(); iStrip++) {
			const Int_t strip = normal.firstStrip + static_cast<Int_t>(iStrip);
			report << "\tStrip " << strip << ": " << normal.stripPeaks[iStrip][iPeak] << std::endl;
		}
	}
	report << std::endl;

	report << "Angled ADC-channels counts are stored in files: " << std::endl;
	for (size_t iPeak = 0; iPeak < peakNames.size() && iPeak < angled.alphaEnergies.size(); iPeak++) {
		TString meanFileName;
		meanFileName.Form("channelsHist_%s_%s_%s.root", peakNames[iPeak].Data(), angledFileBase, detectorName);
		report << " * " << meanFileName << " --- for " << peakNames[iPeak] << " energy peak;" << std::endl;
		for (size_t iStrip = 0; iStrip < angled.stripPeaks.size(); iStrip++) {
			const Int_t strip = angled.firstStrip + static_cast<Int_t>(iStrip);
			report << "\tStrip " << strip << ": " << angled.stripPeaks[iStrip][iPeak] << std::endl;
		}
	}
	report << std::endl;

	report << "Dead layer estimation [um] by strips are stored in file: " << deadFileName << ":" << std::endl;
	Double_t sumDeadLayer = 0.;
	Int_t goodStrips = 0;
	for (size_t i = 0; i < results.size(); i++) {
		report << "\tStrip " << results[i].strip << ": " << results[i].deadLayer;
		if (results[i].status != "ok") {
			report << " (" << results[i].status << ")";
		}
		report << std::endl;
		if (results[i].status == "ok" || results[i].status == "fallback") {
			sumDeadLayer += results[i].deadLayer;
			goodStrips++;
		}
	}
	report << std::endl;
	report << "Effective dead layer thickness is: "
		<< (goodStrips ? sumDeadLayer / goodStrips : 0.)
		<< " [um]." << std::endl;

	report << "Calibration coefficients are stored in files: " << std::endl;
	report << " * " << coeffFileName
		<< " --- linear coefficient (a) in the first column,"
		<< " free coefficient (b) in the second" << std::endl;
	report << " * " << coeffFileNameRev << " --- reversed columns" << std::endl;

	report << std::setw(20) << " b " << std::setw(10) << " a " << std::endl;
	for (size_t i = 0; i < results.size(); i++) {
		report << "\tStrip " << results[i].strip << ": "
			<< results[i].finalB << " "
			<< results[i].finalA
			<< std::endl;
	}
	report.close();
}

} // namespace

void twoAngleFullCalibrationReport(
		const char* normalRawFile = "",
		const char* angledRawFile = "",
		const char* detectorName = "ampDSSD_L_X",
		Double_t angleDeg = 45.,
		const char* treeName = "AnalysisxTree",
		const char* branchPrefix = "NeEvent",
		const char* channelHistRoot = "results",
		const char* outputDirectory = "results/two_angle_full",
		Int_t firstStrip = 0,
		Int_t lastStrip = 15)
{
	gSystem->Load("libTELoss.so");
	gSystem->Load("libAculData.so");

	if (!normalRawFile || !normalRawFile[0] || !angledRawFile || !angledRawFile[0]) {
		Error("twoAngleFullCalibrationReport",
				"Pass two raw ROOT files explicitly: twoAngleFullCalibrationReport(\"normal.root\", \"angled.root\", ...)");
		return;
	}

	gSystem->mkdir(outputDirectory, kTRUE);

	TString normalBase = FileBaseName(normalRawFile);
	TString angledBase = FileBaseName(angledRawFile);
	TString normalDir;
	TString angledDir;
	normalDir.Form("%s/%s/%s", channelHistRoot, normalBase.Data(), detectorName);
	angledDir.Form("%s/%s/%s", channelHistRoot, angledBase.Data(), detectorName);

	AculTwoAnglePeakSet normal = ReadChannelHistsForTwoAngle(
			normalDir.Data(), normalBase.Data(), detectorName, 0., firstStrip, lastStrip);
	AculTwoAnglePeakSet angled = ReadChannelHistsForTwoAngle(
			angledDir.Data(), angledBase.Data(), detectorName, angleDeg, firstStrip, lastStrip);

	normal.alphaEnergies.push_back(4.784);
	normal.alphaEnergies.push_back(6.002);
	normal.alphaEnergies.push_back(7.687);
	angled.alphaEnergies = normal.alphaEnergies;

	AculCalibration calibration;
	std::vector<AculTwoAngleStripResult> results =
		calibration.AnalyzeTwoAnglePreparedPeaks(normal, angled);

	TString textReport;
	TString csvReport;
	TString graphFile;
	textReport.Form("%s/two_angle_dead_layer_report.txt", outputDirectory);
	csvReport.Form("%s/two_angle_dead_layer_report.csv", outputDirectory);
	graphFile.Form("%s/two_angle_dead_layer_graphs.root", outputDirectory);
	calibration.WriteTwoAngleDiagnostics(textReport.Data(), csvReport.Data(), results, normal, angled);
	calibration.WriteTwoAngleDiagnosticGraphs(graphFile.Data(), results);

	WriteCalibrationCoefficientsLikeOneSpectrum(outputDirectory, normalBase.Data(), detectorName, results);
	WriteReportLikeOneSpectrum(
			outputDirectory,
			normalBase.Data(),
			angledBase.Data(),
			detectorName,
			normal,
			angled,
			results);

	TString rawNormalRoot;
	TString rawAngledRoot;
	TString rawNormalPng;
	TString rawAngledPng;
	TString calNormalRoot;
	TString calAngledRoot;
	TString calNormalPng;
	TString calAngledPng;
	rawNormalRoot.Form("%s/raw_%s_%s.root", outputDirectory, normalBase.Data(), detectorName);
	rawAngledRoot.Form("%s/raw_%s_%s.root", outputDirectory, angledBase.Data(), detectorName);
	rawNormalPng.Form("%s/raw_%s_%s.png", outputDirectory, normalBase.Data(), detectorName);
	rawAngledPng.Form("%s/raw_%s_%s.png", outputDirectory, angledBase.Data(), detectorName);
	calNormalRoot.Form("%s/calibrated_%s_%s.root", outputDirectory, normalBase.Data(), detectorName);
	calAngledRoot.Form("%s/calibrated_%s_%s.root", outputDirectory, angledBase.Data(), detectorName);
	calNormalPng.Form("%s/calibrated_%s_%s.png", outputDirectory, normalBase.Data(), detectorName);
	calAngledPng.Form("%s/calibrated_%s_%s.png", outputDirectory, angledBase.Data(), detectorName);

	DrawSpectraWithPeakPositions(
			normalRawFile, treeName, branchPrefix, detectorName, normal, results,
			kFALSE, "normal raw spectra with peaks",
			rawNormalRoot.Data(), rawNormalPng.Data(), 4096, 0., 4096., 100.);
	DrawSpectraWithPeakPositions(
			angledRawFile, treeName, branchPrefix, detectorName, angled, results,
			kFALSE, "angled raw spectra with peaks",
			rawAngledRoot.Data(), rawAngledPng.Data(), 4096, 0., 4096., 100.);
	DrawSpectraWithPeakPositions(
			normalRawFile, treeName, branchPrefix, detectorName, normal, results,
			kTRUE, "normal calibrated spectra with peaks",
			calNormalRoot.Data(), calNormalPng.Data(), 1200, 4., 8.2, 100.);
	DrawSpectraWithPeakPositions(
			angledRawFile, treeName, branchPrefix, detectorName, angled, results,
			kTRUE, "angled calibrated spectra with peaks",
			calAngledRoot.Data(), calAngledPng.Data(), 1200, 4., 8.2, 100.);
}
