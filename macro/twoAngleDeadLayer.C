#include <vector>

#include "../AculData/AculCalibration.h"
#include "../TELoss/TELoss.h"

void twoAngleDeadLayer(
		const char* normalPeakFile,
		const char* angledPeakFile,
		const char* detectorName = "ampDSSD_L_X",
		Double_t angleDeg = 45.,
		const char* textReport = "two_angle_dead_layer_report.txt",
		const char* csvReport = "two_angle_dead_layer_report.csv",
		const char* graphFile = "two_angle_dead_layer_graphs.root")
{
	gSystem->Load("libTELoss.so");
	gSystem->Load("libAculData.so");

	AculCalibration cal;
	std::vector<Double_t> alphaEnergies;
	alphaEnergies.push_back(4.784);
	alphaEnergies.push_back(6.002);
	alphaEnergies.push_back(7.687);

	std::vector<AculTwoAngleStripResult> results = cal.AnalyzeTwoAnglePeakFiles(
			normalPeakFile,
			angledPeakFile,
			detectorName,
			angleDeg,
			alphaEnergies,
			textReport,
			csvReport);

	cal.WriteTwoAngleDiagnosticGraphs(graphFile, results);
}
