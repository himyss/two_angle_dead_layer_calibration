#include <vector>

#include "../AculData/AculCalibration.h"
#include "../TELoss/TELoss.h"

void testTwoAngleDeadLayer()
{
	gSystem->Load("libTELoss.so");
	gSystem->Load("libAculData.so");

	std::vector<Double_t> alphaEnergies;
	alphaEnergies.push_back(4.784);
	alphaEnergies.push_back(6.002);
	alphaEnergies.push_back(7.687);

	const Double_t expectedDeadLayer = 0.8;
	const Double_t expectedA = 0.006;
	const Double_t expectedB = 1.2;
	const Double_t angleDeg = 45.;
	const Double_t cosTheta = TMath::Cos(angleDeg*TMath::DegToRad());

	TELoss silicon;
	silicon.SetEL(1, 2.321);
	silicon.AddEL(14., 28.086, 1);
	silicon.SetZP(2., 4.);
	silicon.SetEtab(100000, 200.);
	silicon.SetDeltaEtab(300);

	std::vector<Double_t> normalChannels;
	std::vector<Double_t> angledChannels;
	for (size_t i = 0; i < alphaEnergies.size(); i++) {
		const Double_t normalEnergy = silicon.GetE(alphaEnergies[i], expectedDeadLayer);
		const Double_t angledEnergy = silicon.GetE(alphaEnergies[i], expectedDeadLayer/cosTheta);
		normalChannels.push_back((normalEnergy - expectedB)/expectedA);
		angledChannels.push_back((angledEnergy - expectedB)/expectedA);
	}

	AculCalibration calibration;
	AculTwoAngleStripResult result = calibration.AnalyzeTwoAngleStrip(
			normalChannels,
			angledChannels,
			alphaEnergies,
			angleDeg,
			0,
			5.,
			kFALSE);

	if (result.status != "ok" ||
			TMath::Abs(result.deadLayer - expectedDeadLayer) > 1.e-3 ||
			TMath::Abs(result.finalA - expectedA) > 1.e-6 ||
			TMath::Abs(result.finalB - expectedB) > 1.e-4) {
		Error("testTwoAngleDeadLayer", "Synthetic recovery failed: status=%s d=%f A=%f B=%f",
				result.status.Data(), result.deadLayer, result.finalA, result.finalB);
		gSystem->Exit(2);
	}

	AculTwoAngleStripResult invalid = calibration.AnalyzeTwoAngleStrip(
			normalChannels,
			angledChannels,
			alphaEnergies,
			0.,
			0,
			5.,
			kFALSE);
	if (invalid.status != "failed" || invalid.finalA != 0.) {
		Error("testTwoAngleDeadLayer", "Invalid input did not fail explicitly.");
		gSystem->Exit(3);
	}

	cout << "testTwoAngleDeadLayer: ok d=" << result.deadLayer
		<< " A=" << result.finalA
		<< " B=" << result.finalB
		<< " rms=" << result.finalRms << endl;
}
