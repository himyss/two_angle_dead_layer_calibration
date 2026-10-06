#include "AculCalibCsI.h"

ClassImp(AculCalibCsI);

AculCalibCsI::AculCalibCsI() : fr("TFile", 5), tr("TTree", 5), cutsCol("TCutG"), fA(16), fB(16) {

	printf("AculCalibCsI::Default constructor called.\n");


}

AculCalibCsI::AculCalibCsI(const char* parfile) : fr("TFile", 5), tr("TTree", 5), cutsCol("TCutG", 10), fA(16), fB(16) {
	printf("AculCalibCsI::Constructor called.\n");

	SetParFile(parfile);
	SetPars();
	cout << "nofiles: " << nofiles << endl;
	OpenTrees();
	LoadCuts();

//	fr.Print();
//	fr.At(0);
	

}

AculCalibCsI::~AculCalibCsI() {

	printf("AculCalibCsI::Destructor called.\n");

}

void AculCalibCsI::OpenTrees() {

	TFile *file;
	for (Int_t i = 0; i < nofiles; i++) {	
		fr[i] = new TFile(fileNames[i], "READ");
		file = (TFile*)fr.At(i);
		cout << file->GetName() << endl;
		tr[i] = (TTree*)file->Get("AnalysisxTree");
	}
}

void AculCalibCsI::LoadCuts() {

	if (cutsFileName == "") {
		printf("AculCalibCsI::LoadCuts: Name of file (*.root) with cuts was not provided.\n");
		printf("AculCalibCsI::LoadCuts: No cuts has been loaded.\n");
		return;
	}
	fCuts = new TFile(cutsFileName.Data(), "READ");
	for (Int_t i = 0; i < nCuts; i++) {
		cutsCol[i] = (TCutG*)fCuts->Get(cutNames[i]);
	}

}

void AculCalibCsI::SetParFile(const char* parfile) {

	fParFileName = parfile;
	return;
}

void AculCalibCsI::PrintParameters(const char* option) {
	
	TString opt = option;

	printf("Energy points: %d\n", energyPoints);
//	printf("Detector: %s, Particle: %s\n", ????);
	for (Int_t i = 0; i < energyPoints; i++) {	
		printf("Peak ranges for energy %f:\n", fE[i]);
		if ( opt.Contains("v") ) {
			for (Int_t j = 0; j < 16; j++) {
				printf("ch: %d\tmin: %d\tmax: %d\n", j, peakMin[i][j], peakMax[i][j]);
			}
		}//if
	}//for

}


void AculCalibCsI::PrintPeakRanges() {
	for (Int_t i = 0; i < 16; i++) {
		printf("%d\t%d\t%d\n", i, peakMin[0][i], peakMax[0][i]);
	}
}

void AculCalibCsI::PrintTrees() {

	TTree *curTree = 0;
	for (Int_t i = 0; i < tr.GetEntries(); i++) {
		curTree = (TTree*)tr.At(i);
		if (curTree) {
			printf("Tree No. %d; File: %s; Name: %s\n", i, curTree->GetDirectory()->GetName(), curTree->GetName());
		} else {
			printf("Tree No. %d was not loaded. Maximal number of trees is %d\n", i, NOCALFILES);
		}
	}

	return;
}

void AculCalibCsI::PrintFiles() {

	printf("Number of loaded files: %d\n", fr.GetEntries());


	TFile *curFile = 0;
	for (Int_t i = 0; i < fr.GetEntries(); i++) {
		curFile = (TFile*)fr.At(i);
		if (curFile) {
			printf("File No. %d: \"%s\"\n", i, curFile->GetName());
		}
//		else {
//			printf("File No. %d was not loaded. Maximal number of files is %d\n", i, NOCALFILES);
//		}
	}

	return;
}

void AculCalibCsI::PrintCuts() {
	
//	printf("AculCalibCsI::PrintCuts: works wrong\n");
//	return;
	if (fCuts) printf("Cuts loaded from file \"%s\"\n", fCuts->GetName());

	TCutG *curCut = 0;
	//cutsCol

	for (Int_t i = 0; i < cutsCol.GetEntries(); i++) {
		curCut = (TCutG*)cutsCol.At(i);
		if (curCut) {
			printf("Cut No. %d; Name: \"%s\"\n", i, curCut->GetName());
		} /*else {
			printf("TOF cut No. %d was not loaded. Maximal number of cuts is %d\n", i, NOCALFILES);
		}*/
	}

	return;
}

Double_t AculCalibCsI::GetA(Int_t i) {
	if (i >= fA.GetSize()) //if i >= number of array element
	{
		return 0.;
	}
	return fA[i];
}

Double_t AculCalibCsI::GetB(Int_t i) {
	if (i >= fB.GetSize()) //if i >= number of array element
	{
		return 0.;
	}
	return fB[i];
}

void AculCalibCsI::DrawVariable(const char* variable, Int_t tree, TCanvas *canvas, Int_t lowRange, Int_t upRange) {

//	if (!canvas) TCanvas *c = new TCanvas();
	if (!canvas) return;

	canvas->Clear();
	canvas->Divide(4,4);
	
	TString canvasTitle;
	TString var;
	TString con;

	TTree *curTree = 0;
	curTree = (TTree*)tr.At(tree);
	if (!curTree) {
		printf("AculCalibCsI::DrawVariable: Tree No. %d was not found.\n", tree);
		return;
	}

	canvasTitle.Form("variable: %s; tree: %d", variable, tree);	
	canvas->SetTitle(canvasTitle.Data());

	for (Int_t i = 0; i < 16; i++) {
		var.Form("%s[%d]", variable, i);
		con.Form("%s[%d]>%d && %s[%d]<%d", variable, i, lowRange, variable, i, upRange);
		canvas->cd(i+1);
		curTree->Draw(var.Data(), con.Data());
		canvas->Update();
	}

}

void AculCalibCsI::DrawBeam(TCanvas *canvas, Int_t files, const char* variable) {

	canvas->SetTitle("Beam");
	canvas->Clear();
//	canvas->Divide(files, 2);
	canvas->Divide(files, 3);

	TTree *curTree = 0;
	TCutG *curCutG = 0;
	TString var;
	TString con;
	
	for (Int_t i = 0; i < files; i++) {
		canvas->cd(i+1);
		curTree = (TTree*)tr.At(i);
		if (!curTree) {
			printf("AculCalibCsI::DrawBeam: Tree No. %d was not found.\n", i);
			continue;
		}
		curTree->Draw("QDC[0]:TDC[0]", "TDC[0]<1000 && QDC[0]<2000", "cont");
//		curTree->Draw("QDC[0]+QDC[1]:(TDC[2]+TDC[3])/2. - (TDC[0]+TDC[1])/2.", "", "cont");
//		cout << "aksjda\t" << energyPoints << endl;
		for (Int_t j = 0; j < energyPoints; j++) {		
			if ( cutsCol.At(j) ) {
				curCutG = (TCutG*)cutsCol.At(j);
				curCutG->Draw("same");
//				printf("AculCalibCsI::DrawBeam: cTOF cut No. %d cannot be drawn, need to repair this function.\n", j);
			}
		}

		canvas->cd(files+1+i);
		curTree->Draw("QDC[0]:QDC[1]", "", "cont");
		for (Int_t j = 0; j < energyPoints; j++) {		
			if ( cutsCol.At(files + j) ) {
				curCutG = (TCutG*)cutsCol.At(files + j);
				curCutG->Draw("same");
//				printf("AculCalibCsI::DrawBeam: cQCD cut No. %d cannot be drawn, need to repair this function.\n", j);
			}
		}

		canvas->cd(2*files+1+i);
		var.Form("%s[5]:TDC[0]", variable);
		con.Form("%s[5]>200", variable);
		curTree->Draw(var.Data(), con.Data(), "cont");
		for (Int_t j = 0; j < energyPoints; j++) {		
			if ( cutsCol.At(files + j) ) {
				curCutG = (TCutG*)cutsCol.At(files + j);
				curCutG->Draw("same");
//				printf("AculCalibCsI::DrawBeam: cQCD cut No. %d cannot be drawn, need to repair this function.\n", j);
			}
		}
		canvas->Update();
	}
}

void AculCalibCsI::DrawVariableCut(const char* variable, Int_t tree, TCanvas *canvas, const char* cut1, const char* cut2, Int_t lowRange) {

//	cout << "kjashbdfjka ajsdbf jakhsdb askjdhb" << endl;

	if (!cutsCol.FindObject(cut1)) {
		printf("Cut %s was not found.\n", cut1);
		return;
	}

	canvas->Clear();
	canvas->Divide(4,4);

	TString canvasTitle;
	TString var;
	TString con;
	TString sVariable = variable;
	sVariable.ToLower();
	Int_t channel = 0;

//	cout << cutsCol.FindObject(cut1)->GetName() << endl;

	canvasTitle.Form("variable: %s; cut1: %s; cut2: %s; tree: %d", variable, cutsCol.FindObject(cut1)->GetName(), cut2 , tree);	
	canvas->SetTitle(canvasTitle.Data());

	TTree *curTree = 0;
	curTree = (TTree*)tr.At(tree);
	if (!curTree) {
		printf("AculCalibCsI::DrawVariableCut: Tree No. %d was not found.\n", tree);
		return;
	}

	for (Int_t i = 0; i<16; i++) {
		if (sVariable.Contains("anc")) {
			channel = i+1;
		} else {
			channel = i;
		}
//		cout << channel << endl;
		var.Form("%s[%d]", variable, channel);
		con.Form("%s[%d]>%d", variable, channel, lowRange);

		canvas->cd(i+1);
		curTree->SetLineColor(1);
		curTree->Draw(var.Data(), con.Data());
		curTree->SetLineColor(3);
		con.Form("%s[%d]>0 && %s", variable, channel, cutsCol.FindObject(cut1)->GetName());
		curTree->Draw(var.Data(), con.Data(), "same");

		TString scut2 = cut2;
		if (scut2.Length() != 0) {
			if (!cutsCol.FindObject(cut2)) {
				printf("Cut %s was not found.\n", cut2);
				continue;
			}
			curTree->SetLineColor(4);
			con.Form("%s[%d]>0 && %s", variable, channel, cutsCol.FindObject(cut2)->GetName());
			curTree->Draw(var.Data(), con.Data(), "same");
		}
		canvas->Update();
	}	

}

void AculCalibCsI::GetPeakMean(const char* variable, Int_t tree, Int_t energy, TCanvas *canvas, const char* beamcut, const Int_t nbins, Int_t lowRange) {

	canvas->Clear();
	canvas->Divide(4,4);

	TString var;
	TString con;
	TString hname;
	TString canvasTitle;

	canvasTitle.Form("variable: %s; tree: %d; cut: %s;", variable, tree, cutsCol.FindObject(beamcut)->GetName());	
	canvas->SetTitle(canvasTitle.Data());

	TTree *curTree = 0;
	curTree = (TTree*)tr.At(tree);
	if (!curTree) {
		printf("AculCalibCsI::GetPeakMean: Tree No. %d was not found.\n", tree);
		return;
	}

	TString sVariable = variable;
	sVariable.ToLower();
	Int_t channel = 0;

	for (Int_t i = 0; i<16; i++) {

		if (sVariable.Contains("anc")) {
			channel = i+1;
		} else {
			channel = i;
		}		
		
		var.Form("%s[%d]>>hfull[%d][%d]", variable, channel, tree, i);
		con.Form("%s[%d]>%d && %s", variable, channel, lowRange, cutsCol.FindObject(beamcut)->GetName());
		canvas->cd(i+1);
		hname.Form("hfull[%d][%d]", tree, i);
		hfull[tree][i] = new TH1I(hname.Data(), "title1", nbins, 0, 4096);
		curTree->SetLineColor(1);
		curTree->Draw(var.Data(), con.Data());

		var.Form("%s[%d]>>hcut[%d][%d]", variable, channel, tree, i);
		con.Form("%s[%d]>%d && %s[%d]<%d && %s", variable, channel, peakMin[energy][i], variable, channel, peakMax[energy][i], cutsCol.FindObject(beamcut)->GetName());
		hname.Form("hcut[%d][%d]", tree, i);
		hcut[tree][i] = new TH1I(hname.Data(), "title2", nbins, 0, 4096);
		hcut[tree][i]->SetLineColor(3);
		curTree->Draw(var.Data(), con.Data(), "same");

		gPad->Update();
		mean[tree][i] = hcut[tree][i]->GetMean();
		meanRMS[tree][i] = hcut[tree][i]->GetRMS();
//		cout << meanRMS[tree][i] << endl;
//		cout << hcut[tree][i]->GetMean() << "\t" << hcut[tree][i]->GetRMS() << endl << endl;
		
		
	}	

	canvas->Update();

}

void AculCalibCsI::Calibrate(TCanvas *canvas, Bool_t savefile, const char* filename, const char* option) {

	canvas->Clear();
	canvas->Divide(4,4);

//	cout << alphas2.GetSize()+1 << endl;
	cout <<" energyPoints+1 "<<energyPoints+1<< endl;

	TString gName;
	TString gTitle;

//	if (savefile) fGraphs->Open(filename, "RECREATE");

	TF1 *fnc;

	for (Int_t i = 0; i<16; i++) {
		canvas->cd(i+1);
		gCal[i] = new TGraphErrors(energyPoints+1);
//		FillGraph(gCal[i], energies.GetSize()+1, energies.GetArray(), i);
		FillGraph(gCal[i], energyPoints+1, fE.GetArray(), i);
//		if (savefile) gCal[i]->Write();
//		gCal[i]->Draw("Al*");
		gCal[i]->Draw("A*");
		gName.Form("g%s%s%d\n", detName.Data(), partName.Data(), i);
		gTitle.Form("%s %s\n", detName.Data(), partName.Data());
//		gCal[i]->SetTitle(gTitle.Data());
		gCal[i]->SetName(gName.Data());
		gCal[i]->Fit("pol1");
		fnc = gCal[i]->GetFunction("pol1");
		fnc->SetLineColor(kRed);
		fA[i] = fnc->GetParameter(1);
		fB[i] = fnc->GetParameter(0);
		canvas->Update();
	}
	if (savefile) SaveClbGraphs(filename, option);

}

void AculCalibCsI::FillGraph(TGraphErrors *g, Int_t npoints, Double_t *energies, Int_t graphNumber, const char* option) {

	TString opt = option;

	//all available energy points and (0,0)
	g->SetPoint(0, 0., 0.);	//(point number, x, y)
	for (Int_t i = 0; i < npoints-1; i++) {
//		g->SetPoint(i+1, energies[i], mean[i][graphNumber]);
//		g->SetPointError(i+1, 0, meanRMS[i][graphNumber]);
		g->SetPoint(i+1, mean[i][graphNumber], energies[i]);
		g->SetPointError(i+1, meanRMS[i][graphNumber], 0);
	}

//	for (Int_t j = 1; j < 4; j++) {
//			g->SetPoint(j, 0., 0.);
//	}
	

}

void AculCalibCsI::WriteClbParameters(const char* filename) {

	std::ofstream outfile(filename);
	if ( !outfile.is_open() ) {
		printf("AculCalibCsI::WriteClbParameters: File %s was not open.\n", filename);
		return;
	}	
	
	outfile << "#detector:\t" << detName << ",\tparticle:\t" << partName << endl;
	outfile << "#channel\tfA\tfB" << endl;
	for (Int_t i = 0; i < 16; i++) {
		outfile << i << "\t" << fA[i] << "\t" << fB[i] << endl;
	}
}

void AculCalibCsI::SaveClbGraphs(const char* filename, const char* option) {

	cout << "asdasd" << endl;	
	cout << fGraphs << endl;
	if (fGraphs) fGraphs->Close();
	cout << "asdasd" << endl;
	fGraphs = new TFile(filename, option);
	cout << fGraphs->IsOpen() << endl;
		cout << fGraphs->GetName() << endl;
	fGraphs->Print();
//	if (!fGraphs->IsOpen()) {
//		printf("AculCalibCsI::SaveClbGraphs: file %s was not open.\n", filename);
//		return;
//	}

	for (Int_t i = 0; i<16; i++) {
		fGraphs->cd();
		gCal[i]->Write();
	}
	fGraphs->Close();
	return;
}

void AculCalibCsI::ReadClbParameters(const char* filename) {

	std::ifstream infile(filename);
	if ( !infile.is_open() ) {
		printf("AculCalibCsI::ReadClbParameters: File %s was not open.\n", filename);
		return;
	}	

	TString line;
	Int_t i;
	Double_t a, b;

	while (!infile.eof()) {
		line.ReadLine(infile);

		if ( line.BeginsWith("#") ) continue;

//		cout << line.Data() << endl;
		sscanf(line.Data(), "%d %lf %lf", &i, &a, &b);
		fA[i] = a;
		fB[i] = b;
//		printf("fA[%d]: %f,\tfB[%d]: %f\n", i, fA[i], i, fB[i]);
//		printf("fA[%d]: %f,\tfB[%d]: %f\n\n", i, a, i, b);
	}

}

void AculCalibCsI::DrawClbGraphs(const char* filename, const char* graphname, TCanvas *canvas) {

	printf("AculCalibCsI::DrawClbGraphs: does not work\n");
	return;

	TFile gfile(filename);
	
	TString gName;
//	TGraph *gr;

	for (Int_t i = 0; i < 16; i++) {
		gName.Form("%s%d", graphname, i);
		gCal[i] = (TGraphErrors*)gfile.Get(gName.Data());
		canvas->cd(i+1);
		gCal[i]->Draw("A*");
	}
}

void AculCalibCsI::DrawEnergyDeposite(const char* variable, TCanvas *canvas, Int_t tree, const char* option) {
	if (!canvas) return;

	canvas->Divide(4,4);
	
	TString opt;
	opt = option;
	opt.ToLower();

	TString canvasTitle;
	TString var;
	TString con;

	TTree *curTree = 0;
	curTree = (TTree*)tr.At(tree);
	if (!curTree) {
		printf("AculCalibCsI::DrawVariable: Tree No. %d was not found.\n", tree);
		return;
	}

	canvasTitle.Form("variable: %s [MeV]; tree: %d", variable, tree);	
	canvas->SetTitle(canvasTitle.Data());

	if (!opt.Contains("same")) canvas->Divide(4,4);
	for (Int_t i = 0; i < 16; i++) {
		var.Form("%s[%d]*%f + %f", variable, i, fB[i], fA[i]);
//		var.Form("%s[%d]*%f + %f", variable, i, fA[i], fB[i]);
		con.Form("%s[%d]>0", variable, i);
		
		if (opt.Contains("same")) {
			if (i==0) curTree->Draw(var.Data(), con.Data());			
			curTree->Draw(var.Data(), con.Data(), "same");
		}
		else {
			canvas->cd(i+1);
			curTree->Draw(var.Data(), con.Data());
		}
		canvas->Update();
	}

	
}

void AculCalibCsI::SetPars() {

	std::ifstream infile(fParFileName.Data());
	if ( !infile.is_open() ) {
		printf("AculCalibCsI::ReadParFile: File %s was not open.\n", fParFileName.Data());
		return;
	}

	energyPoints = 0;

	TString line;

	Int_t i, min, max;
	char det[100], part[100], fname[500], cname[100];

	double en;	//energy

	while (!infile.eof()) {
		line.ReadLine(infile);

		if ( line.BeginsWith("#") ) continue;

		if ( line.BeginsWith("energies") ) {
			sscanf(line.Data(), "%*s %d", &i);
			fE.Set(i);
			continue;
		}

		if ( line.BeginsWith("files") ) {
			sscanf(line.Data(), "%*s %d", &i);
			nofiles = i;
			printf("AculCalibCsI::ReadParFile: %d files will be loaded:\n", nofiles);
			for (Int_t j = 0; j < nofiles; j++) {
				line.ReadLine(infile);
				sscanf(line, "%s", fname);
				fileNames[j] = fname;
				cout << fileNames[j] << endl;
			}
			continue;
		}

		if ( line.BeginsWith("cutFile") ) {
			sscanf(line.Data(), "%*s %s %d", fname, &nCuts);
			cutsFileName = fname;
			for (Int_t j = 0; j < nCuts; j++) {
				line.ReadLine(infile);
				sscanf(line, "%s", cname);
				cutNames[j] = cname;
				cout << cutNames[j] << endl;
			}
			continue;
		}

		if ( line.BeginsWith("detector") ) {
			sscanf(line.Data(), "%*s %s %s", det, part);
			detName = det;
			partName = part;
			printf("%s %s\n", detName.Data(), partName.Data());
			continue;
		}

		if ( line.BeginsWith("energy") ) {
			sscanf(line.Data(), "%*s %lf", &en);
			fE[energyPoints] = en;
			printf("%f\n", fE[energyPoints]);
			energyPoints++;
			continue;
		}

		sscanf(line.Data(), "%d %d %d", &i, &min, &max);
		peakMin[energyPoints-1][i] = min;
		peakMax[energyPoints-1][i] = max;
//		printf("%d %d %d\n", i, peakMin[energyPoints-1][i], peakMax[energyPoints-1][i]);
	}//while

	printf("energyPoints: %d\n", energyPoints);

	infile.close();
	return;
}
