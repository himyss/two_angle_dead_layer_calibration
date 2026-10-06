#pragma once
#include <TROOT.h>
#include <TObject.h>
#include <TSpline.h>
#include <TArrayD.h>
#include <TMath.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
using std::cout;
using std::endl;

class TELoss : public TObject
{
private:
	int mel;
	int nel;
	TArrayD zel, ael, wel;
	double den;
	double zp,ap;
	int ne;
	TArrayD etab, rtab, detab;	//tables: energy, range, deltaE
//	TArrayD vtab;
	double zw,aw;
	double a, b;		//linear interpolation coeficients, y = a*x + b
	int bsearch(int ntab, double *xtab, double x);
	double aitken(int ntab, double *xtab, double *ytab, double x);
	double aitken3(int ntab, double *xtab, double *ytab, double x);
	//====================================================================
	//== Interpolation by 4 points
	//====================================================================

	Double_t linear(int ntab, Double_t *xtab, Double_t *ytab, double x);
	Double_t linear(int ntab, Double_t *xtab, Double_t *ytab, Int_t i0, double x);

public:
	TELoss(void);
	TELoss(int _mel, double _den);
	virtual ~TELoss(void);
	void SetEL(int _mel, double _den);
	void SetDen(double _den) {den = _den;}
	int  AddEL(double _zel, double _ael, double _wel=1.);
	//===============================================================
	//== Add one element
	//===============================================================
	void SetZP(double _zp, double _ap);
	void SetEtab(int _ne, double e2);	//from 0 MeV to e2 MeV, _ne elements
	//===============================================================
	//== Set list of energies and calculate R(E)
	//===============================================================

	void SetDeltaEtab(double r);		//pro kazdou tlostku bude zvlastni tabulka, nebo taky ne
	//===============================================================
	//== Calculate and set list of dE(E)
	//===============================================================

	double GetZ() { return zp; };
	double GetA() { return ap; };
	double GetE(double e0, double r);
	//==================================================================
	//== Calculate new energy using linear interpolation
	//==================================================================

	double GetE0dE(double de, double r);	//de in MeV, r in microns
	//==================================================================
	//== Calculate new energy from deltaE using linear interpolation
	//==================================================================

	double GetE0dE(double de);				//de in MeV
	double GetE0(double e, double r);
	//==================================================================
	//== Calculate new energy using linear interpolation
	//==================================================================

	double GetEold(double e0, double r);
	//==================================================================
	//== Calculate new energy using aitken interpolation
	//==================================================================

	double GetE0old(double e, double r);
	//==================================================================
	//== Calculate new energy
	//==================================================================

	double GetR(double e0, double e);
	//==================================================================
	//== Calculate new energy using linear interpolation
	//==================================================================

	double GetRold(double e0, double e);
	//==================================================================
	//== Calculate new energy
	//==================================================================
	void PrintRE();
	void PrintdEE();
	void PrintREtoFile();
	void PrintdEEtoFile(const char* outfile = "outputdEE.txt");
	ClassDef(TELoss, 1)
};
