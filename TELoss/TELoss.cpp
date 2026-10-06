#include "TELoss.h"

ClassImp(TELoss)

extern "C" void eloss_(int *nel, double *zel, double *ael, double *wel, double *den
		,double *zp, double *ap, int *ne, double *etab, double *rtab
		,double *zw, double *aw);


//extern "C" __declspec( dllimport ) void _stdcall
//	ELOSS(int *nel, double *zel, double *ael, double *wel, double *den
//		,double *zp, double *ap, int *ne, double *etab, double *rtab
//		,double *zw, double *aw);

TELoss::TELoss(void)
{
	ne = 0;
	mel = 0;
} //-----------------------------------------------------------------

TELoss::TELoss(int _mel, double _den)
{
	ne = 0;
	mel = 0;
	SetEL(_mel,_den);
} //-----------------------------------------------------------------

TELoss::~TELoss(void)
{
	//if(RE)   delete RE;
	//if(ER)   delete ER;
} //-----------------------------------------------------------------

void TELoss::SetEL(int _mel, double _den)
{
	den = _den;
	if(_mel != mel)
	{
		mel = _mel;
		zel.Set(mel);
		ael.Set(mel);
		wel.Set(mel);
	}
	nel = 0;
} //-----------------------------------------------------------------

int TELoss::AddEL(double _zel, double _ael, double _wel)
{
	//===============================================================
	//== Add one element
	//===============================================================

	if(nel >= mel) return -1;

	zel[nel] = _zel;
	ael[nel] = _ael;
	wel[nel] = _wel;
	return ++nel;
} //-----------------------------------------------------------------

void TELoss::SetZP(double _zp, double _ap)
{
	zp = _zp;
	ap = _ap;
} //-----------------------------------------------------------------

void TELoss::SetEtab(int _ne, double e2)
{
	//===============================================================
	//== Set list of energies and calculate R(E)
	//===============================================================

	int i;
	double e1 = 0.;
	double es;
	double coef = 10. / den; //== convertion mg/cm^2 -> micron

	if(_ne != ne)
	{
		ne = _ne;		//asi pocet elementu
		etab.Set(ne);
		rtab.Set(ne);
	}

	if ( ((e2-e1)/(static_cast<double>(ne))) > 0.01) { Warning("SetEtab", "Table resolution can be too low for some cases."); }

	es = (e2 - e1) / (double)(ne-1);	//energeticky skok v tabulce
	for(i = 0; i < ne; ++i)
	{
		etab[i] = e1;
		e1 += es;
	}

	eloss_(&nel, zel.GetArray(), ael.GetArray(), wel.GetArray(), &den
		,&zp, &ap
		,&ne, etab.GetArray(), rtab.GetArray(), &zw, &aw);

	for(i = 0; i < ne; ++i) rtab[i] = coef * rtab[i];

	//if(RE) delete RE;
	//if(ER) delete ER;

	//RE = new TSpline3("RE",etab.GetArray(),rtab.GetArray(),ne);
	//ER = new TSpline3("ER",rtab.GetArray(),etab.GetArray(),ne);
} //-----------------------------------------------------------------

void TELoss::SetDeltaEtab(double r)
{
	//===============================================================
	//== Calculate and set list of dE(E)
	//===============================================================

	int i;
//	double es;
	//double coef = 10. / den; //== convertion mg/cm^2 -> micron

	if(ne == 0) {
		Warning("SetDeltaEtab", "Etab was not initialized yet");
		return;
	}

	detab.Set(ne);

//	if ( ((e2-e1)/(static_cast<double>(ne))) > 0.01) { Warning("SetEtab", "Table resolution can be too low for some cases."); }

//	es = (e2 - e1) / (double)(ne-1);	//energeticky skok v tabulce
	for(i = 0; i < ne; ++i)
	{
		if ( GetE(etab[i], r) != 0. ) {
			detab[i] = etab[i] - GetE(etab[i], r);
//			e1 += es;
		}
		else { detab[i] = 0.; }
	}

//	eloss_(&nel, zel.GetArray(), ael.GetArray(), wel.GetArray(), &den
//		,&zp, &ap
//		,&ne, etab.GetArray(), rtab.GetArray(), &zw, &aw);
//
//	for(i = 0; i < ne; ++i) rtab[i] = coef * rtab[i];

	return;
}

double TELoss::GetR(double e0,double e)
{
	//==================================================================
	//== Calculate new energy using linear interpolation
	//==================================================================

	//double r0 = RE->Eval(e0);
	//return ER->Eval(r0 - r);

	double r0 = linear(ne, etab.GetArray(), rtab.GetArray(), e0);
	double r = linear(ne, etab.GetArray(), rtab.GetArray(), e);
//	return aitken3(ne, rtab.GetArray(), etab.GetArray(), r0-r);
	return r0 - r;

} //--------------------------------------------------------------------


double TELoss::GetRold(double e0,double e)
{
	//==================================================================
	//== Calculate new energy
	//==================================================================

	//double r0 = RE->Eval(e0);
	//return ER->Eval(r0 - r);

	double r0 = aitken3(ne, etab.GetArray(), rtab.GetArray(), e0);
	double r = aitken3(ne, etab.GetArray(), rtab.GetArray(), e);
//	return aitken3(ne, rtab.GetArray(), etab.GetArray(), r0-r);
	return r0 - r;

} //--------------------------------------------------------------------


double TELoss::GetE(double e0,double r)
{
	//==================================================================
	//== Calculate new energy using linear interpolation
	//==================================================================

//	Double_t r0 = ylinear(ne, e0);
	Double_t r0 = linear(ne, etab.GetArray(), rtab.GetArray(), e0);
	if (r0 < r) { return 0; }

	return linear(ne, rtab.GetArray(), etab.GetArray(), r0-r);
} //--------------------------------------------------------------------

double TELoss::GetE0dE(double de, double r)
{
	//==================================================================
	//== Calculate new energy from deltaE using linear interpolation
	//==================================================================

	SetDeltaEtab(r);	//problem s interpolaci

	if ( de > TMath::MaxElement(ne, detab.GetArray()) ) { return -1; }

	int i0 = 0;
	while ( detab[i0] == 0 ) {
		i0++;
	}

	//vypis
//	cout << "i0 = \t" << i0 << "\tdetab[i0] = \t" << detab[i0] << "\tetab[i0] = \t" << etab[i0] << endl;

	int i1 = ne-1;
	int ix;
	while(i0 < i1 - 1)
	{
		ix = (i0 + i1) / 2;
		if ( de > detab[ix] ) i1 = ix;
		else             i0 = ix;
	}


//	if (r0 < r) { return 0; }
//	cout << "r0 =\t\t" << r0 << endl;
//	cout << "r0-r =\t\t" << r0-r << endl;
//	cout << "E(r0-r) =\t\t" << xlinear(ne, r0-r) << endl;
//	cout << "E(r0-r) =\t\t" << linear(ne, rtab.GetArray(), etab.GetArray(), r0-r) << endl;
	return linear(ne, detab.GetArray(), etab.GetArray(), i0, de);
} //--------------------------------------------------------------------

double TELoss::GetE0dE(double de)
{
	if ( de > TMath::MaxElement(ne, detab.GetArray()) ) { return -1; }

	int i0 = 0;
	while ( detab[i0] == 0 ) {
		i0++;
	}

//	cout << "detab[i0] je " << detab[i0] << endl;

	int i1 = ne-1;
	int ix;
	while(i0 < i1 - 1)
	{
		ix = (i0 + i1) / 2;
		if ( de > detab[ix] ) i1 = ix;
		else             i0 = ix;
	}

//	cout << "detab[i0] je " << detab[i0] << endl
//		<< "detab[i1] je " << detab[i1] << endl;

//	if (r0 < r) { return 0; }
//	cout << "r0 =\t\t" << r0 << endl;
//	cout << "r0-r =\t\t" << r0-r << endl;
//	cout << "E(r0-r) =\t\t" << xlinear(ne, r0-r) << endl;
//	cout << "E(r0-r) =\t\t" << linear(ne, rtab.GetArray(), etab.GetArray(), r0-r) << endl;
	return linear(ne, detab.GetArray(), etab.GetArray(), i0, de);

	//return linear(ne, detab.GetArray(), etab.GetArray(), de);
}

double TELoss::GetEold(double e0,double r)
{
	//==================================================================
	//== Calculate new energy using aitken interpolation
	//==================================================================

	//double r0 = RE->Eval(e0);
	//return ER->Eval(r0 - r);

	double r0 = aitken3(ne, etab.GetArray(), rtab.GetArray(), e0);
	return aitken3(ne, rtab.GetArray(), etab.GetArray(), r0-r);
} //--------------------------------------------------------------------

double TELoss::GetE0(double e0,double r)
{
	//==================================================================
	//== Calculate new energy using linear interpolation
	//==================================================================

	Double_t r0 = linear(ne, etab.GetArray(), rtab.GetArray(), e0);
	return linear(ne, rtab.GetArray(), etab.GetArray(), r0+r);
} //--------------------------------------------------------------------

double TELoss::GetE0old(double e,double r)
{
	//==================================================================
	//== Calculate new energy
	//==================================================================

	//double r0 = RE->Eval(e);
	//return ER->Eval(r0 + r);

	double r0 = aitken3(ne, etab.GetArray(), rtab.GetArray(), e);
	return aitken3(ne,rtab.GetArray(),etab.GetArray(),r0+r);
} //-----------------------------------------------------------------

void TELoss::PrintRE()
{
	int i;
	for(i = 0; i < ne; ++i)
	{
		printf("%10.2lf %10.3lf\n",etab[i],rtab[i]);
	}
} //-----------------------------------------------------------------

void TELoss::PrintdEE()
{
	int i;
	for(i = 0; i < ne; ++i)
	{
		printf("%10.3lf %10.3lf\n",etab[i],detab[i]);
	}
} //-----------------------------------------------------------------

void TELoss::PrintREtoFile()
{
	FILE *fw;
	fw = fopen("outputRE.txt", "w");

	int i;
	for(i = 0; i < ne; ++i)
	{
		fprintf(fw, "%10.2lf \t%10.3lf\n",etab[i],rtab[i]);
	}

	fclose(fw);
}

void TELoss::PrintdEEtoFile(const char* outfile)
{
	FILE *fw;
	fw = fopen(outfile, "w");

	int i;
	for(i = 0; i < ne; ++i)
	{
		fprintf(fw, "%10.3lf \t%10.3lf\n",etab[i],detab[i]);
	}

	fclose(fw);
}

double TELoss::aitken3(int ntab, double *xtab, double *ytab, double x)
{
	//====================================================================
	//== Interpolation by 4 points
	//====================================================================

//	int ndeg = 4;
	double y[4];

	int i0,i;

	i0 = bsearch(ntab,xtab,x) - 1;
	if(i0 < 0) i0 = 0;
	else if(i0 + 3 >= ntab) i0 = ntab - 4;
	for(i = 0; i < 4; ++i) y[i] = ytab[i0 +i];


	return aitken(3,xtab+i0,y,x);
} //-----------------------------------------------------------------

int TELoss::bsearch(int ntab, double *xtab, double x)
{
	int i0 = 0;
	int i1 = ntab;
	int ix;
	while(i0 < i1 - 1)
	{
		ix = (i0 + i1) / 2;
		if(x < xtab[ix]) i1 = ix;
		else             i0 = ix;
	}
	return i0;
} //-----------------------------------------------------------------

double TELoss::aitken(int ntab, double *xtab, double *ytab, double x)
{
	int i,j;
	for(j = 0; j < ntab; ++j)
	{
		for(i = j+1; i <= ntab; ++i)
		{
			ytab[i] = ((x - xtab[j])*ytab[i] - (x - xtab[i])*ytab[j])
				/(xtab[i] - xtab[j]);
		}
	}
	return ytab[ntab];
} //-----------------------------------------------------------------

Double_t TELoss::linear(int ntab, Double_t *xtab, Double_t *ytab, double x)
{
	Double_t y;
	Double_t	a = 0;
	Double_t	b = 0;

	int i0 = 0;

	i0 = bsearch(ntab, xtab, x) - 1;
	if (i0 < 0) { i0 = 0; }
	else {
		if ( i0+2 >= ntab ) { i0 = i0 - 1; }
	}

//	cout << "xtab[i0] = \t" << xtab[i0] << "\txtab[i0+2] = \t" << xtab[i0+2] << endl;
//	cout << "ytab[i0] = \t" << ytab[i0] << "\tytab[i0+2] = \t" << ytab[i0+2] << endl;

	a = ( ytab[i0] - ytab[i0+2] )/( xtab[i0] - xtab[i0+2] );
	b = ( ytab[i0+2]*xtab[i0] - ytab[i0]*xtab[i0+2] )/( xtab[i0] - xtab[i0+2] );

	Double_t lerror = TMath::Abs(( (1/a)*ytab[i0+1] - (b/a) ) - xtab[i0+1]);
	if (lerror > 0.01) {
		Warning("linear", "Linear interpolation does not give sufficient precision: %f", lerror);
		cout << "xtab[i0+1] \t" << xtab[i0+1] << "\tytab[i0+1] \t" << ytab[i0+1] << endl;
	}

	y = a*x + b;

	return y;
}

Double_t TELoss::linear(int ntab, Double_t *xtab, Double_t *ytab, Int_t i0, double x)
{
	Double_t y;
	Double_t	a = 0;
	Double_t	b = 0;

	//int i0 = 0;

	//i0 = bsearch(ntab, xtab, x) - 1;
	if (i0 < 0) { i0 = 0; }
	else {
		if ( i0+2 >= ntab ) { i0 = i0 - 1; }
	}

//	cout << "xtab[i0] = \t" << xtab[i0] << "\txtab[i0+2] = \t" << xtab[i0+2] << endl;
//	cout << "ytab[i0] = \t" << ytab[i0] << "\tytab[i0+2] = \t" << ytab[i0+2] << endl;

	a = ( ytab[i0] - ytab[i0+2] )/( xtab[i0] - xtab[i0+2] );
	b = ( ytab[i0+2]*xtab[i0] - ytab[i0]*xtab[i0+2] )/( xtab[i0] - xtab[i0+2] );

	Double_t lerror = TMath::Abs(( (1/a)*ytab[i0+1] - (b/a) ) - xtab[i0+1]);
	if (lerror > 0.001) { Warning("linear_ext", "Linear interpolation does not give sufficient precision: %f", lerror); }

//	cout << "a je " << a << "\tb je " << b << endl;

	y = a*x + b;

	return y;
}
