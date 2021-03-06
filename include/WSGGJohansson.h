#ifndef WSGGJohansson_H
#define WSGGJohansson_H

#include <string>
#include <iostream>
#include <vector>

#include "AbsorptionModel.h"

using namespace std;

class WSGGJohansson: public AbsorptionModel
{
    private:

    int NBands = 5;
        
	vector <double> kappa1;
	vector <double> kappa2;

	vector <double>  C1j1;
	vector <double>  C1j2;
	vector <double>  C1j3;
	vector <double>  C2j1;
	vector <double>  C2j2;
	vector <double>  C2j3;
	vector <double>  C3j1;
	vector <double>  C3j2;
	vector <double>  C3j3;

	vector<vector<double>> aOutput;
	vector<vector<double>> kappaOutput;

    public:

	// compute the coefficients for WSGGJohansson model

	double ComputeCoeffs() {return 0;};
	double ComputeCoeffs( int n, double T);

	// read WSGGJohansson parameters
	void ReadData();

	// return coeffs

	double transmissivity(const vector <double>& T, const vector <double>& P, const vector <double>& XH2O, const vector <double>& XCO2, const vector <double>& delta, int iBand, double mu, int i, int j){ return 0;};


	double kappa(int iBand, double pH2O, double pCO2);
	double kappa(int iBand, double T){ return 0;};
	double kappa(double T, double pH2O, double pCO2){return 0;};

	double a();
	double a(int n, double T); // return a coeff for band n at temperature T
	double a(int n, double T, double MR); // return a coeff for band n at temperature T


	double b(int i, int j, double MR);

	vector<vector<double>> solve(vector<double> mu, int Npoints, vector <double> T, vector <double> P, vector <double> XH2O, vector <double> XCO2, vector <double> DeltaX, string sign);


	// wavenumber
	double nu(int iBand) {return 0.;};

	// return number of bands

	int NbBands();

	vector<vector<double>> getAOutput(){return aOutput;};
	vector<vector<double>> getKappaOutput(){return kappaOutput;};

	// write data

	WSGGJohansson();

};

#endif
