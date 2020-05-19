#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>      // std::istringstream
#include <cmath> 
#include <algorithm>

#include "WSGGJohansson.h"
#include "RadiativeFunctions.h"


using namespace std;

WSGGJohansson::WSGGJohansson(){

	kappa1 = {0.055, 0.88, 10.0, 135.0};
	kappa2 = {0.012, -0.021, -1.6, -35.0};

	C1j1 = {0.358, 0.392, 0.142, 0.0798};
	C1j2 = {0.0731, -0.212, -0.0831, -0.0370};
	C1j3 = {-0.0466, 0.0191, 0.0148, 0.0023};
	C2j1 = {-0.165, -0.291, 0.348, 0.0866};
	C2j2 = {-0.0554, 0.644, -0.294, -0.106};
	C2j3 = {0.0930,-0.209, 0.0662, 0.0305};
	C3j1 = {0.0598, 0.0784, -0.122, -0.0127};
	C3j2 = {0.0028,-0.197, 0.118, 0.0169};
	C3j3 = {-0.0256, 0.0662 , -0.0295, -0.0051};

	aOutput = vector<vector<double>>(NBands , vector <double> (0, 0.));  // weights for bands
	kappaOutput = vector<vector<double>>(NBands , vector <double> (0, 0.)); // absorption coefficients for bands

}


/**
 * Read the coefficients for wsgg superposition model of Cassol (2013)
 */

void WSGGJohansson::ReadData(){

}

double WSGGJohansson::a(){

	return 0.;

}

int WSGGJohansson::NbBands(){
	return NBands;
}

// vector<vector<double>> WSGGJohansson::getAOutput(){
// 	return aOutput;
// }

// vector<vector<double>> WSGGJohansson::getKappaOutput(){
// 	return kappaOutput;
// }


double WSGGJohansson::a(int iBand, double T){

	return 0.;
}

double WSGGJohansson::a(int igas, double T, double MR){

	return b(igas,0,MR)*pow(T,0)*1e-0 + b(igas,1,MR)*pow(T,1) + b(igas,2,MR)*pow(T,2);
}

double WSGGJohansson::b(int igas,  int i, double MR){
	switch (i) {
		case 0: //cout << "case 0 " << endl;
			return C1j1[igas] + C2j1[igas]*MR + C3j1[igas]*pow(MR,2);
		case 1: //cout << "case 1 " << endl;
			return C1j2[igas] + C2j2[igas]*MR + C3j2[igas]*pow(MR,2);
		case 2: //cout << "case 2 " << endl;
			return C1j3[igas] + C2j3[igas]*MR + C3j3[igas]*pow(MR,2);	
		default:
			cout << "Error in calculating b!!!! " << endl;
			return 0.;
	}
}

double WSGGJohansson::kappa(int igas, double MolarFracH2O, double MolarFracCO2) {  // Modified by Rui

	double MR;

	if (MolarFracCO2 == 0) {
		return 0.0 ;
	}

	MR = MolarFracH2O / MolarFracCO2;
	MR = min(MR, 2.0);
    MR = max(MR, 0.125);

	return kappa1[igas] + kappa2[igas] * MR;

}

vector<vector<double>> WSGGJohansson::solve(vector<double> mu, int NPoints, vector <double> T, vector <double> P, 
	vector <double> XH2O, vector <double> XCO2, vector <double> DeltaX, string sign)
{

	vector<vector<vector<double>>> I(mu.size(), vector<vector<double>> (NbBands(), vector<double> (NPoints,0)));
	vector<vector<double>> ITot(mu.size(), vector<double> (NPoints,0)); 
	double pH2O, pCO2, MR;

	if (sign == "+") {
		kappaOutput = vector<vector<double>> (NBands , vector <double> (XCO2.size(), 0.));
		aOutput = vector<vector<double>> (NBands , vector <double> (XCO2.size(), 0.));
	}

	for(int muj=0; muj<mu.size(); muj++){

		for(int i=0; i<NPoints; i++){

			// compute a coefficients
			vector <double> acoeffs (5, 0.); 
			vector <double> kappa_coeffs (5, 0.); 
			double sum_of_coeffs = 0;

			for(int iBand = 1; iBand < NbBands(); iBand++){
				if (i == 0) {
					MR = 2.0;
				}
				else {
					if (XCO2[i] != 0.0) {
						MR = XH2O[i] / XCO2[i];
					}
					MR = min(MR, 2.0);
					MR = max(MR, 0.125);
					kappa_coeffs[iBand] = kappa(iBand-1, XH2O[i-1], XCO2[i-1]);
				}
				acoeffs[iBand] = a(iBand-1, T[i]/1200., MR);
				sum_of_coeffs += acoeffs[iBand];
			}

			acoeffs[0] = 1 - sum_of_coeffs;
			kappa_coeffs[0] = 0; //clear gas

			for(int iBand = 0; iBand<NbBands(); iBand++){

				if(i==0){

					I[muj][iBand][i] = acoeffs[iBand]*BlackBodyIntensity(T[i]);
				}
				else{

					I[muj][iBand][i] = acoeffs[iBand]*BlackBodyIntensity(T[i]) + 
					(I[muj][iBand][i-1] - acoeffs[iBand]*BlackBodyIntensity(T[i])) * 
					exp(-kappa_coeffs[iBand] * (XH2O[i-1]+XCO2[i-1]) * 0.986923 * (P[i-1]) * DeltaX[i-1]/mu[muj]);

					if (sign == "+" && i != NPoints - 1) {
						kappaOutput[iBand][i] = kappa_coeffs[iBand];
						aOutput[iBand][i] = acoeffs[iBand];
					}
				}
				ITot[muj][i] = (iBand==0) ? I[muj][iBand][i] : ITot[muj][i] + I[muj][iBand][i]; 
			}
		}

		if(sign == "-"){

			reverse(ITot[muj].begin(), ITot[muj].end());
		}
	}

	return ITot; 

}
