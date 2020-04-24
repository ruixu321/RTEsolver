#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>      // std::istringstream
#include <cmath> 
#include <algorithm>

#include "Grey.h"
#include "RadiativeFunctions.h"

using namespace std;

Grey::Grey(){

}


/**
 */

void Grey::ReadData(){


}

double Grey::a(){

}

int Grey::NbBands(){
    return NBands;
}

double Grey::a(int iBand, double T){

    return 0.;

}


/*!
 * Grey model, kappa is calculate based on http://www.sandia.gov/TNF/radiation.html
 */

double Grey::kappa(double T, double pH2O, double pCO2){

    vector<double> cH2O = {-0.23093,-1.12390,9.41530,-2.99880,0.51582,-1.8684E-5};
    vector<double> cCO2 = {18.741,-121.310,273.5,-194.050,56.310,-5.8169};
    double aH2O = cH2O[0] + cH2O[1]*(1000/T)+ cH2O[2]*pow(1000/T,2) + cH2O[3]*pow(1000/T,3) +cH2O[4]*pow(1000/T,4) + cH2O[5]*pow(1000/T,5);
    double aCO2 = cCO2[0] + cCO2[1]*(1000/T)+ cCO2[2]*pow(1000/T,2) + cCO2[3]*pow(1000/T,3) +cCO2[4]*pow(1000/T,4) + cCO2[5]*pow(1000/T,5);

    return (pCO2*aCO2 + pH2O*aH2O);
}

vector<vector<double>> Grey::solve(vector<double> mu, int NPoints, vector <double> T, vector <double> P, vector <double> XH2O, vector <double> XCO2, vector <double> DeltaX, string sign){

    vector<vector<vector<double>>> I(mu.size(), vector<vector<double>> (NbBands(), vector<double> (NPoints,0)));
    vector<vector<double>> ITot(mu.size(), vector<double> (NPoints,0)); 


    for(int muj=0; muj<mu.size(); muj++){

        for(int i=0; i<NPoints; i++){

            int iBand = 0.;

            double pH2O = XH2O[i]*P[i];
            double pCO2 = XCO2[i]*P[i];
	

            if(i==0){

                I[muj][iBand][i] = BlackBodyIntensity(T[i]); // first approximation cold wall
            }
            else{

                I[muj][iBand][i] = BlackBodyIntensity(T[i]) + (I[muj][iBand][i-1] - BlackBodyIntensity(T[i]))*exp(-kappa(T[i],pH2O,pCO2)*DeltaX[i-1]/mu[muj]); 

            }

            ITot[muj][i] = I[muj][iBand][i]; 

        }

        if(sign == "-"){
            reverse(ITot[muj].begin(), ITot[muj].end());
        }
    }

    return ITot; 

}
