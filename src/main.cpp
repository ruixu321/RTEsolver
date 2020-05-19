#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>      // std::istringstream

#include "RadiationModel.h"

void readData(string input_file, vector<double>& x, vector<double>& XH2O, vector<double>& XCO2, 
                vector <double>& T, vector <double>& P) {

    cout << " start reading " << endl;

    std::ifstream infile(input_file);

    vector <double> xinput;
    vector <double> XH2Oinput;
    vector <double> XCO2input;
    vector <double> Tinput;
    vector <double> Pinput;

    int Nelems = 0;
    if (infile) {
        string line;

        int nline = 0;

        while (getline(infile, line)) {

            istringstream iss(line);

            if(nline!=0){
                iss >> xinput[nline-1];
                iss >> XCO2input[nline-1]; 
                iss >> XH2Oinput[nline-1];
                iss >> Tinput[nline-1];
                iss >> Pinput[nline-1];
                xinput[nline-1] = xinput[nline-1];

                cout << nline-1 << " " << xinput[nline-1] << " " << Tinput[nline-1] 
                    << " " << XCO2input[nline-1] << " " << XH2Oinput[nline-1] << Pinput[nline-1] << endl;
            }
            else{

                iss >> Nelems;

                xinput = vector<double> (Nelems, 0.);
                XH2Oinput = vector<double> (Nelems, 0.);
                XCO2input = vector<double> (Nelems, 0.);
                Tinput = vector<double> (Nelems, 0.);
                Pinput = vector<double> (Nelems, 0.);


		//Tinput[0] = 300.;
		//Tinput[Tinput.size()-1] = 300.;
		//xinput[0] = 0.;
		//xinput[Tinput.size()-1] = xinput[xinput.size()-2] + 0.01;

		//XCO2input[0] = 0.;
		//XCO2input[Nelems+1] = 0.;
		// XH2Oinput[0] = 0.;
		// XH2Oinput[Nelems+1] = 0.;

	    }

	    nline++;
	}
    }

    // interpolation

    auto interpolation = [](double x, double xa, double xb, double ya, double yb){
	    return ya + (x - xa)*(yb-ya)/(xb-xa);
    };

    x = xinput;

    XCO2 = vector <double> (XCO2input.size()-1,0.);

    XH2O = vector <double> (XH2Oinput.size()-1,0.);

    P = vector <double> (XH2Oinput.size()-1,0.);

    T = vector <double> (Tinput.size()+1, 0.);

    T[0] = Tinput[0];

    T[T.size()-1] = Tinput[Tinput.size()-1];

    for(int i = 0; i<XCO2.size(); ++i){

	    double xinterp = (x[i+1] - x[i])/2. + x[i];
	    T[i+1] = interpolation(xinterp, x[i], x[i+1], Tinput[i], Tinput[i+1]);
	    XCO2[i] = interpolation(xinterp, x[i], x[i+1], XCO2input[i], XCO2input[i+1]);
	    XH2O[i] = interpolation(xinterp, x[i], x[i+1], XH2Oinput[i], XH2Oinput[i+1]);
        P[i] = interpolation(xinterp, x[i], x[i+1], Pinput[i], Pinput[i+1]);

	    cout << x[i] << " " << x[i+1] << " " <<  xinterp << " " << T[i+1] << " " << XCO2[i] 
            << " " << XH2O[i] << P[i] << endl;

    }
 
    // constant pressure

    // P = vector <double> (XCO2.size(), 0.);
    // generate(P.begin(), P.end(), [n = 0.] () mutable { return n=1; });


}

void readData_no_interpolation(string input_file, vector<double>& x, vector<double>& XH2O, 
                                vector<double>& XCO2, vector <double>& T, vector <double>& P) {

    cout << " start reading " << endl;

    std::ifstream infile(input_file);

    vector <double> xinput;
    vector <double> XH2Oinput;
    vector <double> XCO2input;
    vector <double> Tinput;
    vector <double> Pinput;

    int Nelems = 0;
    if (infile) {
        string line;

        int nline = 0;

        while (getline(infile, line)) {

            istringstream iss(line);

            if(nline!=0){
                iss >> xinput[nline-1];
                iss >> XCO2input[nline-1]; 
                iss >> XH2Oinput[nline-1];
                iss >> Tinput[nline-1];
                iss >> Pinput[nline-1];
                xinput[nline-1] = xinput[nline-1];

                cout << nline-1 << " " << xinput[nline-1] << " " << Tinput[nline-1] << " " << 
                XCO2input[nline-1] << " " << XH2Oinput[nline-1] << Pinput[nline-1] << endl;
            }
            else{


                iss >> Nelems;

                xinput = vector<double> (Nelems, 0.);
                XH2Oinput = vector<double> (Nelems, 0.);
                XCO2input = vector<double> (Nelems, 0.);
                Tinput = vector<double> (Nelems, 0.);
                Pinput = vector<double> (Nelems, 0.);

		//Tinput[0] = 300.;
		//Tinput[Tinput.size()-1] = 300.;
		//xinput[0] = 0.;
		//xinput[Tinput.size()-1] = xinput[xinput.size()-2] + 0.01;


		//XCO2input[0] = 0.;
		//XCO2input[Nelems+1] = 0.;
		//XH2Oinput[0] = 0.;
		//XH2Oinput[Nelems+1] = 0.;


	    }

	    nline++;
	}
    }

    x = xinput;
    XCO2 = vector <double> (XCO2input.size(),0.);
    XH2O = vector <double> (XH2Oinput.size(),0.);
    P = vector <double> (XH2Oinput.size(),0.);
    T = vector <double> (Tinput.size()+2, 0.);


    T[0] = 300.;
    T[T.size()-1] = 300.;

    for(int i = 0; i<XCO2.size(); ++i){

	    T[i+1] = Tinput[i];
	    XCO2[i] = XCO2input[i];
	    XH2O[i] = XH2Oinput[i];
        P[i] = Pinput[i];

	    cout << x[i] << " " << x[i+1] << " "  << T[i+1] << " " << XCO2[i] << " " << XH2O[i] << P[i] << endl;

    }
 
    // constant pressure

    // P = vector <double> (XCO2.size(), 0.);
    // generate(P.begin(), P.end(), [n = 0.] () mutable { return n=1; });


}


int main(int argc, char** argv){

	// init

	vector <double> x; // no need to modify 
	vector <double> XCO2;
	vector <double> XH2O;
	vector <double> T;
	vector <double> P;

	// input file

	string input_file = string(argv[1]);

	// save

	string directory_save = string(argv[2]);

	// interpolation or not

	string interpolation_bool = string(argv[3]);


	if (interpolation_bool == "True"){

		readData(input_file, x, XH2O, XCO2, T, P);
	}
	else{

		readData_no_interpolation(input_file, x, XH2O, XCO2, T, P);
	}

	// start

	cout << " Start " << endl;

	// creation of the radiation models

	RadiationModel WSGGJ("WSGGJohansson", x, T, XH2O, XCO2, P, 
        "./res/" + directory_save + "/dqr_WSGGJohansson.res", 
        "./res/" + directory_save + "/qr_WSGGJohansson.res", 
        "./res/" + directory_save + "/a_WSGGJohansson.res",
        "./res/" + directory_save + "/kappa_WSGGJohansson.res");

    RadiationModel SNB("SNB", x, T, XH2O, XCO2, P, 
        "./res/" + directory_save + "/dqr_SNB.res", 
        "./res/" + directory_save + "/qr_SNB.res",
        "./res/" + directory_save + "/a_SNB.res",
        "./res/" + directory_save + "/kappa_SNB.res");

	// RadiationModel WSGGB("WSGGBordbar", x, T, XH2O, XCO2, P,
 //        "./res/" + directory_save + "/dqr_WSGGBordbar.res", 
 //        "./res/" + directory_save + "/qr_WSGGBordbar.res",
 //        "./res/" + directory_save + "/a_WSGGBordbar.res",
 //        "./res/" + directory_save + "/kappa_WSGGBordbar.res");

	// RadiationModel WSGGC("WSGG", x, T, XH2O, XCO2, P, "./res/" + directory_save + "/dqr_WSGG.res", "./res/" + directory_save + "/qr_WSGG.res");
	
	// RadiationModel Grey("Grey", x, T, XH2O, XCO2, P, "./res/" + directory_save + "/dqr_Grey.res", "./res/" + directory_save + "/qr_Grey.res");

	// solve RTE

	WSGGJ.Solve();
    SNB.Solve();
	// WSGGB.Solve();
	// WSGGC.Solve();
	
	// Grey.Solve();

	return 0;
}
