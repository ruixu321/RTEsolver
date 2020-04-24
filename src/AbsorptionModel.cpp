#include <string>
#include <iostream>

#include "AbsorptionModel.h"
#include "WSGG.h"
#include "WSGGBordbar.h"
#include "WSGGJohansson.h"
#include "SNB.h"
#include "Grey.h"

using namespace std;

AbsorptionModel *AbsorptionModel::make_absorptionModel(string name){

    if(name == "WSGG"){
        cout << " Creation of the WSGG model " << endl;
        return new WSGG;
    }
    else if(name == "WSGGBordbar"){
        cout << " Creation of the WSGG model with Bordbar model" << endl;
        return new WSGGBordbar;
    }
    else if(name == "WSGGJohansson"){
        cout << " Creation of the WSGG model with Johansson model" << endl;
        return new WSGGJohansson;
    }
    else if(name == "SNB"){
        cout << " Creation of the SNB model " << endl;
        return new SNB;
    }
    else if(name == "Grey"){
        cout << " Creation of the Grey model " << endl;
        return new Grey;
    }
    else{
        cout << "error: no " << name << " Model" << endl;
    }
}

