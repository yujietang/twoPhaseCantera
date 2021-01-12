#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Liquid.h"
#include "AllSpecies.h"

using namespace std;

// NDodecane ndc(101325);
// ICetane ict(101325);
// TransDecalin tdc(101325);
// Toluene tln(101325);

Ethanol eth(101325);

const int tl = 250;
const int th = 500;

int main()
{
    // 250 K ~ 500 K
    vector<double> tem;
    for (int i = tl; i <= th; i++) {
        tem.push_back(static_cast<double>(i));
    }
    // ofstream fndc("NC12H26.csv");
    // ofstream fict("IC16H36.csv");
    // ofstream ftdc("C10H18.csv");
    // ofstream ftln("C7H8.csv");
    ofstream feth("C2H5OH.csv");

    for (double t : tem) {
        // fndc << t << ',' << ndc.rho(t) << ',' << ndc.pv(t) << ',' << ndc.Lv(t) << endl;
        // fict << t << ',' << ict.rho(t) << ',' << ict.pv(t) << ',' << ict.Lv(t) << endl;
        // ftdc << t << ',' << tdc.rho(t) << ',' << tdc.pv(t) << ',' << tdc.Lv(t) << endl;
        // ftln << t << ',' << tln.rho(t) << ',' << tln.pv(t) << ',' << tln.Lv(t) << endl;
        feth << t << ',' << eth.rho(t) << ',' << eth.pv(t) << ',' << eth.Lv(t) 
        << ',' << eth.D(t) << ',' << eth.cpg(t)<< endl;
    }

    return 0;
}

