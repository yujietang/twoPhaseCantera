#include "Inlet1D.h"
#include "cantera/IdealGasMix.h"
#include "cantera/transport.h"
#include "Sim1D.h"
#include "SprayStFlow.h"
#include "Lagrangian.h"

#include <fstream>

using namespace Cantera;
using fmt::print;

int SprayFreeFlame(doublereal flamespeed, bool do_spray)
{
    /************************************************************/
    //                       Injection
    /************************************************************/
    doublereal parcelDiameter(25e-6); // injection droplet diameter [m]
    doublereal injTemperature(300); // injection parcel's temperature [K]
    doublereal injPressure(1.0*OneAtm); // injection pressure [Pa]
    doublereal dtlag = 10e-6; // Parcel Tracking Time Scale(Evaporation Time Scale)
    doublereal Mdot_liquid;
    doublereal dropletNumber;
    /************************************************************/
    //                          Mesh
    /************************************************************/
    const size_t meshPointNumber = 800; // Mesh Point Number
    const doublereal domainLength = 0.1; // Domain Length
    doublereal minGrid = parcelDiameter;// min grid size.
    bool refine_grid = false; // Refined mesh has been turned off
    /************************************************************/
    //                        Gas Flow
    /************************************************************/
    const doublereal temp = 300.0; // inlet gas flow temperature [K]
    const doublereal pressure = 1.0*OneAtm; // inlet gas flow pressure [atm]
    const doublereal phi_over = 1.0;// overall equivalence ratio
    const doublereal phi = 1.0; // gas flow equivalence ratio
    doublereal uin = 0.3; // initial guess of inlet velocity [m/s]
    doublereal Mdot_air;  // initial mass flux of air flow
    doublereal Mdot_fuel; // initial mass flux of fuel flow
    /************************************************************/
    //                        Chemistry
    /************************************************************/
    IdealGasMix gas("Ethanol_31.cti", "gas");//for single component of ethanol fuel:
    doublereal C_atoms = 2.0; // C2H5OH
    doublereal H_atoms = 6.0; // C2H5OH
    doublereal O_atoms = 1.0; // C2H5OH
    size_t nsp = gas.nSpecies(); // number of species
    vector_fp x(nsp, 0.0);
    doublereal ax = C_atoms + H_atoms/4.0 - O_atoms/2.0; // air consumption
    // doublereal fa_stoic = 1.0 / (4.76 * ax); // fuel / air ratio at stoic state
    doublereal fa_stoic = 0.11232;  
    // Species' mole fraction when the condition is fuel/oxidizer flow:
    x[gas.speciesIndex("C2H5OH")] = 1.0; 
    x[gas.speciesIndex("O2")] = 0.21 / phi / fa_stoic;
    x[gas.speciesIndex("N2")] = 0.78 / phi / fa_stoic;
    x[gas.speciesIndex("AR")] = 0.01 / phi / fa_stoic;
    // when the condition is pure air flow  (1 mole air):
    // x[gas.speciesIndex("O2")] = 0.21;
    // x[gas.speciesIndex("N2")] = 0.78;
    // x[gas.speciesIndex("AR")] = 0.01;
    /************************************************************/
    /************************************************************/

    gas.setState_TPX(temp, pressure, x.data());
    // gas.setState_TP(temp, pressure);
    doublereal rho_in = gas.density();
    const doublereal rho_in_fuel = 0.085;
    std::cout << "inlet rhog = " << rho_in << std::endl;
    vector_fp yin(nsp);
    gas.getMassFractions(&yin[0]);
    /**********Evaluate the liquid mass flux**********/
    doublereal uf = flamespeed; //laminar flame speed;
    Mdot_fuel = rho_in_fuel*uf;
    Mdot_air = rho_in*uf - Mdot_fuel;
    Mdot_liquid = fa_stoic*phi_over*Mdot_air - Mdot_fuel; 
    /*************************************************/
    gas.equilibrate("HP");//evaluate the adiabatic flame temperature.
    vector_fp yout(nsp);
    gas.getMassFractions(&yout[0]);
    doublereal rho_out = gas.density();
    doublereal Tad = gas.temperature();
    print("phi = {}, Tad = {}\n", phi, Tad);

    /**********Lagrangian Cloud**********/
    Lagrangian cloud(
        parcelDiameter,
        injTemperature,
        injPressure,
        dtlag,
        Mdot_liquid,
        dropletNumber
    );

    /**********Eulerian Flow**********/
    StFlow gasflow(&gas);
    gasflow.setFreeFlow();
    gasflow.setViscosityFlag(true);
    gasflow.setupCloud(cloud);

    /*************************Create Grid*************************/
    // create an initial grid
    int nz = meshPointNumber;//initial grid point number
    doublereal lz = domainLength;//initial grid length
    vector_fp z(nz);//initial grid point vector
    doublereal dz = lz/((doublereal)(nz - 1));//initial grid size
    //initialize the grid:
    for (int iz = 0; iz < nz; ++iz){
        z[iz] = ((doublereal)iz)*dz;
    }
    
    gasflow.setupGrid(nz, &z[0]);

    // specify the objects to use to compute kinetic rates and transport properties
    std::unique_ptr<Transport> trunityLe(newTransportMgr("UnityLewis", &gas));

    gasflow.setTransport(*trunityLe);
    gasflow.setKinetics(gas);
    gasflow.setPressure(pressure);

    /******Create the inlet******/
    Inlet1D inlet;
    inlet.setMoleFractions(x.data());
    doublereal mdot = uin*rho_in;
    inlet.setMdot(mdot);
    cloud.setMpdot(mdot);
    inlet.setTemperature(temp);
    /******Create the outlet******/
    Outlet1D outlet;
    /******Create the domain******/
    std::vector<Domain1D*> domains { &inlet, &gasflow, &outlet };
    gasflow.setupCloud(cloud);
    cloud.GasFlow(&gasflow);
    std::cout << "Initial grid point number = "
            << gasflow.grid().size() << std::endl;
    
    Sim1D sprayflame(domains);
    sprayflame.Cloud(cloud);
    sprayflame.SprStFlow(gasflow);
    /******Supply the initial guess******/
    vector_fp locs{0.0, 0.3, 0.7, 1.0};
    vector_fp value;
    double uout = inlet.mdot()/rho_out;//TODO: mdot at inlet can not fulfill the continuity
    value = {uin, uin, uout, uout};
    sprayflame.setInitialGuess("u", locs, value);
    value = {temp, temp, Tad, Tad};
    sprayflame.setInitialGuess("T", locs, value);
    
    for(int k = 0; k < nsp; ++k){
        value = {yin[k], yin[k], yout[k], yout[k]};
        sprayflame.setInitialGuess(gas.speciesName(k), locs, value);
    }

    inlet.setMoleFractions(x.data());
    inlet.setMdot(mdot);
    inlet.setTemperature(temp);

    sprayflame.showSolution();
    int flowdomain = 1;
    int loglevel = 1;

    /**********Solve freely propagating flame with spray cloud**********/
    sprayflame.setFixedTemperature(0.5 * (temp + Tad));
    gasflow.solveEnergyEqn();
    sprayflame.solve(loglevel, refine_grid, do_spray); 

    double flameSpeed_mix = sprayflame.value(flowdomain,gasflow.componentIndex("u"),0);
    print("Flame speed with mixture-averaged transport: {} m/s\n",
            flameSpeed_mix);
    const vector_fp& solution = sprayflame.solutionVector();

    vector_fp zvec,Tvec,COvec,CO2vec,O2vec,N2vec,ARvec,C2H5OHvec,Uvec;

    print("\n{:9s}\t{:8s}\t{:5s}\t{:7s}\n",
            "z (m)", "T (K)", "U (m/s)", "Y(CO)");
    for (size_t n = 0; n < gasflow.nPoints(); n++) {
        Tvec.push_back(sprayflame.value(flowdomain,gasflow.componentIndex("T"),n));
        C2H5OHvec.push_back(sprayflame.value(flowdomain,gasflow.componentIndex("C2H5OH"),n));
        O2vec.push_back(sprayflame.value(flowdomain,gasflow.componentIndex("O2"),n));
        N2vec.push_back(sprayflame.value(flowdomain,gasflow.componentIndex("N2"),n));
        ARvec.push_back(sprayflame.value(flowdomain,gasflow.componentIndex("AR"),n));
        COvec.push_back(sprayflame.value(flowdomain,gasflow.componentIndex("CO"),n));
        CO2vec.push_back(sprayflame.value(flowdomain,gasflow.componentIndex("CO2"),n));
        Uvec.push_back(sprayflame.value(flowdomain,gasflow.componentIndex("u"),n));
        zvec.push_back(gasflow.grid(n));
        print("{:9.6f}\t{:8.3f}\t{:5.3f}\t{:7.5f}\t{:7.5f}\t{:7.5f}\t{:7.5f}\t{:7.5f}\n",
                gasflow.grid(n), Tvec[n], Uvec[n], C2H5OHvec[n], O2vec[n], N2vec[n], ARvec[n], COvec[n]);
    }

    print("\nAdiabatic flame temperature from equilibrium is: {}\n", Tad);
    print("Flame speed for phi={} is {} m/s.\n", phi, Uvec[0]);

    std::ofstream outfile("./result/nospray.csv", std::ios::trunc);
    outfile << "  Grid,   Temperature,   Uvec,  C2H5OH, O2, N2, AR,   CO,    CO2\n";
    for (size_t n = 0; n < gasflow.nPoints(); n++) {
        print(outfile, " {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}\n",
                gasflow.grid(n), Tvec[n], Uvec[n], C2H5OHvec[n], O2vec[n], N2vec[n], ARvec[n], COvec[n], CO2vec[n]);
    }

    cloud.write();

    return Uvec[0];
}

int main()
{
    doublereal phi_over;
    print("Enter overall Phi : ");
    std::cin >> phi_over;

    doublereal flamespeed_Old;
    doublereal flamespeed_New;
    bool do_spray;//false; // do_spray = 1: two-way couple; do_spray = 0: pure gas;
    
    double rsd = 1.0;
    int Nloop = 1;
    do{
        if(Nloop == 1){
            do_spray = false;
            flamespeed_Old = 0.3;// pesudo flame speed
            flamespeed_Old = SprayFreeFlame(flamespeed_Old, do_spray);
        }
        else{
            do_spray = true;
            flamespeed_New = SprayFreeFlame(flamespeed_Old, do_spray);

            rsd = abs(flamespeed_New - flamespeed_Old)/(flamespeed_Old + 1e-14);
        }

        Nloop++;
    }while(rsd > 1e-4);

    std::cout << "The spray flame speed is \t" << flamespeed_New << std::endl;
    return 0;
}
