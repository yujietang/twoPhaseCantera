#include "Inlet1D.h"
#include "cantera/IdealGasMix.h"
#include "cantera/transport.h"
#include "Sim1D.h"
#include "SprayStFlow.h"
#include "Lagrangian.h"

#include <fstream>

using namespace Cantera;
using fmt::print;

doublereal SprayFreeFlame(doublereal flamespeed, doublereal phi_over, bool do_spray)
{
    /************************************************************/
    //                       Injection
    /************************************************************/
    doublereal parcelDiameter(20e-6);            // injection droplet diameter [m]
    doublereal injTemperature(300);              // injection parcel's temperature [K]
    doublereal injPressure(1.0*OneAtm);          // injection pressure [Pa]
    doublereal injectionPosition(0.008);
    const doublereal CoLag = 0.025;               //Lagrangian Co number
    /************************************************************/
    //                          Mesh
    /************************************************************/
    const size_t meshPointNumber = 200;          // Mesh Point Number
    const doublereal domainLength = 0.02;        // Domain Length
    bool refine_grid = false;                    // Refined mesh has been turned off
    /************************************************************/
    //                        Gas Flow
    /************************************************************/
    const doublereal temp = 300.0;              // inlet gas flow temperature [K]
    const doublereal pressure = 1.0*OneAtm;     // inlet gas flow pressure [atm]
    const doublereal overallEqRatio = 1.0;      // overall equivalence ratio
    const doublereal phi = 0.8;                 // gas flow equivalence ratio
    /************************************************************/
    //                        Chemistry(for ethanol)
    /************************************************************/
    // IdealGasMix gas("Ethanol_57.cti", "gas");   //for single component of ethanol fuel:
    // const size_t fuelIndex = 48;                //fuel index in ck file
    // const doublereal C_atoms = 2.0;                   // C2H5OH
    // const doublereal H_atoms = 6.0;                   // C2H5OH
    // const doublereal O_atoms = 1.0;                   // C2H5OH
    // size_t nsp = gas.nSpecies();                // number of species
    // vector_fp x(nsp, 0.0);
    // doublereal ax = C_atoms + H_atoms/4.0 - O_atoms/2.0; // air consumption
    // doublereal fa_stoic_mole = 1.0 / (4.77 * ax);        // fuel / air mole ratio at stoic state
    // doublereal fa_stoic_mass = 0.1113;  
    // // Species' mole fraction when the condition is fuel/oxidizer flow:
    // x[gas.speciesIndex("C2H5OH")] = 1.0; 
    // x[gas.speciesIndex("O2")] = 0.21 / phi / fa_stoic_mole;
    // x[gas.speciesIndex("N2")] = 0.78 / phi / fa_stoic_mole;
    // x[gas.speciesIndex("Ar")] = 0.01 / phi / fa_stoic_mole;
    // // when the condition is pure air flow  (1 mole air):
    // // x[gas.speciesIndex("O2")] = 0.21;
    // // x[gas.speciesIndex("N2")] = 0.78;
    // // x[gas.speciesIndex("AR")] = 0.01;
    /************************************************************/
    //                        Chemistry(for N-heptane)
    /************************************************************/
    IdealGasMix gas("nheptane106.cti", "gas");   //for single component of ethanol fuel:
    // IdealGasMix gas("nheptane.cti", "gas");   //for single component of ethanol fuel: 
    const size_t fuelIndex = 45;                //fuel index in ck file
    const doublereal C_atoms = 7.0;                   // C2H5OH
    const doublereal H_atoms = 16.0;                   // C2H5OH
    const doublereal O_atoms = 0.0;                   // C2H5OH
    size_t nsp = gas.nSpecies();                // number of species
    vector_fp x(nsp, 0.0);
    doublereal ax = C_atoms + H_atoms/4.0 - O_atoms/2.0; // air consumption
    doublereal fa_stoic_mole = 1.0 / (4.76 * ax);        // fuel / air mole ratio at stoic state
    doublereal fa_stoic_mass = 0.065;  
    // Species' mole fraction when the condition is fuel/oxidizer flow:
    x[gas.speciesIndex("NC7H16")] = 1.0; 
    x[gas.speciesIndex("O2")] = 0.21 / phi / fa_stoic_mole;
    x[gas.speciesIndex("N2")] = 0.78 / phi / fa_stoic_mole;
    x[gas.speciesIndex("Ar")] = 0.01 / phi / fa_stoic_mole;
    // when the condition is pure air flow  (1 mole air):
    // x[gas.speciesIndex("O2")] = 0.21;
    // x[gas.speciesIndex("N2")] = 0.78;
    // x[gas.speciesIndex("AR")] = 0.01;
    /************************************************************/

    /************************************************************/
    doublereal uin = 0.3; // initial guess of inlet velocity [m/s]
    gas.setState_TPX(temp, pressure, x.data());
    doublereal rho_in = gas.density();
    std::cout << "inlet rhog = " << rho_in << std::endl;
    vector_fp yin(nsp);
    gas.getMassFractions(&yin[0]);

    /***************Evaluate the liquid mass flux***************/
    doublereal uf = flamespeed; //laminar flame speed;
    doublereal Mdot_gas = rho_in*uf;
    doublereal Mdot_air = Mdot_gas/(fa_stoic_mass*phi+1);
    doublereal Mdot_fuel = Mdot_gas - Mdot_air;
    doublereal Mdot_liquid = fa_stoic_mass*phi_over*Mdot_air - Mdot_fuel;

    std::cout << "The intermediate flame speed is " <<  uf << std::endl; 
    std::cout << "gas mass flux = " << Mdot_gas << std::endl;
    std::cout << "fuel vapour mass flux = " << Mdot_fuel << std::endl;
    std::cout << "liquid mass flux = " << Mdot_liquid << std::endl;
    /***********************************************************/

    /*******************Evaluate the tracking time scale******************/
    doublereal dtlag = CoLag*(domainLength/meshPointNumber)/uf;
    std::cout << "The lagrangian time step is " << dtlag << std::endl;
    /*********************************************************************/
    gas.equilibrate("HP");//evaluate the adiabatic flame temperature.
    vector_fp yout(nsp);
    gas.getMassFractions(&yout[0]);
    doublereal rho_out = gas.density();
    doublereal Tad = gas.temperature();
    print("phi = {}, Tad = {}\n", phi, Tad);

    /**********Lagrangian Cloud**********/
    Lagrangian cloud(
        &gas,
        parcelDiameter,
        injTemperature,
        injPressure,
        dtlag,
        Mdot_liquid,
        injectionPosition,
        fuelIndex);

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
    print("Flame speed with UnityLe transport: {} m/s\n",
            flameSpeed_mix);
    const vector_fp& solution = sprayflame.solutionVector();

    // vector_fp zvec,Tvec,COvec,CO2vec,O2vec,N2vec,ARvec,C2H5OHvec,Uvec;

    // print("\n{:9s}\t{:8s}\t{:5s}\t{:7s}\n",
    //         "z (m)", "T (K)", "U (m/s)", "Y(CO)");
    // for (size_t n = 0; n < gasflow.nPoints(); n++) {
    //     Tvec.push_back(sprayflame.value(flowdomain,gasflow.componentIndex("T"),n));
    //     C2H5OHvec.push_back(sprayflame.value(flowdomain,gasflow.componentIndex("C2H5OH"),n));
    //     O2vec.push_back(sprayflame.value(flowdomain,gasflow.componentIndex("O2"),n));
    //     N2vec.push_back(sprayflame.value(flowdomain,gasflow.componentIndex("N2"),n));
    //     ARvec.push_back(sprayflame.value(flowdomain,gasflow.componentIndex("AR"),n));
    //     COvec.push_back(sprayflame.value(flowdomain,gasflow.componentIndex("CO"),n));
    //     CO2vec.push_back(sprayflame.value(flowdomain,gasflow.componentIndex("CO2"),n));
    //     Uvec.push_back(sprayflame.value(flowdomain,gasflow.componentIndex("u"),n));
    //     zvec.push_back(gasflow.grid(n));
    //     print("{:9.6f}\t{:8.3f}\t{:5.3f}\t{:7.5f}\t{:7.5f}\t{:7.5f}\t{:7.5f}\t{:7.5f}\n",
    //             gasflow.grid(n), Tvec[n], Uvec[n], C2H5OHvec[n], O2vec[n], N2vec[n], ARvec[n], COvec[n]);
    // }

    // print("\nAdiabatic flame temperature from equilibrium is: {}\n", Tad);
    // print("Flame speed for phi={} is {} m/s.\n", phi, Uvec[0]);

    // std::ofstream outfile("./result/linj0.006.csv", std::ios::trunc);

    // outfile << "  Grid,   Temperature,   Uvec,  C2H5OH, O2, N2, AR,   CO,    CO2\n";
    // for (size_t n = 0; n < gasflow.nPoints(); n++) {
    //     print(outfile, " {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}\n",
    //             gasflow.grid(n), Tvec[n], Uvec[n], C2H5OHvec[n], O2vec[n], N2vec[n], ARvec[n], COvec[n], CO2vec[n]);
    // }

    /*for C7H16 flame*/
    vector_fp zvec,Tvec,CO2vec,O2vec,N2vec,ARvec,C7H16vec,Uvec;

    for (size_t n = 0; n < gasflow.nPoints(); n++) {
        Tvec.push_back(sprayflame.value(flowdomain,gasflow.componentIndex("T"),n));
        C7H16vec.push_back(sprayflame.value(flowdomain,gasflow.componentIndex("NC7H16"),n));
        O2vec.push_back(sprayflame.value(flowdomain,gasflow.componentIndex("O2"),n));
        CO2vec.push_back(sprayflame.value(flowdomain,gasflow.componentIndex("CO2"),n));
        Uvec.push_back(sprayflame.value(flowdomain,gasflow.componentIndex("u"),n));
        zvec.push_back(gasflow.grid(n));
    }

    print("\nAdiabatic flame temperature from equilibrium is: {}\n", Tad);
    print("Flame speed for phi={} is {} m/s.\n", phi, Uvec[0]);

    std::ofstream outfile("./result/nheptane_test.csv", std::ios::trunc);

    outfile << "  Grid,   Temperature,   Uvec,  C7H16, O2, CO2\n";
    for (size_t n = 0; n < gasflow.nPoints(); n++) {
        print(outfile, " {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}\n",
                gasflow.grid(n), Tvec[n], Uvec[n], C7H16vec[n], O2vec[n], CO2vec[n]);
    }


    cloud.write();
    return flameSpeed_mix;
}

int main()
{
    doublereal phi_over = 0.9;

    // print("Enter overall Phi : ");
    // std::cin >> phi_over;

    // doublereal flamespeed_Old;
    // doublereal flamespeed_New;
    // bool do_spray;//false; // do_spray = 1: two-way couple; do_spray = 0: pure gas;
    
    // double rsd = 1.0;
    // flamespeed_Old = SprayFreeFlame(1.0, phi_over, false);// initial no spray flame speed calculation
    // std::cout << "\nno spray flame speed =\t" << flamespeed_Old << std::endl; 
    // do{
    //     do_spray = true;
    //     flamespeed_New = SprayFreeFlame(flamespeed_Old, phi_over, do_spray);
    //     rsd = abs(flamespeed_New - flamespeed_Old)/(flamespeed_Old + 1e-14);
    //     flamespeed_Old = flamespeed_New;  
    // }while(rsd > 1e-4);

    doublereal flamespeed = 0.3;

    flamespeed = SprayFreeFlame(flamespeed, phi_over, true);

    std::cout << "\nSpray flame speed is \t" << flamespeed << std::endl;
    return 0;
}
