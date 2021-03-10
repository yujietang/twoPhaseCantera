#include "Lagrangian.h"
#include "Liquid.h"
#include "AllSpecies.h"

#include <cmath>

namespace Cantera
{
Lagrangian::Lagrangian(
   const doublereal parcelDiameter,
   const doublereal TInjection,
   const doublereal pinjection,
   const doublereal lagrangianTimeStep,
   const size_t dropletNumber
) :
    fuel(pinjection),
    d_inj(parcelDiameter),
    Tp_inj(TInjection),
    Mdotp_inj(0.0),
    p_inj(pinjection),
    Nd(dropletNumber),
    Np(0),
    dtlag(lagrangianTimeStep),
    gas(0),
    small(1.0e-13)
{
    setupInjection(d_inj, 
                Tp_inj,
                Mdotp_inj);
    std::cout << "====================Parcel Initialization====================" << std::endl;
    std::cout << "injection temperature = " << Tp_inj << std::endl;
    std::cout << "injection pressure = " << p_inj << std::endl;
    std::cout << "liquid mass flow rate = " << Mdotp_inj << std::endl;
    std::cout << "droplet mass = " << Md_inj << std::endl;
    std::cout << "droplet number = " << Nd << std::endl;
    std::cout << "parcel's diameter = " << d_inj << std::endl;
    std::cout << "parcel's mass = " << Nd*Md_inj << std::endl;
}

void Lagrangian::setupInjection(
    doublereal d,
    doublereal Tp,
    doublereal mdotp
)
{
    doublereal Rhop_inj = fuel.rho(Tp_inj);
    doublereal Vd_inj = Pi*std::pow(d, 3.0)/6.0;
    Md_inj = fuel.rho(Tp)*Vd_inj;
    Mdotp_inj = Nd*Md_inj/dtlag;
}

void Lagrangian::evalGasFlow(const vector_fp& solution)
{
    size_t Nvar = c_offset_Y + gas->nsp();
    std::cout << "Number of variables at each grid point: "
            << Nvar << std::endl;
    //pointer to the gas-phase field:
    z = gas->grid();
    p0 = gas->pressure();
    rhog = gas->density();
    mug = gas->viscosity();
    kappag = gas->heatConductivity();
    cpg = gas->cp();
    mw = gas->MW();
    Yg.resize(gas->nsp());  
    for(size_t i = 0; i < solution.size(); ++i)
    {
        if(i%Nvar==0){
            ug.push_back(solution[i]);
        }
        else if(i%Nvar==1){
            vg.push_back(solution[i]);
        }
        else if(i%Nvar==2){
            Tg.push_back(solution[i]);
        }
        else
        {
            for(size_t k=0; k<gas->nsp(); ++k){
                if(i%Nvar==(k+c_offset_Y)){
                    Yg[k].push_back(solution[i]);
                }
            }
            // for(size_t k=0; k<fuelName_.size(); ++k)
            // {
            //     if(i%(gas->nsp()+c_offset_Y)==fuelIndex_[k]){
            //         Yg[k].push_back(solution[i]);
            //     }
            // }
        }
    }
    /******mass and heat transfer******/
    mtf_.resize(z.size(), 0.0);
    htf_.resize(z.size(), 0.0);
    std::cout << "\n========== Check the Eulerian gas field ==========" << std::endl;
    std::cout << "z size = " << z.size() << std::endl;
    std::cout << "rhog size = " << rhog.size() << std::endl;
    std::cout << "mug size = " << mug.size() << std::endl;
    std::cout << "cp size = " << cpg.size() << std::endl;
    std::cout << "kappa size = " << kappag.size() << std::endl;
    std::cout << "mw size = " << mw.size() << std::endl;
    std::cout << "species size = " << fuelIndex_.size() << std::endl;
    std::cout << "Ug size = " << ug.size() << std::endl;
    std::cout << "Vg size = " << vg.size() << std::endl;
    std::cout << "Tg size = " << Tg.size() << std::endl;
    std::cout << "Yg size = " << Yg[0].size() << "\n" << std::endl;
}

void Lagrangian::evalParcelFlow()
{
    rho_.resize(xp.size());
    mu_.resize(xp.size());
    u_.resize(xp.size());
    T_.resize(xp.size());
    cp_.resize(xp.size());
    kappa_.resize(xp.size());
    Y_.resize(gas->nsp());
    for(size_t k = 0; k < Y_.size(); ++k)
    {
        Y_[k].resize(xp.size());
    }

    //interpolate the gas flow field into the parcel's position:
    for(size_t ip = 0; ip < xp.size(); ++ip )
    {
        rho_[ip] = intpfield(rhog, z, xp[ip]);
        mu_[ip] = intpfield(mug, z, xp[ip]);
        u_[ip] = intpfield(ug, z, xp[ip]);
        T_[ip] = intpfield(Tg, z, xp[ip]);
        cp_[ip] = intpfield(cpg, z, xp[ip]);
        kappa_[ip] = intpfield(kappag, z, xp[ip]);
        for(size_t k=0; k<Y_.size(); ++k){
            Y_[k][ip] = intpfield(Yg[k], z, xp[ip]);
        }
    }
}

void Lagrangian::evalTransf()
{
    //resize the vector of quantities transfer:
    htf_.resize(z.size(), 0.0);
    mtf_.resize(z.size(), 0.0);//TODO: multi-component case is a 2d vector.

    doublereal leftz;
    doublereal rightz;
    doublereal leftLength;
    doublereal rightLength;
    doublereal dz;
    
    for(size_t iz = 0; iz < z.size()-1; ++iz)
    {
        dz = z[iz+1] - z[iz];
        leftz = z[iz];
        rightz = z[iz+1];

        for(size_t ip = 1; ip < Np; ++ip)
        {
            if(xp[ip] > leftz && xp[ip] < rightz){
                leftLength = xp[ip] - z[ip];
                rightLength = z[iz+1] - xp[ip];

                mtf_[iz] += mtfp_[ip] 
                            + ((mtfp_[ip] - mtfp_[ip-1])/(xp[ip] - xp[ip-1]))*leftLength;
                mtf_[iz+1] += mtfp_[ip]
                            + ((mtfp_[ip+1] - mtfp_[ip])/(xp[ip+1]-xp[ip]))*rightLength;
                htf_[iz] += htfp_[ip] 
                            + ((htfp_[ip] - htfp_[ip-1])/(xp[ip] - xp[ip-1]))*leftLength;
                htf_[iz+1] += htfp_[ip]
                            + ((htfp_[ip+1] - htfp_[ip])/(xp[ip+1]-xp[ip]))*rightLength;
            }
            else break;
        }
    }
}

void Lagrangian::solve()
{
    clearParcel();
    this->inject();

    doublereal _xp = 0;
    doublereal _ug = ug[0];
    doublereal _Tp;
    doublereal _rhog = rhog[0];
    doublereal _mug = mug[0];
    doublereal _Red;
    //tracking step n:
    size_t marchingStep = 2e+5;

    // size_t marchingStep = 2e+4;
    for(size_t n=1; n<marchingStep; ++n)
    {
        std::cout << "**************** Tracking the parcel [ "
                    << n << " ] ****************\n" << std::endl;
        if((dp.back()<small) || (mp.back()<small) || (_xp > gas->zmax()))
        {
            break;
        }
        _xp += up[n-1]*dtlag;

        //parcel's position:
        xp.push_back(_xp);
        std::cout << "Parcel's position = " << xp[xp.size()-1] << "\n" << std::endl;
        
        evalParcelFlow();

        //parcel's mass:
        mp.push_back(
            mp[n-1] + Nd*dtlag*mddot(n-1)
        );
        std::cout << "Check ... ... mp =                " << mp[mp.size()-1] << std::endl;
        
        //parcel's temperature:
        _Tp = Tp[n-1] + dtlag*Tddot(n-1);
        Tp.push_back(
            (_Tp < fuel.Tb() ? _Tp : fuel.Tb())
        );

        //parcel's density
        rhop.push_back(
            fuel.rho(Tp[n])
        ); 

        //parcel's diameter 
        dp.push_back(
            pow(6*mp[n]/(Pi*rhop[n]+small), 0.33333)
        );

        //parcel's velocity
        //Yuan & Chen (1976):
        _ug = intpfield(ug, z, xp[n-1]);
        std::cout << "LOCAL Parcel VELOCITY AT x_p = " 
                << xp[n-1] 
                << " = "
                << up[n-1] << "\n"<< std::endl;
        _rhog = intpfield(rhog, z, xp[n-1]);
        _mug = intpfield(mug, z, xp[n-1]);
        _Red = _rhog*std::abs(_ug - up[n-1])*dp[n-1]/(_mug+small);
        if(_ug == up[n-1]){
            std::cout << "\n######## Warning: No slip velocity! ########\n" << std::endl;
        }
        std::cout << "Re_d is :" << _Red << "\n" << std::endl;
        std::cout << "Drag coeff. is : " << Cd(_Red) << "\n" << std::endl;
        up.push_back(
            up[n-1] + dtlag*0.75*Cd(_Red)*_rhog*std::abs(_ug-up[n-1])*(_ug-up[n-1])/(dp[n-1]*rhop[n-1]+small)
            // _ug //no drag force
        ); 
        
        /******For droplet******/
        //mass transfer:
        mtfd_.push_back(
            mddot(n)
        );
        //heat transfer:
        htfd_.push_back(
            Tddot(n)
        );
        /******For parcel******/
        //mass transfer:
        mtfp_.push_back(
            mtfp(n)
        );
        //heat transfer:
        htfp_.push_back(
            htfp(n)
        );

        ++Np;

        if(Np != mp.size()){
            std::cerr << "###### Error: parcel's number is not equal to Np! ######";
        }

        std::cout << "******************** End Tracking ********************"
                << "\n" << std::endl;
    }
    // for(size_t ip = 0; ip < mp.size(); ++ip)
    // {
    //     std::cout << "\n******************| Parcel [" 
    //                 << ip
    //                 << "]: xp = "
    //                 << xp[ip] << " |******************\n"
    //                 << std::endl;
    //     std::cout << "\n****** Parcel [" 
    //                 << ip
    //                 << "]: MASS = "
    //                 << mp[ip] << "\n"
    //                 << std::endl;
    //     std::cout << "****** Parcel [" 
    //                 << ip
    //                 << "]: DIAMETER = "
    //                 << dp[ip] << "\n"
    //                 << std::endl;
    //     std::cout << "****** Parcel [" 
    //                 << ip
    //                 << "]: rho = "
    //                 << rhop[ip] << "\n"
    //                 << std::endl;             
    //     std::cout << "****** Parcel [" 
    //                 << ip
    //                 << "]: Tp = "
    //                 << Tp[ip] << "\n"
    //                 << std::endl;        
    //     std::cout << "****** Parcel [" 
    //                 << ip
    //                 << "]: up = "
    //                 << up[ip] << "\n"
    //                 << std::endl;  
    //     std::cout << "****** Parcel [" 
    //                 << ip
    //                 << "]: mtf = "
    //                 << mtfp_[ip] << "\n"
    //                 << std::endl;         
    //     std::cout << "****** Parcel [" 
    //                 << ip
    //                 << "]: htf = "
    //                 << htfp_[ip] << "\n"
    //                 << std::endl;         
    //     std::cout << "****** Droplet [" 
    //                 << ip
    //                 << "]: mtf = "
    //                 << mtfd_[ip] << "\n"
    //                 << std::endl;         
    //     std::cout << "****** Droplet [" 
    //                 << ip
    //                 << "]: htf = "
    //                 << htfd_[ip] << "\n"
    //                 << std::endl;     
    // }
}

void Lagrangian::inject()
{
    xp.push_back(z[0]);
    mp.push_back(Nd*Md_inj);//TODO: evaluate the parcels' and droplets' quantities
    Tp.push_back(Tp_inj);
    up.push_back(ug[0]);
    dp.push_back(d_inj);
    rhop.push_back(fuel.rho(Tp_inj));

    mtfd_.push_back(0.0);
    htfp_.push_back(0.0);
    mtfp_.push_back(0.0);
    htfp_.push_back(0.0);

    Np = 1;
}

void Lagrangian::clearParcel()
{
    //clear parcel:
    dp.clear();
    mp.clear();
    up.clear();
    Tp.clear();
    xp.clear();

    mtfd_.clear();
    htfd_.clear();
    mtfp_.clear();
    htfp_.clear();
    std::cout << "====== Parcels has been cleared up! ======\n"<< std::endl;
}

void Lagrangian::clearGasFlow()
{
    std::cout << "\n================ Starting Clear Gas Flow ================\n"
            << std::endl;
    std::cout << "z size = " << z.size() << std::endl;
    //clear gas flow:
    z.clear();
    ug.clear();
    mug.clear();
    rhog.clear();
    Tg.clear();
    cpg.clear();
    kappag.clear();
    mw.clear();
    Yg.resize(gas->nsp());
    for(size_t k=0; k<Yg.size(); ++k){
        Yg[k].clear();
    }   

    mtf_.clear();
    htf_.clear();
    std::cout << "\n====== Gas Flow has been cleared up! ======"
            << "\n" << std::endl;
}

doublereal Lagrangian::intpfield(
    vector_fp& field,
    vector_fp& grid,
    doublereal xp
) const
{
    if(xp < gas->zmin()){
        std::cerr << "ERROR: Parcel does not reach the inlet!" << std::endl;
    }
    else if(xp > gas->zmax()){
        std::cerr << "ERROR: Parcel moves beyond the domain!" << std::endl;
    }
    else if(xp == gas->zmin()){
        return field[0];
    }
    else if(xp == gas->zmax()){
        return field.back();
    }
    else
    {
        size_t j_L = lower_bound(z.begin(), z.end(), xp) - z.begin() - 1;
        size_t j_R = lower_bound(z.begin(), z.end(), xp) - z.begin();

        doublereal z_L = z[j_L];
        doublereal z_R = z[j_R];

        return field[j_R]*(xp - z_L)/(z_R - z_L) 
            + field[j_L]*(z_R - xp)/(z_R - z_L); 
    }
}

bool Lagrangian::evalRsd(const vector_fp& solution)
{   
    doublereal rsd = 1.0;
    doublereal Phi0 = 0.0;
    doublereal Phi1 = 1.0;

    Told = Tg;
    for(size_t i = 0; i < solution.size(); ++i){
        if(i%(gas->nsp()+c_offset_Y)==2){
            Tnew.push_back(solution[i]);
        }    
    }

    for(size_t ii=0; ii<Told.size(); ++ii){
        Phi0 += sqrt(Told[ii]*Told[ii]);
    }
    for(size_t jj=0; jj<Tnew.size(); ++jj){
        Phi1 += sqrt(Tnew[jj]*Tnew[jj]);
    }

    /*****evaluate the residual*****/
    rsd = (Phi1 -Phi0)/Phi0;

    if(rsd < 1.0e-4){
        return true;
    }
    else{
        return false;
    }
}

void Lagrangian::setMpdot(doublereal mdot)
{
    Mdotp_inj = mdot;
}

doublereal Lagrangian::mddot(size_t n)
{
    const size_t kf = 30;
    const doublereal MWf = mw[kf];//TODO:only for ethanol
    //all parameters needed:
    doublereal Ni;
    doublereal kc;
    doublereal Cs;
    doublereal Cinf;
    doublereal psat;
    doublereal Sh;
    doublereal Red;
    doublereal Sc;
    doublereal rp;
    doublereal Dab;
    
    doublereal pc = p0;
    doublereal Tinf = T_[n];
    doublereal uG = u_[n];
    doublereal muG = mu_[n];
    doublereal rhoG = rho_[n];
    vector_fp Yc(Y_.size(), 0.0);
    for(size_t k=0; k<Yc.size(); ++k){
        Yc[k] = Y_[k][n];
    }

    //calculate the surface (vapour film) values:
    doublereal Ts, rhos, mus, Pr, kappas;
    Ts = (2*Tp[n] + Tinf)/3;
    
    const doublereal TRatio = Tinf/Ts;
    rhos = rhoG*TRatio;
    mus = muG/TRatio;

    // saturation pressure for species i [pa]
    psat = fuel.pv(Tp[n]);

    //vapour diffusivity [m2/s]:
    Dab = fuel.D(Ts);
 
    //Schmidt number:
    Sc = mus/(rhos*Dab + small);

    //droplet Reynold's number:
    Red = rhos*std::abs(uG - up[n])*dp[n]/mus;
    
    //Ranz-Marshall:
    Sh = 2 + 0.522*pow(Red, 2.0)*pow(Sc, 0.33333);

    //species volume fraction in the carrier gas:
    //the fuel index is k = 30:
    vector_fp Xc_(Yc.size(), 0.0);
    for(size_t k=0; k<Yc.size(); ++k){
        Xc_[k] = Yc[k]/mw[k];
    }
    double Xc = 0;
    for(size_t k=0; k<Xc_.size(); ++k){
        Xc += Xc_[k];
    }
    Xc = Xc_[kf]/Xc;

    //vapor concentration at surface [kmol/m3] at film temperature:
    Cs = psat/(RR*Ts);

    //vapor concentration in bulk gas [kmol/m3] at film temperature:
    Cinf = Xc*pc/(RR*Ts);

    //mass transfer coefficient [m/s]:
    kc = Sh*Dab/(dp[n] + small);

    //molar flux of vapour [kmol/m2/s]:
    Ni = std::max(kc*(Cs - Cinf), 0.0);
    
    //Return the mass transfer [kg/s]:
    doublereal massTranfRate = -Pi*Ni*dp[n]*dp[n]*MWf;
    

    return massTranfRate;
}

/***********************Theoretical model of droplet evaporation************************/
doublereal Lagrangian::mddot_th(size_t n)
{
    const size_t kf = 30;
    const doublereal MWf = mw[kf];//TODO:only for ethanol
    //all parameters needed:
    doublereal Ni;
    doublereal kc;
    doublereal Cs;
    doublereal Cinf;
    doublereal psat;
    doublereal Sh;
    doublereal Red;
    doublereal Sc;
    doublereal rp;
    doublereal D;   
    doublereal pc = p0;
    doublereal Tinf = T_[n];
    doublereal Td = Tp[n];
    doublereal uG = u_[n];
    doublereal muG = mu_[n];
    doublereal rhoG = rho_[n];
    vector_fp Yc(Y_.size(), 0.0);
    vector_fp YkMW(Y_.size(), 0.0);
    vector_fp Xc(Yc.size(), 0.0);
    for(size_t k=0; k<Yc.size(); ++k){
        Yc[k] = Y_[k][n];
        YkMW[k] = Yc[k]/mw[k];
    }
    double sumYkMW = std::accumulate(YkMW.begin(), YkMW.end(), 0.0);
    for(size_t k=0; k<Yc.size(); ++k){
        Xc[k] = (Yc[k]/mw[k])/(sumYkMW);
    }
    for(size_t k=0; k<Xc.size(); ++k){
        Xc[k] /= std::accumulate(Xc.begin(), Xc.end(), 0.0);
    }
    //calculate the surface (vapour film) values:
    doublereal Ts, rhos, mus, Pr, kappas;
    Ts = (2*Td + Tinf)/3;
    // saturation pressure for species i [pa]
    psat = fuel.pv(Td);
    //vapour diffusivity [m2/s]:
    // doublereal Dab = fuel.D(Ts);
    doublereal Dab = 10.2e-6;
    //liquid surface fuel molar fraction:
    // Clausius-Clapeyron relation:
    doublereal Tboiling = fuel.Tb();
    doublereal Xsf = (psat/p0)*exp((38.56/(RR))*((1.0/Tboiling)-(1.0/Td)));
    doublereal Ysf = Xsf*mw[kf];
    const double TRatio = Tinf/Ts; 
    rhos = rhoG*TRatio;
    mus = muG/TRatio;
    //Schmidt number:
    Sc = mus/(rhos*Dab + small);
    //droplet Reynold's number:
    Red = rhos*std::abs(uG - up[n])*dp[n]/mus; 
    //Ranz-Marshall:
    Sh = 2 + 0.522*pow(Red, 2.0)*pow(Sc, 0.33333);  
    // Xsf = (psat/p0);
    doublereal Bm = (Xsf - Xc[kf])/(1 - Xsf);
    //Return the mass transfer [kg/s]:
    doublereal massTranfRate = -Pi*dp[n]*Sh*Dab*rhoG*log(1 + Bm);
    return massTranfRate;
}

doublereal Lagrangian::Tddot(size_t n)
{
    //For droplet:
    const size_t kf = 30;
    const doublereal MWf = mw[kf];//TODO:only for ethanol
    doublereal md = mp[n]/Nd;
    doublereal Td = Tp[n];
    //For carrier phase:
    doublereal pc = p0;
    doublereal cpG = cp_[n];
    doublereal Tinf = T_[n];
    doublereal rhoG = rho_[n];
    doublereal uG = u_[n];
    doublereal muG = mu_[n];
    doublereal kG = kappa_[n];
    vector_fp Yc(Y_.size(), 0.0);
    vector_fp YkMW(Y_.size(), 0.0);
    vector_fp Xc(Yc.size(), 0.0);
    for(size_t k=0; k<Yc.size(); ++k){
        Yc[k] = Y_[k][n];
        YkMW[k] = Yc[k]/mw[k];
    }
    double sumYkMW = std::accumulate(YkMW.begin(), YkMW.end(), 0.0);
    for(size_t k=0; k<Yc.size(); ++k){
        Xc[k] = (Yc[k]/mw[k])/(sumYkMW);
    }
    for(size_t k=0; k<Xc.size(); ++k){
        Xc[k] /= std::accumulate(Xc.begin(), Xc.end(), 0.0);
    }

    //surface temperature:
    doublereal Ts = (2*Td + Tinf)/3.0;
    const double TRatio = Tinf/Ts; 
    // saturation pressure for species i [pa]
    doublereal psat = fuel.pv(Ts);

    //liquid surface fuel molar fraction:
    // Clausius-Clapeyron relation:
    doublereal Tboiling = fuel.Tb();
    doublereal Xsf = (psat/p0)*exp((38.56/(RR))*((1.0/Tboiling)-(1.0/Td)));
    doublereal Ysf = Xsf*mw[kf];

    const double sqrtW = sqrt(mw[kf]);
    const double cbrtW = cbrt(mw[kf]);

    doublereal kappas = Ysf*cbrtW*kG;

    doublereal cps = Xsf*cpG;

    // double sumYiSqrtW = Ysf*sqrtW;
    double sumYiCbrtW = Ysf*cbrtW;

    cps = std::max(cps, small);

    kappas /= sumYiCbrtW;
    kappas *= TRatio;
    kappas = std::max(kappas, small);

    doublereal rhos = rhoG*TRatio;

    doublereal mus = muG/TRatio;
    
    doublereal Pr = cps*mus/kappas;

    //mass transfer rate:
    doublereal mddot_ = mddot(n);

    doublereal Red = rhos*std::abs(uG - up[n])*dp[n]/mus;
    
    doublereal Nu = 2.0 + 0.6*std::pow(Red, 0.5)*std::pow(Pr, 0.33333);

    doublereal tddot = (mddot_/(md*fuel.cp(Td)))*fuel.Lv(Td) + (Pi*dp[n]*Nu*kappas/(md*fuel.cp(Td)))*(Tinf-Td);

    return tddot;
}

void Lagrangian::write() const
{
    double t = 0;
    std::ofstream fout1("d=25.csv");
    fout1 << "# t [s], xp [m], d [micron], d^2 [micron^2], mp [mg], Tp [K], Tg [K], ug [m/s]" << std::endl;
    for(size_t ip = 0; ip < xp.size(); ++ip){
        if((ip%5) == 0){
            fout1 << (t+ip*dtlag) << ","
                << xp[ip] << ","
                << dp[ip]*(1e+6) << ","
                << std::pow(dp[ip]*(1e+6), 2.0) << ","
                << mp[ip]*(1e+6) << ","
                << Tp[ip] << ","
                << T_[ip] << ","
                << u_[ip] << std::endl;
        }
    }
}

}
