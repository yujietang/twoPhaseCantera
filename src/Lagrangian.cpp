#include "Lagrangian.h"
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
    doublereal Rhop_inj = rhold(Tp_inj);
    doublereal Vd_inj = Pi*std::pow(d, 3.0)/6.0;
    Md_inj = rhold(Tp)*Vd_inj;
    Mdotp_inj = Nd*Md_inj/dtlag;
}

//double curr_dp, curr_md, curr_up, curr_Td, curr_xp, curr_rhop;
void Lagrangian::RK4_time_advancing()
{
    //xp_curr = xp[n];
    //up_curr = up[n];
    //if(!inertia_parcel) up_curr = intpfield(ug, z, xp_curr);
    // d\ xp/dt = up;  d\ up/dt = acceleration
    double xp_hold = this -> curr_xp, up_hold = this -> curr_up, md_hold = this -> curr_md, Td_hold = this -> curr_Td;
    //double xp_carry = curr_xp, up_carry = curr_up, md_carry = curr_md, Td_carry = curr_Td;
    double xp_res = xp_hold, up_res = up_hold, md_res = md_hold, Td_res = Td_hold;
    double dxp_dt, dup_dt, dmd_dt, dTd_dt;

    //double diameter_p, density_p;
    //density_p = rhold(Tp[n]);
    //diameter_p = dp[n]; //assuming fixed parcel density and diameter during one time step increment, can be further refined

    //RK4 parameters
    int Nstages = 4;
    vector_fp Alpha{0, 0.5, 0.5, 1}; 
    vector_fp Beta{1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
    double alpha, beta;

    for(int j=0; j<Nstages; j++){ // Begin RK4
        dxp_dt = calc_dxp_dt();
        dup_dt = calc_dup_dt();
        dmd_dt = calc_dmd_dt();
        dTd_dt = calc_dTd_dt();
        //dup_dt = parcel_acc(xp_carry, up_carry, diameter_p, density_p);

        if(j != Nstages-1){
            alpha = Alpha[j+1];
            //xp_carry = xp_hold + alpha * dtlag * dxp_dt;
            //up_carry = up_hold + alpha * dtlag * dup_dt;
            //md_carry = md_hold + alpha * dtlag * dmd_dt;
            //Td_carry = Td_hold + alpha * dtlag * dTd_dt;
            this -> curr_xp = xp_hold + alpha * dtlag * dxp_dt;
            this -> curr_up = up_hold + alpha * dtlag * dup_dt;
            this -> curr_md = md_hold + alpha * dtlag * dmd_dt;
            this -> curr_Td = Td_hold + alpha * dtlag * dTd_dt;
            if(!inertia_parcel) this -> curr_up = intpfield(ug, z, this -> curr_xp);
            this -> update_curr();
        }

        beta = Beta[j];

        xp_res += beta * dtlag * dxp_dt;
        up_res += beta * dtlag * dup_dt;
        md_res += beta * dtlag * dmd_dt;
        Td_res += beta * dtlag * dTd_dt;
        if(!inertia_parcel) up_res = intpfield(ug, z, xp_res);
    }

    this -> curr_xp = xp_res;
    this -> curr_up = up_res;
    this -> curr_md = md_res;
    this -> curr_Td = Td_res;
    this -> update_curr();
    return;
}

void Lagrangian::update_curr()
{
    // updated curr_xp, curr_up, curr_md, curr_Td
    this -> curr_rhop = rhold(this -> curr_Td);
    this -> curr_dp = pow(6.0*this->curr_md/(Pi*this->curr_rhop+small), 0.33333);
    return;
}

double Lagrangian::calc_dxp_dt()
{
    return this -> curr_up;
}

double Lagrangian::calc_dup_dt()
{
    double acc = 0.0;
    double _rhog, _Red, _ug, _mug;

    _rhog = intpfield(rhog, z, this -> curr_xp);
    _ug = intpfield(ug, z, this -> curr_xp);
    _mug = intpfield(mug, z, this -> curr_xp);

    _Red = _rhog * std::abs(this -> curr_up - _ug) * this -> curr_dp / (_mug + small);
    acc = 0.75 * _rhog * Cd(_Red) * std::abs(this -> curr_up - _ug) * (_ug - this -> curr_up) /
        (this -> curr_rhop * this -> curr_dp + small);
    return acc;
}

double Lagrangian::calc_dmd_dt()
{
    const double RR = 8314;
    const double MWf = 46.0;//only for ethanol
    //all parameters needed:
    double Ni;
    double kc;
    double Cs;
    double Cinf;
    double psat;
    double Ts;
    double Sh;
    double Red;
    double Sc;
    double rp;
    double D;
    
    double pc = p0;
    double Tinf = intpfield(Tg, z, this -> curr_xp);
    double uG = intpfield(ug, z, this -> curr_xp);
    double muG = intpfield(mug, z, this -> curr_xp);
    double rhoG = intpfield(rhog, z, this -> curr_xp);
    vector_fp Yc(gas->nsp(), 0.0);
    for(size_t k=0; k<Yc.size(); ++k){
        Yc[k] = intpfield(Yg[k], z, this -> curr_xp);
    }
    //1/3 rule to calc the surface temperature:
    Ts = (2*this->curr_Td + Tinf)/3.0;

    //calculate the surface values:
    double TRatio = Tinf/Ts;
    double mus = muG/TRatio;
    double rhos = rhoG*TRatio;
    // saturation pressure for species i [pa]
    // - carrier phase pressure assumed equal to the liquid vapour pressure
    //   close to the surface
    // NOTE: if pSat > pc then particle is superheated
    // calculated evaporation rate will be greater than that of a particle
    // at boiling point, but this is not a boiling model
    psat = pv(Tinf);
    // std::cout << "\nsaturated pressure of fuel pSat = " << psat << std::endl;
    //vapour diffusivity [m2/s]:
    // D = Dab(pc, Ts, 28.0);
    D = 1.5e-7;
    // std::cout << "\nvapour diffusivity of ethanol D = " << D << std::endl;
    //Schmidt number:
    Sc = mus/(rhos*D + small);
    // Sc = 0.75;

    //droplet Reynold's number:
    Red = rhos*std::abs(uG - this->curr_up)*this->curr_dp/mus;
    
    //Ranz-Marshall:
    Sh = 2 + 0.522*pow(Red, 2.0)*pow(Sc, 0.33333);

    //species volume fraction in the carrier gas:
    //the ethanol index is k = 30:
    vector_fp Xc_(Yc.size(), 0.0);
    for(size_t k=0; k<Yc.size(); ++k){
        Xc_[k] = Yc[k]/mw[k];
    }
    double Xc = 0;
    for(size_t k=0; k<Xc_.size(); ++k){
        Xc += Xc_[k];
    }
    Xc = Xc_[30]/Xc;
    // std::cout << "\n$$$ ethanol in carrier phase Xc = " << Xc << std::endl;  
    //vapor concentration at surface [kmol/m3] at film temperature:
    Cs = psat/(RR*Ts);

    //vapor concentration in bulk gas [kmol/m3] at film temperature:
    Cinf = Xc*pc/(RR*Ts);
    // std::cout << "\ncarrier phase vapour pressure = " << Xc*pc << std::endl;
    //mass transfer coefficient [m/s]:
    kc = Sh*D/(this -> curr_dp + small);

    //molar flux of vapour [kmol/m2/s]:
    Ni = std::max(kc*(Cs - Cinf), 0.0);
    
    //Return the mass transfer [kg/s]:
    double massTranfRate = -Pi*Ni*this->curr_dp*this->curr_dp*MWf;
    
    return massTranfRate;
}

double Lagrangian::calc_dTd_dt()
{
    //For droplet:
    double md = this->curr_md/Nd;
    double Td = this->curr_Td;
    //For carrier phase:
    double cp = intpfield(cpg, z, this->curr_xp);
    double TG = intpfield(Tg, z, this->curr_xp);
    double rhoG = intpfield(rhog, z, this->curr_xp);
    double uG = intpfield(ug, z, this->curr_xp);
    double muG = intpfield(mug, z, this->curr_xp);
    double mddot_ = this -> calc_dmd_dt();
    
    double Ts = (2*Td+TG)/3;

    double TRatio = TG/Ts;
    double k = intpfield(kappag, z, this->curr_xp);
    double ks = 0.5*(kappav(Ts) + k/TRatio);
    // doublereal ks = 0.2;
    double mus = muG/TRatio;
    std::cout << "\nvapour viscosity mus = " << mus << std::endl;
    double rhos = rhoG*TRatio;
    double cps = cpv(Ts);
    std::cout << "\nvapour heat capacity cps = " << cps << std::endl;
    double Red = rhos*std::abs(uG - this->curr_up)*this->curr_dp/mus;
    double Pr = std::max(mus*(cps/ks), small);
    // const doublereal Pr = 0.75;
    double Nu = 2.0 + 0.6*std::pow(Red, 0.5)*std::pow(Pr, 0.33333);
    // std::cout << "\n$$$ liquid heat capacity cld = " << cld(Td) << std::endl;
    // std::cout << "\n$$$ surface temperature Ts = " << Ts << std::endl;
    double tddot = -(mddot_/(md*cld(Td)))*Lv(Td) + (Pi*this->curr_dp*Nu*ks/(md*cld(Td)))*(TG-Td);
    
    return tddot;
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
//TODO: refined mesh tansport quantities evaluation is needed.
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

//void Lagrangian::parcel_time_advancing(double & xp_curr, double & up_curr, int n, bool inertia_parcel)
//{
//    xp_curr = xp[n];
//    up_curr = up[n];
//    if(!inertia_parcel) up_curr = intpfield(ug, z, xp_curr);
//    // d\ xp/dt = up;  d\ up/dt = acceleration
//    double xp_hold = xp_curr, up_hold = up_curr;
//    double xp_carry = xp_curr, up_carry = up_curr;
//    double dxp_dt, dup_dt;
//
//    double diameter_p, density_p;
//    density_p = rhold(Tp[n]);
//    diameter_p = dp[n]; //assuming fixed parcel density and diameter during one time step increment, can be further refined
//
//    //RK4 parameters
//    int Nstages = 4;
//    vector_fp Alpha{0, 0.5, 0.5, 1}; 
//    vector_fp Beta{1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
//    double alpha, beta;
//
//    for(int j=0; j<Nstages; j++){ // Begin RK4
//        dxp_dt = up_carry;
//        dup_dt = parcel_acc(xp_carry, up_carry, diameter_p, density_p);
//
//        if(j != Nstages-1){
//            alpha = Alpha[j+1];
//            xp_carry = xp_hold + alpha * dtlag * dxp_dt;
//            up_carry = up_hold + alpha * dtlag * dup_dt;
//            if(!inertia_parcel) up_carry = intpfield(ug, z, xp_carry);
//        }
//
//        beta = Beta[j];
//        xp_curr += beta * dtlag * dxp_dt;
//        up_curr += beta * dtlag * dup_dt;
//        if(!inertia_parcel) up_curr = intpfield(ug, z, xp_curr);
//    }
//
//}

//double Lagrangian::parcel_acc(double & _xp, double & _up, const double & diameter_p, const double & density_p)
//{
//    double acc = 0.0;
//    double _rhog, _Red, _ug, _mug;
//
//    _rhog = intpfield(rhog, z, _xp);
//    _ug = intpfield(ug, z, _xp);
//    _mug = intpfield(mug, z, _xp);
//
//    _Red = _rhog * std::abs(_up - _ug) * diameter_p / (_mug + small);
//    acc = 0.75 * _rhog * Cd(_Red) * std::abs(_up - _ug) * (_ug - _up) /
//        (density_p * diameter_p + small);
//    return acc;
//}

void Lagrangian::solve()
{
    clearParcel();
    this->inject();

    doublereal _xp = 0;
    doublereal _up = 0;
    doublereal _ug = ug[0];
    doublereal _rhog = rhog[0];
    doublereal _mug = mug[0];
    doublereal _Red;
    bool inertia_parcel = false;
    //tracking step n:
    size_t marchingStep = 8000;

    //double curr_dp, curr_md, curr_up, curr_Td, curr_xp, curr_rhop;
    //this -> curr_xp = xp_hold + alpha * dtlag * dxp_dt;
    //this -> curr_up = up_hold + alpha * dtlag * dup_dt;
    //this -> curr_md = md_hold + alpha * dtlag * dmd_dt;
    //this -> curr_Td = Td_hold + alpha * dtlag * dTd_dt;
    //if(!inertia_parcel) this -> curr_up = intpfield(ug, z, this -> curr_xp);
    //this -> update_curr();
    {
        this -> curr_xp = _xp;
        this -> curr_up = _up;
        this -> curr_md = mp[0];
        this -> curr_Td = Tp[0];
        if(!inertia_parcel) this -> curr_up = intpfield(ug, z, this -> curr_xp);
        this -> update_curr();
    }

    for(size_t n=1; n<marchingStep; ++n)
    {
        std::cout << "**************** Tracking the parcel [ "
                    << n << " ] ****************\n" << std::endl;

        //parcel_time_advancing(_xp, _up, n-1, inertia_parcel);
        RK4_time_advancing();

        //parcel's position:
        //xp.push_back(_xp);
        //up.push_back(_up);
        xp.push_back(this -> curr_xp);
        up.push_back(this -> curr_up);
        std::cout << "Parcel's position = " << xp[xp.size()-1] << "\n" << std::endl;
        std::cout << "Parcel's velocity = " << up[up.size()-1] << "\n" << std::endl;

        evalParcelFlow();

        //parcel's mass:
        mp.push_back(
            //mp[n-1] + Nd*dtlag*mddot(n-1)
            this -> curr_md
        );

        //parcel's temperature:
        Tp.push_back(
            //Tp[n-1] + Nd*dtlag*Tddot(n-1)
            this -> curr_Td
        );

        //parcel's density
        rhop.push_back(
            //rhold(Tp[n])
            this -> curr_rhop
        ); 

        //parcel's diameter 
        dp.push_back(
            //pow(6*mp[n]/(Pi*rhop[n]+small), 0.33333)
            this -> curr_dp
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
            std::cerr << "### Error: parcel's number is not equal to Np! ###";
        }

        if(dp[n]<small || mp[n]<small)
        {
            break;
        }

        std::cout << "******************** End Tracking ********************"
                << "\n" << std::endl;
    }
    for(size_t ip = 0; ip < mp.size(); ++ip)
    {
        std::cout << "\n******************| Parcel [" 
                    << ip
                    << "]: xp = "
                    << xp[ip] << " |******************\n"
                    << std::endl;
        std::cout << "\n****** Parcel [" 
                    << ip
                    << "]: MASS = "
                    << mp[ip] << "\n"
                    << std::endl;
        std::cout << "****** Parcel [" 
                    << ip
                    << "]: DIAMETER = "
                    << dp[ip] << "\n"
                    << std::endl;
        std::cout << "****** Parcel [" 
                    << ip
                    << "]: rho = "
                    << rhop[ip] << "\n"
                    << std::endl;             
        std::cout << "****** Parcel [" 
                    << ip
                    << "]: Tp = "
                    << Tp[ip] << "\n"
                    << std::endl;        
        std::cout << "****** Parcel [" 
                    << ip
                    << "]: up = "
                    << up[ip] << "\n"
                    << std::endl;  
        std::cout << "****** Parcel [" 
                    << ip
                    << "]: mtf = "
                    << mtfp_[ip] << "\n"
                    << std::endl;         
        std::cout << "****** Parcel [" 
                    << ip
                    << "]: htf = "
                    << htfp_[ip] << "\n"
                    << std::endl;         
        std::cout << "****** Droplet [" 
                    << ip
                    << "]: mtf = "
                    << mtfd_[ip] << "\n"
                    << std::endl;         
        std::cout << "****** Droplet [" 
                    << ip
                    << "]: htf = "
                    << htfd_[ip] << "\n"
                    << std::endl;     
    }
}

void Lagrangian::_solve()
{
    clearParcel();
    this->inject();

    doublereal _xp = 0;
    doublereal _ug = ug[0];
    doublereal _rhog = rhog[0];
    doublereal _mug = mug[0];
    doublereal _Red;
    //tracking step n:
    size_t marchingStep = 8000;
    for(size_t n=1; n<marchingStep; ++n)
    {
        std::cout << "**************** Tracking the parcel [ "
                    << n << " ] ****************\n" << std::endl;
        _xp += up[n-1]*dtlag;
        //parcel's position:
        xp.push_back(_xp);
        std::cout << "Parcel's position = " << xp[xp.size()-1] << "\n" << std::endl;
        
        evalParcelFlow();

        //parcel's mass:
        mp.push_back(
            // mp[n-1]
            mp[n-1] + Nd*dtlag*mddot(n-1)
        );

        //parcel's temperature:
        Tp.push_back(
            // Tp[n-1]
            Tp[n-1] + Nd*dtlag*Tddot(n-1)
        );

        //parcel's density
        rhop.push_back(
            rhold(Tp[n])
        ); 

        //parcel's diameter 
        dp.push_back(
            pow(6*mp[n]/(Pi*rhop[n]+small), 0.33333)
        );

        //parcel's velocity
        //Yuan & Chen (1976):
        _ug = intpfield(ug, z, xp[n-1]);
        std::cout << "LOCAL GAS VELOCITY AT x_p = " 
                << xp[n-1] 
                << " = "
                << _ug << "\n"<< std::endl;
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
        std::cout << "relative velocity is :" << std::abs(_ug-up[n-1])*(_ug-up[n-1]) << std::endl;
        std::cout << "The drag force is :" << 0.75*Cd(_Red)*_rhog*std::abs(_ug-up[n-1])*(_ug-up[n-1])/(dp[n-1]*rhop[n-1]+small) << std::endl;
        up.push_back(
            // up[n-1] + 0.75*Cd(_Red)*_rhog*std::abs(_ug-up[n-1])*(_ug-up[n-1])/(dp[n-1]*rhop[n-1]+small)
            _ug
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
            std::cerr << "### Error: parcel's number is not equal to Np! ###";
        }

        if(dp[n]<small || mp[n]<small)
        {
            break;
        }

        std::cout << "******************** End Tracking ********************"
                << "\n" << std::endl;
    }
    for(size_t ip = 0; ip < mp.size(); ++ip)
    {
        std::cout << "\n******************| Parcel [" 
                    << ip
                    << "]: xp = "
                    << xp[ip] << " |******************\n"
                    << std::endl;
        std::cout << "\n****** Parcel [" 
                    << ip
                    << "]: MASS = "
                    << mp[ip] << "\n"
                    << std::endl;
        std::cout << "****** Parcel [" 
                    << ip
                    << "]: DIAMETER = "
                    << dp[ip] << "\n"
                    << std::endl;
        std::cout << "****** Parcel [" 
                    << ip
                    << "]: rho = "
                    << rhop[ip] << "\n"
                    << std::endl;             
        std::cout << "****** Parcel [" 
                    << ip
                    << "]: Tp = "
                    << Tp[ip] << "\n"
                    << std::endl;        
        std::cout << "****** Parcel [" 
                    << ip
                    << "]: up = "
                    << up[ip] << "\n"
                    << std::endl;  
        std::cout << "****** Parcel [" 
                    << ip
                    << "]: mtf = "
                    << mtfp_[ip] << "\n"
                    << std::endl;         
        std::cout << "****** Parcel [" 
                    << ip
                    << "]: htf = "
                    << htfp_[ip] << "\n"
                    << std::endl;         
        std::cout << "****** Droplet [" 
                    << ip
                    << "]: mtf = "
                    << mtfd_[ip] << "\n"
                    << std::endl;         
        std::cout << "****** Droplet [" 
                    << ip
                    << "]: htf = "
                    << htfd_[ip] << "\n"
                    << std::endl;     
    }
}

void Lagrangian::inject()
{
    xp.push_back(z[0]);
    mp.push_back(Nd*Md_inj);//TODO: evaluate the parcels' and droplets' quantities
    Tp.push_back(Tp_inj);
    up.push_back(ug[0]);
    dp.push_back(d_inj);
    rhop.push_back(rhold(Tp_inj));

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
    // std::cout << "**********************" << std::endl;
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
    const doublereal RR = 8314;
    const doublereal MWf = 46.0;//only for ethanol
    //all parameters needed:
    doublereal Ni;
    doublereal kc;
    doublereal Cs;
    doublereal Cinf;
    doublereal psat;
    doublereal Ts;
    doublereal Sh;
    doublereal Red;
    doublereal Sc;
    doublereal rp;
    doublereal D;
    
    doublereal pc = p0;
    doublereal Tinf = T_[n];
    doublereal uG = u_[n];
    doublereal muG = mu_[n];
    doublereal rhoG = rho_[n];
    vector_fp Yc(Y_.size(), 0.0);
    for(size_t k=0; k<Yc.size(); ++k){
        Yc[k] = Y_[k][n];
    }
    //1/3 rule to calc the surface temperature:
    Ts = (2*Tp[n] + Tinf)/3;

    //calculate the surface values:
    doublereal TRatio = Tinf/Ts;
    doublereal mus = muG/TRatio;
    doublereal rhos = rhoG*TRatio;
    // std::cout << "T_surf = " << Ts << std::endl;
    // saturation pressure for species i [pa]
    // - carrier phase pressure assumed equal to the liquid vapour pressure
    //   close to the surface
    // NOTE: if pSat > pc then particle is superheated
    // calculated evaporation rate will be greater than that of a particle
    // at boiling point, but this is not a boiling model
    psat = pv(Tinf);
    // std::cout << "\nsaturated pressure of fuel pSat = " << psat << std::endl;
    //vapour diffusivity [m2/s]:
    // D = Dab(pc, Ts, 28.0);
    D = 1.5e-7;
    // std::cout << "\nvapour diffusivity of ethanol D = " << D << std::endl;
    //Schmidt number:
    Sc = mus/(rhos*D + small);
    // Sc = 0.75;

    //droplet Reynold's number:
    Red = rhos*std::abs(uG - up[n])*dp[n]/mus;
    
    //Ranz-Marshall:
    Sh = 2 + 0.522*pow(Red, 2.0)*pow(Sc, 0.33333);

    //species volume fraction in the carrier gas:
    //the ethanol index is k = 30:
    vector_fp Xc_(Yc.size(), 0.0);
    for(size_t k=0; k<Yc.size(); ++k){
        Xc_[k] = Yc[k]/mw[k];
    }
    double Xc = 0;
    for(size_t k=0; k<Xc_.size(); ++k){
        Xc += Xc_[k];
    }
    Xc = Xc_[30]/Xc;
    // std::cout << "\n$$$ ethanol in carrier phase Xc = " << Xc << std::endl;  
    //vapor concentration at surface [kmol/m3] at film temperature:
    Cs = psat/(RR*Ts);

    //vapor concentration in bulk gas [kmol/m3] at film temperature:
    Cinf = Xc*pc/(RR*Ts);
    // std::cout << "\ncarrier phase vapour pressure = " << Xc*pc << std::endl;
    //mass transfer coefficient [m/s]:
    kc = Sh*D/(dp[n] + small);

    //molar flux of vapour [kmol/m2/s]:
    Ni = std::max(kc*(Cs - Cinf), 0.0);
    
    //Return the mass transfer [kg/s]:
    doublereal massTranfRate = -Pi*Ni*dp[n]*dp[n]*MWf;
    
    return massTranfRate;
}

doublereal Lagrangian::Tddot(size_t n)
{
    //For droplet:
    doublereal md = mp[n]/Nd;
    doublereal Td = Tp[n];
    //For carrier phase:
    doublereal cp = cp_[n];
    doublereal TG = T_[n];
    doublereal rhoG = rho_[n];
    doublereal uG = u_[n];
    doublereal muG = mu_[n];
    doublereal mddot_ = mddot(n);
    
    doublereal Ts = (2*Td+TG)/3;

    doublereal TRatio = TG/Ts;
    doublereal k = kappa_[n];
    doublereal ks = 0.5*(kappav(Ts) + k/TRatio);
    // doublereal ks = 0.2;
    doublereal mus = muG/TRatio;
    std::cout << "\nvapour viscosity mus = " << mus << std::endl;
    doublereal rhos = rhoG*TRatio;
    doublereal cps = cpv(Ts);
    std::cout << "\nvapour heat capacity cps = " << cps << std::endl;
    doublereal Red = rhos*std::abs(uG - up[n])*dp[n]/mus;
    doublereal Pr = std::max(mus*(cps/ks), small);
    // const doublereal Pr = 0.75;
    doublereal Nu = 2.0 + 0.6*std::pow(Red, 0.5)*std::pow(Pr, 0.33333);
    // std::cout << "\n$$$ liquid heat capacity cld = " << cld(Td) << std::endl;
    // std::cout << "\n$$$ surface temperature Ts = " << Ts << std::endl;
    doublereal tddot = -(mddot_/(md*cld(Td)))*Lv(Td) + (Pi*dp[n]*Nu*ks/(md*cld(Td)))*(TG-Td);
    
    return tddot;
}

void Lagrangian::write() const
{
    double t = 0;
    std::ofstream fout1("test.csv");
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
