//! @file Lagrangian.h

#ifndef CT_LAGRANGIAN_H
#define CT_LAGRANGIAN_H

#include <iostream>
#include <ostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "cantera/base/global.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "SprayStFlow.h"
#include "Liquid.h"
#include "AllSpecies.h"

namespace Cantera
{
class StFlow;

class Lagrangian
{
    public:
        Lagrangian(IdealGasPhase* ph,
                   const doublereal parcelDiameter,
                   const doublereal TInjection,
                   const doublereal pinjection,
                   const doublereal injectionPosition,
                   const doublereal phi_gas,
                   const doublereal phi_over,
                   const doublereal fa_st_mass,
                   const size_t fuelIndex);

        //set up the injection properties:
        void setupInjection(doublereal d, doublereal Tp, doublereal mdotp);
        
        //set up gas-phase information:
        void evalGasFlow(const vector_fp& solution);
        
        //interpolating the gas flow field into parcel's position:
        void evalParcelFlow();
        
        //set up lagrangian field:
        void GasFlow(StFlow* field){
            gas = field;
        }
                
        //evaluate parcel's transport quantities:
        void evalTransf();
        
        //tracking the parcel:
        void solve();

        //initialize the injection:
        void inject(doublereal injectionPosition);
        
        //clear the field at previous iteration step:
        void clearParcel();
        
        //clear the field of gas-phase information in the last step:
        void clearGasFlow(bool do_spray);

        void recordOldValue();

        //interpolate data from Eulerian field to parcel's position:
        //@ field: Eulerian field i.e. u, T, rho
        //@ grid: Eulerian grid
        //@ xp: lagrangian parcel's position
        doublereal intpfield(vector_fp& field, vector_fp&grid, doublereal xp) const;

        //exterpolate date from lagrangian parcel's position to Eulerian grid point:
        doublereal extpfield(doublereal value, vector_fp& grid, doublereal xp) const;
        
        //Calculate the residual during iteration:
        bool evalRsd(const size_t& Nloop, const vector_fp& solution);
        bool evalResidual(const size_t& Nloop, const bool ifAddSpraySource, const vector_fp& solution);
        //set the mass flux of liquid droplets:
        //@ mdot: gas phase mass flux.
        void setMpdot(doublereal mdot);

        const doublereal getmp(size_t n) const{
            return mp[n];
        }

        //Calculate the dmp/dt using the liquid evaporation model:
        //@ n: parcel's index
        doublereal mddot(size_t n);
        
        //Calculate the dmp/dt using the theoritcal model:
        //@ n: parcel's index
        doublereal mddot_th(size_t n);
        
        //Calculate the dTp/dt using the liquid evaporation model:
        //@ n: parcel's index
        doublereal Tddot(size_t n);

        doublereal hTransRate(size_t n);

        doublereal Ygas(size_t k, size_t j) const
        {
            return Ygtmp[k][j];
        }

        //mass transfer rate for one parcel at position xp_n:
        doublereal mtfp(size_t n)
        {
            return Nd*mddot(n);
        }

        //heat transfer rate for one parcel at position xp_n:
        doublereal htfp(size_t n)
        {
            return Nd*hTransRate(n);
        }

        doublereal stfp(size_t k, size_t n)
        {
            if(k==kf){
                return -Nd*mddot(n);
            }
            else{
                return Nd*mddot(n);
            }
        }

        //mass transfer rate at grid points j [kg/s]:
        doublereal mtf(size_t j) const
        {
            return aa*mtf_[j] + (1-aa)*mtfOld_[j];
        }

        doublereal htf(size_t j) const
        {
            return aa*htf_[j] + (1-aa)*htfOld_[j];
        }

        doublereal stf(size_t k, size_t j) const
        {
            return aa*stf_[k][j] + (1-aa)*stfOld_[k][j];
        }

        void setFuel(const std::vector<std::string>& fuelName)
        {
            fuelName_ = fuelName;
        }

        size_t fuelIndex(size_t k)
        {
            return fuelIndex_[k];
        }

        std::vector<size_t> fuelIndex() const
        {
            return fuelIndex_;
        }

        //Drag force model:
        doublereal Cd(doublereal Red) const
        {
            Red = (Red > small ? Red : small);
            if(Red>1000){
                return 0.424;
            }
            else if(Red<=1000 && Red>0)
            {
                double Cd_ = (24.0/Red)*(1 + (1.0/6.0)*std::pow(Red, 0.66666));
                return Cd_;
            }
            else 
            {
                std::cout << "######## ERROR: Failed Reynolds defination! ########\n" << std::endl;
            }
        }

        void write() const;
    
    private:
        StFlow* gas;
        IdealGasPhase* Thermo;
        // Ethanol fuel;
        NC7H16 fuel;

        doublereal small;
        // //evaluate residual:
        vector_fp TOld_;
        vector_fp TNew_;
        
        //gas phase:
        vector_fp z;
        vector_fp ug;
        vector_fp vg;
        vector_fp mug;
        vector_fp rhog;
        vector_fp Tg;
        vector_fp cpg;
        vector_fp kappag;
        vector_fp mw;
        std::vector<std::string> fuelName_;
        std::vector<size_t> fuelIndex_;
        std::vector<vector_fp> Yg;
        std::vector<vector_fp> Ygtmp;
        doublereal Sl;//laminar flame speed
        doublereal rho_inlet;
        doublereal umax;//max gas velocity
        doublereal fa_st_m;
        doublereal Phi_over;
        doublereal Phi_gas;
        //parcel:
        vector_fp dp;
        vector_fp mp;
        vector_fp up;
        vector_fp Tp;
        vector_fp xp;
        vector_fp rhop;
        //gas phase interpolated to the parcel's position:
        vector_fp rho_;
        vector_fp mu_;
        vector_fp u_;
        vector_fp T_;
        vector_fp cp_;
        vector_fp kappa_;
        std::vector<std::vector<doublereal>> Y_;

        size_t Np; //parcel number.
        size_t Nd; //droplet number per parcel.
        size_t kf; // fuel index in solution vector
        doublereal p0;
        doublereal d_inj; //injection diameter.
        doublereal Mdotp_inj; //injection parcel's mass flow rate.
        doublereal z_inj; // distance between the liquid inlet and gas inlet
        doublereal Tp_inj; //injection temperature.
        doublereal p_inj; //injection pressure.
        doublereal Vold_inj; //injection droplet volume.
        doublereal Md_inj; //injection droplet mass. 
        doublereal dtlag; //evaporation time step.

        //transport quantities of parcel:
        vector_fp mtfd_;//mass transfer of droplet
        vector_fp htfd_;//mass transfer of droplet
        vector_fp mtfp_;//mass transfer of parcel
        vector_fp htfp_;//heat transfer of parcel 
        std::vector<std::vector<doublereal>> stfp_;//species transfer of parcel

        //feedback from lagrangian field to Eulerian field:
        vector_fp mtf_; //mass transfer
        vector_fp htf_; //heat transfer
        std::vector<std::vector<doublereal>> stf_; //species transfer

        //source term in the last iteration step:
        vector_fp mtfOld_;
        vector_fp htfOld_;
        std::vector<std::vector<doublereal>> stfOld_;

        //some const:
        const doublereal RR = 8314.0;
        const doublereal aa = 0.1; //relaxation factor
        const doublereal Co = 0.1; //Corrent number of lagrangian
};
}

#endif