//! @file Lagrangian.h

#ifndef CT_LAGRANGIAN_H
#define CT_LAGRANGIAN_H

#include <iostream>
#include <ostream>
#include <fstream>
#include <string>
#include "cantera/base/global.h"
#include "SprayStFlow.h"

namespace Cantera
{
class StFlow;

class Lagrangian
{
    public:
        Lagrangian(const doublereal parcelDiameter,
                   const doublereal TInjection,
                   const doublereal pinjection,
                   const doublereal lagrangianTimeStep,
                   const size_t dropletNumber);

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
        void _solve();
        void solve();

        //initialize the injection:
        void inject();
        
        //clear the field at previous iteration step:
        void clearParcel();
        
        //clear the field of gas-phase information in the last step:
        void clearGasFlow();
        
        //interpolate data from Eulerian field to parcel's position:
        //@ field: Eulerian field i.e. u, T, rho
        //@ grid: Eulerian grid
        //@ xp: lagrangian parcel's position
        doublereal intpfield(vector_fp& field, vector_fp&grid, doublereal xp) const;

        //exterpolate date from lagrangian parcel's position to Eulerian grid point:
        doublereal extpfield(doublereal value, vector_fp& grid, doublereal xp) const;
        
        //Calculate the residual during iteration:
        bool evalRsd(const vector_fp& solution);

        //set the mass flux of liquid droplets:
        //@ mdot: gas phase mass flux.
        void setMpdot(doublereal mdot);

        const doublereal getmp(size_t n) const{
            return mp[n];
        }

        //Calculate the dmp/dt using the liquid evaporation model:
        //@ n: parcel's index
        doublereal mddot(size_t n);
        
        //Calculate the dTp/dt using the liquid evaporation model:
        //@ n: parcel's index
        doublereal Tddot(size_t n);

        //mass transfer rate for one parcel at position xp_n:
        doublereal mtfp(size_t n)
        {
            return Nd*mddot(n);
        }

        //heat transfer rate for one parcel at position xp_n:
        doublereal htfp(size_t n) 
        {
            return getmp(n)*cp_[n]*Tddot(n);
        }

        //mass transfer rate at grid points j [kg/s]:
        doublereal mtf(size_t j) const
        {
            return mtf_[j];
        }

        doublereal htf(size_t j) const
        {
            return htf_[j];
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

        /******liquid and vapour properties are calculated by NSRD functions******/
        //liquid thermal conductivity:
        doublereal kappa(doublereal Tp) const
        {
            return -0.000281*Tp + 0.253;
        }

        //vapour thermal conductivity:
        doublereal kappav(doublereal Ts) const
        {
            double a = -3.12;
            double b = 0.7152;
            double c = -3550000;
            double d = 0.0;
            return a*pow(Ts, b)/(1+(c/Ts)+d/(Ts*Ts));
        }

        //specific liquid heat capacity:
        doublereal cld(doublereal Tp) const
        {   
            return ((((0.0*Tp + 0.0)*Tp + 5.20523562482363e-5)*Tp 
                    + 0.00714146172046278)*Tp -1.21990926653498)*Tp +  2052.57332;
        }

        //liquid density:
        doublereal rhold(doublereal Tp) const
        {   
            return 70.1308387/std::pow(0.26395, 1 
                    + std::pow(1 - Tp/516.25, 0.2367));
        }
        //liquid vapour dynamic viscousity:
        // doublereal muv(doublereal Tp)

        //latent heat:
        doublereal Lv(doublereal Tp) const
        {
            doublereal Tr = Tp/516.25;
            return  958345.091059064*std::pow(1 - Tr, ((0.0*Tr + 0.0)*Tr + 0.75362)*Tr - 0.4134);
        }

        //@ pc: carrier phase pressure
        //@ Tv: vapor(surface) temperature
        //binary diffusivity derived by API Function:
        doublereal Dab(doublereal pc, doublereal Tv, doublereal Wa)
        {
            const doublereal beta = std::pow((cbrt(147.18)+cbrt(20.1)), 2.0);
            const doublereal alphaBinary = sqrt(1.0/46.069 + 1.0/Wa);
            return (3.6059e-3)*(std::pow(1.8*Tv, 1.75))*alphaBinary/(pc*beta+small);
        }

        //vapor pressure
        doublereal pv(doublereal Tv) const
        {
            return exp(59.796 - 6595/Tv - 5.0474*log(Tv) + 6.3e-7*pow(Tv, 2.0));
        }

        //vapour heat capacity:
        doublereal cpv(doublereal Tv)
        {
            doublereal a1, a2, a3, a4, a5, a6, a7; 
            if(Tv >= 300 && Tv <= 1000){
                a1 = 4.34717120;
                a2 = 1.86288e-2;
                a3 = -6.779467e-6;
                a4 = 8.16592600e-10;
                a5 = 0.0;
                a6 = -3.06615743e+4;
                a7 = 3.24247304;
                return a1*std::pow(Tv, -2.0) 
                        + a2*std::pow(Tv, -1.0) 
                        + a3 
                        + a4*Tv
                        + a5*std::pow(Tv, 2.0)
                        + a6*std::pow(Tv, 3.0)
                        + a7*std::pow(Tv, 4.0);
            }
            else if(Tv >1000 && Tv < 5000){
                a1 = 5.765358e-1;
                a2 = 2.894512e-2;
                a3 = -1.61001e-5;
                a4 = 3.591641e-9;
                a5 = 0.0;
                a6 = -2.963595e+4;
                a7 = 2.270813e+1;
                return a1*std::pow(Tv, -2.0) 
                        + a2*std::pow(Tv, -1.0) 
                        + a3 
                        + a4*Tv
                        + a5*std::pow(Tv, 2.0)
                        + a6*std::pow(Tv, 3.0)
                        + a7*std::pow(Tv, 4.0);            
            }
            else
            {
                std::cerr << "Temperature range overflow!";
            }
            
        }

        void parcel_time_advancing(double & xp_curr, double & up_curr, int n, bool inertia_parcel);
        double parcel_acc(double & _xp, double & _up, const double & diameter_p, const double & density_p);

        void write() const;
    private:
        StFlow* gas;

        doublereal small;

        //evaluate residual:
        vector_fp Told;
        vector_fp Tnew;
        
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
        doublereal p0;
        doublereal d_inj; //injection diameter.
        doublereal Mdotp_inj; //injection parcel's mass flow rate.
        doublereal Tp_inj; //injection temperature.
        doublereal p_inj; //injection pressure.
        doublereal Vold_inj; //injection droplet volume.
        doublereal Md_inj; //injection droplet mass. 
        //TODO: variable evaporation time step:
        doublereal dtlag; //evaporation time step.

        doublereal Ndot; //number of parcels per second. [1/s]

        //transport quantities of parcel:
        vector_fp mtfd_;//mass transfer of droplet
        vector_fp htfd_;//mass transfer of droplet
        vector_fp mtfp_;//mass transfer of parcel
        vector_fp htfp_;//heat transfer of parcel 

        //feedback from lagrangian field to Eulerian field:
        vector_fp mtf_; //mass transfer
        vector_fp htf_; //heat transfer
};
}

#endif
