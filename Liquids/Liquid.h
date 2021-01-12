#ifndef CTF_LIQUID_H_
#define CTF_LIQUID_H_

#include <string>
#include <vector>
#include <cmath>
#include <iostream>

class Liquid
{
public:
    Liquid() {
        std::cout << "Liquid Constructed" << std::endl;
    }
    virtual ~Liquid() = default;

    virtual double Wf() const = 0;

    virtual std::string name() const = 0;

    virtual double Tc() const = 0;

    virtual double Tb() const = 0;

    // Liquid density [kg/m3]
    virtual double rho(const double& T) const = 0;

    // Latent heat of vapourisation [J/kg]
    virtual double Lv(const double& T) const = 0;

    // Liquid heat capacity [J/kg K]
    virtual double cp(const double& T) const = 0;

    // Ideal gas heat capacity [J/kg K]
    virtual double cpg(const double& T) const = 0;

    // Vapour diffussivity [m2/s]
    virtual double D(const double& T) const = 0;

    // Vapor pressure [Pa]
    virtual double pv(const double& T) const  = 0;

    // Liquid thermal conductivity [W/m K]
    virtual double kappa(const double& T) const = 0;

    // Liquid viscousity [Pa s]
    virtual double mu(const double& T) const = 0;
};

#endif  // CTF_LIQUID_H_
