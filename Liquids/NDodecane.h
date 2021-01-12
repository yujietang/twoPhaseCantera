#ifndef CTF_NDODECANE_H_
#define CTF_NDODECANE_H_

#include "Liquid.h"

class NDodecane : public Liquid
{
public:
    NDodecane(const double p0);

    virtual double Wf() const override {
        return Wf_;
    }

    virtual std::string name() const override {
        return name_;
    }

    virtual double Tc() const override {
        return Tc_;
    }

    virtual double Tb() const override {
        return Tb_;
    }

    // Liquid density [kg/m3]
    virtual double rho(const double& T) const override {
        const double a(60.53982858);
        const double b(0.25511);
        const double c(658.0);
        const double d(0.29368);
        return a/std::pow(b, 1 + std::pow(1 - T/c, d));
    }

    // Latent heat of vapourisation [J/kg]
    virtual double Lv(const double& T) const override {
        const double Tc(658.0);
        const double a(454020.829174935);
        const double b(0.40681);
        const double c(0.0);
        const double d(0.0);
        const double e(0.0);
        const double Tr = T/Tc;
        return a*std::pow(1 - Tr, ((e*Tr + d)*Tr + c)*Tr + b);
    }

    // Liquid heat capacity [J/kg K]
    virtual double cp(const double& T) const override {
        const double a(2983.53861146661);
        const double b(-8.0352006011577);
        const double c(0.018207916025784);
        const double d(0.0);
        const double e(0.0);
        const double f(0.0);
        return ((((f*T + e)*T + d)*T + c)*T + b)*T + a;
    }

    // Ideal gas heat capacity [J/kg K]
    virtual double cpg(const double& T) const override {
        const double a(1250.16144371778);
        const double b(3894.02247296552);
        const double c(1715.5);
        const double d(2650.67101879792);
        const double e(777.5);
        return a + b*((c/T)/std::sinh(c/T))*((c/T)/std::sinh(c/T))+ d*((e/T)/std::cosh(e/T))*((e/T)/std::cosh(e/T));
    }

    // Vapour diffussivity [m2/s]
    virtual double D(const double& T) const override {
        const double a(147.18);
        const double b(20.1);
        const double wf(170.338);
        const double wa(28.0);
        const double alpha(std::sqrt(1/wf + 1/wa));
        const double beta((std::cbrt(a) + std::cbrt(b))*(std::cbrt(a) + std::cbrt(b)));
        return 3.6059e-3*(std::pow(1.8*T, 1.75))*alpha/(p0_*beta);
    }

    // Vapor pressure [Pa]
    virtual double pv(const double& T) const override {
        const double a(137.47);
        const double b(-11976.0);
        const double c(-16.698);
        const double d(8.0906e-06);
        const double e(2.0);
        return std::exp(a + b/T + c*std::log(T) + d*std::pow(T, e));
    }

    // Liquid thermal conductivity [W/m K]
    virtual double kappa(const double& T) const override {
        const double a(0.2047);
        const double b(-0.0002326);
        const double c(0.0);
        const double d(0.0);
        const double e(0.0);
        const double f(0.0);
        return ((((f*T + e)*T + d)*T + c)*T + b)*T + a;
    }

    // Liquid viscousity [Pa s]
    virtual double mu(const double& T) const override {
        const double a(-20.607);
        const double b(1943.0);
        const double c(1.3205);
        const double d(0.0);
        const double e(0.0);
        return std::exp(a + b/T + c*std::log(T) + d*std::pow(T, e));
    }

protected:
    double p0_;
    double Wf_;  // [kg/kmol]
    double Tc_;  // critical temperature [K]
    double Tb_;  // boiling temperature [K]
    std::string name_;
};

#endif  // CTF_NDODECANE_H_
