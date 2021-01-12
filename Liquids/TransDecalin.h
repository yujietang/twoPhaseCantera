#ifndef CTF_TRANSDECALIN_H_
#define CTF_TRANSDECALIN_H_

#include "Liquid.h"

class TransDecalin : public Liquid
{
public:
    TransDecalin(const double p0);

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
        const double a(77.40785469999999);
        const double b(2.6991e-001);
        const double c(6.8705e+002);
        const double d(2.9520e-001);
        return a/std::pow(b, 1 + std::pow(1 - T/c, d));
    }

    // Latent heat of vapourisation [J/kg]
    virtual double Lv(const double& T) const override {
        const double Tc(687.05);
        const double a(381980.86117480276);
        const double b(-1.2000e-001);
        const double c(6.0000e-001);
        const double d(0.0);
        const double e(0.0);
        const double Tr = T/Tc;
        return a*std::pow(1 - Tr, ((e*Tr + d)*Tr + c)*Tr + b);
    }

    // Liquid heat capacity [J/kg K]
    virtual double cp(const double& T) const override {
        const double a(611.1983103440793);
        const double b(2.6704664636572084);
        const double c(0.002811512227582765);
        const double d(0.0);
        const double e(0.0);
        const double f(0.0);
        return ((((f*T + e)*T + d)*T + c)*T + b)*T + a;
    }

    // Ideal gas heat capacity [J/kg K]
    virtual double cpg(const double& T) const override {
        const double a(624.5795751267605);
        const double b(4017.2726812438073);
        const double c(-6.6920e+002);
        const double d(-1773.5600674126422);
        const double e(-7.0630e+002);
        return a + b*((c/T)/std::sinh(c/T))*((c/T)/std::sinh(c/T))+ d*((e/T)/std::cosh(e/T))*((e/T)/std::cosh(e/T));
    }

    // Vapour diffussivity [m2/s]
    virtual double D(const double& T) const override {
        const double a(147.18);
        const double b(20.1);
        const double wf(138.253);
        const double wa(28.0);
        const double alpha(std::sqrt(1/wf + 1/wa));
        const double beta((std::cbrt(a) + std::cbrt(b))*(std::cbrt(a) + std::cbrt(b)));
        return 3.6059e-3*(std::pow(1.8*T, 1.75))*alpha/(p0_*beta);
    }

    // Vapor pressure [Pa]
    virtual double pv(const double& T) const override {
        const double a(1.8012e+002);
        const double b(-1.1582e+004);
        const double c(-2.5078e+001);
        const double d(2.2475e-002);
        const double e(1.0);
        return std::exp(a + b/T + c*std::log(T) + d*std::pow(T, e));
    }

    // Liquid thermal conductivity [W/m K]
    virtual double kappa(const double& T) const override {
        const double a(1.3950e-001);
        const double b(-9.0000e-005);
        const double c(0.0);
        const double d(0.0);
        const double e(0.0);
        const double f(0.0);
        return ((((f*T + e)*T + d)*T + c)*T + b)*T + a;
    }

    // Liquid viscousity [Pa s]
    virtual double mu(const double& T) const override {
        const double a(-2.8692e+001);
        const double b(2.3967e+003);
        const double c(2.5297e+000);
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

#endif  // CTF_TRANSDECALIN_H_
