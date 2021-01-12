#ifndef CTF_IOCTANE_H_
#define CTF_IOCTANE_H_

#include "Liquid.h"

class IOctane : public Liquid
{
public:
    IOctane(const double p0);

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
        const double a(67.2363666);
        const double b(0.27373);
        const double c(543.96);
        const double d(0.2846);
        return a/std::pow(b, 1 + std::pow(1 - T/c, d));
    }

    // Latent heat of vapourisation [J/kg]
    virtual double Lv(const double& T) const override {
        const double Tc(543.96);
        const double a(375379.713037617);
        const double b(0.1549);
        const double c(0.138);
        const double d(0.0666);
        const double e(0.0);
        const double Tr = T/Tc;
        return a*std::pow(1 - Tr, ((e*Tr + d)*Tr + c)*Tr + b);
    }

    // Liquid heat capacity [J/kg K]
    virtual double cp(const double& T) const override {
        const double a(1219.89652546156);
        const double b(1.67205049417409);
        const double c(0.00414073237562483);
        const double d(0.0);
        const double e(0.0);
        const double f(0.0);
        return ((((f*T + e)*T + d)*T + c)*T + b)*T + a;
    }

    // Ideal gas heat capacity [J/kg K]
    virtual double cpg(const double& T) const override {
        const double a(997.10236275617);
        const double b(4627.4653990598);
        const double c(1594.0);
        const double d(2933.52942721328);
        const double e(677.94);
        return a + b*((c/T)/std::sinh(c/T))*((c/T)/std::sinh(c/T))+ d*((e/T)/std::cosh(e/T))*((e/T)/std::cosh(e/T));
    }

    // Vapour diffussivity [m2/s]
    virtual double D(const double& T) const override {
        const double a(147.18);
        const double b(20.1);
        const double wf(114.231);
        const double wa(28.0);
        const double alpha(std::sqrt(1/wf + 1/wa));
        const double beta((std::cbrt(a) + std::cbrt(b))*(std::cbrt(a) + std::cbrt(b)));
        return 3.6059e-3*(std::pow(1.8*T, 1.75))*alpha/(p0_*beta);
    }

    // Vapor pressure [Pa]
    virtual double pv(const double& T) const override {
        const double a(120.81);
        const double b(-7550.0);
        const double c(-16.111);
        const double d(0.017099);
        const double e(1.0);
        return std::exp(a + b/T + c*std::log(T) + d*std::pow(T, e));
    }

    // Liquid thermal conductivity [W/m K]
    virtual double kappa(const double& T) const override {
        const double a(0.1508);
        const double b(-0.0001712);
        const double c(0.0);
        const double d(0.0);
        const double e(0.0);
        const double f(0.0);
        return ((((f*T + e)*T + d)*T + c)*T + b)*T + a;
    }

    // Liquid viscousity [Pa s]
    virtual double mu(const double& T) const override {
        const double a(-15.811);
        const double b(1282.5);
        const double c(0.67791);
        const double d(-3.8617e-28);
        const double e(10.0);
        return std::exp(a + b/T + c*std::log(T) + d*std::pow(T, e));
    }

protected:
    double p0_;
    double Wf_;  // [kg/kmol]
    double Tc_;  // critical temperature [K]
    double Tb_;  // boiling temperature [K]
    std::string name_;
};

#endif  // CTF_IOCTANE_H_
