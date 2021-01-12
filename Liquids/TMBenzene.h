#ifndef CTF_TMBENZENE_H_
#define CTF_TMBENZENE_H_

#include "Liquid.h"

class TMBenzene : public Liquid
{
public:
    TMBenzene(const double p0);

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
        const double a(7.25690829e+01);
        const double b(2.59485174e-01);
        const double c(6.49128824e+02);
        const double d(2.77244022e-01);
        return a/std::pow(b, 1 + std::pow(1 - T/c, d));
    }

    // Latent heat of vapourisation [J/kg]
    virtual double Lv(const double& T) const override {
        const double Tc(649.13);
        const double a(4.81080524e+05);
        const double b(3.41999853e-01);
        const double c(0.0);
        const double d(0.0);
        const double e(0.0);
        const double Tr = T/Tc;
        return a*std::pow(1 - Tr, ((e*Tr + d)*Tr + c)*Tr + b);
    }

    // Liquid heat capacity [J/kg K]
    virtual double cp(const double& T) const override {
        const double a(2066.83235476064);
        const double b(-8.14664481609707);
        const double c(0.0322581695445024);
        const double d(-3.01223125427334e-05);
        const double e(0.0);
        const double f(0.0);
        return ((((f*T + e)*T + d)*T + c)*T + b)*T + a;
    }

    // Ideal gas heat capacity [J/kg K]
    virtual double cpg(const double& T) const override {
        const double a(630.989461803106);
        const double b(3107.19440856947);
        const double c(1440.6);
        const double d(2059.88647833212);
        const double e(-650.43);
        return a + b*((c/T)/std::sinh(c/T))*((c/T)/std::sinh(c/T))+ d*((e/T)/std::cosh(e/T))*((e/T)/std::cosh(e/T));
    }

    // Vapour diffussivity [m2/s]
    virtual double D(const double& T) const override {
        const double a(147.18);
        const double b(20.1);
        const double wf(120.19);
        const double wa(28.0);
        const double alpha(std::sqrt(1/wf + 1/wa));
        const double beta((std::cbrt(a) + std::cbrt(b))*(std::cbrt(a) + std::cbrt(b)));
        return 3.6059e-3*(std::pow(1.8*T, 1.75))*alpha/(p0_*beta);
    }

    // Vapor pressure [Pa]
    virtual double pv(const double& T) const override {
        const double a(7.25177292e+01);
        const double b(-7.62519716e+03);
        const double c(-7.28935794e+00);
        const double d(3.28917e-06);
        const double e(2.0);
        return std::exp(a + b/T + c*std::log(T) + d*std::pow(T, e));
    }

    // Liquid thermal conductivity [W/m K]
    virtual double kappa(const double& T) const override {
        const double a(0.2043);
        const double b(-0.000239);
        const double c(0.0);
        const double d(0.0);
        const double e(0.0);
        const double f(0.0);
        return ((((f*T + e)*T + d)*T + c)*T + b)*T + a;
    }

    // Liquid viscousity [Pa s]
    virtual double mu(const double& T) const override {
        const double a(-13.362);
        const double b(1183.0);
        const double c(0.333);
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

#endif  // CTF_TMBENZENE_H_
