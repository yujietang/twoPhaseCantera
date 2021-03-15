#ifndef CTF_NC7H16_H_
#define CTF_NC7H16_H_

#include "Liquid.h"

class NC7H16 : public Liquid
{
public:
    NC7H16(const double p0);

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
        const double a(61.3838396836);
        const double b(0.26211);
        const double c(540.2);
        const double d(0.28141);
        return a/std::pow(b, 1 + std::pow(1 - T/c, d));
    }

    // Latent heat of vapourisation [J/kg]
    virtual double Lv(const double& T) const override {
        const double Tc(540.20);
        const double a(499121.791545248);
        const double b(0.38795);
        const double c(0.0);
        const double d(0.0);
        const double e(0.0);
        const double Tr = T/Tc;
        return a*std::pow(1 - Tr, ((e*Tr + d)*Tr + c)*Tr + b);
    }

    // Liquid heat capacity [J/kg K]
    virtual double cp(const double& T) const override {
        const double Tc(540.20);
        const double a(6.11976102401216);
        const double b(3137.69909384855);
        const double c(182.274175063868);
        const double d(-254.530511150515);
        const double t = 1.0 - (std::min(T, Tc)/Tc);
        return a*a/t + b - t*(2.0*a*c + t*(a*d + t*(c*c/3.0 + t*(0.5*c*d + 0.2*d*d*t))));
    }

    // Ideal gas heat capacity [J/kg K]
    virtual double cpg(const double& T) const override {
        const double a(1199.05392998284);
        const double b(3992.85457666361);
        const double c(1676.6);
        const double d(2734.42177956968);
        const double e(756.4);
        return a + b*((c/T)/std::sinh(c/T))*((c/T)/std::sinh(c/T))+ d*((e/T)/std::cosh(e/T))*((e/T)/std::cosh(e/T));
    }

    // Vapour diffussivity [m2/s]
    virtual double D(const double& T) const override {
        const double a(147.18);
        const double b(20.1);
        const double wf(100.204);
        const double wa(28);
        const double alpha(std::sqrt(1/wf + 1/wa));
        const double beta((std::cbrt(a) + std::cbrt(b))*(std::cbrt(a) + std::cbrt(b)));
        return 3.6059e-3*(std::pow(1.8*T, 1.75))*alpha/(p0_*beta);
    }

    // Vapor pressure [Pa]
    virtual double pv(const double& T) const override {
        const double a(87.829);
        const double b(-6996.4);
        const double c(-9.8802);
        const double d(7.2099e-6);
        const double e(2);
        return std::exp(a + b/T + c*std::log(T) + d*std::pow(T, e));
    }

    // Liquid thermal conductivity [W/m K]
    virtual double kappa(const double& T) const override {
        const double a(0.215);
        const double b(-0.000303);
        const double c(0.0);
        const double d(0.0);
        const double e(0.0);
        const double f(0.0);
        return ((((f*T + e)*T + d)*T + c)*T + b)*T + a;
    }

    // Liquid viscousity [Pa s]
    virtual double mu(const double& T) const override {
        const double a(-24.451);
        const double b(1533.1);
        const double c(2.0087);
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

#endif  // CTF_NC7H16_H_
