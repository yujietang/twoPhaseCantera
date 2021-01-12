#ifndef CTF_ICETANE_H_
#define CTF_ICETANE_H_

#include "Liquid.h"

class ICetane : public Liquid
{
public:
    ICetane(const double p0);

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
        const double a(72.782730504);
        const double b(2.7740e-001);
        const double c(6.9200e+002);
        const double d(2.8571e-001);
        return a/std::pow(b, 1 + std::pow(1 - T/c, d));
    }

    // Latent heat of vapourisation [J/kg]
    virtual double Lv(const double& T) const override {
        const double Tc(692.0);
        const double a(362403.13158559485);
        const double b(4.5022e-001);
        const double c(0.0);
        const double d(0.0);
        const double e(0.0);
        const double Tr = T/Tc;
        return a*std::pow(1 - Tr, ((e*Tr + d)*Tr + c)*Tr + b);
    }

    // Liquid heat capacity [J/kg K]
    virtual double cp(const double& T) const override {
        const double a(1239.792051976407);
        const double b(0.7303882862305975);
        const double c(0.0066887121248253405);
        const double d(0.0);
        const double e(0.0);
        const double f(0.0);
        return ((((f*T + e)*T + d)*T + c)*T + b)*T + a;
    }

    // Ideal gas heat capacity [J/kg K]
    virtual double cpg(const double& T) const override {
        const double a(1235.1109250436757);
        const double b(3838.347438540336);
        const double c(1.7206e+003);
        const double d(2912.6325068052984);
        const double e(7.8550e+002);
        return a + b*((c/T)/std::sinh(c/T))*((c/T)/std::sinh(c/T))+ d*((e/T)/std::cosh(e/T))*((e/T)/std::cosh(e/T));
    }

    // Vapour diffussivity [m2/s]
    virtual double D(const double& T) const override {
        const double a(147.18);
        const double b(20.1);
        const double wf(226.4412);
        const double wa(28.0);
        const double alpha(std::sqrt(1/wf + 1/wa));
        const double beta((std::cbrt(a) + std::cbrt(b))*(std::cbrt(a) + std::cbrt(b)));
        return 3.6059e-3*(std::pow(1.8*T, 1.75))*alpha/(p0_*beta);
    }

    // Vapor pressure [Pa]
    virtual double pv(const double& T) const override {
        const double a(1.2831e+002);
        const double b(-1.1839e+004);
        const double c(-1.5333e+001);
        const double d(6.9676e-006);
        const double e(2.0);
        return std::exp(a + b/T + c*std::log(T) + d*std::pow(T, e));
    }

    // Liquid thermal conductivity [W/m K]
    virtual double kappa(const double& T) const override {
        const double a(1.4278e-001);
        const double b(-1.4671e-004);
        const double c(0.0);
        const double d(0.0);
        const double e(0.0);
        const double f(0.0);
        return ((((f*T + e)*T + d)*T + c)*T + b)*T + a;
    }

    // Liquid viscousity [Pa s]
    virtual double mu(const double& T) const override {
        const double a(-2.4660e+001);
        const double b(2.4888e+003);
        const double c(1.8631e+000);
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

#endif  // CTF_ICETANE_H_
