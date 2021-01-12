#ifndef CTF_MIXTURE_H_
#define CTF_MIXTURE_H_

#include "Liquid.h"
#include <Eigen/Dense>
#include <memory>
#include <cmath>

class Mixture
{
public:
    Mixture(std::vector<std::shared_ptr<Liquid>>&);
    Mixture& operator=(const Mixture&) = delete;

    std::shared_ptr<Liquid> ptr(int i) const {
        return components_[i];
    }

    double TcMin() const {
        return TcMin_;
    }

    double TbMin() const {
        return TbMin_;
    }

    std::string name(int i) const {
        return names_[i];
    }

    void updateProp(Eigen::VectorXd& Td, Eigen::VectorXd& rhod,
        Eigen::VectorXd& kappad, Eigen::VectorXd& cL,
        Eigen::VectorXd& mud, std::vector<Eigen::VectorXd>& Yd) const;

    // Mixture density [kg/m3]
    double rho(const double T, std::vector<double>& Yd) const {
        double res = 0.0;
        std::vector<double> Vf(numComponents_, 0.0);  // volume fractions
        double sumVf = 0.0;
        for (int n = 0; n < numComponents_; n++) {
            Vf[n] = Yd[n] / components_[n]->rho(T);
            sumVf += Vf[n];
        }
        for (int n = 0; n < numComponents_; n++) {
            Vf[n] /= sumVf;
        }
        for (int i = 0; i < numComponents_; i++) {
            res += Vf[i] * components_[i]->rho(T);
        }
        return res;
    }

    // Latent heat of vapourisation [J/kg]
    double Lv(const double T, std::vector<double>& Yd) const {
        double res = 0.0;
        for (int i = 0; i < numComponents_; i++) {
            res += Yd[i] * components_[i]->Lv(T);
        }
        return res;
    }

    // Mixture heat capacity [J/kg K]
    double cp(const double T, std::vector<double>& Yd) const {
        double res = 0.0;
        for (int i = 0; i < numComponents_; i++) {
            res += Yd[i] * components_[i]->cp(T);
        }
        return res;
    }

    // Ideal gas heat capacity [J/kg K]
    double cpg(const double T, std::vector<double>& Yd) const {
        double res = 0.0;
        for (int i = 0; i < numComponents_; i++) {
            res += Yd[i] * components_[i]->cpg(T);
        }
        return res;
    }

    // Vapour diffussivity [m2/s]
    double D(const double T, std::vector<double>& Yd) const {
        double res = 0.0;
        for (int i = 0; i < numComponents_; i++) {
            res += Yd[i] * components_[i]->D(T);
        }
        return res;
    }

    // Mixture thermal conductivity [W/m K]
    double kappa(const double T, std::vector<double>& Yd) const {
        std::vector<double> Xd = calcXd(Yd);
        double res = 0.0;
        for (int i = 0; i < numComponents_; i++) {
            res += Xd[i] * components_[i]->kappa(T);
        }
        return res;
    }

    // Mixture viscousity [Pa s]
    double mu(const double T, std::vector<double>& Yd) const {
        double res = 0.0;
        std::vector<double> Xd = calcXd(Yd);
        for (int i = 0; i < numComponents_; i++) {
            res += Xd[i] * std::log(components_[i]->mu(T));
        }
        return std::exp(res);
    }


private:
    int numComponents_;
    std::vector<std::shared_ptr<Liquid>> components_;
    std::vector<std::string> names_;

    std::vector<double> calcXd(const std::vector<double>& Yd) const;

    double TcMin_;
    double TbMin_;
};

#endif  // CTF_MIXTURE_H_
