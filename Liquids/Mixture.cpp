#include "Mixture.h"

Mixture::Mixture(std::vector<std::shared_ptr<Liquid>>& components)
{
    components_ = components;
    numComponents_ = components.size();
    for (auto liq : components) {
        names_.push_back(liq->name());
    }

    TcMin_ = components[0]->Tc();
    TbMin_ = components[0]->Tb();
    for (int i = 1; i < numComponents_; i++) {
        TcMin_ = std::min(TcMin_, components[i]->Tc());
        TbMin_ = std::max(TbMin_, components[i]->Tb());
    }
}

void Mixture::updateProp(Eigen::VectorXd& Td, Eigen::VectorXd& rhod,
                         Eigen::VectorXd& kappad, Eigen::VectorXd& cL,
                         Eigen::VectorXd& mud, std::vector<Eigen::VectorXd>& Yd) const
{
    for (int j = 0; j < Td.size(); j++) {
        std::vector<double> ydj(numComponents_);
        for (int i = 0; i < numComponents_; i++) {
            ydj[i] = Yd[i](j);
        }

        Td(j) = std::max(0.2*this->TcMin_, std::min(this->TbMin_, Td(j)));
        rhod(j) = this->rho(Td(j), ydj);
        kappad(j) = this->kappa(Td(j), ydj);
        cL(j) = this->cp(Td(j), ydj);
        mud(j) = this->mu(Td(j), ydj);
    }
}

std::vector<double> Mixture::calcXd(const std::vector<double>& Yd) const
{
    std::vector<double> Xd(numComponents_, 0.0);
    double sumxd = 0.0;
    for (int n = 0; n < numComponents_; n++) {
        Xd[n] = Yd[n] / components_[n]->Wf();
        sumxd += Xd[n];
    }
    for (int n = 0; n < numComponents_; n++) {
        Xd[n] /= sumxd;
    }
    return std::move(Xd);
}
