#include "Toluene.h"

Toluene::Toluene(const double p0)
    : Liquid(),
      p0_(p0),
      Wf_(92.141),
      Tc_(591.79),
      name_("C7H8")
{
    std::cout << "Toluene Constructed" << std::endl;
    // Compute boiling temperature using bisection method
    double Thi = Tc_;
    double Tlo = 0.1*Tc_;
    double Tmid = 0.5*(Thi + Tlo);
    while ((Thi - Tlo) > 1e-04) {
        if ((pv(Tmid) - p0_) > 0) Thi = Tmid;
        else Tlo = Tmid;
        Tmid = 0.5*(Thi + Tlo);
    }
    Tb_ = Tlo;  // slightly lower
    std::cout << "Toluene boiling: " << Tb_ << std::endl;
}
