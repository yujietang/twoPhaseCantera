#include "NC7H16.h"

NC7H16::NC7H16(const double p0)
    : Liquid(),
      p0_(p0),
      Wf_(100.204),
      Tc_(540.2),
      name_("NC7H16")
{
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
    std::cout << "N-heptane boiling: " << Tb_ << std::endl;
}
