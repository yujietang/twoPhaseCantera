#include "NDecane.h"

NDecane::NDecane(const double p0)
    : Liquid(),
      p0_(p0),
      Wf_(142.285),
      Tc_(617.70),
      name_("NXC10H22")
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
    std::cout << "NDecane boiling: " << Tb_ << std::endl;
}
