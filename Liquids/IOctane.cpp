#include "IOctane.h"

IOctane::IOctane(const double p0)
    : Liquid(),
      p0_(p0),
      Wf_(114.23),
      Tc_(543.96),
      name_("IC8H18")
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
    std::cout << "IOctane boiling: " << Tb_ << std::endl;
}
