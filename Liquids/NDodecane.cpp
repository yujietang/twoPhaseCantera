#include "NDodecane.h"

NDodecane::NDodecane(const double p0)
    : Liquid(),
      p0_(p0),
      Wf_(170.338),
      Tc_(658.0),
      name_("NC12H26")
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
    std::cout << "NDodecane boiling: " << Tb_ << std::endl;
}
