#include "TMBenzene.h"

TMBenzene::TMBenzene(const double p0)
    : Liquid(),
      p0_(p0),
      Wf_(120.19),
      Tc_(649.13),
      name_("TMBENZ")
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
    std::cout << "TMBenzene boiling: " << Tb_ << std::endl;
}
