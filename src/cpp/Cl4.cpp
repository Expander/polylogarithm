// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "Cl4.hpp"
#include <cmath>

namespace polylogarithm {

/**
 * @brief Clausen function \f$\operatorname{Cl}_4(\theta) = \operatorname{Im}(\operatorname{Li}_4(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\operatorname{Cl}_4(\theta)\f$
 * @author Alexander Voigt
 * @note Implemented as rational function approximation.
 */
double Cl4(double x) noexcept
{
   const double PI = 3.14159265358979324;
   const double PI2 = 2*PI, PIH = PI/2, PI28 = PI*PI/8;
   const double zeta3 = 1.2020569031595943;
   double sgn = 1;

   if (x < 0) {
      x = -x;
      sgn = -1;
   }

   if (x >= PI2) {
      x = std::fmod(x, PI2);
   }

   if (x > PI) {
      const double p0 = 6.28125;
      const double p1 = 0.0019353071795864769253;
      x = (p0 - x) + p1;
      sgn = -sgn;
   }

   if (x == 0 || x == PI) {
      return 0;
   }

   double h = 0;

   if (x < PIH) {
      const double P[] = {
         -3.0555555555555556e-01,  6.0521392328447206e-03,
         -1.9587493942041528e-05, -3.1137343767030358e-08
      };
      const double Q[] = {
         1.0000000000000000e+00, -2.2079728398400851e-02,
         1.0887447112236682e-04, -6.1847621370547954e-08
      };
      const double y = x*x;
      const double y2 = y*y;
      const double p = P[0] + y * P[1] + y2 * (P[2] + y * P[3]);
      const double q = Q[0] + y * Q[1] + y2 * (Q[2] + y * Q[3]);
      h = x*(zeta3 + y*(p/q + 1./6*std::log(x)));
   } else {
      const double P[] = {
         7.6223911686491336e-01, -2.4339587368267260e-01,
         2.8715364937979943e-02, -1.5368612510964667e-03,
         3.6261044225761673e-05, -2.8557977333851308e-07
      };
      const double Q[] = {
         1.0000000000000000e+00, -1.7465715261403233e-01,
         9.5439417991615653e-03, -1.7325070821666274e-04,
         5.9283675098376635e-07,  9.4127575773361230e-10
      };
      const double y = PI - x;
      const double z = y*y - PI28;
      const double z2 = z*z;
      const double z4 = z2*z2;
      const double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]) +
         z4 * (P[4] + z * P[5]);
      const double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]) +
         z4 * (Q[4] + z * Q[5]);
      h = y*p/q;
   }

   return sgn*h;
}

/**
 * @brief Clausen function \f$\operatorname{Cl}_4(\theta) = \operatorname{Im}(\operatorname{Li}_4(e^{i\theta}))\f$ with long double precision
 * @param x real angle
 * @return \f$\operatorname{Cl}_4(\theta)\f$
 * @author Alexander Voigt
 * @note Implemented as rational function approximation.
 */
long double Cl4(long double x) noexcept
{
   const long double PI = 3.14159265358979323846264338327950288L;
   const long double PI2 = 2*PI, PIH = PI/2, PI28 = PI*PI/8;
   const long double zeta3 = 1.2020569031595942853997381615114499908L;
   long double sgn = 1;

   if (x < 0) {
      x = -x;
      sgn = -1;
   }

   if (x >= PI2) {
      x = std::fmod(x, PI2);
   }

   if (x > PI) {
      const long double p0 = 6.28125L;
      const long double p1 = 0.0019353071795864769252867665590057683943L;
      x = (p0 - x) + p1;
      sgn = -sgn;
   }

   if (x == 0 || x == PI) {
      return 0;
   }

   long double h = 0;

   if (x < PIH) {
      const long double P[] = {
        -3.064148293902562247016704418203843255123e-01L,
         2.505569295848759195239609343375835160555e-02L,
        -8.092489304915887233846804452184668671559e-04L,
         1.306550392088097227871332697836673781161e-05L,
        -1.091861195748019042041801413413722149520e-07L,
         4.290976258123337917062729924440640694818e-10L,
        -4.827057875035942039653217415564619613330e-13L,
        -7.817994926896706571532523150572289005814e-16L,
         1.005354194675259498587105145642652708124e-18L
      };
      const long double Q[] = {
         1.000000000000000000000000000000000000000e+00L,
        -8.405033099067044445000733863730421988529e-02L,
         2.827113467192315016976789998041672862170e-03L,
        -4.865778884542162517878400745754972885997e-05L,
         4.544079814149520679837030712405877787979e-07L,
        -2.246626075063535176123041843073598887706e-09L,
         5.278103491358313047352307496085513128432e-12L,
        -4.503548404587830179829476685164670510929e-15L,
         5.683277279497885707083902179699889410498e-19L
      };
      const long double y = x*x;
      const long double z = y - PI28;
      const long double z2 = z*z;
      const long double z4 = z2*z2;
      const long double z8 = z4*z4;
      const long double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]) +
         z4 * (P[4] + z * P[5] + z2 * (P[6] + z * P[7])) + z8 * P[8];
      const long double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]) +
         z4 * (Q[4] + z * Q[5] + z2 * (Q[6] + z * Q[7])) + z8 * Q[8];

      h = x*(zeta3 + y*(p/q + 1.0L/6*std::log(x)));
   } else {
      const long double P[] = {
         7.622391168649133627869594744879702180493e-01L,
        -1.010231891148264580184081150011868846220e+00L,
         5.593311843617106683456367723756675040735e-01L,
        -1.750142535337254911463918175625637736546e-01L,
         3.444639642923636356477676643239523643873e-02L,
        -4.467089630163625626375503936549208950240e-03L,
         3.886778570388940085712534905298979385012e-04L,
        -2.267243475184570369108265113505191140730e-05L,
         8.702335438867849048873191203867941600263e-07L,
        -2.109129915929703527128733018891419256495e-08L,
         2.988467123764362417450870076317332442258e-10L,
        -2.133050389878217363160872116919842900528e-12L,
         5.379728143941033611125864127082270782090e-15L
      };
      const long double Q[] = {
         1.000000000000000000000000000000000000000e+00L,
        -1.180687938097183154430487264715729257263e+00L,
         5.601395009696004947928229952758353602894e-01L,
        -1.452335710129900283499882062220239913044e-01L,
         2.262139281960755856101092904992480836291e-02L,
        -2.192073007443907755835491584063761368816e-03L,
         1.329216132182395210305478468336390325445e-04L,
        -4.948002266894020612513545634482271857510e-06L,
         1.075239424448883617918347989934070034554e-07L,
        -1.231942925472515536287130721129266117605e-09L,
         5.987816416306400663411961902200439505519e-12L,
        -6.217050327376472635914508068068185584426e-15L,
        -3.160114789318241107017689802541627876077e-18L
      };
      const long double y = PI - x;
      const long double z = y*y - PI28;
      const long double z2 = z*z;
      const long double z4 = z2*z2;
      const long double z8 = z4*z4;
      const long double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]) +
         z4 * (P[4] + z * P[5] + z2 * (P[6] + z * P[7])) +
         z8 * (P[8] + z * P[9] + z2 * (P[10] + z * P[11]) + z4 * P[12]);
      const long double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]) +
         z4 * (Q[4] + z * Q[5] + z2 * (Q[6] + z * Q[7])) +
         z8 * (Q[8] + z * Q[9] + z2 * (Q[10] + z * Q[11]) + z4 * Q[12]);

      h = y*p/q;
   }

   return sgn*h;
}

} // namespace polylogarithm
