// ====================================================================
// This file is part of Dilogarithm.
//
// Dilogarithm is licenced under the GNU Lesser General Public License
// (GNU LGPL) version 3.
// ====================================================================

#include "dilog.hpp"
#include <cmath>
#include <limits>

namespace dilogarithm {

namespace {
   template <typename T>
   T sqr(T x) noexcept { return x*x; }

   template <typename T>
   T pow3(T x) noexcept { return x*x*x; }

   template <typename T>
   T pow4(T x) noexcept { return x*x*x*x; }

   // converts -0.0 to 0.0
   std::complex<double> clog(std::complex<double> z) noexcept {
      if (std::real(z) == 0.0) z.real(0.0);
      if (std::imag(z) == 0.0) z.imag(0.0);
      return std::log(z);
   }
} // anonymous namespace

/**
 * @brief Real polylogarithm \f$\mathrm{Li}_0(x) = x/(1-x)\f$
 * @param x real argument
 * @return \f$\mathrm{Li}_0(x)\f$
 */
double Li0(double x)
{
   return x/(1. - x);
}

/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_0(z) = z/(1-z)\f$
 * @param z complex argument
 * @return \f$\mathrm{Li}_0(z)\f$
 */
std::complex<double> Li0(const std::complex<double>& z)
{
   return z/(1. - z);
}

/**
 * @brief Real polylogarithm \f$\mathrm{Li}_1(x) = -\log(1-x)\f$
 * @param x real argument
 * @return \f$\mathrm{Li}_1(x)\f$
 */
double Li1(double x)
{
   return -std::log(1. - x);
}

/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_1(z) = -\log(1-z)\f$
 * @param z complex argument
 * @return \f$\mathrm{Li}_1(z)\f$
 */
std::complex<double> Li1(const std::complex<double>& z)
{
   return -clog(1. - z);
}

/**
 * @brief Real dilogarithm \f$\mathrm{Li}_2(x)\f$
 * @param x real argument
 * @note Implementation translated by R.Brun from CERNLIB DILOG function C332
 * @return \f$\mathrm{Li}_2(x)\f$
 */
double Li2(double x) {
   using std::log;

   const double PI = M_PI;
   const double HF  = 0.5;
   const double PI2 = PI*PI;
   const double PI3 = PI2/3;
   const double PI6 = PI2/6;
   const double PI12 = PI2/12;
   const double C[20] = {0.42996693560813697, 0.40975987533077105,
     -0.01858843665014592, 0.00145751084062268,-0.00014304184442340,
      0.00001588415541880,-0.00000190784959387, 0.00000024195180854,
     -0.00000003193341274, 0.00000000434545063,-0.00000000060578480,
      0.00000000008612098,-0.00000000001244332, 0.00000000000182256,
     -0.00000000000027007, 0.00000000000004042,-0.00000000000000610,
      0.00000000000000093,-0.00000000000000014, 0.00000000000000002};

   double T,H,Y,S,A,ALFA,B1,B2,B0;

   if (x == 1) {
       H = PI6;
   } else if (x == -1) {
       H = -PI12;
   } else {
       T = -x;
       if (T <= -2) {
           Y = -1/(1+T);
           S = 1;
           B1= log(-T);
           B2= log(1+1/T);
           A = -PI3+HF*(B1*B1-B2*B2);
       } else if (T < -1) {
           Y = -1-T;
           S = -1;
           A = log(-T);
           A = -PI6+A*(A+log(1+1/T));
       } else if (T <= -0.5) {
           Y = -(1+T)/T;
           S = 1;
           A = log(-T);
           A = -PI6+A*(-HF*A+log(1+T));
       } else if (T < 0) {
           Y = -T/(1+T);
           S = -1;
           B1= log(1+T);
           A = HF*B1*B1;
       } else if (T <= 1) {
           Y = T;
           S = 1;
           A = 0;
       } else {
           Y = 1/T;
           S = -1;
           B1= log(T);
           A = PI6+HF*B1*B1;
       }
       H    = Y+Y-1;
       ALFA = H+H;
       B1   = 0;
       B2   = 0;
       for (int i=19;i>=0;i--){
          B0 = C[i] + ALFA*B1-B2;
          B2 = B1;
          B1 = B0;
       }
       H = -(S*(B0-H*B2)+A);
    }
    return H;
}

/**
 * @brief Complex dilogarithm \f$\mathrm{Li}_2(z)\f$
 * @param z complex argument
 * @note Implementation translated from SPheno to C++
 * @return \f$\mathrm{Li}_2(z)\f$
 */
std::complex<double> Li2(const std::complex<double>& z) {
   std::complex<double> cy, cz;
   int jsgn, ipi12;
   static const int N = 20;

   // bf[1..N-1] are the even Bernoulli numbers / (2 n + 1)!
   // generated by: Table[BernoulliB[2 n]/(2 n + 1)!, {n, 1, 19}]
   const double bf[N] = {
      - 1./4.,
      + 1./36.,
      - 1./36.e2,
      + 1./21168.e1,
      - 1./108864.e2,
      + 1./52690176.e1,
      - 4.0647616451442255268059093862919666745470571274397078e-11,
      + 8.9216910204564525552179873167527488515142836130490451e-13,
      - 1.9939295860721075687236443477937897056306947496538801e-14,
      + 4.5189800296199181916504765528555932283968190144666184e-16,
      - 1.0356517612181247014483411542218656665960912381686505e-17,
      + 2.3952186210261867457402837430009803816789490019429743e-19,
      - 5.5817858743250093362830745056254199055670546676443981e-21,
      + 1.3091507554183212858123073991865923017498498387833038e-22,
      - 3.0874198024267402932422797648664624315955652561327457e-24,
      + 7.315975652702203420357905609252148591033401063690875e-26,
      - 1.7408456572340007409890551477597025453408414217542713e-27,
      + 4.1576356446138997196178996207752266734882541595115639e-29,
      - 9.9621484882846221031940067024558388498548600173944888e-31,
      + 2.3940344248961653005211679878937495629342791569329158e-32,
   };

   const double rz = std::real(z);
   const double iz = std::imag(z);
   const double az = std::sqrt(sqr(rz) + sqr(iz));

   // special cases
   if (iz == 0.) {
      if (rz <= 1.)
         return {Li2(rz), 0.};
      if (rz > 1.)
         return {Li2(rz), -M_PI*std::log(rz)};
   } else if (az < std::numeric_limits<double>::epsilon()) {
      return z;
   }

   // transformation to |z|<1, Re(z)<=0.5
   if (rz <= 0.5) {
      if (az > 1.) {
         cy = -0.5 * sqr(std::log(-z));
         cz = -std::log(1. - 1. / z);
         jsgn = -1;
         ipi12 = -2;
      } else { // (az <= 1.)
         cy = 0;
         cz = -std::log(1. - z);
         jsgn = 1;
         ipi12 = 0;
      }
   } else { // rz > 0.5
      if (az <= std::sqrt(2*rz)) {
         cz = -std::log(z);
         cy = cz * std::log(1. - z);
         jsgn = -1;
         ipi12 = 2;
      } else { // (az > sqrt(2*rz))
         cy = -0.5 * sqr(std::log(-z));
         cz = -std::log(1. - 1. / z);
         jsgn = -1;
         ipi12 = -2;
      }
   }

   // the dilogarithm
   const std::complex<double> cz2(sqr(cz));
   std::complex<double> sumC;

   for (int i1 = 2; i1 < N; i1++)
      sumC = cz2 * (sumC + bf[N + 1 - i1]);

   // lowest order terms w/ different powers
   sumC = cz + cz2 * (bf[0] + cz * (bf[1] + sumC));

   const std::complex<double> result
      = double(jsgn) * sumC + cy + ipi12 * M_PI * M_PI / 12.;

   return result;
}

/**
 * @brief Complex trilogarithm \f$\mathrm{Li}_3(z)\f$
 * @param z complex argument
 * @return \f$\mathrm{Li}_3(z)\f$
 */
std::complex<double> Li3(const std::complex<double>& z)
{
   const double PI = M_PI;
   const double PI2 = PI*PI;
   const double zeta3 = 1.202056903159594285399738161511449990764986292340498881792271555;
   static const int N = 40;

   const double bf[N] = {
      1., -3./8., 17./216., -5./576.,
      0.00012962962962962962962962962962963, 0.000081018518518518518518518518518519,
      -3.4193571608537594932152755282007e-6, -1.3286564625850340136054421768707e-6,
      8.6608717561098513479465860418241e-8, 2.5260875955320399764844209288654e-8,
      -2.1446944683640647609338850757365e-9, -5.14011062201297891533581769272e-10,
      5.2495821146008294363940888085581e-11, 1.0887754406636318375372971570425e-11,
      -1.2779396094493695305581831754072e-12, -2.3698241773087452099797778810124e-13,
      3.1043578879654622942847532704656e-14, 5.2617586299125060841318392511225e-15,
      -7.5384795499492653659925014322677e-16, -1.1862322577752285253082500951246e-16,
      1.8316979965491383382089273121282e-17, 2.7068171031837350151490734712617e-18,
      -4.4554338978296388264326309921763e-19, -6.2375484922556946503653222473984e-20,
      1.0851521534874534913136560996864e-20, 1.4491174866036081930734904966528e-21,
      -2.6466339754458990334740891186144e-22, -3.3897653488510104721925816586081e-23,
      6.4640477336033108890325309821953e-24, 7.975834489602412424209222725905e-25,
      -1.5809178790287483355921117629383e-25, -1.8861499729622868193110225398853e-26,
      3.8715536638418473303997127188831e-27, 4.4801175002345607304865389832051e-28,
      -9.493033871911836126417536760277e-29, -1.0682813809077381224018214303381e-29,
      2.3304478936103051860078519901928e-30, 2.556077572651975408056356982867e-31,
      -5.7274216061372596844727445803306e-32, -6.1347132137964235825854929689777e-33
   };

   if (z == 0.)
      return 0.;
   if (z == 1.)
      return zeta3;
   if (z == -1.)
      return -0.75*zeta3;
   if (z == 0.5) {
      const double ln2 = std::log(2.);
      const double ln23 = pow3(ln2);
      return (-2*PI2*ln2 + 4*ln23 + 21*zeta3)/24.;
   }

   std::complex<double> u, sum;

   if (std::abs(z) <= 1.) {
      u = -clog(1. - z);
   } else { // az > 1.
      u = -clog(1. - 1./z);
      sum = -pow3(clog(-z))/6. - M_PI*M_PI/6.*clog(-z);
   }

   std::complex<double> p = 1.;

   for (const double b: bf) {
      p *= u;
      sum += b*p;
   }

   return sum;
}

/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_4(z)\f$
 * @param z complex argument
 * @return \f$\mathrm{Li}_4(z)\f$
 */
std::complex<double> Li4(const std::complex<double>& z)
{
   const double PI = M_PI;
   const double PI2 = PI*PI;
   const double PI3 = PI2*PI;
   const double PI4 = PI2*PI2;
   const double zeta4 = 1.0823232337111381915160036965412;
   static const int N = 40;

   const double bf[N] = {
      1., -7./16., 0.11651234567901234567901234567901, -0.019820601851851851851851851851852,
      0.0019279320987654320987654320987654, -0.000031057098765432098765432098765432,
      -0.000015624009114857835298392473643526, 8.4851235467732066371522153835079e-7,
      2.2909616603189711445359383547004e-7, -2.1832614218526916939615352313765e-8,
      -3.8828248791720155722806620380777e-9, 5.4462921032203321182579858808232e-10,
      6.9608052106827254078772334134121e-11, -1.3375737686445215199578072203635e-11,
      -1.2784852685266571604146246361574e-12, 3.2605628580248922428788418178217e-13,
      2.3647571168618257362309504812439e-14, -7.9231351220311617024299900711372e-15,
      -4.3452915709984187250497371626475e-16, 1.9236270062535920116126875526753e-16,
      7.8124143331959546707222938968737e-18, -4.6718038448036555203176282428722e-18,
      -1.3435344329812847856260222675894e-19, 1.1356826851347343244764698375938e-19,
      2.1152756202432586847505983414192e-21, -2.7642026334746517388281729253731e-21,
      -2.7068176608240064256109059558195e-23, 6.7372044828628572143267161265626e-23,
      1.3287265456683822975818009045013e-25, -1.6443773056367826467816763114889e-24,
      8.2836058999339341109829673400349e-27, 4.0190848495069350699709315007621e-26,
      -4.5757138444848790382359734346537e-28, -9.8364109094615127758320974982117e-28,
      1.6900339556037851067729523121903e-29, 2.4104805563059808504664904164902e-29,
      -5.4266127056714182501325034058929e-31, -5.9142429588741767864337599966928e-31,
      1.6232110901087370772711176143968e-32, 1.4527595437740275946132587316158e-32
   };

   if (z == 0.)
      return 0.;
   if (z == 1.)
      return zeta4;
   if (z == -1.)
      return -7.*PI4/720.;

   std::complex<double> u, r;
   double sgn = 1;

   if (std::abs(z) <= 1.) {
      u = -clog(1. - z);
   } else { // az > 1.
      const std::complex<double> lnz  = clog(-z);
      const std::complex<double> lnz2 = sqr(lnz);
      const std::complex<double> lnz4 = pow4(lnz);
      u = -clog(1. - 1./z);
      r = 1./360.*(-7*PI4 - 30.*PI2*lnz2 - 15.*lnz4);
      sgn = -1;
   }

   std::complex<double> p = 1., sum;

   for (const double b: bf) {
      p *= u;
      sum += b*p;
   }

   return sgn*sum + r;
}

/**
 * @brief Clausen function \f$\mathrm{Cl}_2(\theta) = \mathrm{Im}(\mathrm{Li}_2(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_2(\theta)\f$
 */
double Cl2(double x)
{
   using std::exp;
   const std::complex<double> i(0.,1.);

   return std::imag(Li2(exp(i*x)));
}

} // namespace dilogarithm
