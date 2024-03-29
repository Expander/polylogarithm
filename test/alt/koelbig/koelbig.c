#include <math.h>


static double sign(double a, double b) {
   return b >= 0 ? fabs(a) : -fabs(a);
}


/**
 * @brief Clausen function \f$\operatorname{Cl}_2(x)\f$
 * @param x real argument
 * @return \f$\operatorname{Cl}_2(x)\f$
 * @author K.S. Kölbig
 * @note Implementation translated by Alexander Voigt from CERNLIB DCLAUS function C326
 *
 * Journal of Computational and Applied Mathematics 64 (1995) 295-297.
 *
 * Implemented as a truncated series expansion in terms of Chebyshev
 * polynomials, see [Yudell L. Luke: Mathematical functions and their
 * approximations, Academic Press Inc., New York 1975, p.67].
 */
double clausen_2_koelbig(double X)
{
   const double R1 = 1, HF = R1/2;
   const double PI = 3.14159265358979324;
   const double PI2 = 2*PI, PIH = PI/2, RPIH = 2/PI;
   const double A[9] = {
      0.0279528319735756613,
      0.0001763088743898116,
      0.0000012662741461157,
      0.0000000117171818134,
      0.0000000001230064129,
      0.0000000000013952729,
      0.0000000000000166908,
      0.0000000000000002076,
      0.0000000000000000027
   };
   const double B[14] = {
       0.639097088857265341,
      -0.054980569301851716,
      -0.000961261945950606,
      -0.000032054686822550,
      -0.000001329461695426,
      -0.000000062093601824,
      -0.000000003129600656,
      -0.000000000166351954,
      -0.000000000009196527,
      -0.000000000000524004,
      -0.000000000000030580,
      -0.000000000000001820,
      -0.000000000000000110,
      -0.000000000000000007,
   };

   double ALFA, B0, B1, B2, H, U, S, V;

   V = fmod(fabs(X), PI2);
   S = sign(R1, X);

   if (V > PI) {
      V = PI2 - V;
      S = -S;
   }

   if (V == 0 || V == PI) {
      H = 0;
   } else if (V < PIH) {
      U = RPIH*V;
      H = 2*U*U - 1;
      ALFA = H + H;
      B1 = 0;
      B2 = 0;
      for (int i = 8; i >= 0; i--) {
         B0 = A[i] + ALFA*B1 - B2;
         B2 = B1;
         B1 = B0;
      }
      H = V*(1 - log(V) + HF*V*V*(B0 - H*B2));
   } else {
      U = RPIH*V - 2;
      H = 2*U*U - 1;
      ALFA = H + H;
      B1 = 0;
      B2 = 0;
      for (int i = 13; i >= 0; i--) {
         B0 = B[i] + ALFA*B1 - B2;
         B2 = B1;
         B1 = B0;
      }
      H = (PI - V)*(B0 - H*B2);
   }

   return S*H;
}


/**
 * @brief Clausen function \f$\operatorname{Cl}_2(\theta) = \operatorname{Im}(\operatorname{Li}_2(e^{i\theta}))\f$ with long double precision
 * @param x real angle
 * @return \f$\operatorname{Cl}_2(\theta)\f$
 * @author K.S. Kölbig
 * @note Implementation translated from CERNLIB DCLAUS function C326
 * and extended to long double precision by Alexander Voigt.
 *
 * Journal of Computational and Applied Mathematics 64 (1995) 295-297.
 */
long double clausen_2l_koelbig(long double x)
{
   const long double PI = 3.14159265358979323846264338327950288L;
   const long double PI2 = 2*PI, PIH = PI/2, RPIH = 2/PI;
   const long double A[19] = {
      0.0279528319735756613494585924765551791L,
      0.0001763088743898115653057636473920103L,
      0.0000012662741461156530021975187159184L,
      0.0000000117171818134392379295428166866L,
      0.0000000001230064128833746922855709386L,
      0.0000000000013952728970012911958374309L,
      0.0000000000000166907761628567345146740L,
      0.0000000000000002076091315145432983502L,
      0.0000000000000000026609198306058056092L,
      0.0000000000000000000349249563561378275L,
      0.0000000000000000000004673313082962865L,
      0.0000000000000000000000063542322337428L,
      0.0000000000000000000000000875698871820L,
      0.0000000000000000000000000012208003299L,
      0.0000000000000000000000000000171890569L,
      0.0000000000000000000000000000002441331L,
      0.0000000000000000000000000000000034940L,
      0.0000000000000000000000000000000000503L,
      0.0000000000000000000000000000000000007L
   };
   const long double B[30] = {
      0.6390970888572653413071869135953864197L,
     -0.0549805693018517156397035696498958507L,
     -0.0009612619459506064293859076874070709L,
     -0.0000320546868225504765586825318112711L,
     -0.0000013294616954255450141343828695514L,
     -0.0000000620936018243975194590942773212L,
     -0.0000000031296006563911126723262365339L,
     -0.0000000001663519538192669775933926077L,
     -0.0000000000091965272507194254496027281L,
     -0.0000000000005240037738758450093649037L,
     -0.0000000000000305803841873659454134183L,
     -0.0000000000000018196918249487950988000L,
     -0.0000000000000001100398263196261522324L,
     -0.0000000000000000067451775715424687730L,
     -0.0000000000000000004182784651572477035L,
     -0.0000000000000000000261987180876106127L,
     -0.0000000000000000000016553211620337322L,
     -0.0000000000000000000001053943580037580L,
     -0.0000000000000000000000067562822456482L,
     -0.0000000000000000000000004357492971838L,
     -0.0000000000000000000000000282576222954L,
     -0.0000000000000000000000000018415141361L,
     -0.0000000000000000000000000001205472536L,
     -0.0000000000000000000000000000079233873L,
     -0.0000000000000000000000000000005227403L,
     -0.0000000000000000000000000000000346060L,
     -0.0000000000000000000000000000000022982L,
     -0.0000000000000000000000000000000001531L,
     -0.0000000000000000000000000000000000102L,
     -0.0000000000000000000000000000000000007L
   };

   long double sgn = 1;

   if (x < 0) {
      x = -x;
      sgn = -1;
   }

   if (x >= PI2) {
      x = fmodl(x, PI2);
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
      const long double u = RPIH*x;
      h = 2*u*u - 1;
      const long double alfa = h + h;
      long double b0 = 0, b1 = 0, b2 = 0;
      for (int i = 18; i >= 0; i--) {
         b0 = A[i] + alfa*b1 - b2;
         b2 = b1;
         b1 = b0;
      }
      h = x*(1 - logl(x) + x*x*(b0 - h*b2)/2);
   } else {
      const long double u = RPIH*x - 2;
      h = 2*u*u - 1;
      const long double alfa = h + h;
      long double b0 = 0, b1 = 0, b2 = 0;
      for (int i = 29; i >= 0; i--) {
         b0 = B[i] + alfa*b1 - b2;
         b2 = b1;
         b1 = b0;
      }
      h = (PI - x)*(b0 - h*b2);
   }

   return sgn*h;
}


/**
 * @brief Real dilogarithm \f$\operatorname{Li}_2(x)\f$
 * @param x real argument
 * @return \f$\operatorname{Li}_2(x)\f$
 * @author K.S. Kölbig
 * @note Implementation translated by R.Brun from CERNLIB DILOG function C332
 *
 * Implemented as a truncated series expansion in terms of Chebyshev
 * polynomials, see [Yudell L. Luke: Mathematical functions and their
 * approximations, Academic Press Inc., New York 1975, p.67].
 */
double koelbig_dilog(double x)
{
   const double PI  = 3.141592653589793;
   const double HF  = 0.5;
   const double PI2 = PI*PI;
   const double PI3 = PI2/3;
   const double PI6 = PI2/6;
   const double PI12 = PI2/12;
   const double C[20] = {  0.42996693560813697, 0.40975987533077106,
     -0.01858843665014592, 0.00145751084062268,-0.00014304184442340,
      0.00001588415541880,-0.00000190784959387, 0.00000024195180854,
     -0.00000003193341274, 0.00000000434545063,-0.00000000060578480,
      0.00000000008612098,-0.00000000001244332, 0.00000000000182256,
     -0.00000000000027007, 0.00000000000004042,-0.00000000000000610,
      0.00000000000000093,-0.00000000000000014, 0.00000000000000002};

   double T, H, Y, S, A, ALFA, B1, B2, B0;

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
       for (int i = 19; i >= 0; i--) {
          B0 = C[i] + ALFA*B1-B2;
          B2 = B1;
          B1 = B0;
       }
       H = -(S*(B0-H*B2)+A);
    }
    return H;
}


/**
 * @brief Real dilogarithm \f$\operatorname{Li}_2(z)\f$ with long double precision
 * @param x real argument
 * @return \f$\operatorname{Li}_2(z)\f$
 * @author K.S. Kölbig
 * @note Implementation based on translation by R.Brun from CERNLIB
 *    DILOG function C332, extended by Alexander Voigt to quadruple
 *    precision
 *
 * Implemented as a truncated series expansion in terms of Chebyshev
 * polynomials, see [Yudell L. Luke: Mathematical functions and their
 * approximations, Academic Press Inc., New York 1975, p.67].
 */
long double koelbig_dilogl(long double x)
{
   const long double PI  = 3.14159265358979323846264338327950288L;
   const long double HF  = 0.5L;
   const long double PI2 = PI*PI;
   const long double PI3 = PI2/3;
   const long double PI6 = PI2/6;
   const long double PI12 = PI2/12;
   const long double C[] = {
      0.4299669356081369720370336786993879912L,
      0.4097598753307710584682637109252552781L,
     -0.0185884366501459196476416491402122676L,
      0.0014575108406226785536739284164594927L,
     -0.0001430418444234004877488301200908765L,
      0.0000158841554187955323619055047167740L,
     -0.0000019078495938658272271963211420884L,
      0.0000002419518085416474949946146434290L,
     -0.0000000319334127425178346049601414286L,
      0.0000000043454506267691229879571784782L,
     -0.0000000006057848011840744442970533091L,
      0.0000000000861209779935949824428368452L,
     -0.0000000000124433165993886798964242137L,
      0.0000000000018225569623573633006554774L,
     -0.0000000000002700676604911465180857223L,
      0.0000000000000404220926315266464883286L,
     -0.0000000000000061032514526918795037783L,
      0.0000000000000009286297533019575861303L,
     -0.0000000000000001422602085511244683975L,
      0.0000000000000000219263171815395735398L,
     -0.0000000000000000033979732421589786340L,
#if LDBL_DIG > 18
      0.0000000000000000005291954244833147146L,
     -0.0000000000000000000827858081427899765L,
      0.0000000000000000000130037173454556037L,
     -0.0000000000000000000020502222425528249L,
      0.0000000000000000000003243578549148930L,
     -0.0000000000000000000000514779990334321L,
      0.0000000000000000000000081938774771716L,
     -0.0000000000000000000000013077835405713L,
      0.0000000000000000000000002092562930580L,
     -0.0000000000000000000000000335616615054L,
      0.0000000000000000000000000053946577714L,
     -0.0000000000000000000000000008689193209L,
      0.0000000000000000000000000001402281687L,
     -0.0000000000000000000000000000226715578L,
      0.0000000000000000000000000000036717417L,
     -0.0000000000000000000000000000005956152L,
      0.0000000000000000000000000000000967662L,
     -0.0000000000000000000000000000000157439L,
      0.0000000000000000000000000000000025650L,
     -0.0000000000000000000000000000000004185L,
      0.0000000000000000000000000000000000683L,
     -0.0000000000000000000000000000000000112L,
      0.0000000000000000000000000000000000018L,
     -0.0000000000000000000000000000000000003L
#endif
   };

   long double T, H, Y, S, A, ALFA, B1, B2, B0;

   if (x == 1) {
       H = PI6;
   } else if (x == -1) {
       H = -PI12;
   } else {
       T = -x;
       if (T <= -2) {
           Y = -1/(1+T);
           S = 1;
           B1= logl(-T);
           B2= logl(1+1/T);
           A = -PI3+HF*(B1*B1-B2*B2);
       } else if (T < -1) {
           Y = -1-T;
           S = -1;
           A = logl(-T);
           A = -PI6+A*(A+logl(1+1/T));
       } else if (T <= -0.5L) {
           Y = -(1+T)/T;
           S = 1;
           A = logl(-T);
           A = -PI6+A*(-HF*A+logl(1+T));
       } else if (T < 0) {
           Y = -T/(1+T);
           S = -1;
           B1= logl(1+T);
           A = HF*B1*B1;
       } else if (T <= 1) {
           Y = T;
           S = 1;
           A = 0;
       } else {
           Y = 1/T;
           S = -1;
           B1= logl(T);
           A = PI6+HF*B1*B1;
       }
       H    = Y+Y-1;
       ALFA = H+H;
       B1   = 0;
       B2   = 0;
       for (int i = sizeof(C)/sizeof(C[0]) - 1; i >= 0; i--) {
          B0 = C[i] + ALFA*B1-B2;
          B2 = B1;
          B1 = B0;
       }
       H = -(S*(B0-H*B2)+A);
    }
    return H;
}
