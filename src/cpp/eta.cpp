// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "eta.hpp"
#include <limits>

namespace polylogarithm {

namespace {

// Table[PolyLog[n,-1], {n,1,54}]
const double NEG_ETA[54] = {
    -0.69314718055994531, -0.82246703342411322, -0.90154267736969571,
    -0.94703282949724592, -0.97211977044690931, -0.98555109129743510,
    -0.99259381992283028, -0.99623300185264790, -0.99809429754160533,
    -0.99903950759827157, -0.99951714349806075, -0.99975768514385819,
    -0.99987854276326512, -0.99993917034597972, -0.99996955121309924,
    -0.99998476421490611, -0.99999237829204101, -0.99999618786961011,
    -0.99999809350817168, -0.99999904661158152, -0.99999952325821554,
    -0.99999976161323082, -0.99999988080131844, -0.99999994039889239,
    -0.99999997019885696, -0.99999998509923200, -0.99999999254955048,
    -0.99999999627475340, -0.99999999813736942, -0.99999999906868228,
    -0.9999999995343403 , -0.9999999997671699 , -0.9999999998835849 ,
    -0.9999999999417924 , -0.9999999999708962 , -0.9999999999854481 ,
    -0.9999999999927240 , -0.9999999999963620 , -0.9999999999981810 ,
    -0.9999999999990905 , -0.9999999999995453 , -0.9999999999997726 ,
    -0.9999999999998863 , -0.9999999999999432 , -0.9999999999999716 ,
    -0.9999999999999858 , -0.9999999999999929 , -0.9999999999999964 ,
    -0.9999999999999982 , -0.9999999999999991 , -0.9999999999999996 ,
    -0.9999999999999998 , -0.9999999999999999 , -0.9999999999999999
};

constexpr bool is_even(int64_t n) noexcept { return n % 2 == 0; }

} // anonymous namespace

/// negative Dirichlet eta function
double neg_eta(int64_t n) noexcept
{
   if (n < 1) {
      return std::numeric_limits<double>::quiet_NaN();
   } else if (n <= sizeof(NEG_ETA)/sizeof(NEG_ETA[0])) {
      return NEG_ETA[n - 1];
   } else {
      return -1.0;
   }
}

} // namespace polylogarithm
