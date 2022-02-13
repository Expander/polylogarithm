// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "inv_fac.hpp"
#include <limits>

namespace polylogarithm {

namespace {

// 1/n! for integer n = 0, 1, 2, ...
const double INVERSE_FACTORIALS[178] = {
    1.0,    1.0,     0.5   , 1.6666666666666667e-001, 4.1666666666666667e-002,
    8.3333333333333333e-003, 1.3888888888888889e-003, 1.9841269841269841e-004,
    2.4801587301587302e-005, 2.7557319223985891e-006, 2.7557319223985891e-007,
    2.5052108385441719e-008, 2.0876756987868099e-009, 1.6059043836821615e-010,
    1.1470745597729725e-011, 7.6471637318198165e-013, 4.7794773323873853e-014,
    2.8114572543455208e-015, 1.5619206968586226e-016, 8.2206352466243297e-018,
    4.1103176233121649e-019, 1.9572941063391261e-020, 8.8967913924505733e-022,
    3.8681701706306840e-023, 1.6117375710961183e-024, 6.4469502843844734e-026,
    2.4795962632247975e-027, 9.1836898637955461e-029, 3.2798892370698379e-030,
    1.1309962886447717e-031, 3.7699876288159056e-033, 1.2161250415535179e-034,
    3.8003907548547436e-036, 1.1516335620771950e-037, 3.3871575355211618e-039,
    9.6775929586318910e-041, 2.6882202662866364e-042, 7.2654601791530713e-044,
    1.9119632050402819e-045, 4.9024697565135434e-047, 1.2256174391283858e-048,
    2.9893108271424045e-050, 7.1174067312914393e-052, 1.6552108677421952e-053,
    3.7618428812322618e-055, 8.3596508471828040e-057, 1.8173154015614791e-058,
    3.8666285139605939e-060, 8.0554760707512373e-062, 1.6439747083165790e-063,
    3.2879494166331581e-065, 6.4469596404571727e-067, 1.2397999308571486e-068,
    2.3392451525606577e-070, 4.3319354677049217e-072, 7.8762463049180395e-074,
    1.4064725544496499e-075, 2.4674957095607893e-077, 4.2543029475186023e-079,
    7.2106829618959360e-081, 1.2017804936493227e-082, 1.9701319568021683e-084,
    3.1776321883905941e-086, 5.0438606164930064e-088, 7.8810322132703225e-090,
    1.2124664943492804e-091, 1.8370704459837582e-093, 2.7418961880354600e-095,
    4.0322002765227352e-097, 5.8437685166996163e-099, 8.3482407381423090e-101,
    1.1758085546679308e-102, 1.6330674370387928e-104, 2.2370786808750587e-106,
    3.0230792984798090e-108, 4.0307723979730787e-110, 5.3036478920698404e-112,
    6.8878544052855070e-114, 8.8305825708788551e-116, 1.1177952621365639e-117,
    1.3972440776707049e-119, 1.7249926884823518e-121, 2.1036496201004290e-123,
    2.5345176145788301e-125, 3.0172828744986072e-127, 3.5497445582336556e-129,
    4.1276099514344832e-131, 4.7443792545223945e-133, 5.3913400619572665e-135,
    6.0576854628733332e-137, 6.7307616254148146e-139, 7.3964413466096864e-141,
    8.0396101593583548e-143, 8.6447421068369406e-145, 9.1965341562095113e-147,
    9.6805622696942224e-149, 1.0083919030931482e-150, 1.0395792815393280e-152,
    1.0607951852442123e-154, 1.0715102881254669e-156, 1.0715102881254669e-158,
    1.0609012753717494e-160, 1.0400992895801465e-162, 1.0098051355147053e-164,
    9.7096647645644744e-167, 9.2472997757756899e-169, 8.7238677129959339e-171,
    8.1531473953233027e-173, 7.5492105512252803e-175, 6.9258812396562204e-177,
    6.2962556724147458e-179, 5.6723024075808521e-181, 5.0645557210543322e-183,
    4.4819077177471967e-185, 3.9314979980238567e-187, 3.4186939113250928e-189,
    2.9471499235561145e-191, 2.5189315585949697e-193, 2.1346877615211607e-195,
    1.7938552617824880e-197, 1.4948793848187400e-199, 1.2354375081146612e-201,
    1.0126536951759518e-203, 8.2329568713492014e-206, 6.6394813478622592e-208,
    5.3115850782898073e-210, 4.2155437129284185e-212, 3.3193257582113532e-214,
    2.5932232486026197e-216, 2.0102505803121083e-218, 1.5463466002400833e-220,
    1.1804172520916666e-222, 8.9425549400883835e-225, 6.7237255188634463e-227,
    5.0177056110921241e-229, 3.7168189711793512e-231, 2.7329551258671700e-233,
    1.9948577561074233e-235, 1.4455490986285676e-237, 1.0399633803083220e-239,
    7.4283098593451574e-242, 5.2683048647837996e-244, 3.7100738484392955e-246,
    2.5944572366708360e-248, 1.8017064143547472e-250, 1.2425561478308602e-252,
    8.5106585467867134e-255, 5.7895636372698731e-257, 3.9118673224796440e-259,
    2.6254143103890228e-261, 1.7502762069260152e-263, 1.1591233158450432e-265,
    7.6258112884542314e-268, 4.9841903846106088e-270, 3.2364872627341615e-272,
    2.0880562985381687e-274, 1.3384976272680569e-276, 8.5254625940640566e-279,
    5.3958624013063649e-281, 3.3936241517650094e-283, 2.1210150948531309e-285,
    1.3174006800330005e-287, 8.1321029631666700e-290, 4.9890202228016380e-292,
    3.0420855017083159e-294, 1.8436881828535248e-296, 1.1106555318394728e-298,
    6.6506319271824716e-301, 3.9587094804657569e-303, 2.3424316452460100e-305,
    1.3779009677917706e-307, 8.0579003964431028e-310, 4.6848258118855249e-312,
    2.7079917987777601e-314, 1.5563171257343449e-316, 8.8932407184819707e-319,
    5.0529776809556651e-321, 2.8547896502574379e-323
};

} // anonymous namespace

/// negative Dirichlet eta function
double inv_fac(int64_t n) noexcept
{
   if (n < 0) {
      return std::numeric_limits<double>::quiet_NaN();
   } else if (n < sizeof(INVERSE_FACTORIALS)/sizeof(INVERSE_FACTORIALS[0])) {
         return INVERSE_FACTORIALS[n];
      } else {
      return 0.0;
   }
}

} // namespace polylogarithm
