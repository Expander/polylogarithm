// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "zeta.hpp"
#include <cmath>
#include <limits>

namespace polylogarithm {

namespace {

// zeta(n) for n = 2,...,33
const double ZETAS_POS[32] = {
    1.6449340668482264, 1.2020569031595943, 1.0823232337111382,
    1.0369277551433699, 1.0173430619844491, 1.0083492773819228,
    1.0040773561979443, 1.0020083928260822, 1.0009945751278181,
    1.0004941886041195, 1.0002460865533080, 1.0001227133475785,
    1.0000612481350587, 1.0000305882363070, 1.0000152822594087,
    1.0000076371976379, 1.0000038172932650, 1.0000019082127166,
    1.0000009539620339, 1.0000004769329868, 1.0000002384505027,
    1.0000001192199260, 1.0000000596081891, 1.0000000298035035,
    1.0000000149015548, 1.0000000074507118, 1.0000000037253340,
    1.0000000018626597, 1.0000000009313274, 1.0000000004656629,
    1.0000000002328312, 1.0000000001164155
};

// zeta(1 - 2n) for n = 1,...,130, i.e. zeta(-1), zeta(-3), zeta(-5), ...
const double ZETAS_NEG[130] = {
   -8.3333333333333333e-02,  8.3333333333333333e-03, -3.9682539682539683e-03,
    4.1666666666666667e-03, -7.5757575757575758e-03,  2.1092796092796093e-02,
   -8.3333333333333333e-02,  4.4325980392156863e-01, -3.0539543302701197e000,
    2.6456212121212121e001, -2.8146014492753623e002,  3.6075105463980464e003,
   -5.4827583333333333e004,  9.7493682385057471e005, -2.0052695796688079e007,
    4.7238486772162990e008, -1.2635724795916667e010,  3.8087931125245369e011,
   -1.2850850499305083e013,  4.8241448354850170e014, -2.0040310656516253e016,
    9.1677436031953308e017, -4.5979888343656503e019,  2.5180471921451096e021,
   -1.5001733492153929e023,  9.6899578874635941e024, -6.7645882379292821e026,
    5.0890659468662290e028, -4.1147288792557979e030,  3.5666582095375556e032,
   -3.3066089876577577e034,  3.2715634236478716e036, -3.4473782558278054e038,
    3.8614279832705259e040, -4.5892974432454332e042,  5.7775386342770432e044,
   -7.6919858759507135e046,  1.0813635449971655e049, -1.6029364522008965e051,
    2.5019479041560463e053, -4.1067052335810212e055,  7.0798774408494581e057,
   -1.2804546887939509e060,  2.4267340392333524e062, -4.8143218874045769e064,
    9.9875574175727531e066, -2.1645634868435186e069,  4.8962327039620553e071,
   -1.1549023923963520e074,  2.8382249570693707e076, -7.2612008803606716e078,
    1.9323514233419812e081, -5.3450160425288624e083,  1.5356028846422423e086,
   -4.5789872682265798e088,  1.4162025212194809e091, -4.5400652296092655e093,
    1.5076656758807860e096, -5.1830949148264564e098,  1.8435647427256529e101,
   -6.7805554753090959e103,  2.5773326702754605e106, -1.0119112875704598e109,
    4.1016346161542292e111, -1.7155244534032019e114,  7.4003425705269094e116,
   -3.2909225357054443e119,  1.5079831534164771e122, -7.1169879188254549e124,
    3.4580429141577772e127, -1.7290907606676748e130,  8.8936991695032969e132,
   -4.7038470619636015e135,  2.5571938231060206e138, -1.4284067500443528e141,
    8.1952152218313783e143, -4.8276485422727372e146,  2.9189612374770324e149,
   -1.8108932162568904e152,  1.1523577220021169e155, -7.5192311951981770e157,
    5.0294016576411050e160, -3.4473420444477677e163,  2.4207458645868515e166,
   -1.7409465920377677e169,  1.2819489863482243e172, -9.6624121108560918e174,
    7.4526910304300896e177, -5.8808393311674371e180,  4.7462718654907615e183,
   -3.9169132594772825e186,  3.3045071443226032e189, -2.8492890550994583e192,
    2.5103329345077587e195, -2.2593901995475253e198,  2.0769138004287608e201,
   -1.9494732174927259e204,  1.8680731471265914e207, -1.8270752662814577e210,
    1.8235386322595677e213, -1.8568690810125945e216,  1.9287189851195602e219,
   -2.0431170460286448e222,  2.2068411644527846e225, -2.4300821796490274e228,
    2.7274887879083470e231, -3.1197421573755085e234,  3.6358938724282600e237,
   -4.3168300030760883e240,  5.2204244879387200e243, -6.4292606949769305e246,
    8.0623033870130844e249, -1.0292714737903011e253,  1.3375329699780524e256,
   -1.7689480902797380e259,  2.3806479018092397e262, -3.2597127947194185e265,
    4.5404962371601213e268, -6.4328575193147851e271,  9.2687048675749311e274,
   -1.3579619500285181e278,  2.0227839736049322e281, -3.0629906992208336e284,
    4.7143085300742652e287, -7.3741045871355758e290,  1.1720962767050827e294,
   -1.8928866644685657e297,  3.1055517596048927e300, -5.1754977470366798e303,
    8.7601563446229215e306
};

constexpr bool is_even(int64_t n) noexcept { return n % 2 == 0; }

} // anonymous namespace

/// Riemann zeta function for integer arguments
double zeta(int64_t n) noexcept
{
    if (n < 0) {
        if (is_even(n)) {
           return 0.0;
        } else if (-(1 + n)/2 < sizeof(ZETAS_NEG)/sizeof(ZETAS_NEG[0])) {
            return ZETAS_NEG[-(1 + n)/2];
        } else if (is_even((1 - n)/2)) {
           return std::numeric_limits<double>::infinity();
        } else {
           return -std::numeric_limits<double>::infinity();
        }
    } else if (n == 0) {
       return -0.5;
    } else if (n == 1) {
       return std::numeric_limits<double>::infinity();
    } else if ((n - 2) < sizeof(ZETAS_POS)/sizeof(ZETAS_POS[0])) {
       return ZETAS_POS[n - 2];
    }

    return 1.0/(1.0 - std::pow(0.5, n));
}

} // namespace polylogarithm
