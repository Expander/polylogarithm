#pragma once

#ifdef __cplusplus
extern "C" {
#endif

double algorithm_327(double x);

double algorithm_490(double x);
double algorithm_490_2(double x);

double babar_dilog(double x);

double cephes_dilog(double x);

double cephes_dilog_2(double x);

double hassani_dilog(double x);

double koelbig_dilog(double x);

long double koelbig_dilogl(long double x);

void feynhiggs_dilog(long double re, long double im, long double* res_re, long double* res_im);

double morris_dilog(double x);

double pythia_dilog(const double x, const double kmax, const double xerr);

void hdecay_dilog(double re, double im, double* res_re, double* res_im);

void hollik_dilog(double re, double im, double* res_re, double* res_im);

void sherpa_dilog(double re, double im, double* res_re, double* res_im);

void spheno_dilog(double re, double im, double* res_re, double* res_im);

long double tsil_dilog_real(long double x);

void tsil_dilog_complex(long double re, long double im, long double* res_re, long double* res_im);

void tsil_trilog_complex(long double re, long double im, long double* res_re, long double* res_im);

// Clausen functions

double clausen_2_bernoulli(double x);
double clausen_2_koelbig(double x);
long double clausen_2l_koelbig(long double x);
double clausen_2_pade(double x);
double clausen_3_pade(double x);
double clausen_2_wu(double x);
double clausen_3_wu(double x);

#ifdef __cplusplus
}
#endif
