/* Device-side stand-ins for flang runtime / libm calls that amdflang
 * emits inside OpenMP target regions. Host already has these in
 * libflang_rt.runtime / libm; --allow-multiple-definition is set. */
#include <math.h>

#pragma omp begin declare target

double _FortranAModReal8(double a, double p)
{
  return a - trunc(a / p) * p;
}

float _FortranAModReal4(float a, float p)
{
  return a - truncf(a / p) * p;
}

double _FortranAModuloReal8(double a, double p)
{
  return a - floor(a / p) * p;
}

float _FortranAModuloReal4(float a, float p)
{
  return a - floorf(a / p) * p;
}

double lgamma(double xx)
{
  const double cof0 = 76.18009172947146;
  const double cof1 = -86.50532032941677;
  const double cof2 = 24.01409824083091;
  const double cof3 = -1.231739572450155;
  const double cof4 = 1.208650973866179e-3;
  const double cof5 = -5.395239384953e-6;
  double x, y, tmp, ser;
  x = xx;
  y = x;
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  ser = 1.000000000190015;
  y += 1.0; ser += cof0 / y;
  y += 1.0; ser += cof1 / y;
  y += 1.0; ser += cof2 / y;
  y += 1.0; ser += cof3 / y;
  y += 1.0; ser += cof4 / y;
  y += 1.0; ser += cof5 / y;
  return -tmp + log(2.5066282746310005 * ser / x);
}

float lgammaf(float x)
{
  return (float)lgamma((double)x);
}

#pragma omp end declare target

