#include <fem.hpp> // Fortran EMulation library of fable module

#include "ccp4_fftlib_fem.hpp"

namespace ccp4_fftlib_fem {

using namespace fem::major_types;

static const float ccp4_fftlib_twopi = 6.283185f;

void
srfp(
  int const& pts,
  int const& pmax,
  int const& twogrp,
  arr_ref<int> factor,
  arr_ref<int> sym,
  int& psym,
  arr_ref<int> unsym,
  bool& error)
{
  factor(dimension(15));
  sym(dimension(15));
  unsym(dimension(15));
  int nest = fem::int0;
  int n = fem::int0;
  int f = fem::int0;
  int p = fem::int0;
  int q = fem::int0;
  int j = fem::int0;
  arr_1d<14, int> pp;
  arr_1d<7, int> qq;
  int r = fem::int0;
  int jj = fem::int0;
  int ptwo = fem::int0;
  nest = 14;
  n = pts;
  psym = 1;
  f = 2;
  p = 0;
  q = 0;
  statement_10:
  if (n <= 1) {
    goto statement_60;
  }
  else {
    FEM_DO(j, f, pmax) {
      if (n == (n / j) * j) {
        goto statement_30;
      }
    }
    goto statement_40;
    statement_30:
    if (2 * p + q >= nest) {
      goto statement_50;
    }
    else {
      f = j;
      n = n / f;
      if (n == (n / f) * f) {
        n = n / f;
        p++;
        pp(p) = f;
        psym = psym * f;
      }
      else {
        q++;
        qq(q) = f;
      }
      goto statement_10;
    }
  }
  statement_40:
  {
    std::ostringstream o;
    o << "FFTLIB: Largest factor exceeds " << pmax << ".  N = " << pts << ".";
    throw std::runtime_error(o.str());
  }
  statement_50:
  {
    std::ostringstream o;
    o << "FFTLIB: Factor count exceeds " << nest << ".  N = " << pts << ".";
    throw std::runtime_error(o.str());
  }
  statement_60:
  r = 1;
  if (q == 0) {
    r = 0;
  }
  if (p >= 1) {
    FEM_DO(j, 1, p) {
      jj = p + 1 - j;
      sym(j) = pp(jj);
      factor(j) = pp(jj);
      jj = p + q + j;
      factor(jj) = pp(j);
      jj = p + r + j;
      sym(jj) = pp(j);
    }
  }
  if (q >= 1) {
    FEM_DO(j, 1, q) {
      jj = p + j;
      unsym(j) = qq(j);
      factor(jj) = qq(j);
    }
    sym(p + 1) = pts / fem::pow(psym, 2);
  }
  jj = 2 * p + q;
  factor(jj + 1) = 0;
  ptwo = 1;
  j = 0;
  statement_90:
  j++;
  if (factor(j) != 0) {
    if (factor(j) == 2) {
      ptwo = ptwo * 2;
      factor(j) = 1;
      if (ptwo < twogrp) {
        if (factor(j + 1) == 2) {
          goto statement_90;
        }
      }
      factor(j) = ptwo;
      ptwo = 1;
    }
    goto statement_90;
  }
  if (p == 0) {
    r = 0;
  }
  jj = 2 * p + r;
  sym(jj + 1) = 0;
  if (q <= 1) {
    q = 0;
  }
  unsym(q + 1) = 0;
  error = false;
}

void
r2cftk(
  int const& n,
  int const& m,
  arr_ref<float> x0,
  arr_ref<float> y0,
  arr_ref<float> x1,
  arr_ref<float> y1,
  arr_cref<int> dim)
{
  x0(dimension(star));
  y0(dimension(star));
  x1(dimension(star));
  y1(dimension(star));
  dim(dimension(5));
  int nt = fem::int0;
  int sep = fem::int0;
  int l1 = fem::int0;
  int size = fem::int0;
  int k2 = fem::int0;
  int ns = fem::int0;
  int m2 = fem::int0;
  float fm2 = fem::float0;
  int mover2 = fem::int0;
  int mm2 = fem::int0;
  int j = fem::int0;
  bool fold = fem::bool0;
  int k0 = fem::int0;
  bool zero = fem::bool0;
  float angle = fem::float0;
  float c = fem::float0;
  float s = fem::float0;
  int kk = fem::int0;
  int l = fem::int0;
  int k1 = fem::int0;
  int k = fem::int0;
  float rs = fem::float0;
  float is = fem::float0;
  float ru = fem::float0;
  float iu = fem::float0;
  nt = dim(1);
  sep = dim(2);
  l1 = dim(3);
  size = dim(4) - 1;
  k2 = dim(5);
  ns = n * sep;
  m2 = m * 2;
  fm2 = fem::real(m2);
  mover2 = m / 2 + 1;
  mm2 = sep * m2;
  FEM_DO(j, 1, mover2) {
    fold = j > 1 && 2 * j < m + 2;
    k0 = (j - 1) * sep + 1;
    zero = j == 1;
    if (!zero) {
      angle = ccp4_fftlib_twopi * fem::real(j - 1) / fm2;
      c = fem::cos(angle);
      s = fem::sin(angle);
    }
    statement_10:
    FEM_DOSTEP(kk, k0, ns, mm2) {
      FEM_DOSTEP(l, kk, nt, l1) {
        k1 = l + size;
        FEM_DOSTEP(k, l, k1, k2) {
          rs = x0(k) + x1(k);
          is = y0(k) + y1(k);
          ru = x0(k) - x1(k);
          iu = y0(k) - y1(k);
          x0(k) = rs;
          y0(k) = is;
          if (zero) {
            x1(k) = ru;
            y1(k) = iu;
          }
          else {
            x1(k) = ru * c + iu * s;
            y1(k) = iu * c - ru * s;
          }
        }
      }
    }
    if (fold) {
      fold = false;
      k0 = (m + 1 - j) * sep + 1;
      c = -c;
      goto statement_10;
    }
  }
}

void
r3cftk(
  int const& n,
  int const& m,
  arr_ref<float> x0,
  arr_ref<float> y0,
  arr_ref<float> x1,
  arr_ref<float> y1,
  arr_ref<float> x2,
  arr_ref<float> y2,
  arr_cref<int> dim)
{
  x0(dimension(star));
  y0(dimension(star));
  x1(dimension(star));
  y1(dimension(star));
  x2(dimension(star));
  y2(dimension(star));
  dim(dimension(5));
  static const float a = -0.5f;
  static const float b = 0.86602540f;
  int nt = fem::int0;
  int sep = fem::int0;
  int l1 = fem::int0;
  int size = fem::int0;
  int k2 = fem::int0;
  int ns = fem::int0;
  int m3 = fem::int0;
  float fm3 = fem::float0;
  int mm3 = fem::int0;
  int mover2 = fem::int0;
  int j = fem::int0;
  bool fold = fem::bool0;
  int k0 = fem::int0;
  bool zero = fem::bool0;
  float angle = fem::float0;
  float c1 = fem::float0;
  float s1 = fem::float0;
  float c2 = fem::float0;
  float s2 = fem::float0;
  int kk = fem::int0;
  int l = fem::int0;
  int k1 = fem::int0;
  int k = fem::int0;
  float r0 = fem::float0;
  float i0 = fem::float0;
  float rs = fem::float0;
  float is = fem::float0;
  float ra = fem::float0;
  float ia = fem::float0;
  float rb = fem::float0;
  float ib = fem::float0;
  float r1 = fem::float0;
  float i1 = fem::float0;
  float r2 = fem::float0;
  float i2 = fem::float0;
  float t = fem::float0;
  nt = dim(1);
  sep = dim(2);
  l1 = dim(3);
  size = dim(4) - 1;
  k2 = dim(5);
  ns = n * sep;
  m3 = m * 3;
  fm3 = fem::real(m3);
  mm3 = sep * m3;
  mover2 = m / 2 + 1;
  FEM_DO(j, 1, mover2) {
    fold = j > 1 && 2 * j < m + 2;
    k0 = (j - 1) * sep + 1;
    zero = j == 1;
    if (!zero) {
      angle = ccp4_fftlib_twopi * fem::real(j - 1) / fm3;
      c1 = fem::cos(angle);
      s1 = fem::sin(angle);
      c2 = c1 * c1 - s1 * s1;
      s2 = s1 * c1 + c1 * s1;
    }
    statement_10:
    FEM_DOSTEP(kk, k0, ns, mm3) {
      FEM_DOSTEP(l, kk, nt, l1) {
        k1 = l + size;
        FEM_DOSTEP(k, l, k1, k2) {
          r0 = x0(k);
          i0 = y0(k);
          rs = x1(k) + x2(k);
          is = y1(k) + y2(k);
          x0(k) = r0 + rs;
          y0(k) = i0 + is;
          ra = rs * a + r0;
          ia = is * a + i0;
          rb = (x1(k) - x2(k)) * b;
          ib = (y1(k) - y2(k)) * b;
          if (zero) {
            x1(k) = ra + ib;
            y1(k) = ia - rb;
            x2(k) = ra - ib;
            y2(k) = ia + rb;
          }
          else {
            r1 = ra + ib;
            i1 = ia - rb;
            r2 = ra - ib;
            i2 = ia + rb;
            x1(k) = r1 * c1 + i1 * s1;
            y1(k) = i1 * c1 - r1 * s1;
            x2(k) = r2 * c2 + i2 * s2;
            y2(k) = i2 * c2 - r2 * s2;
          }
        }
      }
    }
    if (fold) {
      fold = false;
      k0 = (m + 1 - j) * sep + 1;
      t = c1 * a + s1 * b;
      s1 = c1 * b - s1 * a;
      c1 = t;
      t = c2 * a - s2 * b;
      s2 = -c2 * b - s2 * a;
      c2 = t;
      goto statement_10;
    }
  }
}

void
r4cftk(
  int const& n,
  int const& m,
  arr_ref<float> x0,
  arr_ref<float> y0,
  arr_ref<float> x1,
  arr_ref<float> y1,
  arr_ref<float> x2,
  arr_ref<float> y2,
  arr_ref<float> x3,
  arr_ref<float> y3,
  arr_cref<int> dim)
{
  x0(dimension(star));
  y0(dimension(star));
  x1(dimension(star));
  y1(dimension(star));
  x2(dimension(star));
  y2(dimension(star));
  x3(dimension(star));
  y3(dimension(star));
  dim(dimension(5));
  int nt = fem::int0;
  int sep = fem::int0;
  int l1 = fem::int0;
  int size = fem::int0;
  int k2 = fem::int0;
  int ns = fem::int0;
  int m4 = fem::int0;
  float fm4 = fem::float0;
  int mm4 = fem::int0;
  int mover2 = fem::int0;
  int j = fem::int0;
  bool fold = fem::bool0;
  int k0 = fem::int0;
  bool zero = fem::bool0;
  float angle = fem::float0;
  float c1 = fem::float0;
  float s1 = fem::float0;
  float c2 = fem::float0;
  float s2 = fem::float0;
  float c3 = fem::float0;
  float s3 = fem::float0;
  int kk = fem::int0;
  int l = fem::int0;
  int k1 = fem::int0;
  int k = fem::int0;
  float rs0 = fem::float0;
  float is0 = fem::float0;
  float ru0 = fem::float0;
  float iu0 = fem::float0;
  float rs1 = fem::float0;
  float is1 = fem::float0;
  float ru1 = fem::float0;
  float iu1 = fem::float0;
  float r1 = fem::float0;
  float i1 = fem::float0;
  float r2 = fem::float0;
  float i2 = fem::float0;
  float r3 = fem::float0;
  float i3 = fem::float0;
  float t = fem::float0;
  nt = dim(1);
  sep = dim(2);
  l1 = dim(3);
  size = dim(4) - 1;
  k2 = dim(5);
  ns = n * sep;
  m4 = m * 4;
  fm4 = fem::real(m4);
  mm4 = sep * m4;
  mover2 = m / 2 + 1;
  FEM_DO(j, 1, mover2) {
    fold = j > 1 && 2 * j < m + 2;
    k0 = (j - 1) * sep + 1;
    zero = j == 1;
    if (!zero) {
      angle = ccp4_fftlib_twopi * fem::real(j - 1) / fm4;
      c1 = fem::cos(angle);
      s1 = fem::sin(angle);
      c2 = c1 * c1 - s1 * s1;
      s2 = s1 * c1 + c1 * s1;
      c3 = c2 * c1 - s2 * s1;
      s3 = s2 * c1 + c2 * s1;
    }
    statement_10:
    FEM_DOSTEP(kk, k0, ns, mm4) {
      FEM_DOSTEP(l, kk, nt, l1) {
        k1 = l + size;
        FEM_DOSTEP(k, l, k1, k2) {
          rs0 = x0(k) + x2(k);
          is0 = y0(k) + y2(k);
          ru0 = x0(k) - x2(k);
          iu0 = y0(k) - y2(k);
          rs1 = x1(k) + x3(k);
          is1 = y1(k) + y3(k);
          ru1 = x1(k) - x3(k);
          iu1 = y1(k) - y3(k);
          x0(k) = rs0 + rs1;
          y0(k) = is0 + is1;
          if (zero) {
            x2(k) = ru0 + iu1;
            y2(k) = iu0 - ru1;
            x1(k) = rs0 - rs1;
            y1(k) = is0 - is1;
            x3(k) = ru0 - iu1;
            y3(k) = iu0 + ru1;
          }
          else {
            r1 = ru0 + iu1;
            i1 = iu0 - ru1;
            r2 = rs0 - rs1;
            i2 = is0 - is1;
            r3 = ru0 - iu1;
            i3 = iu0 + ru1;
            x2(k) = r1 * c1 + i1 * s1;
            y2(k) = i1 * c1 - r1 * s1;
            x1(k) = r2 * c2 + i2 * s2;
            y1(k) = i2 * c2 - r2 * s2;
            x3(k) = r3 * c3 + i3 * s3;
            y3(k) = i3 * c3 - r3 * s3;
          }
        }
      }
    }
    if (fold) {
      fold = false;
      k0 = (m + 1 - j) * sep + 1;
      t = c1;
      c1 = s1;
      s1 = t;
      c2 = -c2;
      t = c3;
      c3 = -s3;
      s3 = -t;
      goto statement_10;
    }
  }
}

void
r5cftk(
  int const& n,
  int const& m,
  arr_ref<float> x0,
  arr_ref<float> y0,
  arr_ref<float> x1,
  arr_ref<float> y1,
  arr_ref<float> x2,
  arr_ref<float> y2,
  arr_ref<float> x3,
  arr_ref<float> y3,
  arr_ref<float> x4,
  arr_ref<float> y4,
  arr_cref<int> dim)
{
  x0(dimension(star));
  y0(dimension(star));
  x1(dimension(star));
  y1(dimension(star));
  x2(dimension(star));
  y2(dimension(star));
  x3(dimension(star));
  y3(dimension(star));
  x4(dimension(star));
  y4(dimension(star));
  dim(dimension(5));
  static const float a1 = 0.30901699f;
  static const float b1 = 0.95105652f;
  static const float a2 = -0.80901699f;
  static const float b2 = 0.58778525f;
  int nt = fem::int0;
  int sep = fem::int0;
  int l1 = fem::int0;
  int size = fem::int0;
  int k2 = fem::int0;
  int ns = fem::int0;
  int m5 = fem::int0;
  float fm5 = fem::float0;
  int mm5 = fem::int0;
  int mover2 = fem::int0;
  int j = fem::int0;
  bool fold = fem::bool0;
  int k0 = fem::int0;
  bool zero = fem::bool0;
  float angle = fem::float0;
  float c1 = fem::float0;
  float s1 = fem::float0;
  float c2 = fem::float0;
  float s2 = fem::float0;
  float c3 = fem::float0;
  float s3 = fem::float0;
  float c4 = fem::float0;
  float s4 = fem::float0;
  int kk = fem::int0;
  int l = fem::int0;
  int k1 = fem::int0;
  int k = fem::int0;
  float r0 = fem::float0;
  float i0 = fem::float0;
  float rs1 = fem::float0;
  float is1 = fem::float0;
  float ru1 = fem::float0;
  float iu1 = fem::float0;
  float rs2 = fem::float0;
  float is2 = fem::float0;
  float ru2 = fem::float0;
  float iu2 = fem::float0;
  float ra1 = fem::float0;
  float ia1 = fem::float0;
  float ra2 = fem::float0;
  float ia2 = fem::float0;
  float rb1 = fem::float0;
  float ib1 = fem::float0;
  float rb2 = fem::float0;
  float ib2 = fem::float0;
  float r1 = fem::float0;
  float i1 = fem::float0;
  float r2 = fem::float0;
  float i2 = fem::float0;
  float r3 = fem::float0;
  float i3 = fem::float0;
  float r4 = fem::float0;
  float i4 = fem::float0;
  float t = fem::float0;
  nt = dim(1);
  sep = dim(2);
  l1 = dim(3);
  size = dim(4) - 1;
  k2 = dim(5);
  ns = n * sep;
  m5 = m * 5;
  fm5 = fem::real(m5);
  mm5 = sep * m5;
  mover2 = m / 2 + 1;
  FEM_DO(j, 1, mover2) {
    fold = j > 1 && 2 * j < m + 2;
    k0 = (j - 1) * sep + 1;
    zero = j == 1;
    if (!zero) {
      angle = ccp4_fftlib_twopi * fem::real(j - 1) / fm5;
      c1 = fem::cos(angle);
      s1 = fem::sin(angle);
      c2 = c1 * c1 - s1 * s1;
      s2 = s1 * c1 + c1 * s1;
      c3 = c2 * c1 - s2 * s1;
      s3 = s2 * c1 + c2 * s1;
      c4 = c2 * c2 - s2 * s2;
      s4 = s2 * c2 + c2 * s2;
    }
    statement_10:
    FEM_DOSTEP(kk, k0, ns, mm5) {
      FEM_DOSTEP(l, kk, nt, l1) {
        k1 = l + size;
        FEM_DOSTEP(k, l, k1, k2) {
          r0 = x0(k);
          i0 = y0(k);
          rs1 = x1(k) + x4(k);
          is1 = y1(k) + y4(k);
          ru1 = x1(k) - x4(k);
          iu1 = y1(k) - y4(k);
          rs2 = x2(k) + x3(k);
          is2 = y2(k) + y3(k);
          ru2 = x2(k) - x3(k);
          iu2 = y2(k) - y3(k);
          x0(k) = r0 + rs1 + rs2;
          y0(k) = i0 + is1 + is2;
          ra1 = rs1 * a1 + r0 + rs2 * a2;
          ia1 = is1 * a1 + i0 + is2 * a2;
          ra2 = rs1 * a2 + r0 + rs2 * a1;
          ia2 = is1 * a2 + i0 + is2 * a1;
          rb1 = ru1 * b1 + ru2 * b2;
          ib1 = iu1 * b1 + iu2 * b2;
          rb2 = ru1 * b2 - ru2 * b1;
          ib2 = iu1 * b2 - iu2 * b1;
          if (zero) {
            x1(k) = ra1 + ib1;
            y1(k) = ia1 - rb1;
            x2(k) = ra2 + ib2;
            y2(k) = ia2 - rb2;
            x3(k) = ra2 - ib2;
            y3(k) = ia2 + rb2;
            x4(k) = ra1 - ib1;
            y4(k) = ia1 + rb1;
          }
          else {
            r1 = ra1 + ib1;
            i1 = ia1 - rb1;
            r2 = ra2 + ib2;
            i2 = ia2 - rb2;
            r3 = ra2 - ib2;
            i3 = ia2 + rb2;
            r4 = ra1 - ib1;
            i4 = ia1 + rb1;
            x1(k) = r1 * c1 + i1 * s1;
            y1(k) = i1 * c1 - r1 * s1;
            x2(k) = r2 * c2 + i2 * s2;
            y2(k) = i2 * c2 - r2 * s2;
            x3(k) = r3 * c3 + i3 * s3;
            y3(k) = i3 * c3 - r3 * s3;
            x4(k) = r4 * c4 + i4 * s4;
            y4(k) = i4 * c4 - r4 * s4;
          }
        }
      }
    }
    if (fold) {
      fold = false;
      k0 = (m + 1 - j) * sep + 1;
      t = c1 * a1 + s1 * b1;
      s1 = c1 * b1 - s1 * a1;
      c1 = t;
      t = c2 * a2 + s2 * b2;
      s2 = c2 * b2 - s2 * a2;
      c2 = t;
      t = c3 * a2 - s3 * b2;
      s3 = -c3 * b2 - s3 * a2;
      c3 = t;
      t = c4 * a1 - s4 * b1;
      s4 = -c4 * b1 - s4 * a1;
      c4 = t;
      goto statement_10;
    }
  }
}

void
r8cftk(
  int const& n,
  int const& m,
  arr_ref<float> x0,
  arr_ref<float> y0,
  arr_ref<float> x1,
  arr_ref<float> y1,
  arr_ref<float> x2,
  arr_ref<float> y2,
  arr_ref<float> x3,
  arr_ref<float> y3,
  arr_ref<float> x4,
  arr_ref<float> y4,
  arr_ref<float> x5,
  arr_ref<float> y5,
  arr_ref<float> x6,
  arr_ref<float> y6,
  arr_ref<float> x7,
  arr_ref<float> y7,
  arr_cref<int> dim)
{
  x0(dimension(star));
  y0(dimension(star));
  x1(dimension(star));
  y1(dimension(star));
  x2(dimension(star));
  y2(dimension(star));
  x3(dimension(star));
  y3(dimension(star));
  x4(dimension(star));
  y4(dimension(star));
  x5(dimension(star));
  y5(dimension(star));
  x6(dimension(star));
  y6(dimension(star));
  x7(dimension(star));
  y7(dimension(star));
  dim(dimension(5));
  static const float e = 0.70710678f;
  int nt = fem::int0;
  int sep = fem::int0;
  int l1 = fem::int0;
  int size = fem::int0;
  int k2 = fem::int0;
  int ns = fem::int0;
  int m8 = fem::int0;
  float fm8 = fem::float0;
  int mm8 = fem::int0;
  int mover2 = fem::int0;
  int j = fem::int0;
  bool fold = fem::bool0;
  int k0 = fem::int0;
  bool zero = fem::bool0;
  float angle = fem::float0;
  float c1 = fem::float0;
  float s1 = fem::float0;
  float c2 = fem::float0;
  float s2 = fem::float0;
  float c3 = fem::float0;
  float s3 = fem::float0;
  float c4 = fem::float0;
  float s4 = fem::float0;
  float c5 = fem::float0;
  float s5 = fem::float0;
  float c6 = fem::float0;
  float s6 = fem::float0;
  float c7 = fem::float0;
  float s7 = fem::float0;
  int kk = fem::int0;
  int l = fem::int0;
  int k1 = fem::int0;
  int k = fem::int0;
  float rs0 = fem::float0;
  float is0 = fem::float0;
  float ru0 = fem::float0;
  float iu0 = fem::float0;
  float rs1 = fem::float0;
  float is1 = fem::float0;
  float ru1 = fem::float0;
  float iu1 = fem::float0;
  float rs2 = fem::float0;
  float is2 = fem::float0;
  float ru2 = fem::float0;
  float iu2 = fem::float0;
  float rs3 = fem::float0;
  float is3 = fem::float0;
  float ru3 = fem::float0;
  float iu3 = fem::float0;
  float rss0 = fem::float0;
  float iss0 = fem::float0;
  float rsu0 = fem::float0;
  float isu0 = fem::float0;
  float rss1 = fem::float0;
  float iss1 = fem::float0;
  float rsu1 = fem::float0;
  float isu1 = fem::float0;
  float rus0 = fem::float0;
  float ius0 = fem::float0;
  float ruu0 = fem::float0;
  float iuu0 = fem::float0;
  float rus1 = fem::float0;
  float ius1 = fem::float0;
  float ruu1 = fem::float0;
  float iuu1 = fem::float0;
  float t = fem::float0;
  float r1 = fem::float0;
  float i1 = fem::float0;
  float r2 = fem::float0;
  float i2 = fem::float0;
  float r3 = fem::float0;
  float i3 = fem::float0;
  float r4 = fem::float0;
  float i4 = fem::float0;
  float r5 = fem::float0;
  float i5 = fem::float0;
  float r6 = fem::float0;
  float i6 = fem::float0;
  float r7 = fem::float0;
  float i7 = fem::float0;
  nt = dim(1);
  sep = dim(2);
  l1 = dim(3);
  size = dim(4) - 1;
  k2 = dim(5);
  ns = n * sep;
  m8 = m * 8;
  fm8 = fem::real(m8);
  mm8 = sep * m8;
  mover2 = m / 2 + 1;
  FEM_DO(j, 1, mover2) {
    fold = j > 1 && 2 * j < m + 2;
    k0 = (j - 1) * sep + 1;
    zero = j == 1;
    if (!zero) {
      angle = ccp4_fftlib_twopi * fem::real(j - 1) / fm8;
      c1 = fem::cos(angle);
      s1 = fem::sin(angle);
      c2 = c1 * c1 - s1 * s1;
      s2 = s1 * c1 + c1 * s1;
      c3 = c2 * c1 - s2 * s1;
      s3 = s2 * c1 + c2 * s1;
      c4 = c2 * c2 - s2 * s2;
      s4 = s2 * c2 + c2 * s2;
      c5 = c4 * c1 - s4 * s1;
      s5 = s4 * c1 + c4 * s1;
      c6 = c4 * c2 - s4 * s2;
      s6 = s4 * c2 + c4 * s2;
      c7 = c4 * c3 - s4 * s3;
      s7 = s4 * c3 + c4 * s3;
    }
    statement_10:
    FEM_DOSTEP(kk, k0, ns, mm8) {
      FEM_DOSTEP(l, kk, nt, l1) {
        k1 = l + size;
        FEM_DOSTEP(k, l, k1, k2) {
          rs0 = x0(k) + x4(k);
          is0 = y0(k) + y4(k);
          ru0 = x0(k) - x4(k);
          iu0 = y0(k) - y4(k);
          rs1 = x1(k) + x5(k);
          is1 = y1(k) + y5(k);
          ru1 = x1(k) - x5(k);
          iu1 = y1(k) - y5(k);
          rs2 = x2(k) + x6(k);
          is2 = y2(k) + y6(k);
          ru2 = x2(k) - x6(k);
          iu2 = y2(k) - y6(k);
          rs3 = x3(k) + x7(k);
          is3 = y3(k) + y7(k);
          ru3 = x3(k) - x7(k);
          iu3 = y3(k) - y7(k);
          rss0 = rs0 + rs2;
          iss0 = is0 + is2;
          rsu0 = rs0 - rs2;
          isu0 = is0 - is2;
          rss1 = rs1 + rs3;
          iss1 = is1 + is3;
          rsu1 = rs1 - rs3;
          isu1 = is1 - is3;
          rus0 = ru0 - iu2;
          ius0 = iu0 + ru2;
          ruu0 = ru0 + iu2;
          iuu0 = iu0 - ru2;
          rus1 = ru1 - iu3;
          ius1 = iu1 + ru3;
          ruu1 = ru1 + iu3;
          iuu1 = iu1 - ru3;
          t = (rus1 + ius1) * e;
          ius1 = (ius1 - rus1) * e;
          rus1 = t;
          t = (ruu1 + iuu1) * e;
          iuu1 = (iuu1 - ruu1) * e;
          ruu1 = t;
          x0(k) = rss0 + rss1;
          y0(k) = iss0 + iss1;
          if (zero) {
            x4(k) = ruu0 + ruu1;
            y4(k) = iuu0 + iuu1;
            x2(k) = rsu0 + isu1;
            y2(k) = isu0 - rsu1;
            x6(k) = rus0 + ius1;
            y6(k) = ius0 - rus1;
            x1(k) = rss0 - rss1;
            y1(k) = iss0 - iss1;
            x5(k) = ruu0 - ruu1;
            y5(k) = iuu0 - iuu1;
            x3(k) = rsu0 - isu1;
            y3(k) = isu0 + rsu1;
            x7(k) = rus0 - ius1;
            y7(k) = ius0 + rus1;
          }
          else {
            r1 = ruu0 + ruu1;
            i1 = iuu0 + iuu1;
            r2 = rsu0 + isu1;
            i2 = isu0 - rsu1;
            r3 = rus0 + ius1;
            i3 = ius0 - rus1;
            r4 = rss0 - rss1;
            i4 = iss0 - iss1;
            r5 = ruu0 - ruu1;
            i5 = iuu0 - iuu1;
            r6 = rsu0 - isu1;
            i6 = isu0 + rsu1;
            r7 = rus0 - ius1;
            i7 = ius0 + rus1;
            x4(k) = r1 * c1 + i1 * s1;
            y4(k) = i1 * c1 - r1 * s1;
            x2(k) = r2 * c2 + i2 * s2;
            y2(k) = i2 * c2 - r2 * s2;
            x6(k) = r3 * c3 + i3 * s3;
            y6(k) = i3 * c3 - r3 * s3;
            x1(k) = r4 * c4 + i4 * s4;
            y1(k) = i4 * c4 - r4 * s4;
            x5(k) = r5 * c5 + i5 * s5;
            y5(k) = i5 * c5 - r5 * s5;
            x3(k) = r6 * c6 + i6 * s6;
            y3(k) = i6 * c6 - r6 * s6;
            x7(k) = r7 * c7 + i7 * s7;
            y7(k) = i7 * c7 - r7 * s7;
          }
        }
      }
    }
    if (fold) {
      fold = false;
      k0 = (m + 1 - j) * sep + 1;
      t = (c1 + s1) * e;
      s1 = (c1 - s1) * e;
      c1 = t;
      t = s2;
      s2 = c2;
      c2 = t;
      t = (-c3 + s3) * e;
      s3 = (c3 + s3) * e;
      c3 = t;
      c4 = -c4;
      t = -(c5 + s5) * e;
      s5 = (-c5 + s5) * e;
      c5 = t;
      t = -s6;
      s6 = -c6;
      c6 = t;
      t = (c7 - s7) * e;
      s7 = -(c7 + s7) * e;
      c7 = t;
      goto statement_10;
    }
  }
}

void
rpcftk(
  int const& n,
  int const& m,
  int const& p,
  int const& r,
  arr_ref<float, 2> x,
  arr_ref<float, 2> y,
  arr_cref<int> dim)
{
  x(dimension(r, p));
  y(dimension(r, p));
  dim(dimension(5));
  int nt = fem::int0;
  int sep = fem::int0;
  int l1 = fem::int0;
  int size = fem::int0;
  int k2 = fem::int0;
  int ns = fem::int0;
  int mover2 = fem::int0;
  int mp = fem::int0;
  float fmp = fem::float0;
  int mmp = fem::int0;
  int pp = fem::int0;
  int pm = fem::int0;
  float fp = fem::float0;
  float fu = fem::float0;
  int u = fem::int0;
  float angle = fem::float0;
  int jj = fem::int0;
  arr_1d<18, float> a;
  arr_1d<18, float> b;
  int v = fem::int0;
  arr_2d<9,9, float> aa;
  arr_2d<9,9, float> bb;
  int j = fem::int0;
  bool fold = fem::bool0;
  int k0 = fem::int0;
  bool zero = fem::bool0;
  arr_1d<18, float> c;
  arr_1d<18, float> s;
  int kk = fem::int0;
  int l = fem::int0;
  int k1 = fem::int0;
  int k = fem::int0;
  float xt = fem::float0;
  float yt = fem::float0;
  float rs = fem::float0;
  float is = fem::float0;
  float ru = fem::float0;
  float iu = fem::float0;
  arr_1d<9, float> ra;
  arr_1d<9, float> ia;
  arr_1d<9, float> rb;
  arr_1d<9, float> ib;
  float t = fem::float0;
  nt = dim(1);
  sep = dim(2);
  l1 = dim(3);
  size = dim(4) - 1;
  k2 = dim(5);
  ns = n * sep;
  mover2 = m / 2 + 1;
  mp = m * p;
  fmp = fem::real(mp);
  mmp = sep * mp;
  pp = p / 2;
  pm = p - 1;
  fp = fem::real(p);
  fu = 0.0f;
  FEM_DO(u, 1, pp) {
    fu += 1.0f;
    angle = ccp4_fftlib_twopi * fu / fp;
    jj = p - u;
    a(u) = fem::cos(angle);
    b(u) = fem::sin(angle);
    a(jj) = a(u);
    b(jj) = -b(u);
  }
  FEM_DO(u, 1, pp) {
    FEM_DO(v, 1, pp) {
      jj = u * v - u * v / p * p;
      aa(v, u) = a(jj);
      bb(v, u) = b(jj);
    }
  }
  FEM_DO(j, 1, mover2) {
    fold = j > 1 && 2 * j < m + 2;
    k0 = (j - 1) * sep + 1;
    zero = j == 1;
    if (!zero) {
      angle = ccp4_fftlib_twopi * fem::real(j - 1) / fmp;
      c(1) = fem::cos(angle);
      s(1) = fem::sin(angle);
      FEM_DO(u, 2, pm) {
        c(u) = c(u - 1) * c(1) - s(u - 1) * s(1);
        s(u) = s(u - 1) * c(1) + c(u - 1) * s(1);
      }
    }
    statement_50:
    FEM_DOSTEP(kk, k0, ns, mmp) {
      FEM_DOSTEP(l, kk, nt, l1) {
        k1 = l + size;
        FEM_DOSTEP(k, l, k1, k2) {
          xt = x(k, 1);
          yt = y(k, 1);
          rs = x(k, 2) + x(k, p);
          is = y(k, 2) + y(k, p);
          ru = x(k, 2) - x(k, p);
          iu = y(k, 2) - y(k, p);
          FEM_DO(u, 1, pp) {
            ra(u) = aa(u, 1) * rs + xt;
            ia(u) = aa(u, 1) * is + yt;
            rb(u) = bb(u, 1) * ru;
            ib(u) = bb(u, 1) * iu;
          }
          xt += rs;
          yt += is;
          FEM_DO(u, 2, pp) {
            jj = p - u;
            rs = x(k, u + 1) + x(k, jj + 1);
            is = y(k, u + 1) + y(k, jj + 1);
            ru = x(k, u + 1) - x(k, jj + 1);
            iu = y(k, u + 1) - y(k, jj + 1);
            xt += rs;
            yt += is;
            FEM_DO(v, 1, pp) {
              ra(v) += aa(v, u) * rs;
              ia(v) += aa(v, u) * is;
              rb(v) += bb(v, u) * ru;
              ib(v) += bb(v, u) * iu;
            }
          }
          x(k, 1) = xt;
          y(k, 1) = yt;
          FEM_DO(u, 1, pp) {
            jj = p - u;
            if (zero) {
              x(k, u + 1) = ra(u) + ib(u);
              y(k, u + 1) = ia(u) - rb(u);
              x(k, jj + 1) = ra(u) - ib(u);
              y(k, jj + 1) = ia(u) + rb(u);
            }
            else {
              xt = ra(u) + ib(u);
              yt = ia(u) - rb(u);
              x(k, u + 1) = c(u) * xt + s(u) * yt;
              y(k, u + 1) = c(u) * yt - s(u) * xt;
              xt = ra(u) - ib(u);
              yt = ia(u) + rb(u);
              x(k, jj + 1) = c(jj) * xt + s(jj) * yt;
              y(k, jj + 1) = c(jj) * yt - s(jj) * xt;
            }
          }
        }
      }
    }
    if (fold) {
      fold = false;
      k0 = (m + 1 - j) * sep + 1;
      FEM_DO(u, 1, pm) {
        t = c(u) * a(u) + s(u) * b(u);
        s(u) = -s(u) * a(u) + c(u) * b(u);
        c(u) = t;
      }
      goto statement_50;
    }
  }
}

void
mdftkd(
  int const& n,
  arr_cref<int> factor,
  arr_cref<int> dim,
  arr_ref<float> x,
  arr_ref<float> y)
{
  factor(dimension(15));
  dim(dimension(5));
  x(dimension(star));
  y(dimension(star));
  int s = fem::int0;
  int f = fem::int0;
  int m = fem::int0;
  int p = fem::int0;
  int r = fem::int0;
  s = dim(2);
  f = 0;
  m = n;
  statement_10:
  f++;
  p = factor(f);
  if (p == 0) {
    return;
  }
  else {
    m = m / p;
    r = m * s;
    if (p <= 8) {
      switch (p) {
        case 1: goto statement_10;
        case 2: goto statement_20;
        case 3: goto statement_30;
        case 4: goto statement_40;
        case 5: goto statement_50;
        case 6: goto statement_80;
        case 7: goto statement_70;
        case 8: goto statement_60;
        default: break;
      }
      goto statement_80;
      statement_20:
      r2cftk(n, m, x(1), y(1), x(r + 1), y(r + 1), dim);
      goto statement_10;
      statement_30:
      r3cftk(n, m, x(1), y(1), x(r + 1), y(r + 1), x(2 * r + 1), y(2 * r + 1), dim);
      goto statement_10;
      statement_40:
      r4cftk(n, m, x(1), y(1), x(r + 1), y(r + 1), x(2 * r + 1), y(2 * r + 1), x(3 * r + 1), y(3 * r + 1), dim);
      goto statement_10;
      statement_50:
      r5cftk(n, m, x(1), y(1), x(r + 1), y(r + 1), x(2 * r + 1), y(2 * r + 1), x(3 * r + 1), y(3 * r + 1), x(4 * r + 1), y(4 * r + 1), dim);
      goto statement_10;
      statement_60:
      r8cftk(n, m, x(1), y(1), x(r + 1), y(r + 1), x(2 * r + 1), y(2 * r + 1), x(3 * r + 1), y(3 * r + 1), x(4 * r + 1), y(4 * r + 1), x(5 * r + 1), y(5 * r + 1), x(6 * r + 1), y(6 * r + 1), x(7 * r + 1), y(7 * r + 1), dim);
      goto statement_10;
    }
    statement_70:
    rpcftk(n, m, p, r, x, y, dim);
    goto statement_10;
  }
  statement_80:
  throw std::runtime_error("TRANSFER ERROR DETECTED IN MDFTKD");
}

void
diprp(
  int const& pts,
  arr_cref<int> sym,
  int const& psym,
  arr_cref<int> unsym,
  arr_cref<int> dim,
  arr_ref<float> x,
  arr_ref<float> y)
{
  sym(dimension(15));
  unsym(dimension(15));
  dim(dimension(5));
  x(dimension(star));
  y(dimension(star));
  arr_1d<14, int> u;
  arr_1d<14, int> s;
  int& al = u(1);
  int& bl = u(2);
  int& cl = u(3);
  int& dl = u(4);
  int& el = u(5);
  int& fl = u(6);
  int& gl = u(7);
  int& hl = u(8);
  int& il = u(9);
  int& jl = u(10);
  int& kl = u(11);
  int& ll = u(12);
  int& ml = u(13);
  int& nl = u(14);
  int& bs = s(2);
  int& cs = s(3);
  int& ds = s(4);
  int& es = s(5);
  int& fs = s(6);
  int& gs = s(7);
  int& hs = s(8);
  int& is = s(9);
  int& js = s(10);
  int& ks = s(11);
  int& ls = s(12);
  int& ms = s(13);
  int& ns = s(14);
  int nest = fem::int0;
  int nt = fem::int0;
  int sep = fem::int0;
  int p2 = fem::int0;
  int size = fem::int0;
  int p4 = fem::int0;
  int j = fem::int0;
  int n = fem::int0;
  int jj = fem::int0;
  int a = fem::int0;
  int b = fem::int0;
  int c = fem::int0;
  int d = fem::int0;
  int e = fem::int0;
  int f = fem::int0;
  int g = fem::int0;
  int h = fem::int0;
  int i = fem::int0;
  int k = fem::int0;
  int l = fem::int0;
  int m = fem::int0;
  int delta = fem::int0;
  int p1 = fem::int0;
  int p0 = fem::int0;
  int p3 = fem::int0;
  int p = fem::int0;
  int p5 = fem::int0;
  float t = fem::float0;
  int punsym = fem::int0;
  int mult = fem::int0;
  int test = fem::int0;
  int lk = fem::int0;
  int dk = fem::int0;
  int mods = fem::int0;
  bool onemod = fem::bool0;
  arr_1d<14, int> modulo;
  int kk = fem::int0;
  nest = 14;
  nt = dim(1);
  sep = dim(2);
  p2 = dim(3);
  size = dim(4) - 1;
  p4 = dim(5);
  if (sym(1) != 0) {
    FEM_DO(j, 1, nest) {
      u(j) = 1;
      s(j) = 1;
    }
    n = pts;
    FEM_DO(j, 1, nest) {
      if (sym(j) == 0) {
        goto statement_30;
      }
      else {
        jj = nest + 1 - j;
        u(jj) = n;
        s(jj) = n / sym(j);
        n = n / sym(j);
      }
    }
    statement_30:
    jj = 0;
    FEM_DO(a, 1, al) {
      FEM_DOSTEP(b, a, bl, bs) {
        FEM_DOSTEP(c, b, cl, cs) {
          FEM_DOSTEP(d, c, dl, ds) {
            FEM_DOSTEP(e, d, el, es) {
              FEM_DOSTEP(f, e, fl, fs) {
                FEM_DOSTEP(g, f, gl, gs) {
                  FEM_DOSTEP(h, g, hl, hs) {
                    FEM_DOSTEP(i, h, il, is) {
                      FEM_DOSTEP(j, i, jl, js) {
                        FEM_DOSTEP(k, j, kl, ks) {
                          FEM_DOSTEP(l, k, ll, ls) {
                            FEM_DOSTEP(m, l, ml, ms) {
                              FEM_DOSTEP(n, m, nl, ns) {
                                jj++;
                                if (jj < n) {
                                  delta = (n - jj) * sep;
                                  p1 = (jj - 1) * sep + 1;
                                  FEM_DOSTEP(p0, p1, nt, p2) {
                                    p3 = p0 + size;
                                    FEM_DOSTEP(p, p0, p3, p4) {
                                      p5 = p + delta;
                                      t = x(p);
                                      x(p) = x(p5);
                                      x(p5) = t;
                                      t = y(p);
                                      y(p) = y(p5);
                                      y(p5) = t;
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (unsym(1) != 0) {
    punsym = pts / fem::pow(psym, 2);
    mult = punsym / unsym(1);
    test = (unsym(1) * unsym(2) - 1) * mult * psym;
    lk = mult;
    dk = mult;
    FEM_DO(k, 2, nest) {
      if (unsym(k) == 0) {
        goto statement_210;
      }
      else {
        lk = unsym(k - 1) * lk;
        dk = dk / unsym(k);
        u(k) = (lk - dk) * psym;
        mods = k;
      }
    }
    statement_210:
    onemod = mods < 3;
    if (!onemod) {
      FEM_DO(j, 3, mods) {
        jj = mods + 3 - j;
        modulo(jj) = u(j);
      }
    }
    modulo(2) = u(2);
    jl = (punsym - 3) * psym;
    ms = punsym * psym;
    FEM_DOSTEP(j, psym, jl, psym) {
      k = j;
      statement_230:
      k = k * mult;
      if (!onemod) {
        FEM_DO(i, 3, mods) {
          k = k - (k / modulo(i)) * modulo(i);
        }
      }
      if (k >= test) {
        k = k - (k / modulo(2)) * modulo(2) + modulo(2);
      }
      else {
        k = k - (k / modulo(2)) * modulo(2);
      }
      if (k < j) {
        goto statement_230;
      }
      if (k != j) {
        delta = (k - j) * sep;
        FEM_DO(l, 1, psym) {
          FEM_DOSTEP(m, l, pts, ms) {
            p1 = (m + j - 1) * sep + 1;
            FEM_DOSTEP(p0, p1, nt, p2) {
              p3 = p0 + size;
              FEM_DOSTEP(jj, p0, p3, p4) {
                kk = jj + delta;
                t = x(jj);
                x(jj) = x(kk);
                x(kk) = t;
                t = y(jj);
                y(jj) = y(kk);
                y(kk) = t;
              }
            }
          }
        }
      }
    }
  }
}

void
cmplft(
  arr_ref<float> x,
  arr_ref<float> y,
  int const& n,
  arr_cref<int> d)
{
  x(dimension(star));
  y(dimension(star));
  d(dimension(5));
  int pmax = 19;
  int twogrp = 8;
  arr_1d<15, int> factor;
  arr_1d<15, int> sym;
  arr_1d<15, int> unsym;
  int psym = fem::int0;
  bool error = fem::bool0;
  if (n > 1) {
    srfp(n, pmax, twogrp, factor, sym, psym, unsym, error);
    if (error) {
      std::ostringstream o;
      o << "FFTLIB: invalid number of points for CMPL FT.  N = " << n;
      throw std::runtime_error(o.str());
    }
    else {
      mdftkd(n, factor, d, x, y);
      diprp(n, sym, psym, unsym, d, x, y);
    }
  }
}

void
hermft(
  arr_ref<float> x,
  arr_ref<float> y,
  int const& n,
  arr_cref<int> dim)
{
  x(dimension(star));
  y(dimension(star));
  dim(dimension(5));
  float twon = fem::real(2 * n);
  int nt = dim(1);
  int d2 = dim(2);
  int d3 = dim(3);
  int d4 = dim(4) - 1;
  int d5 = dim(5);
  int i0 = fem::int0;
  int i1 = fem::int0;
  int i = fem::int0;
  float a = fem::float0;
  float b = fem::float0;
  FEM_DOSTEP(i0, 1, nt, d3) {
    i1 = i0 + d4;
    FEM_DOSTEP(i, i0, i1, d5) {
      a = x(i);
      b = y(i);
      x(i) = a + b;
      y(i) = a - b;
    }
  }
  int nover2 = n / 2 + 1;
  float angle = fem::float0;
  float co = fem::float0;
  float si = fem::float0;
  int k = fem::int0;
  int k1 = fem::int0;
  int i2 = fem::int0;
  int j = fem::int0;
  float c = fem::float0;
  float d = fem::float0;
  float e = fem::float0;
  float f = fem::float0;
  if (nover2 >= 2) {
    FEM_DO(i0, 2, nover2) {
      angle = fem::real(i0 - 1) * ccp4_fftlib_twopi / twon;
      co = fem::cos(angle);
      si = fem::sin(angle);
      k = (n + 2 - 2 * i0) * d2;
      k1 = (i0 - 1) * d2 + 1;
      FEM_DOSTEP(i1, k1, nt, d3) {
        i2 = i1 + d4;
        FEM_DOSTEP(i, i1, i2, d5) {
          j = i + k;
          a = x(i) + x(j);
          b = x(i) - x(j);
          c = y(i) + y(j);
          d = y(i) - y(j);
          e = b * co + c * si;
          f = b * si - c * co;
          x(i) = a + f;
          x(j) = a - f;
          y(i) = e + d;
          y(j) = e - d;
        }
      }
    }
    cmplft(x, y, n, dim);
  }
}

void
realft(
  arr_ref<float> even,
  arr_ref<float> odd,
  int const& n,
  arr_cref<int> dim)
{
  even(dimension(star));
  odd(dimension(star));
  dim(dimension(5));
  float twon = fem::real(2 * n);
  cmplft(even, odd, n, dim);
  int nt = dim(1);
  int d2 = dim(2);
  int d3 = dim(3);
  int d4 = dim(4) - 1;
  int d5 = dim(5);
  int nover2 = n / 2 + 1;
  int i = fem::int0;
  float angle = fem::float0;
  float co = fem::float0;
  float si = fem::float0;
  int i0 = fem::int0;
  int j = fem::int0;
  int i1 = fem::int0;
  int i2 = fem::int0;
  int k = fem::int0;
  int l = fem::int0;
  float a = fem::float0;
  float c = fem::float0;
  float b = fem::float0;
  float d = fem::float0;
  float e = fem::float0;
  float f = fem::float0;
  if (nover2 >= 2) {
    FEM_DO(i, 2, nover2) {
      angle = fem::real(i - 1) * ccp4_fftlib_twopi / twon;
      co = fem::cos(angle);
      si = fem::sin(angle);
      i0 = (i - 1) * d2 + 1;
      j = (n + 2 - 2 * i) * d2;
      FEM_DOSTEP(i1, i0, nt, d3) {
        i2 = i1 + d4;
        FEM_DOSTEP(k, i1, i2, d5) {
          l = k + j;
          a = (even(l) + even(k)) / 2.0f;
          c = (even(l) - even(k)) / 2.0f;
          b = (odd(l) + odd(k)) / 2.0f;
          d = (odd(l) - odd(k)) / 2.0f;
          e = c * si + b * co;
          f = c * co - b * si;
          even(k) = a + e;
          even(l) = a - e;
          odd(k) = f - d;
          odd(l) = f + d;
        }
      }
    }
  }
  if (n >= 1) {
    j = n * d2;
    FEM_DOSTEP(i1, 1, nt, d3) {
      i2 = i1 + d4;
      FEM_DOSTEP(k, i1, i2, d5) {
        l = k + j;
        even(l) = even(k) - odd(k);
        odd(l) = 0.0f;
        even(k) += odd(k);
        odd(k) = 0.0f;
      }
    }
  }
}

} // namespace ccp4_fftlib_fem
