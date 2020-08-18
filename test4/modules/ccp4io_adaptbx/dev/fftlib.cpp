// NOT USED -- INCOMPLETE -- DOES NOT COMPILE
// here for future reference

//
//     fftlib.f: fast-fourier transform library routines
//     Copyright (C)  Lynn Ten Eyck
//
//     This library is free software: you can redistribute it and/or
//     modify it under the terms of the GNU Lesser General Public License
//     version 3, modified in accordance with the provisions of the
//     license to address the requirements of UK law.
//
//     You should have received a copy of the modified GNU Lesser General
//     Public License along with this library.  If not, copies may be
//     downloaded from http://www.ccp4.ac.uk/ccp4license.php
//
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU Lesser General Public License for more details.
//
//-- FFT81       F77TRNFM.FOR                        13/09/85    JWC
//
//
//**** FOLLOWING ARE ROUTINES USED BY TEN EYCK'S FFT PROGRAMS***
//

#include <ccp4io_adaptbx/dev/fftlib.hpp>

#include <ObjexxFCL/Fmath.hh>

#include <algorithm>
#include <iostream>

namespace ccp4io_dev {

namespace {

template <typename ElementType>
struct raw_ref1
{
  ElementType* begin;

  inline
  raw_ref1(ElementType& begin_)
  :
    begin(&begin_)
  {}

  inline
  ElementType&
  operator()(int i) const
  {
    return begin[i-1];
  }
};

template <typename ElementType>
struct raw_ref2
{
  ElementType* begin;
  int nr;

  inline
  raw_ref2(ElementType* begin_, int nr_, int /*nc_*/)
  :
    begin(begin_),
    nr(nr_)
  {}

  inline
  ElementType&
  operator()(int i, int j) const
  {
    return begin[(j-1) * nr + (i-1)];
  }
};

void
srfp(
  int const pts,
  int const pmax,
  int const twogrp,
  tiny_ref<int, 15> factor,
  tiny_ref<int, 15> sym,
  int & psym,
  tiny_ref<int, 15> unsym,
  bool & error
)
{
  using std::cout;
//     ==================================================================
//
//---- Symmetrized reordering factoring programme

//     .. Local Scalars ..
  int f, jj, jlast = 0, n, nest, p, ptwo, q, r;

//     .. Local Arrays ..
  tiny_array<int, 14> pp;
  tiny_array<int, 7> qq;
  pp.fill(0); //Objexx Zero initialization added: Seems necessary
  qq.fill(0); //Objexx Zero initialization added: Seems necessary

  nest = 14;
  n = pts;
  psym = 1;
  f = 2;
  p = 0;
  q = 0;
  while ( n > 1 ) {
    for ( int j = f; j <= pmax; ++j ) {
      if ( n == (n/j)*j ) {
        jlast = j;
        goto L30;
      }
    }
    cout << /* FORMAT(" FFTLIB: Largest factor exceeds ",i3,".  N = ",i6,'.') */ boost::format( " FFTLIB: Largest factor exceeds %3d.  N = %6d.\n" ) % pmax % pts;
    error = true;
    goto L100;
L30:
    if ( 2*p + q >= nest ) {
      cout << /* FORMAT(" FFTLIB: Factor count exceeds ",i3,".  N = ",i6,'.') */ boost::format( " FFTLIB: Factor count exceeds %3d.  N = %6d.\n" ) % nest % pts;
      error = true;
      goto L100;
    } else {
      f = jlast;
      n /= f;
      if ( n == (n/f)*f ) {
        n /= f;
        ++p;
        pp(p) = f;
        psym *= f;
      } else {
        ++q;
        qq(q) = f;
      }
    }
  }

  r = 1;
  if ( q == 0 ) r = 0;
  if ( p >= 1 ) {
    for ( int j = 1; j <= p; ++j ) {
      jj = p + 1 - j;
      sym(j) = pp(jj);
      factor(j) = pp(jj);
      jj = p + q + j;
      factor(jj) = pp(j);
      jj = p + r + j;
      sym(jj) = pp(j);
    }
  }
  if ( q >= 1 ) {
    for ( int j = 1; j <= q; ++j ) {
      jj = p + j;
      unsym(j) = qq(j);
      factor(jj) = qq(j);
    }
    sym(p+1) = pts/square( psym );
  }
  jj = 2*p + q;
  factor(jj + 1) = 0;
  ptwo = 1;
  { // Scope for j
    int j = 0;
    while ( true ) {
      ++j;
      if ( factor(j) != 0 ) {
        if ( factor(j) == 2 ) {
          ptwo *= 2;
          factor(j) = 1;
          if ( ptwo < twogrp ) {
            if ( factor(j + 1) == 2 ) continue;
          }
          factor(j) = ptwo;
          ptwo = 1;
        }
        continue;
      }
      if ( p == 0 ) r = 0;
      jj = 2*p + r;
      sym(jj+1) = 0;
      if ( q <= 1 ) q = 0;
      unsym(q+1) = 0;
      error = false;

  //---- Format statements
  //
  //FMT9000    FORMAT (" FFTLIB: Largest factor exceeds ",i3,".  N = ",i6,'.')
  //FMT9010    FORMAT (" FFTLIB: Factor count exceeds ",i3,".  N = ",i6,'.')
      break;
    }
  }
L100:;
}

void
r2cftk(
  int const n,
  int & m,
  raw_ref1<float> x0,
  raw_ref1<float> y0,
  raw_ref1<float> x1,
  raw_ref1<float> y1,
  int5a dim
)
{
//     ==============================================
//
//---- Radix 2 multi-dimensional complex fourier transform kernel

//     .. Local Scalars ..
  static float const twopi = 6.283185f;
  float angle, c, fm2, is, iu, rs, ru, s;
  int k0, k1, k2, l1, m2, mm2, mover2, ns, nt, sep, size;
  bool fold, zero;

  nt = dim(1);
  sep = dim(2);
  l1 = dim(3);
  size = dim(4) - 1;
  k2 = dim(5);
  ns = n*sep;
  m2 = m*2;
  fm2 = float(m2);
  mover2 = m/2 + 1;
  mm2 = sep*m2;

  // avoid "may be used uninitialized" warnings
  c = s = 0.f;
  for ( int j = 1; j <= mover2; ++j ) {
    fold = j > 1 && 2*j < m + 2;
    k0 = (j - 1)*sep + 1;
    zero = j == 1;
    if ( ! zero ) {
      angle = twopi*float(j - 1)/fm2;
      c = std::cos(angle);
      s = std::sin(angle);
    }
    while ( true ) {
      for ( int kk = k0, smm2 = sign( mm2 ); kk * smm2 <= smm2 * ns; kk += mm2 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
        for ( int l = kk, sl1 = sign( l1 ); l * sl1 <= sl1 * nt; l += l1 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
          k1 = l + size;
          for ( int k = l, sk2 = sign( k2 ); k * sk2 <= sk2 * k1; k += k2 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
            rs = x0(k) + x1(k);
            is = y0(k) + y1(k);
            ru = x0(k) - x1(k);
            iu = y0(k) - y1(k);
            x0(k) = rs;
            y0(k) = is;
            if ( zero ) {
              x1(k) = ru;
              y1(k) = iu;
            } else {
              x1(k) = ru*c + iu*s;
              y1(k) = iu*c - ru*s;
            }
          }
        }
      }
      if ( fold ) {
        fold = false;
        k0 = (m + 1 - j)*sep + 1;
        c = -c;
        continue;
      }
      break;
    }
  }

}

void
r3cftk(
  int const n,
  int & m,
  raw_ref1<float> x0,
  raw_ref1<float> y0,
  raw_ref1<float> x1,
  raw_ref1<float> y1,
  raw_ref1<float> x2,
  raw_ref1<float> y2,
  int5a dim
)
{
//     ======================================================
//
//---- Radix 3 multi-dimensional complex fourier transform kernel

//     .. Local Scalars ..
  static const float twopi = 6.283185f;
  static float const a = -0.5f;
  static float const b = 0.86602540f;
  float angle, c1, c2, fm3, i0, i1, i2, ia, ib, is, r0, r1, r2, ra, rb, rs, s1, s2, t;
  int k0, k1, k2, l1, m3, mm3, mover2, ns, nt, sep, size;
  bool fold, zero;

  nt = dim(1);
  sep = dim(2);
  l1 = dim(3);
  size = dim(4) - 1;
  k2 = dim(5);
  ns = n*sep;
  m3 = m*3;
  fm3 = float(m3);
  mm3 = sep*m3;
  mover2 = m/2 + 1;

  // avoid "may be used uninitialized" warnings
  c1 = s1 = c2 = s2 = 0.f;
  for ( int j = 1; j <= mover2; ++j ) {
    fold = j > 1 && 2*j < m + 2;
    k0 = (j - 1)*sep + 1;
    zero = j == 1;
    if ( ! zero ) {
      angle = twopi*float(j - 1)/fm3;
      c1 = std::cos(angle);
      s1 = std::sin(angle);
      c2 = c1*c1 - s1*s1;
      s2 = s1*c1 + c1*s1;
    }
    while ( true ) {
      for ( int kk = k0, smm3 = sign( mm3 ); kk * smm3 <= smm3 * ns; kk += mm3 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
        for ( int l = kk, sl1 = sign( l1 ); l * sl1 <= sl1 * nt; l += l1 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
          k1 = l + size;
          for ( int k = l, sk2 = sign( k2 ); k * sk2 <= sk2 * k1; k += k2 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
            r0 = x0(k);
            i0 = y0(k);
            rs = x1(k) + x2(k);
            is = y1(k) + y2(k);
            x0(k) = r0 + rs;
            y0(k) = i0 + is;
            ra = rs*a + r0;
            ia = is*a + i0;
            rb = (x1(k) - x2(k))*b;
            ib = (y1(k) - y2(k))*b;
            if ( zero ) {
              x1(k) = ra + ib;
              y1(k) = ia - rb;
              x2(k) = ra - ib;
              y2(k) = ia + rb;
            } else {
              r1 = ra + ib;
              i1 = ia - rb;
              r2 = ra - ib;
              i2 = ia + rb;
              x1(k) = r1*c1 + i1*s1;
              y1(k) = i1*c1 - r1*s1;
              x2(k) = r2*c2 + i2*s2;
              y2(k) = i2*c2 - r2*s2;
            }
          }
        }
      }
      if ( fold ) {
        fold = false;
        k0 = (m + 1 - j)*sep + 1;
        t = c1*a + s1*b;
        s1 = c1*b - s1*a;
        c1 = t;
        t = c2*a - s2*b;
        s2 = -c2*b - s2*a;
        c2 = t;
        continue;
      }
      break;
    }
  }

}

void
r4cftk(
  int const n,
  int & m,
  raw_ref1<float> x0,
  raw_ref1<float> y0,
  raw_ref1<float> x1,
  raw_ref1<float> y1,
  raw_ref1<float> x2,
  raw_ref1<float> y2,
  raw_ref1<float> x3,
  raw_ref1<float> y3,
  int5a dim
)
{
//     ==============================================================
//
//---- Radix 4 multi-dimensional complex fourier transform kernel

//     .. Local Scalars ..
  static float const twopi = 6.283185f;
  float angle, c1, c2, c3, fm4, i1, i2, i3, is0, is1, iu0, iu1, r1, r2, r3, rs0, rs1, ru0, ru1, s1, s2, s3, t;
  int k0, k1, k2, l1, m4, mm4, mover2, ns, nt, sep, size;
  bool fold, zero;

  nt = dim(1);
  sep = dim(2);
  l1 = dim(3);
  size = dim(4) - 1;
  k2 = dim(5);
  ns = n*sep;
  m4 = m*4;
  fm4 = float(m4);
  mm4 = sep*m4;
  mover2 = m/2 + 1;

  // avoid "may be used uninitialized" warnings
  c1 = s1 = c2 = s2 = c3 = s3 = 0.f;
  for ( int j = 1; j <= mover2; ++j ) {
    fold = j > 1 && 2*j < m + 2;
    k0 = (j - 1)*sep + 1;
    zero = j == 1;
    if ( ! zero ) {
      angle = twopi*float(j - 1)/fm4;
      c1 = std::cos(angle);
      s1 = std::sin(angle);
      c2 = c1*c1 - s1*s1;
      s2 = s1*c1 + c1*s1;
      c3 = c2*c1 - s2*s1;
      s3 = s2*c1 + c2*s1;
    }
    while ( true ) {
      for ( int kk = k0, smm4 = sign( mm4 ); kk * smm4 <= smm4 * ns; kk += mm4 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
        for ( int l = kk, sl1 = sign( l1 ); l * sl1 <= sl1 * nt; l += l1 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
          k1 = l + size;
          for ( int k = l, sk2 = sign( k2 ); k * sk2 <= sk2 * k1; k += k2 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
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
            if ( zero ) {
              x2(k) = ru0 + iu1;
              y2(k) = iu0 - ru1;
              x1(k) = rs0 - rs1;
              y1(k) = is0 - is1;
              x3(k) = ru0 - iu1;
              y3(k) = iu0 + ru1;
            } else {
              r1 = ru0 + iu1;
              i1 = iu0 - ru1;
              r2 = rs0 - rs1;
              i2 = is0 - is1;
              r3 = ru0 - iu1;
              i3 = iu0 + ru1;
              x2(k) = r1*c1 + i1*s1;
              y2(k) = i1*c1 - r1*s1;
              x1(k) = r2*c2 + i2*s2;
              y1(k) = i2*c2 - r2*s2;
              x3(k) = r3*c3 + i3*s3;
              y3(k) = i3*c3 - r3*s3;
            }
          }
        }
      }
      if ( fold ) {
        fold = false;
        k0 = (m + 1 - j)*sep + 1;
        t = c1;
        c1 = s1;
        s1 = t;
        c2 = -c2;
        t = c3;
        c3 = -s3;
        s3 = -t;
        continue;
      }
      break;
    }
  }

}

void
r5cftk(
  int const n,
  int & m,
  raw_ref1<float> x0,
  raw_ref1<float> y0,
  raw_ref1<float> x1,
  raw_ref1<float> y1,
  raw_ref1<float> x2,
  raw_ref1<float> y2,
  raw_ref1<float> x3,
  raw_ref1<float> y3,
  raw_ref1<float> x4,
  raw_ref1<float> y4,
  int5a dim
)
{
//     =================================================================
//
//---- Radix 5 multi-dimensional complex fourier transform kernel

//     .. Local Scalars ..
  static float const twopi = 6.283185f;
  static float const a1 = 0.30901699f;
  static float const b1 = 0.95105652f;
  static float const a2 = -0.80901699f;
  static float const b2 = 0.58778525f;
  float angle, c1, c2, c3, c4, fm5, i0, i1, i2, i3, i4, ia1, ia2, ib1, ib2, is1, is2, iu1, iu2, r0, r1, r2, r3, r4, ra1, ra2, rb1, rb2, rs1, rs2, ru1, ru2, s1, s2, s3, s4, t;
  int k0, k1, k2, l1, m5, mm5, mover2, ns, nt, sep, size;
  bool fold, zero;

  nt = dim(1);
  sep = dim(2);
  l1 = dim(3);
  size = dim(4) - 1;
  k2 = dim(5);
  ns = n*sep;
  m5 = m*5;
  fm5 = float(m5);
  mm5 = sep*m5;
  mover2 = m/2 + 1;

  // avoid "may be used uninitialized" warnings
  c1 = s1 = c2 = s2 = c3 = s3 = c4 = s4 = 0.f;
  for ( int j = 1; j <= mover2; ++j ) {
    fold = j > 1 && 2*j < m + 2;
    k0 = (j - 1)*sep + 1;
    zero = j == 1;
    if ( ! zero ) {
      angle = twopi*float(j - 1)/fm5;
      c1 = std::cos(angle);
      s1 = std::sin(angle);
      c2 = c1*c1 - s1*s1;
      s2 = s1*c1 + c1*s1;
      c3 = c2*c1 - s2*s1;
      s3 = s2*c1 + c2*s1;
      c4 = c2*c2 - s2*s2;
      s4 = s2*c2 + c2*s2;
    }
    while ( true ) {
      for ( int kk = k0, smm5 = sign( mm5 ); kk * smm5 <= smm5 * ns; kk += mm5 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
        for ( int l = kk, sl1 = sign( l1 ); l * sl1 <= sl1 * nt; l += l1 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
          k1 = l + size;
          for ( int k = l, sk2 = sign( k2 ); k * sk2 <= sk2 * k1; k += k2 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
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
            ra1 = rs1*a1 + r0 + rs2*a2;
            ia1 = is1*a1 + i0 + is2*a2;
            ra2 = rs1*a2 + r0 + rs2*a1;
            ia2 = is1*a2 + i0 + is2*a1;
            rb1 = ru1*b1 + ru2*b2;
            ib1 = iu1*b1 + iu2*b2;
            rb2 = ru1*b2 - ru2*b1;
            ib2 = iu1*b2 - iu2*b1;
            if ( zero ) {
              x1(k) = ra1 + ib1;
              y1(k) = ia1 - rb1;
              x2(k) = ra2 + ib2;
              y2(k) = ia2 - rb2;
              x3(k) = ra2 - ib2;
              y3(k) = ia2 + rb2;
              x4(k) = ra1 - ib1;
              y4(k) = ia1 + rb1;
            } else {
              r1 = ra1 + ib1;
              i1 = ia1 - rb1;
              r2 = ra2 + ib2;
              i2 = ia2 - rb2;
              r3 = ra2 - ib2;
              i3 = ia2 + rb2;
              r4 = ra1 - ib1;
              i4 = ia1 + rb1;
              x1(k) = r1*c1 + i1*s1;
              y1(k) = i1*c1 - r1*s1;
              x2(k) = r2*c2 + i2*s2;
              y2(k) = i2*c2 - r2*s2;
              x3(k) = r3*c3 + i3*s3;
              y3(k) = i3*c3 - r3*s3;
              x4(k) = r4*c4 + i4*s4;
              y4(k) = i4*c4 - r4*s4;
            }
          }
        }
      }
      if ( fold ) {
        fold = false;
        k0 = (m + 1 - j)*sep + 1;
        t = c1*a1 + s1*b1;
        s1 = c1*b1 - s1*a1;
        c1 = t;
        t = c2*a2 + s2*b2;
        s2 = c2*b2 - s2*a2;
        c2 = t;
        t = c3*a2 - s3*b2;
        s3 = -c3*b2 - s3*a2;
        c3 = t;
        t = c4*a1 - s4*b1;
        s4 = -c4*b1 - s4*a1;
        c4 = t;
        continue;
      }
      break;
    }
  }

}

void
r8cftk(
  int const n,
  int & m,
  raw_ref1<float> x0,
  raw_ref1<float> y0,
  raw_ref1<float> x1,
  raw_ref1<float> y1,
  raw_ref1<float> x2,
  raw_ref1<float> y2,
  raw_ref1<float> x3,
  raw_ref1<float> y3,
  raw_ref1<float> x4,
  raw_ref1<float> y4,
  raw_ref1<float> x5,
  raw_ref1<float> y5,
  raw_ref1<float> x6,
  raw_ref1<float> y6,
  raw_ref1<float> x7,
  raw_ref1<float> y7,
  int5a dim
)
{
//      ===============================================================
//
//---- Radix 8 multi-dimensional complex fourier transform kernel

//     .. Local Scalars ..
  static float const twopi = 6.283185f;
  static float const e = 0.70710678f;
  float angle, c1, c2, c3, c4, c5, c6, c7, fm8, i1, i2, i3, i4, i5, i6, i7, is0, is1, is2, is3, iss0, iss1, isu0, isu1, iu0, iu1, iu2, iu3, ius0, ius1, iuu0, iuu1, r1, r2, r3, r4, r5, r6, r7, rs0, rs1, rs2, rs3, rss0, rss1, rsu0, rsu1, ru0, ru1, ru2, ru3, rus0, rus1, ruu0, ruu1, s1, s2, s3, s4, s5, s6, s7, t;
  int k0, k1, k2, l1, m8, mm8, mover2, ns, nt, sep, size;
  bool fold, zero;

  nt = dim(1);
  sep = dim(2);
  l1 = dim(3);
  size = dim(4) - 1;
  k2 = dim(5);
  ns = n*sep;
  m8 = m*8;
  fm8 = float(m8);
  mm8 = sep*m8;
  mover2 = m/2 + 1;

  // avoid "may be used uninitialized" warnings
  c1 = s1 = c2 = s2 = c3 = s3 = c4 = s4 = c5 = s5 = c6 = s6 = c7 = s7 = 0.f;
  for ( int j = 1; j <= mover2; ++j ) {
    fold = j > 1 && 2*j < m + 2;
    k0 = (j - 1)*sep + 1;
    zero = j == 1;
    if ( ! zero ) {
      angle = twopi*float(j - 1)/fm8;
      c1 = std::cos(angle);
      s1 = std::sin(angle);
      c2 = c1*c1 - s1*s1;
      s2 = s1*c1 + c1*s1;
      c3 = c2*c1 - s2*s1;
      s3 = s2*c1 + c2*s1;
      c4 = c2*c2 - s2*s2;
      s4 = s2*c2 + c2*s2;
      c5 = c4*c1 - s4*s1;
      s5 = s4*c1 + c4*s1;
      c6 = c4*c2 - s4*s2;
      s6 = s4*c2 + c4*s2;
      c7 = c4*c3 - s4*s3;
      s7 = s4*c3 + c4*s3;
    }
    while ( true ) {
      for ( int kk = k0, smm8 = sign( mm8 ); kk * smm8 <= smm8 * ns; kk += mm8 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
        for ( int l = kk, sl1 = sign( l1 ); l * sl1 <= sl1 * nt; l += l1 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
          k1 = l + size;
          for ( int k = l, sk2 = sign( k2 ); k * sk2 <= sk2 * k1; k += k2 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
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
            t = (rus1 + ius1)*e;
            ius1 = (ius1 - rus1)*e;
            rus1 = t;
            t = (ruu1 + iuu1)*e;
            iuu1 = (iuu1 - ruu1)*e;
            ruu1 = t;
            x0(k) = rss0 + rss1;
            y0(k) = iss0 + iss1;
            if ( zero ) {
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
            } else {
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
              x4(k) = r1*c1 + i1*s1;
              y4(k) = i1*c1 - r1*s1;
              x2(k) = r2*c2 + i2*s2;
              y2(k) = i2*c2 - r2*s2;
              x6(k) = r3*c3 + i3*s3;
              y6(k) = i3*c3 - r3*s3;
              x1(k) = r4*c4 + i4*s4;
              y1(k) = i4*c4 - r4*s4;
              x5(k) = r5*c5 + i5*s5;
              y5(k) = i5*c5 - r5*s5;
              x3(k) = r6*c6 + i6*s6;
              y3(k) = i6*c6 - r6*s6;
              x7(k) = r7*c7 + i7*s7;
              y7(k) = i7*c7 - r7*s7;
            }
          }
        }
      }
      if ( fold ) {
        fold = false;
        k0 = (m + 1 - j)*sep + 1;
        t = (c1 + s1)*e;
        s1 = (c1 - s1)*e;
        c1 = t;
        t = s2;
        s2 = c2;
        c2 = t;
        t = (-c3 + s3)*e;
        s3 = (c3 + s3)*e;
        c3 = t;
        c4 = -c4;
        t = -(c5 + s5)*e;
        s5 = (-c5 + s5)*e;
        c5 = t;
        t = -s6;
        s6 = -c6;
        c6 = t;
        t = (c7 - s7)*e;
        s7 = -(c7 + s7)*e;
        c7 = t;
        continue;
      }
      break;
    }
  }

}

void
rpcftk(
  int const n,
  int & m,
  int & p,
  int & r,
  raw_ref2<float> const& x,
  raw_ref2<float> const& y,
  int5a dim
)
{
//     ==========================================
//
//---- Radix prime multi-dimensional complex fourier transform kernel
//
//MDW: Note, routine works beyond the nominal bounds of X and Y.
//     We think this is deliberate, so don't be tempted to "fix" it.
//     This routine is therefore incompatible with "bounds-checking"
//     options of compilers.

//     .. Local Scalars ..
  static float const twopi = 6.283185f;
  float angle, fmp, fp, fu, is, iu, rs, ru, t, xt,yt;
  int jj, k0, k1, k2, l1, mmp, mover2, mp, ns, nt, pm, pp, sep, size;
  bool fold, zero;

//     .. Local Arrays ..
  tiny_array<float, 18> a;
  tiny_matrix<float, 9, 9> aa;
  tiny_array<float, 18> b;
  tiny_matrix<float, 9, 9> bb;
  tiny_array<float, 18> c;
  tiny_array<float, 9> ia;
  tiny_array<float, 9> ib;
  tiny_array<float, 9> ra;
  tiny_array<float, 9> rb;
  tiny_array<float, 18> s;

  nt = dim(1);
  sep = dim(2);
  l1 = dim(3);
  size = dim(4) - 1;
  k2 = dim(5);
  ns = n*sep;
  mover2 = m/2 + 1;
  mp = m*p;
  fmp = float(mp);
  mmp = sep*mp;
  pp = p/2;
  pm = p - 1;
  fp = float(p);
  fu = 0.0f;
  for ( int u = 1; u <= pp; ++u ) {
    fu += 1.0f;
    angle = twopi*fu/fp;
    jj = p - u;
    a(u) = std::cos(angle);
    b(u) = std::sin(angle);
    a(jj) = a(u);
    b(jj) = -b(u);
  }
  for ( int u = 1; u <= pp; ++u ) {
    for ( int v = 1; v <= pp; ++v ) {
      jj = u*v - u*v/p*p;
      aa(v,u) = a(jj);
      bb(v,u) = b(jj);
    }
  }

  for ( int j = 1; j <= mover2; ++j ) {
    fold = j > 1 && 2*j < m + 2;
    k0 = (j - 1)*sep + 1;
    zero = j == 1;
    if ( ! zero ) {
      angle = twopi*float(j - 1)/fmp;
      c(1) = std::cos(angle);
      s(1) = std::sin(angle);
      for ( int u = 2; u <= pm; ++u ) {
        c(u) = c(u-1)*c(1) - s(u-1)*s(1);
        s(u) = s(u-1)*c(1) + c(u-1)*s(1);
      }
    }
    while ( true ) {
      for ( int kk = k0, smmp = sign( mmp ); kk * smmp <= smmp * ns; kk += mmp ) { //ObjexxF2Cxx Remove sign() and set <= or >=
        for ( int l = kk, sl1 = sign( l1 ); l * sl1 <= sl1 * nt; l += l1 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
          k1 = l + size;
          for ( int k = l, sk2 = sign( k2 ); k * sk2 <= sk2 * k1; k += k2 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
            xt = x(k,1);
            yt = y(k,1);
            rs = x(k,2) + x(k,p);
            is = y(k,2) + y(k,p);
            ru = x(k,2) - x(k,p);
            iu = y(k,2) - y(k,p);
            for ( int u = 1; u <= pp; ++u ) {
              ra(u) = aa(u,1)*rs + xt;
              ia(u) = aa(u,1)*is + yt;
              rb(u) = bb(u,1)*ru;
              ib(u) = bb(u,1)*iu;
            }
            xt += rs;
            yt += is;
            for ( int u = 2; u <= pp; ++u ) {
              jj = p - u;
              rs = x(k,u+1) + x(k,jj+1);
              is = y(k,u+1) + y(k,jj+1);
              ru = x(k,u+1) - x(k,jj+1);
              iu = y(k,u+1) - y(k,jj+1);
              xt += rs;
              yt += is;
              for ( int v = 1; v <= pp; ++v ) {
                ra(v) = aa(v,u)*rs + ra(v);
                ia(v) = aa(v,u)*is + ia(v);
                rb(v) = bb(v,u)*ru + rb(v);
                ib(v) = bb(v,u)*iu + ib(v);
              }
            }
            x(k,1) = xt;
            y(k,1) = yt;
            for ( int u = 1; u <= pp; ++u ) {
              jj = p - u;
              if ( zero ) {
                x(k,u+1) = ra(u) + ib(u);
                y(k,u+1) = ia(u) - rb(u);
                x(k,jj+1) = ra(u) - ib(u);
                y(k,jj+1) = ia(u) + rb(u);
              } else {
                xt = ra(u) + ib(u);
                yt = ia(u) - rb(u);
                x(k,u+1) = c(u)*xt + s(u)*yt;
                y(k,u+1) = c(u)*yt - s(u)*xt;
                xt = ra(u) - ib(u);
                yt = ia(u) + rb(u);
                x(k,jj+1) = c(jj)*xt + s(jj)*yt;
                y(k,jj+1) = c(jj)*yt - s(jj)*xt;
              }
            }
          }
        }
      }
      if ( fold ) {
        fold = false;
        k0 = (m + 1 - j)*sep + 1;
        for ( int u = 1; u <= pm; ++u ) {
          t = c(u)*a(u) + s(u)*b(u);
          s(u) = -s(u)*a(u) + c(u)*b(u);
          c(u) = t;
        }
        continue;
      }
      break;
    }
  }

}

void
mdftkd(
  int const n,
  tiny_ref<int, 15> factor,
  int5a dim,
  raw_ref1<float> x,
  raw_ref1<float> y
)
{
  using std::cout;
//     ==========================================
//
//---- Multi-dimensional complex fourier transform kernel driver

//     .. Local Scalars ..
  int f, m, p, r, s;

  s = dim(2);
  f = 0;
  m = n;
  while ( true ) {
    ++f;
    p = factor(f);
    if ( p == 0 ) {
      return;
    } else {
      m /= p;
      r = m*s;
      if ( p <= 8 ) {
        switch (p) {
        case 1 :
          break;
        case 2 :
          r2cftk(n,m,x(1),y(1),x(r+1),y(r+1),dim);
          break;
        case 3 :
          r3cftk(n,m,x(1),y(1),x(r+1),y(r+1),x(2*r+1),y(2*r+1),dim);
          break;
        case 4 :
          r4cftk(
            n,m,x(1),y(1),x(r+1),y(r+1),x(2*r+1),y(2*r+1),x(3*r+1),y(3*r+1),
            dim);
          break;
        case 5 :
          r5cftk(
            n,m,x(1),y(1),x(r+1),y(r+1),x(2*r+1),y(2*r+1),x(3*r+1),y(3*r+1),
            x(4*r+1),y(4*r+1),dim);
          break;
        case 7 :
          goto L70;
        case 8 :
          r8cftk(
            n,m,x(1),y(1),x(r+1),y(r+1),x(2*r+1),y(2*r+1),x(3*r+1),y(3*r+1),
            x(4*r+1),y(4*r+1),x(5*r+1),y(5*r+1),x(6*r+1),y(6*r+1),x(7*r+1),
            y(7*r+1),dim);
          break;
        default :
          cout << "TRANSFER ERROR DETECTED IN MDFTKD\n\n\n";
          goto L71;
        }
      } else { // p > 8
L70:
        rpcftk(n,m,p,r,
          raw_ref2<float>(&x(1), r, p),
          raw_ref2<float>(&y(1), r, p),
          dim);
      }
    }
  }
L71:;
}

void
diprp(
  int const pts,
  tiny_ref<int, 15> sym,
  int const psym,
  tiny_ref<int, 15> unsym,
  int5a dim,
  raw_ref1<float> x,
  raw_ref1<float> y
)
{
//     =================================================
//
//---- Double in place reordering programme

//     .. Local Scalars ..
  int delta, dk, kk, lk, mods, mult, nest, nt, p1, p2, p3, p4, p5, punsym, sep, size, test;
  bool onemod;

//     .. Local Arrays ..
  tiny_array<int, 14> modulo;
  tiny_array<int, 14> s;
  tiny_array<int, 14> u;

//     .. Equivalences ..

  int & bs( s( 2 ) );
  int & cs( s( 3 ) );
  int & ds( s( 4 ) );
  int & es( s( 5 ) );
  int & fs( s( 6 ) );
  int & gs( s( 7 ) );
  int & hs( s( 8 ) );
  int & is( s( 9 ) );
  int & js( s( 10 ) );
  int & ks( s( 11 ) );
  int & ls( s( 12 ) );
  int & ms( s( 13 ) );
  int & ns( s( 14 ) );

  int & al( u( 1 ) );
  int & bl( u( 2 ) );
  int & cl( u( 3 ) );
  int & dl( u( 4 ) );
  int & el( u( 5 ) );
  int & fl( u( 6 ) );
  int & gl( u( 7 ) );
  int & hl( u( 8 ) );
  int & il( u( 9 ) );
  int & jl( u( 10 ) );
  int & kl( u( 11 ) );
  int & ll( u( 12 ) );
  int & ml( u( 13 ) );
  int & nl( u( 14 ) );

  nest = 14;
  nt = dim(1);
  sep = dim(2);
  p2 = dim(3);
  size = dim(4) - 1;
  p4 = dim(5);
  assert( p2 > 0 );
  assert( p4 > 0 );
  if ( sym(1) != 0 ) {
    u.fill(1);
    s.fill(1);
    int np = pts;
    for ( int j = 1; j <= nest; ++j ) {
      if ( sym(j) == 0 ) {
        break;
      }
      int const jj = nest + 1 - j;
      u(jj) = np;
      s(jj) = np/sym(j);
      np /= sym(j);
    }
    assert( bs > 0 );
    assert( cs > 0 );
    assert( ds > 0 );
    assert( es > 0 );
    assert( fs > 0 );
    assert( gs > 0 );
    assert( hs > 0 );
    assert( is > 0 );
    assert( js > 0 );
    assert( ks > 0 );
    assert( ls > 0 );
    assert( ms > 0 );
    assert( ns > 0 );
    int jj = 0;
    for ( int a = 1; a <= al; ++a ) {
      for ( int b = a; b <= bl; b += bs ) {
        for ( int c = b; c <= cl; c += cs ) {
          for ( int d = c; d <= dl; d += ds ) {
            for ( int e = d; e <= el; e += es ) {
              for ( int f = e; f <= fl; f += fs ) {
                for ( int g = f; g <= gl; g += gs ) {
                  for ( int h = g; h <= hl; h += hs ) {
                    for ( int i = h; i <= il; i += is ) {
                      for ( int j = i; j <= jl; j += js ) {
                        for ( int k = j; k <= kl; k += ks ) {
                          for ( int l = k; l <= ll; l += ls ) {
                            for ( int m = l; m <= ml; m += ms ) {
                              for ( int n = m; n <= nl; n += ns ) {
                                ++jj;
                                if ( jj < n ) {
                                  delta = ( n - jj ) * sep;
                                  p1 = ( jj - 1 ) * sep + 1;
                                  for ( int p0 = p1; p0 <= nt; p0 += p2 ) {
                                    p3 = p0 + size;
                                    for ( int p = p0; p <= p3; p += p4 ) {
                                      p5 = p + delta;
                                      std::swap( x(p), x(p5) );
                                      std::swap( y(p), y(p5) );
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

  if ( unsym(1) != 0 ) {
    punsym = pts/square( psym );
    mult = punsym/unsym(1);
    test = (unsym(1)*unsym(2) - 1)*mult*psym;
    lk = mult;
    dk = mult;
    mods = -1; // avoid "may be used uninitialized" warning
    for ( int k = 2; k <= nest; ++k ) {
      if ( unsym(k) == 0 ) {
        break;
      }
      lk *= unsym(k - 1);
      dk /= unsym(k);
      u(k) = (lk - dk)*psym;
      mods = k;
    }
    ASSERTBX(mods >= 2);
    onemod = mods < 3;
    if ( ! onemod ) {
      for ( int j = 3; j <= mods; ++j ) {
        int const jj = mods + 3 - j;
        modulo(jj) = u(j);
      }
    }
    modulo(2) = u(2);
    jl = ( punsym - 3 ) * psym;
    ms = punsym * psym;
    for ( int j = psym, spsym = sign( psym ); j * spsym <= spsym * jl; j += psym ) { //ObjexxF2Cxx Remove sign() and set <= or >=
      int k = j;
      while ( true ) {
        k *= mult;
        if ( ! onemod ) {
          for ( int i = 3; i <= mods; ++i ) {
            k = k - (k/modulo(i))*modulo(i);
          }
        }
        if ( k >= test ) {
          k = k - (k/modulo(2))*modulo(2) + modulo(2);
        } else {
          k = k - (k/modulo(2))*modulo(2);
        }
        if ( k < j ) {
        } else {
          if ( k != j ) {
            delta = (k - j)*sep;
            for ( int l = 1; l <= psym; ++l ) {
              for ( int m = l, sms = sign( ms ); m * sms <= sms * pts; m += ms ) { //ObjexxF2Cxx Remove sign() and set <= or >=
                p1 = (m + j - 1)*sep + 1;
                for ( int p0 = p1, sp2 = sign( p2 ); p0 * sp2 <= sp2 * nt; p0 += p2 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
                  p3 = p0 + size;
                  for ( int jj = p0, sp4 = sign( p4 ); jj * sp4 <= sp4 * p3; jj += p4 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
                    kk = jj + delta;
                    std::swap( x(jj), x(kk) );
                    std::swap( y(jj), y(kk) );
                  }
                }
              }
            }
          }
          break;
        }
      }
    }
  }

}

} // namespace <anonymous>

void
cmplft(
  FArray1A_float x,
  FArray1A_float y,
  int const n,
  int5a d
)
{
  using std::cerr;
  x.dim( star );
  y.dim( star );
//     ===============================
//
//     Complex finite discrete fourier transform
//     transforms one dimension of multi-dimensional data
//     modified by L. F. TEN EYCK from a one-dimensional version written
//     by G. T. SANDE, 1969.f
//
//     This program calculates the transform
//               (X(T) + I*Y(T))*(cos(2*PI*T/N) - I*sin(2*PI*T/N))
//
//     INDEXING -- the arrangement of the multi-dimensional data is
//     specified by the integer array D, the values of which are used as
//     control parameters in do loops.  when it is desired to cover all
//     elements of the data for which the subscript being transformed has
//     the value I0, the following is used.
//
//               I1 = (I0 - 1)*D(2) + 1
//               DO 100 I2 = I1, D(1), D(3)
//               I3 = I2 + D(4) - 1
//               DO 100 I = I2, I3, D(5)
//                  .
//                  .
//           100 CONTINUE
//
//     with this indexing it is possible to use a number of arrangements
//     of the data, including normal fortran complex numbers (d(5) = 2)
//     or separate storage of real and imaginary parts.
//
//
//---- PMAX is the largest prime factor that will be tolerated by this
//     program.
//
//---- TWOGRP is the largest power of two that is treated as a special
//     case.

//     .. Local Scalars ..

//     .. Local Arrays ..
  tiny_array<int, 15> factor;
  tiny_array<int, 15> sym;
  tiny_array<int, 15> unsym;

  int const pmax = 19;
  int const twogrp = 8;
  int psym;

  if ( n > 1 ) {
    bool error;
    srfp(n,pmax,twogrp,factor,sym,psym,unsym,error);
    if ( error ) {
      std::ostringstream stream_emess;
      stream_emess << /* FORMAT("FFTLIB: invalid number of points for CMPL FT.  N =",i10) */ boost::format( "FFTLIB: invalid number of points for CMPL FT.  N =%10d\n" ) % n;
      cerr << stream_emess.str() << std::endl; //Objexx Was ccperr(1,stream_emess.str());
    } else {
      assert( d(1) >= 0 );
      assert( d(2) >= 0 );
      assert( d(3) >= 0 );
      assert( d(4) >= 0 );
      assert( d(5) >= 0 );
      mdftkd(n,factor,d,x(1),y(1));
      diprp(n,sym,psym,unsym,d,x(1),y(1));
    }
  }

}

void
hermft(
  FArray1A_float x,
  FArray1A_float y,
  int const n,
  int5a dim
)
{
  x.dim( star );
  y.dim( star );

//     ==================================
//
//---- Hermitian symmetric fourier transform
//
//     Given the unique terms of a hermitian symmetric sequence of length
//     2N this subroutine calculates the 2N real numbers which are its
//     fourier transform.  The even numbered elements of the transform
//     (0, 2, 4, . . ., 2n-2) are returned in X and the odd numbered
//     elements (1, 3, 5, . . ., 2n-1) in Y.
//
//     A finite hermitian sequence of length 2n contains n + 1 unique
//     real numbers and n - 1 unique imaginary numbers.  For convenience
//     the real value for X(n) is stored at Y(0).

//     .. Local Scalars ..
  float a, angle, b, c, co, d, e, f, si;
  int d2, d3, d4, d5, i2, j, k, k1, nover2, nt;

  static float const twopi = 6.283185f;
  float const twon = float(2*n);

  nt = dim(1);
  d2 = dim(2);
  d3 = dim(3);
  d4 = dim(4) - 1;
  d5 = dim(5);

  for ( int i0 = 1, sd3 = sign( d3 ); i0 * sd3 <= sd3 * nt; i0 += d3 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
    int const i1 = i0 + d4;
    for ( int i = i0, sd5 = sign( d5 ); i * sd5 <= sd5 * i1; i += d5 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
      a = x(i);
      b = y(i);
      x(i) = a + b;
      y(i) = a - b;
    }
  }

  nover2 = n/2 + 1;
  if ( nover2 >= 2 ) {
    for ( int i0 = 2; i0 <= nover2; ++i0 ) {
      angle = float(i0 - 1)*twopi/twon;
      co = std::cos(angle);
      si = std::sin(angle);
      k = (n + 2 - 2*i0)*d2;
      k1 = (i0 - 1)*d2 + 1;
      for ( int i1 = k1, sd3 = sign( d3 ); i1 * sd3 <= sd3 * nt; i1 += d3 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
        i2 = i1 + d4;
        for ( int i = i1, sd5 = sign( d5 ); i * sd5 <= sd5 * i2; i += d5 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
          j = i + k;
          a = x(i) + x(j);
          b = x(i) - x(j);
          c = y(i) + y(j);
          d = y(i) - y(j);
          e = b*co + c*si;
          f = b*si - c*co;
          x(i) = a + f;
          x(j) = a - f;
          y(i) = e + d;
          y(j) = e - d;
        }
      }
    }
    cmplft(x,y,n,dim);
  }

}

void
realft(
  FArray1A_float even,
  FArray1A_float odd,
  int const n,
  int5a dim
)
{
  even.dim( star );
  odd.dim( star );

//      ======================================
//
//     REAL FOURIER TRANSFORM
//
//     Given a real sequence of length 2n this subroutine calculates the
//     unique part of the fourier transform.  The fourier transform has
//     n + 1 unique real parts and n - 1 unique imaginary parts.  Since
//     the real part at x(n) is frequently of interest, this subroutine
//     stores it at x(n) rather than in y(0).  Therefore x and y must be
//     of length n + 1 instead of n.  Note that this storage arrangement
//     is different from that employed by the hermitian fourier transform
//     subroutine.
//
//     For convenience the data is presented in two parts, the first
//     containing the even numbered real terms and the second containing
//     the odd numbered terms (numbering starting at 0).  On return the
//     real part of the transform replaces the even terms and the
//     imaginary part of the transform replaces the odd terms.

//     .. Local Scalars ..
  float a, angle, b, c, co, d, e, f, si;
  int d2, d3, d4, d5, i0, i2, j, l, nover2, nt;

  static float const twopi = 6.283185f;
  float const twon = float(2*n);

  cmplft(even,odd,n,dim);

  nt = dim(1);
  d2 = dim(2);
  d3 = dim(3);
  d4 = dim(4) - 1;
  d5 = dim(5);
  nover2 = n/2 + 1;

  if ( nover2 >= 2 ) {
    for ( int i = 2; i <= nover2; ++i ) {
      angle = float(i - 1)*twopi/twon;
      co = std::cos(angle);
      si = std::sin(angle);
      i0 = (i - 1)*d2 + 1;
      j = (n + 2 - 2*i)*d2;
      for ( int i1 = i0, sd3 = sign( d3 ); i1 * sd3 <= sd3 * nt; i1 += d3 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
        i2 = i1 + d4;
        for ( int k = i1, sd5 = sign( d5 ); k * sd5 <= sd5 * i2; k += d5 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
          l = k + j;
          a = (even(l) + even(k))/2.0f;
          c = (even(l) - even(k))/2.0f;
          b = (odd(l) + odd(k))/2.0f;
          d = (odd(l) - odd(k))/2.0f;
          e = c*si + b*co;
          f = c*co - b*si;
          even(k) = a + e;
          even(l) = a - e;
          odd(k) = f - d;
          odd(l) = f + d;
        }
      }
    }
  }

  if ( n >= 1 ) {
    j = n*d2;
    for ( int i1 = 1, sd3 = sign( d3 ); i1 * sd3 <= sd3 * nt; i1 += d3 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
      i2 = i1 + d4;
      for ( int k = i1, sd5 = sign( d5 ); k * sd5 <= sd5 * i2; k += d5 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
        l = k + j;
        even(l) = even(k) - odd(k);
        odd(l) = 0.0f;
        even(k) += odd(k);
        odd(k) = 0.0f;
      }
    }
  }

}

void
rsymft(
  FArray1A_float x,
  int const n,
  int5a dim
)
{
  using std::cerr;
  x.dim( star );

//      ===============================
//
//     REAL SYMMETRIC MULTIDIMENSIONAL FOURIER TRANSFORM
//
//     N must be a multiple of 4.f  The two unique elements are stored at
//     X(1) and X(n+1).
//
//     Decimation in frequency applied to a real symmetric sequence of
//     length 2n gives a real symmetric sequence of length n, the
//     transform of which gives the even numbered fourier coefficients,
//     and a hermitian symmetric sequence of length n, the transform of
//     which gives the odd numbered fourier coefficients.  The sum of
//     the two sequences is a hermitian symmetric sequence of length n,
//     which may be stored in n/2 complex locations.  The transform of
//     this sequence is n real numbers representing the term by term sum
//     of the even and odd numbered fourier coefficients.  This symmetric
//     sequence may be solved if any of the fourier coefficients are
//     known.  For this purpose x0, which is simply the sum of the
//     original sequence, is computed and saved in x(n+1).

//     .. Local Scalars ..
  float a, angle, b, c, co, d, si, twon, twopi;
  int d1, d2, d3, d4, d5, i0, i2, j0, j1, k0, k2, l, m, mj, mk, ml, mm, nn, nover2, nover4, twod2;

  if ( n != 1 ) {
    nover2 = n/2;
    nover4 = n/4;
    if ( 4*nover4 != n ) {
      std::ostringstream stream_emess;
      stream_emess << /* FORMAT("FFTLIB: N not a multiple of 4 in R SYM FT.  N =",i10,//) */ boost::format( "FFTLIB: N not a multiple of 4 in R SYM FT.  N =%10d\n\n\n" ) % n;
      cerr << stream_emess.str() << std::endl; //Objexx Was ccperr(1,stream_emess.str());
    } else {
      d1 = dim(1);
      d2 = dim(2);
      d3 = dim(3);
      d4 = dim(4) - 1;
      d5 = dim(5);
      twopi = 6.283185f;
      twon = float(2*n);
      twod2 = 2*d2;
      k0 = n*d2 + 1;
      for ( int k1 = k0, sd3 = sign( d3 ); k1 * sd3 <= sd3 * d1; k1 += d3 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
        k2 = k1 + d4;
        for ( int k = k1, sd5 = sign( d5 ); k * sd5 <= sd5 * k2; k += d5 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
          x(k) /= 2.0f;
        }
      }

      for ( int i = 2; i <= nover2; ++i ) {
        angle = float(i - 1)*twopi/twon;
        co = std::cos(angle);
        si = std::sin(angle);
        k0 = (i - 1)*d2 + 1;
        j0 = (n + 2 - 2*i)*d2;
        j1 = (n + 1 - i)*d2;
        for ( int k1 = k0, sd3 = sign( d3 ); k1 * sd3 <= sd3 * d1; k1 += d3 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
          k2 = k1 + d4;
          for ( int k = k1, sd5 = sign( d5 ); k * sd5 <= sd5 * k2; k += d5 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
            l = k + j0;
            nn = k + j1;
            a = x(l) + x(k);
            b = x(l) - x(k);
            x(k) = a - b*co;
            x(l) = b*si;
            x(nn) += a;
          }
        }
      }

      if ( nover4 != 1 ) {
        j0 = nover4 - 1;
        for ( int i = 1; i <= j0; ++i ) {
          k0 = (nover2 + i)*d2 + 1;
          j1 = (nover2 - 2*i)*d2;
          for ( int k1 = k0, sd3 = sign( d3 ); k1 * sd3 <= sd3 * d1; k1 += d3 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
            k2 = k1 + d4;
            for ( int k = k1, sd5 = sign( d5 ); k * sd5 <= sd5 * k2; k += d5 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
              l = k + j1;
              a = x(k);
              x(k) = x(l);
              x(l) = a;
            }
          }
        }
      }

      j0 = nover2*d2;
      j1 = n*d2;
      for ( int k1 = 1, sd3 = sign( d3 ); k1 * sd3 <= sd3 * d1; k1 += d3 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
        k2 = k1 + d4;
        for ( int k = k1, sd5 = sign( d5 ); k * sd5 <= sd5 * k2; k += d5 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
          int const i = k + j0;
          l = k + j1;
          x(i) *= 2.0f;
          x(l) = x(k) + x(i) + x(l)*2.0f;
          x(k) *= 2.0f;
        }
      }

      int const k = nover2*d2 + 1;
      hermft(x(1),x(k),nover2,dim);

//---- Solve the equations for all of the sequences
//
      i0 = 1 - d2;
      mk = nover2*d2;
      mj = mk + d2;
      ml = n*d2 + d2;
      mm = ml;
      for ( int ii = 1; ii <= nover4; ++ii ) {
        i0 += d2;
        mj -= twod2;
        ml -= twod2;
        mm -= d2;
        for ( int i1 = i0, sd3 = sign( d3 ); i1 * sd3 <= sd3 * d1; i1 += d3 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
          i2 = i1 + d4;
          for ( int i = i1, sd5 = sign( d5 ); i * sd5 <= sd5 * i2; i += d5 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
            int const j = i + mj;
            int const k = i + mk;
            l = i + ml;
            m = i + mm;
            a = x(i) - x(m);
            b = x(l) - a;
            c = x(k) - b;
            d = x(j) - c;
            x(i) = x(m);
            x(j) = a;
            x(k) = b;
            x(l) = c;
            x(m) = d;
          }
        }
      }

//---- The results are now in a scrambled digit reversed order, i.e.
//     x(1), x(5), x(9), ..., x(10), x(6), x(2), ..., x(3), x(7), x(11),
//     ..., x(12), x(8), x(4).  the following section of program follows
//     the permutation cycles and does the necessary interchanges.
//
      if ( nover4 != 1 ) {
        nn = n - 2;
        for ( int i = 1; i <= nn; ++i ) {
          int k = i;
          while ( true ) {
            k0 = k/4;
            l = k - k0*4;
            if ( l != (l/2)*2 ) k0 = nover4 - 1 - k0;
            k = l*nover4 + k0;
            if ( k < i ) {
            } else {
              if ( k != i ) {
                k0 = i*d2 + 1;
                j0 = (k - i)*d2;
                for ( int k1 = k0, sd3 = sign( d3 ); k1 * sd3 <= sd3 * d1; k1 += d3 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
                  k2 = k1 + d4;
                  for ( int k = k1, sd5 = sign( d5 ); k * sd5 <= sd5 * k2; k += d5 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
                    l = k + j0;
                    a = x(k);
                    x(k) = x(l);
                    x(l) = a;
                  }
                }
              }
              break;
            }
          }
        }
      }
    }
  }

//---- Format statements
//
//FMT9000  FORMAT ("FFTLIB: N not a multiple of 4 in R SYM FT.  N =",i10,//)
}

void
sdiad(
  FArray1A_float x,
  FArray1A_float y,
  int const n,
  int5a dim
)
{
  using std::cerr;
  x.dim( star );
  y.dim( star );
//      ===============================
//
//     This subroutine computes half the fourier synthesis along a screw
//     diad lying along a crystallographic axis given half the fourier
//     coefficients.  That is, it assumes that f(t) = std::conj(f(-t)) for t
//     even and f(t) = -std::conj(f(-t)) for t odd.  n is the length of the
//     desired half of the transform.  The location x(n+1) is required as
//     a scratch location and therefore a value is also returned in
//     x(n+1) and y(n+1).  The value of the second half of the transform
//     may be generated from the first half by the formula x(n+t) = x(t),
//     y(n+t) = -y(t).  In other words, the last half of the transform is
//     the complex conjugate of the first half.
//
//     The transform is calculated by forming the sum of the even terms
//     and the odd terms in place, using the symmetry relations to
//     obtain the values for negative subscripts.  The transform of the
//     resulting sequence may be separated by using the fact that the
//     transform of the even terms is real, while the prodct of the
//     transform of the odd terms and (std::cos(pi*t/n) - i*std::sin(pi*t/n)) is
//     imaginary.  The scratch location is required because the formula
//     for separating the two transforms breaks down when t = n/2.f
//
// Corrections from A.D.MCLACHLAN 1980, put here sep 1985
// errors in original algorithm for the scratch location which
// assumed f(n)=0

//     .. Local Scalars ..
  float a, angle, c, s, twon, twopi;
  int d1, d2, d3, d4, d5, j, k0, k2, l, m, mn, nn, nover2;
  bool fold;

  nover2 = n/2;
  if ( 2*nover2 != n ) {
    std::ostringstream stream_emess;
    stream_emess << /* FORMAT("FFT error: SDIAD: N odd.  N =",i10) */ boost::format( "FFT error: SDIAD: N odd.  N =%10d\n" ) % n;
    cerr << stream_emess.str() << std::endl; //Objexx Was ccperr(1,stream_emess.str());
  } else {
    twon = float(2*n);
    twopi = 6.2831852f;
    d1 = dim(1);
    d2 = dim(2);
    d3 = dim(3);
    d4 = dim(4) - 1;
    d5 = dim(5);

    s = -1.0f;
    if ( nover2 == (2*(nover2/2)) ) s = -s;
    k0 = (n - 1)*d2 + 1;
    for ( int k1 = k0, sd3 = sign( d3 ); k1 * sd3 <= sd3 * d1; k1 += d3 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
      k2 = k1 + d4;
      for ( int k = k1, sd5 = sign( d5 ); k * sd5 <= sd5 * k2; k += d5 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
        l = k + d2;
        y(l) = x(k)*s;
      }
    }
    s = 1.0f;
    nn = n - 2;
    for ( int i = 1; i <= nn; i += 2 ) {
      s = -s;
      mn = (n + 1 - i)*d2;
      k0 = (i - 1)*d2 + 1;
      for ( int k1 = k0, sd3 = sign( d3 ); k1 * sd3 <= sd3 * d1; k1 += d3 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
        k2 = k1 + d4;
        for ( int k = k1, sd5 = sign( d5 ); k * sd5 <= sd5 * k2; k += d5 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
          j = k + d2;
          l = 2*d2 + k;
          m = k + mn;
          y(m) = x(j)*s + y(m);
          x(k) += x(j);
          x(j) = x(l) - x(j);
          y(k) += y(j);
          y(j) -= y(l);
        }
      }
    }
    k0 = (n - 2)*d2 + 1;
    for ( int k1 = k0, sd3 = sign( d3 ); k1 * sd3 <= sd3 * d1; k1 += d3 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
      k2 = k1 + d4;
      for ( int k = k1, sd5 = sign( d5 ); k * sd5 <= sd5 * k2; k += d5 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
        l = k + d2;
        x(k) += x(l);
        y(k) += y(l);
        j = l + d2;
        x(l) = x(j) - x(l);
      }
    }

//---- Reorder scrambled fourier coefficients
//
    for ( int i = 1; i <= nn; ++i ) {
      int k = i;
      while ( true ) {
        k *= 2;
        if ( k > n - 1 ) k = 2*n - 1 - k;
        if ( k < i ) {
        } else {
          if ( k != i ) {
            j = (k - i)*d2;
            k0 = i*d2 + 1;
            for ( int k1 = k0, sd3 = sign( d3 ); k1 * sd3 <= sd3 * d1; k1 += d3 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
              k2 = k1 + d4;
              for ( int k = k1, sd5 = sign( d5 ); k * sd5 <= sd5 * k2; k += d5 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
                l = k + j;
                a = x(k);
                x(k) = x(l);
                x(l) = a;
                a = y(k);
                y(k) = y(l);
                y(l) = a;
              }
            }
          }
          break;
        }
      }
    }

    cmplft(x,y,n,dim);

    m = nover2 - 1;
    for ( int i = 1; i <= m; ++i ) {
      angle = float(i)*twopi/twon;
      c = std::cos(angle);
      s = std::sin(angle);
      k0 = i*d2 + 1;
      fold = true;
      while ( true ) {
        for ( int k1 = k0, sd3 = sign( d3 ); k1 * sd3 <= sd3 * d1; k1 += d3 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
          k2 = k1 + d4;
          for ( int k = k1, sd5 = sign( d5 ); k * sd5 <= sd5 * k2; k += d5 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
            a = y(k)/c;
            x(k) += s*a;
            y(k) = a;
          }
        }
        if ( fold ) {
          c = -c;
          k0 = (n - i)*d2 + 1;
          fold = false;
          continue;
        }
        break;
      }
    }

    m = nover2*d2;
    k0 = m + 1;
    for ( int k1 = k0, sd3 = sign( d3 ); k1 * sd3 <= sd3 * d1; k1 += d3 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
      k2 = k1 + d4;
      for ( int k = k1, sd5 = sign( d5 ); k * sd5 <= sd5 * k2; k += d5 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
        j = k - m;
        l = k + m;
        a = y(l)*2.0f;
        x(k) += a;
        y(k) = a;
        x(l) = x(j);
        y(l) = -y(j);
      }
    }

  }

//---- Format statements
//
//FMT9000  FORMAT ("FFT error: SDIAD: N odd.  N =",i10)
}

void
inv21(
  FArray1A_float x,
  FArray1A_float y,
  int const n,
  int5a d
)
{
  x.dim( star );
  y.dim( star );
//      =========================
//
//---- Inverts fourier transform along a screw
//     diad. the result is scaled by n.

//     .. Local Scalars ..
  float a, b, c, c1, r, s, s1;
  int d1, d2, d3, d4, d5, j3, k, kk, l, ll, m, nover2;

  static float const pi = 3.141593f;

  d1 = d(1);
  d2 = d(2);
  d3 = d(3);
  d4 = d(4) - 1;
  d5 = d(5);

  nover2 = n/2;
  ll = n*d2;
  kk = nover2*d2;
  for ( int j1 = 1, sd3 = sign( d3 ); j1 * sd3 <= sd3 * d1; j1 += d3 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
    int const j2 = j1 + d4;
    for ( int j = j1, sd5 = sign( d5 ); j * sd5 <= sd5 * j2; j += d5 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
      l = ll + j;
      k = kk + j;
      x(l) = x(j) + x(k);
      x(k) += y(k);
      y(l) = 0.0f;
      y(k) = 0.0f;
    }
  }

  c1 = std::cos(pi/float(n));
  s1 = std::sin(pi/float(n));
  c = 1.0f;
  s = 0.0f;
  for ( int i = 2; i <= nover2; ++i ) {
    kk = (n + 2 - 2*i)*d2;
    ll = (n + 1 - i)*d2;
    r = c*c1 - s*s1;
    s = c*s1 + s*c1;
    c = r;
    int const j1 = (i - 1)*d2 + 1;
    for ( int j2 = j1, sd3 = sign( d3 ); j2 * sd3 <= sd3 * d1; j2 += d3 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
      j3 = j2 + d4;
      for ( int j = j2, sd5 = sign( d5 ); j * sd5 <= sd5 * j3; j += d5 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
        l = j + ll;
        k = j + kk;
        x(l) += x(j) + x(k);
        x(j) = y(j)*s + x(j);
        x(k) = y(k)*s + x(k);
        y(j) *= c;
        y(k) = -y(k)*c;
      }
    }
  }

  cmplft(x,y,n,d);

  for ( int i = 1; i <= nover2; ++i ) {
    kk = (n + 1 - 2*i)*d2;
    ll = i*d2 + kk;
    int const j1 = (i - 1)*d2 + 1;
    for ( int j2 = j1, sd3 = sign( d3 ); j2 * sd3 <= sd3 * d1; j2 += d3 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
      j3 = j2 + d4;
      for ( int j = j2, sd5 = sign( d5 ); j * sd5 <= sd5 * j3; j += d5 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
        k = j + kk;
        l = j + ll;
        a = x(j) - x(l);
        b = y(j) + y(l);
        x(j) = x(l);
        y(j) = -y(l);
#if (defined(__GNUC__) \
     && __GNUC__ == 4 && __GNUC_MINOR__ == 1)
        // Problem observed with g++ 4.1.1 (Fedora 6) and 4.1.2 (Fedora 8).
        // Trial to use raw_ref1 for x and y was not successful as a
        // work around.
        throw std::runtime_error(
          "g++ 4.1.x internal compiler error to be resolved.");
#else

        x(l) = x(k) + a;
        y(l) = y(k) - b;
        x(k) = a;
        y(k) = b;
#endif
      }
    }
  }

  m = n - 2;
  for ( int i = 1; i <= m; ++i ) {
    k = i;
    while ( true ) {
      int ko = k;
      k = ko/2;
      if ( 2*k != ko ) k = n - 1 - k;
      if ( k < i ) continue;
      if ( k == i ) {
      } else {
        kk = (k - i)*d2;
        int const j1 = i*d2 + 1;
        for ( int j2 = j1, sd3 = sign( d3 ); j2 * sd3 <= sd3 * d1; j2 += d3 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
          j3 = j2 + d4;
          for ( int j = j2, sd5 = sign( d5 ); j * sd5 <= sd5 * j3; j += d5 ) { //ObjexxF2Cxx Remove sign() and set <= or >=
            k = j + kk;
            a = x(k);
            b = y(k);
            x(k) = x(j);
            y(k) = y(j);
            x(j) = a;
            y(j) = b;
          }
        }
      }
      break;
    }
  }

}

} // namespace ccp4io_dev
