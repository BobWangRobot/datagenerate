# NOT USED -- INCOMPLETE -- DOES NOT RUN
# here for future reference

from __future__ import division
from __future__ import print_function
import ccp4io_dev_ext
import scitbx.fftpack
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.math_utils import prime_factors_of
from libtbx.utils import format_cpu_times
import sys

def compare_fftpack_with_cmplft_1d():
  mt = flex.mersenne_twister(seed=0)
  for n in xrange(1, 101):
    primes = prime_factors_of(n)
    if (n != 1 and max(primes) > 19): continue
    z = (mt.random_double(size=n*2)*2-1).as_float()
    cmplft_xy = z.deep_copy()
    d = flex.int((2*n,2,2*n,2*n,2*n))
    ccp4io_dev_ext.fftlib_cmplft(xy=cmplft_xy, n=n, d=d)
    fft = scitbx.fftpack.complex_to_complex(n)
    fftpack_xy = z.as_double()
    fft.forward(fftpack_xy)
    fftpack_xy = fftpack_xy.as_float()
    if (flex.max_absolute(cmplft_xy-fftpack_xy) > 1e-5):
      assert approx_equal(cmplft_xy, fftpack_xy, eps=1e-5)

def compare_fftpack_with_cmplft_3d():
  mt = flex.mersenne_twister(seed=0)
  for nx,ny,nz in [(30,20,40), (7,19,13), (5,11,4)]:
    z = (mt.random_double(size=2*nx*ny*nz)*2-1).as_float()
    cmplft_xy = z.deep_copy()
    d = flex.int([2*nx*ny*nz, 2*nx*ny, 2*nx*ny*nz, 2*nx*ny, 2])
    ccp4io_dev_ext.fftlib_cmplft(xy=cmplft_xy, n=nz, d=d)
    d = flex.int([2*nx*ny*nz, 2*nx, 2*nx*ny, 2*nx, 2])
    ccp4io_dev_ext.fftlib_cmplft(xy=cmplft_xy, n=ny, d=d)
    d = flex.int([2*nx*ny*nz, 2, 2*nx*ny*nz, 2*nx*ny*nz, 2*nx])
    ccp4io_dev_ext.fftlib_cmplft(xy=cmplft_xy, n=nx, d=d)
    fft = scitbx.fftpack.complex_to_complex_3d((nz,ny,nx))
    fftpack_xy = z.as_double()
    fftpack_xy.reshape(flex.grid(nz,ny,2*nx))
    fft.forward(fftpack_xy)
    fftpack_xy = fftpack_xy.as_float().as_1d()
    if (flex.max_absolute(cmplft_xy-fftpack_xy) > 1e-4):
      assert approx_equal(cmplft_xy, fftpack_xy, eps=1e-4)

def compare_fftpack_with_realft_1d():
  mt = flex.mersenne_twister(seed=0)
  for n_cmpl in xrange(1, 101):
    primes = prime_factors_of(n_cmpl)
    if (n_cmpl != 1 and max(primes) > 19): continue
    n_real = n_cmpl * 2
    m_real = n_real + 2
    z = (mt.random_double(size=m_real)*2-1).as_float()
    realft_xy = z.deep_copy()
    d = flex.int((m_real,2,m_real,m_real,m_real))
    ccp4io_dev_ext.fftlib_realft(xy=realft_xy, n=n_cmpl, d=d)
    fft = scitbx.fftpack.real_to_complex(n_real)
    assert fft.m_real() == m_real
    fftpack_xy = z.as_double()
    fft.forward(fftpack_xy)
    fftpack_xy = fftpack_xy.as_float()
    if (flex.max_absolute(realft_xy-fftpack_xy) > 1e-5):
      assert approx_equal(realft_xy, fftpack_xy, eps=1e-5)

def compare_fftpack_with_hermft_1d():
  mt = flex.mersenne_twister(seed=0)
  for n_cmpl in xrange(1, 101):
    primes = prime_factors_of(n_cmpl)
    if (n_cmpl != 1 and max(primes) > 19): continue
    n_real = n_cmpl * 2
    m_real = n_real + 2
    z = (mt.random_double(size=m_real)*2-1).as_float()
      # The imaginary parts of the first and last complex values should
      # be zero, but z has random values in these places, to prove that
      # they don't change the result.
    hermft_xy = z.deep_copy()
    hermft_xy[1] = hermft_xy[-2]
      # real part of last complex value stored in imaginary part of first
      # (see ccp4/doc/libfft.doc)
    hermft_xy[-2] = 253
      # random value, for consistency check below
    d = flex.int((m_real,2,m_real,m_real,m_real))
    ccp4io_dev_ext.fftlib_hermft(xy=hermft_xy, n=n_cmpl, d=d)
    # consistency check first
    assert hermft_xy[-2] == 253
    assert hermft_xy[-1] == z[-1]
    # reset to zero for comparision with fftpack result further down
    hermft_xy[-2] = 0
    hermft_xy[-1] = 0
    fft = scitbx.fftpack.real_to_complex(n_real)
    assert fft.m_real() == m_real
    fftpack_xy = z.as_double()
    for i in xrange(n_cmpl): fftpack_xy[i*2+1] *= -1 # conjugate
    fft.backward(fftpack_xy)
    fftpack_xy = fftpack_xy.as_float()
    if (flex.max_absolute(hermft_xy-fftpack_xy) > 1e-5):
      assert approx_equal(hermft_xy, fftpack_xy, eps=1e-5)

def show_complete_true_false_cc(nu, nv, nw, recycled, verbose):
  nuvw = nu*nv*nw
  corr = flex.linear_correlation(
    x=recycled[0][:nuvw].as_double(),
    y=recycled[1][:nuvw].as_double())
  if (verbose):
    print("dims:", (nu,nv,nw), "complete=true,false cc:", corr.coefficient())

def exercise_fftlib_real_complex_3d_real_imag_w_given_dims(
      mt, nu, nv, nw, verbose):
  map_size = nu * nv * (nw+2)
  map = (mt.random_double(size=map_size)*2-1).as_float()
  x = map.deep_copy()
  ccp4io_dev_ext.fftlib_real_to_complex_3d_real_imag_w(
    x=x, nu=nu, nv=nv, nw=nw, complete=True)
  f = x.deep_copy()
  ccp4io_dev_ext.hermitian_conjugate_3d_real_imag_w(x=x, nu=nu, nv=nv, nw=nw)
  ccp4io_dev_ext.fftlib_complex_to_real_3d_real_imag_w(
    x=x, nu=nu, nv=nv, nw=nw, complete=True)
  x *= 1/(nu*nv*nw)
  xfocus = x[:nu*nv*nw]
  mfocus = map[:nu*nv*nw]
  if (flex.max_absolute(xfocus-mfocus) > 1e-5):
    assert approx_equal(xfocus, mfocus, eps=1e-5)
  #
  fft = scitbx.fftpack.real_to_complex_3d((nu, nv, nw))
  g = flex.double(flex.grid(fft.m_real()).set_focus(fft.n_real()))
  assert g.size() == f.size()
  for u in xrange(nu):
    for v in xrange(nv):
      for w in xrange(nw+2):
        g[(u*nv+v)*(nw+2)+w] = map[(w*nv+v)*nu+u]
  fft.forward(g)
  for u in xrange(nu):
    for v in xrange(nv):
      for w in xrange(nw//2+1):
        wr = w*2
        wi = wr+1
        for wj in [wr, wi]:
          ge = g[(u*nv+v)*(nw+2)+wj]
          fe = f[(wj*nv+v)*nu+u]
          assert approx_equal(fe, ge, eps=1e-4)
  #
  y = map.deep_copy()
  ccp4io_dev_ext.fftlib_real_to_complex_3d_real_imag_w(
    x=y, nu=nu, nv=nv, nw=nw, complete=False)
  ccp4io_dev_ext.hermitian_conjugate_3d_real_imag_w(x=y, nu=nu, nv=nv, nw=nw)
  ccp4io_dev_ext.fftlib_complex_to_real_3d_real_imag_w(
    x=y, nu=nu, nv=nv, nw=nw, complete=False)
  y *= 1/(nu*nv*nw)
  show_complete_true_false_cc(nu, nv, nw, (x, y), verbose)

def exercise_fftlib_real_complex_3d_real_imag_w(verbose):
  mt = flex.mersenne_twister(seed=0)
  for nu in xrange(1,9):
    for nv in xrange(1,9):
      for nw in xrange(2,9,2):
        exercise_fftlib_real_complex_3d_real_imag_w_given_dims(
          mt, nu, nv, nw, verbose)

def compare_large_3d_real_imag_w_complete_true_false(verbose):
  mt = flex.mersenne_twister(seed=0)
  nu, nv, nw = 50, 52, 54
  map_size = nu * nv * (nw+2)
  map = (mt.random_double(size=map_size)*2-1).as_float()
  recycled = []
  for complete in [True, False]:
    x = map.deep_copy()
    ccp4io_dev_ext.fftlib_real_to_complex_3d_real_imag_w(
      x=x, nu=nu, nv=nv, nw=nw, complete=complete)
    ccp4io_dev_ext.hermitian_conjugate_3d_real_imag_w(x=x, nu=nu, nv=nv, nw=nw)
    ccp4io_dev_ext.fftlib_complex_to_real_3d_real_imag_w(
      x=x, nu=nu, nv=nv, nw=nw, complete=complete)
    recycled.append(x[:nu*nv*nw])
  show_complete_true_false_cc(nu, nv, nw, recycled, verbose)

def run(args):
  assert args in [[], ["--verbose"]]
  verbose = (len(args) != 0)
  compare_fftpack_with_cmplft_1d()
  compare_fftpack_with_cmplft_3d()
  compare_fftpack_with_realft_1d()
  compare_fftpack_with_hermft_1d()
  exercise_fftlib_real_complex_3d_real_imag_w(verbose=verbose)
  compare_large_3d_real_imag_w_complete_true_false(verbose=verbose)
  print(format_cpu_times())

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
