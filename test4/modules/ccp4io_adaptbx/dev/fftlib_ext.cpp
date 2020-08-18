// NOT USED -- INCOMPLETE -- DOES NOT COMPILE
// here for future reference

#include <boost/python.hpp>

#include <ccp4io_adaptbx/dev/fftlib.hpp>

namespace ccp4io_dev { namespace ext {

  void
  fftlib_cmplft(
    af::ref<float> const & xy,
    int const n,
    af::ref<int> const & d)
  {
    ASSERTBX(n > 0);
    ASSERTBX(xy.size() >= 2 * n);
    ASSERTBX(d.size() == 5);
    cmplft(xy[0], xy[1], n, d[0]);
  }

  void
  fftlib_realft(
    af::ref<float> const & xy,
    int const n,
    af::ref<int> const & d)
  {
    ASSERTBX(n > 0);
    ASSERTBX(xy.size() >= 2 * n + 2);
    ASSERTBX(d.size() == 5);
    realft(xy[0], xy[1], n, d[0]);
  }

  void
  fftlib_hermft(
    af::ref<float> const & xy,
    int const n,
    af::ref<int> const & d)
  {
    ASSERTBX(n > 0);
    ASSERTBX(xy.size() >= 2 * n);
    ASSERTBX(d.size() == 5);
    hermft(xy[0], xy[1], n, d[0]);
  }

  void
  fftlib_real_to_complex_3d_real_imag_w(
    af::ref<float> const & x,
    int const nu,
    int const nv,
    int const nw,
    bool complete)
  {
    ASSERTBX(nu > 0);
    ASSERTBX(nv > 0);
    ASSERTBX(nw > 0);
    ASSERTBX(nw % 2 == 0);
    ASSERTBX(x.size() >= nu*nv*(nw+2)); // memory required also if not complete
    int nuv = nu*nv;
    int nuvw2 = (complete ? nuv*(nw+2) : nuv*nw);
    // transform along w
    int5d d;
    d(1)=nuvw2;
    d(2)=2*nuv;
    d(3)=nuvw2;
    d(4)=nuv;
    d(5)=1;
    realft(x[0],x[nuv],nw/2,d);
      // realft always does a complete transform, but
      // if !complete, the layer at nw/2 is ignored below
    // transform along v
    d(1)=nuvw2;
    d(2)=nu;
    d(3)=2*nuv;
    d(4)=nu;
    d(5)=1;
    cmplft(x[0],x[nuv],nv,d);
    // transform along u
    d(1)=nuvw2;
    d(2)=1;
    d(3)=2*nuv;
    d(4)=nuv;
    d(5)=nu;
    cmplft(x[0],x[nuv],nu,d);
  }

  void
  fftlib_complex_to_real_3d_real_imag_w(
    af::ref<float> const & x,
    int const nu,
    int const nv,
    int const nw,
    bool complete)
  {
    ASSERTBX(nu > 0);
    ASSERTBX(nv > 0);
    ASSERTBX(nw > 0);
    ASSERTBX(nw % 2 == 0);
    ASSERTBX(x.size() >= nu*nv*(nw+2)); // memory required also if not complete
    int nuv=nu*nv;
    int nuvw2 = (complete ? nuv*(nw+2) : nuv*nw);
    // transform along u
    int5d d;
    d(1)=nuvw2;
    d(2)=1;
    d(3)=2*nuv;
    d(4)=nuv;
    d(5)=nu;
    cmplft(x[0],x[nuv],nu,d);
    // transform along v
    d(1)=nuvw2;
    d(2)=nu;
    d(3)=2*nuv;
    d(4)=nu;
    d(5)=1;
    cmplft(x[0],x[nuv],nv,d);
    if (complete) {
      // copy real part of last value along w to imag part of first
      for(int v=0;v<nv;v++) {
        for(int u=0;u<nu;u++) {
          int uv1  = ((1 *nv)+v)*nu+u;
          int uvnw = ((nw*nv)+v)*nu+u;
          x[uv1] = x[uvnw];
        }
      }
    }
    // transform along w
    d(1)=nuvw2;
    d(2)=2*nuv;
    d(3)=nuvw2;
    d(4)=nuv;
    d(5)=1;
    hermft(x[0],x[nuv],nw/2,d);
  }

  void
  hermitian_conjugate_3d_real_imag_w(
    af::ref<float> const & x,
    int const nu,
    int const nv,
    int const nw)
  {
    ASSERTBX(nu > 0);
    ASSERTBX(nv > 0);
    ASSERTBX(nw > 0);
    ASSERTBX(nw % 2 == 0);
    ASSERTBX(x.size() >= nu*nv*(nw+2));
    for(int w=0;w<nw+2;w+=2) {
      for(int v=0;v<nv;v++) {
        for(int u=0;u<nu;u++) {
          int uvw_imag = (((w+1)*nv)+v)*nu+u;
          x[uvw_imag] *= -1.0f;
        }
      }
    }
  }

  void
  wrap()
  {
    using namespace boost::python;
    def("fftlib_cmplft", fftlib_cmplft, (arg("xy"), arg("n"), arg("d")));
    def("fftlib_realft", fftlib_realft, (arg("xy"), arg("n"), arg("d")));
    def("fftlib_hermft", fftlib_hermft, (arg("xy"), arg("n"), arg("d")));
    def("fftlib_real_to_complex_3d_real_imag_w",
         fftlib_real_to_complex_3d_real_imag_w, (
      arg("x"), arg("nu"), arg("nv"), arg("nw"), arg("complete")));
    def("fftlib_complex_to_real_3d_real_imag_w",
         fftlib_complex_to_real_3d_real_imag_w, (
      arg("x"), arg("nu"), arg("nv"), arg("nw"), arg("complete")));
    def("hermitian_conjugate_3d_real_imag_w",
         hermitian_conjugate_3d_real_imag_w, (
      arg("x"), arg("nu"), arg("nv"), arg("nw")));
  }

}} // namespace ccp4io_dev::ext

BOOST_PYTHON_MODULE(ccp4io_dev_ext)
{
  ccp4io_dev::ext::wrap();
}
