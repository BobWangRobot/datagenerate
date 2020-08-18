#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/cns/cns_hkl_io.h>

#include <iostream>


using namespace clipper;
using namespace clipper::data32;


int main(int argc, char** argv)
{
  CCP4MTZfile file;

  // import an mtz
  HKL_info mydata;
  HKL_data<F_sigF> myfsig( mydata );
  HKL_data<Phi_fom> myphwt( mydata );
  MTZcrystal xtl;
  MTZdataset set;

  file.open_read( argv[1] );
  file.import_hkl_info( mydata, false );
  file.import_hkl_data( myfsig, set, xtl, "*/*/[FP SIGFP]" );
  file.import_hkl_data( myphwt, set, xtl, "*/*/[PHIB FOM]" );
  file.close_read();

  CNS_HKLfile cns;
  cns.open_write( "1.hkl" );
  cns.export_hkl_info( mydata );
  cns.export_hkl_data( myfsig );
  cns.export_hkl_data( myphwt );
  cns.close_write();

  HKL_info mydata2( mydata.spacegroup(), mydata.cell(), mydata.resolution() );
  HKL_data<F_sigF> myfsig2( mydata2 );
  HKL_data<Phi_fom> myphwt2( mydata2 );
  cns.open_read( "1.hkl" );
  cns.import_hkl_info( mydata2 );
  cns.import_hkl_data( myfsig2 );
  cns.import_hkl_data( myphwt2 );
  cns.close_read();

  cns.open_write( "2.hkl" );
  cns.export_hkl_info( mydata2 );
  cns.export_hkl_data( myfsig2 );
  cns.export_hkl_data( myphwt2 );
  cns.close_write();
}
