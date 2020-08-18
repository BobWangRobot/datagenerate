// Clipper app to perform ffts and map stats
/* Copyright 2003-2004 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-ccp4.h>
#include "ccp4-extras.h"


int main( int argc, char** argv )
{
  CCP4program prog( "cpatterson", "0.1", "$Date: 2004/07/01" );

  // defaults
  clipper::String ipfile = "NONE";
  clipper::String ipcolf = "NONE";
  clipper::String ipcold = "NONE";
  clipper::String ipcola = "NONE";
  clipper::String ipcolh = "NONE";
  clipper::String opfile = "patterson.map";
  bool adiff = false;
  bool oremv = false;
  double limit_e = 1.0e20;
  double weight_e = 0.0;
  clipper::Resolution reso;
  clipper::Grid_sampling grid;
  const int nprm = 12;

  // command input
  CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if ( args[arg] == "-mtzin" ) {
      if ( ++arg < args.size() ) ipfile = args[arg];
    } else if ( args[arg] == "-mapout" ) {
      if ( ++arg < args.size() ) opfile = args[arg];
    } else if ( args[arg] == "-colin-fo" ) {
      if ( ++arg < args.size() ) ipcolf = args[arg];
    } else if ( args[arg] == "-colin-fano" ) {
      if ( ++arg < args.size() ) ipcola = args[arg];
    } else if ( args[arg] == "-colin-fdiff" ) {
      if ( ++arg < args.size() ) ipcold = args[arg];
    } else if ( args[arg] == "-resolution" ) {
      if ( ++arg < args.size() ) {
	reso = clipper::Resolution( clipper::String(args[arg]).f() );
      }
    } else if ( args[arg] == "-grid" ) {
      if ( ++arg < args.size() ) {
	std::vector<clipper::String> g = clipper::String(args[arg]).split(", ");
	grid = clipper::Grid_sampling( g[0].i(), g[1].i(), g[2].i() );
      }
    } else if ( args[arg] == "-e-limit" ) {
      if ( ++arg < args.size() ) limit_e = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-e-weight" ) {
      if ( ++arg < args.size() ) weight_e = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-anomalous" ) {
      adiff = true;
    } else if ( args[arg] == "-origin-removal" ) {
      oremv = true;
    } else {
      std::cout << "Unrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "Usage: cpatterson\n\t-mtzin <filename>\n\t-mapout <filename>\n\t-colin-fo <colpath>\n\t-colin-fano <colpath>\n\t-colin-fdiff <colpath>\n\t-resolution <reso>\n\t-grid <nu>,<nv>,<nw>\n\t-e-limit <limit>\n\t-e-weight <weight>\n\t-anomalous\n\t-origin-removal\nCalculate Patterson from Fobs or Fano.\nFor E-Patterson <weight> = 1, for F-Patterson <weight> = 0, and other values.\n<limit> can reject reflections or differences by E-value.\n";
    exit(1);
  }

  // make data objects
  clipper::CCP4MTZfile mtzin;
  clipper::MTZcrystal cxtl;
  clipper::HKL_info hkls, hklp;
  typedef clipper::HKL_data_base::HKL_reference_index HRI;

  // open file
  mtzin.open_read( ipfile );
  clipper::Spacegroup spgr = mtzin.spacegroup();
  clipper::Cell       cell = mtzin.cell();
  if ( ipcolf != "NONE" ) mtzin.import_crystal( cxtl, ipcolf );
  if ( ipcola != "NONE" ) mtzin.import_crystal( cxtl, ipcola );
  if ( !cxtl.is_null() ) cell = cxtl;
  if ( reso.is_null() ) reso = mtzin.resolution();
  hkls.init( spgr, cell, reso, true );
  clipper::HKL_data<clipper::data32::F_sigF>     fsig( hkls );
  clipper::HKL_data<clipper::data32::F_sigF>     dsig( hkls );
  clipper::HKL_data<clipper::data32::F_sigF_ano> fano( hkls );
  if ( ipcolf != "NONE" ) mtzin.import_hkl_data( fsig, ipcolf );
  if ( ipcola != "NONE" ) mtzin.import_hkl_data( fano, ipcola );
  if ( ipcold != "NONE" ) mtzin.import_hkl_data( dsig, ipcolf );
  mtzin.close_read();

  // make mean/difference F if necessary
  if ( ipcola != "NONE" )
    if ( adiff )
      fsig.compute( fano, clipper::data32::Compute_diff_fsigf_from_fsigfano() );
    else
      fsig.compute( fano, clipper::data32::Compute_mean_fsigf_from_fsigfano() );

  // subtract difference F if necessary
  if ( ipcold != "NONE" ) 
    for ( HRI ih = fsig.first(); !ih.last(); ih.next() )
      if ( !fsig[ih].missing() && !dsig[ih].missing() )
	fsig[ih].f() -= dsig[ih].f();
      else
	fsig[ih].f() = 0.0;

  // calculate E scale factor
  clipper::HKL_data<clipper::data32::E_sigE> esig( hkls );
  esig.compute( fsig, clipper::data32::Compute_EsigE_from_FsigF() );
  std::vector<double> params_init( nprm, 1.0 );
  clipper::BasisFn_spline basis_fo( esig, nprm, 2.0 );
  clipper::TargetFn_scaleEsq<clipper::data32::E_sigE> target_fo( esig );
  clipper::ResolutionFn escale( hkls, basis_fo, target_fo, params_init );

  // reject reflections by E-value
  for ( HRI ih = esig.first(); !ih.last(); ih.next() )
    if ( !esig[ih].missing() )
      if ( esig[ih].E() / sqrt( escale.f(ih) ) > limit_e )
	fsig[ih] = clipper::data32::F_sigF();

  // calculate E-F combination
  for ( HRI ih = fsig.first(); !ih.last(); ih.next() )
    if ( !fsig[ih].missing() )
      fsig[ih].scale( pow( escale.f(ih), 0.5*weight_e ) );

  // get Patterson spacegroup
  clipper::Spacegroup
    pspgr( clipper::Spgr_descr( spgr.generator_ops().patterson_ops() ) );
  hklp.init( pspgr, cell, reso, true );

  // make patterson coeffs
  clipper::HKL_data<clipper::data32::F_phi> fphi( hklp );
  for ( HRI ih = fphi.first(); !ih.last(); ih.next() ) {
    clipper::data32::F_sigF f = fsig[ih.hkl()];
    if ( !f.missing() ) {
      fphi[ih].f() = f.f()*f.f();
      fphi[ih].phi() = 0.0 ;
    }
  }

  // origin removal
  if ( oremv ) {
    clipper::BasisFn_spline basis_fp( fphi, nprm, 2.0 );
    clipper::TargetFn_meanFnth<clipper::data32::F_phi> target_fp( fphi, 1.0 );
    clipper::ResolutionFn oscale( hklp, basis_fp, target_fp, params_init );
    for ( HRI ih = fphi.first(); !ih.last(); ih.next() )
      if ( !fphi[ih].missing() )
	fphi[ih].f() -= oscale.f(ih);
  }

  // make grid if necessary
  if ( grid.is_null() ) grid.init( pspgr, cell, reso );

  // make xmap
  clipper::Xmap<float> xmap( pspgr, cell, grid );
  xmap.fft_from( fphi );

  // write map
  clipper::CCP4MAPfile mapout;
  mapout.open_write( opfile );
  mapout.export_xmap( xmap );
  mapout.close_write();
}
