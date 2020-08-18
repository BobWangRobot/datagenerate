// Clipper app to do sigmaa calc
/* Copyright 2003-2004 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-ccp4.h>
#include "ccp4-extras.h"


int main( int argc, char** argv )
{
  CCP4program prog( "csigmaa", "0.1", "$Date: 2004/06/01" );

  // defaults
  clipper::String ipfile = "NONE";
  clipper::String ipcolfo = "NONE";
  clipper::String ipcolfc = "NONE";
  clipper::String ipcolfree = "NONE";
  clipper::String opfile = "sigmaa.mtz";
  clipper::String opcol = "sigmaa";
  int freeflag = 0;
  int n_refln = 1000;
  int n_param = 20;
  int verbose = 0;

  // command input
  CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if ( args[arg] == "-mtzin" ) {
      if ( ++arg < args.size() ) ipfile = args[arg];
    } else if ( args[arg] == "-mtzout" ) {
      if ( ++arg < args.size() ) opfile = args[arg];
    } else if ( args[arg] == "-colin-fo" ) {
      if ( ++arg < args.size() ) ipcolfo = args[arg];
    } else if ( args[arg] == "-colin-fc" ) {
      if ( ++arg < args.size() ) ipcolfc = args[arg];
    } else if ( args[arg] == "-colin-free" ) {
      if ( ++arg < args.size() ) ipcolfree = args[arg];
    } else if ( args[arg] == "-colout" ) {
      if ( ++arg < args.size() ) opcol = args[arg];
    } else if ( args[arg] == "-free-flag" ) {
      if ( ++arg < args.size() ) freeflag = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-num-reflns" ) {
      if ( ++arg < args.size() ) n_refln = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-num-params" ) {
      if ( ++arg < args.size() ) n_param = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-verbose" ) {
      if ( ++arg < args.size() ) verbose = clipper::String(args[arg]).i();
    } else {
      std::cout << "Unrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "Usage: csigmaa\n\t-mtzin <filename>\n\t-mtzout <filename>\n\t-colin-fo <colpath>\n\t-colin-fc <colpath>\n\t-colin-free <colpath>\n\t-colout <colpath>\n\t-free-flag <free set>\n\t-num-reflns <reflns per spline param>\n\t-num-params <spline params>\nCalculate HL coeffs from Fo and Fc.\n";
    exit(1);
  }

  // make data objects
  clipper::CCP4MTZfile mtzin, mtzout;
  clipper::HKL_info hkls;
  typedef clipper::HKL_data_base::HKL_reference_index HRI;

  mtzin.open_read( ipfile );
  mtzin.import_hkl_info( hkls );
  clipper::HKL_data<clipper::data32::F_sigF> fo( hkls );
  clipper::HKL_data<clipper::data32::F_phi>  fc( hkls );
  clipper::HKL_data<clipper::data32::Flag> free( hkls );
  mtzin.import_hkl_data( fo, ipcolfo );
  mtzin.import_hkl_data( fc, ipcolfc );
  if ( ipcolfree != "NONE" ) mtzin.import_hkl_data( free, ipcolfree );
  if ( opcol[0] != '/' ) opcol = mtzin.assigned_paths()[0].notail()+"/"+opcol;
  mtzin.close_read();

  // flag reflections
  clipper::HKL_data<clipper::data32::F_phi> fb( hkls ), fd( hkls );
  clipper::HKL_data<clipper::data32::Phi_fom> phiw( hkls );
  clipper::HKL_data<clipper::data32::Flag> flag( hkls );
  for ( HRI ih = flag.first(); !ih.last(); ih.next() )
    if ( !fo[ih].missing() && (free[ih].missing()||free[ih].flag()==freeflag) )
      flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
    else
      flag[ih].flag() = clipper::SFweight_spline<float>::NONE;

  // do sigmaa calc
  clipper::SFweight_spline<float> sfw( n_refln, n_param );
  bool fl = sfw( fb, fd, phiw, fo, fc, flag );

  // calc abcd
  clipper::HKL_data<clipper::data32::ABCD> abcd( hkls );
  abcd.compute( phiw, clipper::data32::Compute_abcd_from_phifom() );

  // output data
  mtzout.open_append( ipfile, opfile );
  mtzout.export_hkl_data( abcd, opcol );
  mtzout.export_hkl_data( fb, opcol+"_BEST" );
  mtzout.export_hkl_data( fd, opcol+"_DIFF" );
  mtzout.close_append();

  // DIAGNOSTIC OUTPUT
  if ( verbose > 1 ) {
    std::cout << "\nNumber of spline params: " << sfw.params_scale().size() << "\n";
    clipper::BasisFn_spline basisfn( hkls, sfw.params_scale().size(), 1.0 );
    printf("\n $TABLE: Sigmaa statistics :\n $GRAPHS:scale vs resolution:N:1,2:\n        :lack of closure vs resolution:N:1,3:\n $$\n 1/resol^2   scale   lack_of_closure $$\n $$\n");
    for ( int i = 0; i <= 20.0; i++ ) {
      double s = hkls.resolution().invresolsq_limit()*double(i)/20.0;
      printf( "%6.3f %12.3f %12.3f\n",
	      s, basisfn.f_s(s,sfw.params_scale()), basisfn.f_s(s,sfw.params_error()) );
    }
    printf(" $$\n");
  }
}
