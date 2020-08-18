// Clipper app to combine HL coeffs
/* Copyright 2003-2004 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include "ccp4-extras.h"


int main( int argc, char** argv )
{
  CCP4program prog( "cphasecombine", "0.1", "$Date: 2004/06/01" );

  // defaults
  clipper::String ipfile = "NONE";
  clipper::String ipcolh1 = "NONE";
  clipper::String ipcolh2 = "NONE";
  clipper::String opfile = "phasecombine.mtz";
  clipper::String opcol = "hlcomb";
  float hlwt1 = 1.0;
  float hlwt2 = 1.0;

  // command input
  CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if ( args[arg] == "-mtzin" ) {
      if ( ++arg < args.size() ) ipfile = args[arg];
    } else if ( args[arg] == "-mtzout" ) {
      if ( ++arg < args.size() ) opfile = args[arg];
    } else if ( args[arg] == "-colin-hl-1" ) {
      if ( ++arg < args.size() ) ipcolh1 = args[arg];
    } else if ( args[arg] == "-colin-hl-2" ) {
      if ( ++arg < args.size() ) ipcolh2 = args[arg];
    } else if ( args[arg] == "-weight-hl-1" ) {
      if ( ++arg < args.size() ) hlwt1 = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-weight-hl-2" ) {
      if ( ++arg < args.size() ) hlwt2 = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-colout" ) {
      if ( ++arg < args.size() ) opcol = args[arg];
    } else {
      std::cout << "Unrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "Usage: cphasecombine\n\t-mtzin <filename>\n\t-mtzout <filename>\n\t-colin-hl-1 <colpath>\n\t-colin-hl-2 <colpath>\n\t-colout <colpath>\n\t-weight-hl-1 <weight>\n\t-weight-hl-2 <weight>\nThe specified phase probabilities, given by HL coefficients, are combined.\nIf weights are supplied, these are applied.\nIf the second set of coefficients are omitted, scaling alone occurs.\n";
    exit(1);
  }

  // make data objects
  clipper::CCP4MTZfile mtzin, mtzout;
  clipper::HKL_info hkls;

  // input
  mtzin.open_read( ipfile );
  mtzin.import_hkl_info( hkls );
  clipper::HKL_data<clipper::data32::ABCD>    abcd1( hkls );
  clipper::HKL_data<clipper::data32::ABCD>    abcd2( hkls );
  if ( ipcolh1 != "NONE" ) mtzin.import_hkl_data( abcd1, ipcolh1 );
  if ( ipcolh2 != "NONE" ) mtzin.import_hkl_data( abcd2, ipcolh2 );
  if ( opcol[0] != '/' ) opcol = mtzin.assigned_paths()[0].notail()+"/"+opcol;
  mtzin.close_read();

  // actual work
  clipper::HKL_data<clipper::data32::ABCD> abcd( hkls );
  for ( clipper::HKL_data<clipper::data32::ABCD>::HKL_reference_index
	  ih = abcd.first(); !ih.last(); ih.next() ) {
    abcd[ih].a() = abcd[ih].b() = abcd[ih].c() = abcd[ih].d() = 0.0;
    if ( !abcd1[ih].missing() ) {
      abcd[ih].a() += hlwt1 * abcd1[ih].a();
      abcd[ih].b() += hlwt1 * abcd1[ih].b();
      abcd[ih].c() += hlwt1 * abcd1[ih].c();
      abcd[ih].d() += hlwt1 * abcd1[ih].d();
    }
    if ( !abcd2[ih].missing() ) {
      abcd[ih].a() += hlwt2 * abcd2[ih].a();
      abcd[ih].b() += hlwt2 * abcd2[ih].b();
      abcd[ih].c() += hlwt2 * abcd2[ih].c();
      abcd[ih].d() += hlwt2 * abcd2[ih].d();
    }
  }

  // output
  mtzout.open_append( ipfile, opfile );
  mtzout.export_hkl_data( abcd, opcol );
  mtzout.close_append();
}
