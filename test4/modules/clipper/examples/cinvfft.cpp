// Clipper app to perform inverse ffts
/* Copyright 2003-2004 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include "ccp4-extras.h"


int main( int argc, char** argv )
{
  CCP4program prog( "cinvfft", "0.1", "$Date: 2004/06/01" );

  // defaults
  clipper::String ipmap = "NONE";
  clipper::String ipfile = "NONE";
  clipper::String opfile = "NONE";
  clipper::String opcol = "invfft";

  // command input
  CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if ( args[arg] == "-mapin" ) {
      if ( ++arg < args.size() ) ipmap = args[arg];
    } else if ( args[arg] == "-mtzin" ) {
      if ( ++arg < args.size() ) ipfile = args[arg];
    } else if ( args[arg] == "-mtzout" ) {
      if ( ++arg < args.size() ) opfile = args[arg];
    } else if ( args[arg] == "-colout" ) {
      if ( ++arg < args.size() ) opcol = args[arg];
    } else {
      std::cout << "Unrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "Usage: cinvfft\n\t-mapin <filename>\n\t-mtzin <filename>\n\t-mtzout <filename>\n\t-colout <colpath>\nInverse FFT\n";
    exit(1);
  }

  // make data objects
  clipper::CCP4MTZfile mtzin, mtzout;
  clipper::MTZcrystal cxtl;
  clipper::Xmap<float> xmap;
  clipper::HKL_info hkls;

  // open file
  mtzin.open_read( ipfile );
  if ( opcol[0] != '/' ) opcol = mtzin.column_paths().back().notail()+"/"+opcol;
  mtzin.import_hkl_info( hkls );
  mtzin.import_crystal( cxtl, opcol );
  clipper::HKL_data<clipper::data32::F_phi> fphi( hkls, cxtl );
  mtzin.close_read();

  // import xmap
  clipper::CCP4MAPfile mapin;
  mapin.open_read( ipmap );
  mapin.import_xmap( xmap );
  mapin.close_read();

  // calc fft
  xmap.fft_to( fphi );

  // output data
  mtzout.open_append( ipfile, opfile );
  mtzout.export_hkl_data( fphi, opcol );
  mtzout.close_append();
}
