// Clipper app to convert HL coeffs to phi/fom or back
/* Copyright 2003-2004 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include "ccp4-extras.h"


int main( int argc, char** argv )
{
  CCP4program prog( "chltofom", "0.1", "$Date: 2004/05/01" );

  // defaults
  clipper::String ipfile = "NONE";
  clipper::String ipcolf = "NONE";
  clipper::String ipcolh = "NONE";
  clipper::String ipcolp = "NONE";
  clipper::String opfile = "hltofom.mtz";
  clipper::String opcol = "";

  // command input
  CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if ( args[arg] == "-mtzin" ) {
      if ( ++arg < args.size() ) ipfile = args[arg];
    } else if ( args[arg] == "-mtzout" ) {
      if ( ++arg < args.size() ) opfile = args[arg];
    } else if ( args[arg] == "-colin-hl" ) {
      if ( ++arg < args.size() ) ipcolh = args[arg];
    } else if ( args[arg] == "-colin-phifom" ) {
      if ( ++arg < args.size() ) ipcolp = args[arg];
    } else if ( args[arg] == "-colin-fo" ) {
      if ( ++arg < args.size() ) ipcolf = args[arg];
    } else if ( args[arg] == "-colout" ) {
      if ( ++arg < args.size() ) opcol = args[arg];
    } else {
      std::cout << "Unrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "Usage: chltofom\n\t-mtzin <filename>\n\t-mtzout <filename>\n\t-colin-fo <colpath>\n\t-colin-hl <colpath>\n\t-colin-phifom <colpath>\n\t-colout <colpath>\nIf -colin-hl is specified, conversion is ABCD->phi/fom.\nIf -colin-phifom is specified, conversion is phi/fom->ABCD.\nIf -colin-fo is specified, weighted map coefficients are calculated in addition.\nIf no FOM is present, it is set to 0.99.\n";
    exit(1);
  }

  // make data objects
  clipper::CCP4MTZfile mtzin, mtzout;
  clipper::HKL_info hkls;

  // actual work
  mtzin.open_read( ipfile );
  mtzin.import_hkl_info( hkls );
  clipper::HKL_data<clipper::data32::Phi_fom> phiw( hkls );
  clipper::HKL_data<clipper::data32::ABCD>    abcd( hkls );
  clipper::HKL_data<clipper::data32::F_sigF>  fsig( hkls );
  clipper::HKL_data<clipper::data32::F_phi>   fphi( hkls );
  if ( ipcolh != "NONE" ) mtzin.import_hkl_data( abcd, ipcolh );
  if ( ipcolp != "NONE" ) mtzin.import_hkl_data( phiw, ipcolp );
  if ( ipcolf != "NONE" ) mtzin.import_hkl_data( fsig, ipcolf );
  clipper::String opbase = mtzin.assigned_paths()[0];
  if      ( opcol == "" )     opcol = opbase.substr( 0, opbase.find(".") );
  else if ( opcol[0] != '/' ) opcol = opbase.notail()+"/"+opcol;
  mtzin.close_read();
  mtzout.open_append( ipfile, opfile );
  // compute phi+fom from abcd
  if ( ipcolh != "NONE" ) {
    phiw.compute( abcd, clipper::data32::Compute_phifom_from_abcd() );
    mtzout.export_hkl_data( phiw, opcol );
  }
  // compute abcd from phi+fom
  if ( ipcolp != "NONE" ) {
    // if weights absent, then set to 0.99
    clipper::HKL_data_base::HKL_reference_index ih;
    for ( ih = phiw.first(); !ih.last(); ih.next() )
      if ( !clipper::Util::is_nan( phiw[ih].fom() ) ) break;
    if ( ih.last() )
      for ( ih = phiw.first(); !ih.last(); ih.next() ) phiw[ih].fom() = 0.99;
    abcd.compute( phiw, clipper::data32::Compute_abcd_from_phifom() );
    mtzout.export_hkl_data( abcd, opcol );
  }
  // compute weighted F from F + phi
  if ( ipcolf != "NONE" ) {
    fphi.compute( fsig, phiw, clipper::data32::Compute_fphi_from_fsigf_phifom() );
    mtzout.export_hkl_data( fphi, opcol );
  }
  mtzout.close_append();
}
