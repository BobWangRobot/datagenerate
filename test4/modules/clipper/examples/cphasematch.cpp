// Clipper app to match origins and calculate phase statistics
/* Copyright 2003-2004 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-ccp4.h>
#include "ccp4-extras.h"


int main( int argc, char** argv )
{
  CCP4program prog( "cphasematch", "0.1", "$Date: 2004/06/01" );

  // defaults
  clipper::String ipfile = "NONE";
  clipper::String ipcolfo = "NONE";
  clipper::String ipcolhl1 = "NONE";
  clipper::String ipcolpw1 = "NONE";
  clipper::String ipcolfc1 = "NONE";
  clipper::String ipcolhl2 = "NONE";
  clipper::String ipcolpw2 = "NONE";
  clipper::String ipcolfc2 = "NONE";
  clipper::String opfile = "NONE";
  clipper::String opcol = "phasematch";
  int res_bins = 12;
  int fom_bins = 20;
  bool omatch = true;

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
    } else if ( args[arg] == "-colin-hl-1" ) {
      if ( ++arg < args.size() ) ipcolhl1 = args[arg];
    } else if ( args[arg] == "-colin-phifom-1" ) {
      if ( ++arg < args.size() ) ipcolpw1 = args[arg];
    } else if ( args[arg] == "-colin-fc-1" ) {
      if ( ++arg < args.size() ) ipcolfc1 = args[arg];
    } else if ( args[arg] == "-colin-hl-2" ) {
      if ( ++arg < args.size() ) ipcolhl2 = args[arg];
    } else if ( args[arg] == "-colin-phifom-2" ) {
      if ( ++arg < args.size() ) ipcolpw2 = args[arg];
    } else if ( args[arg] == "-colin-fc-2" ) {
      if ( ++arg < args.size() ) ipcolfc2 = args[arg];
    } else if ( args[arg] == "-colout" ) {
      if ( ++arg < args.size() ) opcol = args[arg];
    } else if ( args[arg] == "-resolution-bins" ) {
      if ( ++arg < args.size() ) res_bins = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-no-origin-hand" ) {
      omatch = false;
    } else {
      std::cout << "Unrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "Usage: cphasematch\n\t-mtzin <filename>\n\t-mtzout <filename>\n\t-colin-fo <colpath>\n\t-colin-hl-1 <colpath>\n\t-colin-phifom-1 <colpath>\n\t-colin-fc-1 <colpath>\n\t-colin-hl-2 <colpath>\n\t-colin-phifom-2 <colpath>\n\t-colin-fc-2 <colpath>\n\t-colout <colpath>\n\t-resolution-bins <number-of-bins>\n\t-no-origin-match\nInput fo and any one of fc/phifom/abcd for each of dataset 1 and 2.\nIf mtzout is given, the second dataset is shifted to match the first.\n";
    exit(1);
  }

  // make data objects
  clipper::CCP4MTZfile mtzin, mtzout;
  clipper::HKL_info hkls;
  clipper::HKL_info::HKL_reference_index ih;

  // actual work
  mtzin.open_read( ipfile );
  mtzin.import_hkl_info( hkls );
  clipper::HKL_data<clipper::data32::F_sigF>  fsig( hkls );
  clipper::HKL_data<clipper::data32::Phi_fom> phiw1( hkls );
  clipper::HKL_data<clipper::data32::ABCD>    abcd1( hkls );
  clipper::HKL_data<clipper::data32::F_phi>   fphi1( hkls );
  clipper::HKL_data<clipper::data32::Phi_fom> phiw2( hkls );
  clipper::HKL_data<clipper::data32::ABCD>    abcd2( hkls );
  clipper::HKL_data<clipper::data32::F_phi>   fphi2( hkls );
  mtzin.import_hkl_data( fsig, ipcolfo );
  if ( ipcolhl1 != "NONE" ) mtzin.import_hkl_data( abcd1, ipcolhl1 );
  if ( ipcolpw1 != "NONE" ) mtzin.import_hkl_data( phiw1, ipcolpw1 );
  if ( ipcolfc1 != "NONE" ) mtzin.import_hkl_data( fphi1, ipcolfc1 );
  if ( ipcolhl2 != "NONE" ) mtzin.import_hkl_data( abcd2, ipcolhl2 );
  if ( ipcolpw2 != "NONE" ) mtzin.import_hkl_data( phiw2, ipcolpw2 );
  if ( ipcolfc2 != "NONE" ) mtzin.import_hkl_data( fphi2, ipcolfc2 );
  if ( opcol[0] != '/' ) opcol = mtzin.assigned_paths()[0].notail()+"/"+opcol;
  mtzin.close_read();

  // fill in data for input 1
  if ( ipcolhl1 != "NONE" ) {
    phiw1.compute( abcd1, clipper::data32::Compute_phifom_from_abcd() );
    fphi1.compute( fsig, phiw1, clipper::data32::Compute_fphi_from_fsigf_phifom() );
  } else if ( ipcolpw1 != "NONE" ) {
    abcd1.compute( phiw1, clipper::data32::Compute_abcd_from_phifom() );
    fphi1.compute( fsig, phiw1, clipper::data32::Compute_fphi_from_fsigf_phifom() );
  } else if ( ipcolfc1 != "NONE" ) {
    for ( ih = hkls.first(); !ih.last(); ih.next() ) {
      phiw1[ih].phi() = fphi1[ih].phi();
      phiw1[ih].fom() = 0.9998;
    }
    abcd1.compute( phiw1, clipper::data32::Compute_abcd_from_phifom() );
  } else {
    clipper::Message::message( clipper::Message_fatal( "Missing input 1" ) );
  }
  // fill in data for input 2
  if ( ipcolhl2 != "NONE" ) {
    phiw2.compute( abcd2, clipper::data32::Compute_phifom_from_abcd() );
    fphi2.compute( fsig, phiw2, clipper::data32::Compute_fphi_from_fsigf_phifom() );
  } else if ( ipcolpw2 != "NONE" ) {
    abcd2.compute( phiw2, clipper::data32::Compute_abcd_from_phifom() );
    fphi2.compute( fsig, phiw2, clipper::data32::Compute_fphi_from_fsigf_phifom() );
  } else if ( ipcolfc2 != "NONE" ) {
    for ( ih = hkls.first(); !ih.last(); ih.next() ) {
      phiw2[ih].phi() = fphi2[ih].phi();
      phiw2[ih].fom() = 0.9998;
    }
    abcd2.compute( phiw2, clipper::data32::Compute_abcd_from_phifom() );
  } else {
    clipper::Message::message( clipper::Message_fatal( "Missing input 2" ) );
  }

  // crude mean phase error
  double sw1 = 0.0, sd1 = 0.0;
  for ( ih = hkls.first(); !ih.last(); ih.next() )
    if ( !fphi1[ih].missing() && !fphi2[ih].missing() ) {
      double w = 1.0 / ih.hkl_class().epsilon();
      sw1 += w;
      sd1 += w * acos( cos( fphi1[ih].phi() - fphi2[ih].phi() ) );
    }

  // calculate phase shift
  bool invert( false );
  clipper::Coord_frac x( 0.0, 0.0, 0.0 );
  if ( omatch ) clipper::OriginMatch<float>( invert, x, fphi1, fphi2 );
  if ( invert ) std::cout << "\n Change of hand  : YES";
  else          std::cout << "\n Change of hand  : NO";
  std::cout << "\n Change or origin: " << x.format() << "\n";

  // shift phases
  for ( ih = hkls.first(); !ih.last(); ih.next() ) {
    clipper::Coord_reci_frac h( ih.hkl() );
    double dphi = clipper::Util::twopi() * ( h * x );
    if ( invert ) fphi2[ih].friedel();
    if ( invert ) phiw2[ih].friedel();
    if ( invert ) abcd2[ih].friedel();
    fphi2[ih].shift_phase( dphi );
    phiw2[ih].shift_phase( dphi );
    abcd2[ih].shift_phase( dphi );
  }

  // crude mean phase error
  double sw2 = 0.0, sd2 = 0.0;
  for ( ih = hkls.first(); !ih.last(); ih.next() )
    if ( !fphi1[ih].missing() && !fphi2[ih].missing() ) {
      double w = 1.0 / ih.hkl_class().epsilon();
      sw2 += w;
      sd2 += w * acos( cos( fphi1[ih].phi() - fphi2[ih].phi() ) );
    }

  // verify the origin shift with crude stats
  std::cout << "\n Mean phase error before origin fixing: " << clipper::Util::rad2d( sd1/sw1 );
  std::cout << "\n Mean phase error after  origin fixing: " << clipper::Util::rad2d( sd2/sw2 ) << "\n";

  // output shifted phases if required
  if ( opfile != "NONE" ) {
    mtzout.open_append( ipfile, opfile );
    if ( ipcolhl2 != "NONE" ) {
      mtzout.export_hkl_data( abcd2, opcol );
    } else if ( ipcolpw2 != "NONE" ) {
      mtzout.export_hkl_data( phiw2, opcol );
    } else if ( ipcolfc2 != "NONE" ) {
      mtzout.export_hkl_data( fphi2, opcol );
    }
    mtzout.close_append();
  }

  // NOW DO DETAILED REFLECTION STATISTICS

  // calculate E scale factor
  const int nprm = 10;
  clipper::HKL_data<clipper::data32::E_sigE> esig( hkls );
  esig.compute( fsig, clipper::data32::Compute_EsigE_from_FsigF() );
  std::vector<double> params_init( nprm, 1.0 );
  clipper::BasisFn_binner basis_fo( hkls, nprm, 2.0 );
  clipper::TargetFn_scaleEsq<clipper::data32::E_sigE> target_fo( esig );
  clipper::ResolutionFn escale( hkls, basis_fo, target_fo, params_init );

  // accumulate stats in bins and overall
  clipper::Resolution_ordinal resord;
  resord.init( hkls, 1.0 );
  int res_binx = res_bins + 1;
  double res_scl = 0.9999 * double( res_bins );
  std::vector<double> snr(res_binx,0.0),
    swr(res_binx,0.0),  sw1r(res_binx,0.0), sw2r(res_binx,0.0),
    sdr(res_binx,0.0), sdw1r(res_binx,0.0), sdw2r(res_binx,0.0),
    scr1(res_binx,0.0), scr2(res_binx,0.0), scr12(res_binx,0.0),
    scer1(res_binx,0.0), scer2(res_binx,0.0), scer12(res_binx,0.0),
    sww1(res_binx,0.0), sww2(res_binx,0.0),
    sdw1(res_binx,0.0), sdw2(res_binx,0.0);
  std::vector<double> nrfl1(fom_bins,0.0), nrfl2(fom_bins,0.0),
    mcosd1(fom_bins,0.0), mcosd2(fom_bins,0.0);
  for ( ih = hkls.first(); !ih.last(); ih.next() )
    if ( !fphi1[ih].missing() && !fphi2[ih].missing() ) {
      int bin = int( res_scl * resord.ordinal( ih.invresolsq() ) );
      double w = 1.0 / ih.hkl_class().epsilon();
      double cosd = cos( phiw1[ih].phi() - phiw2[ih].phi() );
      double d = acos(cosd);
      snr[bin]   += 1.0;
      swr[bin]   += w;
      sw1r[bin]  += w * phiw1[ih].fom();
      sw2r[bin]  += w * phiw2[ih].fom();
      sdr[bin]   += w * d;
      sdw1r[bin] += w * phiw1[ih].fom() * d;
      sdw2r[bin] += w * phiw2[ih].fom() * d;
      scr1[bin]  += w * fphi1[ih].f() * fphi1[ih].f();
      scr2[bin]  += w * fphi2[ih].f() * fphi2[ih].f();
      scr12[bin] += w * fphi1[ih].f() * fphi2[ih].f() * cosd;
      sww1[bin]  += w * phiw1[ih].fom() * phiw1[ih].fom();
      sww2[bin]  += w * phiw2[ih].fom() * phiw2[ih].fom();
      sdw1[bin]  += w * phiw1[ih].fom() * d;
      sdw2[bin]  += w * phiw2[ih].fom() * d;
      snr[res_bins]   += 1.0;
      swr[res_bins]   += w;
      sw1r[res_bins]  += w * phiw1[ih].fom();
      sw2r[res_bins]  += w * phiw2[ih].fom();
      sdr[res_bins]   += w * d;
      sdw1r[res_bins] += w * phiw1[ih].fom() * d;
      sdw2r[res_bins] += w * phiw2[ih].fom() * d;
      scr1[res_bins]  += w * fphi1[ih].f() * fphi1[ih].f();
      scr2[res_bins]  += w * fphi2[ih].f() * fphi2[ih].f();
      scr12[res_bins] += w * fphi1[ih].f() * fphi2[ih].f() * cosd;
      sww1[res_bins]  += w * phiw1[ih].fom() * phiw1[ih].fom();
      sww2[res_bins]  += w * phiw2[ih].fom() * phiw2[ih].fom();
      sdw1[res_bins]  += w * cosd * phiw1[ih].fom();
      sdw2[res_bins]  += w * cosd * phiw2[ih].fom();
      w *= escale.f(ih);
      scer1[bin]  += w * fphi1[ih].f() * fphi1[ih].f();
      scer2[bin]  += w * fphi2[ih].f() * fphi2[ih].f();
      scer12[bin] += w * fphi1[ih].f() * fphi2[ih].f() * cosd;
      scer1[res_bins]  += w * fphi1[ih].f() * fphi1[ih].f();
      scer2[res_bins]  += w * fphi2[ih].f() * fphi2[ih].f();
      scer12[res_bins] += w * fphi1[ih].f() * fphi2[ih].f() * cosd;
      bin = clipper::Util::bound(0, int(double(fom_bins)*phiw1[ih].fom()),
				 fom_bins-1);
      nrfl1[bin] += 1.0;
      mcosd1[bin] += cosd;
      bin = clipper::Util::bound(0, int(double(fom_bins)*phiw2[ih].fom()),
				 fom_bins-1);
      nrfl2[bin] += 1.0;
      mcosd2[bin] += cosd;
    }
  // finalise statistics
  for ( int bin = 0; bin < res_binx; bin++ ) {
    // phase differences
    sdr[bin]   /= clipper::Util::max( swr[bin], 1.0 );
    sdw1r[bin] /= clipper::Util::max( sw1r[bin], 1.0 );
    sdw2r[bin] /= clipper::Util::max( sw2r[bin], 1.0 );
    // mean foms
    sw1r[bin]  /= clipper::Util::max( swr[bin], 1.0 );
    sw2r[bin]  /= clipper::Util::max( swr[bin], 1.0 );
    // correlation
    scr12[bin]  /= sqrt( clipper::Util::max( scr1[bin]*scr2[bin], 1.0 ) );
    scer12[bin] /= sqrt( clipper::Util::max( scer1[bin]*scer2[bin], 1.0 ) );
    // Qfom
    sdw1[bin]   /= sww1[bin];
    sdw2[bin]   /= sww2[bin];
  }
  for ( int bin = 0; bin < fom_bins; bin++ ) {
    if ( nrfl1[bin] > 0.5 ) mcosd1[bin] /= nrfl1[bin];
    if ( nrfl2[bin] > 0.5 ) mcosd2[bin] /= nrfl2[bin];
  }
  // and display
  resord.invert();
  printf("\n $TABLE: Phase statistics with resolution:\n $GRAPHS:<fom> vs resolution:N:1,4,5:\n        :Unweighted and weighted <phase error> vs resolution:N:1,6,7,8:\n        :Reflection correlation vs resolution:N:1,9,10:\n $$\n 1/resol^2 hi   Nrefl <fom1> <fom2>  <dphi> w1<dphi> w2<dphi> wFcorr wEcorr $$\n $$\n");
  for ( int bin = 0; bin < res_bins; bin++ )
    printf( "%6.3f %6.3f %7i %6.3f %6.3f %7.2f  %7.2f  %7.2f %6.3f %6.3f\n",
	    resord.ordinal( double(bin)/res_scl ),
	    resord.ordinal( double(bin+0.999)/res_scl ),
	    int(snr[bin]), sw1r[bin], sw2r[bin],
	    clipper::Util::rad2d(sdr[bin]),
	    clipper::Util::rad2d(sdw1r[bin]),
	    clipper::Util::rad2d(sdw2r[bin]),
	    scr12[bin], scer12[bin] );
  printf(" $$\n");
  printf("\n $TABLE: Phase analysis with FOM:\n $GRAPHS:Histogram of FOM1:N:1,3:\n        :Phase difference vs FOM1:N:1,2,4:\n        :Histogram of FOM2:N:1,5:\n        :Phase difference vs FOM2:N:1,2,6:\n $$\n FOM       acos(FOM)  N(FOM1) acos<cos(dphi)> N(FOM2) acos<cos(dphi)>$$\n $$\n");
  for ( int bin = 0; bin < fom_bins; bin++ )
    printf( "%9.3f %9.3f %9.0f %9.3f %9.0f %9.3f\n",
	    (double(bin)+0.5)/double(fom_bins),
	    clipper::Util::rad2d(acos((double(bin)+0.5)/double(fom_bins))),
	    nrfl1[bin],clipper::Util::rad2d(acos(mcosd1[bin])),
	    nrfl2[bin],clipper::Util::rad2d(acos(mcosd2[bin])) );
  printf(" $$\n");
  printf("\nOverall statistics:\n  Nrefl <fom1> <fom2>  <dphi> w1<dphi> w2<dphi> wFcorr wEcorr    Qfom1  Qfom2\n");
  printf( "%7i %6.3f %6.3f %7.2f  %7.2f  %7.2f %6.3f %6.3f   %6.3f %6.3f\n",
	  int(snr[res_bins]), sw1r[res_bins], sw2r[res_bins],
	  clipper::Util::rad2d(sdr[res_bins]),
	  clipper::Util::rad2d(sdw1r[res_bins]),
	  clipper::Util::rad2d(sdw2r[res_bins]),
	  scr12[res_bins], scer12[res_bins],
	  sdw1[res_bins], sdw2[res_bins] );
}
