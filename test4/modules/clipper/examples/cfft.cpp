// Clipper app to perform ffts and map stats
/* Copyright 2003-2004 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-ccp4.h>
#include "ccp4-extras.h"


double sfunc( double q, int i, int j );
double pfunc( double t, int nt );


int main( int argc, char** argv )
{
  CCP4program prog( "cfft", "0.1", "$Date: 2004/05/01" );

  // defaults
  clipper::String ipfile = "NONE";
  clipper::String ipcolf = "NONE";
  clipper::String ipcold = "NONE";
  clipper::String ipcola = "NONE";
  clipper::String ipcolh = "NONE";
  clipper::String ipcolm = "NONE";
  clipper::String opfile = "NONE";
  bool stats = false;
  bool adiff = false;
  double statsrad = -1.0;
  double uvalue = 0.0;
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
    } else if ( args[arg] == "-colin-hl" ) {
      if ( ++arg < args.size() ) ipcolh = args[arg];
    } else if ( args[arg] == "-colin-fc" ) {
      if ( ++arg < args.size() ) ipcolm = args[arg];
    } else if ( args[arg] == "-u-value" ) {
      if ( ++arg < args.size() ) uvalue = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-b-value" ) {
      if ( ++arg < args.size() ) uvalue = clipper::Util::b2u(clipper::String(args[arg]).f());
    } else if ( args[arg] == "-resolution" ) {
      if ( ++arg < args.size() ) {
	reso = clipper::Resolution( clipper::String(args[arg]).f() );
      }
    } else if ( args[arg] == "-grid" ) {
      if ( ++arg < args.size() ) {
	std::vector<clipper::String> g = clipper::String(args[arg]).split(", ");
	grid = clipper::Grid_sampling( g[0].i(), g[1].i(), g[2].i() );
      }
    } else if ( args[arg] == "-anomalous" ) {
      adiff = true;
    } else if ( args[arg] == "-stats" ) {
      stats = true;
    } else if ( args[arg] == "-stats-radius" ) {
      if ( ++arg < args.size() ) statsrad = clipper::String(args[arg]).f();
    } else {
      std::cout << "Unrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "Usage: cfft\n\t-mtzin <filename>\n\t-mapout <filename>\n\t-colin-fo <colpath>\n\t-colin-fano <colpath>\n\t-colin-fdiff <colpath>\n\t-colin-hl <colpath>\n\t-colin-fc <colpath>\n\t-u-value <U>\n\t-b-value <B>\n\t-resolution <reso>\n\t-grid <nu>,<nv>,<nw>\n\t-anomalous\n\t-stats\n\t-stats-radius <radius>\nIf -colin-hl is specified, conversion is ABCD->phi/fom.\nIf -colin-fo and -colin-hl are specified, they are used to calculate\nmap coefficient. Otherwise the map coefficient may be provided using\n-colin-fc.\n";
    exit(1);
  }

  // make data objects
  clipper::CCP4MTZfile mtzin;
  clipper::MTZcrystal cxtl;
  clipper::HKL_info hkls;
  typedef clipper::HKL_data_base::HKL_reference_index HRI;

  // open file
  mtzin.open_read( ipfile );
  clipper::Spacegroup spgr = mtzin.spacegroup();
  clipper::Cell       cell = mtzin.cell();
  if ( ipcolf != "NONE" ) mtzin.import_crystal( cxtl, ipcolf );
  if ( ipcola != "NONE" ) mtzin.import_crystal( cxtl, ipcola );
  if ( ipcolm != "NONE" ) mtzin.import_crystal( cxtl, ipcolm );
  if ( !cxtl.is_null() ) cell = cxtl;
  if ( reso.is_null() ) reso = mtzin.resolution();
  hkls.init( spgr, cell, reso, true );
  clipper::HKL_data<clipper::data32::ABCD>       abcd( hkls );
  clipper::HKL_data<clipper::data32::F_sigF>     fsig( hkls );
  clipper::HKL_data<clipper::data32::F_sigF>     dsig( hkls );
  clipper::HKL_data<clipper::data32::F_sigF_ano> fano( hkls );
  clipper::HKL_data<clipper::data32::F_phi>      fphi( hkls );
  if ( ipcolf != "NONE" ) mtzin.import_hkl_data( fsig, ipcolf );
  if ( ipcola != "NONE" ) mtzin.import_hkl_data( fano, ipcola );
  if ( ipcold != "NONE" ) mtzin.import_hkl_data( dsig, ipcolf );
  if ( ipcolh != "NONE" ) mtzin.import_hkl_data( abcd, ipcolh );
  if ( ipcolm != "NONE" ) mtzin.import_hkl_data( fphi, ipcolm );
  mtzin.close_read();

  // make anomalous F if necessary
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

  // make map coeffs if necessary
  if ( ipcolm == "NONE" ) {
    clipper::HKL_data<clipper::data32::Phi_fom> phiw( hkls );
    phiw.compute( abcd, clipper::data32::Compute_phifom_from_abcd() );
    fphi.compute( fsig, phiw,clipper::data32::Compute_fphi_from_fsigf_phifom());
  }

  // shift phase for anomalous differences is necessary
  if ( adiff )
    for ( HRI ih = fphi.first(); !ih.last(); ih.next() )
      fphi[ih].shift_phase( 0.5*clipper::Util::pi() );

  // apply U value
  fphi.compute( fphi, clipper::data32::Compute_scale_u_fphi( 1.0, -uvalue ) );

  // make grid if necessary
  if ( grid.is_null() ) grid.init( spgr, cell, reso );

  // make xmap
  clipper::Xmap<float> xmap( spgr, cell, grid );
  xmap.fft_from( fphi );

  // write map
  if ( opfile != "NONE" ) {
    clipper::CCP4MAPfile mapout;
    mapout.open_write( opfile );
    mapout.export_xmap( xmap );
    mapout.close_write();
  }

  /*
    ------------------------------------------------------------------------
    Calculation is done. Now calculate some stats.
    ------------------------------------------------------------------------
  */

  // calculate stats
  if ( stats ) {
    double m1z, m1m, m2z, m2m, m3z, m3m, mn, mi, ma, w, r;
    m1z = m1m = m2z = m2m = m3z = m3m = mn = 0.0;
    mi = 1.0e20; ma = -1.0e20;

    clipper::data32::F_phi f000 = fphi[clipper::HKL(0,0,0)];
    if ( !f000.missing() ) m1z = f000.f() / cell.volume();
    clipper::Xmap<float>::Map_reference_index ix( xmap );
    for ( ix = xmap.first(); !ix.last(); ix.next() ) {
      w = xmap.spacegroup().num_symops() / xmap.multiplicity( ix.coord() );
      r = xmap[ix];
      mn += w;
      m2z += w*r*r;
      m3z += w*r*r*r;
      r -= m1z;
      m1m += w*r;
      m2m += w*r*r;
      m3m += w*r*r*r;
      if ( r < mi ) mi = r;
      if ( r > ma ) ma = r;
    }
    m2z = pow( m2z/mn, 0.50000000 );
    m3z = pow( m3z/mn, 0.33333333 );
    m1m = m1m / mn;
    m2m = pow( m2m/mn, 0.50000000 );
    m3m = pow( m3m/mn, 0.33333333 );
    std::cout << "\nMap statistics:\n Number of points in cell: "
	      << clipper::String(mn,20,12) <<
      "\n 1st moment about zero (mean): " << clipper::String(m1z,12,6)
	      << "    About mean       : " << clipper::String(m1m,12,6) <<
      "\n 2nd moment about zero       : " << clipper::String(m2z,12,6)
	      << "    About mean (rmsd): " << clipper::String(m2m,12,6) <<
      "\n 3rd moment about zero       : " << clipper::String(m3z,12,6)
	      << "    About mean (skew): " << clipper::String(m3m,12,6) <<
      "\n Range:  Min " << clipper::String(mi,12,6)
	      << "  Max " << clipper::String(ma,12,6) <<
      "\n";

    if ( statsrad > 0.0 ) {
      // make squared map
      clipper::Xmap<float> lmom1( xmap );
      clipper::Xmap<float> lmom2( xmap );
      for ( ix = lmom2.first(); !ix.last(); ix.next() )
	lmom2[ix] = clipper::Util::sqr( lmom2[ix] );

      // now calculate local mom1, local mom1 squared
      clipper::MapFilterFn_step fn( statsrad );
      clipper::MapFilter_fft<float> fltr( fn, 1.0, clipper::MapFilter_fft<float>::Relative );
      fltr( lmom1, lmom1 );
      fltr( lmom2, lmom2 );

      // calculate std deviation
      for ( ix = lmom1.first(); !ix.last(); ix.next() )
	lmom2[ix] = sqrt( lmom2[ix] - clipper::Util::sqr( lmom1[ix] ) );

      // now get stats of local standard deviation
      clipper::Map_stats m( lmom2 );
      std::cout << "\nLocal map statistics:\n Mean of local RMSD          : " << clipper::String( m.mean(), 12, 6 ) << "   RMSD of local RMSD: " << clipper::String( m.std_dev(), 12, 6 ) << "\n";

      // next do some handidness checks
      double s, sx, sy, sxx, syy, sxy;
      s = sx = sy = sxx = syy = sxy = 0.0;
      for ( ix = lmom1.first(); !ix.last(); ix.next() ) {
	s += 1.0;
	sx += lmom1[ix];
	sy += lmom2[ix];
	sxx += lmom1[ix]*lmom1[ix];
	syy += lmom2[ix]*lmom2[ix];
	sxy += lmom1[ix]*lmom2[ix];
      }
      int n = fphi.num_obs() - 2;
      double r = (s*sxy-sx*sy)/sqrt((s*sxx-sx*sx)*(s*syy-sy*sy));
      double t = r * sqrt( double(n) / ( 1.0 - r*r ) );
      double p = 1.0 - pfunc( t, n );
      std::cout << "\n Local mean/variance correlation : " << r << "\n";
      std::cout << "\n Local mean/variance significance: " << t << " " << n << " " << p << "\n";
    }
  }
}



double sfunc( double q, int i, int j ) {
  double s, t;
  s = t = 1.0;
  for ( int k = i; k <= j; k+= 2 ) {
    t *= q*double(k)/double(k+1);
    s += t;
  }
  return s;
}

double pfunc( double t, int nt ) {
  double w = t / sqrt(double(nt));
  double th = atan( w );
  double c = cos(th);
  double s = sin(th);
  if ( nt%2 == 1 )
    return 1.0-0.63662*(th+s*c*sfunc(c*c,2,nt-3));
  else
    return 1.0-(s*sfunc(c*c,1,nt-3));
}
