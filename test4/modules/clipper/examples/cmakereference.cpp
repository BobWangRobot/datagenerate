// Clipper app to prepare reference structure
/* Copyright 2003-2004 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-cif.h>
#include "ccp4-extras.h"

extern "C" {
#if defined _MSC_VER
 #include <io.h>
#else
 #include <unistd.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "ftplib.h"               // ftp client lib
extern void decompress(int,int);  // decompress lib
}

int main( int argc, char** argv )
{
  CCP4program prog( "cmakereference", "0.1", "$Date: 2004/09/01" );

  std::cout << "\n  This program includes a modified version of ftplib\n  by Thomas Pfau (see http://nbpfaus.net/~pfau/ftplib/),\n  and a portion of the public domain 'ncompress' code.\n  It is distributed under CCP4 part 0 or LGPL.\n\n";

  // defaults
  clipper::String pdbid = "NONE";
  clipper::String pdbfilez = "NONE";
  clipper::String rflfilez = "NONE";
  clipper::String mtzout = "NONE";
  clipper::String pdbout = "NONE";
  clipper::Resolution reso;
  bool bulk = true;
  int n_refln = 1000;
  int n_param = 20;
  int verbose = 0;

  // command input
  CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if ( args[arg] == "-pdbid" ) {
      if ( ++arg < args.size() ) pdbid = args[arg];
    } else if ( args[arg] == "-pdbin" ) {
      if ( ++arg < args.size() ) pdbfilez = args[arg];
    } else if ( args[arg] == "-cifin" ) {
      if ( ++arg < args.size() ) rflfilez = args[arg];
    } else if ( args[arg] == "-mtzout" ) {
      if ( ++arg < args.size() ) mtzout = args[arg];
    } else if ( args[arg] == "-pdbout" ) {
      if ( ++arg < args.size() ) pdbout = args[arg];
    } else if ( args[arg] == "-resolution" ) {
      if ( ++arg < args.size() ) reso = clipper::Resolution( clipper::String(args[arg]).f() );
    } else if ( args[arg] == "-num-reflns" ) {
      if ( ++arg < args.size() ) n_refln = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-num-params" ) {
      if ( ++arg < args.size() ) n_param = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-no-bulk" ) {
      bulk = false;
    } else if ( args[arg] == "-verbose" ) {
      if ( ++arg < args.size() ) verbose = clipper::String(args[arg]).i();
    } else {
      std::cout << "Unrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "Usage: cmakereference\n\t-pdbid <accession-code>\n\t-pdbin <.ent.Z-file>\n\t-cifin <.ent.Z-file>\n\t-resolution <reso>\n\tGenerate reference structure for pirate.\nIf no .ent.Z files are given, they are fetched by ftp if possible.\n";
    exit(1);
  }

  // post input
  if ( mtzout == "NONE" ) mtzout = "reference-"+pdbid+".mtz";
  if ( pdbout == "NONE" ) pdbout = "reference-"+pdbid+".pdb";

  // ftp settings
  // db specific
  std::string host = "ftp.ebi.ac.uk";
  std::string pdbdir = "/pub/databases/rcsb/pdb/data/structures/all/pdb/";
  std::string rfldir = "/pub/databases/rcsb/pdb/data/structures/all/structure_factors/";
  // generic
  clipper::String user = "anonymous";
  clipper::String pass = "user@host";
  clipper::String tmpdir = "/tmp";
  // from environment
  const char *ptrusr = getenv("USER");
  const char *ptrhst = getenv("HOST");
  const char *ptrscr = getenv("CCP4_SCR");
  if ( ptrusr != NULL && ptrhst != NULL )
    pass = std::string(ptrusr) + "@" + std::string(ptrhst);
  if ( ptrscr != NULL )
    tmpdir = std::string(ptrscr);

  // messages
  clipper::Message_fatal noserver( "Unable to connect to server. Try manual ftp and give filenames." );
  clipper::Message_fatal nologin( "Unable to login to server. Try manual ftp and give filenames." );
  clipper::Message_fatal nopdb( "Unable to fetch coordinate file.\nCheck coordinates are available. Try manual ftp and give filenames." );
  clipper::Message_fatal norfl( "Unable to fetch reflection file.\nCheck reflections are available. Try manual ftp and give filenames." );
  clipper::Message_fatal nordc( "Unable to read compressed file." );
  clipper::Message_fatal nowrc( "Unable to write uncompressed file." );

  // ftp the files, if required
  if ( pdbfilez == "NONE" || rflfilez == "NONE" ) {

    // provide instructions
    std::cout << "\nAttempting to connect to EBI/MSD for file download. If this step fails, use \nyour preferred ftp client to fetch the files, then re-run giving the \nfilenames as arguments, or use the EBI/MSD website.\n\nftp " << host << "\n" << user << "\n" << pass << "\nbinary\ncd " << pdbdir << "\nget pdb" << pdbid << ".ent.Z\ncd " << rfldir << "\nget r" << pdbid << "sf.ent.Z\nquit\n\n";

    netbuf *pbuf;
    int err;
    FtpInit();

    std::cout << "Connecting...\n\n";
    err = FtpConnect( host.c_str(), &pbuf );
    if ( err != 1 ) clipper::Message::message( noserver );
    err = FtpLogin( user.c_str(), pass.c_str(), pbuf );
    if ( err != 1 ) clipper::Message::message( nologin );

    if ( pdbfilez == "NONE" ) {
      std::cout << "Fetching coordinates...\n\n";
      pdbfilez = tmpdir + "/pdb" + pdbid + ".ent.Z";
      err = FtpChdir( pdbdir.c_str(), pbuf );
      if ( err != 1 ) clipper::Message::message( nopdb );
      err = FtpGet(pdbfilez.c_str(),pdbfilez.tail().c_str(),FTPLIB_IMAGE,pbuf);
      if ( err != 1 ) clipper::Message::message( nopdb );
    }

    if ( rflfilez == "NONE" ) {
      std::cout << "Fetching reflections...\n\n";
      rflfilez = tmpdir + "/r" + pdbid + "sf.ent.Z";
      err = FtpChdir( rfldir.c_str(), pbuf );
      if ( err != 1 ) clipper::Message::message( norfl );
      err = FtpGet(rflfilez.c_str(),rflfilez.tail().c_str(),FTPLIB_IMAGE,pbuf);
      if ( err != 1 ) clipper::Message::message( norfl );
    }

    FtpQuit( pbuf );
  }

  // now uncompress the files
  int fdip, fdop;

  std::cout << "Decompressing coordinates...\n\n";
  clipper::String pdbfile = pdbfilez.substr( 0, pdbfilez.length() - 2 );
  fdip = open( pdbfilez.c_str(), O_RDONLY );
  if ( fdip < 0 ) clipper::Message::message( nordc );
  fdop = open( pdbfile.c_str(), O_CREAT|O_WRONLY|O_TRUNC, S_IREAD|S_IWRITE );
  if ( fdop < 0 ) clipper::Message::message( nowrc );
  decompress( fdip, fdop );
  close( fdip );
  close( fdop );

  std::cout << "Decompressing reflections...\n\n";
  clipper::String rflfile = rflfilez.substr( 0, rflfilez.length() - 2 );
  fdip = open( rflfilez.c_str(), O_RDONLY );
  if ( fdip < 0 ) clipper::Message::message( nordc );
  fdop = open( rflfile.c_str(), O_CREAT|O_WRONLY|O_TRUNC, S_IREAD|S_IWRITE );
  if ( fdop < 0 ) clipper::Message::message( nowrc );
  decompress( fdip, fdop );
  close( fdip );
  close( fdop );

  // make data objects
  clipper::CIFfile cifin;
  clipper::CCP4MTZfile mtzfile;
  clipper::HKL_info hkls_in;
  double bulkfrc, bulkscl;
  typedef clipper::HKL_data_base::HKL_reference_index HRI;

  // atomic model
  clipper::MMDBManager mmdb;
  mmdb.SetFlag( MMDBF_AutoSerials | MMDBF_IgnoreDuplSeqNum );
  mmdb.ReadPDBASCII( (char*)pdbfile.c_str() );

  // read reflection info
  cifin.open_read( rflfile );
  clipper::Spacegroup   spgr = mmdb.spacegroup();
  clipper::Cell         cell = mmdb.cell();
  if ( reso.is_null() ) reso = cifin.resolution( cell );
  clipper::HKL_info hkls( spgr, cell, reso, true );
  std::cout << "Number of reflections: " << hkls.num_reflections() << "\n";
  clipper::HKL_data<clipper::data32::F_sigF> fo( hkls );
  cifin.import_hkl_data( fo );
  cifin.close_read();

  // get a list of all the atoms
  clipper::mmdb::PPCAtom psel;
  int hndl, nsel;
  hndl = mmdb.NewSelection();
  mmdb.SelectAtoms( hndl, 0, 0, SKEY_NEW );
  mmdb.GetSelIndex( hndl, psel, nsel );
  clipper::MMDBAtom_list atoms( psel, nsel );
  mmdb.DeleteSelection( hndl );

  // calculate structure factors
  clipper::HKL_data<clipper::data32::F_phi> fc( hkls );
  if ( bulk ) {
    clipper::SFcalc_obs_bulk<float> sfcb;
    sfcb( fc, fo, atoms );
    bulkfrc = sfcb.bulk_frac();
    bulkscl = sfcb.bulk_scale();
  } else {
    clipper::SFcalc_aniso_fft<float> sfc;
    sfc( fc, atoms );
    bulkfrc = bulkscl = 0.0;
  }

  // now do sigmaa calc
  clipper::HKL_data<clipper::data32::F_phi> fb( hkls ), fd( hkls );
  clipper::HKL_data<clipper::data32::Phi_fom> phiw( hkls );
  clipper::HKL_data<clipper::data32::Flag> flag( hkls );
  for ( HRI ih = flag.first(); !ih.last(); ih.next() )
    if ( !fo[ih].missing() )
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
  mtzfile.open_write( mtzout );
  mtzfile.export_hkl_info( hkls );
  mtzfile.export_hkl_data( fo,   "/*/*/FP" );
  mtzfile.export_hkl_data( abcd, "/*/*/FC" );
  mtzfile.export_hkl_data( fb,   "/*/*/FC_BEST" );
  mtzfile.export_hkl_data( fd,   "/*/*/FC_DIFF" );
  mtzfile.close_write();
  mmdb.WritePDBASCII( (char*)pdbout.c_str() );

  // now calc R and R-free
  std::vector<double> params( n_param, 1.0 );
  clipper::BasisFn_spline basisfn( fo, n_param, 1.0 );
  clipper::TargetFn_scaleF1F2<clipper::data32::F_phi,clipper::data32::F_sigF> targetfn( fc, fo );
  clipper::ResolutionFn rfn( hkls, basisfn, targetfn, params );
  double r1w, f1w, r1f, f1f, Fo, Fc;
  r1w = f1w = r1f = f1f = 0.0;
  for ( HRI ih = fo.first(); !ih.last(); ih.next() )
    if ( !fo[ih].missing() ) {
      Fo = fo[ih].f();
      Fc = sqrt( clipper::Util::max( rfn.f(ih), 0.0 ) ) * fc[ih].f();
      r1w += fabs( Fo - Fc );
      f1w += Fo;
    }
  r1f /= clipper::Util::max( f1f, 0.1 );
  r1w /= clipper::Util::max( f1w, 0.1 );
  std::cout << "\n R-factor      : " << r1w << "\n";
}
