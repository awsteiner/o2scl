/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2006-2026, Andrew W. Steiner

  This file is part of O2scl.

  O2scl is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  O2scl is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with O2scl. If not, see <http://www.gnu.org/licenses/>.

  ───────────────────────────────────────────────────────────────────
*/
/*
  This code generates the O2scl HDF files for the
  Brussels-Skyrme-on-a-Grid (BSkG) mass models.

  BSkG3 is parsed from the raw data file distributed on the
  BRUSLIB website, http://www.astro.ulb.ac.be/bruslib/nucdata/ ,
  named 'bskg03-dat' there. See [Grams23]_ for the BSkG3
  reference. This file fills in the extra BSkG-specific fields
  (from \c gamma to \c I3) which the older HFB tables leave
  unused.

  BSkG1, BSkG2, and BSkG4 are parsed from the tables bundled with
  Jerome Margueron's 'nucleardatapy' Python toolkit
  (https://github.com/jeromemargueron/nucleardatapy , files
  '2021-BSkG1.txt', '2022-BSkG2.txt', and '2025-BSkG4.txt' in
  'nucleardatapy/data/nuclei/masses/Theory'), since the BRUSLIB
  site currently only hosts the BSkG3 table. These files use a
  different column layout than the BRUSLIB BSkG3 file, so they
  fill in a different set of extra fields (from \c Erot to \c
  par_n) instead. BSkG1 and BSkG2 do not provide \c beta30 or \c
  beta32. See [Scamps21]_, [Ryssens22]_, and [Grams24]_ for the
  BSkG1, BSkG2, and BSkG4 references, respectively.

  In all cases, the result is stored in an object of type \ref
  o2scl::nucmass_hfb_sp::entry, reusing the same class as the
  HFB17-HFB27 tables (see \ref o2scl_hdf::hfb_sp_load() ). None of
  the BSkG tables provide the 'Deformation and Wigner energies'
  quantity, \c def_wig, so that field is left at zero for all BSkG
  entries.
*/
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// for exit()
#include <cstdlib>

#include <o2scl/string_conv.h>
#include <o2scl/nucmass.h>
#include <o2scl/nucmass_hfb.h>
#include <o2scl/hdf_file.h>

#include <hdf5_hl.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

/** \brief Return \c true if \c x is close to one of the sentinel
    values which the BRUSLIB tables use to mark an unavailable or
    undefined quantity (e.g. 999.99, 99.99, or 99.9, with either
    sign)
*/
bool is_blank(double x, double sentinel) {
  return (fabs(fabs(x)-sentinel)<1.0e-6);
}

/** \brief Return \c true if the raw token \c s is one of the
    sentinel strings ("-" or a numerical value close to \c
    sentinel, with either sign) which the nucleardatapy BSkG1,
    BSkG2, and BSkG4 tables use to mark an unavailable or undefined
    quantity
*/
bool is_blank_tok(std::string s, double sentinel) {
  if (s=="-") return true;
  return is_blank(o2scl::stod(s),sentinel);
}

/** \brief Parse a single-character or signed-integer parity token
    ("+", "-", "+1", "-1", or "0" for an undefined parity) into +1,
    -1, or 99 (blank)
*/
int parse_parity(std::string s) {
  if (s=="+" || s=="+1") return 1;
  if (s=="-" || s=="-1") return -1;
  if (s=="0") return 99;
  O2SCL_ERR("Could not parse parity token in parse_parity().",
            o2scl::exc_einval);
  return 0;
}

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  if (argc<3) {
    cout << "Usage: bskg_parse <dir> <out file>." << endl;
    exit(-1);
  }

  string dir=argv[1];
  string out_fname=argv[2];
  string orig_file;
  string reference;

  // BSkG1, BSkG2, and BSkG4 are parsed from nucleardatapy's tables
  // rather than from BRUSLIB (see the file-level comment above),
  // and fill in a different set of extra fields than BSkG3.
  if (out_fname=="bskg1.o2" || out_fname=="bskg2.o2" ||
      out_fname=="bskg4.o2") {

    bool has_beta3=false;
    if (out_fname=="bskg1.o2") {
      orig_file="2021-BSkG1.txt";
      reference=((string)"G. Scamps, S. Goriely, E. Olsen, ")+
        "M. Bender, and W. Ryssens, Eur. Phys. J. A 57 (2021) 333.";
    } else if (out_fname=="bskg2.o2") {
      orig_file="2022-BSkG2.txt";
      reference=((string)"W. Ryssens, G. Scamps, S. Goriely, ")+
        "and M. Bender, Eur. Phys. J. A 58 (2022) 246.";
    } else {
      orig_file="2025-BSkG4.txt";
      has_beta3=true;
      reference=((string)"G. Grams, W. Ryssens, N. Shchechilin, ")+
        "A. Sanchez-Fernandez, N. Chamel, and S. Goriely, "+
        "arXiv:2411.08007 (2024).";
    }

    string in_fname=dir+"/"+orig_file;
    cout << "Opening file '" << in_fname << "'." << endl;
    ifstream fin(in_fname.c_str());
    if (!fin.is_open()) {
      cout << "Failed to open '" << in_fname << "'." << endl;
      exit(-1);
    }

    // Skip the single header line
    string stemp;
    getline(fin,stemp);
    cout << stemp << endl;

    size_t ncols_exp=has_beta3 ? 20 : 18;

    vector<nucmass_hfb_sp::entry> list;

    while (getline(fin,stemp)) {

      vector<string> sv;
      split_string(stemp,sv);
      if (sv.size()==0) continue;
      if (sv.size()!=ncols_exp) {
        cout << "Expected " << ncols_exp << " columns, got " << sv.size()
             << " in line: " << stemp << endl;
        exit(-1);
      }

      nucmass_hfb_sp::entry he;

      he.Z=o2scl::stoi(sv[0]);
      he.N=o2scl::stoi(sv[1]);
      he.A=he.Z+he.N;

      // Not provided by the nucleardatapy BSkG1/BSkG2/BSkG4 tables;
      // left at zero, as is done for def_wig
      he.def_wig=0.0;
      he.Sn=0.0;
      he.Sp=0.0;
      he.Qbet=0.0;
      he.Jexp=0.0;
      he.Jth=0.0;
      he.Pexp=99;
      he.Pth=99;
      he.gamma=0.0;
      he.S2n=0.0;
      he.S2p=0.0;
      he.delta3n=0.0;
      he.delta3p=0.0;
      he.delta5n=0.0;
      he.delta5p=0.0;
      he.rc4=0.0;
      he.I1=0.0;
      he.I2=0.0;
      he.I3=0.0;

      he.Mcal=o2scl::stod(sv[3]);
      he.Err=is_blank_tok(sv[4],99.9) ? 1.0e99 : o2scl::stod(sv[4]);
      he.beta20=o2scl::stod(sv[6]);
      he.beta22=o2scl::stod(sv[7]);
      // The total beta2 deformation in column 8 (index sv[8]) is
      // stored in the same field used for the "beta 2 deformation"
      // in the older HFB tables
      he.bet2=o2scl::stod(sv[8]);
      // Not provided by any of the BSkG1/BSkG2/BSkG4 tables
      he.bet4=0.0;

      size_t ix=9;
      if (has_beta3) {
        he.beta30=o2scl::stod(sv[9]);
        he.beta32=o2scl::stod(sv[10]);
        ix=11;
      } else {
        // Not provided by BSkG1 or BSkG2
        he.beta30=0.0;
        he.beta32=0.0;
      }

      he.Erot=o2scl::stod(sv[ix]);
      he.avgap_n=o2scl::stod(sv[ix+1]);
      he.avgap_p=o2scl::stod(sv[ix+2]);
      he.Rch=o2scl::stod(sv[ix+3]);
      he.rc_exp=is_blank_tok(sv[ix+4],99.9) ? 1.0e99 :
        o2scl::stod(sv[ix+4]);
      he.rc_err=is_blank_tok(sv[ix+5],99.9) ? 1.0e99 :
        o2scl::stod(sv[ix+5]);
      he.MOI=o2scl::stod(sv[ix+6]);
      he.par_p=parse_parity(sv[ix+7]);
      he.par_n=parse_parity(sv[ix+8]);

      if (list.size()==0) {
        cout << "First line: " << endl;
        cout << he.Z << " " << he.N << " " << he.A << endl;
        cout << "\t" << he.Mcal << " " << he.Err << endl;
        cout << "\t" << he.beta20 << " " << he.beta22 << " "
             << he.bet2 << endl;
        cout << "\t" << he.Erot << " " << he.avgap_n << " "
             << he.avgap_p << endl;
        cout << "\t" << he.Rch << " " << he.rc_exp << " "
             << he.rc_err << " " << he.MOI << endl;
        cout << "\t" << he.par_p << " " << he.par_n << endl;
      }

      list.push_back(he);
    }
    fin.close();

    cout << "Read " << list.size() << " records." << endl;

    // Make HDF table
    size_t offset[38]={HOFFSET(nucmass_hfb_sp::entry,N),
                       HOFFSET(nucmass_hfb_sp::entry,Z),
                       HOFFSET(nucmass_hfb_sp::entry,A),
                       HOFFSET(nucmass_hfb_sp::entry,bet2),
                       HOFFSET(nucmass_hfb_sp::entry,bet4),
                       HOFFSET(nucmass_hfb_sp::entry,Rch),
                       HOFFSET(nucmass_hfb_sp::entry,Sn),
                       HOFFSET(nucmass_hfb_sp::entry,Sp),
                       HOFFSET(nucmass_hfb_sp::entry,Qbet),
                       HOFFSET(nucmass_hfb_sp::entry,Mcal),
                       HOFFSET(nucmass_hfb_sp::entry,Err),
                       HOFFSET(nucmass_hfb_sp::entry,Jexp),
                       HOFFSET(nucmass_hfb_sp::entry,Jth),
                       HOFFSET(nucmass_hfb_sp::entry,Pexp),
                       HOFFSET(nucmass_hfb_sp::entry,Pth),
                       HOFFSET(nucmass_hfb_sp::entry,gamma),
                       HOFFSET(nucmass_hfb_sp::entry,beta20),
                       HOFFSET(nucmass_hfb_sp::entry,beta22),
                       HOFFSET(nucmass_hfb_sp::entry,beta30),
                       HOFFSET(nucmass_hfb_sp::entry,beta32),
                       HOFFSET(nucmass_hfb_sp::entry,S2n),
                       HOFFSET(nucmass_hfb_sp::entry,S2p),
                       HOFFSET(nucmass_hfb_sp::entry,delta3n),
                       HOFFSET(nucmass_hfb_sp::entry,delta3p),
                       HOFFSET(nucmass_hfb_sp::entry,delta5n),
                       HOFFSET(nucmass_hfb_sp::entry,delta5p),
                       HOFFSET(nucmass_hfb_sp::entry,rc4),
                       HOFFSET(nucmass_hfb_sp::entry,I1),
                       HOFFSET(nucmass_hfb_sp::entry,I2),
                       HOFFSET(nucmass_hfb_sp::entry,I3),
                       HOFFSET(nucmass_hfb_sp::entry,Erot),
                       HOFFSET(nucmass_hfb_sp::entry,avgap_n),
                       HOFFSET(nucmass_hfb_sp::entry,avgap_p),
                       HOFFSET(nucmass_hfb_sp::entry,rc_exp),
                       HOFFSET(nucmass_hfb_sp::entry,rc_err),
                       HOFFSET(nucmass_hfb_sp::entry,MOI),
                       HOFFSET(nucmass_hfb_sp::entry,par_p),
                       HOFFSET(nucmass_hfb_sp::entry,par_n)};

    nucmass_hfb_sp::entry he;

    size_t sizes[38]={sizeof(he.N),
                      sizeof(he.Z),
                      sizeof(he.A),
                      sizeof(he.bet2),
                      sizeof(he.bet4),
                      sizeof(he.Rch),
                      sizeof(he.Sn),
                      sizeof(he.Sp),
                      sizeof(he.Qbet),
                      sizeof(he.Mcal),
                      sizeof(he.Err),
                      sizeof(he.Jexp),
                      sizeof(he.Jth),
                      sizeof(he.Pexp),
                      sizeof(he.Pth),
                      sizeof(he.gamma),
                      sizeof(he.beta20),
                      sizeof(he.beta22),
                      sizeof(he.beta30),
                      sizeof(he.beta32),
                      sizeof(he.S2n),
                      sizeof(he.S2p),
                      sizeof(he.delta3n),
                      sizeof(he.delta3p),
                      sizeof(he.delta5n),
                      sizeof(he.delta5p),
                      sizeof(he.rc4),
                      sizeof(he.I1),
                      sizeof(he.I2),
                      sizeof(he.I3),
                      sizeof(he.Erot),
                      sizeof(he.avgap_n),
                      sizeof(he.avgap_p),
                      sizeof(he.rc_exp),
                      sizeof(he.rc_err),
                      sizeof(he.MOI),
                      sizeof(he.par_p),
                      sizeof(he.par_n)};

    const char *names[38]={
      "Neutron number",
      "Proton number",
      "Atomic number",
      "Beta 2 deformation",
      "Beta 4 deformation",
      "RMS charge radius",
      "Neutron separation energy",
      "Proton separation energy",
      "Beta-decay energy",
      "Calculated mass excess",
      "Error between experimental and calculated mass excess",
      "Experimental spin","Theoretical spin",
      "Experimental parity","Theoretical parity",
      "Triaxial deformation angle gamma",
      "Axial quadrupole deformation beta_20",
      "Non-axial quadrupole deformation beta_22",
      "Axial octupole deformation beta_30",
      "Non-axial octupole deformation beta_32",
      "Two-neutron separation energy",
      "Two-proton separation energy",
      "Three-point neutron odd-even mass staggering",
      "Three-point proton odd-even mass staggering",
      "Five-point neutron odd-even mass staggering",
      "Five-point proton odd-even mass staggering",
      "Fourth radial moment of the charge density to the 1/4 power",
      "Moment of inertia about the first axis",
      "Moment of inertia about the second axis",
      "Moment of inertia about the third axis",
      "Rotational correction energy",
      "Average neutron pairing gap",
      "Average proton pairing gap",
      "Experimental RMS charge radius",
      "Error between RMS charge radius and experimental value",
      "Moment of inertia",
      "Parity of the proton subsystem",
      "Parity of the neutron subsystem"};

    hid_t field_type[38]={
      H5T_NATIVE_INT,     // N
      H5T_NATIVE_INT,     // Z
      H5T_NATIVE_INT,     // A
      H5T_NATIVE_DOUBLE,  // bet2
      H5T_NATIVE_DOUBLE,  // bet4
      H5T_NATIVE_DOUBLE,  // Rch
      H5T_NATIVE_DOUBLE,  // Sn
      H5T_NATIVE_DOUBLE,  // Sp
      H5T_NATIVE_DOUBLE,  // Qbet
      H5T_NATIVE_DOUBLE,  // Mcal
      H5T_NATIVE_DOUBLE,  // Err
      H5T_NATIVE_DOUBLE,  // Jexp
      H5T_NATIVE_DOUBLE,  // Jth
      H5T_NATIVE_INT,     // Pexp
      H5T_NATIVE_INT,     // Pth
      H5T_NATIVE_DOUBLE,  // gamma
      H5T_NATIVE_DOUBLE,  // beta20
      H5T_NATIVE_DOUBLE,  // beta22
      H5T_NATIVE_DOUBLE,  // beta30
      H5T_NATIVE_DOUBLE,  // beta32
      H5T_NATIVE_DOUBLE,  // S2n
      H5T_NATIVE_DOUBLE,  // S2p
      H5T_NATIVE_DOUBLE,  // delta3n
      H5T_NATIVE_DOUBLE,  // delta3p
      H5T_NATIVE_DOUBLE,  // delta5n
      H5T_NATIVE_DOUBLE,  // delta5p
      H5T_NATIVE_DOUBLE,  // rc4
      H5T_NATIVE_DOUBLE,  // I1
      H5T_NATIVE_DOUBLE,  // I2
      H5T_NATIVE_DOUBLE,  // I3
      H5T_NATIVE_DOUBLE,  // Erot
      H5T_NATIVE_DOUBLE,  // avgap_n
      H5T_NATIVE_DOUBLE,  // avgap_p
      H5T_NATIVE_DOUBLE,  // rc_exp
      H5T_NATIVE_DOUBLE,  // rc_err
      H5T_NATIVE_DOUBLE,  // MOI
      H5T_NATIVE_INT,     // par_p
      H5T_NATIVE_INT};    // par_n

    // Remove any existing file first, since open_or_create() opens
    // an existing file for read/write rather than truncating it,
    // and also records write access properly (unlike calling
    // set_current_id() directly on a handle from a bare
    // H5Fcreate() call)
    remove(out_fname.c_str());

    hdf_file hf;
    hf.open_or_create(out_fname);
    hid_t file=hf.get_current_id();

    hf.seti("nrecords",list.size());
    cout << "nrecords: " << list.size() << endl;
    hf.sets_fixed("comment",
                  ((string)"HDF5 version of BSkG ")+
                  "mass data created for O2scl. "
                  "See https://awsteiner.org/code/o2scl for details.");

    herr_t status=H5TBmake_table
      (orig_file.c_str(),file,out_fname.c_str(),
       38,list.size(),sizeof(nucmass_hfb_sp::entry),
       names,offset,field_type,100,0,0,&list[0]);

    hf.sets("orig_file",orig_file);
    hf.sets("reference",reference);

    hf.close();

    return 0;
  }

  if (out_fname=="bskg3.o2") {
    orig_file="bskg03-dat";
    reference=((string)"G. Grams, W. Ryssens, G. Scamps, ")+
      "S. Goriely, and N. Chamel, Eur. Phys. J. A 59 (2023) 270.";
  } else {
    cout << "Unknown output file '" << out_fname << "'." << endl;
    exit(-1);
  }

  string in_fname=dir+"/"+orig_file;
  cout << "Opening file '" << in_fname << "'." << endl;
  ifstream fin(in_fname.c_str());
  if (!fin.is_open()) {
    cout << "Failed to open '" << in_fname << "'." << endl;
    exit(-1);
  }

  // Skip the three header/comment lines
  string stemp;
  for(size_t j=0;j<3;j++) {
    getline(fin,stemp);
    cout << stemp << endl;
  }

  vector<nucmass_hfb_sp::entry> list;

  while (getline(fin,stemp)) {

    vector<string> sv;
    split_string(stemp,sv);
    if (sv.size()==0) continue;
    if (sv.size()!=29) {
      cout << "Expected 29 columns, got " << sv.size()
           << " in line: " << stemp << endl;
      exit(-1);
    }

    nucmass_hfb_sp::entry he;

    he.Z=o2scl::stoi(sv[0]);
    he.A=o2scl::stoi(sv[1]);
    he.N=he.A-he.Z;

    he.bet2=o2scl::stod(sv[2]);
    he.bet4=o2scl::stod(sv[3]);
    he.Rch=o2scl::stod(sv[4]);
    he.gamma=o2scl::stod(sv[5]);

    he.Sn=o2scl::stod(sv[6]);
    if (is_blank(he.Sn,999.99)) he.Sn=1.0e99;
    he.Sp=o2scl::stod(sv[7]);
    if (is_blank(he.Sp,999.99)) he.Sp=1.0e99;
    he.Qbet=o2scl::stod(sv[8]);
    if (is_blank(he.Qbet,999.99)) he.Qbet=1.0e99;

    he.Mcal=o2scl::stod(sv[9]);
    he.Err=o2scl::stod(sv[10]);
    if (is_blank(he.Err,99.99)) he.Err=1.0e99;

    he.Jexp=o2scl::stod(sv[11]);
    if (is_blank(he.Jexp,99.9)) he.Jexp=1.0e99;
    he.Jth=o2scl::stod(sv[12]);

    he.Pexp=o2scl::stoi(sv[13]);
    if (he.Pexp==9) he.Pexp=99;
    he.Pth=o2scl::stoi(sv[14]);

    // Not provided by the BSkG3 (BRUSLIB) table; left at zero, as
    // is done for def_wig
    he.def_wig=0.0;
    he.Erot=0.0;
    he.avgap_n=0.0;
    he.avgap_p=0.0;
    he.rc_exp=0.0;
    he.rc_err=0.0;
    he.MOI=0.0;
    he.par_p=0;
    he.par_n=0;

    he.beta20=o2scl::stod(sv[15]);
    he.beta22=o2scl::stod(sv[16]);
    he.beta30=o2scl::stod(sv[17]);
    he.beta32=o2scl::stod(sv[18]);

    he.S2n=o2scl::stod(sv[19]);
    if (is_blank(he.S2n,99.99)) he.S2n=1.0e99;
    he.S2p=o2scl::stod(sv[20]);
    if (is_blank(he.S2p,99.99)) he.S2p=1.0e99;

    he.delta3n=o2scl::stod(sv[21]);
    if (is_blank(he.delta3n,99.99)) he.delta3n=1.0e99;
    he.delta3p=o2scl::stod(sv[22]);
    if (is_blank(he.delta3p,99.99)) he.delta3p=1.0e99;
    he.delta5n=o2scl::stod(sv[23]);
    if (is_blank(he.delta5n,99.99)) he.delta5n=1.0e99;
    he.delta5p=o2scl::stod(sv[24]);
    if (is_blank(he.delta5p,99.99)) he.delta5p=1.0e99;

    he.rc4=o2scl::stod(sv[25]);
    he.I1=o2scl::stod(sv[26]);
    he.I2=o2scl::stod(sv[27]);
    he.I3=o2scl::stod(sv[28]);

    if (list.size()==0) {
      cout << "First line: " << endl;
      cout << he.Z << " " << he.A << " " << he.N << endl;
      cout << "\t" << he.bet2 << " " << he.bet4 << " " << he.Rch
           << " " << he.gamma << endl;
      cout << "\t" << he.Sn << " " << he.Sp << " " << he.Qbet << endl;
      cout << "\t" << he.Mcal << " " << he.Err << endl;
      cout << "\t" << he.Jexp << " " << he.Jth << " "
           << he.Pexp << " " << he.Pth << endl;
    }

    list.push_back(he);
  }
  fin.close();

  cout << "Read " << list.size() << " records." << endl;

  // Make HDF table
  size_t offset[30]={HOFFSET(nucmass_hfb_sp::entry,N),
                     HOFFSET(nucmass_hfb_sp::entry,Z),
                     HOFFSET(nucmass_hfb_sp::entry,A),
                     HOFFSET(nucmass_hfb_sp::entry,bet2),
                     HOFFSET(nucmass_hfb_sp::entry,bet4),
                     HOFFSET(nucmass_hfb_sp::entry,Rch),
                     HOFFSET(nucmass_hfb_sp::entry,Sn),
                     HOFFSET(nucmass_hfb_sp::entry,Sp),
                     HOFFSET(nucmass_hfb_sp::entry,Qbet),
                     HOFFSET(nucmass_hfb_sp::entry,Mcal),
                     HOFFSET(nucmass_hfb_sp::entry,Err),
                     HOFFSET(nucmass_hfb_sp::entry,Jexp),
                     HOFFSET(nucmass_hfb_sp::entry,Jth),
                     HOFFSET(nucmass_hfb_sp::entry,Pexp),
                     HOFFSET(nucmass_hfb_sp::entry,Pth),
                     HOFFSET(nucmass_hfb_sp::entry,gamma),
                     HOFFSET(nucmass_hfb_sp::entry,beta20),
                     HOFFSET(nucmass_hfb_sp::entry,beta22),
                     HOFFSET(nucmass_hfb_sp::entry,beta30),
                     HOFFSET(nucmass_hfb_sp::entry,beta32),
                     HOFFSET(nucmass_hfb_sp::entry,S2n),
                     HOFFSET(nucmass_hfb_sp::entry,S2p),
                     HOFFSET(nucmass_hfb_sp::entry,delta3n),
                     HOFFSET(nucmass_hfb_sp::entry,delta3p),
                     HOFFSET(nucmass_hfb_sp::entry,delta5n),
                     HOFFSET(nucmass_hfb_sp::entry,delta5p),
                     HOFFSET(nucmass_hfb_sp::entry,rc4),
                     HOFFSET(nucmass_hfb_sp::entry,I1),
                     HOFFSET(nucmass_hfb_sp::entry,I2),
                     HOFFSET(nucmass_hfb_sp::entry,I3)};

  nucmass_hfb_sp::entry he;

  size_t sizes[30]={sizeof(he.N),
                    sizeof(he.Z),
                    sizeof(he.A),
                    sizeof(he.bet2),
                    sizeof(he.bet4),
                    sizeof(he.Rch),
                    sizeof(he.Sn),
                    sizeof(he.Sp),
                    sizeof(he.Qbet),
                    sizeof(he.Mcal),
                    sizeof(he.Err),
                    sizeof(he.Jexp),
                    sizeof(he.Jth),
                    sizeof(he.Pexp),
                    sizeof(he.Pth),
                    sizeof(he.gamma),
                    sizeof(he.beta20),
                    sizeof(he.beta22),
                    sizeof(he.beta30),
                    sizeof(he.beta32),
                    sizeof(he.S2n),
                    sizeof(he.S2p),
                    sizeof(he.delta3n),
                    sizeof(he.delta3p),
                    sizeof(he.delta5n),
                    sizeof(he.delta5p),
                    sizeof(he.rc4),
                    sizeof(he.I1),
                    sizeof(he.I2),
                    sizeof(he.I3)};

  const char *names[30]={
    "Neutron number",
    "Proton number",
    "Atomic number",
    "Beta 2 deformation",
    "Beta 4 deformation",
    "RMS charge radius",
    "Neutron separation energy",
    "Proton separation energy",
    "Beta-decay energy",
    "Calculated mass excess",
    "Error between experimental and calculated mass excess",
    "Experimental spin","Theoretical spin",
    "Experimental parity","Theoretical parity",
    "Triaxial deformation angle gamma",
    "Axial quadrupole deformation beta_20",
    "Non-axial quadrupole deformation beta_22",
    "Axial octupole deformation beta_30",
    "Non-axial octupole deformation beta_32",
    "Two-neutron separation energy",
    "Two-proton separation energy",
    "Three-point neutron odd-even mass staggering",
    "Three-point proton odd-even mass staggering",
    "Five-point neutron odd-even mass staggering",
    "Five-point proton odd-even mass staggering",
    "Fourth radial moment of the charge density to the 1/4 power",
    "Moment of inertia about the first axis",
    "Moment of inertia about the second axis",
    "Moment of inertia about the third axis"};

  hid_t field_type[30]={
    H5T_NATIVE_INT,     // N
    H5T_NATIVE_INT,     // Z
    H5T_NATIVE_INT,     // A
    H5T_NATIVE_DOUBLE,  // bet2
    H5T_NATIVE_DOUBLE,  // bet4
    H5T_NATIVE_DOUBLE,  // Rch
    H5T_NATIVE_DOUBLE,  // Sn
    H5T_NATIVE_DOUBLE,  // Sp
    H5T_NATIVE_DOUBLE,  // Qbet
    H5T_NATIVE_DOUBLE,  // Mcal
    H5T_NATIVE_DOUBLE,  // Err
    H5T_NATIVE_DOUBLE,  // Jexp
    H5T_NATIVE_DOUBLE,  // Jth
    H5T_NATIVE_INT,     // Pexp
    H5T_NATIVE_INT,     // Pth
    H5T_NATIVE_DOUBLE,  // gamma
    H5T_NATIVE_DOUBLE,  // beta20
    H5T_NATIVE_DOUBLE,  // beta22
    H5T_NATIVE_DOUBLE,  // beta30
    H5T_NATIVE_DOUBLE,  // beta32
    H5T_NATIVE_DOUBLE,  // S2n
    H5T_NATIVE_DOUBLE,  // S2p
    H5T_NATIVE_DOUBLE,  // delta3n
    H5T_NATIVE_DOUBLE,  // delta3p
    H5T_NATIVE_DOUBLE,  // delta5n
    H5T_NATIVE_DOUBLE,  // delta5p
    H5T_NATIVE_DOUBLE,  // rc4
    H5T_NATIVE_DOUBLE,  // I1
    H5T_NATIVE_DOUBLE,  // I2
    H5T_NATIVE_DOUBLE}; // I3

  // Remove any existing file first, since open_or_create() opens
  // an existing file for read/write rather than truncating it, and
  // also records write access properly (unlike calling
  // set_current_id() directly on a handle from a bare H5Fcreate()
  // call)
  remove(out_fname.c_str());

  hdf_file hf;
  hf.open_or_create(out_fname);
  hid_t file=hf.get_current_id();

  hf.seti("nrecords",list.size());
  cout << "nrecords: " << list.size() << endl;
  hf.sets_fixed("comment",
                ((string)"HDF5 version of BSkG ")+
                "mass data created for O2scl. "
                "See https://awsteiner.org/code/o2scl for details.");

  herr_t status=H5TBmake_table
    (orig_file.c_str(),file,out_fname.c_str(),
     30,list.size(),sizeof(nucmass_hfb_sp::entry),
     names,offset,field_type,100,0,0,&list[0]);

  hf.sets("orig_file",orig_file);
  hf.sets("reference",reference);

  hf.close();

  return 0;
}
