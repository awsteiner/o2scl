/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2021, Andrew W. Steiner

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

  -------------------------------------------------------------------
*/
#include <o2scl/test_mgr.h>
#include <o2scl/hdf_file.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

int main(void) {

  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;

  cout.setf(ios::scientific);
  cout.precision(10);

  test_mgr t;
  t.set_output_level(1);

  system("echo -n 1: ; date");
  
  if (true) {
    
    // This section now checks that gets() now works with
    // both fixed- and variable-length strings
    cout << "Test gets() for fixed- and variable-length strings." << endl;

    hdf_file hf;
    hf.open_or_create("hdf_file_gets.o2");
    hf.sets("str","string");
    hf.sets_fixed("strf","fixed-length string");
    hf.sets("strb","s");
    hf.sets_fixed("strbf","f");
    hf.close();

    hf.open("hdf_file_gets.o2");
    string stmp;
    hf.gets("strb",stmp);
    t.test_gen(stmp=="s","string 1");
    hf.gets("strbf",stmp);
    t.test_gen(stmp=="f","string 2");
    hf.gets("str",stmp);
    t.test_gen(stmp=="string","string 3");
    hf.gets("strf",stmp);
    t.test_gen(stmp=="fixed-length string","string 4");
    hf.close();
    cout << endl;
  }
  
  // Uncompressed section
  {
  
    cout << "Test generic scalars: " << endl;
    {
    
      double d=sqrt(2.0);
      char c='z';
      int i=31;
      float f=sqrt(2.0);
      size_t u=41;
      string s="Te st.";
      string sb="Te st.";
    
      char c2;
      double d2;
      int i2;
      float f2;
      size_t u2;
      string s2;
      string sb2;
    
      hdf_file hf;
    
      // Open a file, set the scalar values
    
      hf.open_or_create("hdf_file.o2");
      hf.setc("testc",c);
      hf.setd("testd",d);
      hf.set_szt("testu",u);
      hf.setf("testf",f);

      hid_t group_id=hf.open_group("Integers");
      hid_t file_id=hf.get_file_id();
      hf.set_current_id(group_id);
      hf.seti("testi",i);
      hf.set_current_id(file_id);
      hf.close_group(group_id);

      hf.sets("tests",s);
      hf.sets_fixed("testsb",sb);
      hf.close();

      // Re-open the file, change the string values

      hf.open("hdf_file.o2",true);
      s="Test two.";
      hf.sets("tests",s);
      hf.sets_fixed("testsb","Texst.");
      hf.close();

      // Re-open the file, get the scalar values
    
      hf.open("hdf_file.o2",true);
      hf.getc("testc",c2);
      hf.getd("testd",d2);
      hf.get_szt("testu",u2);
      hf.getf("testf",f2);

      hid_t group_id2=hf.open_group("Integers");
      hid_t file_id2=hf.get_file_id();
      hf.set_current_id(group_id2);
      hf.geti("testi",i2);
      hf.set_current_id(file_id2);
      hf.close_group(group_id2);

      hf.gets("tests",s2);
      hf.gets_fixed("testsb",sb2);
      hf.close();

      // Compare the old and new values
    
      t.test_gen(c2=='z',"character");
      t.test_rel(d2,sqrt(2.0),1.0e-12,"double");
      t.test_rel(f2,sqrt(2.0f),1.0e-7f,"float");
      t.test_gen(i2==31,"integer");
      t.test_gen(u2==41,"size_t");
      t.test_gen(s2=="Test two.","string");
      t.test_gen(sb2=="Texst.","string fixed");
    
    }
    cout << endl;

    cout << "Test fixed arrays: " << endl;
    {
    
      char c[4]={'z','y','a','b'};
      double d[6]={sqrt(2.0),2.0,4.0,16.0};
      int i[4]={31,41,59,26};
      float f[4]={sqrtf(2.0),2.0,4.0,16.0};
    
      char *c2=new char[4];
      double *d2=new double[4];
      int *i2=new int[4];
      float *f2=new float[4];
    
      hdf_file hf;
    
      // Open a file, set the vector values

      hf.open("hdf_file.o2",true);
      hf.setc_arr_fixed("testca",4,c);
      hf.setd_arr_fixed("testda",4,d);
      hf.setf_arr_fixed("testfa",4,f);

      hid_t group_id=hf.open_group("Integers");
      hid_t file_id=hf.get_file_id();
      hf.set_current_id(group_id);
      hf.seti_arr_fixed("testia",4,i);
      hf.set_current_id(file_id);
      hf.close_group(group_id);

      hf.close();

      // Re-open the file, get the scalar values

      hf.open("hdf_file.o2");
      hf.getc_arr("testca",4,c2);
      hf.getd_arr("testda",4,d2);
      hf.getf_arr("testfa",4,f2);

      hid_t group_id2=hf.open_group("Integers");
      hid_t file_id2=hf.get_file_id();
      hf.set_current_id(group_id2);
      hf.geti_arr("testia",4,i2);
      hf.set_current_id(file_id2);
      hf.close_group(group_id2);

      hf.close();

      // Compare the old and new values
    
      t.test_gen_vec<char *>(4,c,c2,"character");
      t.test_rel_vec<double *>(4,d,d2,1.0e-12,"double");
      t.test_rel_vec<float *>(4,f,f2,1.0e-7f,"float");
      t.test_gen_vec<int *>(4,i,i2,"int");

      delete[] c2;
      delete[] d2;
      delete[] f2;
      delete[] i2;
    
    }
    cout << endl;

    cout << "Test unlimited arrays: " << endl;
    {
    
      char c[6]={'z','y','a','b','c','d'};
      double d[6]={sqrt(2.0),2.0,4.0,16.0,32.0,64.0};
      int i[6]={31,41,59,26,53,58};
      float f[6]={sqrtf(2.0),2.0,4.0,16.0,32.0,64.0};
      vector<string> vs;
      vs.push_back("This");
      vs.push_back("is");
      vs.push_back("a");
      vs.push_back("test.");
    
      char *c2=new char[6];
      double *d2=new double[6];
      int *i2=new int[6];
      float *f2=new float[6];
      vector<string> vs2;
    
      hdf_file hf;
    
      // Open a file, set the vector values

      hf.open("hdf_file.o2",true);
      hf.setc_arr("testca2",4,c);
      hf.setd_arr("testda2",4,d);
      hf.setf_arr("testfa2",4,f);

      hid_t group_id=hf.open_group("Integers");
      hid_t file_id=hf.get_file_id();
      hf.set_current_id(group_id);
      hf.seti_arr("testia2",4,i);
      hf.set_current_id(file_id);
      hf.close_group(group_id);

      hf.sets_vec("testsa2",vs);

      hf.close();

      // Extend the vectors

      hf.open("hdf_file.o2",true);
      hf.setc_arr("testca2",6,c);
      hf.setd_arr("testda2",6,d);
      hf.setf_arr("testfa2",6,f);

      group_id=hf.open_group("Integers");
      file_id=hf.get_file_id();
      hf.set_current_id(group_id);
      hf.seti_arr("testia2",6,i);
      hf.set_current_id(file_id);
      hf.close_group(group_id);

      vs.push_back("Another");
      vs.push_back("test.");
      hf.sets_vec("testsa2",vs);

      hf.close();

      // Re-open the file, get the scalar values

      hf.open("hdf_file.o2");
      hf.getc_arr("testca2",6,c2);
      hf.getd_arr("testda2",6,d2);
      hf.getf_arr("testfa2",6,f2);

      hid_t group_id2=hf.open_group("Integers");
      hid_t file_id2=hf.get_file_id();
      hf.set_current_id(group_id2);
      hf.geti_arr("testia2",6,i2);
      hf.set_current_id(file_id2);
      hf.close_group(group_id2);

      hf.gets_vec("testsa2",vs2);

      hf.close();

      // Compare the old and new values

      t.test_gen_vec<char *>(6,c,c2,"character");
      t.test_rel_vec<double *>(6,d,d2,1.0e-12,"double");
      t.test_rel_vec<float *>(6,f,f2,1.0e-7f,"float");
      t.test_gen_vec<int *>(6,i,i2,"int");

      for(size_t j=0;j<vs.size();j++) {
	t.test_gen(vs[j]==vs2[j],"string");
      }

      delete[] c2;
      delete[] d2;
      delete[] f2;
      delete[] i2;
    
    }
    cout << endl;

    cout << "Test abstract vectors: " << endl;
    {
      std::vector<double> v1a;
      ubvector v2a(4);
      ubvector v3a(4);
      //ubvector_size_t v4a(4);
      for(size_t i=0;i<4;i++) {
	v1a.push_back(2*i);
	v2a[i]=2*i;
	v3a[i]=2*i;
	//v4a[i]=3*i;
      }

      // Create the vector datasets
    
      hdf_file hf;
      hf.open("hdf_file.o2",true);

      hf.setd_vec("vec1",v1a);
      hf.setd_vec_copy("vec2",v2a);
      hf.setd_vec_copy("vec3",v3a);
      //hf.set_szt_vec("vec4",v4a);

      hf.close();

      // Larger vectors

      std::vector<double> v1b;
      ubvector v2b(6);
      ubvector v3b(6);
      //ubvector_size_t v4b(6);
      for(size_t i=0;i<6;i++) {
	v1b.push_back(2*i);
	v2b[i]=2*i;
	v3b[i]=2*i;
	//v4b[i]=3*i;
      }

      // Extend the vector datasets
    
      hf.open("hdf_file.o2",true);

      hf.setd_vec("vec1",v1b);
      hf.setd_vec_copy("vec2",v2b);
      hf.setd_vec_copy("vec3",v3b);
      //hf.set_szt_vec("vec4",v4b);

      hf.close();

      // Get vector data and test

      std::vector<double> v1c;
      ubvector v2c;
      ubvector v3c;
      //ubvector_size_t v4c;

      hf.open("hdf_file.o2");

      hf.getd_vec("vec1",v1c);
      hf.getd_vec_copy("vec2",v2c);
      hf.getd_vec_copy("vec3",v3c);
      //hf.get_szt_vec("vec4",v4c);

      hf.close();
    
      t.test_rel_vec(6,v1b,v1c,1.0e-12,"vector 1");
      t.test_rel_vec(6,v2b,v2c,1.0e-12,"vector 2");
      t.test_rel_vec(6,v3b,v3c,1.0e-12,"vector 3");
      //t.test_rel_vec(6,v4b,v4c,1.0e-12,"vector 4");
    
    }
    cout << endl;

    cout << "Test abstract matrices: " << endl;
    {
      ubmatrix m1a(4,4);
      ubmatrix m2a(4,4);
      for(size_t i=0;i<4;i++) {
	for(size_t j=0;j<4;j++) {
	  m1a(i,j)=((double)(i))-((double)j);
	  m2a(i,j)=((double)(i))-((double)j);
	}
      }

      // Create the vector datasets
    
      hdf_file hf;
      hf.open("hdf_file.o2",true);

      hf.setd_mat_copy("mat1",m1a);
      hf.setd_mat_copy("mat2",m2a);

      hf.close();

      // Larger vectors

      ubmatrix m1b(6,6);
      ubmatrix m2b(6,6);
      for(size_t i=0;i<6;i++) {
	for(size_t j=0;j<6;j++) {
	  m1b(i,j)=((double)(i))-((double)j);
	  m2b(i,j)=((double)(i))-((double)j);
	}
      }

      // Extend the vector datasets
    
      hf.open("hdf_file.o2",true);

      hf.setd_mat_copy("mat1",m1b);
      hf.setd_mat_copy("mat2",m2b);

      hf.close();

      // Get matrix data and test

      ubmatrix m1c;
      ubmatrix m2c;

      hf.open("hdf_file.o2");

      hf.getd_mat_copy("mat1",m1c);
      hf.getd_mat_copy("mat2",m2c);

      hf.close();
    
      t.test_rel_mat(6,6,m1b,m1c,1.0e-12,"matrix 1");
      t.test_rel_mat(6,6,m2b,m2c,1.0e-12,"matrix 2");

    }
    cout << endl;

  }

#ifdef O2SCL_HDF5_COMP

  {
  
    cout << "Test generic scalars: " << endl;
    {
    
      double d=sqrt(2.0);
      char c='z';
      int i=31;
      float f=sqrt(2.0);
      size_t u=41;
      string s="Te st.";
      string sb="Te st.";
    
      char c2;
      double d2;
      int i2;
      float f2;
      size_t u2;
      string s2;
      string sb2;
    
      hdf_file hf;
      hf.compr_type=1;
    
      // Open a file, set the scalar values
    
      hf.open_or_create("hdf_file_comp.o2");
      hf.setc("testc",c);
      hf.setd("testd",d);
      hf.set_szt("testu",u);
      hf.setf("testf",f);

      hid_t group_id=hf.open_group("Integers");
      hid_t file_id=hf.get_file_id();
      hf.set_current_id(group_id);
      hf.seti("testi",i);
      hf.set_current_id(file_id);
      hf.close_group(group_id);

      hf.sets("tests",s);
      hf.sets_fixed("testsb",sb);
      hf.close();

      // Re-open the file, change the string values

      hf.open("hdf_file_comp.o2",true);
      s="Test two.";
      hf.sets("tests",s);
      hf.sets_fixed("testsb","Texst.");
      hf.close();

      // Re-open the file, get the scalar values
    
      hf.open("hdf_file_comp.o2",true);
      hf.getc("testc",c2);
      hf.getd("testd",d2);
      hf.get_szt("testu",u2);
      hf.getf("testf",f2);

      hid_t group_id2=hf.open_group("Integers");
      hid_t file_id2=hf.get_file_id();
      hf.set_current_id(group_id2);
      hf.geti("testi",i2);
      hf.set_current_id(file_id2);
      hf.close_group(group_id2);

      hf.gets("tests",s2);
      hf.gets_fixed("testsb",sb2);
      hf.close();

      // Compare the old and new values
    
      t.test_gen(c2=='z',"character");
      t.test_rel(d2,sqrt(2.0),1.0e-12,"double");
      t.test_rel(f2,sqrt(2.0f),1.0e-7f,"float");
      t.test_gen(i2==31,"integer");
      t.test_gen(u2==41,"size_t");
      t.test_gen(s2=="Test two.","string");
      t.test_gen(sb2=="Texst.","string fixed");
    
    }
    cout << endl;

    cout << "Test fixed arrays: " << endl;
    {
    
      char c[4]={'z','y','a','b'};
      double d[6]={sqrt(2.0),2.0,4.0,16.0};
      int i[4]={31,41,59,26};
      float f[4]={sqrtf(2.0),2.0,4.0,16.0};
    
      char *c2=new char[4];
      double *d2=new double[4];
      int *i2=new int[4];
      float *f2=new float[4];
    
      hdf_file hf;
      hf.compr_type=1;
    
      // Open a file, set the vector values

      hf.open("hdf_file_comp.o2",true);
      hf.setc_arr_fixed("testca",4,c);
      hf.setd_arr_fixed("testda",4,d);
      hf.setf_arr_fixed("testfa",4,f);

      hid_t group_id=hf.open_group("Integers");
      hid_t file_id=hf.get_file_id();
      hf.set_current_id(group_id);
      hf.seti_arr_fixed("testia",4,i);
      hf.set_current_id(file_id);
      hf.close_group(group_id);

      hf.close();

      // Re-open the file, get the scalar values

      hf.open("hdf_file_comp.o2");
      hf.getc_arr("testca",4,c2);
      hf.getd_arr("testda",4,d2);
      hf.getf_arr("testfa",4,f2);

      hid_t group_id2=hf.open_group("Integers");
      hid_t file_id2=hf.get_file_id();
      hf.set_current_id(group_id2);
      hf.geti_arr("testia",4,i2);
      hf.set_current_id(file_id2);
      hf.close_group(group_id2);

      hf.close();

      // Compare the old and new values
    
      t.test_gen_vec<char *>(4,c,c2,"character");
      t.test_rel_vec<double *>(4,d,d2,1.0e-12,"double");
      t.test_rel_vec<float *>(4,f,f2,1.0e-7f,"float");
      t.test_gen_vec<int *>(4,i,i2,"int");

      delete[] c2;
      delete[] d2;
      delete[] f2;
      delete[] i2;
    
    }
    cout << endl;

    cout << "Test unlimited arrays: " << endl;
    {
    
      char c[6]={'z','y','a','b','c','d'};
      double d[6]={sqrt(2.0),2.0,4.0,16.0,32.0,64.0};
      int i[6]={31,41,59,26,53,58};
      float f[6]={sqrtf(2.0),2.0,4.0,16.0,32.0,64.0};
      vector<string> vs;
      vs.push_back("This");
      vs.push_back("is");
      vs.push_back("a");
      vs.push_back("test.");
    
      char *c2=new char[6];
      double *d2=new double[6];
      int *i2=new int[6];
      float *f2=new float[6];
      vector<string> vs2;
    
      hdf_file hf;
      hf.compr_type=1;
    
      // Open a file, set the vector values

      hf.open("hdf_file_comp.o2",true);
      hf.setc_arr("testca2",4,c);
      hf.setd_arr("testda2",4,d);
      hf.setf_arr("testfa2",4,f);

      hid_t group_id=hf.open_group("Integers");
      hid_t file_id=hf.get_file_id();
      hf.set_current_id(group_id);
      hf.seti_arr("testia2",4,i);
      hf.set_current_id(file_id);
      hf.close_group(group_id);

      hf.sets_vec("testsa2",vs);

      hf.close();

      // Extend the vectors

      hf.open("hdf_file_comp.o2",true);
      hf.setc_arr("testca2",6,c);
      hf.setd_arr("testda2",6,d);
      hf.setf_arr("testfa2",6,f);

      group_id=hf.open_group("Integers");
      file_id=hf.get_file_id();
      hf.set_current_id(group_id);
      hf.seti_arr("testia2",6,i);
      hf.set_current_id(file_id);
      hf.close_group(group_id);

      vs.push_back("Another");
      vs.push_back("test.");
      hf.sets_vec("testsa2",vs);

      hf.close();

      // Re-open the file, get the scalar values

      hf.open("hdf_file_comp.o2");
      hf.getc_arr("testca2",6,c2);
      hf.getd_arr("testda2",6,d2);
      hf.getf_arr("testfa2",6,f2);

      hid_t group_id2=hf.open_group("Integers");
      hid_t file_id2=hf.get_file_id();
      hf.set_current_id(group_id2);
      hf.geti_arr("testia2",6,i2);
      hf.set_current_id(file_id2);
      hf.close_group(group_id2);

      hf.gets_vec("testsa2",vs2);

      hf.close();

      // Compare the old and new values

      t.test_gen_vec<char *>(6,c,c2,"character");
      t.test_rel_vec<double *>(6,d,d2,1.0e-12,"double");
      t.test_rel_vec<float *>(6,f,f2,1.0e-7f,"float");
      t.test_gen_vec<int *>(6,i,i2,"int");

      for(size_t j=0;j<vs.size();j++) {
	t.test_gen(vs[j]==vs2[j],"string");
      }

      delete[] c2;
      delete[] d2;
      delete[] f2;
      delete[] i2;
    
    }
    cout << endl;

    cout << "Test abstract vectors: " << endl;
    {
      std::vector<double> v1a;
      ubvector v2a(4);
      ubvector v3a(4);
      //ubvector_size_t v4a(4);
      for(size_t i=0;i<4;i++) {
	v1a.push_back(2*i);
	v2a[i]=2*i;
	v3a[i]=2*i;
	//v4a[i]=3*i;
      }

      // Create the vector datasets
    
      hdf_file hf;
      hf.compr_type=1;
      hf.open("hdf_file_comp.o2",true);

      hf.setd_vec("vec1",v1a);
      hf.setd_vec_copy("vec2",v2a);
      hf.setd_vec_copy("vec3",v3a);
      //hf.set_szt_vec("vec4",v4a);

      hf.close();

      // Larger vectors

      std::vector<double> v1b;
      ubvector v2b(6);
      ubvector v3b(6);
      //ubvector_size_t v4b(6);
      for(size_t i=0;i<6;i++) {
	v1b.push_back(2*i);
	v2b[i]=2*i;
	v3b[i]=2*i;
	//v4b[i]=3*i;
      }

      // Extend the vector datasets
    
      hf.open("hdf_file_comp.o2",true);

      hf.setd_vec("vec1",v1b);
      hf.setd_vec_copy("vec2",v2b);
      hf.setd_vec_copy("vec3",v3b);
      //hf.set_szt_vec("vec4",v4b);

      hf.close();

      // Get vector data and test

      std::vector<double> v1c;
      ubvector v2c;
      ubvector v3c;
      //ubvector_size_t v4c;

      hf.open("hdf_file_comp.o2");

      hf.getd_vec("vec1",v1c);
      hf.getd_vec_copy("vec2",v2c);
      hf.getd_vec_copy("vec3",v3c);
      //hf.get_szt_vec("vec4",v4c);

      hf.close();
    
      t.test_rel_vec(6,v1b,v1c,1.0e-12,"vector 1");
      t.test_rel_vec(6,v2b,v2c,1.0e-12,"vector 2");
      t.test_rel_vec(6,v3b,v3c,1.0e-12,"vector 3");
      //t.test_rel_vec(6,v4b,v4c,1.0e-12,"vector 4");
    
    }
    cout << endl;

    cout << "Test abstract matrices: " << endl;
    {
      ubmatrix m1a(4,4);
      ubmatrix m2a(4,4);
      for(size_t i=0;i<4;i++) {
	for(size_t j=0;j<4;j++) {
	  m1a(i,j)=((double)(i))-((double)j);
	  m2a(i,j)=((double)(i))-((double)j);
	}
      }

      // Create the vector datasets
    
      hdf_file hf;
      hf.compr_type=1;
      hf.open("hdf_file_comp.o2",true);

      hf.setd_mat_copy("mat1",m1a);
      hf.setd_mat_copy("mat2",m2a);

      hf.close();

      // Larger vectors

      ubmatrix m1b(6,6);
      ubmatrix m2b(6,6);
      for(size_t i=0;i<6;i++) {
	for(size_t j=0;j<6;j++) {
	  m1b(i,j)=((double)(i))-((double)j);
	  m2b(i,j)=((double)(i))-((double)j);
	}
      }

      // Extend the vector datasets
    
      hf.open("hdf_file_comp.o2",true);

      hf.setd_mat_copy("mat1",m1b);
      hf.setd_mat_copy("mat2",m2b);

      hf.close();

      // Get matrix data and test

      ubmatrix m1c;
      ubmatrix m2c;

      hf.open("hdf_file_comp.o2");

      hf.getd_mat_copy("mat1",m1c);
      hf.getd_mat_copy("mat2",m2c);

      hf.close();
    
      t.test_rel_mat(6,6,m1b,m1c,1.0e-12,"matrix 1");
      t.test_rel_mat(6,6,m2b,m2c,1.0e-12,"matrix 2");

    }
    cout << endl;

  }

#endif
  
  system("echo -n 2: ; date");
  
  t.report();

  return 0;
}

