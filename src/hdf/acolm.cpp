/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
#include "acolm.h"

#include <o2scl/cloud_file.h>
#include <o2scl/vector_derint.h>
#include <o2scl/cursesw.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_acol;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

/*
  This source file contains the constructor and the functions
  command_add(), get_input(), get_input_one(), command_del(),
  clear_obj(), run(), setup_parameters(), setup_help(), setup_cli(),
  and setup_options().
*/

acol_manager::acol_manager() : cset(this,&acol_manager::comm_set),
			       cng(o2scl_settings.get_convert_units()) {
  
  ff.unicode_mode();
  
  obj_name="acol";
  verbose=1;
  pretty=true;
  names_out=true;
  use_regex=false;
  scientific=true;
  precision=6;
  def_args="";
  ncols=0;
  
  post_interactive=false;

  type="";

  env_var_name="ACOL_DEFAULTS";
  interp_type=1;

#ifdef O2SCL_HDF5_COMP
  compress=1;
#else
  compress=0;
#endif

  cng.err_on_fail=false;

  type_list.push_back("table");
  type_list.push_back("table3d");
  type_list.push_back("hist");
  type_list.push_back("hist_2d");
  type_list.push_back("vector<contour_line>");
  type_list.push_back("int");
  type_list.push_back("double");
  type_list.push_back("char");
  type_list.push_back("string");
  type_list.push_back("int[]");
  type_list.push_back("double[]");
  type_list.push_back("string[]");
  type_list.push_back("size_t");
  type_list.push_back("size_t[]");
  type_list.push_back("uniform_grid<double>");
  type_list.push_back("tensor_grid");
  type_list.push_back("tensor");
  type_list.push_back("tensor<int>");
  type_list.push_back("tensor<size_t>");
  type_list.push_back("prob_dens_mdim_amr");
  type_list.push_back("prob_dens_mdim_gaussian");
  type_list.push_back("prob_dens_mdim_gmm");
  type_list.push_back("vec_vec_string");
  type_list.push_back("vec_vec_double");
  vector_sort<vector<string>,string>(type_list.size(),type_list);
  
  {
    vector<std::string> itmp={"value"};
    type_comm_list.insert(std::make_pair("int",itmp));
    type_comm_list.insert(std::make_pair("double",itmp));
    type_comm_list.insert(std::make_pair("char",itmp));
    type_comm_list.insert(std::make_pair("size_t",itmp));
    type_comm_list.insert(std::make_pair("string",itmp));
  }
  {
    vector<std::string> itmp={"add-vec","average-rows",
      "assign","cat","convert-unit",
      "correl","delete-col",
      "delete-rows","delete-rows-tol",
      "deriv","deriv2","value-grid",
      "find-row","fit","function",
      "get-row","get-unit","value","index",
      "insert","insert-full","integ","interp","interp-table3d",
      "list","max","min","nlines","refine","rename",
      "select","select-rows",
      "set-data","set-unit","sort","stats","sum",
      "to-hist","to-hist-2d","to-table3d","wstats",
      "ser-hist-t3d","to-gaussian","to-pdma","to-gmm"
    };
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("table",itmp));
  }
  {
    vector<std::string> itmp={"cat","contours","deriv-x","deriv-y",
      "function","value","value-grid","get-grid",
      "insert","interp","refine","stats","select",
      "list","max","min","rename","set-data","to-hist",
      "slice","slice-hist","sum","to-hist-2d","to-table",
      "to-tensor-grid","to-tg-fermi","x-name","y-name"};
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("table3d",itmp));
  }
  {
    vector<std::string> itmp={"list","min","max","to-table3d",
      "rearrange"};
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("tensor<int>",itmp));
    type_comm_list.insert(std::make_pair("tensor<size_t>",itmp));
  }
  {
    vector<std::string> itmp={"list","diag","to-table3d","to-table3d-sum",
      "max","min","to-tensor-grid","rearrange",
      "value","function","sum","stats","deriv"};
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("tensor",itmp));
  }
  {
    vector<std::string> itmp={"to-table3d","sample"};
    type_comm_list.insert(std::make_pair("prob_dens_mdim_amr",itmp));
  }
  {
    vector<std::string> itmp={"sample"};
    type_comm_list.insert(std::make_pair("prob_dens_mdim_gaussian",itmp));
  }
  {
    vector<std::string> itmp={"list"};
    type_comm_list.insert(std::make_pair("vector<contour_line>",itmp));
  }
  {
    vector<std::string> itmp={"list"};
    type_comm_list.insert(std::make_pair("vec_vec_double",itmp));
  }
  {
    vector<std::string> itmp={"list"};
    type_comm_list.insert(std::make_pair("vec_vec_string",itmp));
  }
  {
    vector<std::string> itmp={"sample","list"};
    type_comm_list.insert(std::make_pair("prob_dens_mdim_gmm",itmp));
  }
  {
    vector<std::string> itmp={"list","to-table3d","slice","to-table",
      "set-grid","max","min","rearrange",
      "get-grid","interp","value","to-tensor",
      "value-grid","function","sum","stats",
      "binary","deriv"};
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("tensor_grid",itmp));
  }
  {
    vector<std::string> itmp={"max","min","contours","list","to-table3d",
    "refine"};
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("hist_2d",itmp));
  }
  {
    vector<std::string> itmp={"to-table","function","list"};
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("hist",itmp));
  }
  {
    vector<std::string> itmp={"deriv","interp","max","min","sort",
      "to-table","function","sum","find","value","resize"};
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("double[]",itmp));
    type_comm_list.insert(std::make_pair("int[]",itmp));
    type_comm_list.insert(std::make_pair("size_t[]",itmp));
  }
  {
    vector<std::string> itmp={"value","resize"};
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("string[]",itmp));
  }

  // Ensure the RNGs for different types are somewhat uncorrelated
  rng.clock_seed();
  rng_ld.set_seed(rng.get_seed()*2);

  terminal ter;
  command_color=ter.color_from_int(ter.c_cyan+ter.int_high);
  type_color=ter.color_from_int(ter.c_magenta+ter.int_high);
  param_color=ter.color_from_int(ter.c_red+ter.int_high);
  help_color=ter.color_from_int(ter.c_green+ter.int_high);
  exec_color=ter.color_from_int(ter.c_white+ter.int_high);
  url_color=ter.color_from_int(ter.att_underline);
  default_color=ter.default_fgbg();
}

void acol_manager::color_replacements(std::string &s) {
  string_replace(s,"[c]",command_color);
  string_replace(s,"[d]",default_color);
  string_replace(s,"[e]",exec_color);
  string_replace(s,"[h]",help_color);
  string_replace(s,"[p]",param_color);
  string_replace(s,"[t]",type_color);
  string_replace(s,"[u]",url_color);
  return;
}

bool acol_manager::help_found(std::string arg1, std::string arg2) {

  if (arg2=="") {

    // Base commands from cli class
    if (arg1=="alias" || arg1=="commands" || arg1=="get" ||
        arg1=="help" || arg1=="license" || arg1=="no-intro" ||
        arg1=="run" || arg1=="shell" || arg1=="set" ||
        arg1=="warranty") {
      return true;
    }
    
    for(size_t i=0;i<type_list.size();i++) {
      if (cl->string_equal_dash(arg1,type_list[i])) {
        return true;
      }
    }
    for(size_t i=0;i<help_doc_strings.size();i++) {
      if (cl->string_equal_dash(arg1,help_doc_strings[i][0])) {
        return true;
      }
    }
    for(size_t i=0;i<param_doc_strings.size();i++) {
      if (cl->string_equal_dash(arg1,param_doc_strings[i][0])) {
        return true;
      }
    }
    for(size_t i=0;i<cmd_doc_strings.size();i++) {
      if (cl->string_equal_dash(arg1,cmd_doc_strings[i][0])) {
        return true;
      }
    }
  } else {
    for(size_t i=0;i<cmd_doc_strings.size();i++) {
      if (cl->string_equal_dash(arg2,cmd_doc_strings[i][0])) {
        return true;
      }
    }
  }
  
  return false;
}

void acol_manager::update_o2_docs(size_t narr,
                                  o2scl::comm_option_s *options_arr,
                                  std::string new_type) {

  int loc_verbose=0;

  if (verbose>2 || loc_verbose>1) {
    if (new_type.length()>0) {
      cout << "acol_manager::update_o2_docs(): Updating documentation "
           << "for type " << new_type << " with " << narr
           << " options " << endl;
    } else {
      cout << "acol_manager::update_o2_docs(): Updating documentation "
           << "(no current object) with " << narr << " options." << endl;
    }
  }
  
  // Loop through all options, attempting to assign docs for all
  for(size_t j=0;j<narr;j++) {
    
    bool found=false;

    // Loop over all command doc strings, attempting to find a match
    for(size_t k=0;k<cmd_doc_strings.size() && found==false;k++) {

      // We found a match for option j and list of doc strings k
      if (cmd_doc_strings[k][0]==options_arr[j].lng) {
        
        if (loc_verbose>1) {
          cout << "Found documentation for " << options_arr[j].lng << endl;
        }
        found=true;

        // Only proceed if there is at least one line of documentation
        if (cmd_doc_strings[k].size()>=2) {
          
          options_arr[j].desc=cmd_doc_strings[k][1];
          if (loc_verbose>1) {
            cout << "Found brief desc.: " << options_arr[j].desc << endl;
          }

          // If there might be some detailed documentation
          if (cmd_doc_strings[k].size()>=3) {
            
            if (loc_verbose>1) {
              cout << "Found detailed desc." << endl;
            }
            bool generic_docs=true;
            bool no_obj_docs=false;
            
            if (cmd_doc_strings[k][2].substr(0,19)==
                ((string)"For objects of type")) {
              generic_docs=false;
            } else if (cmd_doc_strings[k][2].substr(0,30)==
                       ((string)"If there is no current object:")) {
              generic_docs=false;
              no_obj_docs=true;
            }
            
            if (generic_docs) {
              if (new_type.length()==0) {
                if (loc_verbose>1) {
                  cout << "Found generic docs, no type." << endl;
                }
              } else {
                if (loc_verbose>1) {
                  cout << "Found generic docs, type is " << new_type << endl;
                }
              }
            } else {
              if (new_type.length()==0) {
                if (loc_verbose>1) {
                  cout << "No generic docs, no type." << endl;
                }
              } else {
                if (loc_verbose>1) {
                  cout << "No generic docs, type is " << new_type << endl;
                }
              }
            }
            if (loc_verbose>1 && new_type.length()==0) {
              if (no_obj_docs) {
                cout << "Found no-object docs." << endl;
              } else {
                cout << "No no-object docs." << endl;
              }
            }
            
            if (new_type.length()==0) {
              
              if (no_obj_docs) {
                
                if (loc_verbose>1) {
                  cout << "Reading no-object docs: " << endl;
                }
                
                options_arr[j].desc=cmd_doc_strings[k][3];
                options_arr[j].parm_desc=cmd_doc_strings[k][4];
                bool loop_done=false;
                if (cmd_doc_strings[k].size()>=6) {
                  options_arr[j].help=cmd_doc_strings[k][5];
                }
                for(size_t kk=6;kk<cmd_doc_strings[k].size() &&
                      loop_done==false;kk++) {
                  if (cmd_doc_strings[k][kk].substr(0,19)==
                      ((string)"For objects of type")) {
                    loop_done=true;
                  } else if (cmd_doc_strings[k][kk].substr(0,30)==
                             ((string)"If there is no current object:")) {
                    loop_done=true;
                  }
                  options_arr[j].help+="\n\n"+cmd_doc_strings[k][kk];
                }
                
              } else if (generic_docs==true) {
                
                if (loc_verbose>1) {
                  cout << "Reading generic docs: " << endl;
                }
                options_arr[j].parm_desc=cmd_doc_strings[k][2];
                bool loop_done=false;
                if (cmd_doc_strings[k].size()>=4) {
                  options_arr[j].help=cmd_doc_strings[k][3];
                }
                for(size_t kk=4;kk<cmd_doc_strings[k].size() &&
                      loop_done==false;kk++) {
                  if (cmd_doc_strings[k][kk].substr(0,19)==
                      ((string)"For objects of type")) {
                    loop_done=true;
                  } else if (cmd_doc_strings[k][kk].substr(0,30)==
                             ((string)"If there is no current object:")) {
                    loop_done=true;
                  }
                  options_arr[j].help+="\n\n"+cmd_doc_strings[k][kk];
                }
              } else {
                if (loc_verbose>1) {
                  cout << "No generic docs, no type." << endl;
                }
                options_arr[j].parm_desc="";
                options_arr[j].help="";
              }
              
            } else {
              
              if (loc_verbose>1) {
                cout << "Current type is '" << new_type << "'." << endl;
              }
              options_arr[j].desc="";
              options_arr[j].parm_desc="";
              options_arr[j].help="";
              
              bool loop1_done=false;
              for(size_t kk=2;kk<cmd_doc_strings[k].size() &&
                    loop1_done==false;kk++) {
                
                string s="For objects of type "+new_type+":";
                // (When types contain <>'s, the xml parser manges them
                // and ends up adding a space.)
                string s2="For objects of type "+new_type+" :";
                string s3="For objects of type [t]"+new_type+"[d]:";

                if (cmd_doc_strings[k][kk].substr(0,s.length())==s ||
                    cmd_doc_strings[k][kk].substr(0,s2.length())==s2 ||
                    cmd_doc_strings[k][kk].substr(0,s3.length())==s3) {
                  
                  if (loc_verbose>1) {
                    cout << "Found type-specific docs for type " << new_type
                         << "." << endl;
                  }
                  loop1_done=true;
                  bool loop2_done=false;

                  bool short_set=false;
                  bool pd_set=false;
                  bool long_set=false;
                  for(size_t kl=kk+1;kl<cmd_doc_strings[k].size() &&
                        loop2_done==false;kl++) {
                    
                    if (cmd_doc_strings[k][kl].substr(0,19)==
                        ((string)"For objects of type")) {
                      loop2_done=true;
                    } else if (cmd_doc_strings[k][kl].substr(0,30)==
                               ((string)"If there is no current object:")) {
                      loop2_done=true;
                    } else if (short_set==false) {
                      if (loc_verbose>2) {
                        cout << "Adding desc: "
                             << cmd_doc_strings[k][kl] << endl;
                      }
                      options_arr[j].desc=cmd_doc_strings[k][kl];
                      short_set=true;
                    } else if (pd_set==false) {
                      if (loc_verbose>2) {
                        cout << "Adding parm_desc: "
                             << cmd_doc_strings[k][kl] << endl;
                      }
                      options_arr[j].parm_desc=cmd_doc_strings[k][kl];
                      pd_set=true;
                    } else if (long_set==false) {
                      if (loc_verbose>2) {
                        cout << "Adding first doc string: \n  "
                             << cmd_doc_strings[k][kl] << endl;
                      }
                      options_arr[j].help=cmd_doc_strings[k][kl];
                      long_set=true;
                    } else {
                      if (loc_verbose>2) {
                        cout << "Adding doc string: \n  "
                             << cmd_doc_strings[k][kl] << endl;
                      }
                      options_arr[j].help+="\n\n"+cmd_doc_strings[k][kl];
                    }
                  }
                }
              }
            }
          }
        }
        found=true;
      }

      // End of 'k' loop over command doc strings
    }

    // Afterwards, strip the prefix "Parameters: " from the
    // parameter description
    if (options_arr[j].parm_desc.substr(0,11)==((string)"Arguments: ")) {
      options_arr[j].parm_desc=options_arr[j].parm_desc.substr
        (11,options_arr[j].parm_desc.length()-11);
    }
    
    if (found==true) {
      if (verbose>2 || loc_verbose>1) {
        cout << "acol_manager::update_o2_docs() "
             << "found documentation for command "
             << options_arr[j].lng << " ." << endl;
      }
    } else {
      cout << "acol_manager::update_o2_docs() could not "
           << "find documentation for command "
           << options_arr[j].lng << " ." << endl;
    }
    if (loc_verbose>1) {
      cout << endl;
    }

    color_replacements(options_arr[j].desc);
    color_replacements(options_arr[j].parm_desc);
    color_replacements(options_arr[j].help);
    
    // End of 'j' loop over list of options
  }
    
  return;
}

void acol_manager::command_add(std::string new_type) {

  const int both=cli::comm_option_both;
  
  terminal ter;

  // For now, for each type we just verify that the
  // type_comm_list has the same size as the array here
  
  if (new_type=="int") {
    static const size_t narr=1;
    comm_option_s options_arr[narr]=
      {
        {0,"value","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_value),both}
      };
    if (narr!=type_comm_list["int"].size()) {
      O2SCL_ERR("Type comm list does not match for int",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
  } else if (new_type=="double") {
    static const size_t narr=1;
    comm_option_s options_arr[narr]=
      {
        {0,"value","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_value),both}
      };
    if (narr!=type_comm_list["double"].size()) {
      O2SCL_ERR("Type comm list does not match for double",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
  } else if (new_type=="char") {
    static const size_t narr=1;
    comm_option_s options_arr[narr]=
      {
        {0,"value","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_value),both}
      };
    if (narr!=type_comm_list["char"].size()) {
      O2SCL_ERR("Type comm list does not match for char",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
  } else if (new_type=="size_t") {
    static const size_t narr=1;
    comm_option_s options_arr[narr]=
      {
        {0,"value","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_value),both}
      };
    if (narr!=type_comm_list["size_t"].size()) {
      O2SCL_ERR("Type comm list does not match for size_t",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
  } else if (new_type=="string") {
    static const size_t narr=1;
    comm_option_s options_arr[narr]=
      {
        {0,"value","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_value),both}
      };
    if (narr!=type_comm_list["string"].size()) {
      O2SCL_ERR("Type comm list does not match for string",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
  } else if (new_type=="table") {
    static const size_t narr=45;
    comm_option_s options_arr[narr]=
      {{0,"add-vec","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_add_vec),both},
       {0,"average-rows","",0,3,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_average_rows),both},
       {'a',"assign","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_assign),both},
       {0,"cat","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_cat),both},
       {0,"convert-unit","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_convert_unit),both},
       {0,"correl","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_correl),both},
       {0,"delete-col","",0,1,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_delete_col),both},
       {'d',"delete-rows","",0,1,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_delete_rows),both},
       {0,"delete-rows-tol","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_delete_rows_tol),both},
       {0,"deriv","",0,3,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_deriv),both},
       {0,"deriv2","",0,3,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_deriv2),both},
       {0,"value","",0,3,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_value),both},
       {0,"value-grid","",0,4,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_value_grid),both},
       {0,"find-row","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_find_row),both},
       {0,"fit","",0,7,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_fit),both},
       {0,"function","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_function),both},
       {0,"get-row","",0,1,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_get_row),both},
       {0,"get-unit","",0,1,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_get_unit),both},
       {'N',"index","",0,1,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_index),both},
       {0,"insert","",0,6,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_insert),both},
       {0,"insert-full","",0,4,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_insert_full),both},
       {0,"integ","",0,3,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_integ),both},
       {0,"interp","",0,3,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_interp),both},
       {0,"interp-table3d","",0,-1,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_interp_table3d),both},
       {'l',"list","",0,0,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_list),both},
       {0,"max","",0,1,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_max),both},
       {0,"min","",0,1,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_min),both},
       {0,"nlines","",0,0,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_nlines),both},
       {0,"refine","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_refine),both},
       {0,"rename","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_rename),both},
       {'s',"select","",-1,-1,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_select),both},
       {0,"select-rows","",0,1,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_select_rows),both},
       {0,"ser-hist-t3d","",0,8,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_ser_hist_t3d),both},
       {0,"set-data","",3,4,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_set_data),both},
       {0,"set-unit","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_set_unit),both},
       {'S',"sort","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_sort),both},
       {0,"stats","",0,1,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_stats),both},
       {0,"sum","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_sum),both},
       {0,"to-hist","",0,3,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_to_hist),both},
       {0,"to-hist-2d","",0,5,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_to_hist_2d),both},
       {0,"to-table3d","",0,4,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_to_table3d),both},
       {0,"to-gaussian","",0,4,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_to_gaussian),both},
       {0,"to-gmm","",-1,-1,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_to_gmm),both},
       {0,"to-pdma","",0,4,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_to_pdma),both},
       {0,"wstats","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_wstats),both}
      };
    if (narr!=type_comm_list["table"].size()) {
      std::cout << "1,2: " << narr << " "
                << type_comm_list["table"].size() << std::endl;
      O2SCL_ERR("Type comm list does not match for table",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);

  } else if (new_type=="table3d") {
    
    static const size_t narr=28;
    comm_option_s options_arr[narr]=
      {{0,"to-tensor-grid","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_to_tensor_grid),both},
       {0,"to-table","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_to_table),both},
       {0,"to-hist","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_to_hist),both},
       {0,"to-tg-fermi","",0,5,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_to_tg_fermi),both},
       {0,"get-grid","",0,0,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_get_grid),both},
       {'s',"select","",-1,-1,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_select),both},
       {0,"cat","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_cat),both},
       {0,"contours","",0,5,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_contours),both},
       {0,"deriv-x","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_deriv_x),both},
       {0,"deriv-y","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_deriv_y),both},
       {0,"value","",0,4,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_value),both},
       {0,"value-grid","",0,4,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_value_grid),both},
       {'f',"function","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_function),both},
       {0,"insert","",0,6,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_insert),both},
       {0,"interp","",0,3,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_interp),both},
       {'l',"list","",0,0,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_list),both},
       {0,"max","",0,1,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_max),both},
       {0,"min","",0,1,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_min),both},
       {0,"refine","",0,-1,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_refine),both},
       {0,"rename","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_rename),both},
       {0,"set-data","",0,4,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_set_data),both},
       {0,"slice","Construct a slice.",2,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_slice),both},
       {0,"slice-hist","",1,1,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_slice_hist),both},
       {0,"stats","",0,1,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_stats),both},
       {0,"sum","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_sum),both},
       {0,"to-hist-2d","",0,1,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_to_hist_2d),both},
       {0,"x-name","",0,1,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_x_name),both},
       {0,"y-name","",0,1,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_y_name),both}
      };
    if (narr!=type_comm_list["table3d"].size()) {
      O2SCL_ERR("Type comm list does not match for table3d",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="tensor") {
    
    static const size_t narr=13;
    comm_option_s options_arr[narr]=
      {
        {0,"sum","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_sum),both},
        {0,"deriv","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_deriv),both},
        {0,"stats","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_stats),both},
        {0,"diag","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_diag),both},
        {0,"value","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_value),both},
        {'f',"function","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_function),both},
        {'l',"list","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_list),both},
        {0,"max","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_max),both},
        {0,"min","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_min),both},
        {0,"rearrange","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_rearrange),both},
        {0,"to-table3d","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_to_table3d),both},
        {0,"to-table3d-sum","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_to_table3d_sum),both},
        {0,"to-tensor-grid","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_to_tensor_grid),both}
      };
    if (narr!=type_comm_list["tensor"].size()) {
      O2SCL_ERR("Type comm list does not match for tensor",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="tensor<int>") {
    
    static const size_t narr=5;
    comm_option_s options_arr[narr]=
      {
        {'l',"list","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_list),both},
        {0,"max","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_max),both},
        {0,"min","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_min),both},
        {0,"rearrange","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_rearrange),both},
        {0,"to-table3d","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_to_table3d),both}
      };
    if (narr!=type_comm_list["tensor<int>"].size()) {
      O2SCL_ERR("Type comm list does not match for tensor<int>",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="tensor<size_t>") {
    
    static const size_t narr=5;
    comm_option_s options_arr[narr]=
      {
        {0,"rearrange","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_rearrange),both},
        {'l',"list","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_list),both},
        {0,"to-table3d","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_to_table3d),both},
        {0,"max","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_max),both},
        {0,"min","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_min),both}
      };
    if (narr!=type_comm_list["tensor<size_t>"].size()) {
      O2SCL_ERR("Type comm list does not match for tensor<size_t>",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="tensor_grid") {
    
    static const size_t narr=18;
    comm_option_s options_arr[narr]=
      {
        {0,"sum","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_sum),both},
        {0,"binary","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_binary),both},
        {0,"stats","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_stats),both},
        {0,"value","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_value),both},
        {0,"deriv","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_deriv),both},
        {0,"value-grid","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_value_grid),both},
        {'f',"function","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_function),both},
        {0,"get-grid","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_get_grid),both},
        {0,"interp","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_interp),both},
        {'l',"list","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_list),both},
        {0,"max","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_max),both},
        {0,"min","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_min),both},
        {0,"rearrange","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_rearrange),both},
        {0,"slice","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_slice),both},
        {0,"set-grid","",0,2,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_set_grid),both},
        {0,"to-table3d","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_to_table3d),both},
        {0,"to-table","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_to_table),both},
        {0,"to-tensor","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_to_tensor),both}
      };
    if (narr!=type_comm_list["tensor_grid"].size()) {
      O2SCL_ERR("Type comm list does not match for tensor_grid",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);

  } else if (new_type=="prob_dens_mdim_amr") {
    
    static const size_t narr=2;
    comm_option_s options_arr[narr]=
      {
        {0,"to-table3d","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_to_table3d),both},
        {0,"sample","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_sample),both}
      };
    if (narr!=type_comm_list["prob_dens_mdim_amr"].size()) {
      O2SCL_ERR("Type comm list does not match for prob_dens_mdim_amr",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="prob_dens_mdim_gmm") {
    
    static const size_t narr=2;
    comm_option_s options_arr[narr]=
      {
        {0,"list","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_list),both},
        {0,"sample","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_sample),both}
      };
    if (narr!=type_comm_list["prob_dens_mdim_gmm"].size()) {
      O2SCL_ERR("Type comm list does not match for prob_dens_mdim_amr",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="prob_dens_mdim_gaussian") {
    
    static const size_t narr=1;
    comm_option_s options_arr[narr]=
      {
        {0,"sample","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_sample),both}
      };
    if (narr!=type_comm_list["prob_dens_mdim_gaussian"].size()) {
      O2SCL_ERR("Type comm list does not match for prob_dens_mdim_gaussian",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="hist") {

    static const size_t narr=3;
    comm_option_s options_arr[narr]=
      {
        {0,"function","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_function),both},
        {0,"to-table","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_to_table),both},
        {'l',"list","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_list),both},
      };
    if (narr!=type_comm_list["hist"].size()) {
      O2SCL_ERR("Type comm list does not match for hist",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="double[]") {
    
    static const size_t narr=11;
    comm_option_s options_arr[narr]=
      {
        {0,"deriv","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_deriv),both},
        {0,"value","",0,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_value),both},
        {0,"find","",0,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_find),both},
        {0,"interp","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_interp),both},
        {0,"max","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_max),both},
        {0,"min","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_min),both},
        {0,"resize","",0,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_resize),both},
        {0,"sort","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_sort),both},
        {0,"sum","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_sum),both},
        {0,"to-table","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_to_table),both},      
        {0,"function","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_function),both}      
      };
    if (narr!=type_comm_list["double[]"].size()) {
      O2SCL_ERR("Type comm list does not match for double[]",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="int[]") {

    static const size_t narr=11;
    comm_option_s options_arr[narr]=
      {
        {0,"value","",0,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_value),both},
        {0,"sort","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_sort),both},
        {0,"find","",0,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_find),both},
        {0,"sum","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_sum),both},
        {0,"max","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_max),both},
        {0,"min","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_min),both},
        {0,"deriv","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_deriv),both},
        {0,"interp","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_interp),both},
        {0,"resize","",0,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_resize),both},
        {0,"to-table","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_to_table),both},
        {0,"function","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_function),both}     
      };
    if (narr!=type_comm_list["int[]"].size()) {
      O2SCL_ERR("Type comm list does not match for int[]",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="string[]") {

    static const size_t narr=2;
    comm_option_s options_arr[narr]=
      {
        {0,"value","",0,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_value),both},
        {0,"resize","",0,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_resize),both}
      };
    if (narr!=type_comm_list["string[]"].size()) {
      O2SCL_ERR("Type comm list does not match for string[]",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="size_t[]") {

    static const size_t narr=11;
    comm_option_s options_arr[narr]=
      {
        {0,"value","",0,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_value),both},
        {0,"sort","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_sort),both},
        {0,"find","",0,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_find),both},
        {0,"sum","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_sum),both},
        {0,"max","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_max),both},
        {0,"min","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_min),both},
        {0,"deriv","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_deriv),both},
        {0,"interp","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_interp),both},
        {0,"resize","",0,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_resize),both},
        {0,"to-table","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_to_table),both},
        {0,"function","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_function),both}
      };
    if (narr!=type_comm_list["size_t[]"].size()) {
      O2SCL_ERR("Type comm list does not match for size_t[]",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="vector<contour_line>") {

    static const size_t narr=1;
    comm_option_s options_arr[narr]=
      {
        {'l',"list","",0,0,"","",
         new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_list),both}
      };
    if (narr!=type_comm_list["vector<contour_line>"].size()) {
      O2SCL_ERR("Type comm list does not match for vector<contour_line>",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="vec_vec_double") {

    static const size_t narr=1;
    comm_option_s options_arr[narr]=
      {
        {'l',"list","",0,0,"","",
         new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_list),both}
      };
    if (narr!=type_comm_list["vec_vec_double"].size()) {
      O2SCL_ERR("Type comm list does not match for vec_vec_double",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="vec_vec_string") {

    static const size_t narr=1;
    comm_option_s options_arr[narr]=
      {
        {'l',"list","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_list),both}
      };
    if (narr!=type_comm_list["vec_vec_string"].size()) {
      O2SCL_ERR("Type comm list does not match for vec_vec_string",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="hist_2d") {

    static const size_t narr=6;
    comm_option_s options_arr[narr]=
      {
        {0,"refine","",0,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_refine),both},
        {'l',"list","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_list),both},
        {0,"max","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_max),both},
        {0,"min","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_min),both},
        {0,"to-table3d","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_to_table3d),both},
        {0,"contours","",0,4,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_contours),both},
      };
    if (narr!=type_comm_list["hist_2d"].size()) {
      O2SCL_ERR("Type comm list does not match for hist_2d",
                o2scl::exc_esanity);
    }
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
  }

  return;
}

void acol_manager::command_del(std::string loc_type) {

  std::map<std::string,std::vector<std::string> >::iterator it;
  for(it=type_comm_list.begin();it!=type_comm_list.end();it++) {
    if (it->first==loc_type) {
      std::vector<std::string> &clist=it->second;
      for(size_t j=0;j<clist.size();j++) {
        cl->remove_comm_option(clist[j]);
      }
    }
  }
  
  return;
}

void acol_manager::clear_obj() {

  if (type=="table") {
    table_obj.clear();
  } else if (type=="table3d") {
    table3d_obj.clear();
  } else if (type=="tensor") {
    tensor_obj.clear();
  } else if (type=="tensor<int>") {
    tensor_int_obj.clear();
  } else if (type=="tensor<size_t>") {
    tensor_size_t_obj.clear();
  } else if (type=="hist") {
    hist_obj.clear();
  } else if (type=="hist_2d") {
    hist_2d_obj.clear();
  } else if (type=="vector<contour_line>") {
    cont_obj.clear();
  } else if (type!="string[]") {
    stringv_obj.clear();
  } else if (type!="int[]") {
    intv_obj.clear();
  } else if (type!="double[]") {
    doublev_obj.clear();
  } else if (type!="string") {
    string_obj.clear();
  } else if (type!="size_t[]") {
    size_tv_obj.clear();
  }
  
  type="";
  
  return;
}

int acol_manager::setup_options() {

  cl->cmd_name="acol";

  const int cl_param=cli::comm_option_cl_param;
  const int both=cli::comm_option_both;

  static const int narr=22;

  string type_list_str;
  for(size_t i=0;i<type_list.size()-1;i++) {
    type_list_str+=type_list[i]+", ";
  }
  type_list_str+="or "+type_list[type_list.size()-1]+'.';
  
  terminal ter;

  opts_new.resize(narr);
  opts_new[0]={0,"autocorr","",0,-1,"","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_autocorr),
    both};
  opts_new[1]={0,"calc","",0,2,"","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_calc),
    both};
  opts_new[2]={0,"clear","",0,0,"","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_clear),
    both};
  opts_new[3]={'c',"create","",0,-1,"","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_create),
    both};
  opts_new[4]={0,"docs","",0,1,"","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_docs),
    both};
  opts_new[5]={0,"wdocs","",0,2,"","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_wdocs),
    both};
  opts_new[6]={0,"download","",0,4,"","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_download),
    both};
  opts_new[7]={0,"filelist","",0,2,"","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_filelist),
    both};
  opts_new[8]={'g',"generic","",0,2,"","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_generic),
    both};
  opts_new[9]={0,"convert","",0,9,"","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_convert),
    both};
  opts_new[10]={0,"h5-copy","",-1,-1,"","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_h5_copy),
    both};
  opts_new[11]={0,"constant","",0,-1,"","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_constant),
    both};
  opts_new[12]={'q',"interactive","",0,0,"","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_interactive),
    cl_param};
  opts_new[13]={'i',"internal","",0,1,"","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_internal),
    both};
  opts_new[14]={0,"ninteg","",0,5,"","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_ninteg),
    both};
  opts_new[15]={'o',"output","",0,1,"","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_output),
    both};
  opts_new[16]={'P',"preview","",0,2,"","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_preview),
    both};
  opts_new[17]={'r',"read","",0,2,"","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_read),
    both};
  opts_new[18]={0,"slack","",0,6,"","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_slack),
    both};
  opts_new[19]={0,"type","",0,0,"","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_type),
    both};
  opts_new[20]={'v',"version","",0,0,"","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_version),
    both};
  opts_new[21]={0,"xml-to-o2","",0,0,"","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_xml_to_o2),
    both};

    /*
  // Options, sorted by long name. We allow 0 parameters in many of these
  // options so they can be requested from the user in interactive mode. 
  comm_option_s options_arr[narr]=
    {{0,"autocorr","",0,-1,"","",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_autocorr),
       both},
     {0,"calc","",0,2,"","",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_calc),
      both},
     {0,"clear","",0,0,"","",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_clear),
      both},
     {'c',"create","",0,-1,"","",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_create),
      both},
     {0,"docs","",0,1,"","",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_docs),
      both},
     {0,"wdocs","",0,2,"","",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_wdocs),
      both},
     {0,"download","",0,4,"","",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_download),
      both},
     {0,"filelist","",0,1,"","",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_filelist),
      both},
     {'g',"generic","",0,2,"","",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_generic),
      both},
     {0,"convert","",0,9,"","",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_convert),
      both},
     {0,"h5-copy","",-1,-1,"","",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_h5_copy),
      both},
     {0,"constant","",0,-1,"","",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_constant),
      both},
     {'q',"interactive","",0,0,"","",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_interactive),
      cl_param},
     {'i',"internal","",0,1,"","",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_internal),
      both},
     {0,"ninteg","",0,5,"","",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_ninteg),
      both},
     {'o',"output","",0,1,"","",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_output),
      both},
     {'P',"preview","",0,2,"","",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_preview),
      both},
     {'r',"read","",0,2,"","",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_read),
      both},
     {0,"slack","",0,6,"","",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_slack),
      both},
     {0,"type","",0,0,"","",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_type),
      both},
     {'v',"version","",0,0,"","",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_version),
      both},
     {0,"xml-to-o2","",0,0,"","",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_xml_to_o2),
      both}
    };
    */

  cl->remove_comm_option("xml-to-o2");

  if (true) {
    std::string doc_fn=o2scl::o2scl_settings.get_data_dir()+"/acol_docs.o2";
    if (file_exists(doc_fn)) {
      hdf_file hf;
      hf.open(doc_fn);
      hf.gets_vec_vec_copy("cmd_doc_strings",cmd_doc_strings);
      hf.gets_vec_vec_copy("param_doc_strings",param_doc_strings);
      hf.gets_vec_vec_copy("help_doc_strings",help_doc_strings);
      hf.close();
    } else {
      cout << "Couldn't find file " << doc_fn << endl;
    }
  }
  
  update_o2_docs(narr,&opts_new[0]);
  
  cl->set_comm_option_vec(narr,opts_new);
  
  return 0;
}

int acol_manager::setup_cli() {

  //---------------------------------------------------------------------
  // Get HOME directory and command history
  
  char *hd=getenv("HOME");
  std::string histfile;
  if (hd) {
    histfile=hd;
    histfile+="/.acol_hist";
  }
    
  //--------------------------------------------------------------------
  // Specify command-line option object
  
#ifdef O2SCL_READLINE
  cl=new cli_readline(histfile);
#else
  cl=new cli;
#endif

  //---------------------------------------------------------------------
  // Set colors
  
  char *ac=getenv("ACOL_COLORS");
  if (ac) {
    color_spec=ac;
    
    if (color_spec=="0") {
      cl->set_colors("c:,d:,e:,h:,p:,t:,u:");
    } else if (color_spec=="default") {
      cl->set_colors("c:10006,d:0,e:10015,h:10002,p:10001,t:10005,u:1000");
    } else {
      cl->set_colors(color_spec);
    }
    
    command_color=cl->command_color;
    type_color=cl->type_color;
    param_color=cl->param_color;
    help_color=cl->help_color;
    exec_color=cl->exec_color;
    url_color=cl->url_color;
  }
    
  return 0;
}

int acol_manager::setup_help() {

  cl->cmd_name="acol";

  terminal ter;
  //cl->desc=((string)"acol: A data viewing and processing ")+
  //"program for "+ter.bold()+"O₂scl"+ter.default_fgbg()+".\n";
  cl->desc=((string)"acol: A data viewing and processing ")+
    "program for O₂scl.\n";
  
  return 0;
}

int acol_manager::setup_parameters() {
  
  p_obj_name.str=&obj_name;
  p_def_args.str=&def_args;
  p_color_spec.str=&color_spec;
  p_verbose.i=&verbose;
  p_compress.i=&compress;
  p_precision.i=&precision;
  p_ncols.i=&ncols;
  p_interp_type.i=&interp_type;
  p_scientific.b=&scientific;
  p_pretty.b=&pretty;
  p_names_out.b=&names_out;
  p_use_regex.b=&use_regex;
  
  p_obj_name.help="The current object name.";
  p_color_spec.help="The color specification for terminal output.";
  p_def_args.help=((std::string)"The default arguments from the ")+
    "environment varable "+env_var_name+".";
  p_precision.help="The numerical precision.";
  p_verbose.help="Control the amount of output.";
  p_compress.help=((std::string)"If true, enable compression ")+
    "(defaults to true if compression was enabled in O2scl).";
  p_ncols.help="The number of columns for screen output.";
  p_interp_type.help=((std::string)"The interpolation type ")+
    "(1=linear, 2=cubic spline, 3=periodic cubic spline, 4=Akima, "+
    "5=periodic Akima, 6=monotonic, 7=Steffen's monotonic).";
  p_names_out.help="If true, output column names at top.";
  p_use_regex.help="If true, use regex.";
  p_pretty.help="If true, make the output more readable.";
  p_scientific.help="If true, output in scientific mode.";
  
  cl->par_list.insert(make_pair("obj_name",&p_obj_name));
  cl->par_list.insert(make_pair("def_args",&p_def_args));
  cl->par_list.insert(make_pair("color_spec",&p_color_spec));
  cl->par_list.insert(make_pair("precision",&p_precision));
  cl->par_list.insert(make_pair("verbose",&p_verbose));
  cl->par_list.insert(make_pair("compress",&p_compress));
  cl->par_list.insert(make_pair("ncols",&p_ncols));
  cl->par_list.insert(make_pair("interp_type",&p_interp_type));
  cl->par_list.insert(make_pair("names_out",&p_names_out));
  cl->par_list.insert(make_pair("use_regex",&p_use_regex));
  cl->par_list.insert(make_pair("pretty",&p_pretty));
  cl->par_list.insert(make_pair("scientific",&p_scientific));

  if (true) {
    for(cli::par_t it=cl->par_list.begin();it!=cl->par_list.end();it++) {
      bool found=false;
      for(size_t k=0;k<param_doc_strings.size() && found==false;k++) {
        if (param_doc_strings[k][0]==it->first) {
          it->second->help=param_doc_strings[k][1];
          for(size_t kk=2;kk<param_doc_strings[k].size();kk++) {
            it->second->help+="\n\n"+param_doc_strings[k][kk];
          }
          found=true;
        }
      }
      if (found) {
        if (verbose>2) {
          cout << "Function acol_manager::setup_parameters() found "
               << "documentation for parameter "
               << it->first << " ." << endl;
        }
      } else {
        cout << "Function acol_manager::setup_parameters() could "
             << "not find documentation for parameter "
             << it->first << " ." << endl;
      }
    }
  }
  
  return 0;
}

int acol_manager::run(int argc, char *argv[], bool full_process) {

  //--------------------------------------------------------------------
  // Default to scientific mode

  cout.setf(ios::scientific);

  // 
  if (verbose>2) {
    cout << "Setup cli class." << endl;
  }
  setup_cli();

  // 
  if (verbose>2) {
    cout << "Setup options." << endl;
  }
  setup_options();

  // 
  if (verbose>2) {
    cout << "Setup help." << endl;
  }
  setup_help();

  //-------------------------------------------------------------------

  cl->set_function(cset);

  static const int both=cli::comm_option_both;
  cl->remove_comm_option("help");
  cl->remove_comm_option("commands");
  static const size_t narr2=2;
  comm_option_s options_arr2[narr2]=
    {
      {'h',"help","Show help information.",0,2,
       "[type command] or [command] or [topic]",
       ((std::string)"Show generic help information, or, if an ")+
       "argument is given "+
       "give the documentation for the specified command or topic. "+
       "If two arguments are given, show the help for a type-specific "+
       "command. "+
       "Note that required arguments are typically given inside "+
       "angled brackes <> while optional arguments are given "+
       "inside square brackets [].",
       new o2scl::comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_help),both},
      {0,"commands",
       "List available commands for current or specified type.",
       0,1,"[type or \"all\"]",((string)"Output the commands available for ")+
       "the current type, or, if the optional type argument is given "+
       "then output the commands available for that type. \"commands all\""+
       " outputs all commands for any type.",
       new o2scl::comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_commands),both}
    };
  cl->set_comm_option_vec(narr2,options_arr2);

  //-------------------------------------------------------------------
  // Try to get screen width

  // AWS 3/29/22: we no longer try to do this here but in the
  // individual functions which have complex screen output
  
  // Use curses
  // AWS: this prints out a bunch of crazy characters on Ubuntu
  //get_screen_size_curses(nrow,ncol);
  //#else
  
  // If not, attempt to obtain the result from the environment
  /*
    char *ncstring=getenv("COLUMNS");
    if (ncstring) {
    int nc2;
    int sret=o2scl::stoi_nothrow(ncstring,nc2);
    if (sret==0 && nc2>0) {
    ncol=nc2;
    } else {
    cerr << "Failed to interpret COLUMNS value " << ncstring
    << " as a positive number of columns." << endl;
    }
    }
  */
  
  //#endif
  
  //-------------------------------------------------------------------
  // Setup parameters modified by 'set' and 'get'

  if (verbose>2) {
    cout << "Setup parameters: " << endl;
  }
  setup_parameters();

  //-------------------------------------------------------------------
  // Process default options and call

  std::vector<cmd_line_arg> ca;
  
  char *dc=getenv(env_var_name.c_str());
  if (dc) {
    if (verbose>2) {
      cout << "Process default options: " << dc << endl;
    }
    def_args=dc;
    if (verbose>2) {
      cl->process_args_str(def_args,ca,1,true);
    } else {
      cl->process_args_str(def_args,ca,0,true);
    }
  }
  
  if (full_process) {
    
    //----------------------------------------------------------------
    // Process command-line options
    
    // Note that it's ok that this appears early in the code because it
    // just processes the arguments, it doesn't do any execution based
    // on those arguments until later.
    
    if (verbose>2) {
      cout << "Process command-line options" << endl;
      cl->process_args_c(argc,argv,ca,1,true);
    } else {
      cl->process_args_c(argc,argv,ca,0,true);
    }
    if (argc<2) {
      post_interactive=true;
    }
    
    //------------------------------------------------------------------
    // Post interactive
    
    int ret2=0;
    
    if (verbose>2) {
      cout << "Post_interactive: " << post_interactive << endl;
    }
    
    if (post_interactive) {
      if (verbose>2) {
        cout << "Run interactive mode." << endl;
      }
      ret2=cl->run_interactive();
    }
    
    //--------------------------------------------------------------------
    // Notify user if error occurred
    
    if (ret2!=0) {
      cout << "An error occured." << endl;
    }
    
    delete cl;
  
  }

  return 0;

}

int acol_manager::get_input_one(vector<string> &sv, string directions,
                                string &in, string comm_name,
                                bool itive_com) {

  // If there are enough arguments, then just fill 'in' with the
  // correct values from 'sv'
  if (sv.size()>1) {
    in=sv[1];
    return 0;
  }
  
  if (itive_com) {

    string temp=directions+" (or blank to stop): ";
    in=cl->cli_gets(temp.c_str());
    if (in.length()==0 || o2scl::count_words(in)==0) {
      if (verbose>0) {
        cout << "Command '" << comm_name << "' cancelled." << endl;
      }
      return exc_efailed;
    }
    
  } else {
    cerr << "Not enough arguments to '" << comm_name << "'." << endl;
    return exc_efailed;
  }

  return success;
}

int acol_manager::get_input(vector<string> &sv, vector<string> &directions,
                            vector<string> &in, string comm_name,
                            bool itive_com) {

  size_t ni=directions.size();

  // If there are enough arguments, then just fill the vector
  // 'in' with the correct values from 'sv'
  if (sv.size()>ni) {
    for(size_t i=0;i<sv.size()-1;i++) {
      in.push_back(sv[i+1]);
    }
    return 0;
  }

  // Otherwise, if we're in interactive mode
  if (itive_com) {
    
    // Prompt the user for the correct arguments
    for(size_t i=0;i<ni;i++) {
      string temp=directions[i]+" (or blank to stop): ";
      in.push_back(cl->cli_gets(temp.c_str()));
      // If the user just pressed 'enter', then cancel
      if (in[i].length()==0 || count_words(in[i])==0) {
        cout << "Command '" << comm_name << "' cancelled." << endl;
        return exc_efailed;
      }
    }
    
  } else {
    // We don't have enough arguments and we're not in interactive
    // mode, so we fail
    cerr << "Not enough arguments to '" << comm_name << "'." << endl;
    return exc_efailed;
  }

  return success;
}

int acol_manager::validate_interp_type() {

  if (interp_type<1) {
    cout << "Interpolation type invalid, less than 1. Setting to 1 (linear)."
         << endl;
    interp_type=1;
    return 1;
  } else if (interp_type>10) {
    cout << "Interpolation type invalid, greater than 10. "
         << "Setting to 1 (linear)."
         << endl;
    interp_type=1;
    return 1;
  }

  return 0;
}
