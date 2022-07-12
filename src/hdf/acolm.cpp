/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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

  ffl.html_mode();

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
  type_list.push_back("vec_vec_string");
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
    vector<std::string> itmp={"ac-len","add-vec","average-rows",
      "assign","cat","convert-unit",
      "correl","delete-col",
      "delete-rows","delete-rows-tol",
      "deriv","deriv2","entry-grid",
      "find-row","fit","function",
      "get-row","get-unit","entry","index",
      "insert","insert-full","integ","interp",
      "list","max","min","nlines","refine","rename",
      "select","select-rows",
      "set-data","set-unit","sort","stats","sum",
      "to-hist","to-hist-2d","to-table3d","wstats",
      "ser-hist-t3d",
    };
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("table",itmp));
  }
  {
    vector<std::string> itmp={"cat","contours","deriv-x","deriv-y",
      "function","entry","entry-grid","get-grid",
      "insert","interp","stats","select",
      "list","max","min","rename","set-data",
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
      "entry","function","sum","stats","deriv"};
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("tensor",itmp));
  }
  {
    vector<std::string> itmp={"to-table3d"};
    type_comm_list.insert(std::make_pair("prob_dens_mdim_amr",itmp));
  }
  {
    vector<std::string> itmp={"list","to-table3d","slice","to-table",
      "set-grid","max","min","rearrange",
      "get-grid","interp","entry","to-tensor",
      "entry-grid","function","sum","stats",
      "binary","deriv"};
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("tensor_grid",itmp));
  }
  {
    vector<std::string> itmp={"max","min","contours","list","to-table3d"};
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("hist_2d",itmp));
  }
  {
    vector<std::string> itmp={"to-table","function"};
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("hist",itmp));
  }
  {
    vector<std::string> itmp={"deriv","interp","max","min","sort",
      "to-table","function","sum"};
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("double[]",itmp));
    type_comm_list.insert(std::make_pair("int[]",itmp));
    type_comm_list.insert(std::make_pair("size_t[]",itmp));
  }

  // Ensure the RNGs for different types are somewhat uncorrelated
  rng.clock_seed();
  rng_ld.set_seed(rng.get_seed()*2);
}

void acol_manager::update_o2_docs(size_t narr,
                                  o2scl::comm_option_s *options_arr,
                                  std::string new_type) {

  int loc_verbose=0;

  // Loop through all options, attempting to assign docs for all
  for(size_t j=0;j<narr;j++) {
    
    bool found=false;

    // Loop over all command doc strings, attempting to find a match
    for(size_t k=0;k<cmd_doc_strings.size() && found==false;k++) {

      // We found a match for option j and list of doc strings k
      if (cmd_doc_strings[k][0]==options_arr[j].lng) {

        //if (options_arr[j].lng==(string)"stats") loc_verbose=3;
        //else loc_verbose=0;
        
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
            
            if (cmd_doc_strings[k][2].substr(0,19)==
                ((string)"For objects of type")) {
              generic_docs=false;
            } else if (cmd_doc_strings[k][2].substr(0,30)==
                       ((string)"If there is no current object:")) {
              generic_docs=false;
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
            
            if (new_type.length()==0) {
              if (generic_docs==true) {
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
                cout << "Current type is '" << new_type << "'" << endl;
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

                if (cmd_doc_strings[k][kk].substr(0,s.length())==s ||
                    cmd_doc_strings[k][kk].substr(0,s2.length())==s2) {
                  
                  if (loc_verbose>1) {
                    cout << "Found type-specific docs [" << new_type
                         << "]." << endl;
                  }
                  loop1_done=true;
                  bool loop2_done=false;
                  
                  for(size_t kl=kk+1;kl<cmd_doc_strings[k].size() &&
                        loop2_done==false;kl++) {
                    
                    if (cmd_doc_strings[k][kl].substr(0,19)==
                        ((string)"For objects of type")) {
                      loop2_done=true;
                    } else if (cmd_doc_strings[k][kl].substr(0,30)==
                               ((string)"If there is no current object:")) {
                      loop2_done=true;
                    } else if (options_arr[j].desc=="") {
                      if (loc_verbose>2) {
                        cout << "Adding desc: "
                             << cmd_doc_strings[k][kl] << endl;
                      }
                      options_arr[j].desc=cmd_doc_strings[k][kl];
                    } else if (options_arr[j].parm_desc=="") {
                      if (loc_verbose>2) {
                        cout << "Adding parm_desc: "
                             << cmd_doc_strings[k][kl] << endl;
                      }
                      options_arr[j].parm_desc=cmd_doc_strings[k][kl];
                    } else if (options_arr[j].help=="") {
                      if (loc_verbose>2) {
                        cout << "Adding first doc string: \n  "
                             << cmd_doc_strings[k][kl] << endl;
                      }
                      options_arr[j].help=cmd_doc_strings[k][kl];
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
        cout << "Function acol_manager::update_o2_docs() "
             << "found documentation for command "
             << options_arr[j].lng << " ." << endl;
      }
    } else {
      cout << "Function acol_manager::update_o2_docs() could not "
           << "find documentation for command "
           << options_arr[j].lng << " ." << endl;
    }
    if (loc_verbose>1) {
      cout << endl;
    }

    // End of 'j' loop over list of options
  }
    
  return;
}

void acol_manager::command_add(std::string new_type) {

  const int both=cli::comm_option_both;
  
  terminal ter;
  
  if (new_type=="int") {
    static const size_t narr=1;
    comm_option_s options_arr[narr]=
      {
        {0,"value","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_value),both}
      };
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
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
  } else if (new_type=="table") {
    static const size_t narr=42;
    comm_option_s options_arr[narr]=
      {{0,"ac-len","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_ac_len),both},
       {0,"add-vec","",0,2,"","",
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
       {0,"entry","",0,3,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_entry),both},
       {0,"entry-grid","",0,4,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_entry_grid),both},
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
       {0,"wstats","",0,2,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_wstats),both}
      };
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);

  } else if (new_type=="table3d") {
    
    static const size_t narr=26;
    comm_option_s options_arr[narr]=
      {{0,"to-tensor-grid","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_to_tensor_grid),both},
       {0,"to-table","",0,2,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_to_table),both},
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
       {0,"entry","",0,4,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_entry),both},
       {0,"entry-grid","",0,4,"","",
        new comm_option_mfptr<acol_manager>
        (this,&acol_manager::comm_entry_grid),both},
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
        {0,"entry","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_entry),both},
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
        {0,"entry","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_entry),both},
        {0,"deriv","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_deriv),both},
        {0,"entry-grid","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_entry_grid),both},
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
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);

  } else if (new_type=="prob_dens_mdim_amr") {
    
    static const size_t narr=1;
    comm_option_s options_arr[narr]=
      {
        {0,"to-table3d","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_to_table3d),both}
      };
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="hist") {

    static const size_t narr=2;
    comm_option_s options_arr[narr]=
      {
        {0,"function","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_function),both},
        {0,"to-table","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_to_table),both}
      };
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="double[]") {
    
    static const size_t narr=8;
    comm_option_s options_arr[narr]=
      {
        {0,"deriv","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_deriv),both},
        {0,"interp","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_interp),both},
        {0,"max","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_max),both},
        {0,"min","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_min),both},
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
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="int[]") {

    static const size_t narr=8;
    comm_option_s options_arr[narr]=
      {
        {0,"sort","",0,0,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_sort),both},
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
        {0,"to-table","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_to_table),both},
        {0,"function","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_function),both}     
      };
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="size_t[]") {

    static const size_t narr=8;
    comm_option_s options_arr[narr]=
      {
        {0,"sort","",0,0,"","",
         new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_sort),
         both},
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
        {0,"to-table","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_to_table),both},
        {0,"function","",0,1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_function),both}
      };
    update_o2_docs(narr,&options_arr[0],new_type);
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="vector<contour_line>") {

  } else if (new_type=="hist_2d") {

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
        {0,"to-table3d","",-1,-1,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_to_table3d),both},
        {0,"contours","",0,4,"","",
         new comm_option_mfptr<acol_manager>
         (this,&acol_manager::comm_contours),both},
      };
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

  // Options, sorted by long name. We allow 0 parameters in many of these
  // options so they can be requested from the user in interactive mode. 
  comm_option_s options_arr[narr]=
    {{0,"autocorr","",0,-1,"","",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_autocorr),
       both},
     {0,"calc","",0,1,"","",
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
     {0,"ninteg","",0,4,"","",
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

  cl->remove_comm_option("xml-to-o2");

  if (true) {
    std::string doc_fn=o2scl::o2scl_settings.get_data_dir()+"/acol_docs.o2";
    if (file_exists(doc_fn)) {
      hdf_file hf;
      hf.open(doc_fn);
      hf.gets_vec_vec("cmd_doc_strings",cmd_doc_strings);
      hf.gets_vec_vec("param_doc_strings",param_doc_strings);
      hf.gets_vec_vec("help_doc_strings",help_doc_strings);
      hf.close();
    } else {
      cout << "Couldn't find file " << doc_fn << endl;
    }
  }
  
  /*
    {0,"find-x","Find an entry in the x-grid (3d only)",1,1,"<value>","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_find_x),
    comm_option::both},
    {0,"find-y","Find an entry in the y-grid (3d only)",1,1,"<value>","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_find_y),
    comm_option::both},
    {0,"find-xy",
    "Find the closest grid point to a pair of values (3d only)",
    1,1,"<x value> <y value>","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_find_y),
    comm_option::both},
  */

  update_o2_docs(narr,&options_arr[0]);
  
  cl->set_comm_option_vec(narr,options_arr);
  
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

  return 0;
}

int acol_manager::setup_help() {

  cl->cmd_name="acol";

  terminal ter;
  cl->desc=((string)"acol: A data viewing and processing ")+
    "program for "+ter.bold()+"O₂scl"+ter.default_fg()+".\n";

#ifdef O2SCL_NEVER_DEFINED

  // This code is now moved to acolm_gi.cpp
  
  ostringstream oss;
  oss << ((char)27) << '(' << '0';
  for(size_t i=0;i<78;i++) oss << 'q';
  oss << ((char)27) << '(' << 'B';
  string line=oss.str();
  
  string stemp;
  string dsc=line+"\nNotes:\n\n";
  vector<std::string> sv;
  
  stemp="1. Help for general commands may be obtained with '-help ";
  stemp+="<command>'. Help for type-specific commands can be obtained ";
  stemp+="by '-help <type> <command>'. A list of commands for each type ";
  stemp+="can be obtained with '-commands <type>', or for all commands ";
  stemp+=" use '-commands all'. Required arguments ";
  stemp+="are surrounded by ";
  stemp+="<>'s and optional arguments are surrounded by []'s.\n";
  rewrap_ignore_vt100(stemp,sv,76);
  dsc+=sv[0]+"\n";
  for(size_t j=1;j<sv.size();j++) {
    dsc+="   "+sv[j]+"\n";
  }
  
  stemp="2. Options may also be specified in the environment variable ";
  stemp+=env_var_name+".\n";
  rewrap_ignore_vt100(stemp,sv,76);
  dsc+=sv[0]+"\n";
  for(size_t j=1;j<sv.size();j++) {
    dsc+="   "+sv[j]+"\n";
  }

  stemp="3. Long options may be preceeded by two dashes.\n";
  rewrap_ignore_vt100(stemp,sv,76);
  dsc+=sv[0]+"\n";
  for(size_t j=1;j<sv.size();j++) {
    dsc+="   "+sv[j]+"\n";
  }

  stemp="4. In order to avoid confusion between arguments and functions, ";
  stemp+="use parenthesis and quotes, i.e. \"(-x*2)\" instead of -x*2.\n";
  rewrap_ignore_vt100(stemp,sv,76);
  dsc+=sv[0]+"\n";
  for(size_t j=1;j<sv.size();j++) {
    dsc+="   "+sv[j]+"\n";
  }

  stemp="5. Also, do not use a unary minus next to a binary operator, ";
  stemp+="i.e. use \"a>(-1)\" instead of \"a>-1\".\n\n";
  rewrap_ignore_vt100(stemp,sv,76);
  dsc+=sv[0]+"\n";
  for(size_t j=1;j<sv.size();j++) {
    dsc+="   "+sv[j]+"\n";
  }

  dsc+=line+"\n";
  
  dsc+="List of additional help topics (e.g. \"acol -help <topic>\"): ";
  dsc+="functions, mult-vector-spec, types, value-spec, and vector-spec.\n\n";
  
#ifndef O2SCL_UBUNTU_PKG
  dsc+=((string)"Compiled at ")+((string)__TIME__)+" on "+
    ((string)__DATE__)+" for O₂scl, version "+
    ((string)VERSION)+".\n";
#else
  dsc+=((string)"Compiled for O₂scl, version "+
    ((string)VERSION)+".\n";
#endif
  
  cl->addl_help_cmd=dsc;
  cl->addl_help_cli=dsc;

#endif
  
  return 0;
}

int acol_manager::setup_parameters() {
  
  p_obj_name.str=&obj_name;
  p_def_args.str=&def_args;
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

