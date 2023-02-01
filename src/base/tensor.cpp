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
#include <o2scl/tensor.h>

using namespace std;
using namespace o2scl;

ix_index::ix_index(size_t ix) {
  this->type=index_spec::index;
  this->ix1=ix;
  this->ix2=0;
  this->ix3=0;
  this->val1=0.0;
  this->val2=0.0;
  this->val3=0.0;
}

int o2scl::strings_to_indexes(std::vector<std::string> sv2,
                               std::vector<o2scl::index_spec> &vis,
                               int verbose, bool err_on_fail) {
    
  if (verbose>1) {
    std::cout << "Function strings_to_indexes(): " << std::endl;
    for(size_t j=0;j<sv2.size();j++) {
      std::cout << "  " << j << " " << sv2[j] << std::endl;
    }
  }
      
  for(size_t j=0;j<sv2.size();j++) {
    
    std::vector<std::string> args;
    
    if (sv2[j].find("index(")==0 && sv2[j][sv2[j].size()-1]==')') {
      
      std::string spec=sv2[j].substr(6,sv2[j].length()-7);
      split_string_delim(spec,args,',');
      if (verbose>1) {
        std::cout << "  rearrange, index  : ";
        vector_out(std::cout,args,true);
      }
      vis.push_back(ix_index(o2scl::stoszt(args[0])));
      
    } else if (sv2[j].find("fixed(")==0 && sv2[j][sv2[j].size()-1]==')') {
      
      std::string spec=sv2[j].substr(6,sv2[j].length()-7);
      split_string_delim(spec,args,',');
      if (verbose>1) {
        std::cout << "  rearrange, fixed  : ";
        vector_out(std::cout,args,true);
      }
      if (args.size()<2) {
        if (err_on_fail) {
          O2SCL_ERR("Not enough arguments in fixed().",o2scl::exc_einval);
        }
        return 2;
      }
      vis.push_back(ix_fixed(o2scl::stoszt(args[0]),
                             o2scl::stoszt(args[1])));
      
    } else if (sv2[j].find("sum(")==0 &&
               sv2[j][sv2[j].size()-1]==')') {
      
      std::string spec=sv2[j].substr(4,sv2[j].length()-5);
      split_string_delim(spec,args,',');
      if (verbose>1) {
        std::cout << "  rearrange, sum    : ";
        vector_out(std::cout,args,true);
      }
      
      vis.push_back(ix_sum(o2scl::stoszt(args[0])));
      
    } else if (sv2[j].find("trace(")==0 &&
               sv2[j][sv2[j].size()-1]==')') {
      
      std::string spec=sv2[j].substr(6,sv2[j].length()-7);
      split_string_delim(spec,args,',');
      if (verbose>1) {
        std::cout << "  rearrange, trace  : ";
        vector_out(std::cout,args,true);
      }
      if (args.size()<2) {
        if (err_on_fail) {
          O2SCL_ERR("Not enough arguments in trace().",o2scl::exc_einval);
        }
        return 3;
      }
      vis.push_back(ix_trace(o2scl::stoszt(args[0]),
                             o2scl::stoszt(args[1])));
      
    } else if (sv2[j].find("reverse(")==0 &&
               sv2[j][sv2[j].size()-1]==')') {
      
      std::string spec=sv2[j].substr(8,sv2[j].length()-9);
      split_string_delim(spec,args,',');
      if (verbose>1) {
        std::cout << "  rearrange, reverse: ";
        vector_out(std::cout,args,true);	
      }
      vis.push_back(ix_reverse(o2scl::stoszt(args[0])));
      
    } else if (sv2[j].find("range(")==0 &&
               sv2[j][sv2[j].size()-1]==')') {
      
      std::string spec=sv2[j].substr(6,sv2[j].length()-7);
      split_string_delim(spec,args,',');
      if (verbose>1) {
        std::cout << "  rearrange, range  : ";
        vector_out(std::cout,args,true);
      }
      if (args.size()<3) {
        if (err_on_fail) {
          O2SCL_ERR("Not enough arguments in range().",
                    o2scl::exc_einval);
        }
        return 4;
      }
      vis.push_back(ix_range(o2scl::stoszt(args[0]),
                             o2scl::stoszt(args[1]),
                             o2scl::stoszt(args[2])));
      
    } else if (sv2[j].find("interp(")==0 &&
               sv2[j][sv2[j].size()-1]==')') {
      
      std::string spec=sv2[j].substr(7,sv2[j].length()-8);
      split_string_delim(spec,args,',');
      if (verbose>1) {
        std::cout << "  rearrange, interp : ";
        vector_out(std::cout,args,true);
      }
      if (args.size()<2) {
        if (err_on_fail) {
          O2SCL_ERR("Not enough arguments in interp().",
                    o2scl::exc_einval);
        }
        return 5;
      }
      vis.push_back(ix_interp(o2scl::stoszt(args[0]),
                              o2scl::function_to_double(args[1])));
      
    } else if (sv2[j].find("grid(")==0 && sv2[j][sv2[j].size()-1]==')') {
      
      std::string spec=sv2[j].substr(5,sv2[j].length()-6);
      split_string_delim(spec,args,',');
      if (args.size()<4) {
        if (err_on_fail) {
          O2SCL_ERR("Not enough arguments in grid().",
                    o2scl::exc_einval);
        }
        return 6;
      }
      if (args.size()>4 && o2scl::stob(args[4])==true) {
        if (verbose>1) {
          std::cout << "  rearrange, grid   : ";
          vector_out(std::cout,args,true);
        }
        vis.push_back(ix_grid(o2scl::stoszt(args[0]),
                              o2scl::function_to_double(args[1]),
                              o2scl::function_to_double(args[2]),
                              o2scl::stoszt(args[3]),true));
      } else {
        if (verbose>1) {
          std::cout << "  rearrange, gridw  : ";
          vector_out(std::cout,args,true);
        }
        vis.push_back(ix_grid(o2scl::stoszt(args[0]),
                              o2scl::function_to_double(args[1]),
                              o2scl::function_to_double(args[2]),
                              o2scl::stoszt(args[3]),false));
      }
      
    } else if (sv2[j].find("gridw(")==0 &&
               sv2[j][sv2[j].size()-1]==')') {
      
      std::string spec=sv2[j].substr(6,sv2[j].length()-7);
      split_string_delim(spec,args,',');
      if (args.size()<4) {
        if (err_on_fail) {
          O2SCL_ERR("Not enough arguments in gridw().",
                    o2scl::exc_einval);
        }
        return 7;
      }
      if (args.size()>4 && o2scl::stob(args[4])==true) {
        if (verbose>1) {
          std::cout << "  rearrange, gridw, index, begin, "
                    << "end, bin_width, log: ";
          vector_out(std::cout,args,true);
        }
        vis.push_back(ix_gridw(o2scl::stoszt(args[0]),
                               o2scl::function_to_double(args[1]),
                               o2scl::function_to_double(args[2]),
                               o2scl::function_to_double(args[3]),true));
      } else {
        if (verbose>1) {
          std::cout << "  rearrange, gridw, index, "
                    << "begin, end, bin_width [no log]: ";
          vector_out(std::cout,args,true);
        }
        vis.push_back(ix_gridw(o2scl::stoszt(args[0]),
                               o2scl::function_to_double(args[1]),
                               o2scl::function_to_double(args[2]),
                               o2scl::function_to_double(args[3]),false));
      }
      
    } else {
      if (err_on_fail) {
        O2SCL_ERR2("Function strings_to_indexes() ",
                   "did not understand index specification.",
                   o2scl::exc_einval);
      }
      return 1;
    }
  }
      
  return 0;
}
  
void o2scl::index_spec_preprocess(std::string str,
                                  std::vector<std::string> &sv,
                                  int verbose) {
  
  int paren_count=0;
  std::string entry;
  if (verbose>1) {
    std::cout << "Function index_spec_preprocess(), before: "
              << str << std::endl;
  }
  for (size_t i=0;i<str.length();i++) {
    if (str[i]=='(') {
      entry+=str[i];
      paren_count++;
    } else if (str[i]==')') {
      entry+=str[i];
      paren_count--;
      if (paren_count==0) {
        sv.push_back(entry);
        entry.clear();
        // Skip the character between index specs
        i++;
      }
    } else {
      entry+=str[i];
    }
  }
  if (verbose>1) {
    std::cout << "Function index_spec_preprocess(), after: ";
    o2scl::vector_out(std::cout,sv,true);
  }
  return;
}
