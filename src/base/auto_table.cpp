/*
  -------------------------------------------------------------------
  
  Copyright (C) 2020, Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/auto_table.h>
#include <o2scl/misc.h>
#include <o2scl/string_conv.h>
#include <o2scl/columnify.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_auto_table;

auto_table::auto_table() {
  inside_table=false;
  row_max=100;
}      

void auto_table::done() {
  for(size_t j=0;j<lines.size();j++) {
    cout << lines[j] << endl;
  }
  return;
}

void auto_table::add_string(std::string s, bool include_endl) {

  if (s.length()==0 && include_endl==false) {
    return;
  }
  
  if (s=="\n") {
    if (include_endl==false) {
      s="";
      include_endl=true;
    } else {
      add_string("",true);
      add_string("",true);
      return;
    }
  }

  if (s.find('\n')!=std::string::npos) {
    //cout << "Iere: " << s.length() << " " << s[0] << endl;
    // If there is a carriage return character, then
    // use recursion to add each line in turn
    std::vector<std::string> vs1;
    split_string_delim(s,vs1,'\n');
    // If 's' ends in a \n, then add another empty string to the end.
    // This ensures that include_endl is set to true.
    if (s[s.length()-1]=='\n') {
      //cout << "Jere." << endl;
      vs1.push_back("");
    }
    for(size_t j=0;j<vs1.size();j++) {
      //cout << "split: " << j << ". <" << vs1[j] << ">" << std::endl;
      if (j!=vs1.size()-1) {
	add_string(vs1[j],true);
      } else {
	if (vs1[j].length()>0) {
	  add_string(vs1[j]);
	}
      }
    }
  } else {
    //cout << "Here with string: " << s << " " << include_endl << std::endl;
    if (inside_table==false) {
      if (lines.size()==0) {
	//cout << "0 lines, adding " << s << std::endl;
	lines.push_back(s);
	if (include_endl) {
	  lines.push_back("");
	}
      } else if (lines.size()==1) {
	//cout << "1 line, adding " << s << std::endl;
	if (s.length()>0) {
	  if (lines[0].length()==0 || lines[0][lines[0].length()-1]==' ') {
	    lines[0]+=s;
	  } else {
	    lines[0]+=' '+s;
	  }
	}
	if (include_endl) {
	  lines.push_back("");
	}
      } else if (lines.size()==2) {
	//cout << "2 lines, adding " << s << std::endl;
	std::vector<std::string> vs1, vs2;
	split_string(lines[0],vs1);
	split_string(lines[1],vs2);
	size_t c1=vs1.size();
	size_t c2=vs2.size();
	if (false && c1>0 && c1==c2) {
	  inside_table=true;
	  aligns.resize(c1);
	  for(size_t j=0;j<c1;j++) {
	    if (!is_number(vs2[j])) {
	      aligns[j]=columnify::align_left;
	    } else {
	      aligns[j]=columnify::align_dp;
	    }
	  }
	} else {
	  //cout << "Output line: " << lines[0] << std::endl;
	  cout << lines[0] << endl;
	  string stmp=lines[1];
	  lines.resize(1);
	  lines[0]=stmp;
	  if (lines[0].length()==0 || lines[0][lines[0].length()-1]==' ') {
	    lines[0]+=s;
	  } else {
	    lines[0]+=' '+s;
	  }
	  if (include_endl) {
	    lines.push_back("");
	  }
	}
      } else {
	//cout << "More than two lines." << std::endl;
      }
      //for(size_t k=0;k<lines.size();k++) {
      //cout << k << ". " << lines[k] << std::endl;
      //}
    } else {
      // Inside table
    }
  }
  return;
}

auto_table &o2scl_auto_table::operator<<(auto_table &at, double d) {
  string stmp;
  if (size_of_exponent(d)==3) {
    stmp=o2scl::dtos(d,5);
  } else {
    stmp=o2scl::dtos(d,6);
  }
  at.add_string(stmp);
  return at;
}

auto_table &o2scl_auto_table::operator<<(auto_table &at, float f) {
  string stmp;
  stmp=o2scl::dtos(f,6);
  at.add_string(stmp);
  return at;
}

auto_table &o2scl_auto_table::operator<<(auto_table &at, int i) {
  string stmp=o2scl::itos(i);
  at.add_string(stmp);
  return at;
}

auto_table &o2scl_auto_table::operator<<(auto_table &at, char ch) {
  string stmp;
  stmp+=ch;
  at.add_string(stmp);
  return at;
}

auto_table &o2scl_auto_table::operator<<(auto_table &at, size_t s) {
  string stmp=o2scl::szttos(s);
  at.add_string(stmp);
  return at;
}

auto_table &o2scl_auto_table::operator<<(auto_table &at, std::string s) {
  at.add_string(s);
  return at;
}

/*
auto_table &operator<<(auto_table &at,
		       o2scl_auto_table::auto_table &(*f)(o2scl_auto_table::auto_table&)) {
  (*f)(at);
  return at;
}

auto_table &o2scl_auto_table::endl(auto_table &at) {
at.add_string("",true);
return at;
}
*/
