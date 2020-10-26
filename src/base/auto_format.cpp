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

#include <o2scl/auto_format.h>
#include <o2scl/misc.h>
#include <o2scl/string_conv.h>
#include <o2scl/columnify.h>
#include <o2scl/vector.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_auto_format;

auto_format::auto_format() {
  inside_table=false;
  row_max=1000;
  verbose=0;
  auto_tables=true;
  n_headers=0;
  table_lines=0;
}      

void auto_format::done() {
  // Output all lines
  for(size_t j=0;j<lines.size();j++) {
    cout << lines[j] << endl;
  }
  // Clear the output buffer
  lines.clear();
  return;
}

void auto_format::add_string(std::string s, bool include_endl) {

  // If there is no string and no endline, then there is nothing
  // to do
  if (s.length()==0 && include_endl==false) {
    return;
  }

  // Easily handle a single carriage return
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

  // If there is a carriage return mid-string, then use recursion to
  // handle each line separately
  if (s.find('\n')!=std::string::npos) {
    
    // Split into lines
    std::vector<std::string> vs1;
    split_string_delim(s,vs1,'\n');
    // If the original string, s, ends in a '\n', then add another
    // empty string to the end. This ensures that include_endl is set
    // to true.
    if (s[s.length()-1]=='\n') {
      vs1.push_back("");
    }
    // Iterate through each line
    for(size_t j=0;j<vs1.size();j++) {
      // Perform the recursion, setting include_endl to true if there
      // was a '\n' character.
      if (j!=vs1.size()-1) {
	add_string(vs1[j],true);
      } else {
	if (vs1[j].length()>0) {
	  add_string(vs1[j]);
	}
      }
    }
    
  } else {
    
    if (verbose>0) {
      cout << "Start add_string(): \"" << s << "\", "
	   << include_endl << std::endl;
    }
    if (inside_table==false) {
      if (lines.size()==0) {
	if (verbose>1) {
	  cout << "0 lines, adding " << s << std::endl;
	}
	lines.push_back(s);
	if (include_endl) {
	  lines.push_back("");
	}

      } else if (lines.size()==1) {

	if (verbose>1) {
	  cout << "1 line, adding " << s << std::endl;
	}
	// We don't consider adding 's' to the current line unless
	// it is non-empty
	if (s.length()>0) {
	  // If the current line is empty or ends with a space, then
	  // we don't need a space
	  if (lines[0].length()==0 || lines[0][lines[0].length()-1]==' ') {
	    lines[0]+=s;
	  } else {
	    // Add a space
	    lines[0]+=' '+s;
	  }
	}
	if (include_endl) {
	  lines.push_back("");
	}
	
      } else if (lines.size()==2) {

	if (verbose>1) {
	  cout << "2 lines, adding \"" << s << "\"" << std::endl;
	}
	
	// We don't consider adding 's' to the current line unless
	// it is non-empty
	if (s.length()>0) {
	  // If the current line is empty or ends with a space, then
	  // we don't need a space
	  if (lines[1].length()==0 || lines[1][lines[1].length()-1]==' ') {
	    lines[1]+=s;
	  } else {
	    // Add a space
	    lines[1]+=' '+s;
	  }
	}

	// If an endl is requested, then we need to check if we're
	// starting a new table
	if (include_endl) {

	  // Count the number of words in the two lines
	  std::vector<std::string> vs1, vs2;
	  split_string(lines[0],vs1);
	  split_string(lines[1],vs2);
	  size_t c1=vs1.size();
	  size_t c2=vs2.size();
	  if (verbose>0) {
	    std::cout << "First and second line count: " << c1 << " "
		      << c2 << std::endl;
	  }

	  // If they match, then presume we are in a new table
	  if (c1>0 && c1==c2) {
	    
	    if (verbose>0) {
	      std::cout << "Setting inside_table to true." << std::endl;
	    }
	    inside_table=true;

	    // Number of columns
	    size_t ncol=c1;

	    // Determine alignment specifications
	    aligns.resize(ncol);
	    for(size_t j=0;j<ncol;j++) {
	      if (!is_number(vs2[j])) {
		aligns[j]=columnify::align_left;
	      } else {
		aligns[j]=columnify::align_dp;
	      }
	    }
	    if (verbose>0) {
	      cout << "Aligns: ";
	      vector_out(cout,aligns,true);
	    }

	    // Now fill the columns data 
	    columns.clear();
	    columns.resize(ncol);
	    for(size_t j=0;j<ncol;j++) {
	      columns[j].push_back(vs1[j]);
	      columns[j].push_back(vs2[j]);
	    }

	    // Initialize the next_column counter
	    next_column=0;
	    
	  } else {
	    
	    // If we're not in a table, then output the first line in
	    // the cache and move the second line to the first line
	    if (verbose>0) {
	      cout << "Output line: " << lines[0] << std::endl;
	    } else {
	      cout << lines[0] << endl;
	    }
	    lines[0]=lines[1];
	    // Variable include_endl is true, so we're ready to
	    // start a new line
	    lines[1]="";
	    
	  }
	}
      } else {
	O2SCL_ERR("More than two lines in auto_format::add_string().",
		  o2scl::exc_esanity);
      }
      if (verbose>0) {
	std::cout << "Function add_string(), lines object: " << std::endl;
	for(size_t k=0;k<lines.size();k++) {
	  cout << "  " << k << ". \"" << lines[k] << "\"" << std::endl;
	}
	cout << endl;
      }
    } else {
      
      if (verbose>0) {
	cout << "Function add_string(), inside_table: \"" << s << "\", "
	     << include_endl << " " << next_column << std::endl;
	if (false) {
	  for(size_t k=0;k<columns.size();k++) {
	    cout << columns[k].size() << " ";
	  }
	  cout << endl;
	}
      }

      // There are a couple situations to handle:
      
      // 1. We were given a new string in the table, so just add it
      // to the table,
      
      // 2. There is an extra string (s.length()>0 &&
      // next_column>=columns.size()), and thus we have to end the
      // table and go back to non-table mode, and
      
      // 3. include_endl is true even though there were not enough
      // strings for a row (next_column<columns.size()-1),
      // so we have to go back to non-table mode.

      if ((s.length()>0 && next_column>=columns.size()) ||
	  (include_endl && next_column<columns.size()-1)) {

	// Cases 2 and 3 above
	
	if (verbose>0) {
	  cout << "Ending table." << endl;
	}
	inside_table=false;
	if (verbose>0) {
	  cout << "Setting inside_table to false." << endl;
	}
	
	// First, if there is any non-table data, add it to
	// the non-table output buffer
	lines.clear();
	lines.push_back("");
	if (columns[0].size()==0) {
	  O2SCL_ERR("Column 0 empty.",o2scl::exc_einval);
	}
	size_t last_row=columns[0].size()-1;
	for(size_t k=0;k<next_column;k++) {
	  if (columns[k].size()>last_row) {
	    if (verbose>0) {
	      cout << "Adding " << columns[k][last_row]
		   << " to lines[0]." << endl;
	    }
	    if (k!=next_column-1) {
	      lines[0]+=columns[k][last_row]+' ';
	    } else {
	      lines[0]+=columns[k][last_row];
	    }
	  }
	}

	// If there was a string s, then we need to add it as well
	if (s.length()>0) {
	  if (verbose>0) {
	    cout << "Adding " << s << " to lines[0]." << endl;
	  }
	  if (lines[0].length()>0 && lines[0][lines[0].length()-1]!=' ') {
	    lines[0]+=' '+s;
	  } else {
	    lines[0]+=s;
	  }
	}

	if (verbose>0) {
	  cout << "Running columnify::align() " << columns.size() << " "
	       << columns[0].size()-1 << endl;
	  cout << "Column sizes: ";
	  for(size_t k=0;k<columns.size();k++) {
	    cout << columns[k].size() << " ";
	  }
	}
	
	// Now output the table
	columnify c;
	vector<string> tab_out(columns[0].size()-1);
	c.align(columns,columns.size(),columns[0].size()-1,
		tab_out,aligns);
	for(size_t j=0;j<tab_out.size();j++) {
	  cout << tab_out[j] << endl;
	}
	
      } else {
	
	// Inside table, case 1 above
	
	if (s.length()>0) {
	  if (verbose>1) {
	    cout << "Adding \"" << s << "\" to table." << endl;
	  }
	  columns[next_column].push_back(s);
	  next_column++;
	}
	if (include_endl && next_column==columns.size()) {
	  if (verbose>0) {
	    cout << "Resetting next_column." << endl;
	  }
	  next_column=0;
	}
	
      }

      if (verbose>0) {
	cout << endl;
      }
    }
  }
  return;
}

auto_format &o2scl_auto_format::operator<<(auto_format &at, double d) {
  string stmp;
  if (size_of_exponent(d)==3) {
    stmp=o2scl::dtos(d,5);
  } else {
    stmp=o2scl::dtos(d,6);
  }
  at.add_string(stmp);
  return at;
}

auto_format &o2scl_auto_format::operator<<(auto_format &at, float f) {
  string stmp;
  stmp=o2scl::dtos(f,6);
  at.add_string(stmp);
  return at;
}

auto_format &o2scl_auto_format::operator<<(auto_format &at, int i) {
  string stmp=o2scl::itos(i);
  at.add_string(stmp);
  return at;
}

auto_format &o2scl_auto_format::operator<<(auto_format &at, char ch) {
  string stmp;
  stmp+=ch;
  at.add_string(stmp);
  return at;
}

auto_format &o2scl_auto_format::operator<<(auto_format &at, size_t s) {
  string stmp=o2scl::szttos(s);
  at.add_string(stmp);
  return at;
}

auto_format &o2scl_auto_format::operator<<(auto_format &at, std::string s) {
  at.add_string(s);
  return at;
}

/*
  auto_format &operator<<(auto_format &at,
  o2scl_auto_format::auto_format &(*f)(o2scl_auto_format::auto_format&)) {
  (*f)(at);
  return at;
  }

  auto_format &o2scl_auto_format::endl(auto_format &at) {
  at.add_string("",true);
  return at;
  }
*/
