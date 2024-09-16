/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2020-2024, Andrew W. Steiner
  
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
  enabled=true;
  precision_=6;
  outs=&std::cout;
  align_matrices=true;
}      

void auto_format::attach(std::ostream &out) {
  outs=&out;
  return;
}

void auto_format::unattach() {
  done();
  outs=&std::cout;
  return;
}

void auto_format::on() {
  enabled=true;
  return;
}

void auto_format::precision(size_t p) {
  if (p==0) precision_=6;
  else precision_=p;
  return;
}

void auto_format::off() {
  // Make sure that nothing is left in the cache before turning off
  done();
  enabled=false;
  return;
}

void auto_format::done() {
  // Output all lines
  for(size_t j=0;j<lines.size();j++) {
    (*outs) << lines[j] << endl;
  }
  // Clear the output buffer
  lines.clear();
  return;
}

void auto_format::debug_table() {
  (*outs) << "Headers:" << endl;
  for(size_t j=0;j<headers.size();j++) {
    for(size_t k=0;k<headers[j].size();k++) {
      (*outs) << "\"" << headers[j][k] << "\" ";
    }
    (*outs) << endl;
  }
  if (columns.size()>0) {
    (*outs) << "Columns " << columns.size() << endl;
    (*outs) << "Column sizes: ";
    for(size_t j=0;j<columns.size();j++) {
      (*outs) << columns[j].size() << " ";
    }
    (*outs) << endl;
    for(size_t j=0;j<columns[0].size();j++) {
      for(size_t k=0;k<columns.size();k++) {
	(*outs) << "\"" << columns[k][j] << "\" ";
      }
      (*outs) << endl;
    }
  }
  return;
}

void auto_format::start_table() {
  if (enabled==false) return;
  if (lines.size()>0) {
    // If there is still data left in the buffer then output it,
    // which forces an endl before the table even if not explicitly
    // given.
    done();
  }

  n_headers=1;
  inside_table=true;
  
  vector<string> empty;
  headers.clear();
  headers.push_back(empty);

  return;
}

void auto_format::end_table() {
  if (enabled==false) return;

  if (verbose>0) {
    (*outs) << "Running columnify::align() " << columns.size() << " "
	 << columns[0].size() << endl;
    (*outs) << "Column sizes: ";
    for(size_t k=0;k<columns.size();k++) {
      (*outs) << columns[k].size() << " ";
    }
    (*outs) << endl;
  }

  // Determine alignments
  size_t row=n_headers;
  aligns.resize(columns.size()-1);
  for(size_t i=0;i<columns.size()-1;i++) {
    if (!is_number(columns[i][row])) {
      aligns[i]=columnify::align_left;
    } else {
      aligns[i]=columnify::align_dp;
    }
  }
  
  // Now output the table
  columnify c;
  c.table_lines=table_lines;
  vector<string> tab_out(columns[0].size());
  c.align(columns,columns.size()-1,columns[0].size(),
	  tab_out,aligns);

  for(size_t j=0;j<tab_out.size();j++) {
    if (verbose>0) {
      (*outs) << "Output: ";
    }
    if (columns[columns.size()-1][j].length()>0) {
      (*outs) << tab_out[j] << " "
	   << columns[columns.size()-1][j] << endl;
    } else {
      (*outs) << tab_out[j] << endl;
    }
  }
  
  inside_table=false;
  n_headers=0;
  headers.clear();
  columns.clear();
  next_column=0;
  return;
}

void auto_format::add_string(std::string s) {

  if (enabled==false) {
    (*outs) << s;
    return;
  }
  
  bool include_endl=false;
  
  // If there is no string then there is nothing to do
  if (s.length()==0) {
    return;
  }

  // Easily handle a single carriage return
  if (s=="\n") {
    endline();
  }

  // If there is a carriage return mid-string, then use recursion to
  // handle each line separately
  if (s.find('\n')!=std::string::npos) {
    
    // Split into lines
    std::vector<std::string> vs1;
    split_string_delim(s,vs1,'\n');
    
    // Iterate through each line
    for(size_t j=0;j<vs1.size();j++) {
      // Perform the recursion, calling endline() if there was a '\n'
      // character.
      if (j!=vs1.size()-1) {
	add_string(vs1[j]);
	endline();
      } else {
	add_string(vs1[j]);
      }
    }

    return;
  }

  if (inside_table==false) {
    
    if (verbose>0) {
      (*outs) << "Start add_string(): \"" << s << "\"" << std::endl;
    }
  
    if (lines.size()==0) {
      lines.push_back("");
    }

    // If the current line is empty or ends with a space, then
    // we don't need a space
    size_t next_line=lines.size()-1;
    if (lines[next_line].length()==0 ||
	lines[next_line][lines[next_line].length()-1]==' ') {
      lines[next_line]+=s;
    } else {
      // Otherwise, add a space
      lines[next_line]+=' '+s;
    }
    
  } else {

    // Inside table=true section

    if (verbose>0) {
      (*outs) << "Start add_string(), inside_table: \"" << s << "\""
	   << std::endl;
    }

    if (headers.size()==0) {
      O2SCL_ERR("Headers zero in add_string().",o2scl::exc_efailed);
    }

    // If columns is empty, then we're in the header section still
    if (columns.size()==0) {
    
      // Add the string to the appropriate header row
      headers[headers.size()-1].push_back(s);

    } else {

      // Main column section

      if (verbose>0) {
	(*outs) << "Adding \"" << s << "\" to column " << next_column
	     << std::endl;
      }

      size_t n_cols=headers[0].size();

      if (next_column>n_cols) {
	size_t n=columns[n_cols].size();
	columns[n_cols][n-1]+=" "+s;
      } else {
	columns[next_column].push_back(s);
      }
      next_column++;

      if (verbose>0) {
	(*outs) << "Next column is " << next_column << std::endl;
      }
      
    }

#ifdef O2SCL_NEVER_DEFINED

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
	(*outs) << "Ending table." << endl;
      }
      inside_table=false;
      if (verbose>0) {
	(*outs) << "Setting inside_table to false." << endl;
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
	    (*outs) << "Adding " << columns[k][last_row]
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
	  (*outs) << "Adding " << s << " to lines[0]." << endl;
	}
	if (lines[0].length()>0 && lines[0][lines[0].length()-1]!=' ') {
	  lines[0]+=' '+s;
	} else {
	  lines[0]+=s;
	}
      }
      
      if (verbose>0) {
	(*outs) << "Running columnify::align() " << columns.size() << " "
	     << columns[0].size()-1 << endl;
	(*outs) << "Column sizes: ";
	for(size_t k=0;k<columns.size();k++) {
	  (*outs) << columns[k].size() << " ";
	}
      }
      
      // Now output the table
      columnify c;
      vector<string> tab_out(columns[0].size()-1);
      c.align(columns,columns.size(),columns[0].size()-1,
	      tab_out,aligns);
      for(size_t j=0;j<tab_out.size();j++) {
	(*outs) << tab_out[j] << endl;
      }
      
    } else {
      
      // Inside table, case 1 above
      
      if (s.length()>0) {
	if (verbose>1) {
	  (*outs) << "Adding \"" << s << "\" to table." << endl;
	}
	columns[next_column].push_back(s);
	next_column++;
      }
      if (include_endl && next_column==columns.size()) {
	if (verbose>0) {
	  (*outs) << "Resetting next_column." << endl;
	}
	next_column=0;
      }
      
    }
    
    if (verbose>0) {
      (*outs) << endl;
    }

#endif    
    
  }
  
  return;
}

void auto_format::endline() {

  if (enabled==false) {
    (*outs) << endl;
    return;
  }
  
  if (inside_table==false) {

    if (verbose>0) {
      (*outs) << "Output: ";
    }
    if (lines.size()>0) {
      (*outs) << lines[0] << endl;
    } else {
      (*outs) << endl;
    }
    lines.clear();

  } else {

    // Inside the table section

    if (columns.size()==0) {

      // Header section

      // If necessary, add a new header
      if (headers.size()<n_headers) {
	if (verbose>0) {
	  (*outs) << "Adding a new header row." << std::endl;
	}
	vector<string> empty;
	headers.push_back(empty);
	
      } else {
	
	// Otherwise, move to the body section
	size_t n_cols=headers[0].size();
	if (verbose>0) {
	  (*outs) << "Moving to body section, n_cols="
		    << n_cols << std::endl;
	}
	vector<string> empty;
	for (size_t j=0;j<n_cols+1;j++) {
	  columns.push_back(empty);
	}
	next_column=0;

	// Copy headers to the columns
	for(size_t j=0;j<headers.size();j++) {
	  //cout << "Here: " << j << " " << n_cols+1 << " "
	  //<< headers[j].size() << endl;
	  for(size_t k=0;k<n_cols+1 || k<headers[j].size();k++) {
	    if (k<n_cols+1) {
	      if (k<headers[j].size()) {
		//cout << "Iere: " << k << " " << j << " "
		//<< headers[j][k] << endl;
		columns[k].push_back(headers[j][k]);
	      } else {
		//cout << "Jere: " << k << endl;
		columns[k].push_back("");
	      }
	    } else {
	      columns[n_cols][j]+=" "+headers[j][k];
	    }
	  }
	}
      }
      
    } else {

      if (columns.size()==row_max) {
	if (verbose>0) {
	  (*outs) << "Ending table." << std::endl;
	}
	end_table();
      }

      if (next_column!=0) {
	// If some columns were not filled, then fill them
	// before proceeding to the next line
	size_t n_cols=headers[0].size();
	while (next_column<=n_cols) {
	  columns[next_column].push_back("");
	  next_column++;
	}
	
	next_column=0;
      }
      
    }
  }
    
#ifdef NEVER_DEFINED
  
    // Count the number of words in the two lines
    std::vector<std::string> vs1, vs2;
    split_string(lines[0],vs1);
    split_string(lines[1],vs2);
    size_t c1=vs1.size();
    size_t c2=vs2.size();
    if (verbose>0) {
      (*outs) << "First and second line count: " << c1 << " "
		<< c2 << std::endl;
    }
    
    // If they match, then presume we are in a new table
    if (c1>0 && c1==c2) {
      
      if (verbose>0) {
	(*outs) << "Setting inside_table to true." << std::endl;
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
      
    }
  }

#endif
  
  return;
}

auto_format &o2scl_auto_format::operator<<(auto_format &at, double d) {
  string stmp;
  if (size_of_exponent(d)==3 && at.precision_>1) {
    stmp=o2scl::dtos(d,at.precision_-1);
  } else {
    stmp=o2scl::dtos(d,at.precision_);
  }
  at.add_string(stmp);
  return at;
}

auto_format &o2scl_auto_format::operator<<(auto_format &at,
					   long double d) {
  string stmp;
  stmp=o2scl::dtos(d,at.precision_);
  at.add_string(stmp);
  return at;
}

#ifndef O2SCL_NO_BOOST_MULTIPRECISION
auto_format &o2scl_auto_format::operator<<(auto_format &at,
					   const cpp_dec_float_35 &d) {
  string stmp;
  stmp=o2scl::dtos(d,at.precision_);
  at.add_string(stmp);
  return at;
}

auto_format &o2scl_auto_format::operator<<(auto_format &at,
					   const cpp_dec_float_50 &d) {
  string stmp;
  stmp=o2scl::dtos(d,at.precision_);
  at.add_string(stmp);
  return at;
}

auto_format &o2scl_auto_format::operator<<(auto_format &at,
					   const cpp_dec_float_100 &d) {
  string stmp;
  stmp=o2scl::dtos(d,at.precision_);
  at.add_string(stmp);
  return at;
}
#endif

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

auto_format &o2scl_auto_format::operator<<(auto_format &at,
					   const std::vector<double> &vd) {
  for(size_t i=0;i<vd.size();i++) {
    at << vd[i];
  }
  return at;
}

auto_format &o2scl_auto_format::operator<<(auto_format &at,
					   const std::vector<int> &vi) {
  for(size_t i=0;i<vi.size();i++) {
    at << vi[i];
  }
  return at;
}

auto_format &o2scl_auto_format::operator<<(auto_format &at,
					   const std::vector<size_t> &vi) {
  for(size_t i=0;i<vi.size();i++) {
    at << vi[i];
  }
  return at;
}

auto_format &o2scl_auto_format::operator<<(auto_format &at,
					   const std::vector<char> &vi) {
  for(size_t i=0;i<vi.size();i++) {
    at << vi[i];
  }
  return at;
}

auto_format &o2scl_auto_format::operator<<(auto_format &at,
					   const std::vector<std::string> &vi) {
  for(size_t i=0;i<vi.size();i++) {
    at << vi[i];
  }
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
