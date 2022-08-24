/*
  -------------------------------------------------------------------
  
  Copyright (C) 2020-2022, Andrew W. Steiner
  
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
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <o2scl/string_conv.h>
#include <o2scl/vector.h>

using namespace std;
using namespace o2scl;

/** \brief Convert all non-alphanumeric characters to underscores
 */
std::string underscoreify(std::string s) {
  std::string s2;
  bool last_char_is_underscore=false;
  for(size_t i=0;i<s.length();i++) {
    if (std::isalnum(s[i])==false) {
      if (i==0 || last_char_is_underscore==false) {
        s2+='_';
      }
      last_char_is_underscore=true;
    } else {
      s2+=s[i];
      last_char_is_underscore=false;
    }
  }
  return s2;
}

/** \brief Read the next line from \c fin, skip blank and comment 
    lines, split by spaces into \c vs, and set \c done if done
*/
void next_line(ifstream &fin, std::string &line,
               std::vector<std::string> &vs, bool &done) {
  if (fin.eof()) {
    done=true;
    return;
  }
  getline(fin,line);
  while (line.length()==0 || line[0]=='#') {
    if (fin.eof()) {
      done=true;
      return;
    } else {
      getline(fin,line);
    }
  }
  vs.clear();
  split_string(line,vs);
  return;
}

/** \brief Parse a potentially multi-line string in the interface
    file

    This is used for rst_header, py_class_doc and extra_py. 
*/
void parse_vector_string(ifstream &fin, std::string type,
                         std::string &line,
                         std::vector<std::string> &vs, bool &done,
                         std::vector<std::string> &output) {

  bool debug=false;
  if (debug) {
    cout << "parse_vector_string() start: " << line << endl;
  }
  if (vs.size()==1) {
    if (debug) {
      cout << "Clearing " << type << "." << endl;
    }
    output.clear();
    next_line(fin,line,vs,done);
  } else if ((vs[0]==type && vs[1]=="|") ||
             (vs.size()>=3 && vs[1]==type && vs[2]=="|")) {
    if (debug) {
      cout << "Here " << line << endl;
    }
    output.clear();
    next_line(fin,line,vs,done);
    while (line.length()>0 && line[0]=='|') {
      if (line[1]==' ' && line.length()>2) {
        line=line.substr(2,line.length()-2);
      } else {
        line=line.substr(1,line.length()-1);
      }
      output.push_back(line);
      next_line(fin,line,vs,done);
    }
  } else {
    if (debug) {
      cout << "Here2 " << line << endl;
    }
    output.clear();
    output.push_back("");
    for(size_t k=1;k<vs.size();k++) {
      output[0]+=vs[k];
      if (k!=vs.size()-1) output[0]+=" ";
    }
    next_line(fin,line,vs,done);
  }
  /*
    Remove surrounding quotes
    
    if (output[0]=='\"' &&
    output[output.length()-
    1]=='\"') {
    output=output.substr
    (1,output.length()-2);
    }
  */
  cout << "Setting " << type << " to:" << endl;
  for(size_t i=0;i<output.size();i++) {
    cout << i << ": " << output[i] << endl;
  }
  
  return;
}

/** \brief Interface base
 */
class if_base {
  
public:

  /// Object name
  std::string name;

  /// Python name for the object
  std::string py_name;
  
};

/** \brief Interface for shared pointer
 */
class if_shared_ptr : public if_base {

public:

  /// Namespace
  std::string ns;
  
};

/** \brief Interface type

    \note The value of \c name from the parent class is used to store
    the type itself, \c prefix is used to store prefixes like "const"
    or "static" and \c suffix is used to store suffixes like "*", or
    "&".
 */
class if_type : public if_base {
  
public:
  
  /** \brief Type qualifiers

      Includes either
      - const
      - static 
      - short
      - long
      - io (yanic-specific)
      - out (yanic-specific)

      and 

      The shared_ptr is also used as a prefix to indicate a 
      shared pointer of a particular type.
  */
  std::string prefix;
  
  /// Ampersands and asterisks, *, *&, &, **, etc.
  std::string suffix;

  /** \brief Convert the type to a string for output
      
      Prefix and suffix are surrounded in parenthesis to make
      them more clear
  */
  std::string to_string() {
    std::string s=this->name;
    if (prefix.length()>0) s=((std::string)"(")+prefix+") "+s;
    if (suffix.length()>0) s+=" ("+suffix+")";
    return s;
  }
  
  /** \brief Parse a vector string object as a type with a prefix,
      a suffix, and a name
   */
  void parse(std::vector<std::string> &vs, size_t start, size_t end) {

    if (start>=end) {
      O2SCL_ERR("if_type::parse() called with start>=end.",
                o2scl::exc_einval);
    }
    if (vs[end-1][0]=='[' &&
        vs[end-1][vs[end-1].length()-1]==']') {
    }
    
    // Create a new vector string object for convenience
    std::vector<std::string> vs2;
    for(size_t i=start;i<end;i++) {
      vs2.push_back(vs[i]);
    }

    // Parse according to the number of arguments in the list.
    // We will handle *'s or &'s in the 'name' field later.
    if (vs2.size()==1) {
      
      prefix="";
      suffix="";
      // If there is only one string in the list with no spaces, then
      // store the type in name
      name=vs2[0];
      
    } else if (vs2.size()==2) {

      // If's a two-string list, it's either "prefix" "name" or
      // "name" "suffix", where the suffix is made up only of
      // *'s and &'s. 
      
      bool symb_only=true;
      for(size_t i=0;i<vs2[1].length();i++) {
        if (vs2[1][i]!='*' && vs2[1][i]!='&') symb_only=false;
      }
      if (symb_only) {
        name=vs2[0];
        suffix=vs2[1];
      } else {
        prefix=vs2[0];
        name=vs2[1];
      }
      
    } else if (vs2.size()==3) {

      // If it's a three-string list, then it could be a shared
      // pointer, a type of the form "prefix" "name" "suffix", or a
      // type "prefix 1" "prefix2" "name". Note that this function
      // cannot handle, e.g. "double * *".
      
      if (vs2[0]=="shared_ptr" || vs2[0]=="std::shared_ptr") {
        prefix=vs2[0];
        name=vs2[1];
        suffix=vs2[2];
      } else {
        bool symb_only=true;
        for(size_t i=0;i<vs2[2].length();i++) {
          if (vs2[2][i]!='*' && vs2[2][i]!='&') symb_only=false;
        }
        if (symb_only) {
          prefix=vs2[0];
          name=vs2[1];
          suffix=vs2[2];
        } else {
          prefix=vs2[0]+" "+vs2[1];
          name=vs2[2];
        }
      }
      
    } else if (vs2.size()==4) {
      
      bool symb_only=true;
      for(size_t i=0;i<vs2[3].length();i++) {
        if (vs2[3][i]!='*' && vs2[3][i]!='&') symb_only=false;
      }
      if (symb_only) {
        prefix=vs2[0]+" "+vs2[1];
        name=vs2[2];
        suffix=vs2[3];
      } else {
        prefix=vs2[0]+" "+vs2[1]+" "+vs2[2];
        name=vs2[3];
      }
      
    } else {
      
      cerr << "Unsupported number of type arguments, " << vs2.size()
           << ", in if_type::parse()." << endl;
      cout << "vs: ";
      vector_out(cout,vs,true);
      cout << "vs2: ";
      vector_out(cout,vs2,true);
      O2SCL_ERR("Unsupported type format in if_type::parse().",
                o2scl::exc_eunimpl);
      
    }
    
    // Take *'s and &'s in name and add them to suffix
    char last_ch=name[name.length()-1];
    while (last_ch=='*' && last_ch=='&') {
      string stemp;
      stemp+=last_ch;
      suffix=stemp+suffix;
      if (name.length()==1) {
        last_ch=' ';
      } else {
        name.substr(0,name.length()-1);
        last_ch=name[name.length()-1];
      }
    }

    return;
  }

  /// Return true if the type is a pointer
  bool is_pointer() {
    if (suffix=="*") return true;
    return false;
  }
  
  /// Return true if the type is a reference
  bool is_reference() {
    if (suffix=="&") return true;
    return false;
  }
  
  /// Return true if the type is const
  bool is_const() {
    if (prefix.find("const")!=std::string::npos) return true;
    return false;
  }
  
  /// Return true if the type is const
  bool is_shared_ptr() {
    if (prefix.find("shared_ptr")!=std::string::npos ||
        prefix.find("std::shared_ptr")!=std::string::npos) {
      return true;
    }
    return false;
  }
  
  /// Return true if the type is static
  bool is_static() {
    if (prefix.find("static")!=std::string::npos) return true;
    return false;
  }

  /// Return true if the type is long
  bool is_long() {
    if (prefix.find("long")!=std::string::npos) return true;
    return false;
  }

  /// Return true if the type is an input/output reference
  bool is_io() {
    if (prefix.find("io")!=std::string::npos) return true;
    return false;
  }

  /// Return true if the type is an output reference
  bool is_out() {
    if (prefix.find("out")!=std::string::npos) return true;
    return false;
  }

  /// Return true if the type is a standard C type
  bool is_ctype() {
    if (name=="bool" || name=="double" || name=="int" ||
        name=="char" || name=="size_t" || name=="float") {
      return true;
    }
    return false;
  }

  // End of class if_type
};

/** \brief Interface variable

    The name of the variable is stored in \c name, the type is
    stored in \c type, and a default value (if present) is 
    stored in \c value.
 */
class if_var : public if_base {
  
public:

  /// The variable type 
  if_type ift;

  /// The value (used for default function values)
  std::string value;

  /// The python name for the variable
  std::string py_name;
  
  /** \brief Parse a list of strings to this variable

      Note that this does not set the python name because that 
      is specified on a different line in the interface file.
   */
  void parse(vector<string> &vs) {

    if (vs.size()==0) {
      O2SCL_ERR("if_var.parse() called with empty list.",
                o2scl::exc_esanity);
    }
    
    // Check if ampersands or asterisks are on the LHS
    // of the variable name. If so, include them in
    // the variable type instead

    // Temporarily store the last string in the list
    string last_string=vs[vs.size()-1];

    // This will be the last index in the list which does not
    // contain a default value
    size_t last_non_def=vs.size()-1;

    // See if we can find a default value, if so, set the "value" field
    if (last_string[0]=='[' && last_string[last_string.length()-1]==']') {
      value=last_string.substr(1,last_string.length()-2);
      if (vs.size()==0) {
        O2SCL_ERR("if_var.parse() called with default value but no type.",
                  o2scl::exc_esanity);
      }
      last_string=vs[vs.size()-2];
      last_non_def--;
    }
    
    if (last_string[0]=='&' || last_string[0]=='*') {
      
      // Clear the last list entry so we can fill it with the
      // ampersands and asterisks we find in 'last_string'
      vs[last_non_def]="";

      // Progressively move ampersands and asterisks from last_string
      // to the last string in the list
      while (last_string[0]=='&' || last_string[0]=='*') {
        vs[last_non_def]=last_string[0]+vs[last_non_def];
        last_string=last_string.substr(1,last_string.length()-1);
      }
      
      name=last_string;
      ift.parse(vs,1,last_non_def+1);
      
    } else {
      
      name=last_string;
      ift.parse(vs,1,last_non_def);
    }
    
    return;
  }
  
};

/** \brief Interface function
 */
class if_func : public if_base {
  
public:
  
  /// Return value 
  if_type ret;

  /// Function arguments
  std::vector<if_var> args;
  
  /// Namespace
  std::string ns;

  /// If true, then another function has the same name
  bool overloaded;

  if_func() {
    overloaded=false;
  }
};

/** \brief Interface class
 */
class if_class : public if_base {
  
public:

  /// Pattern for the class documentation in python
  std::vector<std::string> py_class_doc;
  
  /// True if the class is abstract
  bool is_abstract;

  /// If true, the standard copy constructors are included
  bool std_cc;

  /// If true, then the default constructor is included (default true)
  bool def_cons;
  
  /// Members
  std::vector<if_var> members;
  
  /// Methods
  std::vector<if_func> methods;
  
  /// Constructors
  std::vector<if_func> cons;
  
  /// List of parent classes
  std::vector<std::string> parents;
  
  /// Lines of extra python code for the python class
  std::vector<std::string> extra_py;
  
  /// Namespace
  std::string ns;

  if_class() {
    is_abstract=false;
    std_cc=false;
    def_cons=true;
  }
  
};

int main(int argc, char *argv[]) {

  if (argc<5) {
    cerr << "Args: <interface file> <c++ output prefix> "
         << "<python prefix> <rst prefix> [header file]" << endl;
    exit(-1);
  }
  
  std::string fname=argv[1];
  std::string cpp_prefix=argv[2];
  std::string py_prefix=argv[3];
  std::string rst_prefix=argv[4];

  std::string header_file="";

  cout << "Reading interface file " << fname << " ." << endl;
  cout << "Setting C++ output prefix to " << cpp_prefix << " ." << endl;
  cout << "Setting Python output prefix to " << py_prefix << " ." << endl;
  if (argc>=6) {
    header_file=argv[5];
    cout << "Setting header file to " << header_file << " ." << endl;
  }
  cout << endl;

  // Read header file
  vector<string> header_strings;
  ifstream fin2;
  string line_temp;
  fin2.open(header_file.c_str());
  while (getline(fin2,line_temp)) {
    header_strings.push_back(line_temp);
    cout << "Line: " << line_temp << endl;
  }
  fin2.close();
  
  cout << "Parsing interface " << fname << " ." << endl;

  // The internal version of the C++/python interface is stored in
  // these data structures. The parser reads the interface file into
  // this data
  
  // Current namespace
  std::string ns;
  // Current dll_name
  std::string dll_name;
  // Current rst_header
  std::vector<std::string> rst_header;
  /// Python documentation pattern
  std::vector<std::string> py_class_doc;
  // Current list of includes
  std::vector<std::string> h_includes;
  // Current list of includes
  std::vector<std::string> cpp_includes;
  // Current list of includes
  std::vector<std::string> py_headers;
  // Current list of includes
  std::vector<std::string> cpp_using;
  // Current list of classes
  std::vector<if_class> classes;
  // Current list of functions
  std::vector<if_func> functions;
  // Current list of shared pointers
  std::vector<if_shared_ptr> sps;

  // Mapping classes to python names for internal classes
  std::map<std::string,std::string,std::less<std::string>> class_py_names;
  typedef
    std::map<std::string,std::string,std::less<std::string>>::iterator
    cpn_it;
  
  class_py_names.insert(std::make_pair("std::string","std_string"));
  class_py_names.insert(std::make_pair("std::vector<double>",
                                       "std_vector"));
  class_py_names.insert(std::make_pair("std::vector<int>",
                                       "std_vector_int"));
  class_py_names.insert(std::make_pair("std::vector<size_t>",
                                       "std_vector_size_t"));
  class_py_names.insert(std::make_pair("std::vector<std::string>",
                                       "std_vector_string"));
  class_py_names.insert(std::make_pair("std::vector<std::vector<std::string>>",
                                       "vec_vec_string"));
  class_py_names.insert(std::make_pair("std::vector<std::vector<double>>",
                                       "std_vector_vector"));
  class_py_names.insert(std::make_pair("tensor<>","tensor"));
  class_py_names.insert(std::make_pair
                        ("boost::numeric::ublas::vector<double>",
                         "ublas_vector"));
  class_py_names.insert(std::make_pair
                        ("boost::numeric::ublas::matrix<double>",
                         "ublas_matrix"));
  class_py_names.insert(std::make_pair
                        ("boost::numeric::ublas::vector<int>",
                         "ublas_vector_int"));
  class_py_names.insert(std::make_pair
                        ("boost::numeric::ublas::matrix<int>",
                         "ublas_matrix_int"));
  class_py_names.insert(std::make_pair("std::vector<contour_line>",
                                       "vector_contour_line"));
  
  // Open file
  ifstream fin;
  fin.open(fname.c_str());
  
  // Parse file
  std::string line;
  std::vector<std::string> vs;
  bool done=false;

  // If true, import numpy
  bool import_numpy=false;

  // Read the first line
  next_line(fin,line,vs,done);
  if (done==true) {
    cerr << "Empty interface file." << endl;
    exit(-1);
  }

  // Continue parsing, line by line, until 'done' is true
  while (done==false) {
    
    if (vs[0]=="namespace") {
      
      if (vs.size()==1) {
        ns="";
        cout << "Clearing namespace." << endl;
      } else {
        ns=vs[1];
        cout << "Setting namespace to " << ns << "." << endl;
      }
      
      next_line(fin,line,vs,done);
      
    } else if (vs[0]=="dll_name") {
      
      if (vs.size()==1) {
        dll_name="";
        cout << "Clearing dll_name." << endl;
      } else {
        dll_name=vs[1];
        cout << "Setting dll_name to " << dll_name << "." << endl;
      }
      
      next_line(fin,line,vs,done);
      
    } else if (vs[0]=="rst_header") {
      
      parse_vector_string(fin,"rst_header",line,vs,done,rst_header);
      
    } else if (vs[0]=="py_class_doc") {

      parse_vector_string(fin,"py_class_doc",line,vs,done,py_class_doc);
      
    } else if (vs[0]=="h_include") {
      
      if (vs.size()<2) {
        cerr << "No argument for h_include." << endl;
        exit(-1);
      }
      h_includes.push_back(vs[1]);
      
      next_line(fin,line,vs,done);
      
    } else if (vs[0]=="py_header") {
      
      if (vs.size()<2) {
        cerr << "No argument for py_header." << endl;
        exit(-1);
      }
      std::string pyh=vs[1];
      for(size_t j=2;j<vs.size();j++) {
        pyh+=" "+vs[j];
      }
      py_headers.push_back(pyh);
      
      next_line(fin,line,vs,done);
      
    } else if (vs[0]=="cpp_include") {
      
      if (vs.size()<2) {
        cerr << "No argument for cpp_include." << endl;
        exit(-1);
      }
      cpp_includes.push_back(vs[1]);
      
      next_line(fin,line,vs,done);
      
    } else if (vs[0]=="shared_ptr") {
      
      if (vs.size()<2) {
        cerr << "No argument for shared_ptr." << endl;
        exit(-1);
      }
      
      if_shared_ptr ifss;
      ifss.ns=ns;
      ifss.name=vs[1];
      
      cout << "Shared pointer for type " << vs[1] << endl;
      
      next_line(fin,line,vs,done);

      // If present, record the python name of the shared pointer
      if (line[0]=='-' && line[1]==' ' && vs[1]=="py_name" &&
          vs.size()==3) {
        ifss.py_name=vs[2];
        cout << "  with python name " << vs[2] << endl;
        next_line(fin,line,vs,done);
      }
      
      sps.push_back(ifss);
      cout << endl;
      
    } else if (vs[0]=="cpp_using") {
      
      if (vs.size()<2) {
        cerr << "No argument for cpp_using." << endl;
        exit(-1);
      }
      cpp_using.push_back(vs[1]);
      
      next_line(fin,line,vs,done);
      
    } else if (vs[0]=="class") {
      
      if (vs.size()<2) {
        cerr << "No name for class." << endl;
        exit(-1);
      }
      
      if_class ifc;
      ifc.name=vs[1];
      ifc.ns=ns;
      ifc.py_class_doc=py_class_doc;

      // Determine if the class is abstract
      if (vs.size()>=3 && vs[2]=="abstract") {
        ifc.is_abstract=true;
        cout << "Starting abstract class " << ifc.name << " in namespace "
             << ifc.ns << endl;
      } else {
        cout << "Starting class " << ifc.name << " in namespace "
             << ifc.ns << endl;
      }

      // Read the next line
      next_line(fin,line,vs,done);

      // Loop until the class definition is finished
      bool class_done=false;
      while (class_done==false) {

        if (vs.size()>=3 && vs[0]=="-" && vs[1]=="function") {
          
          // Section for a member function
          
          if_func iff;
          iff.name=vs[2];
          cout << "  Starting member function " << vs[2] << endl;

          // Read the return value
          
          next_line(fin,line,vs,done);
          
          if (done==true || vs.size()<2 || vs[0]!="-" || line[0]!=' ' ||
              line[1]!=' ') {
            cerr << "Could not get return value for function "
                 << iff.name << " in class " << ifc.name << endl;
            exit(-1);
          }
          
          iff.ret.parse(vs,1,vs.size());

          cout << "    Member function " << iff.name
               << " has return type "
               << iff.ret.to_string() << endl;

          // If it returns a vector<double>, then we need to
          // import numpy so we can return a numpy array
          if ((iff.ret.name=="vector<double>" ||
               iff.ret.name=="std::vector<double>") &&
              iff.ret.suffix=="&") {
            import_numpy=true;
          }

          // Next line after return type
          next_line(fin,line,vs,done);

          // Check for a python name for this member function
          if (vs[1]=="py_name" && vs.size()>=3) {
            
            iff.py_name=vs[2];

            cout << "    Member function " << iff.name
                 << " has py_name " << iff.py_name << endl;

            next_line(fin,line,vs,done);
            if (done) class_done=true;
            
          }

          // Read the member function arguments if present
          while (vs.size()>=2 && line[0]==' ' && line[1]==' ' &&
                 vs[0]=="-") {
            
            if_var ifv;

            ifv.parse(vs);

            cout << "    Member function " << iff.name
                 << " has argument " << ifv.name << " with type "
                 << ifv.ift.to_string();
            if (ifv.value.length()!=0) {
              cout << " and default value " << ifv.value;
            }
            cout << " [";
            if (ifv.ift.is_pointer()) cout << "p";
            if (ifv.ift.is_reference()) cout << "r";
            if (ifv.ift.is_const()) cout << "t";
            if (ifv.ift.is_static()) cout << "s";
            if (ifv.ift.is_shared_ptr()) cout << "h";
            if (ifv.ift.is_io()) cout << "i";
            if (ifv.ift.is_out()) cout << "o";
            if (ifv.ift.is_out()) cout << "o";
            if (ifv.ift.is_ctype()) cout << "c";
            cout << "]." << endl;

            // Check that C type references which are not return values
            // are marked 'out' or 'io'
            if (iff.name!="operator[]" && iff.name!="operator()" &&
                ifv.ift.is_ctype() && ifv.ift.is_reference()) {
              if (!ifv.ift.is_out() && !ifv.ift.is_io()) {
                O2SCL_ERR("C type reference must be marked 'out' or 'io'",
                          o2scl::exc_einval);
              }
            }
            
            iff.args.push_back(ifv);
            
            next_line(fin,line,vs,done);
            if (done) class_done=true;

            // Here, check if there is a python name to this member
            // function argument

            if (vs.size()>=3 && vs[0]=="-" && vs[1]=="py_name" &&
                line[0]==' ' && line[1]==' ' && line[2]==' ' &&
                line[3]==' ') {
              ifv.py_name=vs[2];
              cout << "      Member function argument "
                   << ifv.name << " has python name "
                   << ifv.py_name << " ." << endl;
            }
          }

          ifc.methods.push_back(iff);

          // End of member function section
          
        } else if (vs.size()>=3 && vs[0]=="-" && vs[1]=="cons") {

          // Section for a constructor
          
          if_func iff;
          iff.name=vs[2];
          cout << "  Starting constructor to be named " << iff.name
               << " ." << endl;

          // Read the next line
          next_line(fin,line,vs,done);

          // Check for a python name for this member function
          if (vs[1]=="py_name" && vs.size()>=3) {
            
            iff.py_name=vs[2];

            cout << "    Member function " << iff.name
                 << " has py_name " << iff.py_name << endl;

            next_line(fin,line,vs,done);
            if (done) class_done=true;
            
          }

          // Constructors don't have return values, so just
          // look for constructor arguments
          while (vs.size()>=2 && line[0]==' ' && line[1]==' ' &&
                 vs[0]=="-") {
            
            if_var ifv;

            ifv.parse(vs);
            
            cout << "    Constructor with name " << iff.name
                 << " has argument " << ifv.name << " with type "
                 << ifv.ift.to_string() << endl;
            
            iff.args.push_back(ifv);
            
            next_line(fin,line,vs,done);
            if (done) class_done=true;
          }

          ifc.cons.push_back(iff);

          // End of constructor section
          
        } else if (vs.size()>=3 && vs[0]=="-" && vs[1]=="parent") {

          // Parent class 
          
          ifc.parents.push_back(vs[2]);
          cout << "  Parent class " << vs[2] << endl;
          
          next_line(fin,line,vs,done);
          if (done) class_done=true;
          
        } else if (vs.size()>=2 && vs[0]=="-" && vs[1]=="no_def_cons") {

          // No default constructor flag
          
          if (ifc.std_cc==true) {
            O2SCL_ERR2("Standard copy constructor without default ",
                       "constructor not implemented.",o2scl::exc_eunimpl);
          }
          
          ifc.def_cons=false;
          cout << "  Parent class does not have a default constructor."
               << endl;
          
          next_line(fin,line,vs,done);
          if (done) class_done=true;
          
        } else if (vs.size()>=2 && vs[0]=="-" && vs[1]=="std_cc") {

          // Standard copy constructor flag

          if (ifc.def_cons==false) {
            O2SCL_ERR2("Standard copy constructor without default ",
                       "constructor not implemented.",o2scl::exc_eunimpl);
          }
          
          ifc.std_cc=true;
          cout << "  Class " << ifc.name << " has the standard "
               << "copy constructors." << endl;
          
          next_line(fin,line,vs,done);
          if (done) class_done=true;
          
        } else if (vs.size()>=3 && vs[0]=="-" && vs[1]=="py_name") {

          // Python name for the class
          
          ifc.py_name=vs[2];
          cout << "  Class " << ifc.name
               << " has python name " << vs[2] << endl;
          
          next_line(fin,line,vs,done);
          if (done) class_done=true;
          
        } else if (vs.size()>=3 && vs[0]=="-" &&
                   vs[1]=="py_class_doc") {

          // Python documentation for this class
          
          parse_vector_string(fin,"py_class_doc",line,vs,done,
                              ifc.py_class_doc);
          
          if (done) class_done=true;
          
        } else if (vs.size()>=3 && vs[0]=="-" &&
                   vs[1]=="extra_py") {

          // Extra python code for this class

          parse_vector_string(fin,"extra_py",line,vs,done,ifc.extra_py);
          
          if (done) class_done=true;
          
        } else if (vs.size()>=3 && vs[0]=="-") {

          // Member data section
          
          if_var ifv;
          
          ifv.parse(vs);
          
          cout << "  Member " << ifv.name << " with type "
               << ifv.ift.to_string() << " [";
          if (ifv.ift.is_pointer()) cout << "p";
          if (ifv.ift.is_reference()) cout << "r";
          if (ifv.ift.is_const()) cout << "t";
          if (ifv.ift.is_static()) cout << "s";
          if (ifv.ift.is_shared_ptr()) cout << "h";
          if (ifv.ift.is_io()) cout << "i";
          if (ifv.ift.is_out()) cout << "o";
          if (ifv.ift.is_out()) cout << "o";
          if (ifv.ift.is_ctype()) cout << "c";
          cout << "]." << endl;

          next_line(fin,line,vs,done);
          
          if (done) {
            class_done=true;
          } else {

            // Check if the member data has a py_name
            if (line[0]==' ' && line[1]==' ' && line[2]=='-' &&
                line[3]==' ' && vs[1]=="py_name" && vs.size()>=3) {
              
              ifv.py_name=vs[2];
              cout << "    Member " << ifv.name << " has py_name "
                   << ifv.py_name << " ." << endl;
              
              next_line(fin,line,vs,done);
              if (done) {
                class_done=true;
              }
            }
          }
          
          ifc.members.push_back(ifv);

          // End of member data section

        } else {

          // If it failed to match a part of the class, then the
          // class is done
          class_done=true;
          
        }

        if (class_done==true) {
          cout << "Class " << ifc.name << " done." << endl;
          cout << endl;
          classes.push_back(ifc);
        }

        // End of 'while (class_done==false)'
      }
      
    } else if (vs[0]=="function") {

      if_func iff;
      iff.name=vs[1];
      iff.ns=ns;
      cout << "Starting function " << vs[1] << endl;

      // Read and parse function return type
      
      next_line(fin,line,vs,done);
          
      if (done==true || vs.size()<2 || vs[0]!="-" || line[0]!='-' ||
          line[1]!=' ') {
        cerr << "Could not get return value for function "
             << iff.name << endl;
        exit(-1);
      }
          
      iff.ret.parse(vs,1,vs.size());
      
      cout << "  Function " << iff.name
           << " has return type "
           << iff.ret.to_string() << endl;

      // Read next line
      
      next_line(fin,line,vs,done);

      // Determine if a python name is specified
      
      if (vs[1]=="py_name" && vs.size()>=3) {
        
        iff.py_name=vs[2];
        
        cout << "  Function " << iff.name
             << " has py_name " << iff.py_name << endl;
        next_line(fin,line,vs,done);
        
      }

      // Read function arguments
      
      bool function_done=false;
      
      while (vs.size()>=2 && line[0]=='-' && line[1]==' ' &&
             vs[0]=="-") {

        if_var ifv;

        ifv.parse(vs);

        if (ifv.value.length()==0) {
          cout << "  Function " << iff.name
               << " has argument " << ifv.name << " with type "
               << ifv.ift.to_string() << endl;
        } else {
          cout << "  Function " << iff.name
               << " has argument " << ifv.name << " with type "
               << ifv.ift.to_string() << " with default value "
               << ifv.value << endl;
        }
        
        iff.args.push_back(ifv);
        
        next_line(fin,line,vs,done);
        if (done) function_done=true;
      }
      
      cout << "Function " << iff.name << " done." << endl;
      cout << endl;
      functions.push_back(iff);
    }

    if (false) {
      cout << "Done with section." << endl;
      char ch;
      cin >> ch;
    }
    
  }

  // Close file
  fin.close();

  // ----------------------------------------------------------------
  // Post-process to handle function overloading

  // For each class, manually examine each pair of member functions to
  // see if they have the same name
  for(size_t i=0;i<classes.size();i++) {
    if_class &ifc=classes[i];
    for(size_t j=0;j<ifc.methods.size();j++) {
      if_func &iff=ifc.methods[j];
      for(size_t k=j+1;k<ifc.methods.size();k++) {
        if_func &iff2=ifc.methods[k];
        if (iff.name==iff2.name) {
          // They must have different python names
          if (iff.py_name==iff2.py_name || iff.py_name.length()==0 ||
              iff2.py_name.length()==0) {
            cout << j << " " << k << " " << ifc.methods.size() << endl;
            cout << "In class: " << ifc.name << endl;
            cout << "Function names: '" << iff.name << "' and '"
                 << iff2.name << "'." << endl;
            cout << "Function python names: '"
                 << iff.py_name << "' and '" << iff2.py_name
                 << "'." << endl;
            O2SCL_ERR2("Member functions with same name and same ",
                       "py_name or with empty py_name.",
                       o2scl::exc_einval);
          }
          // Set the overloaded flag
          iff.overloaded=true;
          iff2.overloaded=true;
          cout << "In class " << ifc.name << ", member functions "
               << iff.name << " are overloaded and will use names\n  "
               << iff.py_name << " and " << iff2.py_name
               << "." << endl;
        }
      }
    }
  }

  // Manually examine each pair of global functions to see if they
  // have the same name
  for(size_t i=0;i<functions.size();i++) {
    if_func &iff=functions[i];
    for(size_t j=i+1;j<functions.size();j++) {
      if_func &iff2=functions[j];
      if (iff.name==iff2.name) {
        // They must have different python names
        if (iff.py_name==iff2.py_name || iff.py_name.length()==0 ||
            iff2.py_name.length()==0) {
          std::cout << "iff.name: " << iff.name << std::endl;
          std::cout << "iff2.name: " << iff2.name << std::endl;
          std::cout << "iff.py_name: " << iff.py_name << std::endl;
          std::cout << "iff2.py_name: " << iff2.py_name << std::endl;
          O2SCL_ERR2("Functions with same name and same py_name or ",
                     "empty py_name.",o2scl::exc_einval);
        }
        iff.overloaded=true;
        iff2.overloaded=true;
        cout << "Functions "
             << iff.name << " are overloaded and will use names\n  "
             << iff.py_name << " and " << iff2.py_name
             << "." << endl;
      }
    }
  }
  
  // End of interface file parsing code
  // ----------------------------------------------------------------
  // Create C++ header/source

  for(size_t kk=0;kk<2;kk++) {
    
    bool header=true, source=false;
    if (kk==1) {
      source=true;
      header=false;
    }
    
    ofstream fout;
    if (header) {
      fout.open((cpp_prefix+".h").c_str());
    } else {
      fout.open((cpp_prefix+".cpp").c_str());
    }
    
    fout << "/*" << endl;
    for(size_t i=0;i<header_strings.size();i++) {
      fout << header_strings[i] << endl;
    }
    fout << "*/" << endl;
    fout << endl;

    if (header) {
      
      for(size_t i=0;i<h_includes.size();i++) {
        fout << "#include " << h_includes[i] << endl;
      }
      fout << endl;
      
      fout << "extern \"C\" {" << endl;
      fout << endl;
      
    } else {
      
      for(size_t i=0;i<cpp_includes.size();i++) {
        fout << "#include " << cpp_includes[i] << endl;
      }
      fout << endl;
      
      for(size_t i=0;i<cpp_using.size();i++) {
        fout << "using namespace " << cpp_using[i] << ";" << endl;
      }
      fout << endl;
      
    }

    for(size_t i=0;i<classes.size();i++) {

      if_class &ifc=classes[i];

      // Generate an object creation function as long as the class is
      // not abstract and has a default constructor
      
      if (ifc.is_abstract==false && ifc.def_cons==true) {
        fout << "void *" << underscoreify(ifc.ns) << "_create_"
             << underscoreify(ifc.name) << "()";
        if (header) {
          fout << ";" << endl;
        } else {
          fout << " {" << endl;
          fout << "  " << ifc.name << " *ptr=new " << ifc.name
               << ";" << endl;
          fout << "  return ptr;" << endl;
          fout << "}" << endl;
        }
        fout << endl;
      }
      
      // Generate a destructor function as long as the class is not
      // abstract
      
      if (ifc.is_abstract==false) {
        fout << "void " << underscoreify(ifc.ns) << "_free_"
             << underscoreify(ifc.name) << "(void *vptr)";
        if (header) {
          fout << ";" << endl;
        } else {
          fout << " {" << endl;
          fout << "  " << ifc.name << " *ptr=(" << ifc.name
               << " *)vptr;" << endl;
          fout << "  delete ptr;" << endl;
          fout << "  return;" << endl;
          fout << "}" << endl;
        }
        fout << endl;
      }

      // If the standard copy constructors are included, create a
      // _copy_ method on the C side which we can use to create a
      // deep_copy method in python. This function presumes that the
      // source and destination pointers already point to a valid
      // object.
      
      if (ifc.std_cc) {
        fout << "void " << underscoreify(ifc.ns) << "_copy_"
             << underscoreify(ifc.name) << "(void *vsrc, void *vdest)";
        if (header) {
          fout << ";" << endl;
        } else {
          fout << " {" << endl;
          fout << "  " << ifc.name << " *src=(" << ifc.name
               << " *)vsrc;" << endl;
          fout << "  " << ifc.name << " *dest=(" << ifc.name
               << " *)vdest;" << endl;
          fout << "  *dest=*src;" << endl;
          fout << "}" << endl;
        }
        fout << endl;
      }

      for(size_t j=0;j<ifc.members.size();j++) {

        if_var &ifv=ifc.members[j];

        // --------------------------------------------------------------
        // Get functions for class data members

        if (ifv.ift.is_ctype()) {
          
          // Get function for a C data type
          if (ifv.ift.is_reference()) {
            fout << "*";
          }
          fout << ifv.ift.name << " " << underscoreify(ifc.ns) << "_"
               << underscoreify(ifc.name) << "_get_" << ifv.name
               << "(void *vptr)";
          if (header) {
            fout << ";" << endl;
          } else {
            fout << " {" << endl;
            fout << "  " << ifc.name << " *ptr=(" << ifc.name
                 << " *)vptr;" << endl;
            if (ifv.ift.is_reference()) {
              fout << "  return &ptr->" << ifv.name << ";" << endl;
            } else {
              fout << "  return ptr->" << ifv.name << ";" << endl;
            }
            fout << "}" << endl;
          }
          
        } else if (ifv.ift.name=="std::string" ||
                   ifv.ift.name=="string") {
          
          // Get function for string data
          fout << "void *" << underscoreify(ifc.ns) << "_"
               << underscoreify(ifc.name) << "_get_" << ifv.name
               << "(void *vptr)";
          if (header) {
            fout << ";" << endl;
          } else {
            fout << " {" << endl;
            fout << "  " << ifc.name << " *ptr=(" << ifc.name
                 << " *)vptr;" << endl;
            fout << "  std::string *sptr=new std::string;" << endl;
            fout << "  *sptr=ptr->" << ifv.name << ";" << endl;
            fout << "  return sptr;" << endl;
            fout << "}" << endl;
          }
          
        } else if (ifv.ift.is_shared_ptr()) {
          
          // Get function for shared pointers
          fout << "void " << underscoreify(ifc.ns) << "_"
               << underscoreify(ifc.name) << "_get_" << ifv.name
               << "(void *vptr, void *p_v)";
          if (header) {
            fout << ";" << endl;
          } else {
            fout << " {" << endl;            
            fout << "  " << ifc.name << " *ptr=(" << ifc.name
                 << " *)vptr;" << endl;
            fout << "  std::shared_ptr<" << ifv.ift.name
                 << " > *p_tgsp=(std::shared_ptr<"
                 << ifv.ift.name << " > *)p_v;" << endl;
            fout << "  *(p_tgsp)=ptr->" << ifv.name << ";" << endl;
            fout << "  return;" << endl;
            fout << "}" << endl;
          }
          
        } else {
          
          // Get function for other types
          fout << "void *" << underscoreify(ifc.ns) << "_"
               << underscoreify(ifc.name) << "_get_" << ifv.name
               << "(void *vptr)";
          if (header) {
            fout << ";" << endl;
          } else {
            fout << " {" << endl;            
            fout << "  " << ifc.name << " *ptr=(" << ifc.name
                 << " *)vptr;" << endl;
            //fout << "  " << ifv.ift.name << " *p_tgot=("
            //<< ifv.ift.name << " *)p_v;" << endl;
            //fout << "  *(p_tgot)=ptr->" << ifv.name << ";" << endl;
            if (ifv.ift.is_pointer()) {
              fout << "  return (void *)((ptr->" << ifv.name << "));" << endl;
            } else {
              fout << "  return (void *)(&(ptr->" << ifv.name << "));" << endl;
            }
            fout << "}" << endl;
          }
          
        }
        fout << endl;
      
        // Set functions for class data
        if (ifv.ift.is_ctype()) {
          if (!ifv.ift.is_const() &&
              !ifv.ift.is_reference()) {
            // Set function for a C data type
            fout << "void " << underscoreify(ifc.ns) << "_"
                 << underscoreify(ifc.name) << "_set_"
                 << ifv.name << "(void *vptr, "
                 << ifv.ift.name << " v)";
            if (header) {
              fout << ";" << endl;
            } else {
              fout << " {" << endl;
              fout << "  " << ifc.name << " *ptr=(" << ifc.name
                   << " *)vptr;" << endl;
              fout << "  ptr->" << ifv.name << "=v;" << endl;
              fout << "  return;" << endl;
              fout << "}" << endl;
            }
          }
        } else {
          // Set function for other types
          fout << "void " << underscoreify(ifc.ns) << "_"
               << underscoreify(ifc.name) << "_set_" << ifv.name
               << "(void *vptr, void *p_v)";
          if (header) {
            fout << ";" << endl;
          } else {
            fout << " {" << endl;
            fout << "  " << ifc.name << " *ptr=(" << ifc.name
                 << " *)vptr;" << endl;
            if (ifv.ift.is_shared_ptr()) {
              // Shared pointers
              fout << "  std::shared_ptr<" << ifv.ift.name
                   << " > *p_tssp=(std::shared_ptr<"
                   << ifv.ift.name << " > *)p_v;" << endl;
              fout << "  ptr->" << ifv.name << "=*(p_tssp);" << endl;
              fout << "  return;" << endl;
              fout << "}" << endl;
            } else if (ifv.ift.is_pointer()) {
              // Other types
              fout << "  " << ifv.ift.name << " *p_tsptr=("
                   << ifv.ift.name << " *)p_v;" << endl;
              fout << "  ptr->" << ifv.name << "=p_tsptr;" << endl;
              fout << "  return;" << endl;
              fout << "}" << endl;
            } else {
              // Other types
              fout << "  " << ifv.ift.name << " *p_tsot=("
                   << ifv.ift.name << " *)p_v;" << endl;
              fout << "  ptr->" << ifv.name << "=*(p_tsot);" << endl;
              fout << "  return;" << endl;
              fout << "}" << endl;
            }
          }
        }
        fout << endl;

        // End of loop over class data members
      }
    
      for(size_t j=0;j<ifc.methods.size();j++) {
      
        if_func &iff=ifc.methods[j];

        // The C code for the return type 
        string ret_type;
        // Extra arguments for the current method
        string extra_args;
        
        // Function header, first determine return type and
        // any extra arguments
        if (iff.ret.name=="std::string") {
          ret_type="void *";
        } else if ((iff.ret.name=="vector<double>" ||
                    iff.ret.name=="std::vector<double>") &&
                   iff.ret.suffix=="&") {
          // For a std::vector<double> &, we return a void, but add
          // double pointer and integer output parameters
          ret_type="void ";
          extra_args=", double **dptr, int *n_";
        } else if (iff.ret.name=="void" || iff.ret.is_ctype()) {
          if (iff.ret.suffix=="*") {
            ret_type=iff.ret.name+" *";
          } else {
            ret_type=iff.ret.name+" ";
          }
          if (iff.ret.prefix.length()>0) {
            ret_type=iff.ret.prefix+" "+ret_type;
          }
        } else {
          ret_type="void *";
        }

        // Generate the code for the function declaration part
        if (iff.overloaded) {
          fout << ret_type << underscoreify(ifc.ns) << "_"
               << underscoreify(ifc.name) << "_" << iff.py_name
               << "(void *vptr";
        } else if (iff.name=="operator[]" || iff.name=="operator()") {
          fout << ret_type << underscoreify(ifc.ns) << "_"
               << underscoreify(ifc.name) << "_getitem"
               << "(void *vptr";
        } else {
          fout << ret_type << underscoreify(ifc.ns) << "_"
               << underscoreify(ifc.name) << "_" << iff.name
               << "(void *vptr";
        }
        if (iff.args.size()>0) {
          fout << ", ";
        }
        for(size_t k=0;k<iff.args.size();k++) {
          if (iff.args[k].ift.suffix=="") {
            if (iff.args[k].ift.name=="std::string") {
              fout << "char *" << iff.args[k].name;
            } else {
              fout << iff.args[k].ift.name << " " << iff.args[k].name;
            }
          } else if (iff.args[k].ift.suffix=="&") {
            if (iff.args[k].ift.is_ctype()) {
              fout << iff.args[k].ift.name << " *"
                   << iff.args[k].name;
            } else {
              fout << "void *ptr_" << iff.args[k].name;
            }
          }
          // Output default value if we're in the header file
          if (header) {
            // String arguments are converted to char *'s, so
            // they don't need a default value
            if (iff.args[k].ift.name!="std::string") {
              if (iff.args[k].value.length()>0) {
                if (iff.args[k].value=="True") {
                  fout << "=true";
                } else if (iff.args[k].value=="False") {
                  fout << "=false";
                } else {
                  fout << "=" << iff.args[k].value;
                }
              }
            }
          }
          if (k!=iff.args.size()-1) {
            fout << ", ";
          }
        }
        
        if (header) {
          
          fout << extra_args << ");" << endl;
          
        } else {
          
          fout << extra_args << ") {" << endl;
          
          // Pointer assignment for class
          fout << "  " << ifc.name << " *ptr=("
               << ifc.name << " *)vptr;" << endl;
          
          // If the argument is a reference and not a standard C type,
          // then we'll need to convert from a void *
          for(size_t k=0;k<iff.args.size();k++) {
            if (iff.args[k].ift.suffix=="&" && !iff.args[k].ift.is_ctype()) {
              std::string type_temp=iff.args[k].ift.name;
              if (type_temp=="std_vector") {
                type_temp="std::vector<double>";
              }
              fout << "  " << type_temp << " *"
                   << iff.args[k].name << "=("
                   << type_temp << " *)ptr_" << iff.args[k].name
                   << ";" << endl;
            }
          }
          
          // Now generate code to call the member function from
          // inside the wrapper function. Also generate the code
          // for the return
          
          if (iff.ret.prefix.find("shared_ptr")!=std::string::npos ||
              iff.ret.prefix.find("std::shared_ptr")!=std::string::npos) {
            // If we're returning a shared pointer, then we actually
            // need to create a pointer to the shared pointer first
            fout << "  std::shared_ptr<" << iff.ret.name << " > *ret=new"
                 << " std::shared_ptr<" << iff.ret.name << " >;" << endl;
            // Then we set the shared pointer using the class
            // member function
            fout << "  *ret=ptr->" << iff.name << "(";
          } else if ((iff.ret.name=="vector<double>" ||
                      iff.ret.name=="std::vector<double>") &&
                     iff.ret.suffix=="&") {
            // In case it's const, we have to explicitly typecast
            fout << "  const std::vector<double> &r=ptr->" << iff.name << "(";
          } else if (iff.name=="operator[]" || iff.name=="operator()") {
            if (iff.ret.name=="std::string") {
              fout << "  std::string *sptr=new std::string;" << endl;
              fout << "  *sptr=ptr->" << iff.name << "(";
            } else if (iff.ret.name=="std::vector<std::string>") {
              fout << "  std::vector<std::string> *vsptr="
                   << "new std::vector<std::string>;" << endl;
              fout << "  *vsptr=ptr->" << iff.name << "(";
            } else if (iff.ret.name=="contour_line") {
              fout << "  " << iff.ret.name
                   << " *ret=&(ptr->" << iff.name << "(";
            } else if (iff.ret.name=="void") {
              fout << "  ptr->" << iff.name << "(";
            } else {
              fout << "  " << iff.ret.name
                   << " ret=ptr->" << iff.name << "(";
            }
          } else if (iff.ret.name=="void") {
            fout << "  ptr->" << iff.name << "(";
          } else if (iff.ret.name=="std::string") {
            fout << "  std::string *sptr=new std::string;" << endl;
            fout << "  *sptr=ptr->" << iff.name << "(";
          } else {
            if (iff.ret.suffix=="&") {
              // If the function returns a reference, return a pointer
              // instead
              if (iff.ret.prefix.length()>0) {
                fout << "  " << iff.ret.prefix << " ";
              } else {
                fout << "  ";
              }
              fout << iff.ret.name << " *ret=&ptr->" << iff.name << "(";
            } else if (iff.ret.suffix=="*") {
              // If it's a pointer, then we can just return the pointer
              if (iff.ret.prefix.length()>0) {
                fout << "  " << iff.ret.prefix << " ";
              } else {
                fout << "  ";
              }
              fout << iff.ret.name << " *ret=ptr->" << iff.name << "(";
            } else if (iff.ret.is_ctype()) {
              // If it's a C type, then we can just return by value
              if (iff.ret.prefix.length()>0) {
                fout << "  " << iff.ret.prefix << " ";
              } else {
                fout << "  ";
              }
              fout << iff.ret.name << " ret=ptr->" << iff.name << "(";
            } else {
              // If it's a class, and we're returning by value,
              // then we have to return a pointer to a new object
              fout << "  " << iff.ret.name << " *ret=new "
                   << iff.ret.name << ";" << endl;
              if (iff.ret.prefix.length()>0) {
                fout << "  " << iff.ret.prefix << " ";
              } else {
                fout << "  ";
              }
              fout << "*ret=ptr->" << iff.name << "(";
            }
          }
          
          for(size_t k=0;k<iff.args.size();k++) {
            if (iff.args[k].ift.suffix=="") {
              fout << iff.args[k].name;
            } else if (iff.args[k].ift.suffix=="&") {
              fout << "*" << iff.args[k].name;
            }
            if (k!=iff.args.size()-1) {
              fout << ",";
            }
          }
          if (iff.name=="operator[]" && iff.ret.name=="contour_line") {
            fout << "));" << endl;
          } else {
            fout << ");" << endl;
          }

          // Construct the function return statement
          
          if ((iff.ret.name=="vector<double>" ||
               iff.ret.name=="std::vector<double>") &&
              iff.ret.suffix=="&") {
            
            // Extra code for array types
            fout << "  *dptr=(double *)(&(r[0]));" << endl;
            fout << "  *n_=r.size();" << endl;
            
            fout << "  return;" << endl;
            
          } else {
            
            if (iff.ret.name=="std::string") {
              fout << "  return sptr;" << endl;
            } else if (iff.ret.name=="std::vector<std::string>") {
              fout << "  return vsptr;" << endl;
            } else if (iff.ret.name=="void") {
              fout << "  return;" << endl;
            } else if (iff.ret.is_ctype() || iff.ret.is_reference() ||
                       iff.ret.is_shared_ptr()) {
              if ((iff.ret.name=="std::vector<double>" ||
                   iff.ret.name=="std::vector<int>" ||
                   iff.ret.name=="boost::numeric::ublas::vector<double>" ||
                   iff.ret.name=="std::vector<size_t>") &&
                  iff.ret.is_const()) {
                // Cast away const for conversions between const pointers
                // and void *
                fout << "  return (void *)ret;" << endl;
              } else {
                fout << "  return ret;" << endl;
              }
            } else if (iff.ret.is_ctype()) {
              fout << "  return &ret;" << endl;
            } else {
              // If it's not a reference or a pointer, then if it's
              // not a ctype we need to return a pointer to a new
              // object.
              fout << "  return ret;" << endl;
            }
          }
          
          // Ending function brace
          fout << "}" << endl;          
        }
        fout << endl;

        // Generate setitem code for operator[] if it returns a
        // non-const reference. Presume ifc.name is something like
        // std_vector and iff.name is "operator[]" and iff.ret.name is
        // "double"
        if (iff.name=="operator[]" && !iff.ret.is_const() &&
            iff.ret.suffix=="&") {
          if (iff.ret.name=="std::string") {
            fout << "void " << ifc.ns << "_" << underscoreify(ifc.name)
                 << "_setitem(void *vptr, size_t i, std::string *val)";
          } else if (iff.ret.name=="std::vector<double>") {
            fout << "void " << ifc.ns << "_" << underscoreify(ifc.name)
                 << "_setitem(void *vptr, size_t i, "
                 << "void *valptr)";
          } else if (iff.ret.name=="contour_line") {
            fout << "void " << ifc.ns << "_" << underscoreify(ifc.name)
                 << "_setitem(void *vptr, size_t i, "
                 << "void *valptr)";
          } else {
            fout << "void " << ifc.ns << "_" << underscoreify(ifc.name)
                 << "_setitem(void *vptr, size_t i, " << iff.ret.name
                 << " val)";
          }
          if (header) {
            fout << ";" << endl;
          } else {
            fout << " {" << endl;
            fout << "  " << ifc.name << " *ptr=(" << ifc.name
                 << " *)vptr;" << endl;
            if (iff.ret.name=="std::vector<double>") {
              fout << "  std::vector<double> *valptr2="
                   << "(std::vector<double> *)valptr;" << endl;
              fout << "  (*ptr)[i]=*valptr2;" << endl;
            } else if (iff.ret.name=="contour_line") {
              fout << "  contour_line *valptr2="
                   << "(contour_line *)valptr;" << endl;
              fout << "  (*ptr)[i]=*valptr2;" << endl;
            } else if (iff.ret.name=="std::string") {
              fout << "  (*ptr)[i]=*val;" << endl;
            } else {
              fout << "  (*ptr)[i]=val;" << endl;
            }
            fout << "  return;" << endl;
            fout << "}" << endl;
          }
          fout << endl;
        }
        
        // Generate setitem code for operator() if it returns a
        // non-const reference. Presume ifc.name is something like
        // std_vector and iff.name is "operator()" and iff.ret.name is
        // "double"
        if (iff.name=="operator()" && !iff.ret.is_const() &&
            iff.ret.suffix=="&") {
          fout << "void " << ifc.ns << "_" << underscoreify(ifc.name)
               << "_setitem(void *vptr, size_t i, size_t j, " << iff.ret.name
               << " val)";
          if (header) {
            fout << ";" << endl;
          } else {
            fout << " {" << endl;
            fout << "  " << ifc.name << " *ptr=(" << ifc.name
                 << " *)vptr;" << endl;
            fout << "  (*ptr)(i,j)=val;" << endl;
            fout << "  return;" << endl;
            fout << "}" << endl;
          }
          fout << endl;
        }
        
      }

      for(size_t j=0;j<ifc.cons.size();j++) {
      
        if_func &iff=ifc.cons[j];

        // Now generate the actual code for the constructor

        fout << "void *" << underscoreify(ifc.ns) << "_"
             << underscoreify(ifc.name) << "_" << iff.name
             << "(";
        for(size_t k=0;k<iff.args.size();k++) {
          if (iff.args[k].ift.suffix=="") {
            if (iff.args[k].ift.name=="std::string") {
              fout << "char *" << iff.args[k].name;
            } else {
              fout << iff.args[k].ift.name << " " << iff.args[k].name;
            }
          } else if (iff.args[k].ift.suffix=="&") {
            if (iff.args[k].ift.name=="std::string") {
              fout << "void *ptr_" << iff.args[k].name;
            } else if (iff.args[k].ift.name=="vector<size_t>" ||
                       iff.args[k].ift.name=="std::vector<size_t>") {
              fout << "void *ptr_" << iff.args[k].name;
            } else {
              cout << "class: " << ifc.name << " constructor name: "
                   << iff.name << endl;
              O2SCL_ERR2("Cannot yet handle this kind of reference in a ",
                         "constructor.",o2scl::exc_eunimpl);
            }
          }
          if (k!=iff.args.size()-1) {
            fout << ", ";
          }
        }
        
        if (header) {
          fout << ");" << endl;
          
        } else {
          
          fout << ") {" << endl;
          
          // If the argument is a reference and not a standard C type,
          // then we'll need to convert from a void *
          for(size_t k=0;k<iff.args.size();k++) {
            if (iff.args[k].ift.suffix=="&" && !iff.args[k].ift.is_ctype()) {
              std::string type_temp=iff.args[k].ift.name;
              if (type_temp=="std_vector") {
                type_temp="std::vector<double>";
              }
              fout << "  " << type_temp << " *"
                   << iff.args[k].name << "=("
                   << type_temp << " *)ptr_" << iff.args[k].name
                   << ";" << endl;
            }
          }
          
          // Pointer assignment for class
          fout << "  " << ifc.name << " *ptr=new "
               << ifc.name << "(";
          
          for(size_t k=0;k<iff.args.size();k++) {
            if (iff.args[k].ift.suffix=="") {
              fout << iff.args[k].name;
            } else if (iff.args[k].ift.suffix=="&") {
              fout << "*" << iff.args[k].name;
            }
            if (k!=iff.args.size()-1) {
              fout << ",";
            }
          }
          
          fout << ");" << endl;

          fout << "  return ptr;" << endl;
          
          // Ending function brace
          fout << "}" << endl;          
        }
        fout << endl;
        
      }

    }
  
    for(size_t i=0;i<sps.size();i++) {
    
      if_shared_ptr &ifsp=sps[i];

      // Generate code for the create pointer function
      // Always create a void *
      fout << "void *" << underscoreify(ifsp.ns) << "_create_shared_ptr_"
           << underscoreify(ifsp.name) << "()";
      if (header) {
        fout << ";" << endl;
      } else {
        fout << " {" << endl;
        fout << "  std::shared_ptr<" << ifsp.name
             << " > *ptr=new std::shared_ptr<" << ifsp.name
             << " >(new " << ifsp.name << ");" << endl;
        fout << "  return ptr;" << endl;
        fout << "}" << endl;
      }
      fout << endl;
      
      // Generate code for the free pointer function
      // Always free a void *
      fout << "void " << underscoreify(ifsp.ns) << "_free_shared_ptr_"
           << underscoreify(ifsp.name) << "(void *vptr)";
      if (header) {
        fout << ";" << endl;
      } else {
        fout << " {" << endl;
        fout << "  std::shared_ptr<" << ifsp.name
             << " > *ptr=(std::shared_ptr<" << ifsp.name
             << " > *)vptr;" << endl;
        fout << "  delete ptr;" << endl;
        fout << "}" << endl;
      }
      fout << endl;
      
      // Function to allow python to create a raw pointer from the
      // shared pointer
      fout << "void *" << underscoreify(ifsp.ns) << "_shared_ptr_"
           << underscoreify(ifsp.name) << "_ptr(void *vp)";
      if (header) {
        fout << ";" << endl;
      } else {
        fout << " {" << endl;
        fout << "  std::shared_ptr<" << ifsp.name
             << " > *p=(std::shared_ptr<" << ifsp.name << " > *)vp;" << endl;
        fout << "  " << ifsp.name << " *ref=p->get();" << endl;
        fout << "  return ref;" << endl;
        fout << "}" << endl;        
      }
      fout << endl;
      
    }    
  
    for(size_t i=0;i<functions.size();i++) {
    
      if_func &iff=functions[i];
    
      // Function header
      
      string func_name=iff.name;
      if (iff.overloaded) {
        func_name=iff.py_name;
      } else {
        func_name=underscoreify(func_name);
      }
      
      if (iff.ret.name=="std::string") {
        fout << "void *" << underscoreify(iff.ns) << "_"
             << func_name << "(";
      } else if (iff.ret.name=="void" || iff.ret.is_ctype()) {
        fout << iff.ret.name << " " << underscoreify(iff.ns) << "_"
             << func_name << "_wrapper(";
      } else {
        fout << "void *" << underscoreify(iff.ns) << "_"
             << func_name << "_wrapper(";
      }
    
      for(size_t k=0;k<iff.args.size();k++) {
        
        if (iff.args[k].ift.suffix=="") {
          if (iff.args[k].ift.name=="std::string") {
            fout << "char *" << iff.args[k].name;
          } else {
            fout << iff.args[k].ift.name << " " << iff.args[k].name;
          }
        } else if (iff.args[k].ift.suffix=="&") {
          if (iff.args[k].ift.name=="std::string") {
            fout << "void *&ptr_" << iff.args[k].name;
          } else {
            fout << "void *ptr_" << iff.args[k].name;
          }
        }
        
        // Output default value
        if (header) {
          if (iff.args[k].value.length()>0) {
            if (iff.args[k].value=="True") {
              fout << "=true";
            } else if (iff.args[k].value=="False") {
              fout << "=false";
            } else {
              fout << "=" << iff.args[k].value;
            }
          }
        }

        // Comma before next function parameter
        if (k!=iff.args.size()-1) {
          fout << ", ";
        }
      }
    
      fout << ")";
      
      if (header) {
        fout << ";" << endl;
      } else {
        fout << " {" << endl;

        // Pointer assignments for arguments
        for(size_t k=0;k<iff.args.size();k++) {
          if (iff.args[k].ift.suffix=="&") {
            if (iff.args[k].ift.name=="std::string") {
              fout << "  std::string *"
                   << iff.args[k].name
                   << "=new std::string;" << endl;
            } else {
              fout << "  " << iff.args[k].ift.name << " *"
                   << iff.args[k].name << "=("
                   << iff.args[k].ift.name << " *)ptr_" << iff.args[k].name
                   << ";" << endl;
            }
          }
        }
        
        // Now generate code for actual function call and the
        // return statement
        if (iff.ret.name=="void") {
          
          fout << "  " << iff.name << "(";
          
          vector<std::string> addl_code;

          for(size_t k=0;k<iff.args.size();k++) {
            if (iff.args[k].ift.suffix=="") {
              fout << iff.args[k].name;
            } else if (iff.args[k].ift.suffix=="&") {
              fout << "*" << iff.args[k].name;
              if (iff.args[k].ift.name=="std::string") {
                addl_code.push_back(((string)"ptr_")+iff.args[k].name+
                                    "=(void *)"+iff.args[k].name+";");
              }
            }
            if (k!=iff.args.size()-1) {
              fout << ",";
            }
          }
          fout << ");" << endl;
          
          for(size_t k=0;k<addl_code.size();k++) {
            fout << "  " << addl_code[k] << endl;
          }

          fout << "  return;" << endl;
          
        } else {
          
          if (iff.ret.is_ctype() || iff.ret.is_reference()) {
            fout << "  " << iff.ret.name << " ret=" << iff.name << "(";
          } else {
            fout << "  " << iff.ret.name << " *ret=new "
                 << iff.ret.name << ";" << endl;
            fout << "  *ret=" << iff.name << "(";
          }

          for(size_t k=0;k<iff.args.size();k++) {
            if (iff.args[k].ift.suffix=="") {
              fout << iff.args[k].name;
            } else if (iff.args[k].ift.suffix=="&") {
              fout << "*" << iff.args[k].name;
            }
            if (k!=iff.args.size()-1) {
              fout << ",";
            }
          }
          fout << ");" << endl;
          
          fout << "  return ret;" << endl;
        }
        fout << "}" << endl;
        
      }
      fout << endl;
      
    }

    if (header) {
      // End of the extern C section
      fout << "}" << endl;
    }
  
    fout.close();

  }
  
  // ----------------------------------------------------------------
  // Create python source code

  ofstream fout;
  fout.open((py_prefix+".py").c_str());

  fout << "\"\"\"" << endl;
  for(size_t i=0;i<header_strings.size();i++) {
    fout << header_strings[i] << endl;
  }
  fout << "\"\"\"" << endl;
  fout << endl;
  
  fout << "import ctypes" << endl;
  fout << "from abc import abstractmethod" << endl;
  fout << "from o2sclpy.utils import force_bytes" << endl;
  if (import_numpy) {
    fout << "import numpy" << endl;
  }
  fout << endl;

  if (py_headers.size()>0) {
    for(size_t i=0;i<py_headers.size();i++) {
      fout << py_headers[i] << endl;
    }
    fout << endl;
  }
  
  for(size_t i=0;i<classes.size();i++) {
    
    if_class &ifc=classes[i];

    // The first line of the class definition, including the class
    // and the parent class
    
    if (ifc.parents.size()==0) {
      if (ifc.py_name!="") {
        fout << "class " << ifc.py_name << ":" << endl;
      } else {
        fout << "class " << ifc.name << ":" << endl;
      }
    } else {
      // Go through all the classes to see if the
      // parent class also has a py_name
      string parent_py_name=ifc.parents[0];
      for(size_t j=0;j<classes.size();j++) {
        if_class &jfc=classes[j];
        if (jfc.name==ifc.parents[0] && jfc.py_name!="") {
          parent_py_name=jfc.py_name;
        }
      }
      if (ifc.py_name!="") {
        fout << "class " << ifc.py_name << "(" << parent_py_name << "):"
             << endl;
      } else {
        fout << "class " << ifc.name << "(" << parent_py_name << "):"
             << endl;
      }
    }

    // The class documentation
    fout << "    \"\"\"" << endl;
    if (ifc.py_class_doc.size()>0) {
      for(size_t j=0;j<ifc.py_class_doc.size();j++) {
        string s=ifc.py_class_doc[j];
        while (s.find("%name%")!=std::string::npos) {
          s.replace(s.find("%name%"),6,ifc.name);
        }
        fout << "    " << s << endl;
      }
    } else if (py_class_doc.size()>0) {
      for(size_t j=0;j<py_class_doc.size();j++) {
        string s=py_class_doc[j];
        while (s.find("%name%")!=std::string::npos) {
          s.replace(s.find("%name%"),6,ifc.name);
        }
        fout << "    " << s << endl;
      }
    }
    fout << "    \"\"\"" << endl;
    fout << endl;
    
    // List and initialize class data members
    if (ifc.parents.size()==0) {
      fout << "    _ptr=0" << endl;
      fout << "    _link=0" << endl;
      fout << "    _owner=True" << endl;
      fout << endl;
    }
    
    // Define the __init__() function
    if (ifc.is_abstract || ifc.def_cons==false) {
      fout << "    @abstractmethod" << endl;
    }
    fout << "    def __init__(self,link,pointer=0):" << endl;
    fout << "        \"\"\"" << endl;
    fout << "        Init function for class ";
    if (ifc.py_name!="") {
      fout << ifc.py_name << endl;
    } else {
      fout << ifc.name << endl;
    }
    fout << endl;
    fout << "        | Parameters:" << endl;
    fout << "        | *link* :class:`linker` object" << endl;
    fout << "        | *pointer* ``ctypes.c_void_p`` pointer" << endl;
    fout << endl;
    fout << "        \"\"\"" << endl;
    fout << endl;
    fout << "        if pointer==0:" << endl;
    fout << "            f=link." << dll_name << "." << ifc.ns
         << "_create_" << underscoreify(ifc.name) << endl;
    fout << "            f.restype=ctypes.c_void_p" << endl;
    fout << "            f.argtypes=[]" << endl;
    fout << "            self._ptr=f()" << endl;
    fout << "        else:" << endl;
    fout << "            self._ptr=pointer" << endl;
    fout << "            self._owner=False" << endl;
    fout << "        self._link=link" << endl;
    fout << "        return" << endl;
    fout << endl;
    
    // Define the __del__() function
    fout << "    def __del__(self):" << endl;
    fout << "        \"\"\"" << endl;
    fout << "        Delete function for class ";
    if (ifc.py_name!="") {
      fout << ifc.py_name << endl;
    } else {
      fout << ifc.name << endl;
    }
    fout << "        \"\"\"" << endl;
    fout << endl;
    fout << "        if self._owner==True:" << endl;
    fout << "            f=self._link." << dll_name << "." << ifc.ns
         << "_free_" << underscoreify(ifc.name) << endl;
    fout << "            f.argtypes=[ctypes.c_void_p]" << endl;
    fout << "            f(self._ptr)" << endl;
    fout << "            self._owner=False" << endl;
    fout << "            self._ptr=0" << endl;
    fout << "        return" << endl;
    fout << endl;

    // Define the copy() function
    fout << "    def __copy__(self):" << endl;
    fout << "        \"\"\"" << endl;
    fout << "        Shallow copy function for class ";
    if (ifc.py_name!="") {
      fout << ifc.py_name << endl;
    } else {
      fout << ifc.name << endl;
    }
    fout << "        " << endl;
    fout << "        Returns: ";
    if (ifc.py_name!="") {
      fout << ifc.py_name << " object" <<  endl;
    } else {
      fout << ifc.name << " object" << endl;
    }
    
    fout << "        \"\"\"" << endl;
    fout << endl;
    fout << "        new_obj=type(self)(self._link,self._ptr)" << endl;
    fout << "        return new_obj" << endl;
    fout << endl;

    // Define the deepcopy() function if it has a standard
    // copy constructor
    
    if (ifc.std_cc) {
      fout << "    def __deepcopy__(self,memo):" << endl;
      fout << "        \"\"\"" << endl;
      fout << "        Deep copy function for class ";
      if (ifc.py_name!="") {
        fout << ifc.py_name << endl;
      } else {
        fout << ifc.name << endl;
      }
      fout << "        " << endl;
      fout << "        Returns: new copy of the ";
      if (ifc.py_name!="") {
        fout << ifc.py_name << " object" <<  endl;
      } else {
        fout << ifc.name << " object" << endl;
      }
      fout << "        \"\"\"" << endl;
      fout << endl;
      // Create the new object
      fout << "        new_obj=type(self)(self._link)" << endl;
      fout << "        f2=self._link." << dll_name << "." << ifc.ns
           << "_copy_" << underscoreify(ifc.name) << endl;
      fout << "        f2.argtypes=[ctypes.c_void_p,ctypes.c_void_p]"
           << endl;
      fout << "        f2(self._ptr,new_obj._ptr)" << endl;
      fout << "        return new_obj" << endl;
      fout << endl;
    }

    // Define member get and set properties
    for(size_t j=0;j<ifc.members.size();j++) {
      
      if_var &ifv=ifc.members[j];

      // Hack, because del has a special meaning in python
      if (ifv.name=="del") ifv.name="delta";

      // Set up get functions for member data

      if (ifv.ift.is_ctype()) {

        // For raw C-types, set up a property and return a value
        fout << "    @property" << endl;
        fout << "    def " << ifv.name << "(self):" << endl;
        fout << "        \"\"\"" << endl;
        fout << "        Property of type ``ctypes.c_" << ifv.ift.name
             << "``" << endl;
        fout << "        \"\"\"" << endl;
        fout << "        func=self._link." << dll_name << "."
             << ifc.ns << "_" << underscoreify(ifc.name)
             << "_get_" << ifv.name << endl;
        fout << "        func.restype=ctypes.c_" << ifv.ift.name << endl;
        fout << "        func.argtypes=[ctypes.c_void_p]" << endl;
        fout << "        return func(self._ptr)" << endl;
        
      } else if (ifv.ift.prefix.find("shared_ptr")!=std::string::npos ||
                 ifv.ift.prefix.find("std::shared_ptr")!=std::string::npos) {

        // For shared pointers, return a new shared pointer object
        
        fout << "    def get_" << ifv.name << "(self," << ifv.name
             << "):" << endl;
        fout << "        \"\"\"" << endl;
        fout << "        Object of type :class:`"
             << ifv.ift.name << "`" << endl;
        fout << "        \"\"\"" << endl;

        // Create a new shared pointer object
        size_t len=ifv.ift.name.length();
        std::string tmps=ifv.ift.name;
        // Manually remove '<>' from the typename if necessary
        if (len>2 && ifv.ift.name[len-2]=='<' &&
            ifv.ift.name[len-1]=='>') {
          tmps=ifv.ift.name.substr(0,len-2);
        }
        fout << "        sp=shared_ptr_"+tmps+
          "(self._link)" << endl;
        
        fout << "        func=self._link." << dll_name << "." << ifc.ns << "_"
             << underscoreify(ifc.name)
             << "_get_" << ifv.name << endl;
        fout << "        func.argtypes=[ctypes.c_void_p," 
             << "ctypes.c_void_p]" << endl;
        fout << "        func(self._ptr,sp._s_ptr)" << endl;
        fout << "        return" << endl;

      } else if (ifv.ift.name=="std::string") {

        // Get strings by converting to a bytes object using the
        // std_string class

        fout << "    def get_" << ifv.name << "(self):" << endl;
        fout << "        \"\"\"" << endl;
        fout << "        Get object of type :class:`"
             << ifv.ift.name << "`" << endl;
        fout << "        \"\"\"" << endl;
        fout << "        func=self._link." << dll_name << "." << ifc.ns << "_"
             << underscoreify(ifc.name)
             << "_get_" << ifv.name << endl;

        fout << "        func.restype=ctypes.c_void_p" << endl;
        
        fout << "        func.argtypes=[ctypes.c_void_p]" << endl;
        fout << "        s=std_string(self._link)" << endl;
        fout << "        s._ptr=func(self._ptr)" << endl;
        fout << "        return s.to_bytes()" << endl;
        
      } else {

        // Get a reference to other types
        
        fout << "    def get_" << ifv.name << "(self):" << endl;
        fout << "        \"\"\"" << endl;
        fout << "        Get object of type :class:`"
             << ifv.ift.name << "`" << endl;
        fout << "        \"\"\"" << endl;
        fout << "        func1=self._link." << dll_name << "." << ifc.ns << "_"
             << underscoreify(ifc.name)
             << "_get_" << ifv.name << endl;
        
        fout << "        func1.restype=ctypes.c_void_p" << endl;
        
        fout << "        func1.argtypes=[ctypes.c_void_p]" << endl;
        if (true) {
          fout << "        ptr=func1(self._ptr)" << endl;
          std::string type_temp=ifv.ift.name;
          if (ifv.ift.py_name.length()!=0) {
            type_temp=ifv.ift.py_name;
          }
          if (type_temp.substr(type_temp.length()-2,2)==
              (string)"<>") {
            type_temp=type_temp.substr(0,type_temp.length()-2);
          }
          for(cpn_it it=class_py_names.begin();
              it!=class_py_names.end();it++) {
            if (it->first==type_temp) type_temp=it->second;
          }
          if (type_temp.substr(0,7)==(string)"o2scl::") {
            type_temp=type_temp.substr(7,type_temp.length()-7);
          }
          fout << "        obj=" << type_temp << "(self._link,ptr)"
               << endl;
          fout << "        return obj" << endl;
        } else {
          fout << "        " << ifv.name << "._ptr=func1(self._ptr)" << endl;
          fout << "        " << ifv.name << "._owner=False" << endl;
          fout << "        return" << endl;
        }
        
      }
      fout << endl;

      // Setter
      if (ifv.ift.is_ctype()) {
        
        if (!ifv.ift.is_const() &&
            !ifv.ift.is_reference()) {
          
          fout << "    @" << ifv.name << ".setter" << endl;
          fout << "    def " << ifv.name << "(self,value):" << endl;
          fout << "        \"\"\"" << endl;
          fout << "        Setter function for " << ifc.name << "::"
               << ifv.name << " ." << endl;
          fout << "        \"\"\"" << endl;
          fout << "        func=self._link." << dll_name << "."
               << ifc.ns << "_" << underscoreify(ifc.name)
               << "_set_" << ifv.name << endl;
          fout << "        func.argtypes=[ctypes.c_void_p,ctypes.c_"
               << ifv.ift.name << "]" << endl;
          fout << "        func(self._ptr,value)" << endl;
          fout << "        return" << endl;
          
        }
        
      } else {
        
        fout << "    def set_" << ifv.name << "(self,value):" << endl;
        fout << "        \"\"\"" << endl;
        fout << "        Set object of type :class:`"
             << ifv.ift.name << "`" << endl;
        fout << "        \"\"\"" << endl;
        fout << "        func=self._link." << dll_name << "." << ifc.ns << "_"
             << underscoreify(ifc.name) << "_set_" << ifv.name << endl;
        fout << "        func.argtypes=[ctypes.c_void_p,"
             << "ctypes.c_void_p]" << endl;
        if (ifv.ift.prefix.find("shared_ptr")!=std::string::npos ||
            ifv.ift.prefix.find("std::shared_ptr")!=std::string::npos) {
          fout << "        func(self._ptr,value._s_ptr)" << endl;
        } else {
          fout << "        func(self._ptr,value._ptr)" << endl;
        }
        fout << "        return" << endl;
        
      }
      fout << endl;
    }

    // Define methods
    for(size_t j=0;j<ifc.methods.size();j++) {

      if_func &iff=ifc.methods[j];

      // Function header
      if (iff.name=="operator[]" || iff.name=="operator()") {
        fout << "    def __getitem__(self";
      } else if (iff.py_name!="") {
        fout << "    def " << iff.py_name << "(self";
      } else {
        fout << "    def " << iff.name << "(self";
      }
      for(size_t k=0;k<iff.args.size();k++) {
        if (k==0 && iff.name=="operator()" &&
            (ifc.py_name=="ublas_matrix" ||
             ifc.py_name=="ublas_matrix_int")) {
          fout << ",matrix_tuple";
          k++;
        } else if (!iff.args[k].ift.is_ctype() ||
                   !iff.args[k].ift.is_reference() ||
                   !iff.args[k].ift.is_out()) {
          fout << "," << iff.args[k].name;
          if (iff.args[k].value.length()>0) {
            if (iff.args[k].value=="true") {
              fout << "=True";
            } else if (iff.args[k].value=="false") {
              fout << "=False";
            } else {
              fout << "=" << iff.args[k].value;
            }
          }
        }
      }
      fout << "):" << endl;

      // Begin documentation for member function
      fout << "        \"\"\"" << endl;

      // Parameter list
      if (iff.args.size()>0) {
        fout << "        | Parameters:" << endl;
      }
      for(size_t k=0;k<iff.args.size();k++) {
        if (iff.args[k].ift.is_ctype()) {
          if (iff.args[k].ift.suffix=="&") {
            if (iff.args[k].ift.is_io()) {
              fout << "        | *" << iff.args[k].name << "*";
              if (iff.args[k].value.length()>0) {
                fout << " =" << iff.args[k].value;
              }
              fout << ": ``ctypes.c_" << iff.args[k].ift.name
                   << "``" << endl;
            }
          } else {
            fout << "        | *" << iff.args[k].name << "*";
            if (iff.args[k].value.length()>0) {
              fout << " =" << iff.args[k].value;
            }
            fout << ": ``" << iff.args[k].ift.name << "``" << endl;
          }
        } else if (iff.args[k].ift.suffix=="&") {
          std::string type_temp=iff.args[k].ift.name;
          for(cpn_it it=class_py_names.begin();
              it!=class_py_names.end();it++) {
            if (it->first==type_temp) type_temp=it->second;
          }
          fout << "        | *" << iff.args[k].name << "*";
          if (iff.args[k].value.length()>0) {
            fout << "=" << iff.args[k].value;
          }
          fout << ": :class:`" << type_temp << "` object"
               << endl;
        } else if (iff.args[k].ift.name=="std::string") {
          fout << "        | *" << iff.args[k].name << "*";
          if (iff.args[k].value.length()>0) {
            fout << " =" << iff.args[k].value;
          }
          fout << ": string" << endl;
        }
      }

      // Reformat return type
      size_t len=iff.ret.name.length();
      std::string reformat_ret_type=iff.ret.name;
      // Manually remove '<>' from the typename if necessary
      if (len>2 && iff.ret.name[len-2]=='<' &&
          iff.ret.name[len-1]=='>') {
        reformat_ret_type=iff.ret.name.substr(0,len-2);
      }
      for(cpn_it it=class_py_names.begin();
          it!=class_py_names.end();it++) {
        if (it->first==reformat_ret_type) reformat_ret_type=it->second;
      }
      
      //if (reformat_ret_type=="boost::numeric::ublas::matrix<double>") {
      //reformat_ret_type="ublas_matrix";
      //}
      
      // Depending on return type, set return document string
      // and return python code
      std::string return_docs, restype_string;
      if (iff.name=="operator[]" || iff.name=="operator()") {
        if (iff.ret.name=="std::string") {
          //return_docs="std_string object";
          return_docs="Python bytes object";
          restype_string="ctypes.c_void_p";
        } else if (iff.ret.name=="std::vector<std::string>") {
          return_docs="std_vector_string object";
          restype_string="ctypes.c_void_p";
        } else if ((iff.ret.name!="vector<double>" &&
             iff.ret.name!="std::vector<double>") ||
            iff.ret.suffix!="&") {
          return_docs="";
          restype_string="ctypes.c_"+iff.ret.name;
        } else if (iff.ret.suffix=="&") {
          return_docs=((string)":class:`")+reformat_ret_type+"` object";
          restype_string="ctypes.c_void_p";
        }
      } else if ((iff.ret.name=="vector<double>" ||
                  iff.ret.name=="std::vector<double>") &&
                 iff.ret.suffix=="&") {
        return_docs="``numpy`` array";
        restype_string="";
      } else if (iff.ret.name=="boost::numeric::ublas::vector<double>" &&
                 iff.ret.suffix=="&") {
        return_docs="ublas_vector object";
        restype_string="ctypes.c_void_p";
      } else if (iff.ret.name=="boost::numeric::ublas::vector<int>" &&
                 iff.ret.suffix=="&") {
        return_docs="ublas_vector_int object";
        restype_string="";
      } else if (iff.ret.prefix.find("shared_ptr")!=std::string::npos ||
                 iff.ret.prefix.find("std::shared_ptr")!=std::string::npos) {
        return_docs=((string)":class:`shared_ptr_")+reformat_ret_type+"`.";
        restype_string="ctypes.c_void_p";
      } else if (iff.ret.suffix=="&") {
        return_docs=((string)":class:`")+reformat_ret_type+"` object";
        restype_string="ctypes.c_void_p";
      } else if (iff.ret.name=="std::string") {
        if (iff.ret.is_reference()) {
          return_docs="std_string object";
          restype_string="ctypes.c_void_p";
        } else {
          return_docs="Python bytes object";
          restype_string="ctypes.c_void_p";
        }
      } else if (iff.ret.name=="size_t" || iff.ret.name=="int") {
        return_docs="a Python int";
        restype_string=((string)"ctypes.c_")+reformat_ret_type;
      } else if (iff.ret.name=="float" || iff.ret.name=="double") {
        return_docs="a Python float";
        restype_string=((string)"ctypes.c_")+reformat_ret_type;
      } else if (iff.ret.name=="bool") {
        return_docs="a Python boolean";
        restype_string=((string)"ctypes.c_")+reformat_ret_type;
      } else if (iff.ret.name!="void") {
        return_docs=((string)":class:`")+reformat_ret_type+"` object";
        restype_string="ctypes.c_void_p";
      }
      for(size_t k=0;k<iff.args.size();k++) {
        if (iff.args[k].ift.is_ctype()) {
          if (iff.args[k].ift.suffix=="&") {
            if (iff.args[k].ift.name=="double" ||
                iff.args[k].ift.name=="float") {
              return_docs+=", a Python float";
            } else if (iff.args[k].ift.name=="size_t" ||
                       iff.args[k].ift.name=="int") {
              return_docs+=", a Python int";
            } else {
              return_docs+=", a Python obj";
            }
          }
        }
      }
                            
      // Output return value documentation
      if (return_docs.length()>0) {
        fout << "        | Returns: " << return_docs << endl;
      }

      // End of member function documentation
      fout << "        \"\"\"" << endl;

      // Perform necessary conversions
      for(size_t k=0;k<iff.args.size();k++) {
        if (iff.args[k].ift.name=="std::string") {
          fout << "        " << iff.args[k].name
               << "_=ctypes.c_char_p(force_bytes("
               << iff.args[k].name << "))" << endl;
        }
      }

      // Ctypes function alias
      if (iff.overloaded) {
        fout << "        func=self._link." << dll_name << "."
             << ifc.ns << "_" << underscoreify(ifc.name) << "_"
             << iff.py_name << endl;
      } else if (iff.name=="operator[]" || iff.name=="operator()") {
        fout << "        func=self._link." << dll_name << "."
             << ifc.ns << "_" << underscoreify(ifc.name) << "_getitem"
             << endl;
      } else {
        fout << "        func=self._link." << dll_name << "."
             << ifc.ns << "_" << underscoreify(ifc.name) << "_"
             << iff.name << endl;
      }
      
      // Additional code for returning vector<double> objects
      if ((iff.ret.name=="vector<double>" ||
           iff.ret.name=="std::vector<double>") &&
          iff.ret.suffix=="&") {
        fout << "        n_=ctypes.c_int(0)" << endl;
        fout << "        ptr_=ctypes.POINTER(ctypes.c_double)()" << endl;
      }

      // Specify C interface return type
      if (restype_string.length()>0) {
        fout << "        func.restype=" << restype_string << endl;
      }

      // Extra code before and after the ctypes function call
      vector<string> pre_func_code, post_func_code;
      // Extra return values (for C-type references)
      string addl_ret;
      
      // Ctypes function argument types
      fout << "        func.argtypes=[ctypes.c_void_p";

      // This is set to true if we need to pass the argument by
      // reference and convert it to a python object
      vector<bool> needs_conv(iff.args.size());
      
      for(size_t k=0;k<iff.args.size();k++) {
        needs_conv[k]=false;
        if (iff.args[k].ift.suffix=="&") {
          if (iff.args[k].ift.is_ctype()) {
            // If it's a C type, then we will convert from a python
            // object to a ctypes object and then convert back
            // afterwards
            if (iff.args[k].ift.is_io()) {
              pre_func_code.push_back(iff.args[k].name+"_conv=ctypes.c_"+
                                      iff.args[k].ift.name+
                                      "("+iff.args[k].name+")");
            } else {
              pre_func_code.push_back(iff.args[k].name+"_conv=ctypes.c_"+
                                      iff.args[k].ift.name+
                                      "(0)");
            }
            //post_func_code.push_back(iff.args[k].name+"="+
            //iff.args[k].name+"_conv.value");
            if (addl_ret.length()==0) {
              addl_ret=iff.args[k].name+"_conv.value";
            } else {
              addl_ret+=","+iff.args[k].name+"_conv.value";
            }
            needs_conv[k]=true;
            fout << ",ctypes.POINTER(ctypes.c_"
                 << iff.args[k].ift.name << ")";
          } else {
            fout << ",ctypes.c_void_p";
          }
        } else if (iff.args[k].ift.name=="std::string") {
          fout << ",ctypes.c_char_p";
        } else {
          fout << ",ctypes.c_" << iff.args[k].ift.name;
        }
      }

      if (iff.name=="operator()" &&
          (ifc.py_name=="ublas_matrix" ||
           ifc.py_name=="ublas_matrix_int")) {
        pre_func_code.push_back("m,n=matrix_tuple");
      }
      
      // Instead of returning an array, we send a pointer to the C
      // wrapper
      if ((iff.ret.name=="vector<double>" ||
           iff.ret.name=="std::vector<double>") &&
          iff.ret.suffix=="&") {
        fout << ",ctypes.POINTER(ctypes.POINTER(ctypes.c_double))";
        fout << ",ctypes.POINTER(ctypes.c_int)";
      }
      fout << "]" << endl;

      // Set up the variables which formulate the ctypes function call
      string function_start, function_end;
      
      if ((iff.ret.name=="vector<double>" ||
           iff.ret.name=="std::vector<double>") &&
          iff.ret.suffix=="&") {
        
        function_start="func(self._ptr";
        function_end=",ctypes.byref(ptr_),ctypes.byref(n_))";
        
      } else if (iff.ret.prefix.find("shared_ptr")!=std::string::npos ||
                 iff.ret.prefix.find("std::shared_ptr")!=std::string::npos) {

        function_start="sp=shared_ptr_"+reformat_ret_type+
          "(self._link,func(self._ptr)";
        function_end=")";
        
      } else if (iff.ret.name=="std::string" || iff.ret.name=="string") {
        
        function_start="ret=func(self._ptr";
        function_end=")";
        post_func_code.push_back("strt=std_string(self._link,ret)");
        post_func_code.push_back("strt._owner=True");
        
      } else if (iff.ret.name=="std::vector<std::string>") {
        
        function_start="ret=func(self._ptr";
        function_end=")";
        post_func_code.push_back("vstrt=std_vector_string(self._link,ret)");
        post_func_code.push_back("vstrt._owner=True");
        
      } else if (iff.ret.name=="void") {
        function_start="func(self._ptr";
        function_end=")";
      } else if (iff.ret.is_ctype() || iff.ret.is_reference()) {
        function_start="ret=func(self._ptr";
        function_end=")";
      } else {
        function_start="ret2=func(self._ptr";
        function_end=")";

        std::string ret_temp=iff.ret.name;
        for(cpn_it it=class_py_names.begin();
            it!=class_py_names.end();it++) {
          if (it->first==ret_temp) ret_temp=it->second;
        }
        size_t len=ret_temp.length();
        // Manually remove '<>' from the return type if necessary
        if (len>2 && ret_temp[len-2]=='<' && ret_temp[len-1]=='>') {
          ret_temp=ret_temp.substr(0,len-2);
        }
        
        post_func_code.push_back(((string)"ret=")+ret_temp+
                                 "(self._link,ret2)");
        post_func_code.push_back("ret.owner=True");
      }

      // Write the code for the actual ctypes function call
      
      for(size_t iii=0;iii<pre_func_code.size();iii++) {
        fout << "        " << pre_func_code[iii] << endl;
      }
      
      fout << "        " << function_start;
      
      for(size_t k=0;k<iff.args.size();k++) {
        if (iff.args[k].ift.suffix=="&") {
          if (iff.args[k].ift.is_ctype()) {
            if (needs_conv[k]) {
              fout << ",ctypes.byref(" << iff.args[k].name << "_conv)";
            } else {
              fout << "," << iff.args[k].name;
            }
          } else {
            fout << "," << iff.args[k].name << "._ptr";
          }
        } else if (iff.args[k].ift.name=="std::string") {
          fout << "," << iff.args[k].name << "_";
        } else {
          fout << "," << iff.args[k].name;
        }
      }

      fout << function_end << endl;
      
      for(size_t iii=0;iii<post_func_code.size();iii++) {
        fout << "        " << post_func_code[iii] << endl;
      }
      
      // Python code to set up the return value
      if ((iff.ret.name=="vector<double>" ||
           iff.ret.name=="std::vector<double>") &&
          iff.ret.suffix=="&") {
        fout << "        ret=numpy.ctypeslib.as_array"
             << "(ptr_,shape=(n_.value,))" << endl;
        fout << "        return ret";
        if (addl_ret.length()>0) fout << "," << addl_ret << endl;
        else fout << endl;
      } else if (iff.ret.name=="std::vector<std::string>") {
        fout << "        return vstrt" << endl;
      } else if (iff.ret.name=="std::string" || iff.ret.name=="string") {
        if (iff.name=="operator[]") {
          fout << "        return strt.to_bytes()" << endl;
        } else if (iff.ret.is_reference()) {
          fout << "        return strt";
          if (addl_ret.length()>0) fout << "," << addl_ret << endl;
          else fout << endl;
        } else {
          fout << "        return strt.to_bytes()";
          if (addl_ret.length()>0) fout << "," << addl_ret << endl;
          else fout << endl;
        }
      } else if (iff.ret.suffix=="&" && iff.name!="operator[]" &&
                 iff.name!="operator()") {
        fout << "        ret2=" 
             << underscoreify(reformat_ret_type) << "(self._link,ret)"
             << endl;
        fout << "        return ret2";
        if (addl_ret.length()>0) fout << "," << addl_ret << endl;
        else fout << endl;
      } else if (iff.ret.prefix.find("shared_ptr")!=std::string::npos ||
                 iff.ret.prefix.find("std::shared_ptr")!=std::string::npos) {
        fout << "        return sp";
        if (addl_ret.length()>0) fout << "," << addl_ret << endl;
        else fout << endl;
      } else if (iff.ret.name=="void") {
        fout << "        return";
        if (addl_ret.length()>0) fout << " " << addl_ret << endl;
        else fout << endl;
      } else {
        fout << "        return ret";
        if (addl_ret.length()>0) fout << "," << addl_ret << endl;
        else fout << endl;
      }
      fout << endl;

      // For operator[] functions, __getitem__ python code was already
      // taken care of. Here, we take care of the __setitem__ python
      // code.
      if (iff.name=="operator[]" && !iff.ret.is_const() &&
          iff.ret.suffix=="&") {
        fout << "    def __setitem__(self,i,value):" << endl;
        fout << "        \"\"\"" << endl;
        fout << "        | Parameters:" << endl;
        fout << "        | *i*: ``size_t``" << endl;
        if (iff.ret.name=="std::string") {
          fout << "        | *value*: Python bytes object" << endl;
        } else if (iff.ret.name=="std::vector<double>") {
          fout << "        | *value*: Python array" << endl;
        } else if (iff.ret.name=="std::vector<std::string>") {
          fout << "        | *value*: std_vector_string object" << endl;
        } else {
          fout << "        | *value*: ``" << iff.ret.name << "``" << endl;
        }
        fout << "        \"\"\"" << endl;
        fout << "        func=self._link." << dll_name << "."
             << ifc.ns << "_" << underscoreify(ifc.name) << "_setitem"
             << endl;
        if (iff.ret.name=="std::string") {
          fout << "        func.argtypes=[ctypes.c_void_p,"
               << "ctypes.c_size_t,ctypes.c_void_p]"
               << endl;
          fout << "        s=std_string(self._link)" << endl;
          fout << "        s.init_bytes(value)" << endl;
          fout << "        func(self._ptr,i,s._ptr)" << endl;
        } else if (iff.ret.name=="std::vector<std::string>") {
          fout << "        func.argtypes=[ctypes.c_void_p,"
               << "ctypes.c_size_t,ctypes.c_void_p]"
               << endl;
          fout << "        func(self._ptr,i,value._ptr)" << endl;
        } else if (iff.ret.name=="std::vector<double>") {
          fout << "        sv=std_vector(self._link)" << endl;
          fout << "        sv.resize(len(value))" << endl;
          fout << "        for j in range(0,len(value)):" << endl;
          fout << "            sv[j]=value[j]" << endl;
          fout << "        func.argtypes=[ctypes.c_void_p,"
               << "ctypes.c_size_t,ctypes.c_void_p]"
               << endl;
          fout << "        func(self._ptr,i,sv._ptr)" << endl;
        } else {
          fout << "        func.argtypes=[ctypes.c_void_p,"
               << "ctypes.c_size_t,ctypes.c_" << iff.ret.name << "]"
               << endl;
          fout << "        func(self._ptr,i,value)" << endl;
        }
        fout << "        return" << endl;
        fout << endl;
      }

      // For operator() functions, __getitem__ python code was already
      // taken care of. Here, we take care of the __setitem__ python
      // code.
      if (iff.name=="operator()" && !iff.ret.is_const() &&
          iff.ret.suffix=="&") {
        fout << "    def __setitem__(self,matrix_tuple,value):" << endl;
        fout << "        m,n=matrix_tuple" << endl;
        fout << "        func=self._link." << dll_name << "."
             << ifc.ns << "_" << underscoreify(ifc.name) << "_setitem"
             << endl;
        fout << "        func.argtypes=[ctypes.c_void_p,"
             << "ctypes.c_size_t,ctypes.c_size_t,ctypes.c_"
             << iff.ret.name << "]" << endl;
        fout << "        func(self._ptr,m,n,value)" << endl;
        fout << "        return" << endl;
        fout << endl;
      }

      // End of python code generation for this method
    }

    // Python code generation for additional constructors
    
    for(size_t j=0;j<ifc.cons.size();j++) {
      
      if_func &iff=ifc.cons[j];

      fout << "    @classmethod" << endl;

      // The function header
      fout << "    def " << iff.name << "(cls,link,";
      for(size_t k=0;k<iff.args.size();k++) {
        fout << iff.args[k].name;
        if (k!=iff.args.size()-1) {
          fout << ",";
        }
      }
      fout << "):" << endl;

      // The function documentation
      fout << "        \"\"\"" << endl;
      fout << "        Constructor-like class method for "
           << ifc.name << " ." << endl;
      fout << endl;
      fout << "        | Parameters:" << endl;
      fout << endl;
      fout << "        \"\"\"" << endl;
      fout << endl;

      // The C wrapper function from the DLL
      fout << "        f=link." << dll_name << "." << ifc.ns << "_"
           << underscoreify(ifc.name) << "_" << iff.name << endl;

      // Output the constructor return type
      fout << "        f.restype=ctypes.c_void_p" << endl;

      // Output the constructor argument types
      fout << "        f.argtypes=[";
      for(size_t k=0;k<iff.args.size();k++) {
        if (iff.args[k].ift.suffix=="&") {
          fout << "ctypes.c_void_p";
        } else if (iff.args[k].ift.name=="std::string") {
          fout << "ctypes.c_char_p";
        } else {
          fout << "ctypes.c_" << iff.args[k].ift.name;
        }
        if (k!=iff.args.size()-1) {
          fout << ",";
        }
      }
      fout << "]" << endl;

      // Output the constructor function call
      fout << "        return cls(link,f(";
      // Arguments 
      for(size_t k=0;k<iff.args.size();k++) {
        if (iff.args[k].ift.suffix=="&") {
          fout << iff.args[k].name << "._ptr";
        } else if (iff.args[k].ift.name=="std::string") {
          fout << iff.args[k].name << "_";
        } else {
          fout << iff.args[k].name;
        }
        if (k!=iff.args.size()-1) fout << ",";
      }
      fout << "))" << endl;
      fout << endl;
      
    }    

    // Output any additional python code for this class
    for(size_t j=0;j<ifc.extra_py.size();j++) {
      fout << "    " << ifc.extra_py[j] << endl;
    }
    fout << endl;

    // End of loop over classes
  }

  // Python code for shared pointers
  for(size_t i=0;i<sps.size();i++) {
    
    if_shared_ptr &ifsp=sps[i];
    
    if (ifsp.py_name!="") {
      fout << "class shared_ptr_" << ifsp.py_name << "("
           << ifsp.py_name << "):" << endl;
    } else {
      fout << "class shared_ptr_" << ifsp.name << "("
           << ifsp.name << "):" << endl;
    }
    fout << "    \"\"\"" << endl;
    fout << "    Python interface for a shared pointer to a class of "
         << "type ``" << ifsp.name << "``" << endl;
    fout << "    \"\"\"" << endl;
    
    fout << endl;
    
    // Initialize pointer
    fout << "    _s_ptr=0" << endl;
    fout << "    _ptr=0" << endl;
    fout << "    _link=0" << endl;
    fout << "    _owner=True" << endl;
    fout << endl;
    
    // Define __init__() function
    fout << "    def __init__(self,link,shared_ptr=0):" << endl;
    fout << "        \"\"\"" << endl;
    if (ifsp.py_name!="") {
      fout << "        Init function for shared_ptr_"
           << ifsp.py_name << " ." << endl;
    } else {
      fout << "        Init function for shared_ptr_"
           << ifsp.name << " ." << endl;
    }
    fout << "        \"\"\"" << endl;
    fout << endl;
    fout << "        self._link=link" << endl;
    fout << "        if shared_ptr==0:" << endl;
    fout << "            f2=self._link." << dll_name << "." << ifsp.ns
         << "_create_shared_ptr_" << underscoreify(ifsp.name) << endl;
    fout << "            f2.restype=ctypes.c_void_p" << endl;
    fout << "            self._s_ptr=f2()" << endl;
    fout << "        else:" << endl;
    fout << "            self._s_ptr=shared_ptr" << endl;
    fout << endl;
    fout << "        f=self._link." << dll_name << "." << ifsp.ns
         << "_shared_ptr_" << underscoreify(ifsp.name) << "_ptr" << endl;
    fout << "        f.argtypes=[ctypes.c_void_p]" << endl;
    fout << "        f.restype=ctypes.c_void_p" << endl;
    fout << "        self._ptr=f(self._s_ptr)" << endl;
    fout << "        return" << endl;
    fout << endl;
    
    // Define __del__() function
    fout << "    def __del__(self):" << endl;
    fout << "        \"\"\"" << endl;
    if (ifsp.py_name!="") {
      fout << "        Delete function for shared_ptr_"
           << ifsp.py_name << " ." << endl;
    } else {
      fout << "        Delete function for shared_ptr_"
           << ifsp.name << " ." << endl;
    }
    fout << "        \"\"\"" << endl;
    fout << endl;
    fout << "        f=self._link." << dll_name << "." << ifsp.ns
         << "_free_shared_ptr_"
         << underscoreify(ifsp.name) << endl;
    fout << "        f.argtypes=[ctypes.c_void_p]" << endl;
    fout << "        f(self._s_ptr)" << endl;
    fout << "        return" << endl;
    fout << endl;

  }
  
  // Define functions
  for(size_t j=0;j<functions.size();j++) {
    
    if_func &iff=functions[j];
    
    // Function header
    if (iff.overloaded || iff.py_name.length()>0) {
      fout << "def " << iff.py_name << "(link,";
    } else {
      fout << "def " << iff.name << "(link,";
    }
    for(size_t k=0;k<iff.args.size();k++) {
      fout << iff.args[k].name;
      if (iff.args[k].value.length()>0) {
        if (iff.args[k].value=="true") {
          fout << "=True";
        } else if (iff.args[k].value=="false") {
          fout << "=False";
        } else {
          fout << "=" << iff.args[k].value;
        }
      }
      if (k!=iff.args.size()-1) {
        fout << ",";
      }
    }
    fout << "):" << endl;

    // ----------------------------------------------------
    // Documentation for python code for a function
    
    fout << "    \"\"\"" << endl;
    
    fout << "        | Parameters:" << endl;
    fout << "        | *link* :class:`linker` object" << endl;
    for(size_t k=0;k<iff.args.size();k++) {
      if (iff.args[k].ift.is_ctype()) {
        fout << "        | *" << iff.args[k].name
             << "*: ``" << iff.args[k].ift.name << "``" << endl;
      } else if (iff.args[k].ift.suffix=="&") {
        fout << "        | *" << iff.args[k].name
             << "*: :class:`" << iff.args[k].ift.name << "` object"
             << endl;
      } else if (iff.args[k].ift.name=="std::string") {
        fout << "        | *" << iff.args[k].name
             << "*: string" << endl;
      }
    }
    
    if ((iff.ret.name=="vector<double>" ||
         iff.ret.name=="std::vector<double>") &&
        iff.ret.suffix=="&") {
      fout << "        | Returns: ``numpy`` array" << endl;
    } else if (iff.ret.prefix.find("shared_ptr")!=std::string::npos ||
               iff.ret.prefix.find("std::shared_ptr")!=std::string::npos) {
      size_t len=iff.ret.name.length();
      std::string tmps=iff.ret.name;
      // Manually remove '<>' from the typename if necessary
      if (len>2 && iff.ret.name[len-2]=='<' &&
          iff.ret.name[len-1]=='>') {
        tmps=iff.ret.name.substr(0,len-2);
      }
      fout << "        | Returns: :class:`"
           << "shared_ptr_" << tmps << "`." << endl;
    } else if (iff.ret.suffix=="&") {
      size_t len=iff.ret.name.length();
      std::string tmps=iff.ret.name;
      // Manually remove '<>' from the typename if necessary
      if (len>2 && iff.ret.name[len-2]=='<' &&
          iff.ret.name[len-1]=='>') {
        tmps=iff.ret.name.substr(0,len-2);
      }
      fout << "        | Returns: :class:`"
           << tmps << "` object" << endl;
    } else if (iff.ret.name=="std::string") {
      fout << "        | Returns: python bytes object"
           << endl;
    } else if (iff.ret.is_ctype()) {
      fout << "        | Returns: ``ctypes.c_"
           << iff.ret.name << "`` object" << endl;
    } else if (iff.ret.name!="void") {

      size_t len=iff.ret.name.length();
      std::string reformat_ret_type=iff.ret.name;
      // Manually remove '<>' from the typename if necessary
      if (len>2 && iff.ret.name[len-2]=='<' &&
          iff.ret.name[len-1]=='>') {
        reformat_ret_type=iff.ret.name.substr(0,len-2);
      }
      for(cpn_it it=class_py_names.begin();
          it!=class_py_names.end();it++) {
        if (it->first==reformat_ret_type) reformat_ret_type=it->second;
      }
      
      fout << "        | Returns: ``"
           << reformat_ret_type << "`` object" << endl;
    }
      
    fout << "    \"\"\"" << endl;

    // End of documention of python code for a function
    
    // ----------------------------------------------------
    // Perform necessary conversions
    
    for(size_t k=0;k<iff.args.size();k++) {
      if (iff.args[k].ift.name=="std::string") {
        if (iff.args[k].ift.suffix=="&") {
          fout << "    " << iff.args[k].name << ".__del__()" << endl;
          fout << "    " << iff.args[k].name
               << "._ptr=ctypes.c_void_p()" << endl;
        } else {
          fout << "    " << iff.args[k].name
               << "_=ctypes.c_char_p(force_bytes("
               << iff.args[k].name << "))" << endl;
        }
      }
    }
      
    // Ctypes interface for function
    
    if (iff.overloaded) {
      fout << "    func=link." << dll_name << "." << iff.ns << "_"
           << iff.py_name << "_wrapper" << endl;
    } else {
      fout << "    func=link." << dll_name << "." << iff.ns << "_"
           << underscoreify(iff.name) << "_wrapper" << endl;
    }
    if (iff.ret.name!="void") {
      if (iff.ret.name=="std::string") {
        fout << "    func.restype=ctypes.c_char_p" << endl;
      } else {

        if (iff.ret.is_ctype()) {
          fout << "    func.restype=ctypes.c_" << iff.ret.name << endl;
        } else {
          size_t len=iff.ret.name.length();
          std::string reformat_ret_type=iff.ret.name;
          // Manually remove '<>' from the typename if necessary
          if (len>2 && iff.ret.name[len-2]=='<' &&
              iff.ret.name[len-1]=='>') {
            reformat_ret_type=iff.ret.name.substr(0,len-2);
          }
          for(cpn_it it=class_py_names.begin();
              it!=class_py_names.end();it++) {
            if (it->first==reformat_ret_type) reformat_ret_type=it->second;
          }
          
          fout << "    func.restype=" << reformat_ret_type << endl;
        }
      }
    }
    fout << "    func.argtypes=[";
    for(size_t k=0;k<iff.args.size();k++) {
      if (iff.args[k].ift.suffix=="&") {
        if (iff.args[k].ift.name=="std::string") {
          fout << "ctypes.POINTER(ctypes.c_void_p)";
        } else {
          fout << "ctypes.c_void_p";
        }
      } else if (iff.args[k].ift.name=="std::string") {
        fout << "ctypes.c_char_p";
      } else {
        fout << "ctypes.c_" << iff.args[k].ift.name;
      }
      if (k!=iff.args.size()-1) {
        fout << ",";
      }
    }
    fout << "]" << endl;
    
    // Call C++ wrapper function
    if (iff.ret.name=="void") {
      fout << "    func(";
    } else {
      fout << "    ret=func(";
    }
    for(size_t k=0;k<iff.args.size();k++) {
      if (iff.args[k].ift.suffix=="&") {
        if (iff.args[k].ift.name=="std::string") {
          fout << "ctypes.byref(" << iff.args[k].name << "._ptr)";
        } else {
          fout << iff.args[k].name << "._ptr";
        }
      } else if (iff.args[k].ift.name=="std::string") {
        fout << iff.args[k].name << "_";
      } else {
        fout << iff.args[k].name;
      }
      if (k!=iff.args.size()-1) fout << ",";
    }
    fout << ")" << endl;
    
    for(size_t k=0;k<iff.args.size();k++) {
      if (iff.args[k].ift.suffix=="&") {
        if (iff.args[k].ift.name=="std::string") {
          fout << "    " << iff.args[k].name << "._owner=True" << endl;
        }
      }
    }
    
    // Return
    if (iff.ret.name=="void") {
      fout << "    return" << endl;
    } else {
      fout << "    return ret" << endl;
    }
    fout << endl;
    
  }

  fout.close();

  // ----------------------------------------------------------------
  // Create rst file for python documentation

  ofstream fout2;
  fout2.open((rst_prefix+".rst").c_str());

  for(size_t j=0;j<rst_header.size();j++) {
    fout2 << rst_header[j] << endl;
  }
  fout2 << endl;

  // First, begin with a "table of contents" which lists
  // all classes and functions
  
  for(size_t i=0;i<classes.size();i++) {

    if_class &ifc=classes[i];

    fout2 << "* :ref:`Class ";
    if (ifc.py_name=="") {
      fout2 << ifc.name << "`" << endl;
    } else {
      fout2 << ifc.py_name << "`" << endl;
    }
    
  }

  for(size_t i=0;i<sps.size();i++) {

    if_shared_ptr &ifsp=sps[i];

    fout2 << "* :ref:`Class shared_ptr_";
    
    if (ifsp.py_name=="") {
      fout2 << ifsp.name << "`" << endl;
    } else {
      fout2 << ifsp.py_name << "`" << endl;
    }

  }
  
  for(size_t j=0;j<functions.size();j++) {
    
    if_func &iff=functions[j];

    fout2 << "* :ref:`Function ";
    
    if (iff.overloaded || iff.py_name.length()>0) {
      fout2 << iff.py_name << "`" << endl;
    } else {
      fout2 << iff.name << "`" << endl;
    }
  }

  fout2 << endl;

  // Now go through each object invidiually
  
  for(size_t i=0;i<classes.size();i++) {

    if_class &ifc=classes[i];

    // This is the section header for the class. The variable 'len'
    // counts the number of characters so that we can match with and
    // equal number of dashes.

    size_t len=6;
    fout2 << "Class ";
    if (ifc.py_name=="") {
      fout2 << ifc.name << endl;
      len+=ifc.name.length();
    } else {
      fout2 << ifc.py_name << endl;
      len+=ifc.py_name.length();
    }
    for(size_t kk=0;kk<len;kk++) {
      fout2 << "-";
    }
    fout2 << "\n" << endl;

    // Construct the 'autoclass' documentation
    if (ifc.py_name=="") {
      fout2 << ".. autoclass:: o2sclpy." << ifc.name << endl;
    } else {
      fout2 << ".. autoclass:: o2sclpy." << ifc.py_name << endl;
    }
    fout2 << "        :members:" << endl;
    fout2 << "        :undoc-members:" << endl;
    fout2 << endl;
    fout2 << "        .. automethod:: __init__" << endl;
    fout2 << "        .. automethod:: __del__" << endl;
    fout2 << "        .. automethod:: __copy__" << endl;
    if (ifc.std_cc) {
      fout2 << "        .. automethod:: __deepcopy__" << endl;
    }
    
    // If the class contains an operator[] or an operator(), then ensure
    // that the __getitem__ and __setitem__ methods are documented.
    for(size_t k=0;k<ifc.methods.size();k++) {
      if (ifc.methods[k].name=="operator[]" ||
          ifc.methods[k].name=="operator()") {
        fout2 << "        .. automethod:: __getitem__" << endl;
        if (!ifc.methods[k].ret.is_const() &&
            ifc.methods[k].ret.suffix=="&") {
          fout2 << "        .. automethod:: __setitem__" << endl;
        }
        k=ifc.methods.size();
      }
    }
    fout2 << endl;
  }

  for(size_t i=0;i<sps.size();i++) {

    if_shared_ptr &ifsp=sps[i];

    // This is the section header for the shared pointer. The variable
    // 'len' counts the number of characters so that we can match with
    // and equal number of dashes.
    
    size_t len=17;
    fout2 << "Class shared_ptr_";
    
    if (ifsp.py_name=="") {
      fout2 << ifsp.name << endl;
      len+=ifsp.name.length();
    } else {
      fout2 << ifsp.py_name << endl;
      len+=ifsp.py_name.length();
    }
    for(size_t kk=0;kk<len;kk++) {
      fout2 << "-";
    }
    fout2 << "\n" << endl;

    // Output the 'autoclass' command
    
    if (ifsp.py_name!="") {
      
      // Manually remove '<>' from the typename if necessary
      size_t len=ifsp.name.length();
      std::string tmps=ifsp.name;
      if (len>2 && ifsp.name[len-2]=='<' &&
          ifsp.name[len-1]=='>') {
        tmps=ifsp.name.substr(0,len-2);
      }
      
      fout2 << ".. autoclass:: o2sclpy.shared_ptr_" << tmps << endl;
      fout2 << "        :members:" << endl;
      fout2 << "        :undoc-members:\n" << endl;
      fout2 << "        .. automethod:: __init__" << endl;
      fout2 << "        .. automethod:: __del__" << endl;
      fout2 << endl;
      
    } else {
      
      fout2 << ".. autoclass:: o2sclpy.shared_ptr_" << ifsp.py_name << endl;
      fout2 << "        :members:" << endl;
      fout2 << "        :undoc-members:\n" << endl;
    }
    
  }

  // RST help for functions
  for(size_t j=0;j<functions.size();j++) {
    
    if_func &iff=functions[j];

    // This is the section header for the function. The variable 'len'
    // counts the number of characters so that we can match with and
    // equal number of dashes.
    
    size_t len=9;
    fout2 << "Function ";
    
    if (iff.overloaded || iff.py_name.length()>0) {
      fout2 << iff.py_name << endl;
      len+=iff.py_name.length();
    } else {
      fout2 << iff.name << endl;
      len+=iff.name.length();
    }
    for(size_t kk=0;kk<len;kk++) {
      fout2 << "-";
    }
    fout2 << "\n" << endl;
    
    // Now construct the 'autofunction' command
    
    if (iff.overloaded || iff.py_name.length()>0) {
      fout2 << ".. autofunction:: o2sclpy." << iff.py_name << "(link,";
    } else {
      fout2 << ".. autofunction:: o2sclpy." << iff.name << "(link,";
    }      
    for(size_t k=0;k<iff.args.size();k++) {
      fout2 << iff.args[k].name;
      if (k!=iff.args.size()-1) {
        fout2 << ",";
      }
    }
    fout2 << ")\n" << endl;
    
  }    
  
  fout2.close();

  cout << "Yanic completed successfully." << endl;
  
  return 0;
}
