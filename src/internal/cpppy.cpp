/*
  -------------------------------------------------------------------
  
  Copyright (C) 2020-2021, Andrew W. Steiner
  
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

using namespace std;
using namespace o2scl;

/*
  Get and set methods for C++ class members always imply a copy. Right
  now cpppy does not yet support, for example, obtaining a reference
  or a pointer to a c++ class member.

  Todos: need to fix function names in case where there is no namespace.
*/

/** \brief Convert all non-alphanumeric characters to underscores
 */
std::string underscoreify(std::string s) {
  std::string s2=s;
  for(size_t i=0;i<s2.length();i++) {
    if (std::isalnum(s2[i])==false) s2[i]='_';
  }
  return s2;
}

/** \brief Read the next line from \c fin, skipping blank
    lines, splitting by spaces into \c vs, and setting \c done
    if done
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
 */
void parse_vector_string(ifstream &fin, std::string type, std::string &line,
                         std::vector<std::string> &vs, bool &done,
                         std::vector<std::string> &output) {

  //cout << "Here0 " << line << endl;
  if (vs.size()==1) {
    cout << "Clearing " << type << "." << endl;
    output.clear();
    next_line(fin,line,vs,done);
  } else if ((vs[0]==type && vs[1]=="|") ||
             (vs.size()>=3 && vs[1]==type && vs[2]=="|")) {
    //cout << "Here " << line << endl;
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
    //cout << "Here2 " << line << endl;
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
  
  std::string name;
  
};

/** \brief Interface for shared pointer
 */
class if_shared_ptr : public if_base {

public:

  /// Namespace
  std::string ns;
  
  /// Python name for the shared pointer type
  std::string py_name;
  
};

/** \brief Interface type
 */
class if_type : public if_base {
  
public:
  
  /** \brief Type qualifiers

      For example
      - const
      - static 
      - static const

      and 

      The shared_ptr is also used as a prefix to indicate a 
      shared pointer of a particular type.
  */
  std::string prefix;
  
  /// Ampersands and asterisks, *, *&, &, etc.
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
  
  /** \brief Parse a vector string object into the type parts
   */
  void parse(std::vector<std::string> &vs, size_t start, size_t end) {
    
    // Create a new vector string object for convenience
    std::vector<std::string> vs2;
    for(size_t i=start;i<end;i++) {
      vs2.push_back(vs[i]);
    }
    
    // Handle according to the number of arguments
    if (vs2.size()==1) {
      prefix="";
      suffix="";
      name=vs2[0];
    } else if (vs2.size()==2) {
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
    } else {
      cerr << "Unsupported number of type arguments." << endl;
      exit(-1);
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
  }
  
  /// Return true if the type is a reference
  bool is_reference() {
    if (suffix=="&") return true;
  }
  
  /// Return true if the type is const
  bool is_const() {
    if (prefix.find("const")!=std::string::npos) return true;
  }
  
  /// Return true if the type is static
  bool is_static() {
    if (prefix.find("const")!=std::string::npos) return true;
  }

  // End of class if_type
};

/** \brief A variable with a type and a name
 */
class if_var : public if_base {
  
public:

  /// Python name for the variable
  std::string py_name;

  /// The variable type 
  if_type ift;
  
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
  
};

/** \brief Interface class
 */
class if_class : public if_base {
  
public:

  /// Python name for the class
  std::string py_name;

  /// Pattern for the class documentation in python
  std::vector<std::string> py_class_doc;
  
  /// True if the class is abstract
  bool is_abstract;
  
  /// Members
  std::vector<if_var> members;
  
  /// Methods
  std::vector<if_func> methods;
  
  /// Constructors
  std::vector<if_func> cons;
  
  /// List of parent classes
  std::vector<std::string> parents;
  
  /// Namespace
  std::string ns;

  if_class() {
    is_abstract=false;
  }
  
};

int main(int argc, char *argv[]) {

  if (argc<5) {
    cerr << "Args: <interface file> <c++ output prefix> "
         << "<python prefix> <rst prefix>" << endl;
    exit(-1);
  }
  std::string fname=argv[1];
  std::string cpp_prefix=argv[2];
  std::string py_prefix=argv[3];
  std::string rst_prefix=argv[4];

  cout << "Reading interface file " << fname << " ." << endl;
  cout << "Setting C++ output prefix to " << cpp_prefix << " ." << endl;
  cout << "Setting Python output prefix to " << py_prefix << " ." << endl;
  cout << endl;

  vector<string> header_strings=
    {"  -------------------------------------------------------------------",
     "",
     "  Copyright (C) 2020-2021, Andrew W. Steiner",
     "",
     "  This file is part of O2scl.",
     "",
     "  O2scl is free software; you can redistribute it and/or modify",
     "  it under the terms of the GNU General Public License as published by",
     "  the Free Software Foundation; either version 3 of the License, or",
     "  (at your option) any later version.",
     "",
     "  O2scl is distributed in the hope that it will be useful,",
     "  but WITHOUT ANY WARRANTY; without even the implied warranty of",
     "  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the",
     "  GNU General Public License for more details.",
     "",
     "  You should have received a copy of the GNU General Public License",
     "  along with O2scl. If not, see <http://www.gnu.org/licenses/>."
     "",
     "  -------------------------------------------------------------------"
    };
  
  cout << "Parsing interface " << fname << " ." << endl;
  
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
  
  // Open file
  ifstream fin;
  fin.open(fname.c_str());
  
  // Parse file
  std::string line;
  std::vector<std::string> vs;
  bool done=false;

  bool import_numpy=false;
  
  next_line(fin,line,vs,done);
  if (done==true) {
    cerr << "Empty interface file." << endl;
    exit(-1);
  }
  
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
      
      if (line[0]=='-' && line[1]==' ' && vs[1]=="py_name" && vs.size()==3) {
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
      if (vs.size()>=3 && vs[2]=="abstract") {
        ifc.is_abstract=true;
        cout << "Starting abstract class " << ifc.name << " in namespace "
             << ifc.ns << endl;
      } else {
        cout << "Starting class " << ifc.name << " in namespace "
             << ifc.ns << endl;
      }
      
      next_line(fin,line,vs,done);
      
      bool class_done=false;
      
      while (class_done==false) {

        if (vs.size()>=3 && vs[0]=="-" && vs[1]=="function") {
          
          if_func iff;
          iff.name=vs[2];
          if (iff.name=="operator[]") {
            iff.name="index_operator";
          }
          cout << "  Starting member function " << vs[2] << endl;

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
          if ((iff.ret.name=="vector<double>" ||
               iff.ret.name=="std::vector<double>") &&
              iff.ret.suffix=="&") {
            import_numpy=true;
          }

          next_line(fin,line,vs,done);
          
          while (vs.size()>=2 && line[0]==' ' && line[1]==' ' &&
                 vs[0]=="-") {
            
            if_var ifv;
            string last_string=vs[vs.size()-1];
            if (last_string[0]=='&' || last_string[0]=='*') {
              vs[vs.size()-1]="";
              while (last_string[0]=='&' || last_string[0]=='*') {
                vs[vs.size()-1]=last_string[0]+vs[vs.size()-1];
                last_string=last_string.substr(1,last_string.length()-1);
              }
              ifv.name=last_string;
              ifv.ift.parse(vs,1,vs.size());
            } else {
              ifv.name=last_string;
              ifv.ift.parse(vs,1,vs.size()-1);
            }
            cout << "    Member function " << iff.name
                 << " has argument " << ifv.name << " with type "
                 << ifv.ift.to_string() << endl;

            iff.args.push_back(ifv);
            
            next_line(fin,line,vs,done);
            if (done) class_done=true;
          }

          ifc.methods.push_back(iff);
          
        } else if (vs.size()>=3 && vs[0]=="-" && vs[1]=="parent") {
          
          ifc.parents.push_back(vs[2]);
          cout << "  Parent class " << vs[2] << endl;
          
          next_line(fin,line,vs,done);
          if (done) class_done=true;
          
        } else if (vs.size()>=3 && vs[0]=="-" && vs[1]=="py_name") {
          
          ifc.py_name=vs[2];
          cout << "  Class has python name " << vs[2] << endl;
          
          next_line(fin,line,vs,done);
          if (done) class_done=true;
          
        } else if (vs.size()>=3 && vs[0]=="-" &&
                   vs[1]=="py_class_doc") {

          parse_vector_string(fin,"py_class_doc",line,vs,done,ifc.py_class_doc);
          
          if (done) class_done=true;
          
        } else if (vs.size()>=3 && vs[0]=="-") {
          
          if_var ifv;
          // First, check if there are any & or *'s before the
          // variable name
          if (vs[2][0]=='&' || vs[2][0]=='*') {
            ifv.name=vs[2].substr(1,vs[2].length()-1);
            vs[2]=vs[2].substr(0,1);
            ifv.ift.parse(vs,1,3);
          } else {
            ifv.name=vs[2];
            ifv.ift.parse(vs,1,2);
          }
          
          cout << "  Member " << ifv.name << " with type "
               << ifv.ift.to_string() << " ." << endl;

          next_line(fin,line,vs,done);
          if (done) {
            class_done=true;
          } else {
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

        } else {
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
      
      next_line(fin,line,vs,done);
          
      if (done==true || vs.size()<2 || vs[0]!="-" || line[0]!='-' ||
          line[1]!=' ') {
        cerr << "Could not get return value for function "
             << iff.name << endl;
        exit(-1);
      }
          
      iff.ret.parse(vs,1,vs.size());
      
      cout << "  Member function " << iff.name
           << " has return type "
           << iff.ret.to_string() << endl;
      
      next_line(fin,line,vs,done);
      
      bool function_done=false;
      
      while (vs.size()>=2 && line[0]=='-' && line[1]==' ' &&
             vs[0]=="-") {
        
        if_var ifv;
        string last_string=vs[vs.size()-1];
        if (last_string[0]=='&' || last_string[0]=='*') {
          vs[vs.size()-1]="";
          while (last_string[0]=='&' || last_string[0]=='*') {
            vs[vs.size()-1]=last_string[0]+vs[vs.size()-1];
            last_string=last_string.substr(1,last_string.length()-1);
          }
          ifv.name=last_string;
          ifv.ift.parse(vs,1,vs.size());
        } else {
          ifv.name=last_string;
          ifv.ift.parse(vs,1,vs.size()-1);
        }
        cout << "  Function " << iff.name
             << " has argument " << ifv.name << " with type "
             << ifv.ift.to_string() << endl;
        
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

      if (ifc.is_abstract==false) {
        // The create pointer function
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
      
        // The free pointer function
        fout << "void " << underscoreify(ifc.ns) << "_free_"
             << underscoreify(ifc.name) << "(void *vptr)";
        if (header) {
          fout << ";" << endl;
        } else {
          fout << " {" << endl;
          fout << "  " << ifc.name << " *ptr=(" << ifc.name
               << " *)vptr;" << endl;
          fout << "  delete ptr;" << endl;
          fout << "}" << endl;
        }
        fout << endl;
      }

      for(size_t j=0;j<ifc.members.size();j++) {

        if_var &ifv=ifc.members[j];

        // Get function
        
        if (ifv.ift.name=="bool" ||
            ifv.ift.name=="double" ||
            ifv.ift.name=="int" ||
            ifv.ift.name=="size_t") {
          fout << ifv.ift.name << " " << underscoreify(ifc.ns) << "_"
               << underscoreify(ifc.name) << "_get_" << ifv.name
               << "(void *vptr)";
          if (header) {
            fout << ";" << endl;
          } else {
            fout << " {" << endl;
            fout << "  " << ifc.name << " *ptr=(" << ifc.name
                 << " *)vptr;" << endl;
            fout << "  return ptr->" << ifv.name << ";" << endl;
            fout << "}" << endl;
          }
        } else if (ifv.ift.name=="std::string") {
          fout << "const char *" << underscoreify(ifc.ns) << "_"
               << underscoreify(ifc.name) << "_get_" << ifv.name
               << "(void *vptr)";
          if (header) {
            fout << ";" << endl;
          } else {
            fout << " {" << endl;
            fout << "  " << ifc.name << " *ptr=(" << ifc.name
                 << " *)vptr;" << endl;
            fout << "  python_temp_string=ptr->" << ifv.name
                 << ";" << endl;
            fout << "  return python_temp_string.c_str();" << endl;
            fout << "}" << endl;
          }
          
        } else {
          fout << "void " << underscoreify(ifc.ns) << "_"
               << underscoreify(ifc.name) << "_get_" << ifv.name
               << "(void *vptr, void *p_v)";
          if (header) {
            fout << ";" << endl;
          } else {
            fout << " {" << endl;
            fout << "  " << ifc.name << " *ptr=(" << ifc.name
                 << " *)vptr;" << endl;
            fout << "  " << ifv.ift.name << " *p_t=("
                 << ifv.ift.name << " *)p_v;" << endl;
            fout << "  *(p_t)=ptr->" << ifv.name << ";" << endl;
            fout << "  return;" << endl;
            fout << "}" << endl;
          }
        }
        fout << endl;
      
        // Set function declaration
        if (ifv.ift.name=="bool" ||
            ifv.ift.name=="double" ||
            ifv.ift.name=="int" ||
            ifv.ift.name=="size_t") {
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
        } else {
          fout << "void " << underscoreify(ifc.ns) << "_"
               << underscoreify(ifc.name) << "_set_" << ifv.name
               << "(void *vptr, void *p_v)";
          if (header) {
            fout << ";" << endl;
          } else {
            fout << " {" << endl;
            fout << "  " << ifc.name << " *ptr=(" << ifc.name
                 << " *)vptr;" << endl;
            fout << "  " << ifv.ift.name << " *p_t=("
                 << ifv.ift.name << " *)p_v;" << endl;
            fout << "  ptr->" << ifv.name << "=*(p_t);" << endl;
            fout << "  return;" << endl;
            fout << "}" << endl;
          }
        }
        fout << endl;
      }
    
      for(size_t j=0;j<ifc.methods.size();j++) {
      
        if_func &iff=ifc.methods[j];

        string ret_type;
        string extra_args;
        
        // Function header, first determine return type and
        // any extra arguments
        if (iff.ret.name=="std::string") {
          ret_type="const char *";
        } else if ((iff.ret.name=="vector<double>" ||
                    iff.ret.name=="std::vector<double>") &&
                   iff.ret.suffix=="&") {
          // For a std::vector<double> &, we return a void, but add
          // double pointer and integer output parameters
          ret_type="void ";
          extra_args=", double **dptr, int *n";
        } else if (iff.ret.name=="void" ||
                   iff.ret.name=="bool" ||
                   iff.ret.name=="double" ||
                   iff.ret.name=="int" ||
                   iff.ret.name=="size_t") {
          ret_type=iff.ret.name+" ";
        } else {
          ret_type="void *";
        }

        // Now generate the actual code
        fout << ret_type << underscoreify(ifc.ns) << "_"
             << underscoreify(ifc.name) << "_" << iff.name << "(void *vptr";
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
            fout << "void *ptr_" << iff.args[k].name;
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
          
          // Pointer assignments for arguments
          for(size_t k=0;k<iff.args.size();k++) {
            if (iff.args[k].ift.suffix=="&") {
              fout << "  " << iff.args[k].ift.name << " *"
                   << iff.args[k].name << "=("
                   << iff.args[k].ift.name << " *)ptr_" << iff.args[k].name
                   << ";" << endl;
            }
          }
          
          // Now generate code for actual function call and the
          // return statement
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
            // For a std::vector<double> &, just return a pointer
            if (iff.name=="index_operator") {
              // In case it's const, we have to explicitly typecast
              fout << "  *dptr=(double *)(&(ptr->operator[](";
            } else {
              fout << "  *dptr=ptr->" << iff.name << "(";
            }
          } else if (iff.ret.name=="void") {
            fout << "  ptr->" << iff.name << "(";
          } else {
            // If the function returns a reference, return a pointer
            // instead
            if (iff.ret.suffix=="&") {
              fout << "  " << iff.ret.name << " *ret=&ptr->" << iff.name << "(";
            } else {
              fout << "  " << iff.ret.name << " ret=ptr->" << iff.name << "(";
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
          
          // Extra code for array types
          if ((iff.ret.name=="vector<double>" ||
               iff.ret.name=="std::vector<double>") &&
              iff.ret.suffix=="&") {
            fout << ")[0]));" << endl;
            fout << "  *n=ptr->get_nlines();" << endl;
            fout << "  return;" << endl;
          } else {
            fout << ");" << endl;
            
            if (iff.ret.name=="std::string") {
              fout << "  python_temp_string=ret;" << endl;
              fout << "  return python_temp_string.c_str();" << endl;
            } else if (iff.ret.name=="void") {
              fout << "  return;" << endl;
            } else {
              fout << "  return ret;" << endl;
            }
          }
          
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
        fout << "  " << ifsp.name << " *ptr=new " << ifsp.name << ";" << endl;
        fout << "  return ptr;" << endl;
        fout << "}" << endl;
        fout << endl;        
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
        fout << "  " << ifsp.name << " *ptr=(" << ifsp.name
             << " *)vptr;" << endl;
        fout << "  delete ptr;" << endl;
        fout << "}" << endl;
        fout << endl;        
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
      if (iff.ret.name=="std::string") {
        fout << "const char *" << underscoreify(iff.ns) << "_"
             << iff.name << "(";
      } else if (iff.ret.name=="void" ||
                 iff.ret.name=="bool" ||
                 iff.ret.name=="double" ||
                 iff.ret.name=="int" ||
                 iff.ret.name=="size_t") {
        fout << iff.ret.name << " " << underscoreify(iff.ns) << "_"
             << iff.name << "_wrapper(";
      } else {
        fout << "void *" << underscoreify(iff.ns) << "_"
             << iff.name << "_wrapper(";
      }
    
      for(size_t k=0;k<iff.args.size();k++) {
        if (iff.args[k].ift.suffix=="") {
          if (iff.args[k].ift.name=="std::string") {
            fout << "char *" << iff.args[k].name;
          } else {
            fout << iff.args[k].ift.name << " " << iff.args[k].name;
          }
        } else if (iff.args[k].ift.suffix=="&") {
          fout << "void *ptr_" << iff.args[k].name;
        }
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
            fout << "  " << iff.args[k].ift.name << " *"
                 << iff.args[k].name << "=("
                 << iff.args[k].ift.name << " *)ptr_" << iff.args[k].name
                 << ";" << endl;
          }
        }
        
        // Now generate code for actual function call and the
        // return statement
        if (iff.ret.name=="void") {
          
          fout << "  " << iff.name << "(";
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
          
          fout << "  return;" << endl;
          
        } else {
          
          fout << "  " << iff.ret.name << " ret=ptr->" << iff.name << "(";
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
          
          if (iff.ret.name=="std::string") {
            fout << "  python_temp_string=ret;" << endl;
            fout << "  return python_temp_string.c_str();" << endl;
          } else {
            fout << "  return ret;" << endl;
          }
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
  
  for(size_t i=0;i<py_headers.size();i++) {
    fout << py_headers[i] << endl;
  }
  fout << endl;
  
  for(size_t i=0;i<classes.size();i++) {
    
    if_class &ifc=classes[i];
    
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
    
    // Initialize pointer
    if (ifc.parents.size()==0) {
      fout << "    _ptr=0" << endl;
      fout << "    _link=0" << endl;
      fout << "    _owner=True" << endl;
      fout << endl;
    }
    
    // Define __init__() function
    if (ifc.is_abstract) {
      fout << "    @abstractmethod" << endl;
    }
    fout << "    def __init__(self,link,pointer=0):" << endl;
    fout << "        \"\"\"" << endl;
    fout << "        Init function for class " << ifc.name << " ." << endl;
    fout << "        \"\"\"" << endl;
    fout << endl;
    fout << "        if pointer==0:" << endl;
    fout << "            f=link." << dll_name << "." << ifc.ns << "_create_"
         << underscoreify(ifc.name) << endl;
    fout << "            f.restype=ctypes.c_void_p" << endl;
    fout << "            f.argtypes=[]" << endl;
    fout << "            self._ptr=f()" << endl;
    fout << "        else:" << endl;
    fout << "            self._ptr=pointer" << endl;
    fout << "            self._owner=False" << endl;
    fout << "        self._link=link" << endl;
    fout << "        return" << endl;
    fout << endl;
    
    // Define __del__() function
    fout << "    def __del__(self):" << endl;
    fout << "        \"\"\"" << endl;
    fout << "        Delete function for class " << ifc.name << " ." << endl;
    fout << "        \"\"\"" << endl;
    fout << endl;
    fout << "        if self._owner==True:" << endl;
    fout << "            f=self._link." << dll_name << "." << ifc.ns << "_free_"
         << underscoreify(ifc.name) << endl;
    fout << "            f.argtypes=[ctypes.c_void_p]" << endl;
    fout << "            f(self._ptr)" << endl;
    fout << "        return" << endl;
    fout << endl;

    // Define member get and set properties
    for(size_t j=0;j<ifc.members.size();j++) {
      if_var &ifv=ifc.members[j];

      // Hack, because del has a special meaning in python
      if (ifv.name=="del") ifv.name="delta";

      // Getter
      if (ifv.ift.name=="bool" ||
          ifv.ift.name=="double" ||
          ifv.ift.name=="int" ||
          ifv.ift.name=="size_t") {
        fout << "    @property" << endl;
        fout << "    def " << ifv.name << "(self):" << endl;
        fout << "        \"\"\"" << endl;
        fout << "        Property of type ctypes.c_" << ifv.ift.name << endl;
        fout << "        \"\"\"" << endl;
        fout << "        func=self._link." << dll_name << "." << ifc.ns << "_"
             << underscoreify(ifc.name)
             << "_get_" << ifv.name << endl;
        fout << "        func.restype=ctypes.c_" << ifv.ift.name << endl;
        fout << "        func.argtypes=[ctypes.c_void_p]" << endl;
        fout << "        return func(self._ptr)" << endl;
      } else {
        fout << "    def get_" << ifv.name << "(self," << ifv.name
             << "):" << endl;
        fout << "        \"\"\"" << endl;
        fout << "        Object of type :class:`"
             << ifv.ift.name << "`" << endl;
        fout << "        \"\"\"" << endl;
        fout << "        func=self._link." << dll_name << "." << ifc.ns << "_"
             << underscoreify(ifc.name)
             << "_get_" << ifv.name << endl;
        
        // If it's not a string, then we don't need a return type
        if (ifv.ift.name=="std::string") {
          fout << "        func.restype=ctypes.c_char_p" << endl;
        }
        
        fout << "        func.argtypes=[ctypes.c_void_p," 
             << "ctypes.c_void_p]" << endl;
        fout << "        func(self._ptr," << ifv.name
             << "._ptr)" << endl;
        fout << "        return" << endl;
      }
      fout << endl;

      // Setter
      if (ifv.ift.name=="bool" ||
          ifv.ift.name=="double" ||
          ifv.ift.name=="int" ||
          ifv.ift.name=="size_t") {
        fout << "    @" << ifv.name << ".setter" << endl;
        fout << "    def " << ifv.name << "(self,value):" << endl;
        fout << "        \"\"\"" << endl;
        fout << "        Setter function for " << ifc.name << "::"
             << ifv.name << " ." << endl;
        fout << "        \"\"\"" << endl;
        fout << "        func=self._link." << dll_name << "." << ifc.ns << "_"
             << underscoreify(ifc.name)
             << "_set_" << ifv.name << endl;
        fout << "        func.argtypes=[ctypes.c_void_p,ctypes.c_"
             << ifv.ift.name << "]" << endl;
        fout << "        func(self._ptr,value)" << endl;
        fout << "        return" << endl;
      } else {
        fout << "    def set_" << ifv.name << "(self,value):" << endl;
        fout << "        \"\"\"" << endl;
        fout << "        Setter function for " << ifc.name << "::"
             << ifv.name << " ." << endl;
        fout << "        \"\"\"" << endl;
        fout << "        func=self._link." << dll_name << "." << ifc.ns << "_"
             << underscoreify(ifc.name) << "_set_" << ifv.name << endl;
        fout << "        func.argtypes=[ctypes.c_void_p,"
             << "ctypes.c_void_p]" << endl;
        fout << "        func(self._ptr,value._ptr)" << endl;
        fout << "        return" << endl;
      }
      fout << endl;
    }

    // Define methods
    for(size_t j=0;j<ifc.methods.size();j++) {

      if_func &iff=ifc.methods[j];

      // Function header
      if (iff.name=="index_operator") {
        fout << "    def __getitem__(self";
      } else {
        fout << "    def " << iff.name << "(self";
      }
      if (iff.args.size()>0) {
        fout << ",";
      }
      for(size_t k=0;k<iff.args.size();k++) {
        fout << iff.args[k].name;
        if (k!=iff.args.size()-1) {
          fout << ",";
        }
      }
      fout << "):" << endl;

      // Comment
      fout << "        \"\"\"" << endl;
      
      //fout << "        Wrapper for " << ifc.name << "::"
      //<< iff.name << "() ." << endl;
      if (iff.args.size()>0) {
        fout << "        | Parameters:" << endl;
      }
      for(size_t k=0;k<iff.args.size();k++) {
        if (iff.args[k].ift.name=="bool" ||
            iff.args[k].ift.name=="int" ||
            iff.args[k].ift.name=="size_t" ||
            iff.args[k].ift.name=="double") {
          fout << "        | *" << iff.args[k].name
               << "*: ``" << iff.args[k].ift.name << "``" << endl;
        } else if (iff.args[k].ift.suffix=="&") {
          fout << "        | *" << iff.args[k].name
               << "*: :ref:`" << iff.args[k].ift.name << "` object"
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
             << tmps << "`" << endl;
      } else if (iff.ret.name=="std::string") {
        fout << "        | Returns: python bytes object"
             << endl;
      } else if (iff.ret.name!="void") {
        fout << "        | Returns: ``ctypes.c_"
             << iff.ret.name << "`` object" << endl;
      }
      
      fout << "        \"\"\"" << endl;

      // Perform necessary conversions
      for(size_t k=0;k<iff.args.size();k++) {
        if (iff.args[k].ift.name=="std::string") {
          fout << "        " << iff.args[k].name
               << "_=ctypes.c_char_p(force_bytes("
               << iff.args[k].name << "))" << endl;
        }
      }

      // Ctypes interface for function
      fout << "        func=self._link." << dll_name << "." << ifc.ns << "_"
           << underscoreify(ifc.name) << "_"
           << iff.name << endl;
      if ((iff.ret.name=="vector<double>" ||
           iff.ret.name=="std::vector<double>") &&
          iff.ret.suffix=="&") {
        fout << "        n_=ctypes.c_int(0)" << endl;
        fout << "        ptr_=ctypes.POINTER(ctypes.c_double)()" << endl;
      } else if (iff.ret.prefix.find("shared_ptr")!=std::string::npos ||
                 iff.ret.prefix.find("std::shared_ptr")!=std::string::npos) {
        fout << "        func.restype=ctypes.c_void_p" << endl;
      } else if (iff.ret.name!="void") {
        if (iff.ret.name=="std::string") {
          fout << "        func.restype=ctypes.c_char_p" << endl;
        } else if (iff.ret.suffix=="&") {
          fout << "        func.restype=ctypes.c_void_p" << endl;
        } else {
          fout << "        func.restype=ctypes.c_" << iff.ret.name << endl;
        }
      }
      fout << "        func.argtypes=[ctypes.c_void_p";
      for(size_t k=0;k<iff.args.size();k++) {
        if (iff.args[k].ift.suffix=="&") {
          fout << ",ctypes.c_void_p";
        } else if (iff.args[k].ift.name=="std::string") {
          fout << ",ctypes.c_char_p";
        } else {
          fout << ",ctypes.c_" << iff.args[k].ift.name;
        }
      }
      if ((iff.ret.name=="vector<double>" ||
           iff.ret.name=="std::vector<double>") &&
          iff.ret.suffix=="&") {
        fout << ",ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),";
        fout << "ctypes.POINTER(ctypes.c_int)";
      }
      fout << "]" << endl;

      // Call C++ wrapper function
      if ((iff.ret.name=="vector<double>" ||
           iff.ret.name=="std::vector<double>") &&
          iff.ret.suffix=="&") {
        fout << "        func(self._ptr";
        for(size_t k=0;k<iff.args.size();k++) {
          if (iff.args[k].ift.suffix=="&") {
            fout << "," << iff.args[k].name << "._ptr";
          } else if (iff.args[k].ift.name=="std::string") {
            fout << "," << iff.args[k].name << "_";
          } else {
            fout << "," << iff.args[k].name;
          }
        }
        fout << ",ctypes.byref(ptr_),ctypes.byref(n_))" << endl;;
        
      } else if (iff.ret.prefix.find("shared_ptr")!=std::string::npos ||
                 iff.ret.prefix.find("std::shared_ptr")!=std::string::npos) {
        
        size_t len=iff.ret.name.length();
        std::string tmps=iff.ret.name;
        // Manually remove '<>' from the typename if necessary
        if (len>2 && iff.ret.name[len-2]=='<' &&
            iff.ret.name[len-1]=='>') {
          tmps=iff.ret.name.substr(0,len-2);
        }
        fout << "        sp=shared_ptr_"+tmps+
          "(self._link,func(self._ptr))" << endl;
        
      } else {
        
        if (iff.ret.name=="void") {
          fout << "        func(self._ptr";
        } else {
          fout << "        ret=func(self._ptr";
        }
        for(size_t k=0;k<iff.args.size();k++) {
          if (iff.args[k].ift.suffix=="&") {
            fout << "," << iff.args[k].name << "._ptr";
          } else if (iff.args[k].ift.name=="std::string") {
            fout << "," << iff.args[k].name << "_";
          } else {
            fout << "," << iff.args[k].name;
          }
        }
        fout << ")" << endl;
      }
      
      if ((iff.ret.name=="vector<double>" ||
           iff.ret.name=="std::vector<double>") &&
          iff.ret.suffix=="&") {
        fout << "        ret=numpy.ctypeslib.as_array(ptr_,shape=(n_.value,))"
             << endl;
        fout << "        return ret" << endl;
      } else if (iff.ret.suffix=="&") {
        std::string tmps=iff.ret.name;
        size_t len=iff.ret.name.length();
        // Manually remove '<>' from the typename if necessary
        if (len>2 && iff.ret.name[len-2]=='<' &&
            iff.ret.name[len-1]=='>') {
          tmps=iff.ret.name.substr(0,len-2);
        }
        fout << "        ret2=" 
             << underscoreify(tmps) << "(self._link,ret)"
             << endl;
        fout << "        return ret2" << endl;
      } else {
        // Return
        if (iff.ret.prefix.find("shared_ptr")!=std::string::npos ||
            iff.ret.prefix.find("std::shared_ptr")!=std::string::npos) {
          fout << "        return sp" << endl;
        } else if (iff.ret.name=="void") {
          fout << "        return" << endl;
        } else {
          fout << "        return ret" << endl;
        }
      }
      fout << endl;
      
    }
    
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
    
    fout << endl;
    
    // Initialize pointer
    fout << "    _s_ptr=0" << endl;
    fout << "    _ptr=0" << endl;
    fout << "    _link=0" << endl;
    fout << "    _owner=True" << endl;
    fout << endl;
    
    // Define __init__() function
    fout << "    def __init__(self,link,shared_ptr):" << endl;
    fout << "        \"\"\"" << endl;
    fout << "        Init function for shared_ptr_"
         << ifsp.name << " ." << endl;
    fout << "        \"\"\"" << endl;
    fout << endl;
    fout << "        self._link=link" << endl;
    fout << "        self._s_ptr=shared_ptr" << endl;
    fout << endl;
    fout << "        f=self._link." << dll_name << "." << ifsp.ns << "_shared_ptr_"
         << underscoreify(ifsp.name) << "_ptr" << endl;
    fout << "        f.argtypes=[ctypes.c_void_p]" << endl;
    fout << "        f.restype=ctypes.c_void_p" << endl;
    fout << "        self._ptr=f(self._s_ptr)" << endl;
    fout << "        return" << endl;
    fout << endl;
    
    // Define __del__() function
    fout << "    def __del__(self):" << endl;
    fout << "        \"\"\"" << endl;
    fout << "        Delete function for shared_ptr_"
         << ifsp.name << " ." << endl;
    fout << "        \"\"\"" << endl;
    fout << endl;
    fout << "        f=self._link." << dll_name << "." << ifsp.ns << "_free_shared_ptr_"
         << underscoreify(ifsp.name) << endl;
    fout << "        f.argtypes=[ctypes.c_void_p]" << endl;
    fout << "        f(self._ptr)" << endl;
    fout << "        return" << endl;
    fout << endl;

  }
  
  // Define functions
  for(size_t j=0;j<functions.size();j++) {
    
    if_func &iff=functions[j];
    
    // Function header
    fout << "def " << iff.name << "(link,";
    for(size_t k=0;k<iff.args.size();k++) {
      fout << iff.args[k].name;
      if (k!=iff.args.size()-1) {
        fout << ",";
      }
    }
    fout << "):" << endl;
    
    // Comment
    fout << "    \"\"\"" << endl;
    fout << "    Wrapper for " 
         << iff.name << "() ." << endl;
    fout << "    \"\"\"" << endl;
    
    // Perform necessary conversions
    for(size_t k=0;k<iff.args.size();k++) {
      if (iff.args[k].ift.name=="std::string") {
        fout << "    " << iff.args[k].name << "_=ctypes.c_char_p(force_bytes("
             << iff.args[k].name << "))" << endl;
      }
    }
      
    // Ctypes interface for function
    fout << "    func=link." << dll_name << "." << iff.ns << "_" << iff.name
         << "_wrapper" << endl;
    if (iff.ret.name!="void") {
      if (iff.ret.name=="std::string") {
        fout << "    func.restype=ctypes.c_char_p" << endl;
      } else {
        fout << "    func.restype=ctypes.c_" << iff.ret.name << endl;
      }
    }
    fout << "    func.argtypes=[";
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
    
    // Call C++ wrapper function
    if (iff.ret.name=="void") {
      fout << "    func(";
    } else {
      fout << "    ret=func(";
    }
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
    fout << ")" << endl;
    
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
  
  for(size_t i=0;i<classes.size();i++) {

    if_class &ifc=classes[i];

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
    
    if (ifc.py_name=="") {
      fout2 << ".. autoclass:: o2sclpy." << ifc.name << endl;
    } else {
      fout2 << ".. autoclass:: o2sclpy." << ifc.py_name << endl;
    }
    fout2 << "        :members:" << endl;
    fout2 << "        :undoc-members:\n" << endl;
  }

  for(size_t i=0;i<sps.size();i++) {
    if_shared_ptr &ifsp=sps[i];
    
    if (ifsp.py_name!="") {
      size_t len=ifsp.name.length();
      std::string tmps=ifsp.name;
      // Manually remove '<>' from the typename if necessary
      if (len>2 && ifsp.name[len-2]=='<' &&
          ifsp.name[len-1]=='>') {
        tmps=ifsp.name.substr(0,len-2);
      }
      fout2 << ".. autoclass:: o2sclpy.shared_ptr_" << tmps << endl;
      fout2 << "        :members:" << endl;
      fout2 << "        :undoc-members:\n" << endl;
    } else {
      /*
      size_t len=ifsp.py_name.length();
      std::string tmps=ifsp.py_name;
      // Manually remove '<>' from the typename if necessary
      if (len>2 && ifsp.py_name[len-2]=='<' &&
          ifsp.py_name[len-1]=='>') {
        tmps=ifsp.py_name.substr(0,len-2);
      }
      */
      fout2 << ".. autoclass:: o2sclpy.shared_ptr_" << ifsp.py_name << endl;
      fout2 << "        :members:" << endl;
      fout2 << "        :undoc-members:\n" << endl;
    }
    
  }
  
  fout2.close();
  
  return 0;
}
