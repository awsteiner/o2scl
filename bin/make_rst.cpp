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
#include <vector>
#include <set>
#include <map>

#ifdef O2SCL_PUGIXML
#include "pugixml.hpp"
#endif

using namespace std;

std::string remove_template_params(std::string in) {
  std::string stmp;
  int n_bracket=0;
  for(size_t ij=0;ij<in.length();ij++) {
    if (in[ij]=='<') {
      n_bracket++;
    } else if (in[ij]=='>') {
      n_bracket--;
    } else if (n_bracket==0) {
      stmp+=in[ij];
    }
  }
  return stmp;
}

/** \brief Convert <,>, ,: to _
 */
std::string underscoreify(std::string s) {
  std::string s2=s;
  for(size_t i=0;i<s2.length();i++) {
    if (s2[i]==' ' || s2[i]=='<' || s2[i]==':' || s2[i]=='>') {
      s2[i]='_';
    }
  }
  while (s2.find("__")!=std::string::npos) {
    s2.replace(s2.find("__"),2,"_");
  }
  return s2;
}

/** \brief Inside string \c s, extract the element inside tag \c tag
    and store the element in \c result
*/
int xml_get_tag_element(string s, string tag, string &result) {
  // Remove enclosing <> in tag
  if (tag[0]=='<') {
    if (tag[tag.length()-1]!='>') {
      return 1;
    }
    tag=tag.substr(1,tag.length()-2);
  }
  size_t tag_len=tag.length();
  string tag_start=((string)"<")+tag+">";
  string tag_end=((string)"</")+tag+">";

  // Exit early if the tag is not found
  size_t istart=s.find(tag_start);
  if (istart==std::string::npos) {
    return 2;
  }
  size_t iend=s.find(tag_end);
  if (iend==std::string::npos) {
    return 3;
  }
  if (istart+tag_len+2>iend) return 4;
  
  // Make sure the <name> tag is also present in the string somewhere,
  // return early if it is not
  size_t loc1=s.find("<name>");
  size_t loc2=s.find("</name>");
  size_t loct=loc1+6;
  if (loc2<loct) {
    return 5;
  }
  
  // Finally, return the result, adjusting by 2 for the < and >
  size_t len=iend-istart-tag_len-2;
  result=s.substr(istart+tag_len+2,len);
  return 0;
}

int main(int argc, char *argv[]) {

#ifdef O2SCL_PUGIXML
  
  if (argc<2) {
    cerr << "Requires input file." << endl;
    exit(-1);
  }

  std::string in_dir=argv[1];
  std::string in_file=in_dir+"/index.xml";
  cout << "Using input file " << in_file << endl;

  // -----------------------------------------------------------------
  
  typedef struct generic_s {
    std::string ns;
    std::string name;
  } generic;

  typedef struct function_s {
    std::string ns;
    std::string name;
    std::vector<std::string> args;
  } function;

  std::vector<generic> class_list, function_list;
  
  std::vector<function> overloaded_list;
  
  // -----------------------------------------------------------------
  // Step 1: Parse the doxygen XML to find all the classes and
  // functions

  int xml_verbose=0;
    
  { 
    pugi::xml_document doc;
    
    pugi::xml_parse_result result=doc.load_file(in_file.c_str());
    if (!result) {
      cout << result.description() << endl;
      cout << "Failed to read file." << endl;
      exit(-1);
    }
    
    pugi::xml_node dindex=doc.first_child();
    
    if (xml_verbose>0) {
      std::cout << "1: " << dindex.name() << endl;
    }
    int i=0;
    
    for (pugi::xml_node_iterator it=dindex.begin();it!=dindex.end();++it) {
      if (xml_verbose>0) {
        cout << i << " " << it->name() << " "
             << it->child_value("name") << endl;
        cout << "  kind: " << it->attribute("kind").value() << std::endl;
      }
      
      if (((string)it->attribute("kind").value())==((string)"class")) {
        generic o;
        o.name=it->child_value("name");
        size_t loc=o.name.find_last_of(':');
        if (loc!=std::string::npos && loc>1) {
          o.ns=o.name.substr(0,loc-1);
          o.name=o.name.substr(loc+1,o.name.length()-loc);
        }
        class_list.push_back(o);
      }
      
      if (((string)it->attribute("kind").value())==((string)"namespace")) {
        
        int j=0;
        
        for (pugi::xml_node_iterator it2=it->begin();
             it2!=it->end();++it2) {
          
          if (xml_verbose>1) {
            cout << j << " " << it2->name() << " "
                 << it2->child_value("name") << endl;
            cout << "  kind: " << it2->attribute("kind").value()
                 << std::endl;
          }
          
          if (((string)it2->attribute("kind").value())==((string)"function")) {
            generic o;
            o.ns=it->child_value("name");
            o.name=it2->child_value("name");

            // Remove template parameters from function name
            o.name=remove_template_params(o.name);
            
            //
            bool found_in_list=false;
            for(size_t k=0;k<function_list.size();k++) {
              if (function_list[k].ns==o.ns &&
                  function_list[k].name==o.name) {
                found_in_list=true;
              }
            }
            if (found_in_list) {
              function f;
              f.ns=o.ns;
              f.name=o.name;
              bool found_in_dups=false;
              for(size_t k=0;k<overloaded_list.size();k++) {
                if (overloaded_list[k].ns==o.ns &&
                    overloaded_list[k].name==o.name) {
                  found_in_dups=true;
                }
              }
              if (found_in_dups==false) {
                overloaded_list.push_back(f);
              }
            } else {
              function_list.push_back(o);
            }
          }
          j++;
        }  
      }
      i++;
    }

    // End of reading index.xml
  }

  vector<std::string> ns_list;

  for(size_t k=0;k<overloaded_list.size();k++) {
    
    bool found_in_ns=false;
    for(size_t i=0;i<ns_list.size();i++) {
      if (ns_list[i]==overloaded_list[k].ns) {
        found_in_ns=true;
      }
    }
    if (found_in_ns==false) {
      ns_list.push_back(overloaded_list[k].ns);
    }

  }

  for(size_t k=0;k<ns_list.size();k++) {

    // Doxygen doubles namespaces for the namespace.xml file
    // so we do that here
    std::string stmp=ns_list[k];
    for(int ik=0;ik<((int)stmp.length());ik++) {
      if (stmp[ik]=='_') {
        stmp=stmp.substr(0,ik)+"__"+stmp.substr(ik+1,stmp.length()-ik-1);
        ik+=2;
        if (ik==100) {
          cout << "Failed to rename namespace file: " << stmp << endl;
          exit(-1);
        }
      }
    }
    
    if (ns_list[k]!="std") {
      
      pugi::xml_document ns_doc;
      
      std::string ns_file=in_dir+"/namespace"+stmp+".xml";

      cout << "Reading namespace file: " << ns_file << endl;
      pugi::xml_parse_result result=ns_doc.load_file(ns_file.c_str());
      if (!result) {
        cout << result.description() << endl;
        cout << "Failed to read namespace file " << ns_file << endl;
        exit(-1);
      }

      // Parse through <doxygen><compounddef>
      pugi::xml_node dindex=ns_doc.first_child().first_child();
      
      for (pugi::xml_node_iterator it=dindex.begin();it!=dindex.end();++it) {

        // Parse through the namespace
        
        //std::cout << "2: " << it->name() << std::endl;

        // The namespace name is in a <compoundname> object, classes
        // are in <innerclass> objects, and some functions are 
        // stored in sections, <sectiondef> objects
        
        if (it->name()==((string)"sectiondef")) {

          cout << "Section: " << it->child("header").child_value() << endl;
          
          // In each section, look for a <memberdef> object with a
          // kind attribute of "function"
          
          for (pugi::xml_node_iterator it2=it->begin();
               it2!=it->end();++it2) {
            
            //std::cout << "3: " << it2->name() << " "
            //<< it2->attribute("kind").value() << std::endl;
            
            if (it2->name()==((string)"memberdef") &&
                it2->attribute("kind").value()==((string)"function")) {
              
              //cout << "  " << it2->child("name").child_value() << endl;
              //cout << "  " << it2->child("argsstring").child_value() << endl;

              // If we found a function, see if its in overloaded_list
              // so we can set the argstrings
              
              bool found=false;
              std::string func_name=it2->child("name").child_value();
              func_name=remove_template_params(func_name);
              
              for(size_t j=0;j<overloaded_list.size() && found==false;j++) {
                if (overloaded_list[j].ns==ns_list[k] &&
                    overloaded_list[j].name==func_name) {
                  
                  found=true;

                  pugi::xml_node nt=it2->child("argsstring");
                  overloaded_list[j].args.push_back(nt.child_value());

                  /*
                    cout << "Found overloaded function in namespace "
                    << overloaded_list[j].ns << " with name "
                    << overloaded_list[j].name << endl;
                    for(size_t ik=0;ik<overloaded_list[j].args.size();ik++) {
                    cout << "  " << ik << " "
                    << overloaded_list[j].args[ik].substr(0,74)
                    << endl;
                    }
                    cout << endl;
                  */
                }
              }

              /*
                if (found==false) {
                cout << "This function is not in overloaded list." << endl;
                }
              */
              
            }
          }
          
        }
        
      }
      
    }
    
  }

  for(size_t i=0;i<overloaded_list.size();i++) {
    bool found=false;
    size_t jfound;
    for(size_t j=0;j<function_list.size();j++) {
      if (function_list[j].ns==overloaded_list[i].ns &&
          function_list[j].name==overloaded_list[i].name) {
        found=true;
        jfound=j;
      }
    }
    if (found==true) {
      std::vector<generic> function_list_new;
      for(size_t k=0;k<function_list.size();k++) {
        if (k!=jfound) {
          generic g;
          g.ns=function_list[k].ns;
          g.name=function_list[k].name;
          function_list_new.push_back(g);
        }
      }
      function_list=function_list_new;
    }
  }
  
  if (true || xml_verbose>1) {
    cout << "Class list: " << endl;
    for(size_t i=0;i<class_list.size();i++) {
      cout << i << " " << class_list[i].ns << " "
           << class_list[i].name << endl;
    }
    cout << endl;
    cout << "Function list: " << endl;
    for(size_t i=0;i<function_list.size();i++) {
      cout << i << " " << function_list[i].ns << " "
           << function_list[i].name << endl;
    }
    cout << endl;
    cout << "Overloaded list: " << endl;
    for(size_t i=0;i<overloaded_list.size();i++) {
      cout << i << " " << overloaded_list[i].ns << " "
           << overloaded_list[i].name << endl;
      for(size_t ik=0;ik<overloaded_list[i].args.size();ik++) {
        cout << "  " << ik << " ";
        if (overloaded_list[i].args[ik].length()<70) {
          cout << overloaded_list[i].args[ik];
        } else {
          cout << overloaded_list[i].args[ik].substr(0,70) << "..";
        }
        cout << endl;
      }
    }
    cout << endl;
    cout << "Namespace list: " << endl;
    for(size_t i=0;i<ns_list.size();i++) {
      cout << i << " " << ns_list[i] << endl;
    }
    cout << endl;
  }

  // -------------------------------------------------------
  
  cout << "---------------------------------------------------" << endl;
  cout << "Creating rst files for classes:" << endl;
  cout << endl;
  
  // Proceed through the list
  for(size_t j=0;j<class_list.size();j++) {
    
    // Open the rst file
    string fname_out="class/";
    fname_out+=class_list[j].name+".rst";
    if (class_list[j].name.length()<4 ||
        class_list[j].name[0]!='o' || class_list[j].name[1]!='p' ||
        class_list[j].name[2]!='e' || class_list[j].name[3]!='r') {
      fname_out=underscoreify(fname_out);
    }
    ofstream fout(fname_out);

    // Title
    string head="Class "+class_list[j].name;
    // If we're in "class mode" and the namespace is not empty,
    // then add the namespace to the header
    head+=" ("+class_list[j].ns+")";

    fout << head << endl;
    for(size_t i=0;i<head.length();i++) {
      fout << "=";
    }
    fout << endl;
    fout << endl;
    
    fout << ":ref:`O2scl <o2scl>` : :ref:`Class List`\n" << endl;
    fout << ".. _" << class_list[j].name << ":\n" << endl;
    
    // Output the class or function directive
    if (class_list[j].ns.length()>0) {
      fout << ".. doxygenclass:: " << class_list[j].ns << "::"
           << class_list[j].name << endl;
    } else {
      fout << ".. doxygenclass:: " << class_list[j].name << endl;
    }
    
    // Close the file
    fout.close();

    cout << "Wrote class file: " << fname_out << endl;
  }
  cout << endl;

  cout << "---------------------------------------------------" << endl;
  cout << "Creating rst files for functions not overloaded:" << endl;
  cout << endl;
  
  // Proceed through the list
  for(size_t j=0;j<function_list.size();j++) {
    
    // Open the rst file
    string fname_out="function/";
    fname_out+=function_list[j].name+".rst";
    if (function_list[j].name.length()<4 ||
        function_list[j].name[0]!='o' || function_list[j].name[1]!='p' ||
        function_list[j].name[2]!='e' || function_list[j].name[3]!='r') {
      fname_out=underscoreify(fname_out);
    }
    ofstream fout(fname_out);

    // Title
    string head="Function "+function_list[j].name;
    // If we're in "function mode" and the namespace is not empty,
    // then add the namespace to the header
    head+=" ("+function_list[j].ns+")";
	
    fout << head << endl;
    for(size_t i=0;i<head.length();i++) {
      fout << "=";
    }
    fout << endl;
    fout << endl;
    
    fout << ":ref:`O2scl <o2scl>` : :ref:`Function List`\n" << endl;
    
    // Output the function or function directive
    fout << ".. doxygenfunction:: " << function_list[j].ns << "::"
         << function_list[j].name << endl;
    
    // Close the file
    fout.close();

    cout << "Wrote function file: " << fname_out << endl;
    
  }
  cout << endl;

  cout << "---------------------------------------------------" << endl;
  cout << "Creating rst files for overloaded functions:" << endl;
  cout << endl;
  
  // Proceed through the list
  for(size_t j=0;j<overloaded_list.size();j++) {
    
    // Open the rst file
    string fname_out="function/";
    fname_out+=overloaded_list[j].name+".rst";
    if (overloaded_list[j].name.length()<4 ||
        overloaded_list[j].name[0]!='o' || overloaded_list[j].name[1]!='p' ||
        overloaded_list[j].name[2]!='e' || overloaded_list[j].name[3]!='r') {
      fname_out=underscoreify(fname_out);
    }
    ofstream fout(fname_out);

    // Title
    string head="Function "+overloaded_list[j].name;
    // If we're in "function mode" and the namespace is not empty,
    // then add the namespace to the header
    head+=" ("+overloaded_list[j].ns+")";
	
    fout << head << endl;
    for(size_t i=0;i<head.length();i++) {
      fout << "=";
    }
    fout << endl;
    fout << endl;
    
    fout << ":ref:`O2scl <o2scl>` : :ref:`Function List`\n" << endl;
    
    for(size_t k=0;k<overloaded_list[j].args.size();k++) {
      // Output the function or function directive
      fout << ".. doxygenfunction:: " << overloaded_list[j].ns << "::"
           << overloaded_list[j].name
           << overloaded_list[j].args[k] << "\n" << endl;
    }
    
    // Close the file
    fout.close();

    cout << "Wrote function file: " << fname_out << endl;
    
  }
  cout << endl;

#endif
  
  return 0;
}
