/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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

#include <o2scl/xml.h>
#include <o2scl/set_pugixml.h>

#ifdef O2SCL_SET_PUGIXML
    
pugi::xml_node o2scl::doxygen_xml_get(std::string fname,
                                      std::string func_name,
                                      std::string node_name,
                                      pugi::xml_document &doc,
                                      int verbose) {

  if (verbose>1) {
    std::cout << "Looking for " << func_name << " in file "
              << fname << std::endl;
  }
    
  pugi::xml_parse_result result=doc.load_file(fname.c_str());
  if (!result) {
    std::cout << result.description() << std::endl;
    std::cout << "Failed to read namespace file " << fname << std::endl;
    exit(-1);
  }
    
  // Parse through <doxygen><compounddef>
  pugi::xml_node dindex=doc.first_child().first_child();
    
  for (pugi::xml_node_iterator it=dindex.begin();it!=dindex.end();++it) {
      
    // Parse through the namespace

    if (verbose>1) {
      std::cout << "2: " << it->name() << std::endl;
    }
      
    // The namespace name is in a <compoundname> object, classes
    // are in <innerclass> objects, and some functions are 
    // stored in sections, <sectiondef> objects
      
    if (it->name()==((std::string)"sectiondef")) {
        
      std::string section_name=it->child("header").child_value();
      if (verbose>0) {
        if (section_name.length()==0) {
          std::cout << "Section: <no name>" << std::endl;
        } else {
          std::cout << "Section: "
                    << it->child("header").child_value() << std::endl;
        }
      }
        
      // In each section, look for a <memberdef> object with a
      // kind attribute of "function"
        
      for (pugi::xml_node_iterator it2=it->begin();
           it2!=it->end();++it2) {

        if (verbose>1) {
          std::cout << "3: " << it2->name() << " "
                    << it2->attribute("kind").value() << std::endl;
        }
          
        if (it2->name()==((std::string)"memberdef") &&
            it2->attribute("kind").value()==((std::string)"function")) {

          if (verbose>1) {
            std::cout << "Found function named: "
                      << it2->child("name").child_value() << std::endl;
          }
            
          // If we found a function, see if its in overloaded_list
          // so we can set the argstrings
            
          std::string func_name2=it2->child("name").child_value();
          
          //std::string tlate_parms;
          //separate_template_params(func_name2,func_name2,
          //tlate_parms);

          if (func_name==func_name2) {
            return it2->child(node_name.c_str());
          }
            
          // End of if statement we're dealing with a function
        }

        // End of loop over XML nodes in this section
      }

      // End of 'if (it->name()==((std::string)"sectiondef"))'
    }

    // Main loop over the XML nodes
  }

  pugi::xml_node empty;
  return empty;
}

pugi::xml_node o2scl::doxygen_xml_member_get
(std::string fname, std::string class_name, std::string func_name, 
 std::string node_name,
 pugi::xml_document &doc,
 int verbose) {
  
  if (verbose>1) {
    std::cout << "Looking for member " << func_name
              << " of class " << class_name << " in file "
              << fname << std::endl;
  }
    
  pugi::xml_parse_result result=doc.load_file(fname.c_str());
  if (!result) {
    std::cout << result.description() << std::endl;
    std::cout << "Failed to read namespace file " << fname << std::endl;
    exit(-1);
  }
    
  // Parse through <doxygen><compounddef>
  pugi::xml_node dindex=doc.first_child().first_child();
    
  for (pugi::xml_node_iterator it=dindex.begin();it!=dindex.end();++it) {
      
    // Parse through the namespace

    if (verbose>1) {
      std::cout << "2: " << it->name() << std::endl;
    }
      
    // The namespace name is in a <compoundname> object, classes
    // are in <innerclass> objects, and some functions are 
    // stored in sections, <sectiondef> objects
      
    if (it->name()==((std::string)"sectiondef")) {
        
      std::string section_name=it->child("header").child_value();
      if (verbose>0) {
        if (section_name.length()==0) {
          std::cout << "Section: <no name>" << std::endl;
        } else {
          std::cout << "Section: "
                    << it->child("header").child_value() << std::endl;
        }
      }
        
      // In each section, look for a <memberdef> object with a
      // kind attribute of "function" or "variable"
        
      for (pugi::xml_node_iterator it2=it->begin();
           it2!=it->end();++it2) {

        if (verbose>1) {
          std::cout << "3: " << it2->name() << " "
                    << it2->attribute("kind").value() << std::endl;
        }
          
        if (it2->name()==((std::string)"memberdef") &&
            (it2->attribute("kind").value()==((std::string)"function") ||
             it2->attribute("kind").value()==((std::string)"variable"))) {

          if (verbose>1) {
            std::cout << "Found function named: "
                      << it2->child("name").child_value() << std::endl;
          }
            
          // If we found a function, see if its in overloaded_list
          // so we can set the argstrings
            
          std::string func_name2=it2->child("name").child_value();
          //std::string tlate_parms;
          //separate_template_params(func_name2,func_name2,
          //tlate_parms);

          if (func_name==func_name2) {
            return it2->child(node_name.c_str());
          }
            
          // End of if statement we're dealing with a function
        }

        // End of loop over XML nodes in this section
      }

      // End of 'if (it->name()==((std::string)"sectiondef"))'
    }

    // Main loop over the XML nodes
  }

  pugi::xml_node empty;
  return empty;
}
  
#endif
