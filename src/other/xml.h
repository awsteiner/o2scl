/*
  -------------------------------------------------------------------

  Copyright (C) 2022, Andrew W. Steiner

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
/** \file xml.h
    \brief File defining \ref o2scl::xml
*/
#ifndef O2SCL_XML_H
#define O2SCL_XML_H

#include <string>
#include <vector>

#include <o2scl/err_hnd.h>

#ifdef O2SCL_PUGIXML
#include "pugixml.hpp"
#endif

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

#ifdef O2SCL_PUGIXML

  /** \brief A base class for PugiXML walkers
   */
  class walker_base  : public pugi::xml_tree_walker {
    
  public:

    int last_depth;

    int verbose;
    
    std::vector<std::string> names;

    virtual bool begin(pugi::xml_node &node) {
      last_depth=-1;
      verbose=0;
      names.clear();
      return true;
    }
      
    walker_base() {
      last_depth=-1;
      verbose=0;
    }
    
  };
  
  /** \brief A PugiXML walker which outputs the full contents of the 
      node (no attributes yet) to an ostream
   */
  class ostream_walker : public walker_base {
    
  public:
    
    std::ostream *outs;

    ostream_walker() {
      outs=&std::cout;
    }

    virtual bool for_each(pugi::xml_node &node) {

      if (verbose>0) {
        (*outs) << " d: "
                  << depth() << " l: " << last_depth << std::endl;
      }
      
      if (depth()<last_depth) {
        int n=last_depth-depth();
        for(int i=0;i<n;i++) {
          for (int j = 0; j < depth()-i+n-1; j++) {
            (*outs) << "  ";
          }
          (*outs) << "</" << names[names.size()-1] << ">"
                    << std::endl;
          names.pop_back();
        }
      }
      
      if (((std::string)node.name()).length()>0) {
        names.push_back(node.name());
      }
      
      for (int i = 0; i < depth(); i++) {
        (*outs) << "  ";
      }
      
      if (((std::string)node.name()).length()>0) {
        (*outs) << "<" << node.name()
                  << ">" << node.value();
      } else {
        (*outs) << node.value();
      }
      if (verbose>0) {
        for(size_t k=0;k<names.size();k++) {
          (*outs) << "." << names[k] << ". ";
        }
      }
      (*outs) << std::endl;
      
      last_depth=depth();
      
      //names=node.name();
      return true;
    }

    virtual bool end(pugi::xml_node &node) {
      int n=last_depth;
      for(int i=0;i<n;i++) {
        if (names.size()>0) {
          for (int j = 0; j < depth()-i+n; j++) {
            (*outs) << "  ";
          }
          (*outs) << "</" << names[names.size()-1] << ">"
                  << std::endl;
          names.pop_back();
        }
      }
      return true;
    }
    
  };
  
  /** \brief A PugiXML walker which outputs the full contents of the 
      node (no attributes yet) to a vector<string> object
   */
  class vec_string_walker : public walker_base {
    
  public:

    /** \brief The traversal output
     */
    std::vector<std::string> output;

    /** \brief Before traversing
     */
    virtual bool begin(pugi::xml_node &node) {
      output.clear();
      walker_base::begin(node);
      return true;
    }

    /** \brief For each subnode
     */
    virtual bool for_each(pugi::xml_node &node) {

      std::string stmp;
      
      if (verbose>0) {
        std::cout << " d: "
                  << depth() << " l: " << last_depth << std::endl;
      }
      
      if (depth()<last_depth) {
        int n=last_depth-depth();
        for(int i=0;i<n;i++) {
          for (int j = 0; j < depth()-i+n-1; j++) {
            stmp+=((std::string)"  ");
          }
          stmp+=((std::string)"</")+names[names.size()-1]+">";
          output.push_back(stmp);
          stmp.clear();
          names.pop_back();
        }
      }
      
      if (((std::string)node.name()).length()>0) {
        names.push_back(node.name());
      }
      
      for (int i = 0; i < depth(); i++) {
        stmp+=((std::string)"  ");
      }
      
      if (((std::string)node.name()).length()>0) {
        stmp+=((std::string)"<")+node.name()+">"+node.value();
      } else {
        stmp+=node.value();
      }
      if (verbose>0) {
        for(size_t k=0;k<names.size();k++) {
          std::cout << "." << names[k] << ". ";
        }
      }
      output.push_back(stmp);
      stmp.clear();
      
      last_depth=depth();
      
      //names=node.name();
      return true;
    }

    /** \brief After traversing
     */
    virtual bool end(pugi::xml_node &node) {
      int n=last_depth;
      for(int i=0;i<n;i++) {
        std::string stmp;
        if (names.size()>0) {
          for (int j = 0; j < depth()-i+n; j++) {
            stmp+=((std::string)"  ");
          }
          stmp+=((std::string)"</")+names[names.size()-1]+">";
          output.push_back(stmp);
          stmp.clear();
          names.pop_back();
        }
      }
      return true;
    }
    
  };
  
  /** \brief Extract XML node named \c node_name in the doxygen
      documentation for a global function named \c func_name from a
      file named \c fname
  */
  pugi::xml_node doxygen_xml_get
  (std::string fname, std::string func_name,
   std::string node_name, pugi::xml_document &doc,
   int verbose=0);
  
  /** \brief Extract XML node named \c node_name in the doxygen
      documentation for a member function named \c func_name from a
      class named \c class_name from a file named \c fname
  */
  pugi::xml_node doxygen_xml_member_get
  (std::string fname, std::string class_name, std::string func_name, 
   std::string node_name, pugi::xml_document &doc, int verbose=0);
   
#endif
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif

