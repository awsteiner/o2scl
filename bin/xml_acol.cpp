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
#include <iostream>
#include <fstream>
#include <string>

#include <o2scl/string_conv.h>
#include <o2scl/hdf_io.h>
#include <o2scl/hdf_file.h>
#include <o2scl/xml.h>
#include <o2scl/acolm.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_acol;
using namespace o2scl_hdf;

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);
  
#ifdef O2SCL_PUGIXML

  vector<string> data;
  
  acol_manager am;
  am.run(0,0,false);

  vector<string> cli_list={"run","shell","quit","exit",
    "set","get","no-intro","license","warranty","alias"};
  vector<string> clist=am.cl->get_option_list();

  for(size_t i=0;i<clist.size();i++) {

    std::string cmd_name=clist[i], fn_name="comm_";
    for(size_t k=0;k<cmd_name.length();k++) {
      if (cmd_name[k]=='-') {
        fn_name+='_';
      } else {
        fn_name+=cmd_name[k];
      }
    }
    cout << cmd_name << " " << fn_name << endl;
      
    if (std::find(cli_list.begin(),cli_list.end(),cmd_name)==
        cli_list.end()) {
      
      std::string doc_fn=o2scl_settings.get_doc_dir()+
        "xml/classo2scl__acol_1_1acol__manager.xml";
      
      pugi::xml_document doc1, doc2;
      
      pugi::xml_node n1=doxygen_xml_member_get(doc_fn,"acol_manager",fn_name,
                                               "briefdescription",doc1);
      if (n1!=0) {
        pugi::xml_node n2=doxygen_xml_member_get(doc_fn,"acol_manager",fn_name,
                                                 "detaileddescription",doc2);
        if (n2!=0) {
          cout << "dxg: " << n2.name() << endl;
          ostream_walker walker;
          n2.traverse(walker);
        }
      }
      
    }

  }
  
#endif
  
  return 0;
}
