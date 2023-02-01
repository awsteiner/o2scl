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
#include <climits>
#include <o2scl/misc.h>
#include <o2scl/test_mgr.h>
#include <o2scl/xml.h>

using namespace std;
using namespace o2scl;

int main(void) {

  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(2);

#ifdef O2SCL_PUGIXML
  
  pugi::xml_document doc1, doc2;
  
  std::string doc_fn=((string)"../../doc/o2scl/")+
    "xml/namespaceo2scl.xml";
  
  pugi::xml_node n=doxygen_xml_get(doc_fn,"glob_wrapper",
                                   "detaileddescription",doc1);
  cout << "dxg: " << n.name() << endl;
  ostream_walker walker;
  n.traverse(walker);

  /*
    AWS, 7/21/22: this file is no longer part of the installation so
    this test no longer works
    
    std::string doc2_fn=o2scl_settings.get_doc_dir()+
    "xml/classo2scl_1_1gen__test__number.xml";
    
    pugi::xml_node n2=doxygen_xml_member_get
    (doc2_fn,"gen_test_number","set_radix","detaileddescription",doc2);
    
    cout << "dxmg: " << n2.name() << endl;
    n2.traverse(walker);
  */
  
#endif
  
  t.report();

  return 0;
}

