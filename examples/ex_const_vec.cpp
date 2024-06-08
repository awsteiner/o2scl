/*
  -------------------------------------------------------------------
  
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

  -------------------------------------------------------------------
*/

/* Example: ex_const_vec.cpp
   -------------------------------------------------------------------
   This demonstrates const-correctness in vectors

   (This is still under construction.)

*/

#include <o2scl/test_mgr.h>
#include <o2scl/ovector_tlate.h>

using namespace std;
using namespace o2scl;

int main(void) {
  
  cout.setf(ios::scientific);
  
  test_mgr t;
  // Only print something out if one of the tests fails
  t.set_output_level(1);

  {
    // Make several starting objects
    // ovector
    ovector vec(2);
    vec[0]=2.0;
    vec[1]=1.0;
    // const ovector
    const ovector cvec=vec;
    // vector_view
    ovector_view view=vec;
    // const vector_view
    const ovector_view cview=vec;
    // vector_base
    ovector_base bview=vec;
    // const vector_base
    const ovector_base cbview=vec;
    // vector_const_view
    ovector_const_view ocview=cvec;
    // const vector_const_view
    const ovector_const_view cocview=cvec;
    
    {
      // These are all deep copies

      //ovector copies of const ovectors
      ovector vec1a=cvec;
      const ovector vec2a=cvec;
      // Show we're const
      // vec2a[0]=3.0;
      
      //ovector copies of non-const ovectors
      ovector vec1b=vec;
      const ovector vec2b=vec;
      // Show we're const
      // vec2b[0]=3.0;

      // const_views of const ovector_views
      ovector vec1c=cview;
      const ovector vec2c=cview;
      // Show we're const
      // vec2c[0]=3.0;

      // const_views of non-const ovector_views
      ovector vec1d=view;
      const ovector vec2d=view;
      // Show we're const
      // vec2d[0]=3.0;

      // const_views of const ovector_base
      ovector vec1e=cbview;
      const ovector vec2e=cbview;
      // Show we're const
      // vec2e[0]=3.0;

      // const_views of non-const ovector_base
      ovector vec1f=bview;
      const ovector vec2f=bview;
      // Show we're const
      // vec2f[0]=3.0;

      // const_views of const ovector_const_views
      ovector vec1g=cocview;
      const ovector vec2g=cocview;
      // Show we're const
      // vec2g[0]=3.0;

      // const_views of non-const ovector_const_views
      ovector vec1h=ocview;
      const ovector vec2h=ocview;
      // Show we're const
      // vec2h[0]=3.0;
    }

    {
      // const_views of const ovectors
      ovector_const_view ocview1a=cvec;
      const ovector_const_view ocview2a=cvec;
      // Show we're const
      // ocview1a[0]=3.0;
      // ocview2a[0]=3.0;

      // const_views of non-const ovectors
      ovector_const_view ocview1b=vec;
      const ovector_const_view ocview2b=vec;
      // Show we're const
      // ocview1b[0]=3.0;
      // ocview2b[0]=3.0;

      // const_views of const ovector_views
      ovector_const_view ocview1c=cview;
      const ovector_const_view ocview2c=cview;
      // Show we're const
      // ocview1c[0]=3.0;
      // ocview2c[0]=3.0;

      // const_views of non-const ovector_views
      ovector_const_view ocview1d=view;
      const ovector_const_view ocview2d=view;
      // Show we're const
      // ocview1d[0]=3.0;
      // ocview2d[0]=3.0;

      // const_views of const ovector_base
      ovector_const_view ocview1e=cbview;
      const ovector_const_view ocview2e=cbview;
      // Show we're const
      // ocview1e[0]=3.0;
      // ocview2e[0]=3.0;

      // const_views of non-const ovector_base
      ovector_const_view ocview1f=bview;
      const ovector_const_view ocview2f=bview;
      // Show we're const
      // ocview1f[0]=3.0;
      // ocview2f[0]=3.0;

      // const_views of const ovector_const_views
      ovector_const_view ocview1g=cocview;
      const ovector_const_view ocview2g=cocview;
      // Show we're const
      // ocview1g[0]=3.0;
      // ocview2g[0]=3.0;

      // const_views of non-const ovector_const_views
      ovector_const_view ocview1h=ocview;
      const ovector_const_view ocview2h=ocview;
      // Show we're const
      // ocview1h[0]=3.0;
      // ocview2h[0]=3.0;
    }

    {

      // base views of const ovectors 
      // (Not allowed because ovector_base data is not const.)
      // ovector_base ocbase1a=cvec;
      // const ovector_base ocbase2a=cvec;

      // Reference is allowed, as long as it's const
      // ovector_base &ocbase1a=cvec;
      const ovector_base &ocbase2a=cvec;

      // base views of non-const ovectors
      ovector_base ocbase1b=vec;
      const ovector_base ocbase2b=vec;

      // base views of const ovector_views
      // (Not allowed to prevent non-const copying of const ovectors.)
      // ovector_base ocbase1c=cview;
      // const ovector_base ocbase2c=cview;

      // Reference is allowed, as long as it's const
      // ovector_base &ocbase1c=cview;
      const ovector_base &ocbase2c=cview;

      // base views of non-const ovector_views
      ovector_base ocbase1d=view;
      const ovector_base ocbase2d=view;
      
      // base views of const ovector_base
      // (Not allowed to prevent non-const copying of const ovectors.)
      // ovector_base ocbase1e=cbview;
      // const ovector_base ocbase2e=cbview;

      // Reference is allowed, as long as it's const
      // ovector_base &ocbase1e=cbview;
      const ovector_base &ocbase2e=cbview;

      // base views of non-const ovector_base
      ovector_base ocbase1f=bview;
      const ovector_base ocbase2f=bview;
      
      // base views of const ovector_const_views
      // (Not allowed because ovector_base data is not const.)
      // ovector_base ocbase1g=cocview;
      // const ovector_base ocbase2g=cocview;
      // ovector_base &ocbase1g=cocview;
      // const ovector_base &ocbase2g=cocview;

      // base views of non-const ovector_const_views
      // (Not allowed because ovector_base data is not const.)
      // ovector_base ocbase1h=ocview;
      // const ovector_base ocbase2h=ocview;
      // ovector_base &ocbase1h=ocview;
      // const ovector_base &ocbase2h=ocview;
    }

    {

      // normal views of const ovectors
      // (Not allowed because ovector_view data is not const.)
      // ovector_view ocview1a=cvec;
      // const ovector_view ocview2a=cvec;
      // ovector_view &ocview1a=cvec;
      // const ovector_view &ocview2a=cvec;

      // normal views of non-const ovectors
      ovector_view ocview1b=vec;
      const ovector_view ocview2b=vec;

      // normal views of const ovector_views
      ovector_view ocview1c=cview;
      const ovector_view ocview2c=cview;

      // normal views of non-const ovector_views
      ovector_view ocview1d=view;
      const ovector_view ocview2d=view;

      // normal views of const ovector_base
      // (Not allowed to prevent non-const copying of const ovectors.)
      // ovector_view ocview1e=cbview;
      // const ovector_view ocview2e=cbview;
      // ovector_view &ocview1e=cbview;
      // const ovector_view &ocview2e=cbview;
      
      // normal views of non-const ovector_base
      ovector_view ocview1f=bview;
      const ovector_view ocview2f=bview;
      
      // normal views of const ovector_const_views
      // (Not allowed because ovector_view data is not const.)
      // ovector_view ocview1g=cocview;
      // const ovector_view ocview2g=cocview;
      // ovector_view &ocview1g=cocview;
      // const ovector_view &ocview2g=cocview;

      // normal views of non-const ovector_const_views
      // (Not allowed because ovector_view data is not const.)
      // ovector_view ocview1h=ocview;
      // const ovector_view ocview2h=ocview;
      // ovector_view &ocview1h=ocview;
      // const ovector_view &ocview2h=ocview;
    }

    {
      // Check that an const ovector_base reference still has const data.
      const ovector_base &ref=cvec;
      // const ovector_view viewx=ref;
      // const ovector_view &viewx=ref;
      // const ovector_base base=ref;
    }

  }

  {
    // Make several starting objects
    // uvector
    uvector vec(2);
    vec[0]=2.0;
    vec[1]=1.0;
    // const uvector
    const uvector cvec=vec;
    // vector_view
    uvector_view view=vec;
    // const vector_view
    const uvector_view cview=vec;
    // vector_base
    uvector_base bview=vec;
    // const vector_base
    const uvector_base cbview=vec;
    // vector_const_view
    uvector_const_view ocview=cvec;
    // const vector_const_view
    const uvector_const_view cocview=cvec;
    
    {
      // These are all deep copies

      //uvector copies of const uvectors
      uvector vec1a=cvec;
      const uvector vec2a=cvec;
      // Show we're const
      // vec2a[0]=3.0;
      
      //uvector copies of non-const uvectors
      uvector vec1b=vec;
      const uvector vec2b=vec;
      // Show we're const
      // vec2b[0]=3.0;

      // const_views of const uvector_views
      uvector vec1c=cview;
      const uvector vec2c=cview;
      // Show we're const
      // vec2c[0]=3.0;

      // const_views of non-const uvector_views
      uvector vec1d=view;
      const uvector vec2d=view;
      // Show we're const
      // vec2d[0]=3.0;

      // const_views of const uvector_base
      uvector vec1e=cbview;
      const uvector vec2e=cbview;
      // Show we're const
      // vec2e[0]=3.0;

      // const_views of non-const uvector_base
      uvector vec1f=bview;
      const uvector vec2f=bview;
      // Show we're const
      // vec2f[0]=3.0;

      // const_views of const uvector_const_views
      uvector vec1g=cocview;
      const uvector vec2g=cocview;
      // Show we're const
      // vec2g[0]=3.0;

      // const_views of non-const uvector_const_views
      uvector vec1h=ocview;
      const uvector vec2h=ocview;
      // Show we're const
      // vec2h[0]=3.0;
    }

    {
      // const_views of const uvectors
      uvector_const_view ocview1a=cvec;
      const uvector_const_view ocview2a=cvec;
      // Show we're const
      // ocview1a[0]=3.0;
      // ocview2a[0]=3.0;

      // const_views of non-const uvectors
      uvector_const_view ocview1b=vec;
      const uvector_const_view ocview2b=vec;
      // Show we're const
      // ocview1b[0]=3.0;
      // ocview2b[0]=3.0;

      // const_views of const uvector_views
      uvector_const_view ocview1c=cview;
      const uvector_const_view ocview2c=cview;
      // Show we're const
      // ocview1c[0]=3.0;
      // ocview2c[0]=3.0;

      // const_views of non-const uvector_views
      uvector_const_view ocview1d=view;
      const uvector_const_view ocview2d=view;
      // Show we're const
      // ocview1d[0]=3.0;
      // ocview2d[0]=3.0;

      // const_views of const uvector_base
      uvector_const_view ocview1e=cbview;
      const uvector_const_view ocview2e=cbview;
      // Show we're const
      // ocview1e[0]=3.0;
      // ocview2e[0]=3.0;

      // const_views of non-const uvector_base
      uvector_const_view ocview1f=bview;
      const uvector_const_view ocview2f=bview;
      // Show we're const
      // ocview1f[0]=3.0;
      // ocview2f[0]=3.0;

      // const_views of const uvector_const_views
      uvector_const_view ocview1g=cocview;
      const uvector_const_view ocview2g=cocview;
      // Show we're const
      // ocview1g[0]=3.0;
      // ocview2g[0]=3.0;

      // const_views of non-const uvector_const_views
      uvector_const_view ocview1h=ocview;
      const uvector_const_view ocview2h=ocview;
      // Show we're const
      // ocview1h[0]=3.0;
      // ocview2h[0]=3.0;
    }

    {

      // base views of const uvectors 
      // (Not allowed because uvector_base data is not const.)
      // uvector_base ocbase1a=cvec;
      // const uvector_base ocbase2a=cvec;

      // Reference is allowed, as long as it's const
      // uvector_base &ocbase1a=cvec;
      const uvector_base &ocbase2a=cvec;

      // base views of non-const uvectors
      uvector_base ocbase1b=vec;
      const uvector_base ocbase2b=vec;

      // base views of const uvector_views
      // (Not allowed to prevent non-const copying of const uvectors.)
      // uvector_base ocbase1c=cview;
      // const uvector_base ocbase2c=cview;

      // Reference is allowed, as long as it's const
      // uvector_base &ocbase1c=cview;
      const uvector_base &ocbase2c=cview;

      // base views of non-const uvector_views
      uvector_base ocbase1d=view;
      const uvector_base ocbase2d=view;
      
      // base views of const uvector_base
      // (Not allowed to prevent non-const copying of const uvectors.)
      // uvector_base ocbase1e=cbview;
      // const uvector_base ocbase2e=cbview;

      // Reference is allowed, as long as it's const
      // uvector_base &ocbase1e=cbview;
      const uvector_base &ocbase2e=cbview;

      // base views of non-const uvector_base
      uvector_base ocbase1f=bview;
      const uvector_base ocbase2f=bview;
      
      // base views of const uvector_const_views
      // (Not allowed because uvector_base data is not const.)
      // uvector_base ocbase1g=cocview;
      // const uvector_base ocbase2g=cocview;
      // uvector_base &ocbase1g=cocview;
      // const uvector_base &ocbase2g=cocview;

      // base views of non-const uvector_const_views
      // (Not allowed because uvector_base data is not const.)
      // uvector_base ocbase1h=ocview;
      // const uvector_base ocbase2h=ocview;
      // uvector_base &ocbase1h=ocview;
      // const uvector_base &ocbase2h=ocview;
    }

    {

      // normal views of const uvectors
      // (Not allowed because uvector_view data is not const.)
      // uvector_view ocview1a=cvec;
      // const uvector_view ocview2a=cvec;
      // uvector_view &ocview1a=cvec;
      // const uvector_view &ocview2a=cvec;

      // normal views of non-const uvectors
      uvector_view ocview1b=vec;
      const uvector_view ocview2b=vec;

      // normal views of const uvector_views
      uvector_view ocview1c=cview;
      const uvector_view ocview2c=cview;

      // normal views of non-const uvector_views
      uvector_view ocview1d=view;
      const uvector_view ocview2d=view;

      // normal views of const uvector_base
      // (Not allowed to prevent non-const copying of const uvectors.)
      // uvector_view ocview1e=cbview;
      // const uvector_view ocview2e=cbview;
      // uvector_view &ocview1e=cbview;
      // const uvector_view &ocview2e=cbview;
      
      // normal views of non-const uvector_base
      uvector_view ocview1f=bview;
      const uvector_view ocview2f=bview;
      
      // normal views of const uvector_const_views
      // (Not allowed because uvector_view data is not const.)
      // uvector_view ocview1g=cocview;
      // const uvector_view ocview2g=cocview;
      // uvector_view &ocview1g=cocview;
      // const uvector_view &ocview2g=cocview;

      // normal views of non-const uvector_const_views
      // (Not allowed because uvector_view data is not const.)
      // uvector_view ocview1h=ocview;
      // const uvector_view ocview2h=ocview;
      // uvector_view &ocview1h=ocview;
      // const uvector_view &ocview2h=ocview;
    }

    {
      // Check that an const uvector_base reference still has const data.
      const uvector_base &ref=cvec;
      // const uvector_view viewx=ref;
      // const uvector_view &viewx=ref;
      // const uvector_base base=ref;
    }

  }

  t.report();
  return 0;
}
// End of example

