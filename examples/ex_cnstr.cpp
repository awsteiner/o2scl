#include <iostream>
#include <hal/table.h>
#include <hal/minimize.h>
#include <hal/base_ioc.h>

using namespace std;
using namespace hal;

int main(void) {
  table t;
  t.line_of_names("x lb clb co cco");
  cout.setf(ios::scientific);
  for(double x=3.0;x<=7.0;x+=0.02) {
    double line[5]={x,lower_bound(x,5.0,1.0,10.0),
		    cont_lower_bound(x,5.0,1.0,10.0),
		    constraint(x,5.0,1.0,10.0),
		    cont_constraint(x,5.0,1.0,10.0)};
    t.line_of_data(5,line);
  }
  
  base_ioc bio;
  collection co;
  text_out_file *tof=new text_out_file("ex_cnstr.out");
  co.out_one(tof,"table","t",&t);
  delete tof;
  return 0;
}
  
