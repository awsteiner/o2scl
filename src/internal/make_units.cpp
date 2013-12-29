#include <iostream>
#include <o2scl/convert_units_gnu.h>

using namespace std;
using namespace o2scl;

int main(void) {

  convert_units_gnu cug;
  cug.make_units_dat("units.dat",false,false,false);
  cug.make_units_dat("units_hck.dat",true,true,true);

  return 0;
}
