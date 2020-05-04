/*
  -------------------------------------------------------------------
  
  Copyright (C) 2020, Andrew W. Steiner
  
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
#include <o2scl/cursesw.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

int main(int argc, char *argv[]) {

  test_mgr t;
  t.set_output_level(1);

#ifdef O2SCL_READLINE
  
  cursesw cw;

  if (argc>=2 && ((string)argv[1])==((string)"1")) {

    cout << "C++ mode, character and enter." << endl;
    char ch1;
    cin >> ch1;
    cout << "Character was: " << ch1 << endl;
    
    cw.init();
    cw.cw_printw("Ten characters without enter.\n");
    for(size_t i=0;i<10;i++) {
      int ch2=cw.cw_getch();
      cw.cw_printw(cw.identify(ch2)+"\n");
    }
    cw.finalize();
    
    cout << "C++ mode, character and enter." << endl;
    char ch2;
    cin >> ch2;
    cout << "Character was: " << ch2 << endl;
  }

  int row, col;
  o2scl::get_screen_size(row,col);
  cout << row << " " << col << endl;
  
#endif
  
  t.report();
  return 0;
}
