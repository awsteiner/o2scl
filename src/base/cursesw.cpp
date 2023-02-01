/*
  -------------------------------------------------------------------
  
  Copyright (C) 2020-2023, Andrew W. Steiner
  
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
#include <cctype>
// for ioctl
#include <sys/ioctl.h>
//for STDOUT_FILENO
#include <unistd.h>

#include <o2scl/cursesw.h>
#include <o2scl/string_conv.h>

using namespace std;
using namespace o2scl;

#ifdef O2SCL_NCURSES

#include "curses.h"

cursesw::cursesw() {

  std::vector<key_def> kd={{KEY_BREAK,"break"},
			   {KEY_DOWN,"down arrow"},
			   {KEY_UP,"up arrow"},
			   {KEY_LEFT,"left arrow"},
			   {KEY_RIGHT,"right arrow"},
			   {KEY_HOME,"home"},
			   {KEY_BACKSPACE,"backspace"},
			   {KEY_F0,"F0"},
			   {KEY_F(1),"F1"},
			   {KEY_F(2),"F2"},
			   {KEY_F(3),"F3"},
			   {KEY_F(4),"F4"},
			   {KEY_F(5),"F5"},
			   {KEY_F(6),"F6"},
			   {KEY_F(7),"F7"},
			   {KEY_F(8),"F8"},
			   {KEY_F(9),"F9"},
			   {KEY_F(10),"F10"},
			   {KEY_F(11),"F11"},
			   {KEY_F(12),"F12"},
			   {KEY_DL,"delete line"},
			   {KEY_IL,"insert line"},
			   {KEY_DC,"delete character"},
			   {KEY_IC,"insert char or enter insert mode"},
			   {KEY_EIC,"exit insert char mode"},
			   {KEY_CLEAR,"clear screen"},
			   {KEY_EOS,"clear to end of screen"},
			   {KEY_EOL,"clear to end of line"},
			   {KEY_SF,"scroll 1 line forward"},
			   {KEY_SR,"scroll 1 line backward (reverse)"},
			   {KEY_NPAGE,"next page"},
			   {KEY_PPAGE,"previous page"},
			   {KEY_STAB,"set tab"},
			   {KEY_CTAB,"clear tab"},
			   {KEY_CATAB,"clear all tabs"},
			   {KEY_ENTER,"enter or send"},
			   {KEY_SRESET,"soft (partial) reset"},
			   {KEY_RESET,"reset or hard reset"},
			   {KEY_PRINT,"print or copy"},
			   {KEY_LL,"home down or bottom (lower left)"},
			   {KEY_A1,"upper left of keypad"},
			   {KEY_A3,"upper right of keypad"},
			   {KEY_B2,"center of keypad"},
			   {KEY_C1,"lower left of keypad"},
			   {KEY_C3,"lower right of keypad"},
			   {KEY_BTAB,"back tab"},
			   {KEY_BEG,"beginning"},
			   {KEY_CANCEL,"cancel"},
			   {KEY_CLOSE,"close"},
			   {KEY_COMMAND,"command"},
			   {KEY_COPY,"copy"},
			   {KEY_CREATE,"create"},
			   {KEY_END,"end"},
			   {KEY_EXIT,"exit"},
			   {KEY_FIND,"find"},
			   {KEY_HELP,"help"},
			   {KEY_MARK,"mark"},
			   {KEY_MESSAGE,"message"},
			   {KEY_MOUSE,"mouse event read"},
			   {KEY_MOVE,"move"},
			   {KEY_NEXT,"next object"},
			   {KEY_OPEN,"open"},
			   {KEY_OPTIONS,"options"},
			   {KEY_PREVIOUS,"previous object"},
			   {KEY_REDO,"redo"},
			   {KEY_REFERENCE,"reference"},
			   {KEY_REFRESH,"refresh"},
			   {KEY_REPLACE,"replace"},
			   {KEY_RESIZE,"screen resized"},
			   {KEY_RESTART,"restart"},
			   {KEY_RESUME,"resume"},
			   {KEY_SAVE,"save"},
			   {KEY_SBEG,"shifted beginning"},
			   {KEY_SCANCEL,"shifted cancel"},
			   {KEY_SCOMMAND,"shifted command"},
			   {KEY_SCOPY,"shifted copy"},
			   {KEY_SCREATE,"shifted create"},
			   {KEY_SDC,"shifted delete char"},
			   {KEY_SDL,"shifted delete line"},
			   {KEY_SELECT,"Select"},
			   {KEY_SEND,"shifted end"},
			   {KEY_SEOL,"shifted clear line"},
			   {KEY_SEXIT,"shifted exit"},
			   {KEY_SFIND,"shifted find"},
			   {KEY_SHELP,"shifted help"},
			   {KEY_SHOME,"shifted home"},
			   {KEY_SIC,"shifted input"},
			   {KEY_SLEFT,"shifted left arrow"},
			   {KEY_SMESSAGE,"shifted message"},
			   {KEY_SMOVE,"shifted move"},
			   {KEY_SNEXT,"shifted next"},
			   {KEY_SOPTIONS,"shifted options"},
			   {KEY_SPREVIOUS,"shifted prev"},
			   {KEY_SPRINT,"shifted print"},
			   {KEY_SREDO,"shifted redo"},
			   {KEY_SREPLACE,"shifted replace"},
			   {KEY_SRIGHT,"shifted right arrow"},
			   {KEY_SRSUME,"shifted resume"},
			   {KEY_SSAVE,"shifted save"},
			   {KEY_SSUSPEND,"shifted suspend"},
			   {KEY_SUNDO,"shifted undo"},
			   {KEY_SUSPEND,"suspend"},
			   {KEY_UNDO,"undo"}};
  key_list=kd;
};

void cursesw::init() {

  // Initialize curses session and return terminal to normal status
  initscr();
  
  // Ensure characters are not output to the terminal
  noecho();
  
  // Enable keypad. This appears to be critical to getting
  // the arrows to work correctly on OSX.
  keypad(stdscr,TRUE);

  int row, col;
  getmaxyx(stdscr,row,col);
  
}

int cursesw::cw_getch() {
  return getch();
}

void cursesw::finalize() {
  
  // Return to normal "cooked" mode with line buffering
  noraw();
  
  // Finalize curses session and return terminal to normal status
  endwin();
  
}

std::string cursesw::identify(int ch) {

  if (isprint(ch)) {
    char ch2=ch;
    std::string ret;
    ret+=ch2;
    return ret;
  }
  
  for(size_t i=0;i<key_list.size();i++) {
    if (ch==key_list[i].key) {
      return key_list[i].desc;
    }
  }
  
  return "not found";
}

void cursesw::cw_printw(std::string s) {
  printw(s.c_str());
  return;
}

void o2scl::get_screen_size_curses(int &row, int &col) {
  // Use the curses function
  initscr();
  getmaxyx(stdscr,row,col);
  endwin();
  return;
}

#endif

void o2scl::get_screen_size_tput(int &row, int &col) {
  string s_row, s_col;
  pipe_cmd_string("tput lines",s_row);
  pipe_cmd_string("tput cols",s_col);
  row=o2scl::stoi(s_row);
  col=o2scl::stoi(s_col);
  return;
}

int o2scl::get_screen_size_ioctl(int &row, int &col) {
  struct winsize w;
  int ret=ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
  
  row=w.ws_row;
  col=w.ws_col;
  
  return ret;
}
