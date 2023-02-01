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
#ifndef O2SCL_CURSESW_H
#define O2SCL_CURSESW_H

/** \file cursesw.h
    \brief Desc
*/

#include <iostream>
#include <string>
#include <vector>

namespace o2scl {

#if defined(O2SCL_NCURSES) || defined(DOXYGEN)
  
  /** \brief Simple interface for curses
  */
  class cursesw {

  protected:
    
    /** \brief Desc
     */
    typedef struct key_def_s {
      int key;
      std::string desc;
    } key_def;
    
    /** \brief Desc
     */
    std::vector<key_def> key_list;

    /** \brief Desc
     */
    int row, col;
    
  public:

    /** \brief Desc
     */
    cursesw();

    /** \brief Desc
     */
    void init();

    /** \brief Desc
     */
    int cw_getch();
    
    /** \brief Desc
     */
    void finalize();

    /** \brief Desc
     */
    std::string identify(int ch);

    /** \brief Desc
     */
    void cw_printw(std::string);
    
  };

  /** \brief Use curses to determine window size
   */
  void get_screen_size_curses(int &row, int &col);

#endif

  /** \brief Use tput to determine window size
   */
  void get_screen_size_tput(int &row, int &col);

  /** \brief Use ioctl to determine window size

      If this returns non-zero, then the attempt to determine the 
      size failed.
   */
  int get_screen_size_ioctl(int &row, int &col);
  
}

#endif
