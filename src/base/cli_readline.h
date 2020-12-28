/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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
/** \file cli_readline.h
    \brief File defining \ref o2scl::cli_readline
*/
#ifndef O2SCL_CLI_READLINE_H
#define O2SCL_CLI_READLINE_H

#include <readline/readline.h>
#include <readline/history.h>

#include <o2scl/cli.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief An extension to \ref o2scl::cli which uses readline
      
      This header-only class requires the GNU <tt>readline</tt>
      library for use, but is not referenced by \o2 code at the
      moment to make the library usable without <tt>readline</tt>.

  */
  class cli_readline : public cli {
  
  protected:

    /// Buffer for readline
    char *line_read;

    /// String containing filename
    std::string histfile;

    /// Maximum history file size
    size_t msize;

  public:
    
    cli_readline(std::string fname="", size_t max_size=100) {
      line_read=0;
      msize=max_size;
      
      histfile=fname;
      if (histfile.size()>0) {
	read_history(histfile.c_str());
      }
    }
    
    ~cli_readline() {
      if (histfile.size()>0) {
	stifle_history(((int)msize));
	write_history(histfile.c_str());
      }
    }

    /** \brief Set history file
     */
    void set_histfile(std::string fname) {
      histfile=fname;
      if (histfile.size()>0) {
	read_history(histfile.c_str());
      }
      return;
    }
    
    /** \brief Function to get a string from the user
	\nothing
    */
    virtual char *cli_gets(const char *c) {
    
      /* If the buffer has already been allocated, return the memory to
	 the free pool.
      */
      if (line_read) {
	free(line_read);
	line_read=(char *)0;
      }
    
      line_read=readline(c);
    
      /* If the line has any text in it, save it on the history.
       */
      if (line_read && *line_read && histfile.size()>0) {
	add_history(line_read);
      }
    
      return(line_read);
    }

    //#endif
  
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
