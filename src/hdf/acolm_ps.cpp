/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
#include "acolm.h"

#include <regex>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/cloud_file.h>
#include <o2scl/vector_derint.h>
#include <o2scl/cursesw.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_acol;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

int acol_manager::comm_ser_hist_t3d(std::vector<std::string> &sv,
                                    bool itive_com) {

  vector<string> in, pr;
  pr.push_back("Vector spec. for grid");
  pr.push_back("Direction (x or y)");
  pr.push_back("Grid name");
  pr.push_back("Vector spec. for bin edges (or \"auto\")");
  pr.push_back("Vector spec. for bin grid (or \"auto\" or <size>)");
  pr.push_back("Bin name");
  pr.push_back("Pattern");
  pr.push_back("New slice name");
  int ret=get_input(sv,pr,in,"ser-hist-t3d",itive_com);
  if (ret!=0) return ret;
  
  std::vector<double> grid, bin_edges, bin_grid;
  int ret2=vector_spec(in[0],grid,0,false);

  if (in[3]=="auto") {
    size_t n_bins=o2scl::stoszt(in[4]);
    table3d_obj.create_table_hist_set(grid,in[1],in[2],n_bins,
                                      in[5],table_obj,in[6],in[7]);
  } else if (in[4]=="auto") {
    int ret3=vector_spec(in[3],bin_edges,0,false);
    table3d_obj.create_table_hist_set(grid,in[1],in[2],bin_edges,
                                      in[5],table_obj,in[6],in[7]);
  } else {
    int ret3=vector_spec(in[3],bin_edges,0,false);
    int ret4=vector_spec(in[4],bin_grid,0,false);
    table3d_obj.create_table_hist_set(grid,in[1],in[2],bin_edges,
                                      bin_grid,in[5],table_obj,in[6],in[7]);
  }

  command_del(type);
  clear_obj();
  command_add("table3d");
  type="table3d";

  return 0;
}

int acol_manager::comm_refine(std::vector<std::string> &sv, bool itive_com) {

  vector<string> in, pr;
  pr.push_back("Index column");
  pr.push_back("Refinement factor");
  int ret=get_input(sv,pr,in,"slice",itive_com);
  if (ret!=0) return ret;

  string index=in[0];
  size_t factor=o2scl::stoszt(in[1]);

  if (type=="table") {

    o2scl::table_units<> table_new;

    // Copy over column names and units
    for(size_t j=0;j<table_obj.get_ncolumns();j++) {
      cout << "New column: " << table_obj.get_column_name(j) << endl;
      table_new.new_column(table_obj.get_column_name(j));
      table_new.set_unit(table_obj.get_column_name(j),
			 table_obj.get_unit(table_obj.get_column_name(j)));
    }

    table_new.set_nlines((table_obj.get_nlines()-1)*factor+1);

    for(size_t k=0;k<table_obj.get_nlines()-1;k++) {
      for(size_t j=0;j<factor;j++) {
	table_new.set(index,k*factor+j,table_obj.get(index,k)+
		      (table_obj.get(index,k+1)-table_obj.get(index,k))/
		      ((double)factor)*((double)j));
      }
    }
    table_new.set(index,table_new.get_nlines()-1,
		  table_obj.get(index,table_obj.get_nlines()-1));

    table_new.set_interp_type(table_obj.get_interp_type());
    table_new.insert_table(table_obj,index);
    table_obj=table_new;
    
  } else {
    cerr << "Refine does not work with " << type << " objects." << endl;
  }
    
  return 0;
}

int acol_manager::comm_slack(std::vector<std::string> &sv, bool itive_com) {

  if (smess.url.length()==0) {
    if (smess.set_url_from_env("O2SCL_SLACK_URL")==false) {
      cerr << "Slack webhook URL not specified." << endl;
      return 1;
    }
    cout << "Set Slack URL to " << smess.url << endl;
  }

  string channel, message, image_url, alt_text;
  smess.verbose=verbose;

  if (sv.size()>=2 && sv[1]=="image") {
    if (sv.size()>=6) {
      channel=sv[2];
      message=sv[3];
      image_url=sv[4];
      alt_text=sv[5];
    } else if (sv.size()>=5) {
      message=sv[2];
      image_url=sv[3];
      alt_text=sv[4];
    } else {
      cerr << "Image specified, but not enough arguments given." << endl;
      return 10;
    }
  } else if (sv.size()>=3) {
    channel=sv[1];
    message=sv[2];
  } else if (sv.size()>=2) {
    message=sv[1];
  } else {
    cerr << "No slack message given." << endl;
    return 4;
  }

  // Take care of slack channel
  if (channel.length()>0) {
    smess.channel=channel;
  } else {
    if (smess.channel.length()==0) {
      if (smess.set_channel_from_env("O2SCL_SLACK_CHANNEL")==false) {
	cerr << "Slack channel not specified." << endl;
	return 2;
      }
    }
  }
  if (smess.channel[0]!='#') {
    smess.channel=((string)("#"))+smess.channel;
  }
  cout << "Set Slack channel to " << smess.channel << endl;

  // Take care of slack username
  if (smess.username.length()==0) {
    if (smess.set_username_from_env("O2SCL_SLACK_USERNAME")==false) {
      cerr << "Slack username not specified." << endl;
      return 3;
    }
    cout << "Set Slack username to " << smess.username << endl;
  }

  // Construct array of strings from string specification
  std::vector<std::string> slist;
  int ss_ret;
  if (sv.size()>=3) {
    ss_ret=strings_spec(message,slist,3,true);
  } else {
    ss_ret=strings_spec(message,slist,3,true);
  }
  if (ss_ret!=0) {
    cerr << "String specification failed." << endl;
  }

  // Collect array of strings into single string
  std::string stmp;
  for(size_t j=0;j<slist.size();j++) {
    stmp+=slist[j];
    if (j!=slist.size()-1) stmp+="\n";
  }

  // Now send the message
  int sret;
  if (image_url.length()>0) {
    if (verbose>=2) {
      cout << "message: " << stmp << endl;
      cout << "image_url: " << image_url << endl;
      cout << "alt_text: " << alt_text << endl;
    }
    sret=smess.send_image(stmp,image_url,alt_text,false);
  } else {
    if (verbose>=2) {
      cout << "message: " << stmp << endl;
    }
    sret=smess.send(stmp,false);
  }
  if (sret!=0) {
    cout << "Class slack_messenger failed to send message (error "
         << sret << ")." << endl;
    return sret;
  }
  
  return 0;
}

int acol_manager::comm_preview(std::vector<std::string> &sv, bool itive_com) {

  if (type.length()==0) {
    cerr << "No object to preview." << endl;
    return 3;
  }
  
  if (scientific) cout.setf(ios::scientific);
  else cout.unsetf(ios::scientific);
  
  cout.precision(precision);

  int ncols_loc;
  if (ncols<=0) {
    int srow, scol;
    int iret=get_screen_size_ioctl(srow,scol);
    if (scol>10 || iret!=0) ncols_loc=scol;
    else ncols_loc=80;
  } else {
    ncols_loc=ncols;
  }

  cout << "Type: " << type << " Name: " << obj_name << endl;
  
  if (type=="table3d") {
    
    int lmar=table3d_obj.get_y_name().length()+1;
    
    size_t nx, ny;
    table3d_obj.get_size(nx,ny);
    if (nx==0 || ny==0) {
      cout << "No size set. Blank table3d." << endl;
    } else {
      int nrows, ncls;
      if (sv.size()>=2) {
	int sret=o2scl::stoi_nothrow(sv[1],nrows);
	if (sret!=0 || nrows<=0) {
	  cerr << "Failed to interpret " << sv[1]
	       << "as a positive number of rows." << endl;
	  return exc_einval;
	}
      } else {
	nrows=10;
      }
      if (((size_t)nrows)>nx) nrows=nx;
      if (sv.size()>=3) {
	int sret=o2scl::stoi_nothrow(sv[2],ncls);
	if (sret!=0 || nrows<=0) {
	  cerr << "Failed to interpret " << sv[2]
	       << "as a positive number of rows." << endl;
	  return exc_einval;
	}
      } else {
	// 8+prec for the grid point, 4 for extra spacing,
	// and lmar for the left margin which has the x label
	if (ncols_loc<=precision+lmar+12) ncls=1;
	else ncls=(ncols_loc-precision-12-lmar)/(precision+8);
	if (verbose>1) {
	  std::cout << "Screen width: " << ncols_loc << " prec: " << precision
		    << " lmar: " << lmar << " flag: "
                    << (ncols_loc<=precision+lmar+12)
		    << " ncols_loc: " << ncls << endl;
	}
      }
      if (((size_t)ncls)>ny) ncls=ny;
      int dx=nx/ncls;
      int dy=ny/nrows;
      if (dx==0) dx=1;
      if (dy==0) dy=1;

      cout << "x: " << table3d_obj.get_x_name() << " [";
      if (nx<4) {
	for(size_t i=0;i<nx-1;i++) {
	  cout << table3d_obj.get_grid_x(i) << " ";
	}
	cout << table3d_obj.get_grid_x(nx-1) << "] ";
      } else {
	cout << table3d_obj.get_grid_x(0) << " ";
	cout << table3d_obj.get_grid_x(1) << " ... ";
	cout << table3d_obj.get_grid_x(nx-2) << " ";
	cout << table3d_obj.get_grid_x(nx-1) << "] ";
      }
      cout << endl;
      cout << "y: " << table3d_obj.get_y_name() << " [";
      if (ny<4) {
	for(size_t i=0;i<ny-1;i++) {
	  cout << table3d_obj.get_grid_y(i) << " ";
	}
	cout << table3d_obj.get_grid_y(ny-1) << "] ";
      } else {
	cout << table3d_obj.get_grid_y(0) << " ";
	cout << table3d_obj.get_grid_y(1) << " ... ";
	cout << table3d_obj.get_grid_y(ny-2) << " ";
	cout << table3d_obj.get_grid_y(ny-1) << "] ";
      }
      cout << endl;
      
      size_t nt=table3d_obj.get_nslices(); 
      if (nt!=0) {
	for(size_t k=0;k<nt;k++) {
	  // Slice name
	  cout << "Slice " << k << ": "
	       << table3d_obj.get_slice_name(k) << endl;

	  // Label row
	  for(int i=0;i<lmar+14;i++) cout << " ";
	  cout << "   " << table3d_obj.get_x_name() << endl;

	  // Set showpos
	  cout.setf(ios::showpos);

	  // Grid row
	  for(size_t i=0;i<((size_t)precision)+8+lmar;i++) cout << " ";
	  cout << "| ";
	  
	  for(size_t i=0;i<((size_t)ncls);i++) {
	    cout << table3d_obj.get_grid_x(i*dx) << " ";
	  }
	  cout << endl;

	  // Divider row
	  for(int i=0;i<lmar;i++) cout << " ";
	  for(size_t i=0;i<((size_t)precision)+8;i++) cout << "-";
	  cout << "|";
	  for(size_t i=0;i<((size_t)ncls)*(precision+8);i++) {
	    cout << "-";
	  }
	  cout << endl;

	  // Data output
	  for(int j=((int)nrows)-1;j>=0;j--) {
	    if (j==((int)nrows)-1) {
	      cout << table3d_obj.get_y_name() << " ";
	    } else {
	      for(int i=0;i<lmar;i++) cout << " ";
	    }
	    cout << table3d_obj.get_grid_y(j*dy) << " | ";
	    for(size_t i=0;i<((size_t)ncls);i++) {
	      cout << table3d_obj.get(i*dx,j*dy,k) << " ";
	    }
	    cout << endl;
	  }

	  // Unset showpos
	  cout.unsetf(ios::showpos);

	  // Newline between slices
	  cout << endl;
	}
      }
    }
    return 0;

  } else if (type=="hist") {

    int nrows=10;
    if (sv.size()>=2) {
      int sret=o2scl::stoi_nothrow(sv[1],nrows);
      if (sret!=0 || nrows<=0) {
	cerr << "Failed to interpret " << sv[1]
	     << "as a positive number of rows." << endl;
	return exc_einval;
      }
    }
    int inr=(hist_obj.size()+(nrows-1))/(nrows);
    if (inr<1) inr=1;

    cout.precision(precision);

    cout.setf(ios::left);
    cout.width(precision+8);
    cout << " low";
    cout.width(precision+8);
    cout << " high";
    cout.width(precision+8);
    cout << " weight" << endl;
    cout.unsetf(ios::left);

    for(int i=0;i<((int)hist_obj.size());i+=inr) {
      double val=hist_obj.get_bin_low_i(i);
      if (val<0.0) {
	cout << val << " ";
      } else {
	cout << " " << val << " ";
      }
      val=hist_obj.get_bin_high_i(i);
      if (val<0.0) {
	cout << val << " ";
      } else {
	cout << " " << val << " ";
      }
      val=hist_obj.get_wgt_i(i);
      if (val<0.0) {
	cout << val << endl;
      } else {
	cout << " " << val << endl;
      }
    }
    
    return 0;

  } else if (type=="hist_2d") {

    size_t nx=hist_2d_obj.size_x(), ny=hist_2d_obj.size_y();
    if (nx==0 || ny==0) {
      cout << "No size set. Blank hist_2d." << endl;
    } else {

      int nrows, ncls;
      if (sv.size()>=2) {
	int sret=o2scl::stoi_nothrow(sv[1],nrows);
	if (sret!=0 || nrows<=0) {
	  cerr << "Failed to interpret " << sv[1]
	       << "as a positive number of rows." << endl;
	  return exc_einval;
	}
      } else {
	nrows=10;
      }
      if (((size_t)nrows)>nx) nrows=nx;
      if (sv.size()>=3) {
	int sret=o2scl::stoi_nothrow(sv[2],ncls);
	if (sret!=0 || ncls<=0) {
	  cerr << "Failed to interpret " << sv[2]
	       << "as a positive number of columns.." << endl;
	  return exc_einval;
	}
      } else {
	// 8+prec for the grid point, 3 for extra spacing,
	if (ncols_loc<=precision+11) ncls=1;
	else ncls=(ncols_loc-precision-11)/(precision+8);
	if (verbose>1) {
	  std::cout << "Screen width: " << ncols_loc << " prec: " << precision
		    << " flag: " << (ncols_loc<=precision+11)
		    << " ncols_loc: " << ncls << endl;
	}
      }
      if (((size_t)ncls)>ny) ncls=ny;
      size_t dx=nx/ncls;
      size_t dy=ny/nrows;
      if (dx==0) dx=1;
      if (dy==0) dy=1;

      cout << "nx: " << nx << " ny: " << ny << endl;
      if (nx<=4) {
	cout << "x: [";
	for(size_t i=0;i<nx-1;i++) {
	  cout << hist_2d_obj.get_x_low_i(i) << " ";
	}
	cout << hist_2d_obj.get_x_low_i(nx-1) << "]" << endl;
	cout << "   [";
	for(size_t i=0;i<nx-1;i++) {
	  cout << hist_2d_obj.get_x_rep_i(i) << " ";
	}
	cout << hist_2d_obj.get_x_rep_i(nx-1) << "]" << endl;
	cout << "   [";
	for(size_t i=0;i<nx-1;i++) {
	  cout << hist_2d_obj.get_x_high_i(i) << " ";
	}
	cout << hist_2d_obj.get_x_high_i(nx-1) << "]" << endl;
      } else {
	cout << "x: [";
	cout << hist_2d_obj.get_x_low_i(0) << " ";
	cout << hist_2d_obj.get_x_low_i(1) << " ... ";
	cout << hist_2d_obj.get_x_low_i(nx-2) << " ";
	cout << hist_2d_obj.get_x_low_i(nx-1) << "]" << endl;
	cout << "   [";
	cout << hist_2d_obj.get_x_rep_i(0) << " ";
	cout << hist_2d_obj.get_x_rep_i(1) << " ... ";
	cout << hist_2d_obj.get_x_rep_i(nx-2) << " ";
	cout << hist_2d_obj.get_x_rep_i(nx-1) << "]" << endl;
	cout << "   [";
	cout << hist_2d_obj.get_x_high_i(0) << " ";
	cout << hist_2d_obj.get_x_high_i(1) << " ... ";
	cout << hist_2d_obj.get_x_high_i(nx-2) << " ";
	cout << hist_2d_obj.get_x_high_i(nx-1) << "]" << endl;
      }
      cout << endl;
      if (ny<=4) {
	cout << "y: [";
	for(size_t i=0;i<ny-1;i++) {
	  cout << hist_2d_obj.get_y_low_i(i) << " ";
	}
	cout << hist_2d_obj.get_y_low_i(ny-1) << "]" << endl;
	cout << "   [";
	for(size_t i=0;i<ny-1;i++) {
	  cout << hist_2d_obj.get_y_rep_i(i) << " ";
	}
	cout << hist_2d_obj.get_y_rep_i(ny-1) << "]" << endl;
	cout << "   [";
	for(size_t i=0;i<ny-1;i++) {
	  cout << hist_2d_obj.get_y_high_i(i) << " ";
	}
	cout << hist_2d_obj.get_y_high_i(ny-1) << "]" << endl;
      } else {
	cout << "y: [";
	cout << hist_2d_obj.get_y_low_i(0) << " ";
	cout << hist_2d_obj.get_y_low_i(1) << " ... ";
	cout << hist_2d_obj.get_y_low_i(ny-2) << " ";
	cout << hist_2d_obj.get_y_low_i(ny-1) << "]" << endl;
	cout << "   [";
	cout << hist_2d_obj.get_y_rep_i(0) << " ";
	cout << hist_2d_obj.get_y_rep_i(1) << " ... ";
	cout << hist_2d_obj.get_y_rep_i(ny-2) << " ";
	cout << hist_2d_obj.get_y_rep_i(ny-1) << "]" << endl;
	cout << "   [";
	cout << hist_2d_obj.get_y_high_i(0) << " ";
	cout << hist_2d_obj.get_y_high_i(1) << " ... ";
	cout << hist_2d_obj.get_y_high_i(ny-2) << " ";
	cout << hist_2d_obj.get_y_high_i(ny-1) << "]" << endl;
      }
      cout << endl;
      
      // Set showpos
      cout.setf(ios::showpos);
      
      // Grid row
      for(size_t i=0;i<((size_t)precision)+8;i++) cout << " ";
      cout << "| ";
      
      for(size_t i=0;i<((size_t)ncls);i++) {
	cout << hist_2d_obj.get_x_rep_i(i*dx) << " ";
      }
      cout << endl;
      
      // Divider row
      for(size_t i=0;i<((size_t)precision)+8;i++) cout << "-";
      cout << "|";
      for(size_t i=0;i<((size_t)ncls)*(precision+8);i++) {
	cout << "-";
      }
      cout << endl;
      
      // Data output
      for(int j=((int)nrows)-1;j>=0;j--) {
	cout << hist_2d_obj.get_y_rep_i(j*dy) << " | ";
	for(size_t i=0;i<((size_t)ncls);i++) {
	  cout << hist_2d_obj.get_wgt_i(i*dx,j*dy) << " ";
	}
	cout << endl;
      }
      
      // Unset showpos
      cout.unsetf(ios::showpos);
      
    }
    
    return 0;
    
  } else if (type=="table") {
    
    if (table_obj.get_nlines()==0) {
      cerr << "No table to preview." << endl;
      return exc_efailed;
    }
    
    if (table_obj.get_ncolumns()>0) {
      
      string nlast;
      int inr;
      int inc;
      
      //----------------------------------------------------------------------
      // Compute number of columns which will fit
      
      size_t max_cols=(ncols_loc)/(8+precision);
      if (max_cols>table_obj.get_ncolumns()) max_cols=table_obj.get_ncolumns();
      
      //--------------------------------------------------------------------
      // Compute column and row increment
      
      if (sv.size()==2) {
	int nrows;
	int ret=o2scl::stoi_nothrow(sv[1],nrows);
	if (ret!=0 || nrows<=0) {
	  std::cerr << "Could not interpret " << sv[1]
		    << " as a positive and nonzero number of rows." << endl;
	  return exc_efailed;
	}
	inr=(table_obj.get_nlines()+(nrows-1))/(nrows);
	if (inr<1) inr=1;
      } else {
	inr=(table_obj.get_nlines()+9)/10;
	if (inr<1) inr=1;
      }
      inc=(table_obj.get_ncolumns()+(1))/max_cols;
      if (inc<1) inc=1;
      
      //--------------------------------------------------------------------
      // Get last row number if necessary
      
      if (pretty==true) {
	nlast=itos(table_obj.get_nlines()-1);
      }
      
      //--------------------------------------------------------------------
      // Output column names
      
      if (names_out==true) {
	
	for(size_t ki=0;ki<max_cols;ki++) {
	  
	  size_t i=ki*inc;
	  if (i>=table_obj.get_ncolumns()) i=table_obj.get_ncolumns()-1;
	  
	  // Preceeding space
	  if (pretty==true) {
	    cout << ' ';
	  }
	  
	  // Column name
	  cout << table_obj.get_column_name(i) << " ";
	  
	  // Trailing spaces
	  if (pretty==true) {
	    int nsp=precision+6-((int)(table_obj.get_column_name(i).size()));
	    for(int j=0;j<nsp;j++) cout << ' ';
	  } else {
	    for(size_t kk=1;kk<nlast.length();kk++) cout << ' ';
	  }
	  
	}
	cout << endl;
      }
      
      //--------------------------------------------------------------------
      // Output units
      
      if (names_out==true && table_obj.get_nunits()>0) {
	
	for(size_t ki=0;ki<max_cols;ki++) {
	  
	  size_t i=ki*inc;
	  if (i>=table_obj.get_ncolumns()) i=table_obj.get_ncolumns()-1;
	  
	  // Preceeding space
	  if (pretty==true) {
	    cout << ' ';
	  }
	  
	  // Column name
	  string cunit=table_obj.get_unit(table_obj.get_column_name(i));
	  cout << '[' << cunit << "] ";
	  
	  // Trailing spaces
	  if (pretty==true) {
	    int nsp=precision+6-cunit.size()-2;
	    if (nsp<0) nsp=0;
	    for(int j=0;j<nsp;j++) cout << ' ';
	  } else {
	    for(size_t kk=1;kk<nlast.length();kk++) cout << ' ';
	  }
	  
	}
	cout << endl;
      }
      
      //--------------------------------------------------------------------
      // Output data
      
      for(size_t i=0;i<table_obj.get_nlines();i+=inr) {
	
	for(size_t kj=0;kj<max_cols;kj++) {
	  
	  size_t j=kj*inc;
	  if (j>=table_obj.get_ncolumns()) j=table_obj.get_ncolumns()-1;
	  
	  if (pretty==true) {
	    double d=table_obj.get(j,i);
	    if (!has_minus_sign(&d)) {
	      cout << ' ';
	    }
	  }
	  cout << table_obj.get(j,i) << ' ';
	  if (pretty==true) {
	    for(int kk=0;kk<((int)(table_obj.get_column_name(j).size()-
				   precision-6));kk++) {
	      cout << ' ';
	    }
	  }
	  
	}
	cout << endl;
	
	//--------------------------------------------------------------------
	// Continue to next row
      }
      
    }
    
    return 0;

  } else if (type=="vector<contour_line>") {

    for(size_t i=0;i<cont_obj.size();i++) {
      cout << "Value: " << cont_obj[i].level << endl;
      for(size_t j=0;j<cont_obj[i].x.size();j++) {
	cout << cont_obj[i].x[j] << " ";
	cout << cont_obj[i].y[j] << endl;
      }
      cout << endl;
    }
    
    return 0;
  } else if (type=="int") {
    cout << "Value is " << int_obj << endl;
    return 0;
  } else if (type=="double") {
    cout << "Value is " << double_obj << endl;
    return 0;
  } else if (type=="char") {
    cout << "Value is " << char_obj << endl;
    return 0;
  } else if (type=="string") {
    cout << "Value is " << string_obj << endl;
    return 0;
  } else if (type=="size_t") {
    cout << "Value is " << size_t_obj << endl;
    return 0;
  } else if (type=="int[]") {
    
    vector<string> inc, outc;
    for(size_t i=0;i<intv_obj.size();i++) {
      string tmp=o2scl::szttos(i)+". "+o2scl::itos(intv_obj[i]);
      inc.push_back(tmp);
    }
    screenify(inc.size(),inc,outc);
    
    int inr;
    if (sv.size()==2) {
      int nrows;
      int ret=o2scl::stoi_nothrow(sv[1],nrows);
      if (ret!=0 || nrows<=0) {
	std::cerr << "Could not interpret " << sv[1]
		  << " as a positive and nonzero number of rows." << endl;
	return exc_efailed;
      }
      inr=(outc.size()+(nrows-1))/(nrows);
      if (inr<1) inr=1;
    } else {
      inr=(outc.size()+9)/10;
      if (inr<1) inr=1;
    }

    for(size_t i=0;i<outc.size();i+=inr) {
      cout << outc[i] << endl;
    }
    return 0;

  } else if (type=="double[]") {

    vector<string> inc, outc;
    for(size_t i=0;i<doublev_obj.size();i++) {
      string tmp;
      if (has_minus_sign(&doublev_obj[i])) {
	tmp=o2scl::szttos(i)+". "+o2scl::dtos(doublev_obj[i]);
      } else {
	tmp=o2scl::szttos(i)+".  "+o2scl::dtos(doublev_obj[i]);
      }
      inc.push_back(tmp);
    }
    screenify(inc.size(),inc,outc);

    int inr;
    if (sv.size()==2) {
      int nrows;
      int ret=o2scl::stoi_nothrow(sv[1],nrows);
      if (ret!=0 || nrows<=0) {
	std::cerr << "Could not interpret " << sv[1]
		  << " as a positive and nonzero number of rows." << endl;
	return exc_efailed;
      }
      inr=(outc.size()+(nrows-1))/(nrows);
      if (inr<1) inr=1;
    } else {
      inr=(outc.size()+9)/10;
      if (inr<1) inr=1;
    }

    for(size_t i=0;i<outc.size();i+=inr) {
      cout << outc[i] << endl;
    }
    return 0;

  } else if (type=="size_t[]") {

    vector<string> inc, outc;
    for(size_t i=0;i<size_tv_obj.size();i++) {
      string tmp=o2scl::szttos(i)+". "+o2scl::szttos(size_tv_obj[i]);
      inc.push_back(tmp);
    }
    screenify(inc.size(),inc,outc);
    
    int inr;
    if (sv.size()==2) {
      int nrows;
      int ret=o2scl::stoi_nothrow(sv[1],nrows);
      if (ret!=0 || nrows<=0) {
	std::cerr << "Could not interpret " << sv[1]
		  << " as a positive and nonzero number of rows." << endl;
	return exc_efailed;
      }
      inr=(outc.size()+(nrows-1))/(nrows);
      if (inr<1) inr=1;
    } else {
      inr=(outc.size()+9)/10;
      if (inr<1) inr=1;
    }

    for(size_t i=0;i<outc.size();i+=inr) {
      cout << outc[i] << endl;
    }
    return 0;

  } else if (type=="string[]") {

    vector<string> inc, outc;
    for(size_t i=0;i<stringv_obj.size();i++) {
      string tmp=o2scl::szttos(i)+". \""+stringv_obj[i]+"\"";
      inc.push_back(tmp);
    }
    screenify(inc.size(),inc,outc);
    
    int inr;
    if (sv.size()==2) {
      int nrows;
      int ret=o2scl::stoi_nothrow(sv[1],nrows);
      if (ret!=0 || nrows<=0) {
	std::cerr << "Could not interpret " << sv[1]
		  << " as a positive and nonzero number of rows." << endl;
	return exc_efailed;
      }
      inr=(outc.size()+(nrows-1))/(nrows);
      if (inr<1) inr=1;
    } else {
      inr=(outc.size()+9)/10;
      if (inr<1) inr=1;
    }

    for(size_t i=0;i<outc.size();i+=inr) {
      cout << outc[i] << endl;
    }
    return 0;

  } else if (type=="uniform_grid<double>") {
    
    cout << "Uniform grid " << obj_name << endl;
    cout << "Number of bins: " << ug_obj.get_nbins() << endl;
    cout << "Start: " << ug_obj.get_start() << endl;
    cout << "End: " << ug_obj.get_end() << endl;
    cout << "Width: " << ug_obj.get_width() << endl;
    
    return 0;
    
  } else if (type=="tensor_grid") {
    
    size_t rk=tensor_grid_obj.get_rank();
    cout << "Rank: " << rk << endl;
    for(size_t i=0;i<rk;i++) {
      size_t sz=tensor_grid_obj.get_size(i);
      cout << i << " : ";
      vector<double> grid(sz);
      tensor_grid_obj.copy_grid(i,grid);
      if (grid.size()==4) {
	cout << grid[0] << " " << grid[1] << " " << grid[2] << " "
	     << grid[3] << endl;
      } if (grid.size()==3) {
	cout << grid[0] << " " << grid[1] << " " << grid[2] << endl;
      } else if (grid.size()==2) {
	cout << grid[0] << " " << grid[1] << endl;
      } else if (grid.size()==1) {
	cout << grid[0] << endl;
      } else {
	cout << grid[0] << " " << grid[1] << " ... "
	     << grid[sz-2] << " " << grid[sz-1] << endl;
      }
    }

    int nrows=10;
    if (sv.size()>=2) {
      int sret=o2scl::stoi_nothrow(sv[1],nrows);
      if (sret!=0 || nrows<=0) {
	cerr << "Failed to interpret " << sv[1]
	     << "as a positive number of rows." << endl;
	return exc_einval;
      }
    }
    
    size_t total_size=tensor_grid_obj.total_size();

    double x=tensor_grid_obj.get_data()[total_size-1];
    vector<size_t> ix(rk);
    tensor_grid_obj.unpack_index(total_size-1,ix);
    string test="(";
    for(size_t i=0;i<rk-1;i++) {
      test+=o2scl::dtos(tensor_grid_obj.get_grid(i,ix[i]))+",";
    }
    test+=o2scl::dtos(tensor_grid_obj.get_grid(rk-1,ix[rk-1]))+"): ";
    test+=o2scl::dtos(x);
    size_t maxwid=test.length();
    if (!has_minus_sign(&x)) maxwid++;

    // Number of 'columns' is equal to the number of columns
    // in the screen 'ncols_loc' divided by the maximum width
    // of one column
    size_t nct=ncols_loc/maxwid;
    size_t step=total_size/nrows/nct;
    vector<string> svin, svout;
    for(size_t i=0;i<total_size;i+=step) {
      tensor_grid_obj.unpack_index(i,ix);
      string stemp="(";
      for(size_t j=0;j<rk-1;j++) {
	stemp+=o2scl::dtos(tensor_grid_obj.get_grid(j,ix[j]))+",";
      }
      stemp+=o2scl::dtos(tensor_grid_obj.get_grid(rk-1,ix[rk-1]))+"): ";
      stemp+=o2scl::dtos(tensor_grid_obj.get_data()[i]);
      svin.push_back(stemp);
    }
    screenify(svin.size(),svin,svout);
    for(size_t i=0;i<svout.size();i++) {
      cout << svout[i] << endl;
    }
    
    return 0;
    
  } else if (type=="tensor") {
    
    int nrows=10;
    if (sv.size()>=2) {
      int sret=o2scl::stoi_nothrow(sv[1],nrows);
      if (sret!=0 || nrows<=0) {
	cerr << "Failed to interpret " << sv[1]
	     << "as a positive number of rows." << endl;
	return exc_einval;
      }
    }
    
    size_t rk=tensor_obj.get_rank();
    cout << "Rank: " << rk << " sizes: (";
    for(size_t i=0;i<rk-1;i++) {
      cout << tensor_obj.get_size(i) << ",";
    }
    cout << tensor_obj.get_size(rk-1) << ")" << endl;
    
    size_t total_size=tensor_obj.total_size();

    double x=tensor_obj.get_data()[total_size-1];
    vector<size_t> ix(rk);
    tensor_obj.unpack_index(total_size-1,ix);
    string test="(";
    for(size_t i=0;i<rk-1;i++) {
      test+=o2scl::szttos(ix[i])+",";
    }
    test+=o2scl::szttos(ix[rk-1])+"): ";
    test+=o2scl::dtos(x);
    size_t maxwid=test.length();
    if (!has_minus_sign(&x)) maxwid++;

    // Number of 'columns' is equal to the number of columns
    // in the screen 'ncols_loc' divided by the maximum width
    // of one column
    size_t nct=ncols_loc/maxwid;
    size_t step=total_size/nrows/nct;
    vector<string> svin, svout;
    for(size_t i=0;i<total_size;i+=step) {
      tensor_obj.unpack_index(i,ix);
      string stemp="(";
      for(size_t j=0;j<rk-1;j++) {
	stemp+=o2scl::szttos(ix[j])+",";
      }
      stemp+=o2scl::szttos(ix[rk-1])+"): ";
      stemp+=o2scl::dtos(tensor_obj.get_data()[i]);
      svin.push_back(stemp);
    }
    screenify(svin.size(),svin,svout);
    for(size_t i=0;i<svout.size();i++) {
      cout << svout[i] << endl;
    }
    
    return 0;

  } else if (type=="tensor<int>") {
    
    int nrows=10;
    if (sv.size()>=2) {
      int sret=o2scl::stoi_nothrow(sv[1],nrows);
      if (sret!=0 || nrows<=0) {
	cerr << "Failed to interpret " << sv[1]
	     << "as a positive number of rows." << endl;
	return exc_einval;
      }
    }
    
    size_t rk=tensor_obj.get_rank();
    cout << "Rank: " << rk << " sizes: (";
    for(size_t i=0;i<rk-1;i++) {
      cout << tensor_obj.get_size(i) << ",";
    }
    cout << tensor_obj.get_size(rk-1) << ")" << endl;
    
    size_t total_size=tensor_obj.total_size();

    int x=tensor_obj.get_data()[total_size-1];
    vector<size_t> ix(rk);
    tensor_obj.unpack_index(total_size-1,ix);
    string test="(";
    for(size_t i=0;i<rk-1;i++) {
      test+=o2scl::szttos(ix[i])+",";
    }
    test+=o2scl::szttos(ix[rk-1])+"): ";
    test+=o2scl::itos(x);
    size_t maxwid=test.length();
    if (x>=0) maxwid++;

    // Number of 'columns' is equal to the number of columns
    // in the screen 'ncols_loc' divided by the maximum width
    // of one column
    size_t nct=ncols_loc/maxwid;
    size_t step=total_size/nrows/nct;
    vector<string> svin, svout;
    for(size_t i=0;i<total_size;i+=step) {
      tensor_obj.unpack_index(i,ix);
      string stemp="(";
      for(size_t j=0;j<rk-1;j++) {
	stemp+=o2scl::szttos(ix[j])+",";
      }
      stemp+=o2scl::szttos(ix[rk-1])+"): ";
      stemp+=o2scl::itos(tensor_obj.get_data()[i]);
      svin.push_back(stemp);
    }
    screenify(svin.size(),svin,svout);
    for(size_t i=0;i<svout.size();i++) {
      cout << svout[i] << endl;
    }
    
    return 0;

  } else if (type=="tensor<size_t>") {
    
    int nrows=10;
    if (sv.size()>=2) {
      int sret=o2scl::stoi_nothrow(sv[1],nrows);
      if (sret!=0 || nrows<=0) {
	cerr << "Failed to interpret " << sv[1]
	     << "as a positive number of rows." << endl;
	return exc_einval;
      }
    }
    
    size_t rk=tensor_obj.get_rank();
    cout << "Rank: " << rk << " sizes: (";
    for(size_t i=0;i<rk-1;i++) {
      cout << tensor_obj.get_size(i) << ",";
    }
    cout << tensor_obj.get_size(rk-1) << ")" << endl;
    
    size_t total_size=tensor_obj.total_size();

    size_t x=tensor_obj.get_data()[total_size-1];
    vector<size_t> ix(rk);
    tensor_obj.unpack_index(total_size-1,ix);
    string test="(";
    for(size_t i=0;i<rk-1;i++) {
      test+=o2scl::szttos(ix[i])+",";
    }
    test+=o2scl::szttos(ix[rk-1])+"): ";
    test+=o2scl::szttos(x);
    size_t maxwid=test.length();
    if (x>=0) maxwid++;

    // Number of 'columns' is equal to the number of columns
    // in the screen 'ncols_loc' divided by the maximum width
    // of one column
    size_t nct=ncols_loc/maxwid;
    size_t step=total_size/nrows/nct;
    vector<string> svin, svout;
    for(size_t i=0;i<total_size;i+=step) {
      tensor_obj.unpack_index(i,ix);
      string stemp="(";
      for(size_t j=0;j<rk-1;j++) {
	stemp+=o2scl::szttos(ix[j])+",";
      }
      stemp+=o2scl::szttos(ix[rk-1])+"): ";
      stemp+=o2scl::szttos(tensor_obj.get_data()[i]);
      svin.push_back(stemp);
    }
    screenify(svin.size(),svin,svout);
    for(size_t i=0;i<svout.size();i++) {
      cout << svout[i] << endl;
    }
    
    return 0;
  }

  cerr << "Cannot preview type " << type << "." << endl;
  return 0;
}

int acol_manager::comm_slice(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table3d") {
  
    if (table3d_obj.get_nslices()==0) {
      cerr << "No data in current table3d object." << endl;
      return exc_efailed;
    }
    
    vector<string> in, pr;
    pr.push_back("Slice in 'x' or 'y' direction");
    pr.push_back("Value to interpolate in for slice");
    int ret=get_input(sv,pr,in,"slice",itive_com);
    if (ret!=0) return ret;
    
    if (sv[1]=="x") {
      table3d_obj.extract_x(std::stod(sv[2]),table_obj);
      command_del(type);
      clear_obj();
      command_add("table");
      type="table";
    } else if (sv[1]=="y") {
      table3d_obj.extract_y(std::stod(sv[2]),table_obj);
      command_del(type);
      clear_obj();
      command_add("table");
      type="table";
    } else {
      cerr << "Invalid first argument to 'slice' must be either 'x' or 'y'."
	   << endl;
    }

  } else if (type=="tensor_grid") {

    size_t rank=tensor_grid_obj.get_rank();
    size_t nfix=(sv.size()-1)/2;

    vector<size_t> ifix;
    vector<double> vals;
    
    if (nfix==0) {
      std::string i1;
      int ret=get_input_one(sv,"Number of indices to fix",i1,"slice",
			    itive_com);
      if (ret!=0) return ret;

      nfix=o2scl::stoszt(i1);
      if (nfix==0) {
	cerr << "User specified zero indices to fix." << endl;
	return 1;
      }

      vector<string> sv2;
      vector<string> in(nfix*2), pr(nfix*2);
      for(size_t j=0;j<nfix;j++) {
	pr.push_back("Index to fix");
	pr.push_back("Value to fix it at");
      }
      ret=get_input(sv2,pr,in,"slice",itive_com);
      if (ret!=0) return ret;

      for(size_t j=0;j<nfix;j++) {
	ifix.push_back(o2scl::stoszt(in[2*j+0]));
	vals.push_back(o2scl::stod(in[2*j+1]));
      }
      
    } else {

      for(size_t j=0;j<nfix;j++) {
	ifix.push_back(o2scl::stoszt(sv[2*j+1]));
	vals.push_back(o2scl::stod(sv[2*j+2]));
      }

    }

    tensor_grid<> tg_old=tensor_grid_obj;
    tensor_grid_obj=tg_old.copy_slice_interp(ifix,vals);

    cout << "Old rank is " << rank << " and new rank is "
	 << tensor_grid_obj.get_rank() << endl;
    
  } else {
    cerr << "Slice does not work with " << type << " objects." << endl;
  }
    
  
  return 0;
}

int acol_manager::comm_slice_hist(std::vector<std::string> &sv,
				  bool itive_com) {

  if (type=="table3d") {

    if (sv.size()<2) {
      cerr << "Command 'slice-hist' needs a slice argument." << endl;
      return 1;
    }

    double min=matrix_min_value<ubmatrix,double>
      (table3d_obj.get_slice(sv[1]));
    double max=matrix_max_value<ubmatrix,double>
      (table3d_obj.get_slice(sv[1]));
    uniform_grid<double> grid=uniform_grid_end<double>(min,max,40);
    hist_obj.set_bin_edges(grid);
    for(size_t i=0;i<table3d_obj.get_nx();i++) {
      for(size_t j=0;j<table3d_obj.get_ny();j++) {
	hist_obj.update(table3d_obj.get(i,j,sv[1]));
      }
    }

    command_del(type);
    clear_obj();
    command_add("hist");
    type="hist";
    
  } else {
    cerr << "Command 'slice-hist' does not work with "
	 << type << " objects." << endl;
  }
    
  
  return 0;
}

int acol_manager::comm_set_data(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table3d") {

    if (type!="table3d") {
      cerr << "No table3d with data to set." << endl;
      return exc_efailed;
    }
    
    vector<string> in, pr;
    pr.push_back(table3d_obj.get_x_name()+" value of point to set");
    pr.push_back(table3d_obj.get_y_name()+" value of point to set");
    pr.push_back("Slice name to set");
    pr.push_back("New value");
    int ret=get_input(sv,pr,in,"set-data",itive_com);
    if (ret!=0) return ret;
    
    table3d_obj.set_val(o2scl::stod(in[0]),o2scl::stod(in[1]),in[2],
			o2scl::stod(in[3]));
    return 0;
  }

  if (table_obj.get_nlines()==0) {
    cerr << "No table to insert columns into." << endl;
    return exc_efailed;
  }

  vector<string> in, pr;
  pr.push_back("Function describing rows to modify");
  pr.push_back("Name of column to modify");
  pr.push_back("Function describing new values");
  int ret=get_input(sv,pr,in,"set-data",itive_com);
  if (ret!=0) return ret;
  
  if (table_obj.is_column(in[1])==false) {
    cerr << "Could not find column named '" << in[1] << "'." << endl;
    return exc_efailed;
  }

  for(size_t i=0;i<table_obj.get_nlines();i++) {
    if (table_obj.row_function(in[0],i)>0.5) {
      table_obj.set(in[1],i,table_obj.row_function(in[2],i));
    }
  }

  return 0;
}

int acol_manager::comm_set_unit(std::vector<std::string> &sv, bool itive_com) {

  if (type!="table") {
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }

  if (table_obj.get_nlines()==0) {
    cerr << "No table to set units of." << endl;
    return exc_efailed;
  }

  vector<string> in, pr;
  pr.push_back("Column to set units of");
  pr.push_back("New unit");
  int ret=get_input(sv,pr,in,"set-unit",itive_com);
  if (ret!=0) return ret;
  
  if (table_obj.is_column(in[0])==false) {
    cerr << "Could not find column named '" << in[0] << "'." << endl;
    return exc_efailed;
  }

  table_obj.set_unit(in[0],in[1]);

  return 0;
}

int acol_manager::comm_rename(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table3d") {

    vector<string> pr, in;
    pr.push_back("Enter slice to be renamed");
    pr.push_back("Enter new name");
    
    int ret=get_input(sv,pr,in,"rename",itive_com);
    if (ret!=0) return ret;
    
    table3d_obj.rename_slice(in[0],in[1]);
    
    return 0;

  } else if (type=="table") {

    vector<string> pr, in;
    pr.push_back("Enter column to be renamed");
    pr.push_back("Enter new name");
    
    int ret=get_input(sv,pr,in,"rename",itive_com);
    if (ret!=0) return ret;
    
    if (table_obj.is_column(in[0])==false) {
      cerr << "Could not find column named '" << in[0] << "'." << endl;
      return exc_efailed;
    }
    
    table_obj.new_column(in[1]);
    
    table_obj.copy_column(in[0],in[1]);
    table_obj.delete_column(in[0]);

  } else {
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }    

  return 0;
}

int acol_manager::comm_sort(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table") {
  
    std::string i1, i2;
    
    if (table_obj.get_nlines()==0) {
      cerr << "No table to sort." << endl;
      return exc_efailed;
    }
    
    if (sv.size()>2) {
      i1=sv[1];
      i2=sv[2];
    } else if (sv.size()>1) {
      i1=sv[1];
    } else {
      if (itive_com) {
	i1=cl->cli_gets("Enter column to sort by (or blank to stop): ");
	if (i1.length()==0) {
	  cout << "Command 'sort' cancelled." << endl;
	  return 0;
	}
      } else {
	cerr << "Not enough arguments for 'sort'." << endl;
	return exc_efailed;
      }
    }
    
    bool unique=false;
    if (i2==((std::string)"unique")) unique=true;
    
    if (table_obj.is_column(i1)==false) {
      cerr << "Could not find column named '" << i1 << "'." << endl;
      return exc_efailed;
    }
    
    if (verbose>1) {
      cout << "Sorting by column " << i1 << endl; 
    }
    table_obj.sort_table(i1);
    
    if (unique) {
      if (verbose>0) {
	std::cout << "Making sort unique by deleting identical rows. "
		  << "Starting with " << table_obj.get_nlines() << " rows."
		  << std::endl;
      }
      table_obj.delete_idadj_rows();
      if (verbose>0) {
	std::cout << "Done. Now " << table_obj.get_nlines() << " rows."
		  << std::endl;
      }
    }

  } else if (type=="double[]") {

    vector_sort<vector<double>,double>(doublev_obj.size(),doublev_obj);
    if (verbose>0) {
      cout << "Object of type double[] sorted." << endl;
    }
    
  } else if (type=="int[]") {

    vector_sort<vector<int>,int>(intv_obj.size(),intv_obj);
    if (verbose>0) {
      cout << "Object of type double[] sorted." << endl;
    }
    
  } else if (type=="size_t[]") {

    vector_sort<vector<size_t>,size_t>(size_tv_obj.size(),size_tv_obj);
    if (verbose>0) {
      cout << "Object of type double[] sorted." << endl;
    }

  } else {

    cout << "Not implemented for type " << type << endl;
    return 1;
    
  }
  
  return 0;
}

int acol_manager::comm_stats(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table") {
    
    if (table_obj.get_nlines()==0) {
      cerr << "No table to analyze." << endl;
      return exc_efailed;
    }
    
    std::string i1;
    int ret=get_input_one(sv,"Enter column to get info on",i1,"stats",
			  itive_com);
    if (ret!=0) return ret;
    
    if (table_obj.is_column(i1)==false) {
      cerr << "Could not find column named '" << i1 << "'." << endl;
      return exc_efailed;
    }
    
    const vector<double> &cref=table_obj.get_column(i1);
    cout << "N        : " << table_obj.get_nlines() << endl;
    cout << "Sum      : " << vector_mean(table_obj.get_nlines(),cref)*
      table_obj.get_nlines() << endl;
    cout << "Mean     : " << vector_mean(table_obj.get_nlines(),cref) << endl;
    cout << "Std. dev.: " << vector_stddev(table_obj.get_nlines(),cref) << endl;
    size_t ix;
    double val;
    vector_min(table_obj.get_nlines(),cref,ix,val);
    cout << "Min      : " << val << " at index: " << ix << endl;
    vector_max(table_obj.get_nlines(),cref,ix,val);
    cout << "Max      : " << val << " at index: " << ix << endl;
    
    size_t dup=0, inc=0, dec=0;
    for(size_t i=0;i<table_obj.get_nlines()-1;i++) {
      if (cref[i+1]==cref[i]) dup++;
      if (cref[i]<cref[i+1]) inc++;
      if (cref[i]>cref[i+1]) dec++;
    }

    size_t ninf=0, nnan=0;
    for(size_t i=0;i<table_obj.get_nlines();i++) {
      if (std::isinf(cref[i])) ninf++;
      if (std::isnan(cref[i])) nnan++;
    }
    
    if (inc>0 && dec==0) {
      if (dup>0) {
	cout << "Increasing (" << dup << " duplicates)." << endl;
      } else {
	cout << "Strictly increasing. No duplicates." << endl;
      }
    } else if (dec>0 && inc==0) {
      if (dup>0) {
	cout << "Decreasing (" << dup << " duplicates)." << endl;
      } else {
	cout << "Strictly decreasing. No duplicates." << endl;
      }
    } else if (dec==0 && inc==0) {
      cout << "Constant (" << dup << " duplicates)." << endl;
    } else {
      cout << "Non-monotonic (" << inc << " increasing, " << dec 
	   << " decreasing, and " << dup << " duplicates)." << endl;
    }
    if (ninf>0) {
      cout << ninf << " infinite values." << endl;
    }
    if (nnan>0) {
      cout << nnan << " NaN values." << endl;
    }
    if ((dup+inc+dec)!=(table_obj.get_nlines()-1)) {
      cout << "Counting mismatch from non-finite values or signed zeros."
	   << endl;
    }
    
  } else if (type=="table3d") {
    
    if (table3d_obj.get_nx()==0) {
      cerr << "No table3d object to analyze." << endl;
      return exc_efailed;
    }
    
    std::string i1;
    int ret=get_input_one(sv,"Enter slice to get info on",i1,"stats",
			  itive_com);
    if (ret!=0) return ret;

    size_t sl_index;
    if (table3d_obj.is_slice(i1,sl_index)==false) {
      cerr << "Could not find slice named '" << i1 << "'." << endl;
      return exc_efailed;
    }
    
    const ubmatrix &cmref=table3d_obj.get_slice(i1);
    cout << "Nx       : " << table3d_obj.get_nx() << endl;
    cout << "Ny       : " << table3d_obj.get_ny() << endl;
    cout << "Sum      : "
	 << matrix_sum<ubmatrix,double>(table3d_obj.get_nx(),
					table3d_obj.get_ny(),
					cmref) << endl;
    double min, max;
    size_t i_min, i_max, j_min, j_max;
    matrix_minmax_index(table3d_obj.get_nx(),
			table3d_obj.get_ny(),cmref,i_min,j_min,min,
			i_max,j_max,max);
    cout << "Min      : " << min << " at (" << i_min << ","
	 << j_min << ")" << endl;
    cout << "Max      : " << max << " at (" << i_max << ","
	 << j_max << ")" << endl;
    
    size_t ninf=0, nnan=0;
    for(size_t i=0;i<table3d_obj.get_nx();i++) {
      for(size_t j=0;j<table3d_obj.get_ny();j++) {
	if (std::isinf(cmref(i,j))) ninf++;
	if (std::isnan(cmref(i,j))) nnan++;
      }
    }
    
    if (ninf>0) {
      cout << ninf << " infinite values." << endl;
    }
    if (nnan>0) {
      cout << nnan << " NaN values." << endl;
    }
    
  } else if (type=="tensor") {
    
    const std::vector<double> &data=tensor_obj.get_data();
    size_t N=data.size();
    cout << "N        : " << N << endl;
    cout << "Sum      : " << vector_mean(N,data)*N << endl;
    cout << "Mean     : " << vector_mean(N,data) << endl;
    cout << "Std. dev.: " << vector_stddev(N,data) << endl;
    size_t ix;
    vector<size_t> ix2(tensor_obj.get_rank());
    double val;
    vector_min(N,data,ix,val);
    tensor_obj.unpack_index(ix,ix2);
    cout << "Min      : " << val << " at indices: ";
    vector_out(cout,ix2,true);
    vector_max(N,data,ix,val);
    tensor_obj.unpack_index(ix,ix2);
    cout << "Max      : " << val << " at indices: ";
    vector_out(cout,ix2,true);
    
    size_t ninf=0, nnan=0;
    for(size_t i=0;i<N;i++) {
      if (std::isinf(data[i])) ninf++;
      if (std::isnan(data[i])) nnan++;
    }
    if (ninf>0) {
      cout << ninf << " infinite values." << endl;
    }
    if (nnan>0) {
      cout << nnan << " NaN values." << endl;
    }
  
  } else if (type=="tensor_grid") {
    
    const std::vector<double> &data=tensor_grid_obj.get_data();
    size_t N=data.size();
    cout << "N        : " << N << endl;
    cout << "Sum      : " << vector_mean(N,data)*N << endl;
    cout << "Mean     : " << vector_mean(N,data) << endl;
    cout << "Std. dev.: " << vector_stddev(N,data) << endl;
    size_t ix;
    vector<size_t> ix2(tensor_grid_obj.get_rank());
    double val;
    vector_min(N,data,ix,val);
    tensor_grid_obj.unpack_index(ix,ix2);
    cout << "Min      : " << val << " at indices: ";
    vector_out(cout,ix2,true);
    vector_max(N,data,ix,val);
    tensor_grid_obj.unpack_index(ix,ix2);
    cout << "Max      : " << val << " at indices: ";
    vector_out(cout,ix2,true);
    
    size_t ninf=0, nnan=0;
    for(size_t i=0;i<N;i++) {
      if (std::isinf(data[i])) ninf++;
      if (std::isnan(data[i])) nnan++;
    }
    if (ninf>0) {
      cout << ninf << " infinite values." << endl;
      if (ninf<10) {
	for(size_t i=0;i<N;i++) {
	  if (std::isinf(data[i])) {
	    tensor_grid_obj.unpack_index(i,ix2);
	    cout << "Value at (";
	    for(size_t j=0;j<tensor_grid_obj.get_rank();j++) {
	      if (j!=tensor_grid_obj.get_rank()-1) {
		cout << ix2[j] << ",";
	      } else {
		cout << ix2[j];
	      }
	    }
	    cout << ") (";
	    for(size_t j=0;j<tensor_grid_obj.get_rank();j++) {
	      if (j!=tensor_grid_obj.get_rank()-1) {
		cout << tensor_grid_obj.get_grid(j,ix2[j]) << ",";
	      } else {
		cout << tensor_grid_obj.get_grid(j,ix2[j]);
	      }
	    }
	    cout << ") is " << data[i] << endl;
	  }
	}
      }
    }
    if (nnan>0) {
      cout << nnan << " NaN values." << endl;
      if (nnan<10) {
	for(size_t i=0;i<N;i++) {
	  if (std::isnan(data[i])) {
	    tensor_grid_obj.unpack_index(i,ix2);
	    cout << "Value at (";
	    for(size_t j=0;j<tensor_grid_obj.get_rank();j++) {
	      if (j!=tensor_grid_obj.get_rank()-1) {
		cout << ix2[j] << ",";
	      } else {
		cout << ix2[j];
	      }
	    }
	    cout << ") (";
	    for(size_t j=0;j<tensor_grid_obj.get_rank();j++) {
	      if (j!=tensor_grid_obj.get_rank()-1) {
		cout << tensor_grid_obj.get_grid(j,ix2[j]) << ",";
	      } else {
		cout << tensor_grid_obj.get_grid(j,ix2[j]);
	      }
	    }
	    cout << ") is " << data[i] << endl;
	  }
	}
      }
    }
  
  } else {
    
    cout << "Not implemented for type " << type << endl;
    return 1;
    
  }
  
  return 0;
}

int acol_manager::comm_set(std::vector<std::string> &sv, bool itive_com) {

  // Make sure the object interpolation types coincide with the
  // variable setting
  if (type=="table") {
    table_obj.set_interp_type(interp_type);
  } else if (type=="table3d") {
    table3d_obj.set_interp_type(interp_type);
  } else if (type=="hist") {
    hist_obj.set_interp_type(interp_type);
  }

  return 0;
}

int acol_manager::comm_select(std::vector<std::string> &sv, bool itive_com) {

  if (type!="table" && type!="table3d") {
    cout << "Not implemented for type " << type << endl;
    return 0;
  }
  
  // Remove outermost single or double quotes from all arguments
  for(size_t i=0;i<sv.size();i++) {
    size_t svs=sv[i].size();
    if (svs>2 && ((sv[i][0]=='\'' && sv[i][svs-1]=='\'') ||
		  (sv[i][0]=='\"' && sv[i][svs-1]=='\"'))) {
      sv[i]=sv[i].substr(1,sv[i].size()-2);
    }
  }

  // ----------------------------------------------------------------
  // Parse arguments into names and values
  // ----------------------------------------------------------------

  // The number of column arguments
  int nargs=((int)sv.size())-1;
  // The name part of the argument
  std::vector<std::string> names;
  // The function part of the argument
  std::vector<std::string> args;
  // Each element is true if the corresponding argument is a pattern
  std::vector<bool> is_pattern;
  
  // ----------------------------------------------------------------
  
  if (nargs==0) {
    cerr << "No arguments to 'select'. Command 'select' aborted." << endl;
    return exc_efailed;
  }
  
  // Reserve enough space
  names.reserve(nargs);
  args.reserve(nargs);
  is_pattern.reserve(nargs);
  
  // Go through each column
  for(int i=0;i<nargs;i++) {

    args.push_back(sv[i+1]);
    
    if (args[i][0]==':') {
      
      // If it has a colon, add to the patterns list
      args[i]=args[i].substr(1,args[i].length()-1);
      is_pattern.push_back(true);
      names.push_back("");
      
    } else {
      
      // If there's no colon, look for an '=' sign
      int tmpindx;
      tmpindx=args[i].find('=',0);
      
      if (tmpindx!=((int)string::npos)) {
	
	// Parse column name and function
	names.push_back(args[i].substr(0,tmpindx));
	args[i]=args[i].substr(tmpindx+1,args[i].size()-tmpindx-1);
	
      } else {
	
	// There's no equals sign, just assume its a column name
	names.push_back(args[i]);
      }
      
      is_pattern.push_back(false);
      
    }
  }
  
  // Output the results of the parsing
  if (verbose>2) {
    cout << "Parsed arguments: " << endl;
    for(int i=0;i<nargs;i++) {
      cout.width(2);
      cout << i << " ";
      if (is_pattern[i]==true) {
	cout << "Pattern: " << args[i] << endl;
      } else {
	cout << names[i] << " = " << args[i] << endl;
      }
    }
  }

  if (type=="table3d") {

    if (type!="table3d") {
      cerr << "No table3d to select from." << endl;
      return exc_efailed;
    }
    
    // ---------------------------------------------------------------------
    // Create new table3d and copy grid over
    // ---------------------------------------------------------------------

    table3d *new_table3d=new table3d;
    size_t nx, ny;
    table3d_obj.get_size(nx,ny);
    new_table3d->set_xy(table3d_obj.get_x_name(),nx,table3d_obj.get_x_data(),
			table3d_obj.get_y_name(),ny,table3d_obj.get_y_data());
	
    // ---------------------------------------------------------------------
    // Copy constants from old to new table3d
    // ---------------------------------------------------------------------

    for(size_t i=0;i<table3d_obj.get_nconsts();i++) {
      string tnam;
      double tval;
      table3d_obj.get_constant(i,tnam,tval);
      if (verbose>2) {
	cout << "Adding constant " << tnam << " = " << tval << endl;
      }
      new_table3d->add_constant(tnam,tval);
    }
  
    // ---------------------------------------------------------------------
    // Add slides and data to new table3d
    // ---------------------------------------------------------------------

    std::vector<bool> matched(table3d_obj.get_nslices());
    for(size_t i=0;i<table3d_obj.get_nslices();i++) {
      matched[i]=false;
    }

    // In this next loop, we need fix the code to ensure that when a
    // non-pattern column is given, it sets matched to true
    // appropriately.

    // Copy the data over
    for(int i=0;i<nargs;i++) {

      if (is_pattern[i]==false) {
	
        // Return an error if the slice doesn't exist
	size_t ix;
        if (names[i]==args[i] && table3d_obj.is_slice(args[i],ix)==false) {
          cerr << "Slice '" << args[i] << "' is not in the table." << endl;
          return exc_einval;
        }

        // Add the new slice to the new table
        new_table3d->new_slice(names[i]);

        // Fill slice with the new data
        ubmatrix mat(nx,ny);
        table3d_obj.function_matrix(args[i],mat,false);
        new_table3d->copy_to_slice(mat,names[i]);

      } else {

        if (use_regex==true) {
          
          regex r(args[i]);
        
          // Find the matching slices
          for(size_t j=0;j<table3d_obj.get_nslices();j++) {
	  
            if (matched[j]==false &&
                regex_search(table3d_obj.get_slice_name(j),r)) {
	    
              // If we've found a match, add it to the new table
              matched[j]=true;
	    
              // Add the new slice to the new table
              new_table3d->new_slice(table3d_obj.get_slice_name(j));
	    
              // Fill it with the new data
              ubmatrix mat(nx,ny);
              table3d_obj.function_matrix(args[i],mat,false);
              new_table3d->copy_to_slice(mat,table3d_obj.get_slice_name(j));
            }
          }
        
        } else {
          
          // Find the matching slices
          for(size_t j=0;j<table3d_obj.get_nslices();j++) {
	  
            if (matched[j]==false &&
                fnmatch(args[i].c_str(),
                        table3d_obj.get_slice_name(j).c_str(),0)==0) {
	    
              // If we've found a match, add it to the new table
              matched[j]=true;
	    
              // Add the new slice to the new table
              new_table3d->new_slice(table3d_obj.get_slice_name(j));
	    
              // Fill it with the new data
              ubmatrix mat(nx,ny);
              table3d_obj.function_matrix(args[i],mat,false);
              new_table3d->copy_to_slice(mat,table3d_obj.get_slice_name(j));
            }
          }
        
        }
      }
    }
    
    // Todo: Replace this copy with std::swap
    table3d_obj=*new_table3d;
    
  } else {

    if (table_obj.get_nlines()==0) {
      cerr << "No table to select from." << endl;
      return exc_efailed;
    }

    // ---------------------------------------------------------------------
    // Create new table
    // ---------------------------------------------------------------------
  
    table_units<> new_table;
  
    new_table.set_nlines(table_obj.get_nlines());

    // ---------------------------------------------------------------------
    // Copy constants from old to new table
    // ---------------------------------------------------------------------

    for(size_t i=0;i<table_obj.get_nconsts();i++) {
      string tnam;
      double tval;
      table_obj.get_constant(i,tnam,tval);
      new_table.add_constant(tnam,tval);
    }
  
    // ---------------------------------------------------------------------
    // Add columns and data to new table
    // ---------------------------------------------------------------------

    std::vector<bool> matched(table_obj.get_ncolumns());
    for(size_t i=0;i<table_obj.get_ncolumns();i++) {
      matched[i]=false;
    }

    // In this next loop, we need fix the code to ensure that when a
    // non-pattern column is given, it sets matched to true
    // appropriately.

    for(int i=0;i<nargs;i++) {

      if (is_pattern[i]==false) {
      
	// Return an error if the column doesn't exist
	if (names[i]==args[i] && table_obj.is_column(args[i])==false) {
	  cerr << "Column '" << args[i] << "' is not in the table." << endl;
	  return exc_einval;
	}

	// Add the new column to the new table
	new_table.new_column(names[i]);

	// If necessary, set units
	if (names[i]==args[i] && table_obj.get_unit(args[i]).length()!=0) {
	  new_table.set_unit(names[i],table_obj.get_unit(args[i]));
	}

	// Fill column with the new data
	ubvector vec(table_obj.get_nlines());
	table_obj.function_vector(args[i],vec,false);
	new_table.copy_to_column(vec,names[i]);

      } else {

        if (use_regex==true) {
          
          regex r(args[i]);
        
          // Find the matching columns
          for(size_t j=0;j<table_obj.get_ncolumns();j++) {

            if (matched[j]==false &&
                regex_search(table_obj.get_column_name(j),r)) {

              // If we've found a match, add it to the new table
              matched[j]=true;

              // Add the new column to the new table
              new_table.new_column(table_obj.get_column_name(j));

              // If necessary, set units
              string tmp=table_obj.get_column_name(j);
              if (table_obj.get_unit(tmp).length()!=0) {
                new_table.set_unit(tmp,table_obj.get_unit(tmp));
              }

              // Fill it with the new data
              ubvector vec(table_obj.get_nlines());
              table_obj.function_vector(table_obj.get_column_name(j),vec,false);
              new_table.copy_to_column(vec,table_obj.get_column_name(j));
            }
          }

        } else {
          
          // Find the matching columns
          for(size_t j=0;j<table_obj.get_ncolumns();j++) {

            if (matched[j]==false &&
                fnmatch(args[i].c_str(),
                        table_obj.get_column_name(j).c_str(),0)==0) {

              // If we've found a match, add it to the new table
              matched[j]=true;

              // Add the new column to the new table
              new_table.new_column(table_obj.get_column_name(j));

              // If necessary, set units
              string tmp=table_obj.get_column_name(j);
              if (table_obj.get_unit(tmp).length()!=0) {
                new_table.set_unit(tmp,table_obj.get_unit(tmp));
              }

              // Fill it with the new data
              ubvector vec(table_obj.get_nlines());
              table_obj.function_vector(table_obj.get_column_name(j),vec,false);
              new_table.copy_to_column(vec,table_obj.get_column_name(j));
            }
          }

        }
      }
    }

    // Todo: Replace this copy with std::swap
    std::swap(table_obj,new_table);

  }

  // Call list command
  if (verbose>0) {
    comm_list(sv,itive_com);
  }

  return 0;
}

int acol_manager::comm_select_rows(std::vector<std::string> &sv, 
                                   bool itive_com) {

  if (type!="table") {
    cout << "Not implemented for type " << type << endl;
    return 0;
  }

  if (table_obj.get_nlines()==0) {
    cerr << "No table to select rows from." << endl;
    return exc_efailed;
  }

  std::string i1;
  if (sv.size()>=2) {
    i1=sv[1];
  } else if (itive_com==true) {
    i1=cl->cli_gets("Function to specify rows (or blank to quit): ");
    if (i1.length()==0) {
      if (verbose>0) cout << "Command 'select_rows' cancelled." << endl;
      return 0;
    }
  } else {
    cerr << "No rows to delete." << endl;
    return exc_efailed;
  }
    
  // ---------------------------------------------------------------------
  // Create new table
  // ---------------------------------------------------------------------
  
  table_units<> new_table;
  
  // ---------------------------------------------------------------------
  // Copy constants from old to new table
  // ---------------------------------------------------------------------

  for(size_t i=0;i<table_obj.get_nconsts();i++) {
    string tnam;
    double tval;
    table_obj.get_constant(i,tnam,tval);
    new_table.add_constant(tnam,tval);
    if (verbose>=2) {
      cout << "Copying constant " << tnam << " " << tval
           << " to new table." << endl;
    }
  }
  
  // ---------------------------------------------------------------------
  // Add column names to new table
  // ---------------------------------------------------------------------

  for(int i=0;i<((int)table_obj.get_ncolumns());i++) {
    new_table.new_column(table_obj.get_column_name(i));
  }

  // ---------------------------------------------------------------------
  // Copy data from selected rows
  // ---------------------------------------------------------------------

  calc_utf8<> calc;
  calc.compile(i1.c_str(),0);
  vector<std::u32string> cols32=calc.get_var_list();
  vector<std::string> cols(cols32.size());
  for(size_t ij=0;ij<cols32.size();ij++) {
    char32_to_utf8(cols32[ij],cols[ij]);
  }
  if (verbose>=2) {
    cout << "Calculating expression: " << i1 << endl;
  }
  
  std::map<std::string,double> vars;

  int new_lines=0;
  for(int i=0;i<((int)table_obj.get_nlines());i++) {
    
    if (verbose>0 && i%10000==0) {
      std::cout << "Finished " << i << " of "
		<< table_obj.get_nlines() << " lines." << endl;
    }
    
    for(size_t j=0;j<cols.size();j++) {
      vars[cols[j]]=table_obj.get(cols[j],i);
      if (verbose>=2 && i==0) {
        cout << "At row 0, setting variable " << cols[j] << " to "
             << table_obj.get(cols[j],i) << endl;
      }
    }
    
    //calc.compile(i1.c_str(),&vars);
    if (calc.eval(&vars)>0.5) {
      
      // It is important to use set_nlines_auto() here because it
      // increases the table size fast enough to avoid poor scaling
      new_table.set_nlines_auto(new_lines+1);

      // Trivially parallize the assignment over all columns
#ifdef O2SCL_OPENMP
#pragma omp parallel
#endif
      {
#ifdef O2SCL_OPENMP
#pragma omp for
#endif
	for(int j=0;j<((int)table_obj.get_ncolumns());j++) {
	  new_table.set(j,new_lines,table_obj.get(j,i));
	}
	
	// End of parallel region
      }
      
      new_lines++;
    }
  }
  
  // Swap the old table with the new one
  std::swap(table_obj,new_table);
  
  return 0;
}

int acol_manager::comm_set_grid(std::vector<std::string> &sv, bool itive_com) {

  if (type=="tensor_grid") {

    size_t rank=tensor_grid_obj.get_rank();
    
    vector<string> pr, in;
    pr.push_back("Index of grid to set");
    pr.push_back("Function defining grid");
    int ret=get_input(sv,pr,in,"set-grid",itive_com);
    if (ret!=0) return ret;

    vector<double> grid;

    size_t k=o2scl::stoszt(in[0]);
    
    if (in[1].find(':')==std::string::npos) {
      
      calc_utf8<> calc;
      std::map<std::string,double> vars;
      for(size_t i=0;i<tensor_grid_obj.get_size(k);i++) {
	vars["i"]=((double)i);
	vars["x"]=tensor_grid_obj.get_grid(k,i);
	calc.compile(in[1].c_str(),&vars);
	double val=calc.eval(&vars);
	tensor_grid_obj.set_grid(k,i,val);
      }
      
    } else {
      
      std::vector<double> vtemp;
      int ret=vector_spec(in[1],vtemp,3,false);
      if (ret!=0) {
	cerr << "Interpretation of vector specification failed."
	     << endl;
	return 3;
      }
      
      if (vtemp.size()<tensor_grid_obj.get_size(k)) {
	cerr << "Vector specification results in vector "
	     << "smaller than tensor grid." << endl;
	return 2;
      }
      
      for(size_t i=0;i<tensor_grid_obj.get_size(k);i++) {
	tensor_grid_obj.set_grid(k,i,vtemp[i]);
      }
      
    }

    return 0;
  }

  cout << "Not implemented for type " << type << endl;
  return 1;
}

int acol_manager::comm_read(std::vector<std::string> &sv, 
			    bool itive_com) {

  vector<string> in, pr;
  if (sv.size()<2) {
    //|| (sv.size()<3 && itive_com)) {
    pr.push_back("Enter filename");
    pr.push_back("Enter object name");
    int ret=get_input(sv,pr,in,"table",itive_com);
    if (ret!=0) return ret;
  } else {
    in.resize(2);
    in[0]=sv[1];
    if (sv.size()>=3) {
      in[1]=sv[2];
    } 
  }
  
  // Delete previous object
  command_del(type);
  clear_obj();

  // Use hdf_file to open the file
  hdf_file hf;
  string type2;
  int ret;

  string fname_old=in[0];
  std::vector<std::string> matches;
  int wret=wordexp_wrapper(fname_old,matches);
  if (matches.size()>1 || matches.size()==0 || wret!=0) {
    cerr << "Function wordexp_wrapper() returned non-zero value. "
	 << "Bad filename?" << endl;
    return 1;
  }
  string fname=matches[0];
  if (verbose>1) {
    cout << "Function wordexp() converted "
	 << fname_old << " to " << fname << endl;
  }

  ret=hf.open(fname.c_str(),false,false);
  if (ret!=0) {
    cerr << "Could not find file named '" << fname << "'. Wrong file name?" 
	 << endl;
    return exc_efailed;
  }

  if (in[1].length()!=0) {

    if (verbose>1) {
      cout << "Command read looking for object with name " << in[1] << endl;
    }

    hdf_file::iterate_parms ip={in[1],&hf,false,type,verbose,
      hdf_file::ip_type_from_name,false};
    
    H5Literate(hf.get_current_id(),H5_INDEX_NAME,H5_ITER_NATIVE,
	       0,hdf_file::iterate_func,&ip);
    
    if (ip.found==false) {
      cerr << "Could not find object named " << in[1]
	   << " in file " << fname << endl;
      return 1;
    }
    
    if (verbose>1) {
      cout << "Command read found object with type " << ip.type << endl;
    }

    if (ip.type=="table") {
      if (verbose>2) {
	cout << "Reading table." << endl;
      }
      hdf_input(hf,table_obj,in[1]);
      obj_name=in[1];
      interp_type=table_obj.get_interp_type();
      command_add("table");
      type="table";
      return 0;
    } else if (ip.type=="table3d") {
      if (verbose>2) {
	cout << "Reading table3d." << endl;
      }
      hdf_input(hf,table3d_obj,in[1]);
      obj_name=in[1];
      interp_type=table3d_obj.get_interp_type();
      command_add("table3d");
      type="table3d";
      return 0;
    } else if (ip.type=="tensor_grid") {
      if (verbose>2) {
	cout << "Reading tensor_grid." << endl;
      }
      hdf_input(hf,tensor_grid_obj,in[1]);
      obj_name=in[1];
      command_add("tensor_grid");
      type="tensor_grid";
      return 0;
    } else if (ip.type=="prob_dens_mdim_amr") {
      if (verbose>2) {
	cout << "Reading prob_dens_mdim_amr." << endl;
      }
      hdf_input(hf,pdma_obj,in[1]);
      obj_name=in[1];
      command_add("prob_dens_mdim_amr");
      type="prob_dens_mdim_amr";
      return 0;
    } else if (ip.type.substr(0,10)==((string)"double[][]").substr(0,10)) {
      if (verbose>2) {
	cout << "Reading tensor." << endl;
      }
      hf.getd_ten(in[1],tensor_obj);
      obj_name=in[1];
      command_add("tensor");
      type="tensor";
      return 0;
    } else if (ip.type.substr(0,10)==((string)"int[][]").substr(0,10)) {
      if (verbose>2) {
	cout << "Reading tensor<int>." << endl;
      }
      hf.geti_ten(in[1],tensor_int_obj);
      obj_name=in[1];
      command_add("tensor<int>");
      type="tensor<int>";
      return 0;
    } else if (ip.type.substr(0,10)==((string)"size_t[][]").substr(0,10)) {
      if (verbose>2) {
	cout << "Reading tensor<size_t>." << endl;
      }
      hf.get_szt_ten(in[1],tensor_size_t_obj);
      obj_name=in[1];
      command_add("tensor<size_t>");
      type="tensor<size_t>";
      return 0;
    } else if (ip.type=="hist") {
      if (verbose>2) {
	cout << "Reading hist." << endl;
      }
      hdf_input(hf,hist_obj,in[1]);
      obj_name=in[1];
      command_add("hist");
      type="hist";
      return 0;
    } else if (ip.type=="hist_2d") {
      if (verbose>2) {
	cout << "Reading hist_2d." << endl;
      }
      hdf_input(hf,hist_2d_obj,in[1]);
      obj_name=in[1];
      command_add("hist_2d");
      type="hist_2d";
      return 0;
    } else if (ip.type=="vector<contour_line>") {
      if (verbose>2) {
	cout << "Reading vector<contour_line>." << endl;
      }
      hdf_input(hf,cont_obj,in[1]);
      obj_name=in[1];
      command_add("vector<contour_line>");
      type="vector<contour_line>";
      return 0;
    } else if (ip.type=="uniform_grid<double>") {
      if (verbose>2) {
	cout << "Reading uniform_grid<double>." << endl;
      }
      hdf_input(hf,ug_obj,in[1]);
      obj_name=in[1];
      command_add("uniform_grid<double>");
      type="uniform_grid<double>";
      return 0;
    } else if (ip.type=="string[]") {
      if (verbose>2) {
	cout << "Reading string[]." << endl;
      }
      hf.gets_vec(in[1],stringv_obj);
      obj_name=in[1];
      command_add("string[]");
      type="string[]";
      return 0;
    } else if (ip.type=="vec_vec_string") {
      if (verbose>2) {
	cout << "Reading vec_vec_string." << endl;
      }
      hf.gets_vec_vec(in[1],vvstring_obj);
      obj_name=in[1];
      command_add("vec_vec_string");
      type="vec_vec_string";
      return 0;
    } else if (ip.type=="int") {
      if (verbose>2) {
	cout << "Reading int." << endl;
      }
      hf.geti(in[1],int_obj);
      obj_name=in[1];
      command_add("int");
      type="int";
      return 0;
    } else if (ip.type=="char") {
      if (verbose>2) {
	cout << "Reading char." << endl;
      }
      hf.getc(in[1],char_obj);
      obj_name=in[1];
      command_add("char");
      type="char";
      return 0;
    } else if (ip.type=="string") {
      if (verbose>2) {
	cout << "Reading string." << endl;
      }
      hf.gets(in[1],string_obj);
      obj_name=in[1];
      command_add("string");
      type="string";
      return 0;
    } else if (ip.type=="char[fixed]") {
      if (true || verbose>2) {
	cout << "Reading char[fixed] and storing as string." << endl;
      }
      hf.gets_fixed(in[1],string_obj);
      obj_name=in[1];
      command_add("string");
      type="string";
      return 0;
    } else if (ip.type=="double") {
      if (verbose>2) {
	cout << "Reading double." << endl;
      }
      hf.getd(in[1],double_obj);
      obj_name=in[1];
      command_add("double");
      type="double";
      return 0;
    } else if (ip.type=="size_t") {
      if (verbose>2) {
	cout << "Reading size_t." << endl;
      }
      hf.get_szt(in[1],size_t_obj);
      obj_name=in[1];
      command_add("size_t");
      type="size_t";
      return 0;
    } else if (ip.type=="int[]") {
      if (verbose>2) {
	cout << "Reading int[]." << endl;
      }
      hf.geti_vec(in[1],intv_obj);
      obj_name=in[1];
      command_add("int[]");
      type="int[]";
      return 0;
    } else if (ip.type=="double[]") {
      if (verbose>2) {
	cout << "Reading double[]." << endl;
      }
      hf.getd_vec(in[1],doublev_obj);
      obj_name=in[1];
      command_add("double[]");
      type="double[]";
      return 0;
    } else if (ip.type=="size_t[]") {
      if (verbose>2) {
	cout << "Reading size_t[]." << endl;
      }
      hf.get_szt_vec(in[1],size_tv_obj);
      obj_name=in[1];
      command_add("size_t[]");
      type="size_t[]";
      return 0;
    }

    cerr << "Found object with name " << in[1]
	 << " in file " << fname << " but type " << ip.type
	 << " is not readable." << endl;
    return 2;
  }

  ret=hf.find_object_by_type("table",in[1],use_regex,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, found first table object named '"
	   << in[1] << "'." << endl;
    }
    hdf_input(hf,table_obj,in[1]);
    obj_name=in[1];
    command_add("table");
    type="table";
    interp_type=table_obj.get_interp_type();
    return 0;
  }
    
  ret=hf.find_object_by_type("table3d",in[1],use_regex,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, found first table3d object named '"
	   << in[1] << "'." << endl;
    }
    hdf_input(hf,table3d_obj,in[1]);
    obj_name=in[1];
    interp_type=table3d_obj.get_interp_type();
      
    command_add("table3d");
    type="table3d";
      
    return 0;
  }
    
  ret=hf.find_object_by_type("hist",in[1],use_regex,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, found first hist object named '"
	   << in[1] << "'." << endl;
    }
    hdf_input(hf,hist_obj,in[1]);
    obj_name=in[1];
    command_add("hist");
    type="hist";
    return 0;
  }
    
  ret=hf.find_object_by_type("hist_2d",in[1],use_regex,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, found first hist_2d object named '"
	   << in[1] << "'." << endl;
    }
    hdf_input(hf,hist_2d_obj,in[1]);
    obj_name=in[1];
    command_add("hist_2d");
    type="hist_2d";
    return 0;
  }
  
  ret=hf.find_object_by_type("tensor_grid",in[1],use_regex,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, found first tensor_grid object named '"
	   << in[1] << "'." << endl;
    }
    hdf_input(hf,tensor_grid_obj,in[1]);
    obj_name=in[1];
    command_add("tensor_grid");
    type="tensor_grid";
    return 0;
  }
  
  ret=hf.find_object_by_type("tensor",in[1],use_regex,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, found first tensor object named '"
	   << in[1] << "'." << endl;
    }
    hf.getd_ten(in[1],tensor_obj);
    obj_name=in[1];
    command_add("tensor");
    type="tensor";
    return 0;
  }
  
  ret=hf.find_object_by_type("tensor<size_t>",in[1],use_regex,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, found first tensor<size_t> object named '"
	   << in[1] << "'." << endl;
    }
    hf.get_szt_ten(in[1],tensor_size_t_obj);
    obj_name=in[1];
    command_add("tensor<size_t>");
    type="tensor<size_t>";
    return 0;
  }
  
  ret=hf.find_object_by_type("tensor<int>",in[1],use_regex,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, found first tensor<int> object named '"
	   << in[1] << "'." << endl;
    }
    hf.geti_ten(in[1],tensor_int_obj);
    obj_name=in[1];
    command_add("tensor<int>");
    type="tensor<int>";
    return 0;
  }
  
  ret=hf.find_object_by_type("vector<contour_line>",in[1],use_regex,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, found first vector<contour_line> "
	   << "object named '"
	   << in[1] << "'." << endl;
    }
    hdf_input(hf,cont_obj,in[1]);
    obj_name=in[1];
    command_add("vector<contour_line>");
    type="vector<contour_line>";
    return 0;
  }
  
  ret=hf.find_object_by_type("uniform_grid<double>",in[1],use_regex,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, found first uniform_grid<double> "
	   << "object named '"
	   << in[1] << "'." << endl;
    }
    hdf_input(hf,ug_obj,in[1]);
    obj_name=in[1];
    command_add("uniform_grid<double>");
    type="uniform_grid<double>";
    return 0;
  }
  
  ret=hf.find_object_by_type("prob_dens_mdim_amr",in[1],use_regex,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, found first prob_dens_mdim_amr "
	   << "object named '"
	   << in[1] << "'." << endl;
    }
    hdf_input(hf,pdma_obj,in[1]);
    obj_name=in[1];
    command_add("prob_dens_mdim_amr");
    type="prob_dens_mdim_amr";
    return 0;
  }

  ret=hf.find_object_by_type("double[]",in[1],use_regex,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, found first double[] "
	   << "object named '"
	   << in[1] << "'." << endl;
    }
    hf.getd_vec(in[1],doublev_obj);
    obj_name=in[1];
    command_add("double[]");
    type="double[]";
    return 0;
  }
  
  ret=hf.find_object_by_type("int[]",in[1],use_regex,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, found first int[] "
	   << "object named '"
	   << in[1] << "'." << endl;
    }
    hf.geti_vec(in[1],intv_obj);
    obj_name=in[1];
    command_add("int[]");
    type="int[]";
    return 0;
  }
  
  ret=hf.find_object_by_type("string[]",in[1],use_regex,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, found first string[] "
	   << "object named '"
	   << in[1] << "'." << endl;
    }
    hf.gets_vec(in[1],stringv_obj);
    obj_name=in[1];
    command_add("string[]");
    type="string[]";
    return 0;
  }
  
  ret=hf.find_object_by_type("size_t[]",in[1],use_regex,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, found first size_t[] "
	   << "object named '"
	   << in[1] << "'." << endl;
    }
    hf.get_szt_vec(in[1],size_tv_obj);
    obj_name=in[1];
    command_add("size_t[]");
    type="size_t[]";
    return 0;
  }
  
  ret=hf.find_object_by_type("double",in[1],use_regex,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, found first double "
	   << "object named '"
	   << in[1] << "'." << endl;
    }
    hf.getd(in[1],double_obj);
    obj_name=in[1];
    command_add("double");
    type="double";
    return 0;
  }
  
  ret=hf.find_object_by_type("int",in[1],use_regex,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, found first int "
	   << "object named '"
	   << in[1] << "'." << endl;
    }
    hf.geti(in[1],int_obj);
    obj_name=in[1];
    command_add("int");
    type="int";
    return 0;
  }
  
  ret=hf.find_object_by_type("string",in[1],use_regex,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, found first string "
	   << "object named '"
	   << in[1] << "'." << endl;
    }
    hf.gets(in[1],string_obj);
    obj_name=in[1];
    command_add("string");
    type="string";
    return 0;
  }
  
  ret=hf.find_object_by_type("size_t",in[1],use_regex,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, found first size_t "
	   << "object named '"
	   << in[1] << "'." << endl;
    }
    hf.get_szt(in[1],size_t_obj);
    obj_name=in[1];
    command_add("size_t");
    type="size_t";
    return 0;
  }
  
  ret=hf.find_object_by_type("char",in[1],use_regex,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, found first char "
	   << "object named '"
	   << in[1] << "'." << endl;
    }
    hf.getc(in[1],char_obj);
    obj_name=in[1];
    command_add("char");
    type="char";
    return 0;
  }
  
  cout << "Could not find object of any readable type in file '" << fname
       << "'." << endl;
  
  return exc_efailed;
}

int acol_manager::comm_rearrange(std::vector<std::string> &sv,
				 bool itive_com) {

  if (type=="tensor") {

    vector<string> sv2;
    for(size_t j=1;j<sv.size();j++) {
      tensor_obj.index_spec_preprocess(sv[j],sv2);
    }
    
    vector<o2scl::index_spec> vis;
    tensor_obj.strings_to_indexes(sv2,vis,verbose);
    
    tensor<> t;
    t=tensor_obj.rearrange_and_copy(vis,verbose,false);
    if (t.total_size()==0) {
      cerr << "Function rearrange_and_copy() failed." << endl;
      return 1;
    }
    tensor_obj=t;
    
  } else if (type=="tensor<int>") {
    
    vector<string> sv2;
    for(size_t j=1;j<sv.size();j++) {
      tensor_int_obj.index_spec_preprocess(sv[j],sv2);
    }
    
    vector<o2scl::index_spec> vis;
    tensor_int_obj.strings_to_indexes(sv2,vis,verbose);
    
    tensor<int> t;
    t=tensor_int_obj.rearrange_and_copy(vis,verbose,false);
    if (t.total_size()==0) {
      cerr << "Function rearrange_and_copy() failed." << endl;
      return 1;
    }
    tensor_int_obj=t;
    
  } else if (type=="tensor<size_t>") {

    vector<string> sv2;
    for(size_t j=1;j<sv.size();j++) {
      tensor_size_t_obj.index_spec_preprocess(sv[j],sv2);
    }
    
    vector<o2scl::index_spec> vis;
    tensor_size_t_obj.strings_to_indexes(sv2,vis,verbose);
    
    tensor<size_t> t;
    t=tensor_size_t_obj.rearrange_and_copy(vis,verbose,false);
    if (t.total_size()==0) {
      cerr << "Function rearrange_and_copy() failed." << endl;
      return 1;
    }
    tensor_size_t_obj=t;

  } else if (type=="tensor_grid") {

    vector<string> sv2;
    for(size_t j=1;j<sv.size();j++) {
      tensor_grid_obj.index_spec_preprocess(sv[j],sv2);
    }

    vector<o2scl::index_spec> vis;
    tensor_grid_obj.strings_to_indexes(sv2,vis,verbose);
    
    tensor_grid<> t;
    t=tensor_grid_obj.rearrange_and_copy(vis,verbose,false);
    if (t.total_size()==0) {
      cerr << "Function rearrange_and_copy() failed." << endl;
      return 1;
    }
    tensor_grid_obj=t;
    
  } else {
    cerr << "Rearrange does not work with " << type << " objects." << endl;
  }
    
  return 0;
}

int acol_manager::comm_sum(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table3d") {
    
    if (sv.size()<2) {
      cerr << "Not enough arguments to sum." << endl;
      return exc_efailed;
    }

    string s2=sv[1];
    string name2;
    if (sv.size()>=3) name2=sv[2];

    table3d t2;

    hdf_file hf;
    int hfret=hf.open(s2,false,false);
    if (hfret!=0) {
      cerr << "Failed to read file named " << s2 << endl;
      return exc_efailed;
    }
    hdf_input(hf,t2,name2);
    hf.close();
  
    size_t nx=t2.get_nx();
    size_t ny=t2.get_ny();
    const ubvector &xg=t2.get_x_data();
    const ubvector &yg=t2.get_y_data();
  
    for(size_t k=0;k<t2.get_nslices();k++) {
      string slname=t2.get_slice_name(k);
      size_t slix;
      if (table3d_obj.is_slice(slname,slix)==false) {
	table3d_obj.new_slice(slname);
	table3d_obj.set_slice_all(slname,0.0);
      }
      for(size_t i=0;i<nx;i++) {
	for(size_t j=0;j<ny;j++) {
	  double x=xg[i];
	  double y=yg[j];
	  double value=table3d_obj.get_val(x,y,slname)+t2.get(i,j,slname);
	  table3d_obj.set_val(x,y,slname,value);
	}
      }
    }

  } else if (type=="table") {

    if (sv.size()<2) {
      cerr << "Not enough arguments to add." << endl;
      return exc_efailed;
    }
    if (table_obj.get_nlines()==0) {
      cerr << "No table to add to." << endl;
      return exc_efailed;
    }

    string file2=sv[1];
    std::string name2;
    if (sv.size()>=3) name2=sv[2];
    table_units<> tab2;

    hdf_file hf;
    int hfret=hf.open(file2,false,false);
    if (hfret!=0) {
      cerr << "Failed to read file named " << file2 << endl;
      return exc_efailed;
    }
    hdf_input(hf,tab2,name2);
    hf.close();

    size_t n1=table_obj.get_nlines();
    size_t n2=tab2.get_nlines();
    if (n2>n1) {
      table_obj.set_nlines(n1+n2);
      for(size_t j=0;j<table_obj.get_ncolumns();j++) {
	for(size_t i=n1;i<n1+n2;i++) {
	  table_obj.set(j,i,0.0);
	}
      }
    }
    for(size_t j=0;j<tab2.get_ncolumns();j++) {
      std::string col_name=tab2.get_column_name(j);
      if (!table_obj.is_column(col_name)) {
	table_obj.new_column(col_name);
	for(size_t i=0;i<table_obj.get_nlines();i++) {
	  table_obj.set(col_name,i,0.0);
	}
      }
      for(size_t i=0;i<n2;i++) {
	table_obj.set(col_name,i,tab2.get(col_name,i)+
		      table_obj.get(col_name,i));
      }
    }
    
  } else if (type=="tensor_grid") {
    const std::vector<double> &data=tensor_grid_obj.get_data();
    cout << "Sum of all " << data.size() << " entries in the tensor is: "
	 << o2scl::vector_sum<const vector<double>,double>(data) << endl;
  } else if (type=="tensor") {
    const std::vector<double> &data=tensor_obj.get_data();
    cout << "Sum of all " << data.size() << " entries in the tensor is: "
	 << o2scl::vector_sum<const vector<double>,double>(data) << endl;
  } else if (type=="double[]") {
    cout << "Sum of " << doublev_obj.size() << " elements is: "
	 << o2scl::vector_sum<vector<double>,double>(doublev_obj) << endl;
  } else if (type=="int[]") {
    cout << "Sum of " << intv_obj.size() << " elements is: "
	 << o2scl::vector_sum<vector<int>,int>(intv_obj) << endl;
  } else if (type=="size_t[]") {
    cout << "Sum of " << size_tv_obj.size() << " elements is: "
	 << o2scl::vector_sum<vector<size_t>,size_t>(size_tv_obj) << endl;
  } else {

    cerr << "Cannot 'sum' with object of type " << type << endl;
    return exc_efailed;
    
  }
  
  return 0;
}

