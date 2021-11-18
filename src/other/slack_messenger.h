/*
  -------------------------------------------------------------------

  Copyright (C) 2019-2021, Andrew W. Steiner

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
/** \file slack_messenger.h
    \brief File defining \ref o2scl::slack_messenger
*/
#ifndef O2SCL_SLACK_MESSENGER_H
#define O2SCL_SLACK_MESSENGER_H

#include <string>
#include <vector>

#include <o2scl/err_hnd.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Object to send messages to Slack using curl

      \note Experimental
   */
  class slack_messenger {

  protected:

    /** \brief Time (in seconds) the last message was sent
     */
    double time_last_message;

    /** \brief If true, use MPI to determine time
     */
    bool mpi_time;
  
  public:

    /** \brief The URL for the Slack webhook
     */
    std::string url;

    /** \brief The destination channel
     */
    std::string channel;

    /** \brief Minimum time between messages in seconds (default 300)
     */
    double min_time_between;

    /** \brief Icon to use (without colons; default "computer")
     */
    std::string icon;

    /** \brief Slack username 
     */
    std::string username;

    /** \brief Verbosity parameter (default 1)
     */
    int verbose;

    /** \brief Create a messenger object with specified channel,
	username, URL and time method

	If the URL is not specified, this constructor tries to determine
	it from the environment variable \c O2SCL_SLACK_URL .
    */
    slack_messenger(std::string p_channel="", std::string p_username="", 
		    std::string p_url="", bool p_mpi_time=false) {
    
      if (p_url.length()==0) {
	set_url_from_env("O2SCL_SLACK_URL");
      } else {
	url=p_url;
      } 
      channel=p_channel;
      username=p_username;
      min_time_between=300.0;
      icon="computer";
      verbose=1;
      mpi_time=p_mpi_time;
      if (mpi_time) {
#ifdef O2SCL_MPI
	time_last_message=MPI_Wtime()-min_time_between-1.0;
#else
	O2SCL_ERR2("Value mpi_time is true but O2SCL_MPI not defined ",
		   "in slack_messenger::slack_messenger().",
		   o2scl::exc_einval);
#endif
      } else {
	time_last_message=time(0)-min_time_between-1.0;
      }
    }

    /** \brief Set the time mode (normal or MPI)
     */
    void set_time_mode(bool loc_mpi_time) {
      if (mpi_time) {
#ifdef O2SCL_MPI
	time_last_message=MPI_Wtime();
#else
	O2SCL_ERR2("Value mpi_time is true but O2SCL_MPI not defined ",
		   "in slack_message::slack_message().",
		   o2scl::exc_einval);
#endif
      } else {
	time_last_message=time(0);
      }
      mpi_time=loc_mpi_time;
      return;
    }
  
    /** \brief Set the Slack webhook URL from the environment variable
        \c env_var
     */
    bool set_url_from_env(std::string env_var) {
      char *cstring=getenv(env_var.c_str());
      if (cstring) {
	url=cstring;
	return true;
      }
      return false;
    }
  
    /** \brief Set the channel from the environment variable
        \c env_var
     */
    bool set_channel_from_env(std::string env_var) {
      char *cstring=getenv(env_var.c_str());
      if (cstring) {
	channel=cstring;
	return true;
      }
      return false;
    }

    /** \brief Set the username from the environment variable
        \c env_var
     */
    bool set_username_from_env(std::string env_var) {
      char *cstring=getenv(env_var.c_str());
      if (cstring) {
	username=cstring;
	return true;
      }
      return false;
    }
  
    /** \brief Send a message
     */
    int send(std::string message, bool err_on_fail=true) {
      
      int iret=0;

      if (url.length()>0) {

	if (channel.length()==0) {
	  O2SCL_ERR2("No slack channel specified in ",
		     "slack_messenger::send().",o2scl::exc_einval);
	}
	if (username.length()==0) {
	  O2SCL_ERR2("No slack username specified in ",
		     "slack_messenger::send().",o2scl::exc_einval);
	}

	double time_now=time_last_message;
	if (mpi_time) {
#ifdef O2SCL_MPI
	  time_now=MPI_Wtime();
#else
	  O2SCL_ERR2("Value mpi_time is true but O2SCL_MPI not defined ",
		     "in slack_message::slack_message().",
		    o2scl::exc_einval);
#endif
	} else {
	  time_now=time(0);
	}
      
	if (time_now-time_last_message>min_time_between) {
	
	  std::string scr;
	  if (icon.length()>0) {
	    scr=((std::string)"curl -X POST --data-urlencode ")+
	      "\"payload={\\\"channel\\\": \\\""+channel+"\\\", "+
	      "\\\"username\\\": \\\""+username+"\\\", "+
	      "\\\"text\\\": \\\""+message+"\\\", "+
	      "\\\"icon_emoji\\\": \\\":"+icon+":\\\"}\" "+url;
	  } else {
	    scr=((std::string)"curl -X POST --data-urlencode ")+
	      "\"payload={\\\"channel\\\": \\\""+channel+"\\\", "+
	      "\\\"username\\\": \\\""+username+"\\\", "+
	      "\\\"text\\\": \\\""+message+"\\\"}\" "+url;
	  }

	  if (verbose>1) {
	    std::cout << "Executing: " << scr << std::endl;
	  }
	
	  iret=system(scr.c_str());
	
	  if (iret!=0 && err_on_fail) {
	    O2SCL_ERR2("System command failed in ",
		       "slack_messenger::send().",o2scl::exc_efailed);
	  }
	
	  time_last_message=time_now;
	
	}
      
      } else {
	O2SCL_ERR2("No slack URL specified in ",
		   "slack_messenger::send().",o2scl::exc_einval);
      }
    
      return iret;
    }

    /** \brief Send a message
     */
    int send_image(std::string message, std::string image_url,
                   std::string alt_text, bool err_on_fail=true) {
      
      int iret=0;

      if (url.length()>0) {

	if (channel.length()==0) {
	  O2SCL_ERR2("No slack channel specified in ",
		     "slack_messenger::send().",o2scl::exc_einval);
	}
	if (username.length()==0) {
	  O2SCL_ERR2("No slack username specified in ",
		     "slack_messenger::send().",o2scl::exc_einval);
	}

	double time_now=time_last_message;
	if (mpi_time) {
#ifdef O2SCL_MPI
	  time_now=MPI_Wtime();
#else
	  O2SCL_ERR2("Value mpi_time is true but O2SCL_MPI not defined ",
		     "in slack_message::slack_message().",
		    o2scl::exc_einval);
#endif
	} else {
	  time_now=time(0);
	}
      
	if (time_now-time_last_message>min_time_between) {
	
	  std::string scr;
	  if (icon.length()>0) {
	    scr=((std::string)"curl -X POST --data-urlencode ")+
	      "\"payload={"+
              "\\\"channel\\\": \\\""+channel+"\\\", "+
	      "\\\"username\\\": \\\""+username+"\\\", "+
	      "\\\"icon_emoji\\\": \\\":"+icon+":\\\", "+ 
	      "\\\"blocks\\\": [ { \\\"type\\\": \\\"image\\\", "+
              "\\\"title\\\": { \\\"type\\\": \\\"plain_text\\\", "+
              "\\\"text\\\": \\\""+message+"\\\" }, "+
              "\\\"image_url\\\": \\\""+image_url+"\\\", "+
              "\\\"alt_text\\\": \\\""+alt_text+"\\\" "+
              "} ] }\" "+url;
          } else {
	    scr=((std::string)"curl -X POST --data-urlencode ")+
	      "\"payload={"+
              "\\\"channel\\\": \\\""+channel+"\\\", "+
	      "\\\"username\\\": \\\""+username+"\\\", "+
	      "\\\"blocks\\\": [ { \\\"type\\\": \\\"image\\\", "+
              "\\\"title\\\": { \\\"type\\\": \\\"plain_text\\\", "+
              "\\\"text\\\": \\\""+message+"\\\" }, "+
              "\\\"image_url\\\": \\\""+image_url+"\\\", "+
              "\\\"alt_text\\\": \\\""+alt_text+"\\\" "+
              "} ] }\" "+url;
	  }

	  if (verbose>1) {
	    std::cout << "Executing: " << scr << std::endl;
	  }
	
	  iret=system(scr.c_str());
	
	  if (iret!=0 && err_on_fail) {
	    O2SCL_ERR2("System command failed in ",
		       "slack_messenger::send().",o2scl::exc_efailed);
	  }
	
	  time_last_message=time_now;
	
	}
      
      } else {
	O2SCL_ERR2("No slack URL specified in ",
		   "slack_messenger::send().",o2scl::exc_einval);
      }
    
      return iret;
    }
    
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif

