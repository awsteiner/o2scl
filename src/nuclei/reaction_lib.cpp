/*
  -------------------------------------------------------------------
  
  Copyright (C) 2008-2019, Andrew W. Steiner
  
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
#include <o2scl/reaction_lib.h>

using namespace std;
using namespace o2scl;

bool reaction_lib::matches(size_t ul, size_t ri) {
  if (fZ[ul]==((int)lib[fi].Z[ri]) &&
      fA[ul]==((int)lib[fi].A[ri])) {
    return true;
  }
  return false;
}

int reaction_lib::find_in_chap(std::vector<nuclear_reaction> &nrl,
			       size_t chap, std::string nuc1, std::string nuc2, 
			       std::string nuc3, std::string nuc4, 
			       std::string nuc5, std::string nuc6) {
  
  nucmass_info nmi;
  size_t nspec=1;
  nmi.parse_elstring(nuc1,fZ[0],fN[0],fA[0]);
  if (nuc2.length()>0) {
    nmi.parse_elstring(nuc2,fZ[1],fN[1],fA[1]);
    nspec++;
    if (nuc3.length()>0) {
      nmi.parse_elstring(nuc3,fZ[2],fN[2],fA[2]);
      nspec++;
      if (nuc4.length()>0) {
	nmi.parse_elstring(nuc4,fZ[3],fN[3],fA[3]);
	nspec++;
	if (nuc5.length()>0) {
	  nmi.parse_elstring(nuc5,fZ[4],fN[4],fA[4]);
	  nspec++;
	  if (nuc6.length()>0) {
	    nmi.parse_elstring(nuc6,fZ[5],fN[5],fA[5]);
	    nspec++;
	  }
	}
      }
    }
  }

  if (chap==1) {
    // nuc1 -> nuc2

    if (nspec==1) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && matches(0,0)) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else if (nspec==2) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && matches(0,0) && matches(1,1)) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else {
      O2SCL_ERR2("Too many nuclei specified for chapter in ",
		     "reaction_lib::find_in_chap().",exc_efailed);
    }

  } else if (chap==2) {
    // nuc1 -> nuc2 + nuc3

    if (nspec==1) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && matches(0,0)) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else if (nspec==2) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && 
	    ((matches(0,0) && matches(1,1)) ||
	     (matches(0,0) && matches(1,2)))) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else if (nspec==3) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && 
	    ((matches(0,0) && matches(1,1) && matches(2,2)) ||
	     (matches(0,0) && matches(1,2) && matches(2,1)))) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else {
      O2SCL_ERR2("Too many nuclei specified for chapter in ",
		     "reaction_lib::find_in_chap().",exc_efailed);
    }

  } else if (chap==3) {
    // nuc1 -> nuc2 + nuc3 + nuc4

    if (nspec==1) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && matches(0,0)) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else if (nspec==2) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && 
	    ((matches(0,0) && matches(1,1)) ||
	     (matches(0,0) && matches(1,2)) ||
	     (matches(0,0) && matches(1,3)))) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else if (nspec==3) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && 
	    ((matches(0,0) && matches(1,1) && matches(2,2)) ||
	     (matches(0,0) && matches(1,2) && matches(2,3)) ||
	     (matches(0,0) && matches(1,1) && matches(2,3)) ||
	     (matches(0,0) && matches(1,2) && matches(2,1)) ||
	     (matches(0,0) && matches(1,3) && matches(2,2)) ||
	     (matches(0,0) && matches(1,3) && matches(2,1)))) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else if (nspec==4) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && 
	    ((matches(0,0) && matches(1,1) && 
	      matches(2,2) && matches(3,3)) ||
	     (matches(0,0) && matches(1,1) && 
	      matches(2,3) && matches(3,2)) ||
	     (matches(0,0) && matches(1,2) && 
	      matches(2,1) && matches(3,3)) ||
	     (matches(0,0) && matches(1,2) && 
	      matches(2,3) && matches(3,1)) ||
	     (matches(0,0) && matches(1,3) && 
	      matches(2,1) && matches(3,2)) ||
	     (matches(0,0) && matches(1,3) && 
	      matches(2,2) && matches(3,1)))) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else {
      O2SCL_ERR2("Too many nuclei specified for chapter in ",
		     "reaction_lib::find_in_chap().",exc_efailed);
    }

  } else if (chap==4) {
    // nuc1 + nuc2 -> nuc3

    if (nspec==1) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && (matches(0,0) || matches(0,1))) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else if (nspec==2) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && 
	    ((matches(0,0) && matches(1,1)) ||
	     (matches(0,1) && matches(1,0)))) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else if (nspec==3) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && 
	    ((matches(0,0) && matches(1,1) && matches(2,2)) ||
	     (matches(0,1) && matches(1,0) && matches(2,2)))) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else {
      O2SCL_ERR2("Too many nuclei specified for chapter in ",
		     "reaction_lib::find_in_chap().",exc_efailed);
    }

  } else if (chap==5) {
    // nuc1 + nuc2 -> nuc3 + nuc4

    if (nspec==1) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && (matches(0,0) || matches(0,1))) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else if (nspec==2) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && 
	    ((matches(0,0) && matches(1,1)) ||
	     (matches(0,1) && matches(1,0)))) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else if (nspec==3) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && 
	    ((matches(0,0) && matches(1,1) && matches(2,2)) ||
	     (matches(0,0) && matches(1,1) && matches(2,3)) ||
	     (matches(0,1) && matches(1,0) && matches(2,2)) ||
	     (matches(0,1) && matches(1,0) && matches(2,3)))) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else if (nspec==4) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && 
	    ((matches(0,0) && matches(1,1) && 
	      matches(2,2) && matches(3,3)) ||
	     (matches(0,0) && matches(1,1) && 
	      matches(2,3) && matches(3,2)) ||
	     (matches(0,1) && matches(1,0) && 
	      matches(2,2) && matches(3,3)) ||
	     (matches(0,1) && matches(1,0) && 
	      matches(2,3) && matches(3,2)))) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else {
      O2SCL_ERR2("Too many nuclei specified for chapter in ",
		     "reaction_lib::find_in_chap().",exc_efailed);
    }

  } else if (chap==6) {
    // nuc1 + nuc2 -> nuc3 + nuc4 + nuc5 

    if (nspec==1) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && (matches(0,0) || matches(0,1))) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else if (nspec==2) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && 
	    ((matches(0,0) && matches(1,1)) ||
	     (matches(0,1) && matches(1,0)))) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else if (nspec==3) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap) {
	  if ((matches(0,0) && matches(1,1)) ||
	      (matches(1,0) && matches(0,1))) {
	    if (matches(2,2) || matches(2,3) || matches(2,4)) {
	      nrl.push_back(lib[fi]);
	    }
	  }
	}
      }
    } else if (nspec==4) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap) {
	  if ((matches(0,0) && matches(1,1)) ||
	      (matches(1,0) && matches(0,1))) {
	    if ((matches(2,2) && matches(3,3)) ||
		(matches(2,2) && matches(3,4)) ||
		(matches(2,3) && matches(3,2)) ||
		(matches(2,3) && matches(3,4)) ||
		(matches(2,4) && matches(3,2)) ||
		(matches(2,4) && matches(3,3))) {
	      nrl.push_back(lib[fi]);
	    }
	  }
	}
      }
    } else if (nspec==5) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap) {
	  if ((matches(0,0) && matches(1,1)) ||
	      (matches(1,0) && matches(0,1))) {
	    if ((matches(2,2) && matches(3,3) && matches(4,4)) ||
		(matches(2,2) && matches(3,4) && matches(4,3)) ||
		(matches(2,3) && matches(3,2) && matches(4,4)) ||
		(matches(2,3) && matches(3,4) && matches(4,2)) ||
		(matches(2,4) && matches(3,2) && matches(4,3)) ||
		(matches(2,4) && matches(3,3) && matches(4,2))) {
	      nrl.push_back(lib[fi]);
	    }
	  }
	}
      }
    } else {
      O2SCL_ERR2("Too many nuclei specified for chapter in ",
		     "reaction_lib::find_in_chap().",exc_efailed);
    }

  } else if (chap==7) {
    // nuc1 + nuc2 -> nuc3 + nuc4 + nuc5 + nuc6

    if (nspec==1) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && (matches(0,0) || matches(0,1))) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else if (nspec==2) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && 
	    ((matches(0,0) && matches(1,1)) ||
	     (matches(0,1) && matches(1,0)))) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else if (nspec==3) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap) {
	  if ((matches(0,0) && matches(1,1)) ||
	      (matches(1,0) && matches(0,1))) {
	    if (matches(2,2) || matches(2,3) || 
		matches(2,4) || matches(2,5)) {
	      nrl.push_back(lib[fi]);
	    }
	  }
	}
      }
    } else if (nspec==4) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap) {
	  if ((matches(0,0) && matches(1,1)) ||
	      (matches(1,0) && matches(0,1))) {
	    if ((matches(2,2) && matches(3,3)) ||
		(matches(2,2) && matches(3,4)) ||
		(matches(2,2) && matches(3,5)) ||
		(matches(2,3) && matches(3,2)) ||
		(matches(2,3) && matches(3,4)) ||
		(matches(2,3) && matches(3,5)) ||
		(matches(2,4) && matches(3,2)) ||
		(matches(2,4) && matches(3,3)) ||
		(matches(2,4) && matches(3,5)) ||
		(matches(2,5) && matches(3,2)) ||
		(matches(2,5) && matches(3,3)) ||
		(matches(2,5) && matches(3,4))) {
	      nrl.push_back(lib[fi]);
	    }
	  }
	}
      }
    } else if (nspec==5) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap) {
	  if ((matches(0,0) && matches(1,1)) ||
	      (matches(1,0) && matches(0,1))) {
	    if (matches(2,2)) {
	      if ((matches(3,3) && matches(4,4)) ||
		  (matches(3,3) && matches(4,5)) ||
		  (matches(3,5) && matches(4,4)) ||
		  (matches(3,4) && matches(4,3)) ||
		  (matches(3,4) && matches(4,5)) ||
		  (matches(3,5) && matches(4,3))) {
		nrl.push_back(lib[fi]);
	      }
	    } else if (matches(2,3)) {
	      if ((matches(3,2) && matches(4,4)) ||
		  (matches(3,2) && matches(4,5)) ||
		  (matches(3,5) && matches(4,4)) ||
		  (matches(3,4) && matches(4,2)) ||
		  (matches(3,4) && matches(4,5)) ||
		  (matches(3,5) && matches(4,2))) {
		nrl.push_back(lib[fi]);
	      }
	    } else if (matches(2,4)) {
	      if ((matches(3,3) && matches(4,2)) ||
		  (matches(3,3) && matches(4,5)) ||
		  (matches(3,5) && matches(4,2)) ||
		  (matches(3,2) && matches(4,3)) ||
		  (matches(3,2) && matches(4,5)) ||
		  (matches(3,5) && matches(4,3))) {
		nrl.push_back(lib[fi]);
	      }
	    } else if (matches(2,5)) {
	      if ((matches(3,3) && matches(4,4)) ||
		  (matches(3,3) && matches(4,2)) ||
		  (matches(3,2) && matches(4,4)) ||
		  (matches(3,4) && matches(4,3)) ||
		  (matches(3,4) && matches(4,2)) ||
		  (matches(3,2) && matches(4,3))) {
		nrl.push_back(lib[fi]);
	      }
	    }
	  }
	}
      }
    } else if (nspec==6) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap) {
	  if ((matches(0,0) && matches(1,1)) ||
	      (matches(1,0) && matches(0,1))) {
	    if (matches(2,2)) {
	      if (matches(3,3))  {
		if ((matches(4,4) && matches(5,5)) ||
		    (matches(4,5) && matches(5,4))) {
		  nrl.push_back(lib[fi]);
		}
	      } else if (matches(3,4)) {
		if ((matches(4,3) && matches(5,5)) ||
		    (matches(4,5) && matches(5,3))) {
		  nrl.push_back(lib[fi]);
		}
	      } else if (matches(3,5)) {
		if ((matches(4,3) && matches(5,4)) ||
		    (matches(4,4) && matches(5,3))) {
		  nrl.push_back(lib[fi]);
		}
	      }
	    } else if (matches(2,3)) {
	      if (matches(3,2))  {
		if ((matches(4,4) && matches(5,5)) ||
		    (matches(4,5) && matches(5,4))) {
		  nrl.push_back(lib[fi]);
		}
	      } else if (matches(3,4)) {
		if ((matches(4,2) && matches(5,5)) ||
		    (matches(4,5) && matches(5,2))) {
		  nrl.push_back(lib[fi]);
		}
	      } else if (matches(3,5)) {
		if ((matches(4,4) && matches(5,2)) ||
		    (matches(4,2) && matches(5,4))) {
		  nrl.push_back(lib[fi]);
		}
	      }
	    } else if (matches(2,4)) {
	      if (matches(3,2))  {
		if ((matches(4,3) && matches(5,5)) ||
		    (matches(4,5) && matches(5,3))) {
		  nrl.push_back(lib[fi]);
		}
	      } else if (matches(3,3)) {
		if ((matches(4,2) && matches(5,5)) ||
		    (matches(4,5) && matches(5,2))) {
		  nrl.push_back(lib[fi]);
		}
	      } else if (matches(3,5)) {
		if ((matches(4,2) && matches(5,2)) ||
		    (matches(4,3) && matches(5,3))) {
		  nrl.push_back(lib[fi]);
		}
	      }
	    } else if (matches(2,5)) {
	      if (matches(3,2))  {
		if ((matches(4,3) && matches(5,4)) ||
		    (matches(4,4) && matches(5,3))) {
		  nrl.push_back(lib[fi]);
		}
	      } else if (matches(3,3)) {
		if ((matches(4,2) && matches(5,4)) ||
		    (matches(4,4) && matches(5,2))) {
		  nrl.push_back(lib[fi]);
		}
	      } else if (matches(3,4)) {
		if ((matches(4,2) && matches(5,3)) ||
		    (matches(4,3) && matches(5,2))) {
		  nrl.push_back(lib[fi]);
		}
	      }
	    }
	  }
	}
      }
    } else {
      O2SCL_ERR2("Too many nuclei specified for chapter in ",
		     "reaction_lib::find_in_chap().",exc_efailed);
    }
    
  } else if (chap==8) {
    // nuc1 + nuc2 + nuc3 -> nuc4
    
    if (nspec==1) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && (matches(0,0) || matches(0,1) || 
				   matches(0,2))) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else if (nspec==2) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && 
	    ((matches(0,0) && matches(1,1)) ||
	     (matches(0,0) && matches(1,2)) ||
	     (matches(0,1) && matches(1,0)) ||
	     (matches(0,1) && matches(1,2)) ||
	     (matches(0,2) && matches(1,0)) ||
	     (matches(0,2) && matches(1,1)))) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else if (nspec==3) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap) {
	  if ((matches(0,0) && matches(1,1) && matches(2,2)) ||
	      (matches(0,0) && matches(1,2) && matches(2,1)) ||
	      (matches(0,1) && matches(1,0) && matches(2,2)) ||
	      (matches(0,1) && matches(1,2) && matches(2,0)) ||
	      (matches(0,2) && matches(1,0) && matches(2,1)) ||
	      (matches(0,2) && matches(1,1) && matches(2,0))) {
	    nrl.push_back(lib[fi]);
	  }
	}
      }
    } else if (nspec==4) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && matches(3,3)) {
	  if ((matches(0,0) && matches(1,1) && matches(2,2)) ||
	      (matches(0,0) && matches(1,2) && matches(2,1)) ||
	      (matches(0,1) && matches(1,0) && matches(2,2)) ||
	      (matches(0,1) && matches(1,2) && matches(2,0)) ||
	      (matches(0,2) && matches(1,0) && matches(2,1)) ||
	      (matches(0,2) && matches(1,1) && matches(2,0))) {
	    nrl.push_back(lib[fi]);
	  }
	}
      }
    } else {
      O2SCL_ERR2("Too many nuclei specified for chapter in ",
		     "reaction_lib::find_in_chap().",exc_efailed);
    }

  } else if (chap==9) {
    // nuc1 + nuc2 + nuc3 -> nuc4 + nuc5
    
    if (nspec==1) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && (matches(0,0) || matches(0,1) || 
				   matches(0,2))) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else if (nspec==2) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && 
	    ((matches(0,0) && matches(1,1)) ||
	     (matches(0,0) && matches(1,2)) ||
	     (matches(0,1) && matches(1,0)) ||
	     (matches(0,1) && matches(1,2)) ||
	     (matches(0,2) && matches(1,0)) ||
	     (matches(0,2) && matches(1,1)))) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else if (nspec==3) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap) {
	  if ((matches(0,0) && matches(1,1) && matches(2,2)) ||
	      (matches(0,0) && matches(1,2) && matches(2,1)) ||
	      (matches(0,1) && matches(1,0) && matches(2,2)) ||
	      (matches(0,1) && matches(1,2) && matches(2,0)) ||
	      (matches(0,2) && matches(1,0) && matches(2,1)) ||
	      (matches(0,2) && matches(1,1) && matches(2,0))) {
	    nrl.push_back(lib[fi]);
	  }
	}
      }
    } else if (nspec==4) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap) {
	  if (matches(3,3) || matches(3,4)) {
	    if ((matches(0,0) && matches(1,1) && matches(2,2)) ||
		(matches(0,0) && matches(1,2) && matches(2,1)) ||
		(matches(0,1) && matches(1,0) && matches(2,2)) ||
		(matches(0,1) && matches(1,2) && matches(2,0)) ||
		(matches(0,2) && matches(1,0) && matches(2,1)) ||
		(matches(0,2) && matches(1,1) && matches(2,0))) {
	      nrl.push_back(lib[fi]);
	    }
	  }
	}
      }
    } else if (nspec==5) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap) {
	  if ((matches(3,3) && matches(4,4)) ||
	      (matches(3,4) && matches(4,3))) {
	    if ((matches(0,0) && matches(1,1) && matches(2,2)) ||
		(matches(0,0) && matches(1,2) && matches(2,1)) ||
		(matches(0,1) && matches(1,0) && matches(2,2)) ||
		(matches(0,1) && matches(1,2) && matches(2,0)) ||
		(matches(0,2) && matches(1,0) && matches(2,1)) ||
		(matches(0,2) && matches(1,1) && matches(2,0))) {
	      nrl.push_back(lib[fi]);
	    }
	  }
	}
      }
    } else {
      O2SCL_ERR2("Too many nuclei specified for chapter in ",
		     "reaction_lib::find_in_chap().",exc_efailed);
    }

  } else if (chap==10) {
    // nuc1 + nuc2 + nuc3 + nuc4 -> nuc5 + nuc6
    
    if (nspec==1) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && 
	    (matches(0,0) || matches(0,1) || matches(0,2) || matches(0,3))) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else if (nspec==2) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap) {
	  if ((matches(0,0) && matches(1,1)) ||
	      (matches(0,0) && matches(1,2)) ||
	      (matches(0,0) && matches(1,3)) ||
	      (matches(0,1) && matches(1,0)) ||
	      (matches(0,1) && matches(1,2)) ||
	      (matches(0,1) && matches(1,3)) ||
	      (matches(0,2) && matches(1,0)) ||
	      (matches(0,2) && matches(1,1)) ||
	      (matches(0,2) && matches(1,3)) ||
	      (matches(0,3) && matches(1,0)) ||
	      (matches(0,3) && matches(1,1)) ||
	      (matches(0,3) && matches(1,2))) {
	    nrl.push_back(lib[fi]);
	  }
	}
      }
    } else if (nspec==3) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap) {
	  if (matches(2,3)) {
	    if ((matches(0,0) && matches(0,1)) ||
		(matches(0,0) && matches(0,2)) ||
		(matches(0,1) && matches(0,0)) ||
		(matches(0,1) && matches(0,2)) ||
		(matches(0,2) && matches(0,0)) ||
		(matches(0,2) && matches(0,1))) {
	      nrl.push_back(lib[fi]);
	    }
	  } else if (matches(2,2)) {
	    if ((matches(0,0) && matches(0,1)) ||
		(matches(0,0) && matches(0,3)) ||
		(matches(0,1) && matches(0,0)) ||
		(matches(0,1) && matches(0,3)) ||
		(matches(0,3) && matches(0,0)) ||
		(matches(0,3) && matches(0,1))) {
	      nrl.push_back(lib[fi]);
	    }
	  } else if (matches(2,1)) {
	    if ((matches(0,0) && matches(0,3)) ||
		(matches(0,0) && matches(0,2)) ||
		(matches(0,3) && matches(0,0)) ||
		(matches(0,3) && matches(0,2)) ||
		(matches(0,2) && matches(0,0)) ||
		(matches(0,2) && matches(0,3))) {
	      nrl.push_back(lib[fi]);
	    }
	  } else if (matches(2,0)) {
	    if ((matches(0,3) && matches(0,1)) ||
		(matches(0,3) && matches(0,2)) ||
		(matches(0,1) && matches(0,3)) ||
		(matches(0,1) && matches(0,2)) ||
		(matches(0,2) && matches(0,3)) ||
		(matches(0,2) && matches(0,1))) {
	      nrl.push_back(lib[fi]);
	    }
	  }
	}
      }
    } else if (nspec==4) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap) {
	  if (matches(0,0)) {
	    if (matches(1,1)) {
	      if ((matches(2,2) && matches(3,3)) ||
		  (matches(2,3) && matches(3,2))) {
		nrl.push_back(lib[fi]);
	      }
	    } else if (matches(1,2)) {
	      if ((matches(2,1) && matches(3,3)) ||
		  (matches(2,3) && matches(3,1))) {
		nrl.push_back(lib[fi]);
	      }
	    } else if (matches(1,3)) {
	      if ((matches(2,1) && matches(3,2)) ||
		  (matches(2,2) && matches(3,1))) {
		nrl.push_back(lib[fi]);
	      }
	    }
	  } else if (matches(0,1)) {
	    if (matches(1,0)) {
	      if ((matches(2,2) && matches(3,3)) ||
		  (matches(2,3) && matches(3,2))) {
		nrl.push_back(lib[fi]);
	      }
	    } else if (matches(1,2)) {
	      if ((matches(2,0) && matches(3,3)) ||
		  (matches(2,3) && matches(3,0))) {
		nrl.push_back(lib[fi]);
	      }
	    } else if (matches(1,3)) {
	      if ((matches(2,0) && matches(3,2)) ||
		  (matches(2,2) && matches(3,0))) {
		nrl.push_back(lib[fi]);
	      }
	    }
	  } else if (matches(0,2)) {
	    if (matches(1,0)) {
	      if ((matches(2,1) && matches(3,3)) ||
		  (matches(2,3) && matches(3,1))) {
		nrl.push_back(lib[fi]);
	      }
	    } else if (matches(1,1)) {
	      if ((matches(2,0) && matches(3,3)) ||
		  (matches(2,3) && matches(3,0))) {
		nrl.push_back(lib[fi]);
	      }
	    } else if (matches(1,3)) {
	      if ((matches(2,0) && matches(3,1)) ||
		  (matches(2,1) && matches(3,0))) {
		nrl.push_back(lib[fi]);
	      }
	    }
	  } else if (matches(0,3)) {
	    if (matches(1,0)) {
	      if ((matches(2,1) && matches(3,2)) ||
		  (matches(2,2) && matches(3,1))) {
		nrl.push_back(lib[fi]);
	      }
	    } else if (matches(1,1)) {
	      if ((matches(2,0) && matches(3,2)) ||
		  (matches(2,2) && matches(3,0))) {
		nrl.push_back(lib[fi]);
	      }
	    } else if (matches(1,2)) {
	      if ((matches(2,0) && matches(3,1)) ||
		  (matches(2,1) && matches(3,0))) {
		nrl.push_back(lib[fi]);
	      }
	    }
	  }
	}
      }
    } else if (nspec==5) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && (matches(4,4) || matches(4,5))) {
	  if (matches(0,0)) {
	    if (matches(1,1)) {
	      if ((matches(2,2) && matches(3,3)) ||
		  (matches(2,3) && matches(3,2))) {
		nrl.push_back(lib[fi]);
	      }
	    } else if (matches(1,2)) {
	      if ((matches(2,1) && matches(3,3)) ||
		  (matches(2,3) && matches(3,1))) {
		nrl.push_back(lib[fi]);
	      }
	    } else if (matches(1,3)) {
	      if ((matches(2,1) && matches(3,2)) ||
		  (matches(2,2) && matches(3,1))) {
		nrl.push_back(lib[fi]);
	      }
	    }
	  } else if (matches(0,1)) {
	    if (matches(1,0)) {
	      if ((matches(2,2) && matches(3,3)) ||
		  (matches(2,3) && matches(3,2))) {
		nrl.push_back(lib[fi]);
	      }
	    } else if (matches(1,2)) {
	      if ((matches(2,0) && matches(3,3)) ||
		  (matches(2,3) && matches(3,0))) {
		nrl.push_back(lib[fi]);
	      }
	    } else if (matches(1,3)) {
	      if ((matches(2,0) && matches(3,2)) ||
		  (matches(2,2) && matches(3,0))) {
		nrl.push_back(lib[fi]);
	      }
	    }
	  } else if (matches(0,2)) {
	    if (matches(1,0)) {
	      if ((matches(2,1) && matches(3,3)) ||
		  (matches(2,3) && matches(3,1))) {
		nrl.push_back(lib[fi]);
	      }
	    } else if (matches(1,1)) {
	      if ((matches(2,0) && matches(3,3)) ||
		  (matches(2,3) && matches(3,0))) {
		nrl.push_back(lib[fi]);
	      }
	    } else if (matches(1,3)) {
	      if ((matches(2,0) && matches(3,1)) ||
		  (matches(2,1) && matches(3,0))) {
		nrl.push_back(lib[fi]);
	      }
	    }
	  } else if (matches(0,3)) {
	    if (matches(1,0)) {
	      if ((matches(2,1) && matches(3,2)) ||
		  (matches(2,2) && matches(3,1))) {
		nrl.push_back(lib[fi]);
	      }
	    } else if (matches(1,1)) {
	      if ((matches(2,0) && matches(3,2)) ||
		  (matches(2,2) && matches(3,0))) {
		nrl.push_back(lib[fi]);
	      }
	    } else if (matches(1,2)) {
	      if ((matches(2,0) && matches(3,1)) ||
		  (matches(2,1) && matches(3,0))) {
		nrl.push_back(lib[fi]);
	      }
	    }
	  }
	}
      }
    } else if (nspec==6) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap) {
	  if ((matches(4,4) && matches(5,5)) ||
	      (matches(4,5) && matches(5,4))) {
	    if (matches(0,0)) {
	      if (matches(1,1)) {
		if ((matches(2,2) && matches(3,3)) ||
		    (matches(2,3) && matches(3,2))) {
		  nrl.push_back(lib[fi]);
		}
	      } else if (matches(1,2)) {
		if ((matches(2,1) && matches(3,3)) ||
		    (matches(2,3) && matches(3,1))) {
		  nrl.push_back(lib[fi]);
		}
	      } else if (matches(1,3)) {
		if ((matches(2,1) && matches(3,2)) ||
		    (matches(2,2) && matches(3,1))) {
		  nrl.push_back(lib[fi]);
		}
	      }
	    } else if (matches(0,1)) {
	      if (matches(1,0)) {
		if ((matches(2,2) && matches(3,3)) ||
		    (matches(2,3) && matches(3,2))) {
		  nrl.push_back(lib[fi]);
		}
	      } else if (matches(1,2)) {
		if ((matches(2,0) && matches(3,3)) ||
		    (matches(2,3) && matches(3,0))) {
		  nrl.push_back(lib[fi]);
		}
	      } else if (matches(1,3)) {
		if ((matches(2,0) && matches(3,2)) ||
		    (matches(2,2) && matches(3,0))) {
		  nrl.push_back(lib[fi]);
		}
	      }
	    } else if (matches(0,2)) {
	      if (matches(1,0)) {
		if ((matches(2,1) && matches(3,3)) ||
		    (matches(2,3) && matches(3,1))) {
		  nrl.push_back(lib[fi]);
		}
	      } else if (matches(1,1)) {
		if ((matches(2,0) && matches(3,3)) ||
		    (matches(2,3) && matches(3,0))) {
		  nrl.push_back(lib[fi]);
		}
	      } else if (matches(1,3)) {
		if ((matches(2,0) && matches(3,1)) ||
		    (matches(2,1) && matches(3,0))) {
		  nrl.push_back(lib[fi]);
		}
	      }
	    } else if (matches(0,3)) {
	      if (matches(1,0)) {
		if ((matches(2,1) && matches(3,2)) ||
		    (matches(2,2) && matches(3,1))) {
		  nrl.push_back(lib[fi]);
		}
	      } else if (matches(1,1)) {
		if ((matches(2,0) && matches(3,2)) ||
		    (matches(2,2) && matches(3,0))) {
		  nrl.push_back(lib[fi]);
		}
	      } else if (matches(1,2)) {
		if ((matches(2,0) && matches(3,1)) ||
		    (matches(2,1) && matches(3,0))) {
		  nrl.push_back(lib[fi]);
		}
	      }
	    }
	  }
	}
      }
    } else {
      O2SCL_ERR2("Too many nuclei specified for chapter in ",
		     "reaction_lib::find_in_chap().",exc_efailed);
    }

  } else if (chap==11) {
    // nuc1 -> nuc2 + nuc3 + nuc4 + nuc5
    
    if (nspec==1) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && matches(0,0)) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else if (nspec==2) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && 
	    ((matches(0,0) && matches(1,1)) ||
	     (matches(0,0) && matches(1,2)) ||
	     (matches(0,0) && matches(1,3)) ||
	     (matches(0,0) && matches(1,4)))) {
	  nrl.push_back(lib[fi]);
	}
      }
    } else if (nspec==3) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && matches(0,0)) {
	  if (matches(1,1)) {
	    if (matches(2,2)) {
	      nrl.push_back(lib[fi]);
	    } else if (matches(2,3)) {
	      nrl.push_back(lib[fi]);
	    } else if (matches(2,4)) {
	      nrl.push_back(lib[fi]);
	    }
	  } else if (matches(1,2)) {
	    if (matches(2,1)) {
	      nrl.push_back(lib[fi]);
	    } else if (matches(2,3)) {
	      nrl.push_back(lib[fi]);
	    } else if (matches(2,4)) {
	      nrl.push_back(lib[fi]);
	    }
	  } else if (matches(1,3)) {
	    if (matches(2,1)) {
	      nrl.push_back(lib[fi]);
	    } else if (matches(2,2)) {
	      nrl.push_back(lib[fi]);
	    } else if (matches(2,4)) {
	      nrl.push_back(lib[fi]);
	    }
	  } else if (matches(1,4)) {
	    if (matches(2,1)) {
	      nrl.push_back(lib[fi]);
	    } else if (matches(2,2)) {
	      nrl.push_back(lib[fi]);
	    } else if (matches(2,3)) {
	      nrl.push_back(lib[fi]);
	    }
	  }
	}
      }
    } else if (nspec==4) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && matches(0,0)) {
	  if (matches(1,1)) {
	    if ((matches(2,2) && matches(3,3)) ||
		(matches(2,2) && matches(3,4)) ||
		(matches(2,3) && matches(3,2)) ||
		(matches(2,3) && matches(3,4)) ||
		(matches(2,4) && matches(3,2)) ||
		(matches(2,4) && matches(3,3))) {
	      nrl.push_back(lib[fi]);
	    }
	  } else if (matches(1,2)) {
	    if ((matches(2,1) && matches(3,3)) ||
		(matches(2,1) && matches(3,4)) ||
		(matches(2,3) && matches(3,1)) ||
		(matches(2,3) && matches(3,4)) ||
		(matches(2,4) && matches(3,1)) ||
		(matches(2,4) && matches(3,3))) {
	      nrl.push_back(lib[fi]);
	    }
	  } else if (matches(1,3)) {
	    if ((matches(2,2) && matches(3,1)) ||
		(matches(2,2) && matches(3,4)) ||
		(matches(2,1) && matches(3,2)) ||
		(matches(2,1) && matches(3,4)) ||
		(matches(2,4) && matches(3,2)) ||
		(matches(2,4) && matches(3,1))) {
	      nrl.push_back(lib[fi]);
	    }
	  } else if (matches(1,4)) {
	    if ((matches(2,2) && matches(3,3)) ||
		(matches(2,2) && matches(3,1)) ||
		(matches(2,3) && matches(3,2)) ||
		(matches(2,3) && matches(3,1)) ||
		(matches(2,1) && matches(3,2)) ||
		(matches(2,1) && matches(3,3))) {
	      nrl.push_back(lib[fi]);
	    }
	  }
	}
      }
    } else if (nspec==5) {
      for(fi=0;fi<lib.size();fi++) {
	if (lib[fi].chap==chap && matches(0,0)) {
	  if (matches(1,1)) {
	    if ((matches(2,2) && matches(3,3) && matches(4,4)) ||
		(matches(2,2) && matches(3,4) && matches(4,3)) ||
		(matches(2,3) && matches(3,2) && matches(4,4)) ||
		(matches(2,3) && matches(3,4) && matches(4,2)) ||
		(matches(2,4) && matches(3,2) && matches(4,3)) ||
		(matches(2,4) && matches(3,3) && matches(4,2))) {
	      nrl.push_back(lib[fi]);
	    }
	  } else if (matches(1,2)) {
	    if ((matches(2,1) && matches(3,3) && matches(4,4)) ||
		(matches(2,1) && matches(3,4) && matches(4,3)) ||
		(matches(2,3) && matches(3,1) && matches(4,4)) ||
		(matches(2,3) && matches(3,4) && matches(4,1)) ||
		(matches(2,4) && matches(3,1) && matches(4,3)) ||
		(matches(2,4) && matches(3,3) && matches(4,1))) {
	      nrl.push_back(lib[fi]);
	    }
	  } else if (matches(1,3)) {
	    if ((matches(2,2) && matches(3,1) && matches(4,4)) ||
		(matches(2,2) && matches(3,4) && matches(4,1)) ||
		(matches(2,1) && matches(3,2) && matches(4,4)) ||
		(matches(2,1) && matches(3,4) && matches(4,2)) ||
		(matches(2,4) && matches(3,2) && matches(4,1)) ||
		(matches(2,4) && matches(3,1) && matches(4,2))) {
	      nrl.push_back(lib[fi]);
	    }
	  } else if (matches(1,4)) {
	    if ((matches(2,2) && matches(3,3) && matches(4,1)) ||
		(matches(2,2) && matches(3,1) && matches(4,3)) ||
		(matches(2,3) && matches(3,2) && matches(4,1)) ||
		(matches(2,3) && matches(3,1) && matches(4,2)) ||
		(matches(2,1) && matches(3,2) && matches(4,3)) ||
		(matches(2,1) && matches(3,3) && matches(4,2))) {
	      nrl.push_back(lib[fi]);
	    }
	  }
	}
      }
    } else {
      O2SCL_ERR2("Too many nuclei specified for chapter in ",
		     "reaction_lib::find_in_chap().",exc_efailed);
    }
    
  } else {
    O2SCL_ERR("Invalid chapter in reaction_lib::find_in_chap().",
		  exc_efailed);
  }

  return 0;
}

/// Read from a file
int reaction_lib::read_file_reaclib2(std::string fname) {

  nuclear_reaction r;
  ifstream ins(fname.c_str());
  // for parse_elstring()
  nucmass_info nm;
  
  bool eof=false;

  do {

    size_t chap;
    string s0,s1,s2,s3;

    if (!(ins >> chap)) {
      eof=true;
    }
    if (!eof) {
      getline(ins,s0);

      getline(ins,s1);
      getline(ins,s2);
      getline(ins,s3);
  
      r.chap=chap;
      if (s1.substr(0,5)!="     ") {
	O2SCL_ERR("Failed first five spaces.",exc_efailed);
      }
      for(size_t i=0;i<6;i++) {
	r.name[i]=s1.substr(5+i*5,5);

	if (r.name[i]=="     ") {
	  r.Z[i]=0;
	  r.A[i]=0;
	  r.isomer[i]=0;
	} else {
	  // Remove whitespace from string
	  istringstream iss(r.name[i]);
	  iss >> r.name[i];
	  // Parse the aluminum isomer separately
	  if (r.name[i]=="Al*6") {
	    r.Z[i]=13;
	    r.A[i]=26;
	    r.isomer[i]=1;
	  } else {
	    // Normal nuclei
	    int N, Z, A;
	    int ret=nm.parse_elstring(r.name[i],Z,N,A);
	    if (ret!=0) {
	      O2SCL_ERR("Failed to parse.",exc_efailed);
	    }
	    r.Z[i]=(size_t)Z;
	    r.A[i]=(size_t)A;
	    r.isomer[i]=0;
	  }
	}
      }
      if (s1.substr(35,8)!="        ") {
	O2SCL_ERR("Failed eight spaces.",exc_efailed);
      }
      r.ref=s1.substr(43,4);
      r.type=s1[47];
      r.rev=s1[48];
      if (s1.substr(49,3)!="   ") {
	O2SCL_ERR("Failed three spaces.",exc_efailed);
      }
      if (s1.substr(52,12)[8]!='e') {
	O2SCL_ERR("Number failed.",exc_efailed);
      }
      r.Q=o2scl::stod(s1.substr(52,12));
  
      for(size_t i=0;i<4;i++) {
	if (s2.substr(i*13,13)[9]!='e') {
	  O2SCL_ERR("Number failed.",exc_efailed);
	}
	r.a[i]=o2scl::stod(s2.substr(i*13,13));
      }

      for(size_t i=0;i<3;i++) {
	if (s3.substr(i*13,13)[9]!='e') {
	  O2SCL_ERR("Number failed.",exc_efailed);
	}
	r.a[i+4]=o2scl::stod(s3.substr(i*13,13));
      }
      
      lib.push_back(r);
    }

  } while (eof==false);
    
  return 0;
}

