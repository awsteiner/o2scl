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
#include <o2scl/find_constants.h>
#include <o2scl/lib_settings.h>
#include <o2scl/convert_units.h>

#include <fnmatch.h>

#include <boost/algorithm/string.hpp>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

find_constants::find_constants() {
    
  list={{{"Schwarzchild radius"},
	 "m",o2scl_const::o2scl_mks,o2scl_mks::schwarzchild_radius},
	{{"Schwarzchild radius"},
	 "cm",o2scl_const::o2scl_cgs,o2scl_cgs::schwarzchild_radius},
	{{"Schwarzchild radius"},
	 "km",o2scl_const::o2scl_mks,o2scl_mks::schwarzchild_radius/1.0e3},
	{{"speed of light","c"},
	 "m/s",o2scl_const::o2scl_mks,
	 o2scl_const::speed_of_light_f<double>(o2scl_const::o2scl_mks)},
	{{"speed of light","c"},
	 "cm/s",o2scl_const::o2scl_cgs,
	 o2scl_const::speed_of_light_f<double>(o2scl_const::o2scl_cgs)},
	{{"Gravitational constant","G","Newton's constant"},
	 "m^3/kg/s^2",o2scl_const::o2scl_mks,
	 o2scl_mks::gravitational_constant},
	{{"Gravitational constant","G","Newton's constant"},
	 "cm^3/g/s^2",o2scl_const::o2scl_cgs,
	 o2scl_cgs::gravitational_constant},
	{{"Boltzmann constant","kB","k_B"},
	 "m^2/kg/s^2/K",o2scl_const::o2scl_mks,
	 o2scl_mks::boltzmann},
	{{"Boltzmann constant","kB","k_B"},
	 "cm^2/g/s^2/K",o2scl_const::o2scl_cgs,
	 o2scl_cgs::boltzmann},
	{{"Stefan-Boltzmann constant","sigma_SB"},
	 "kg/s^3/K^4",o2scl_const::o2scl_mks,
	 o2scl_mks::stefan_boltzmann_constant},
	{{"Stefan-Boltzmann constant","sigma_SB"},
	 "g/s^3/K^4",o2scl_const::o2scl_cgs,
	 o2scl_cgs::stefan_boltzmann_constant},
	{{"Planck's constant","h"},
	 "kg*m^2/s",o2scl_const::o2scl_mks,
	 o2scl_const::planck_f<double>(o2scl_const::o2scl_mks)},
	{{"Planck's constant","h"},
	 "g*cm^2/s",o2scl_const::o2scl_cgs,
	 o2scl_const::planck_f<double>(o2scl_const::o2scl_cgs)},
	{{"reduced Planck's constant","hbar"},
	 "kg*m^2/s",o2scl_const::o2scl_mks,
	 o2scl_const::hbar_f<double>(o2scl_const::o2scl_mks)},
	{{"reduced Planck's constant","hbar"},
	 "g*cm^2/s",o2scl_const::o2scl_cgs,
	 o2scl_const::hbar_f<double>(o2scl_const::o2scl_cgs)},
	{{"Avogadro's number","N_A","NA"},
	 "",0,o2scl_const::avogadro},
	{{"pi"},"",0,o2scl_const::pi},
	{{"pi^2","pi squared","pi2"},"",0,o2scl_const::pi2},
	{{"sin^2(theta_W)"},"",0,o2scl_const::sin2_theta_weak},
	{{"GFermi","G_Fermi"},"s^4/m^4/kg^2",o2scl_const::o2scl_mks,
	 o2scl_mks::Gfermi},
	{{"GFermi","G_Fermi"},"s^4/cm^4/g^2",o2scl_const::o2scl_cgs,
	 o2scl_cgs::Gfermi},
	{{"elementary charge","electron charge","e"},"C",
	 o2scl_const::o2scl_mks,o2scl_const::elem_charge_f<double>()},
	{{"hbar c","hbar*c","hbarc"},"MeV*fm",
	 0,o2scl_const::hc_mev_fm_f<double>()},
	{{"mass electron","electron mass"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_electron},
	{{"mass muon","muon mass"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_muon},
	{{"mass tau","tau mass"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_tau},
	{{"mass neutron","neutron mass"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_neutron},
	{{"mass proton","proton mass"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_proton},
	{{"mass deuteron","deuteron mass"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_deuteron},
	{{"mass triton","triton mass"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_triton},
	{{"mass helion","helion mass"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_helion},
	{{"mass alpha","alpha mass"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_alpha},
	{{"mass electron","electron mass"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_electron},
	{{"mass muon","muon mass"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_muon},
	{{"mass tau","tau mass"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_tau},
	{{"mass neutron","neutron mass"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_neutron},
	{{"mass proton","proton mass"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_proton},
	{{"mass deuteron","deuteron mass"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_deuteron},
	{{"mass triton","triton mass"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_triton},
	{{"mass helion","helion mass"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_helion},
	{{"mass alpha","alpha mass"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_alpha},
	{{"mass lambda","lambda mass"},"MeV",0,
	 o2scl_const::mass_lambda_MeV},
	{{"mass sigma minus","sigma minus mass"},"MeV",0,
	 o2scl_const::mass_sigma_minus_MeV},
	{{"mass sigma zero","sigma zero mass"},"MeV",0,
	 o2scl_const::mass_sigma_zero_MeV},
	{{"mass sigma plus","sigma plus mass"},"MeV",0,
	 o2scl_const::mass_sigma_plus_MeV},
	{{"mass cascade zero","cascade zero mass"},"MeV",0,
	 o2scl_const::mass_cascade_zero_MeV},
	{{"mass cascade minus","cascade minus mass"},"MeV",0,
	 o2scl_const::mass_cascade_minus_MeV},
	{{"mass up","up mass"},"MeV",0,
	 o2scl_const::mass_up_MeV},
	{{"mass down","down mass"},"MeV",0,
	 o2scl_const::mass_down_MeV},
	{{"mass strange","strange mass"},"MeV",0,
	 o2scl_const::mass_strange_MeV}
  };

}

int find_constants::find_nothrow(std::string name, std::string unit,
				 vector<find_constants_list> &matches) {
  
  o2scl::convert_units<> &cu=o2scl_settings.get_convert_units();
  
  // Simplify by setting dashes and hyphens to spaces
  for(size_t i=0;i<name.length();i++) {
    if (name[i]=='-' || name[i]=='_') name[i]=' ';
  }
    
  // Start with a fresh list
  matches.clear();
    
  // Temporarily store matching indexes
  vector<size_t> indexes;

  int match_type=0, match_exact=1, match_pattern=2;
    
  // Initial pass, exact matches
  for(size_t i=0;i<list.size();i++) {
    for(size_t j=0;j<list[i].names.size();j++) {
      if (boost::iequals(name,list[i].names[j])) {
	indexes.push_back(i);
	// Now that we've found a match, don't look in the
	// other names for this list entry
	j=list[i].names.size();
	match_type=match_exact;
      }
    }
  }

  // No matches, so try wildcard matches
  if (indexes.size()==0) {
      
    string name_wc=((string)"*")+name+"*";
    for(size_t i=0;i<list.size();i++) {
      for(size_t j=0;j<list[i].names.size();j++) {
	if (fnmatch(name_wc.c_str(),(list[i].names[j]).c_str(),
		    FNM_CASEFOLD)==0) {
	  indexes.push_back(i);
	  // Now that we've found a match, don't look in the
	  // other names for this list entry
	  j=list[i].names.size();
	  match_type=match_pattern;
	}
      }
    }
      
  }

  // There was only one match
  if (indexes.size()==1) {

    // Add to 'matches' list
    matches.push_back(list[indexes[0]]);
      
    // Unit unspecified or matching
    if (unit.length()==0 ||
	(unit=="mks" &&
	 list[indexes[0]].unit_flag==o2scl_const::o2scl_mks) ||
	(unit=="cgs" &&
	 list[indexes[0]].unit_flag==o2scl_const::o2scl_cgs) ||
	boost::iequals(unit,list[indexes[0]].unit)) {
      if (match_type==match_exact) {
	return one_exact_match_unit_match;
      } else {
	return one_pattern_match_unit_match;
      }
    }
      
    // Try to convert units
    if (unit.length()>0) {
      double val2;
      cout << "Trying to convert from "
	   << list[indexes[0]].unit << " to "
	   << unit << endl;
      int cret=cu.convert_ret(list[indexes[0]].unit,unit,
			      list[indexes[0]].val,val2);
      if (cret==0) {
	// Update the value with the unit conversion and
	// the unit with the new unit
	matches[0].val=val2;
	matches[0].unit=unit;
	if (match_type==match_exact) {
	  return one_exact_match_unit_match;
	} else {
	  return one_pattern_match_unit_match;
	}
      }
    }

    if (match_type==match_exact) {
      return one_exact_match_unit_diff;
    } else {
      return one_pattern_match_unit_diff;
    }
  }

  if (indexes.size()>0 && unit=="") {
    // No unit string was given, so just return
    for(size_t i=0;i<indexes.size();i++) {
      matches.push_back(list[indexes[i]]);
    }
    if (match_type==match_exact) {
      return exact_matches_no_unit;
    } else {
      return pattern_matches_no_unit;
    }
  }
    
  if (indexes.size()>0) {

    // We found at least one match, check unit
      
    vector<size_t> indexes2;
      
    // Look for entries with matching unit
    for(size_t i=0;i<indexes.size();i++) {
      if ((unit=="cgs" &&
	   list[indexes[i]].unit_flag==o2scl_const::o2scl_cgs) ||
	  (unit=="mks" &&
	   list[indexes[i]].unit_flag==o2scl_const::o2scl_mks) ||
	  boost::iequals(list[indexes[i]].unit,unit)) {
	indexes2.push_back(indexes[i]);
      }
    }
      
    if (indexes2.size()==0) {
      // No matching unit, try to convert
      for(size_t i=0;i<indexes.size();i++) {
	double val2;
	cout << "Trying to convert from "
	     << list[indexes[i]].unit << " to "
	     << unit << endl;
	int cret=cu.convert_ret(list[indexes[i]].unit,unit,
				list[indexes[i]].val,val2);
	if (cret==0) {
	  matches.push_back(list[indexes[i]]);
	  // Update the value with the unit conversion and
	  // the unit with the new unit
	  matches[0].val=val2;
	  matches[0].unit=unit;
	  if (match_type==match_exact) {
	    return one_exact_match_unit_match;
	  } else {
	    return one_pattern_match_unit_match;
	  }
	}
      }

      // If no matching unit, just return the list of name matches
      for(size_t i=0;i<indexes.size();i++) {
	matches.push_back(list[indexes[i]]);
      }
      if (match_type==match_exact) {
	return exact_matches_unit_diff;
      } else {
	return pattern_matches_unit_diff;
      }
    }

    // Return only the entries with matching units
    for(size_t i=0;i<indexes2.size();i++) {
      matches.push_back(list[indexes2[i]]);
    }
    if (match_type==match_exact) {
      return exact_matches_unit_match;
    } else {
      return pattern_matches_unit_match;
    }
  }

  return no_matches;
}

void find_constants::find_print(std::string name, std::string unit,
				size_t prec) {

  cout.precision(prec);
    
  vector<find_constants_list> matches;
  int ret=find_nothrow(name,unit,matches);
  if (ret==no_matches) {
    cout << "No matches found for name " << name << endl;
    return;
  }

  if (unit.length()==0) {
    cout << "Matches for " << name;
  } else {
    cout << "Matches for " << name << " in " << unit;
  }
  if (ret==one_exact_match_unit_diff ||
      ret==exact_matches_unit_diff) {
    cout << " (no matching units)" << endl;
  }
  cout << ": " << endl;
  for(size_t i=0;i<matches.size();i++) {
    cout << "(" << i+1 << "/" << matches.size() << ") "
	 << matches[i].names[0] << " = "
	 << matches[i].val << " "
	 << matches[i].unit << endl;
  }
  return;
}
  
double find_constants::find_unique(std::string name, std::string unit) {
  vector<find_constants_list> matches;
  int ret=find_nothrow(name,unit,matches);
  if (ret!=one_exact_match_unit_match &&
      ret!=one_pattern_match_unit_match) {
    O2SCL_ERR2("Failed to find unique match in ",
	       "find_constants::find_unique().",o2scl::exc_einval);
  }
  return matches[0].val;
}
  
