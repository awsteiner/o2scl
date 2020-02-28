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
    
  list={{{"schwarzchildradius","rschwarz"},
	 "m",o2scl_const::o2scl_mks,o2scl_mks::schwarzchild_radius},
	{{"schwarzchildradius","rschwarz"},
	 "cm",o2scl_const::o2scl_cgs,o2scl_cgs::schwarzchild_radius},
	{{"schwarzchildradius","rschwarz"},
	 "km",o2scl_const::o2scl_mks,o2scl_mks::schwarzchild_radius/1.0e3},
	{{"speedoflight","c","lightspeed"},
	 "m/s",o2scl_const::o2scl_mks,
	 o2scl_const::speed_of_light_f<double>(o2scl_const::o2scl_mks)},
	{{"speedoflight","c","lightspeed"},
	 "cm/s",o2scl_const::o2scl_cgs,
	 o2scl_const::speed_of_light_f<double>(o2scl_const::o2scl_cgs)},
	{{"gravitationalconstant","g","newtonsconstant",
	  "newtonconstant"},
	 "m^3/kg/s^2",o2scl_const::o2scl_mks,
	 o2scl_mks::gravitational_constant},
	{{"gravitationalconstant","g","newtonsconstant",
	  "newtonconstant"},
	 "cm^3/g/s^2",o2scl_const::o2scl_cgs,
	 o2scl_cgs::gravitational_constant},
	{{"boltzmannconstant","kb",
	  "boltzmannsconstant"},
	 "m^2/kg/s^2/K",o2scl_const::o2scl_mks,
	 o2scl_mks::boltzmann},
	{{"boltzmannconstant","kb",
	  "boltzmannsconstant"},
	 "cm^2/g/s^2/K",o2scl_const::o2scl_cgs,
	 o2scl_cgs::boltzmann},
	{{"stefanboltzmannconstant","sigmasb","stefanboltzmann","ssb"},
	 "kg/s^3/K^4",o2scl_const::o2scl_mks,
	 o2scl_mks::stefan_boltzmann_constant},
	{{"stefanboltzmannconstant","sigmasb","stefanboltzmann","ssb"},
	 "g/s^3/K^4",o2scl_const::o2scl_cgs,
	 o2scl_cgs::stefan_boltzmann_constant},
	{{"plancksconstant","h","planckconstant"},
	 "kg*m^2/s",o2scl_const::o2scl_mks,
	 o2scl_const::planck_f<double>(o2scl_const::o2scl_mks)},
	{{"plancksconstant","h","planckconstant"},
	 "g*cm^2/s",o2scl_const::o2scl_cgs,
	 o2scl_const::planck_f<double>(o2scl_const::o2scl_cgs)},
	{{"reducedplancksconstant","hbar"},
	 "kg*m^2/s",o2scl_const::o2scl_mks,
	 o2scl_const::hbar_f<double>(o2scl_const::o2scl_mks)},
	{{"reducedplancksconstant","hbar"},
	 "g*cm^2/s",o2scl_const::o2scl_cgs,
	 o2scl_const::hbar_f<double>(o2scl_const::o2scl_cgs)},
	{{"avogadrosnumber","NA","avogadro"},
	 "",0,o2scl_const::avogadro},
	{{"alphaem","finestructure","alpha"},"",0,
	 o2scl_const::fine_structure},
	{{"pi"},"",0,o2scl_const::pi},
	{{"zeta32"},"",0,o2scl_const::zeta32},
	{{"zeta2"},"",0,o2scl_const::zeta2},
	{{"zeta52"},"",0,o2scl_const::zeta52},
	{{"zeta3"},"",0,o2scl_const::zeta3},
	{{"zeta5"},"",0,o2scl_const::zeta5},
	{{"zeta7"},"",0,o2scl_const::zeta7},
	{{"pi2","pisquared"},"",0,o2scl_const::pi2},
	{{"sin2thetaW"},"",0,o2scl_const::sin2_theta_weak},
	{{"GFermi"},"s^4/m^4/kg^2",o2scl_const::o2scl_mks,
	 o2scl_mks::Gfermi},
	{{"GFermi"},"s^4/cm^4/g^2",o2scl_const::o2scl_cgs,
	 o2scl_cgs::Gfermi},
	{{"elementarycharge","electroncharge","e","chargeelectron",
	  "qelectron"},"C",
	 o2scl_const::o2scl_mks,o2scl_const::elem_charge_f<double>()},
	{{"hbarc"},"MeV*fm",
	 0,o2scl_const::hc_mev_fm_f<double>()},
	{{"masselectron","electronmass"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_electron},
	{{"massmuon","muonmass"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_muon},
	{{"masstau","taumass"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_tau},
	{{"massneutron","neutronmass"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_neutron},
	{{"massproton","protonmass"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_proton},
	{{"massdeuteron","deuteronmass"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_deuteron},
	{{"masstriton","tritonmass"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_triton},
	{{"masshelion","helionmass"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_helion},
	{{"massalpha","alphamass"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_alpha},
	{{"masselectron","electronmass"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_electron},
	{{"massmuon","muonmass"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_muon},
	{{"masstau","taumass"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_tau},
	{{"massneutron","neutronmass"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_neutron},
	{{"massproton","protonmass"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_proton},
	{{"massdeuteron","deuteronmass"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_deuteron},
	{{"masstriton","tritonmass"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_triton},
	{{"masshelion","helionmass"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_helion},
	{{"massalpha","alphamass"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_alpha},
	{{"masslambda","lambdamass"},"MeV",0,
	 o2scl_const::mass_lambda_MeV},
	{{"masssigma minus","sigma minusmass"},"MeV",0,
	 o2scl_const::mass_sigma_minus_MeV},
	{{"masssigma zero","sigma zeromass"},"MeV",0,
	 o2scl_const::mass_sigma_zero_MeV},
	{{"masssigma plus","sigma plusmass"},"MeV",0,
	 o2scl_const::mass_sigma_plus_MeV},
	{{"masscascade zero","cascade zeromass"},"MeV",0,
	 o2scl_const::mass_cascade_zero_MeV},
	{{"masscascade minus","cascade minusmass"},"MeV",0,
	 o2scl_const::mass_cascade_minus_MeV},
	{{"massup","upmass"},"MeV",0,
	 o2scl_const::mass_up_MeV},
	{{"massdown","downmass"},"MeV",0,
	 o2scl_const::mass_down_MeV},
	{{"massstrange","strangemass"},"MeV",0,
	 o2scl_const::mass_strange_MeV},
	{{"masssolar","solarmass","masssun","sunmass","msun"},
	 "kg",o2scl_const::o2scl_mks,
	 o2scl_mks::solar_mass},
	{{"massmercury","mercurymass","mmercury"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mercury_mass},
	{{"massvenus","venusmass","mvenus"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::venus_mass},
	{{"massearth","earthmass","mearth"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::earth_mass},
	{{"massmars","marsmass","mmars"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mars_mass},
	{{"massjupiter","jupitermass","mjupiter","mjup"},
	 "kg",o2scl_const::o2scl_mks,
	 o2scl_mks::jupiter_mass},
	{{"masssaturn","saturnmass","msaturn"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::saturn_mass},
	{{"massuranus","uranusmass"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::uranus_mass},
	{{"massneptune","neptunemass","mneptune"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::neptune_mass},
	{{"masspluto","plutomass","mpluto"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::pluto_mass},
	{{"radiussolar","solarradius","radiussun","sunradius","rsun"},
	 "m",o2scl_const::o2scl_mks,o2scl_mks::solar_radius},
	{{"radiusmercury","mercuryradius","rmercury"},"m",
	 o2scl_const::o2scl_mks,
	 o2scl_mks::mercury_radius},
	{{"radiusvenus","venusradius","rvenus"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::venus_radius},
	{{"radiusearthequatorial","earthequatorialradius",
	  "earthradiusequatorial"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::earth_radius_equatorial},
	{{"radiusearthpolar","earthpolarradius",
	  "earthradiuspolar"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::earth_radius_polar},
	{{"radiusmarsequatorial","marsequatorialradius",
	  "marsradiusequatorial"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::mars_radius_equatorial},
	{{"radiusmarspolar","marspolarradius",
	  "marsradiuspolar"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::mars_radius_polar},
	{{"radiusjupiterequatorial","jupiterequatorialradius",
	  "jupiterradiusequatorial"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::jupiter_radius_equatorial},
	{{"radiusjupiterpolar","jupiterpolarradius",
	  "jupiterradiuspolar"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::jupiter_radius_polar},
	{{"radiussaturnequatorial","saturnequatorialradius",
	  "saturnradiusequatorial"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::saturn_radius_equatorial},
	{{"radiussaturnpolar","saturnpolarradius",
	  "saturnradiuspolar"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::saturn_radius_polar},
	{{"radiusuranusequatorial","uranusequatorialradius",
	  "uranusradiusequatorial"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::uranus_radius_equatorial},
	{{"radiusuranuspolar","uranuspolarradius",
	  "uranusradiuspolar"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::uranus_radius_polar},
	{{"radiusneptuneequatorial","neptuneequatorialradius",
	  "neptuneradiusequatorial"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::neptune_radius_equatorial},
	{{"radiusneptunepolar","neptunepolarradius",
	  "neptuneradiuspolar"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::neptune_radius_polar},
	{{"radiuspluto","plutoradius","rpluto"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::pluto_radius},
	{{"rydberg"},"g*cm^2/s^2",o2scl_const::o2scl_cgs,
	 o2scl_cgs::rydberg}
  };

}

int find_constants::find_nothrow(std::string name, std::string unit,
				 vector<find_constants_list> &matches) {
  
  o2scl::convert_units<> &cu=o2scl_settings.get_convert_units();
  
  // Simplify by setting dashes and hyphens to spaces
  cout << "Before: " << name << endl;
  for(size_t i=0;i<name.length();i++) {
    if (!isalnum(name[i])) {
      name.erase(i,1);
      i=0;
    }
  }
  cout << "After: " << name << endl;
    
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
      cu.verbose=1;
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
	cu.verbose=1;
	int cret=cu.convert_ret(list[indexes[i]].unit,unit,
				list[indexes[i]].val,val2);
	if (cret==0) {
	  matches.push_back(list[indexes[i]]);
	  // Update the value with the unit conversion and
	  // the unit with the new unit
	  matches[matches.size()-1].val=val2;
	  matches[matches.size()-1].unit=unit;
	}
      }

      if (matches.size()>0) {
	if (matches.size()==1) {
	  if (match_type==match_exact) {
	    return one_exact_match_unit_match;
	  } else {
	    return one_pattern_match_unit_match;
	  }
	} else {
	  if (match_type==match_exact) {
	    return exact_matches_unit_match;
	  } else {
	    return pattern_matches_unit_match;
	  }
	}
      }

      // If no matching unit conversions, just return the list of name
      // matches
      for(size_t i=0;i<indexes.size();i++) {
	matches.push_back(list[indexes[i]]);
      }
      if (match_type==match_exact) {
	return exact_matches_unit_diff;
      } else {
	return pattern_matches_unit_diff;
      }
    }

    /*
    // Return only the entries with matching units
    for(size_t i=0;i<indexes2.size();i++) {
    matches.push_back(list[indexes2[i]]);
    }
    if (match_type==match_exact) {
    return exact_matches_unit_match;
    } else {
    return pattern_matches_unit_match;
    }
    */
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
  
