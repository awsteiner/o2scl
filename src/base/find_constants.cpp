/*
  -------------------------------------------------------------------
  
  Copyright (C) 2020-2021, Andrew W. Steiner
  
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
#include <o2scl/vector.h>

#include <fnmatch.h>

#include <boost/algorithm/string.hpp>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

find_constants::find_constants() {

  /*
    The code is written in such a way that multiple entries for the
    same constant must be grouped together in this list and have
    exactly the same 'name' array so that the find_nothrow() function
    doesn't report the same constant multiple times.
  */
  list={{{"schwarzchildradius","rschwarz"},
	 "m",o2scl_const::o2scl_mks,o2scl_mks::schwarzchild_radius,
         "derived from the IAU 2015 nominal solar mass parameter"},
	{{"schwarzchildradius","rschwarz"},
	 "cm",o2scl_const::o2scl_cgs,o2scl_cgs::schwarzchild_radius,
         "derived from the IAU 2015 nominal solar mass parameter"},
	{{"schwarzchildradius","rschwarz"},
	 "km",o2scl_const::o2scl_mks,o2scl_mks::schwarzchild_radius/1.0e3,
         "derived from the IAU 2015 nominal solar mass parameter"},
	{{"speedoflight","c","lightspeed"},
	 "m/s",o2scl_const::o2scl_mks,
	 o2scl_const::speed_of_light_f<double>(o2scl_const::o2scl_mks),
         "exact"},
	{{"speedoflight","c","lightspeed"},
	 "cm/s",o2scl_const::o2scl_cgs,
	 o2scl_const::speed_of_light_f<double>(o2scl_const::o2scl_cgs),
         "exact"},
	{{"gravitationalconstant","g","newtonsconstant",
	  "newtonconstant"},"m^3/kg/s^2",o2scl_const::o2scl_mks,
	 o2scl_mks::gravitational_constant,"CODATA 2018"},
	{{"gravitationalconstant","g","newtonsconstant",
	  "newtonconstant"},"cm^3/g/s^2",o2scl_const::o2scl_cgs,
	 o2scl_cgs::gravitational_constant,"CODATA 2018"},
	{{"boltzmannconstant","kb","boltzmannsconstant","boltzmann"},
	 "m^2/kg/s^2/K",o2scl_const::o2scl_mks,o2scl_mks::boltzmann,"exact"},
	{{"boltzmannconstant","kb","boltzmannsconstant","boltzmann"},
	 "cm^2/g/s^2/K",o2scl_const::o2scl_cgs,o2scl_cgs::boltzmann,"exact"},
	{{"stefanboltzmannconstant","sigmasb","stefanboltzmann","ssb"},
	 "kg/s^3/K^4",o2scl_const::o2scl_mks,
	 o2scl_mks::stefan_boltzmann_constant,
         "exact; derived from k_B, c, and h bar"},
	{{"stefanboltzmannconstant","sigmasb","stefanboltzmann","ssb"},
	 "g/s^3/K^4",o2scl_const::o2scl_cgs,
	 o2scl_cgs::stefan_boltzmann_constant,
         "exact; derived from k_B, c, and h bar"},
	{{"plancksconstant","h","planckconstant","planck"},
	 "kg*m^2/s",o2scl_const::o2scl_mks,
	 o2scl_const::planck_f<double>(o2scl_const::o2scl_mks),"exact"},
	{{"plancksconstant","h","planckconstant","planck"},
	 "g*cm^2/s",o2scl_const::o2scl_cgs,
	 o2scl_const::planck_f<double>(o2scl_const::o2scl_cgs),"exact"},
	{{"reducedplancksconstant","hbar"},
	 "kg*m^2/s",o2scl_const::o2scl_mks,
	 o2scl_const::hbar_f<double>(o2scl_const::o2scl_mks),
         "exact; derived from the Planck constant"},
	{{"reducedplancksconstant","hbar"},
	 "g*cm^2/s",o2scl_const::o2scl_cgs,
	 o2scl_const::hbar_f<double>(o2scl_const::o2scl_cgs),
         "exact; derived from the Planck constant"},
	{{"avogadrosnumber","NA","avogadro"},
	 "",0,o2scl_const::avogadro,"exact"},
	{{"alphaem","finestructure","alpha"},"",0,
	 o2scl_const::fine_structure,"CODATA 2018"},
	{{"pi"},"",0,o2scl_const::pi,"exact"},
	{{"zeta32"},"",0,o2scl_const::zeta32,"exact"},
	{{"zeta2"},"",0,o2scl_const::zeta2,"exact"},
	{{"zeta52"},"",0,o2scl_const::zeta52,"exact"},
	{{"zeta3"},"",0,o2scl_const::zeta3,"exact"},
	{{"zeta5"},"",0,o2scl_const::zeta5,"exact"},
	{{"zeta7"},"",0,o2scl_const::zeta7,"exact"},
	{{"pi2","pisquared"},"",0,o2scl_const::pi2,"exact"},
	{{"rootpi","squarerootpi"},"",0,o2scl_const::root_pi,"exact"},
	{{"sin2thetaw"},"",0,o2scl_const::sin2_theta_weak,"PDG 2020 value"},
	{{"gfermi"},"s^4/m^4/kg^2",o2scl_const::o2scl_mks,
	 o2scl_mks::gfermi,
         ((string)"derived from CODATA 2018 value for G_Fermi (identical to ")+
         "PDG 2020 value) and CODATA 2018 value of electron volt"},
	{{"gfermi"},"s^4/cm^4/g^2",o2scl_const::o2scl_cgs,
	 o2scl_cgs::gfermi,
         ((string)"derived from CODATA 2018 value for G_Fermi (identical to ")+
         "PDG 2020 value) and CODATA 2018 value of electron volt"},
	{{"gfermi"},"1/GeV^2",0,o2scl_const::gfermi_gev2,
         "CODATA 2018 (identical to PDG 2020 value)"},
	{{"elementarycharge","electroncharge","e","chargeelectron",
	  "qelectron"},"C",
	 o2scl_const::o2scl_mks,o2scl_const::elem_charge_f<double>(),
         "exact"},
	{{"hbarc"},"MeV*fm",
	 0,o2scl_const::hc_mev_fm_f<double>(),"derived from Plack constant"},
	{{"hbarc"},"J*m",o2scl_const::o2scl_mks,
	 o2scl_const::hbarc_f<double>(o2scl_const::o2scl_mks),
         "derived from Plack constant"},
	{{"hbarc"},"erg*cm",o2scl_const::o2scl_cgs,
	 o2scl_const::hbarc_f<double>(o2scl_const::o2scl_cgs),
         "derived from Plack constant"},
	{{"masselectron","electronmass","melectron","melec"},"kg",
	 o2scl_const::o2scl_mks,o2scl_mks::mass_electron,"CODATA 2018"},
	{{"masselectron","electronmass","melectron","melec"},"g",
	 o2scl_const::o2scl_cgs,o2scl_cgs::mass_electron,"CODATA 2018"},
	{{"massmuon","muonmass","mmuon"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_muon,"CODATA 2018"},
	{{"massmuon","muonmass","mmuon"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_muon,"CODATA 2018"},
	{{"masstau","taumass","mtau"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_tau,"CODATA 2018"},
	{{"masstau","taumass","mtau"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_tau,"CODATA 2018"},
	{{"massneutron","neutronmass","mneutron","mneut"},"kg",
	 o2scl_const::o2scl_mks,o2scl_mks::mass_neutron,"CODATA 2018"},
	{{"massneutron","neutronmass","mneutron","mneut"},"g",
	 o2scl_const::o2scl_cgs,o2scl_cgs::mass_neutron,"CODATA 2018"},
	{{"massproton","protonmass","mproton","mprot"},"kg",
	 o2scl_const::o2scl_mks,o2scl_mks::mass_proton,"CODATA 2018"},
	{{"massproton","protonmass","mproton","mprot"},"g",
	 o2scl_const::o2scl_cgs,o2scl_cgs::mass_proton,"CODATA 2018"},
	{{"massdeuteron","deuteronmass","mdeuteron","mdeut"},"kg",
	 o2scl_const::o2scl_mks,o2scl_mks::mass_deuteron,"CODATA 2018"},
	{{"massdeuteron","deuteronmass","mdeuteron","mdeut"},"g",
	 o2scl_const::o2scl_cgs,o2scl_cgs::mass_deuteron,"CODATA 2018"},
	{{"masstriton","tritonmass","mtriton"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_triton,"CODATA 2018"},
	{{"masstriton","tritonmass","mtriton"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_triton,"CODATA 2018"},
	{{"masshelion","helionmass","mhelion"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_helion,"CODATA 2018"},
	{{"masshelion","helionmass","mhelion"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_helion,"CODATA 2018"},
	{{"massalpha","alphamass","malpha"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_alpha,"CODATA 2018"},
	{{"massalpha","alphamass","malpha"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_alpha,"CODATA 2018"},
	{{"masslambda","lambdamass","mlambda"},"MeV",0,
	 o2scl_const::mass_lambda_MeV,"\"OUR FIT\" value from PDG 2020"},
	{{"masssigmaminus","sigmaminusmass","msigma-"},"MeV",0,
	 o2scl_const::mass_sigma_minus_MeV,"\"OUR FIT\" value from PDG 2020"},
	{{"masssigmazero","sigmazeromass","msigma0"},"MeV",0,
	 o2scl_const::mass_sigma_zero_MeV,"\"OUR FIT\" value from PDG 2020"},
	{{"masssigmaplus","sigmaplusmass","msigma+"},"MeV",0,
	 o2scl_const::mass_sigma_plus_MeV,"\"OUR FIT\" value from PDG 2020"},
	{{"masscascadezero","cascadezeromass","mcascade0","mxi0"},
	 "MeV",0,o2scl_const::mass_cascade_zero_MeV,
         "\"OUR FIT\" value from PDG 2020"},
	{{"masscascademinus","cascademinusmass","mcascade-","mxi-"},
	 "MeV",0,o2scl_const::mass_cascade_minus_MeV,
         "\"OUR FIT\" value from PDG 2020"},
	{{"massup","upmass","mup"},"MeV",0,
	 o2scl_const::mass_up_MeV,"\"OUR EVALUATION\" value from PDG 2020"},
	{{"massdown","downmass","mdown"},"MeV",0,
	 o2scl_const::mass_down_MeV,"\"OUR EVALUATION\" value from PDG 2020"},
	{{"massstrange","strangemass","mstrange"},"MeV",0,
	 o2scl_const::mass_strange_MeV,
         "\"OUR EVALUATION\" value from PDG 2020"},
	{{"masssolar","solarmass","masssun","sunmass","msun","modot"},
	 "kg",o2scl_const::o2scl_mks,o2scl_mks::solar_mass,
         ((string)"derived from IAU's 2015 nominal value of the solar ")
         +"mass parameter divided by the CODATA 2018 value of the "+
         "gravitational constant"},
	{{"masssolar","solarmass","masssun","sunmass","msun","modot"},
	 "g",o2scl_const::o2scl_cgs,o2scl_cgs::solar_mass,
         ((string)"derived from IAU's 2015 nominal value of the solar ")
         +"mass parameter divided by the CODATA 2018 value of the "+
         "gravitational constant"},
	{{"massmercury","mercurymass","mmercury"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mercury_mass,""},
	{{"massmercury","mercurymass","mmercury"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mercury_mass,""},
	{{"massvenus","venusmass","mvenus"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::venus_mass,""},
	{{"massvenus","venusmass","mvenus"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::venus_mass,""},
	{{"massearth","earthmass","mearth"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::earth_mass,""},
	{{"massearth","earthmass","mearth"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::earth_mass,""},
	{{"massmars","marsmass","mmars"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mars_mass,""},
	{{"massmars","marsmass","mmars"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mars_mass,""},
	{{"massjupiter","jupitermass","mjupiter","mjup"},
	 "kg",o2scl_const::o2scl_mks,o2scl_mks::jupiter_mass,""},
	{{"massjupiter","jupitermass","mjupiter","mjup"},
	 "g",o2scl_const::o2scl_cgs,o2scl_cgs::jupiter_mass,""},
	{{"masssaturn","saturnmass","msaturn"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::saturn_mass,""},
	{{"masssaturn","saturnmass","msaturn"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::saturn_mass,""},
	{{"massuranus","uranusmass","muranus"},"kg",
	 o2scl_const::o2scl_mks,o2scl_mks::uranus_mass,""},
	{{"massuranus","uranusmass","muranus"},"g",
	 o2scl_const::o2scl_cgs,o2scl_cgs::uranus_mass,""},
	{{"massneptune","neptunemass","mneptune"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::neptune_mass,""},
	{{"massneptune","neptunemass","mneptune"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::neptune_mass,""},
	{{"masspluto","plutomass","mpluto"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::pluto_mass,""},
	{{"masspluto","plutomass","mpluto"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::pluto_mass,""},
	{{"radiussolar","solarradius","radiussun","sunradius","rsun"},
	 "m",o2scl_const::o2scl_mks,o2scl_mks::solar_radius,""},
	{{"radiussolar","solarradius","radiussun","sunradius","rsun"},
	 "cm",o2scl_const::o2scl_cgs,o2scl_cgs::solar_radius,""},
	{{"radiusmercury","mercuryradius","rmercury"},"m",
	 o2scl_const::o2scl_mks,o2scl_mks::mercury_radius,""},
	{{"radiusmercury","mercuryradius","rmercury"},"cm",
	 o2scl_const::o2scl_cgs,o2scl_cgs::mercury_radius,""},
	{{"radiusvenus","venusradius","rvenus"},"m",
	 o2scl_const::o2scl_mks,o2scl_mks::venus_radius,""},
	{{"radiusvenus","venusradius","rvenus"},"cm",
	 o2scl_const::o2scl_cgs,o2scl_cgs::venus_radius,""},
	{{"radiusearthequatorial","earthequatorialradius",
	  "earthradiusequatorial"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::earth_radius_equatorial,""},
	{{"radiusearthequatorial","earthequatorialradius",
	  "earthradiusequatorial"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::earth_radius_equatorial,""},
	{{"radiusearthpolar","earthpolarradius",
	  "earthradiuspolar"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::earth_radius_polar,""},
	{{"radiusearthpolar","earthpolarradius",
	  "earthradiuspolar"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::earth_radius_polar,""},
	{{"radiusmarsequatorial","marsequatorialradius",
	  "marsradiusequatorial"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::mars_radius_equatorial,""},
	{{"radiusmarsequatorial","marsequatorialradius",
	  "marsradiusequatorial"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mars_radius_equatorial,""},
	{{"radiusmarspolar","marspolarradius",
	  "marsradiuspolar"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::mars_radius_polar,""},
	{{"radiusmarspolar","marspolarradius",
	  "marsradiuspolar"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mars_radius_polar,""},
	{{"radiusjupiterequatorial","jupiterequatorialradius",
	  "jupiterradiusequatorial"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::jupiter_radius_equatorial,""},
	{{"radiusjupiterequatorial","jupiterequatorialradius",
	  "jupiterradiusequatorial"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::jupiter_radius_equatorial,""},
	{{"radiusjupiterpolar","jupiterpolarradius",
	  "jupiterradiuspolar"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::jupiter_radius_polar,""},
	{{"radiusjupiterpolar","jupiterpolarradius",
	  "jupiterradiuspolar"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::jupiter_radius_polar,""},
	{{"radiussaturnequatorial","saturnequatorialradius",
	  "saturnradiusequatorial"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::saturn_radius_equatorial,""},
	{{"radiussaturnequatorial","saturnequatorialradius",
	  "saturnradiusequatorial"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::saturn_radius_equatorial,""},
	{{"radiussaturnpolar","saturnpolarradius",
	  "saturnradiuspolar"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::saturn_radius_polar,""},
	{{"radiussaturnpolar","saturnpolarradius",
	  "saturnradiuspolar"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::saturn_radius_polar,""},
	{{"radiusuranusequatorial","uranusequatorialradius",
	  "uranusradiusequatorial"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::uranus_radius_equatorial,""},
	{{"radiusuranusequatorial","uranusequatorialradius",
	  "uranusradiusequatorial"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::uranus_radius_equatorial,""},
	{{"radiusuranuspolar","uranuspolarradius",
	  "uranusradiuspolar"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::uranus_radius_polar,""},
	{{"radiusuranuspolar","uranuspolarradius",
	  "uranusradiuspolar"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::uranus_radius_polar,""},
	{{"radiusneptuneequatorial","neptuneequatorialradius",
	  "neptuneradiusequatorial"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::neptune_radius_equatorial,""},
	{{"radiusneptuneequatorial","neptuneequatorialradius",
	  "neptuneradiusequatorial"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::neptune_radius_equatorial,""},
	{{"radiusneptunepolar","neptunepolarradius",
	  "neptuneradiuspolar"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::neptune_radius_polar,""},
	{{"radiusneptunepolar","neptunepolarradius",
	  "neptuneradiuspolar"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::neptune_radius_polar,""},
	{{"radiuspluto","plutoradius","rpluto"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::pluto_radius,""},
	{{"radiuspluto","plutoradius","rpluto"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::pluto_radius,""},
	{{"rydberg"},"kg*m^2/s^2",o2scl_const::o2scl_mks,
	 o2scl_mks::rydberg,"CODATA 2018"},
	{{"rydberg"},"g*cm^2/s^2",o2scl_const::o2scl_cgs,
	 o2scl_cgs::rydberg,"CODATA 2018"}
  };

}

int find_constants::find_nothrow(std::string name, std::string unit,
				 vector<find_constants_list> &matches,
				 int verbose) {
  
  o2scl::convert_units<> &cu=o2scl_settings.get_convert_units();
  
  // Simplify by removing alphanumerics except + and -,
  // which we need to distinguish between positive and negative
  // particle masses
  if (verbose>1) {
    std::cout << "find_constants::find_nothrow(): "
	      << "before simplify: " << name << endl;
  }
  for(size_t i=0;i<name.length();i++) {
    if (!isalnum(name[i]) && name[i]!='+' && name[i]!='-') {
      name.erase(i,1);
      i=0;
    }
  }
  if (verbose>1) {
    std::cout << "find_constants::find_nothrow(): "
	      << "after simplify: " << name << endl;
  }
    
  // Start with a fresh list
  matches.clear();
    
  // Temporarily store matching indexes
  vector<size_t> indexes;

  int match_type=0, match_exact=1, match_pattern=2;
    
  // Initial pass, exact matches
  for(size_t i=0;i<list.size();i++) {
    for(size_t j=0;j<list[i].names.size();j++) {
      if (verbose>2) {
	std::cout << "find_constants::find_nothrow(): "
		  << name << " " << i << " " << j << " "
		  << list[i].names[j]
		  << " " << boost::iequals(name,list[i].names[j]) << endl;
      }
      if (boost::iequals(name,list[i].names[j])) {
	indexes.push_back(i);
	// Now that we've found a match, don't look in the
	// other names for this list entry
	j=list[i].names.size();
	match_type=match_exact;
      }
    }
  }

  if (verbose>1) {
    std::cout << "find_constants::find_nothrow(): "
	      << "pass 1 indexes: ";
    vector_out(std::cout,indexes,true);
  }
  
  // No matches, so try wildcard matches
  if (indexes.size()==0) {
      
    string name_wc=((string)"*")+name+"*";
    for(size_t i=0;i<list.size();i++) {
      for(size_t j=0;j<list[i].names.size();j++) {
	int fn_ret=fnmatch(name_wc.c_str(),(list[i].names[j]).c_str(),
			   FNM_CASEFOLD);
	if (verbose>2) {
	  std::cout << "find_constants::find_nothrow(): "
		    << name_wc << " " << i << " " << j << " "
		    << list[i].names[j]
		    << " " << fn_ret << endl;
	}
	if (fn_ret==0) {
	  indexes.push_back(i);
	  // Now that we've found a match, don't look in the
	  // other names for this list entry
	  j=list[i].names.size();
	  match_type=match_pattern;
	}
      }
    }
      
  }

  if (verbose>1) {
    std::cout << "find_constants::find_nothrow(): "
	      << "pass 2 indexes: ";
    vector_out(std::cout,indexes,true);
  }
  
  // There was only one match
  if (indexes.size()==1) {

    // Add to 'matches' list
    matches.push_back(list[indexes[0]]);

    if (verbose>1) {
      std::cout << "find_constants::find_nothrow(): "
		<< "one match unit: " << unit << " "
		<< list[indexes[0]].unit_flag << " "
		<< list[indexes[0]].unit << std::endl;
    }
  
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
      if (verbose>0) {
	cout << "find_constant::find_nothrow(): Trying to convert from "
	     << list[indexes[0]].unit << " to "
	     << unit << endl;
      }
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
    
    if (verbose>1) {
      std::cout << "find_constants::find_nothrow(): "
		<< "Multiple matches found. No unit given." << std::endl;
    }

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

    if (verbose>1) {
      std::cout << "find_constants::find_nothrow(): "
		<< "Multiple name matches found. Checking units." << std::endl;
    }
    
    // We found at least one match, check unit
      
    vector<size_t> indexes2;
    
    // Look for entries with matching unit
    for(size_t i=0;i<indexes.size();i++) {

      if (verbose>1) {
	std::cout << "find_constants::find_nothrow(): "
		  << "many name matches unit: " << unit << " "
		  << list[indexes[i]].unit_flag << " "
		  << list[indexes[i]].unit << std::endl;
      }
      
      if ((unit=="cgs" &&
	   list[indexes[i]].unit_flag==o2scl_const::o2scl_cgs) ||
	  (unit=="mks" &&
	   list[indexes[i]].unit_flag==o2scl_const::o2scl_mks) ||
	  boost::iequals(list[indexes[i]].unit,unit)) {
	indexes2.push_back(indexes[i]);
	if (verbose>2) {
	  std::cout << "find_constants::find_nothrow(): Added."
		    << std::endl;
	}
      }
    }
      
    if (indexes2.size()==0) {

      if (verbose>1) {
	std::cout << "find_constants::find_nothrow(): "
		  << "many name matches and unit " << unit
                  << " specified, "
                  << "but no unit matches." << std::endl;
      }
      
      // No matching unit, try to convert
      for(size_t i=0;i<indexes.size();i++) {
	double val2;
	cout << "Trying to convert from "
	     << list[indexes[i]].unit << " to "
	     << unit << endl;
	int cret=cu.convert_ret(list[indexes[i]].unit,unit,
				list[indexes[i]].val,val2);
	if (cret==0 &&
	    (matches.size()==0 ||
	     list[indexes[i]].names!=matches[matches.size()-1].names)) {
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
	if (i==0 ||
	    list[indexes[i]].names!=matches[matches.size()-1].names) {
	  matches.push_back(list[indexes[i]]);
	}
      }
      if (match_type==match_exact) {
	return exact_matches_unit_diff;
      } else {
	return pattern_matches_unit_diff;
      }

    } else {

      if (verbose>1) {
        std::cout << "At least one exact unit match was found." << std::endl;
      }

      // There were exact unit matches, so set up the matches list
      for(size_t i=0;i<indexes2.size();i++) {
	if (i==0 ||
	    list[indexes2[i]].names!=matches[matches.size()-1].names) {
	  matches.push_back(list[indexes2[i]]);
	}
      }
      if (match_type==match_exact) {
        if (matches.size()==1) {
          return one_exact_match_unit_match;
        } else {
          return exact_matches_unit_match;
        }
      } else {
        if (matches.size()==1) {
          return one_pattern_match_unit_match;
        } else {
          return pattern_matches_unit_match;
        }
      }
    }
  }

  return no_matches;
}

void find_constants::find_print(std::string name, std::string unit,
				size_t prec, int verbose) {

  cout.precision(prec);
    
  vector<find_constants_list> matches;
  int ret=find_nothrow(name,unit,matches,verbose);
  if (ret==no_matches) {
    cout << "find_constant::find_print(): No matches found for name "
	 << name << endl;
    return;
  }
  
  cout << "find_constant::find_print(): Matches for " << name;
  if (ret==one_exact_match_unit_diff ||
      ret==exact_matches_unit_diff) {
    cout << " (no matching units)" << endl;
  } else if (unit.length()>0) {
    cout << " in " << unit;
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
    std::string err=((std::string)"Failed to find unique match for name ")+
      name+" and unit "+unit+" in find_constants::find_unique(). "+
      "Returned "+o2scl::itos(ret)+".";
    O2SCL_ERR(err.c_str(),o2scl::exc_einval);
  }
  return matches[0].val;
}

void find_constants::output_list(std::ostream &os) {
  os << "name unit flag value source" << endl;
  os << "  alternate names" << endl;
  for(size_t i=0;i<list.size();i++) {
    os << list[i].names[0] << " ";
    os << list[i].unit << " ";
    os << list[i].unit_flag << " ";
    os << list[i].val << " ";
    os << list[i].source << endl;
    if (list[i].names.size()>0) {
      os << "  ";
      for(size_t j=1;j<list[i].names.size();j++) {
        os << list[i].names[j] << " ";
      }
      os << endl;
    } else {
      os << "  (no alternate names)" << endl;
    }
  }
  return;
}
