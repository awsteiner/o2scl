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
#ifdef O2SCL_REGEX
#include <regex>
#else
#include <fnmatch.h>
#endif

#include <boost/algorithm/string.hpp>

#include <o2scl/find_constants.h>
#include <o2scl/lib_settings.h>
#include <o2scl/convert_units.h>
#include <o2scl/vector.h>

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
  list={{{"Schwarzchild radius","rschwarz"},
	 "m",o2scl_const::o2scl_mks,o2scl_mks::schwarzchild_radius,
         "derived from the IAU 2015 nominal solar mass parameter",
         1,0,0,0,0,0,0},
	{{"Schwarzchild radius","rschwarz"},
	 "cm",o2scl_const::o2scl_cgs,o2scl_cgs::schwarzchild_radius,
         "derived from the IAU 2015 nominal solar mass parameter",
         0,0,0,0,0,0,0},
	{{"Schwarzchild radius","rschwarz"},
	 "km",o2scl_const::o2scl_mks,o2scl_mks::schwarzchild_radius/1.0e3,
         "derived from the IAU 2015 nominal solar mass parameter",
         1,0,0,0,0,0,0},
	{{"speed of light","c","lightspeed"},
	 "m/s",o2scl_const::o2scl_mks,
	 o2scl_const::speed_of_light_f<double>(o2scl_const::o2scl_mks),
         "exact",1,0,-1,0,0,0,0},
	{{"speed of light","c","lightspeed"},
	 "cm/s",o2scl_const::o2scl_cgs,
	 o2scl_const::speed_of_light_f<double>(o2scl_const::o2scl_cgs),
         "exact",0,0,0,0,0,0,0},
	{{"gravitational","g","newtonsconstant","gn",
	  "newtonconstant"},"m^3/kg/s^2",o2scl_const::o2scl_mks,
	 o2scl_mks::gravitational_constant,"CODATA 2018",3,-1,-2,0,0,0,0},
	{{"gravitational","g","newtonsconstant","gn",
	  "newtonconstant"},"cm^3/g/s^2",o2scl_const::o2scl_cgs,
	 o2scl_cgs::gravitational_constant,"CODATA 2018",0,0,0,0,0,0,0},
	{{"Boltzmann's","kb","boltzmannsconstant","boltzmann"},
	 "m^2/kg/s^2/K",o2scl_const::o2scl_mks,o2scl_mks::boltzmann,
         "exact",2,-1,-2,-1,0,0,0},
	{{"Boltzmann's","kb","boltzmannsconstant","boltzmann"},
	 "cm^2/g/s^2/K",o2scl_const::o2scl_cgs,o2scl_cgs::boltzmann,
         "exact",0,0,0,0,0,0,0},
	{{"Stefan-Boltzmann","sigmasb","stefanboltzmann","ssb","σsb"},
	 "kg/s^3/K^4",o2scl_const::o2scl_mks,
	 o2scl_mks::stefan_boltzmann_constant,
         "exact; derived from k_B, c, and h bar",0,1,-3,-4,0,0,0},
	{{"Stefan-Boltzmann","sigmasb","stefanboltzmann","ssb","σsb"},
	 "g/s^3/K^4",o2scl_const::o2scl_cgs,
	 o2scl_cgs::stefan_boltzmann_constant,
         "exact; derived from k_B, c, and h bar",0,0,0,0,0,0,0},
	{{"Planck","h","planck","plancks"},
	 "kg*m^2/s",o2scl_const::o2scl_mks,
	 o2scl_const::planck_f<double>(o2scl_const::o2scl_mks),
         "exact",2,1,-1,0,0,0,0},
	{{"Planck","h","planckconstant","plancks"},
	 "g*cm^2/s",o2scl_const::o2scl_cgs,
	 o2scl_const::planck_f<double>(o2scl_const::o2scl_cgs),
         "exact",0,0,0,0,0,0,0},
	{{"reduced Planck","hbar","ħ","reducedplancks"},
	 "kg*m^2/s",o2scl_const::o2scl_mks,
	 o2scl_const::hbar_f<double>(o2scl_const::o2scl_mks),
         "exact; derived from the Planck constant",2,1,-1,0,0,0,0},
	{{"reduced Planck","hbar","ħ","reducedplancks"},
	 "g*cm^2/s",o2scl_const::o2scl_cgs,
	 o2scl_const::hbar_f<double>(o2scl_const::o2scl_cgs),
         "exact; derived from the Planck constant",0,0,0,0,0,0,0},
	{{"Avogadro's number","na","avogadro"},
	 "",fc_none,o2scl_const::avogadro,"exact",0,0,0,0,0,0,0},
	{{"Fine structure","alphaem","alpha","αem"},"",0,
	 o2scl_const::fine_structure,"CODATA 2018",0,0,0,0,0,0,0},
	{{"pi","π"},"",fc_none,o2scl_const::pi,"exact",0,0,0,0,0,0,0},
	{{"zeta32","zeta(3/2)","ζ(3/2)"},"",fc_none,o2scl_const::zeta32,
         "exact",0,0,0,0,0,0,0},
	{{"zeta2","zeta(2)","ζ(2)"},"",fc_none,o2scl_const::zeta2,
         "exact",0,0,0,0,0,0,0},
	{{"zeta52","zeta(5/2)","ζ(5/2)"},"",fc_none,o2scl_const::zeta52,
         "exact",0,0,0,0,0,0,0},
	{{"zeta3","zeta(3)","ζ(3)"},"",fc_none,o2scl_const::zeta3,
         "exact",0,0,0,0,0,0,0},
	{{"zeta5","zeta(5)","ζ(5)"},"",fc_none,o2scl_const::zeta5,
         "exact",0,0,0,0,0,0,0},
	{{"zeta7","zeta(7)","ζ(7)"},"",fc_none,o2scl_const::zeta7,
         "exact",0,0,0,0,0,0,0},
	{{"pi2","pisquared","π²"},"",fc_none,o2scl_const::pi2,
         "exact",0,0,0,0,0,0,0},
	{{"pi2","pisquared","π²"},"",fc_none,o2scl_const::pi2,
         "exact",0,0,0,0,0,0,0},
	{{"pi2","pisquared","π²"},"",fc_none,o2scl_const::pi2,
         "exact",0,0,0,0,0,0,0},
	{{"pi3","picubed","π³"},"",fc_none,o2scl_const::pi2*o2scl_const::pi,
         "exact",0,0,0,0,0,0,0},
	{{"pi4","pifourth","π⁴"},"",fc_none,o2scl_const::pi2*o2scl_const::pi2,
         "exact",0,0,0,0,0,0,0},
	{{"rootpi","squarerootpi","√π"},"",fc_none,o2scl_const::root_pi,
         "exact",0,0,0,0,0,0,0},
	{{"Euler-Mascheroni","euler"},"",fc_none,
         o2scl_const::euler_mascheroni,"exact",0,0,0,0,0,0,0},
	{{"sin2thetaw","sin2θW","sin²θW"},"",fc_none,
         o2scl_const::sin2_theta_weak,"PDG 2020 value",0,0,0,0,0,0,0},
	{{"gfermi","gf"},"s^4/m^4/kg^2",o2scl_const::o2scl_mks,
	 o2scl_mks::gfermi,
         ((string)"derived from CODATA 2018 value for G_Fermi (identical to ")+
         "PDG 2020 value) and CODATA 2018 value of electron volt",
         -4,-2,4,0,0,0,0},
	{{"gfermi","gf"},"s^4/cm^4/g^2",o2scl_const::o2scl_cgs,
	 o2scl_cgs::gfermi,
         ((string)"derived from CODATA 2018 value for G_Fermi (identical to ")+
         "PDG 2020 value) and CODATA 2018 value of electron volt",
         0,0,0,0,0,0,0},
	{{"gfermi","gf"},"1/GeV^2",0,o2scl_const::gfermi_gev2,
         "CODATA 2018 (identical to PDG 2020 value)",0,0,0,0,0,0,0},
	{{"elementarycharge","electroncharge","e","chargeelectron",
	  "qelectron"},"C",
	 o2scl_const::o2scl_mks,o2scl_const::elem_charge_f<double>(),
         "exact",0,0,1,0,1,0,0},
	{{"hbarc","ħc"},"MeV*fm",
	 0,o2scl_const::hc_mev_fm_f<double>(),
         "derived from Plack constant",0,0,0,0,0,0,0},
	{{"hbarc","ħc"},"J*m",o2scl_const::o2scl_mks,
	 o2scl_const::hbarc_f<double>(o2scl_const::o2scl_mks),
         "derived from Plack constant",3,1,-2,0,0,0,0},
	{{"hbarc","ħc"},"erg*cm",o2scl_const::o2scl_cgs,
	 o2scl_const::hbarc_f<double>(o2scl_const::o2scl_cgs),
         "derived from Plack constant",0,0,0,0,0,0,0},
	{{"mass electron","electronmass","melectron","melec"},"kg",
	 o2scl_const::o2scl_mks,o2scl_mks::mass_electron,
         "CODATA 2018",0,1,0,0,0,0,0},
	{{"mass electron","electronmass","melectron","melec"},"g",
	 o2scl_const::o2scl_cgs,o2scl_cgs::mass_electron,
         "CODATA 2018",0,0,0,0,0,0,0},
	{{"mass muon","muonmass","mmuon"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_muon,"CODATA 2018",0,1,0,0,0,0,0},
	{{"mass muon","muonmass","mmuon"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_muon,"CODATA 2018",0,0,0,0,0,0,0},
	{{"mass tau","taumass","mtau"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_tau,"CODATA 2018",0,1,0,0,0,0,0},
	{{"mass tau","taumass","mtau"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_tau,"CODATA 2018",0,0,0,0,0,0,0},
	{{"mass neutron","neutronmass","mneutron","mneut"},"kg",
	 o2scl_const::o2scl_mks,o2scl_mks::mass_neutron,
         "CODATA 2018",0,1,0,0,0,0,0},
	{{"mass neutron","neutronmass","mneutron","mneut"},"g",
	 o2scl_const::o2scl_cgs,o2scl_cgs::mass_neutron,
         "CODATA 2018",0,0,0,0,0,0,0},
	{{"mass proton","protonmass","mproton","mprot"},"kg",
	 o2scl_const::o2scl_mks,o2scl_mks::mass_proton,
         "CODATA 2018",0,1,0,0,0,0,0},
	{{"mass proton","protonmass","mproton","mprot"},"g",
	 o2scl_const::o2scl_cgs,o2scl_cgs::mass_proton,
         "CODATA 2018",0,0,0,0,0,0,0},
	{{"mass deuteron","deuteronmass","mdeuteron","mdeut"},"kg",
	 o2scl_const::o2scl_mks,o2scl_mks::mass_deuteron,
         "CODATA 2018",0,1,0,0,0,0,0},
	{{"mass deuteron","deuteronmass","mdeuteron","mdeut"},"g",
	 o2scl_const::o2scl_cgs,o2scl_cgs::mass_deuteron,
         "CODATA 2018",0,0,0,0,0,0,0},
	{{"mass triton","tritonmass","mtriton"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_triton,"CODATA 2018",0,1,0,0,0,0,0},
	{{"mass triton","tritonmass","mtriton"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_triton,"CODATA 2018",0,0,0,0,0,0,0},
	{{"mass helion","helionmass","mhelion"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_helion,"CODATA 2018",0,1,0,0,0,0,0},
	{{"mass helion","helionmass","mhelion"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_helion,"CODATA 2018",0,0,0,0,0,0,0},
	{{"mass alpha","alphamass","malpha","mα"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mass_alpha,"CODATA 2018",0,1,0,0,0,0,0},
	{{"mass alpha","alphamass","malpha","mα"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mass_alpha,"CODATA 2018",0,0,0,0,0,0,0},
	{{"mass lambda","lambdamass","mlambda","mΛ"},"MeV",0,
	 o2scl_const::mass_lambda_MeV,
         "\"OUR FIT\" value from PDG 2020",0,0,0,0,0,0,0},
	{{"mass sigma minus","sigmaminusmass","msigma-","mΣ-"},"MeV",0,
	 o2scl_const::mass_sigma_minus_MeV,
         "\"OUR FIT\" value from PDG 2020",0,0,0,0,0,0,0},
	{{"mass sigma zero","sigmazeromass","msigma0","mΣ0"},"MeV",0,
	 o2scl_const::mass_sigma_zero_MeV,
         "\"OUR FIT\" value from PDG 2020",0,0,0,0,0,0,0},
	{{"mass sigma plus","sigmaplusmass","msigma+","mΣ+"},"MeV",0,
	 o2scl_const::mass_sigma_plus_MeV,
         "\"OUR FIT\" value from PDG 2020",0,0,0,0,0,0,0},
	{{"mass cascade zero","cascadezeromass","mcascade0","mxi0","mΞ0"},
	 "MeV",0,o2scl_const::mass_cascade_zero_MeV,
         "\"OUR FIT\" value from PDG 2020",0,0,0,0,0,0,0},
	{{"mass cascade minus","cascademinusmass","mcascade-","mxi-","mΞ-"},
	 "MeV",0,o2scl_const::mass_cascade_minus_MeV,
         "\"OUR FIT\" value from PDG 2020",0,0,0,0,0,0,0},
	{{"mass up","upmass","mup"},"MeV",0,
	 o2scl_const::mass_up_MeV,
         "\"OUR EVALUATION\" value from PDG 2020",0,0,0,0,0,0,0},
	{{"mass down","downmass","mdown"},"MeV",0,
	 o2scl_const::mass_down_MeV,
         "\"OUR EVALUATION\" value from PDG 2020",0,0,0,0,0,0,0},
	{{"mass strange","strangemass","mstrange"},"MeV",0,
	 o2scl_const::mass_strange_MeV,
         "\"OUR EVALUATION\" value from PDG 2020",0,0,0,0,0,0,0},
	{{"mass solar","solarmass","masssun","sunmass","msun","modot","m☉"},
	 "kg",o2scl_const::o2scl_mks,o2scl_mks::solar_mass,
         ((string)"derived from IAU's 2015 nominal value of the solar ")
         +"mass parameter divided by the CODATA 2018 value of the "+
         "gravitational constant",0,1,0,0,0,0,0},
	{{"mass solar","solarmass","masssun","sunmass","msun","modot","m☉"},
	 "g",o2scl_const::o2scl_cgs,o2scl_cgs::solar_mass,
         ((string)"derived from IAU's 2015 nominal value of the solar ")
         +"mass parameter divided by the CODATA 2018 value of the "+
         "gravitational constant",0,0,0,0,0,0,0},
	{{"mass mercury","mercurymass","mmercury","m☿"},
         "kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mercury_mass,"",0,1,0,0,0,0,0},
	{{"mass mercury","mercurymass","mmercury","m☿"},
         "g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mercury_mass,"",0,0,0,0,0,0,0},
	{{"mass venus","venusmass","mvenus","m♀"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::venus_mass,"",0,1,0,0,0,0,0},
	{{"mass venus","venusmass","mvenus","m♀"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::venus_mass,"",0,0,0,0,0,0,0},
	{{"mass earth","earthmass","mearth","m♁","m⊕"},
         "kg",o2scl_const::o2scl_mks,
	 o2scl_mks::earth_mass,"IAU 2015 nominal value",0,1,0,0,0,0,0},
	{{"mass earth","earthmass","mearth","m♁","m⊕"},
         "g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::earth_mass,"IAU 2015 nominal value",0,0,0,0,0,0,0},
	{{"mass mars","marsmass","mmars","m♂"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::mars_mass,"",0,1,0,0,0,0,0},
	{{"mass mars","marsmass","mmars","m♂"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mars_mass,"",0,0,0,0,0,0,0},
	{{"mass jupiter","jupitermass","mjupiter","mjup","m♃"},
	 "kg",o2scl_const::o2scl_mks,o2scl_mks::jupiter_mass,
         "IAU 2015 nominal value",0,1,0,0,0,0,0},
	{{"mass jupiter","jupitermass","mjupiter","mjup","m♃"},
	 "g",o2scl_const::o2scl_cgs,o2scl_cgs::jupiter_mass,
         "IAU 2015 nominal value",0,0,0,0,0,0,0},
	{{"mass saturn","saturnmass","msaturn","m♄"},"kg",
         o2scl_const::o2scl_mks,
	 o2scl_mks::saturn_mass,"",0,1,0,0,0,0,0},
	{{"mass saturn","saturnmass","msaturn","m♄"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::saturn_mass,"",0,0,0,0,0,0,0},
	{{"mass uranus","uranusmass","muranus","m♅"},"kg",
	 o2scl_const::o2scl_mks,o2scl_mks::uranus_mass,"",0,1,0,0,0,0,0},
	{{"mass uranus","uranusmass","muranus","m♅"},"g",
	 o2scl_const::o2scl_cgs,o2scl_cgs::uranus_mass,"",0,0,0,0,0,0,0},
	{{"mass neptune","neptunemass","mneptune","m♆"},"kg",
         o2scl_const::o2scl_mks,
	 o2scl_mks::neptune_mass,"",0,1,0,0,0,0,0},
	{{"mass neptune","neptunemass","mneptune","m♆"},"g",
         o2scl_const::o2scl_cgs,
	 o2scl_cgs::neptune_mass,"",0,0,0,0,0,0,0},
	{{"mass pluto","plutomass","mpluto","m♇"},"kg",o2scl_const::o2scl_mks,
	 o2scl_mks::pluto_mass,"",0,1,0,0,0,0,0},
	{{"masspluto","plutomass","mpluto","m♇"},"g",o2scl_const::o2scl_cgs,
	 o2scl_cgs::pluto_mass,"",0,0,0,0,0,0,0},
	{{"radius solar","solarradius","radiussun","sunradius","rsun","r☉"},
	 "m",o2scl_const::o2scl_mks,o2scl_mks::solar_radius,"",1,0,0,0,0,0,0},
	{{"radius solar","solarradius","radiussun","sunradius","rsun","r☉"},
	 "cm",o2scl_const::o2scl_cgs,o2scl_cgs::solar_radius,"",0,0,0,0,0,0,0},
	{{"radius mercury","mercuryradius","rmercury","r☿"},"m",
	 o2scl_const::o2scl_mks,o2scl_mks::mercury_radius,"",1,0,0,0,0,0,0},
	{{"radius mercury","mercuryradius","rmercury","r☿"},"cm",
	 o2scl_const::o2scl_cgs,o2scl_cgs::mercury_radius,"",0,0,0,0,0,0,0},
	{{"radius venus","venusradius","rvenus","r♀"},"m",
	 o2scl_const::o2scl_mks,o2scl_mks::venus_radius,"",1,0,0,0,0,0,0},
	{{"radius venus","venusradius","rvenus","r♀"},"cm",
	 o2scl_const::o2scl_cgs,o2scl_cgs::venus_radius,"",0,0,0,0,0,0,0},
	{{"radius earthequatorial","earthequatorialradius",
	  "earthradiusequatorial","r♁eq","r⊕eq"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::earth_radius_equatorial,"",1,0,0,0,0,0,0},
	{{"radius earthequatorial","earthequatorialradius",
	  "earthradiusequatorial","r♁eq","r⊕eq"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::earth_radius_equatorial,"",0,0,0,0,0,0,0},
	{{"radius earthpolar","earthpolarradius",
	  "earthradiuspolar","r♁pol","r⊕pol"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::earth_radius_polar,"",1,0,0,0,0,0,0},
	{{"radius earthpolar","earthpolarradius",
	  "earthradiuspolar","r♁pol","r⊕pol"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::earth_radius_polar,"",0,0,0,0,0,0,0},
	{{"radius marsequatorial","marsequatorialradius",
	  "marsradiusequatorial","r♂eq"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::mars_radius_equatorial,"",1,0,0,0,0,0,0},
	{{"radius marsequatorial","marsequatorialradius",
	  "marsradiusequatorial","r♂eq"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mars_radius_equatorial,"",0,0,0,0,0,0,0},
	{{"radius marspolar","marspolarradius",
	  "marsradiuspolar","r♂pol"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::mars_radius_polar,"",1,0,0,0,0,0,0},
	{{"radius marspolar","marspolarradius",
	  "marsradiuspolar","r♂pol"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::mars_radius_polar,"",0,0,0,0,0,0,0},
	{{"radius jupiterequatorial","jupiterequatorialradius",
	  "jupiterradiusequatorial","r♃eq"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::jupiter_radius_equatorial,"",1,0,0,0,0,0,0},
	{{"radius jupiterequatorial","jupiterequatorialradius",
	  "jupiterradiusequatorial","r♃eq"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::jupiter_radius_equatorial,"",0,0,0,0,0,0,0},
	{{"radius jupiterpolar","jupiterpolarradius",
	  "jupiterradiuspolar","r♃pol"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::jupiter_radius_polar,"",1,0,0,0,0,0,0},
	{{"radius jupiterpolar","jupiterpolarradius",
	  "jupiterradiuspolar","r♃pol"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::jupiter_radius_polar,"",0,0,0,0,0,0,0},
	{{"radius saturnequatorial","saturnequatorialradius",
	  "saturnradiusequatorial","r♄eq"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::saturn_radius_equatorial,"",1,0,0,0,0,0,0},
	{{"radius saturnequatorial","saturnequatorialradius",
	  "saturnradiusequatorial","r♄eq"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::saturn_radius_equatorial,"",0,0,0,0,0,0,0},
	{{"radius saturnpolar","saturnpolarradius",
	  "saturnradiuspolar","r♄pol"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::saturn_radius_polar,"",1,0,0,0,0,0,0},
	{{"radius saturnpolar","saturnpolarradius",
	  "saturnradiuspolar","r♄pol"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::saturn_radius_polar,"",0,0,0,0,0,0,0},
	{{"radius uranusequatorial","uranusequatorialradius",
	  "uranusradiusequatorial","r♅eq"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::uranus_radius_equatorial,"",1,0,0,0,0,0,0},
	{{"radius uranusequatorial","uranusequatorialradius",
	  "uranusradiusequatorial","r♅eq"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::uranus_radius_equatorial,"",0,0,0,0,0,0,0},
	{{"radius uranuspolar","uranuspolarradius",
	  "uranusradiuspolar","r♅pol"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::uranus_radius_polar,"",1,0,0,0,0,0,0},
	{{"radius uranuspolar","uranuspolarradius",
	  "uranusradiuspolar","r♅pol"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::uranus_radius_polar,"",0,0,0,0,0,0,0},
	{{"radius neptuneequatorial","neptuneequatorialradius",
	  "neptuneradiusequatorial","r♆eq"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::neptune_radius_equatorial,"",1,0,0,0,0,0,0},
	{{"radius neptuneequatorial","neptuneequatorialradius",
	  "neptuneradiusequatorial","r♆eq"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::neptune_radius_equatorial,"",0,0,0,0,0,0,0},
	{{"radius neptunepolar","neptunepolarradius",
	  "neptuneradiuspolar","r♆pol"},"m",o2scl_const::o2scl_mks,
	 o2scl_mks::neptune_radius_polar,"",1,0,0,0,0,0,0},
	{{"radius neptunepolar","neptunepolarradius",
	  "neptuneradiuspolar","r♆pol"},"cm",o2scl_const::o2scl_cgs,
	 o2scl_cgs::neptune_radius_polar,"",0,0,0,0,0,0,0},
	{{"radius pluto","plutoradius","rpluto","r♇"},"m",
         o2scl_const::o2scl_mks,o2scl_mks::pluto_radius,"",1,0,0,0,0,0,0},
	{{"radius pluto","plutoradius","rpluto","r♇"},"cm",
         o2scl_const::o2scl_cgs,o2scl_cgs::pluto_radius,"",0,0,0,0,0,0,0},
	{{"Rydberg"},"kg*m^2/s^2",o2scl_const::o2scl_mks,
	 o2scl_mks::rydberg,"CODATA 2018",1,0,0,0,0,0,0},
	{{"Rydberg"},"g*cm^2/s^2",o2scl_const::o2scl_cgs,
	 o2scl_cgs::rydberg,"CODATA 2018",0,0,0,0,0,0,0}
  };

}

bool find_constants::unit_match_logic(std::string unit,
                                      const find_constants_list &f) {
  if (unit.length()==0 || f.unit_flag==fc_unknown ||
      (boost::iequals(unit,"mks") &&
       (f.unit_flag==o2scl_const::o2scl_mks ||
        f.unit_flag==fc_none)) ||
      (boost::iequals(unit,"cgs") &&
       (f.unit_flag==o2scl_const::o2scl_cgs ||
        f.unit_flag==fc_none)) ||
      (boost::iequals(unit,"none") && f.unit_flag==fc_none) ||
      boost::iequals(unit,f.unit) ||
      boost::iequals(unit,"any")) {
    return true;
  }
  return false;
}
                                      

int find_constants::find_nothrow(std::string name, std::string unit,
				 vector<find_constants_list> &matches,
				 int verbose) {
  
  o2scl::convert_units<> &cu=o2scl_settings.get_convert_units();
  
  if (verbose>1) {
    std::cout << "find_constants::find_nothrow(): "
	      << "before simplify: " << name << endl;
  }

  // Remove whitespace and punctuation. We need + and - to distinguish
  // between positive and negative particle masses.
  remove_ws_punct(name);
  
  if (verbose>1) {
    std::cout << "find_constants::find_nothrow(): "
	      << "after simplify: " << name << endl;
  }
    
  // Start with a fresh list
  matches.clear();
    
  // Temporarily store matching indexes
  vector<size_t> indexes;

  int match_type=0, match_exact=1, match_pattern=2;
    
  // Initial pass, exact name matches
  for(size_t i=0;i<list.size();i++) {
    for(size_t j=0;j<list[i].names.size();j++) {
      if (verbose>2) {
	std::cout << "find_constants::find_nothrow(): "
		  << name << " " << i << " " << j << " "
		  << list[i].names[j]
		  << " " << boost::iequals(name,list[i].names[j]) << endl;
      }
      string temp=list[i].names[j];
      remove_ws_punct(temp);
      if (boost::iequals(name,temp)) {
        if (verbose>2) {
          cout << "find_constants::find_nothrow(): Found match."
               << endl;
        }
	indexes.push_back(i);
	// Now that we've found a match, don't look in the
	// other names for this list entry
	j=list[i].names.size();
	match_type=match_exact;
      }
    }
  }

#ifdef O2SCL_REGEX    
  std::regex r(name);
#else
  string fn_pat=((string)"*")+name+"*";
#endif
  
  if (verbose>1) {
    std::cout << "find_constants::find_nothrow(): "
	      << "pass 1 indexes: ";
    vector_out(std::cout,indexes,true);
#ifdef O2SCL_REGEX    
    std::cout << "find_constants::find_nothrow(): Using regex "
              << name << std::endl;
#else
    std::cout << "find_constants::find_nothrow(): Using fnmatch() "
              << "with pattern " << fn_pat << std::endl;
#endif
  }
  
  // No matches, so try wildcard matches
  if (indexes.size()==0) {
      
    for(size_t i=0;i<list.size();i++) {
      for(size_t j=0;j<list[i].names.size();j++) {
        string temp=list[i].names[j];
        remove_ws_punct(temp);
#ifdef O2SCL_REGEX
        bool fn_ret=std::regex_search(temp,r);
#else
        bool fn_ret=(fnmatch(fn_pat.c_str(),temp.c_str(),0)==0);
#endif
	if (verbose>2) {
	  std::cout << "find_constants::find_nothrow(): "
		    << name << " " << i << " " << j << " "
		    << list[i].names[j]
		    << " " << fn_ret << endl;
	}
	if (fn_ret==true) {
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
    if (unit_match_logic(unit,list[indexes[0]])) {
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
      
      if (unit_match_logic(unit,list[indexes[i]])) {
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
  cout << "List.size(): " << list.size() << endl;
  os << "name unit flag value units source" << endl;
  os << "  alternate names" << endl;
  for(size_t i=0;i<list.size();i++) {
    os << list[i].names[0] << " ";
    os << list[i].unit << " ";
    if (list[i].unit_flag==o2scl_const::o2scl_mks) {
      os << "MKS ";
    } else if (list[i].unit_flag==o2scl_const::o2scl_cgs) {
      os << "CGS ";
    } else {
      os << "other ";
    }
    os << list[i].val << " ";
    os << "(m:" << list[i].m << ",kg:" << list[i].k
       << ",s:" << list[i].s << ",K:" << list[i].K
       << ",A:" << list[i].A << ",mol:" << list[i].mol
       << ",cd:" << list[i].cd << ")" << endl;
    vector<string> sv;
    rewrap(list[i].source,sv,77);
    for(size_t j=0;j<sv.size();j++) {
      os << "  " << sv[j] << endl;
    }
    if (list[i].names.size()>1) {
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

void find_constants::add_constant(const find_constants_list &f,
                                  int verbose) {

  if (verbose>1) {
    cout << "find_constants::add_constant() attempting to add constant "
         << f.names[0] << " with value " << f.val << endl;
  }
  
  if (f.names.size()==0) {
    O2SCL_ERR2("No names specified in ",
               "find_constants::add_constant().",o2scl::exc_einval);
  }
  
  // Double check that there are no name duplicates before we add
  size_t n_matches=0;
  for(size_t i=0;i<list.size();i++) {
    for(size_t j=0;j<list[i].names.size();j++) {
      for(size_t k=0;k<f.names.size();k++) {
        if (list[i].names[j]==f.names[k]) {
          n_matches++;
        }
      }
    }
  }
  if (n_matches>0) {
    O2SCL_ERR2("Name already found in ",
               "find_constants::add_constant().",o2scl::exc_einval);
  }

  if (verbose>0) {
    cout << "find_constants::add_constant() adding constant "
         << f.names[0] << " with value " << f.val << endl;
  }
  cout << "List.size(): " << list.size() << endl;
  list.push_back(f);
  cout << "List.size(): " << list.size() << endl;
  
  return;
}

void find_constants::del_constant(std::string &name, int verbose) {

  if (verbose>1) {
    cout << "find_constants::add_constant() attempting to remove "
         << "constant named " << name << endl;
  }
  
  size_t n_matches=0, i_match;
  for(size_t i=0;i<list.size();i++) {
    for(size_t j=0;j<list[i].names.size();j++) {
      if (list[i].names[j]==name) {
        n_matches++;
        i_match=i;
      }
    }
  }
  if (n_matches==1) {
    std::vector<find_constants_list>::iterator it=list.begin();
    it+=i_match;
    if (verbose>1) {
      cout << "find_constants::add_constant() Removing "
           << "constant named " << name << " with value "
           << it->val << endl;
    }
    list.erase(it);
    return;
  }
  if (n_matches>1) {
    O2SCL_ERR2("More than one match found in attempting to delete in ",
               "find_constants::del_constant",o2scl::exc_einval);
  }
  O2SCL_ERR2("No matches in ",
             "find_constants::del_constant",o2scl::exc_einval);
  return;
}
