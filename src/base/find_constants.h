/*
  ───────────────────────────────────────────────────────────────────
  
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

  ───────────────────────────────────────────────────────────────────
*/
#ifndef O2SCL_FIND_CONSTANTS_H
#define O2SCL_FIND_CONSTANTS_H
#include <iostream>
#include <vector>

#include <boost/algorithm/string.hpp>

#include <o2scl/constants.h>
#include <o2scl/string_conv.h>

namespace o2scl {

  /** \brief A searchable database of constants with units

      \note This class cannot handle conversions to F and C since
      those are more complicated than a simply multiplication or
      division. 
   */
  template<class fp_t=double> class find_constants {

  public:

    find_constants()  {

      /*
        The code is written in such a way that multiple entries for the
        same constant must be grouped together in this list and have
        exactly the same 'name' array so that the find_nothrow() function
        doesn't report the same constant multiple times.
      */
      list={{{"Schwarzchild radius","rschwarz"},
          "m",o2scl_const::o2scl_mks,
          o2scl_const::schwarzchild_radius_f<fp_t>(),
          "derived from the IAU 2015 nominal solar mass parameter",
          1,0,0,0,0,0,0},
            {{"vacuum permittivity","vacuum electric permittivity",
               "permittivity of free space","epsilon0","ε0"},
             "F/m",o2scl_const::o2scl_mks,
             o2scl_const::vacuum_permittivity_f<fp_t>(),
             "CODATA 2018",-3,-1,4,0,2,0,0},
            {{"vacuum permeability","vacuum electric permeability",
               "permeability of free space","mu0","μ0","magnetic constant"},
             "N/A^2",o2scl_const::o2scl_mks,
             o2scl_const::vacuum_permeability_f<fp_t>(),
             "CODATA 2018",1,1,-2,0,-2,0,0},
            {{"Bohr radius","rbohr"},
             "m",o2scl_const::o2scl_mks,
             o2scl_const::bohr_radius_f<fp_t>(),
             "CODATA 2018",1,0,0,0,0,0,0},
            {{"Thomson cross section","σThomson"},
             "m^2",o2scl_const::o2scl_mks,
             o2scl_const::thomson_csec_f<fp_t>(),
             "CODATA 2018",2,0,0,0,0,0,0},
            {{"classical electron radius","electron radius",
               "relectron","re"},
             "m",o2scl_const::o2scl_mks,2.8179403262e-15,
             "CODATA 2018",1,0,0,0,0,0,0},
            {{"Wien frequency displacement law","b'","bprime","b′"},
             "Hz/K",o2scl_const::o2scl_mks,5.878925757e10,
             "CODATA 2018",0,0,-1,-1,0,0,0},
            {{"Wien wavelength displacement law","b"},
             "m/K",o2scl_const::o2scl_mks,2.897771955e-3,
             "CODATA 2018",1,0,0,-1,0,0,0},
            {{"Planck length"},"m",o2scl_const::o2scl_mks,
             sqrt(o2scl_const::gravitational_constant_f<fp_t>()*
                  o2scl_const::hbar_f<fp_t>()/
                  o2scl_const::speed_of_light_f<fp_t>()/
                  o2scl_const::speed_of_light_f<fp_t>()/
                  o2scl_const::speed_of_light_f<fp_t>()),
             "derived",1,0,0,0,0,0,0},
            {{"Planck mass"},"kg",o2scl_const::o2scl_mks,
             sqrt(o2scl_const::hbar_f<fp_t>()/
                  o2scl_const::gravitational_constant_f<fp_t>()*
                  o2scl_const::speed_of_light_f<fp_t>()),
             "derived",0,1,0,0,0,0,0},
            {{"Planck time"},"s",o2scl_const::o2scl_mks,
             sqrt(o2scl_const::gravitational_constant_f<fp_t>()*
                  o2scl_const::hbar_f<fp_t>()/
                  o2scl_const::speed_of_light_f<fp_t>()/
                  o2scl_const::speed_of_light_f<fp_t>()/
                  o2scl_const::speed_of_light_f<fp_t>()/
                  o2scl_const::speed_of_light_f<fp_t>()/
                  o2scl_const::speed_of_light_f<fp_t>()),
             "derived",0,0,1,0,0,0,0},
            {{"Planck temperature"},"K",o2scl_const::o2scl_mks,
             sqrt(o2scl_const::hbar_f<fp_t>()*
                  o2scl_const::speed_of_light_f<fp_t>()*
                  o2scl_const::speed_of_light_f<fp_t>()*
                  o2scl_const::speed_of_light_f<fp_t>()*
                  o2scl_const::speed_of_light_f<fp_t>()*
                  o2scl_const::speed_of_light_f<fp_t>()/
                  o2scl_const::gravitational_constant_f<fp_t>()/
                  o2scl_const::boltzmann_f<fp_t>()/
                  o2scl_const::boltzmann_f<fp_t>()),
             "derived",0,0,0,1,0,0,0},
            /*
              Things to add in the future:
              Astrophysical/cosmological constants from PDG
            */
            {{"elementary charge","e"},"C",
             o2scl_const::o2scl_mks,o2scl_const::elem_charge_f<fp_t>(),
             "exact",0,0,1,0,1,0,0},
            {{"Bohr magneton"},"J/T",
             o2scl_const::o2scl_mks,o2scl_const::bohr_magneton_f<fp_t>(),
             "CODATA 2018",1,1,0,0,1,0,0},
            {{"nuclear magneton"},"J/T",
             o2scl_const::o2scl_mks,o2scl_const::nuclear_magneton_f<fp_t>(),
             "CODATA 2018",1,1,0,0,1,0,0},
            {{"strong coupling constant at the Z mass"},"",
             fc_none,0.1179,
             ((std::string)"https://pdg.lbl.gov/2021/reviews/")+
             "contents_sports.html",
             0,0,0,0,0,0,0},
            {{"Schwarzchild radius","rschwarz"},
             "cm",o2scl_const::o2scl_cgs,
             o2scl_const::schwarzchild_radius_f<fp_t>(o2scl_const::o2scl_cgs),
             "derived from the IAU 2015 nominal solar mass parameter",
             0,0,0,0,0,0,0},
            {{"Schwarzchild radius","rschwarz"},
             "km",o2scl_const::o2scl_mks,
             o2scl_const::schwarzchild_radius_f<fp_t>()/1.0e3,
             "derived from the IAU 2015 nominal solar mass parameter",
             1,0,0,0,0,0,0},
            {{"speed of light","c","lightspeed"},
             "m/s",o2scl_const::o2scl_mks,
             o2scl_const::speed_of_light_f<fp_t>(o2scl_const::o2scl_mks),
             "exact",1,0,-1,0,0,0,0},
            {{"speed of light","c","lightspeed"},
             "cm/s",o2scl_const::o2scl_cgs,
             o2scl_const::speed_of_light_f<fp_t>(o2scl_const::o2scl_cgs),
             "exact",0,0,0,0,0,0,0},
            {{"gravitational","g","gnewton"},"m^3/kg/s^2",
             o2scl_const::o2scl_mks,
             o2scl_const::gravitational_constant_f<fp_t>(),"CODATA 2018",
             3,-1,-2,0,0,0,0},
            {{"gravitational","g","gnewton"},"cm^3/g/s^2",
             o2scl_const::o2scl_cgs,
             o2scl_const::gravitational_constant_f<fp_t>
             (o2scl_const::o2scl_cgs),"CODATA 2018",0,0,0,0,0,0,0},
            {{"Boltzmann's","kb","boltzmann"},
             "m^2/kg/s^2/K",o2scl_const::o2scl_mks,
             o2scl_const::boltzmann_f<fp_t>(),
             "exact",2,-1,-2,-1,0,0,0},
            {{"Boltzmann's","kb","boltzmann"},
             "cm^2/g/s^2/K",o2scl_const::o2scl_cgs,
             o2scl_const::boltzmann_f<fp_t>(o2scl_const::o2scl_cgs),
             "exact",0,0,0,0,0,0,0},
            {{"Stefan-Boltzmann","sigmasb","stefanboltzmann","ssb","σsb"},
             "kg/s^3/K^4",o2scl_const::o2scl_mks,
             o2scl_const::stefan_boltz_cons_f<fp_t>(),
             "exact; derived from k_B, c, and ħ",0,1,-3,-4,0,0,0},
            {{"Stefan-Boltzmann","sigmasb","stefanboltzmann",
               "ssb","σsb"},
             "g/s^3/K^4",o2scl_const::o2scl_cgs,
             o2scl_const::stefan_boltz_cons_f<fp_t>(o2scl_const::o2scl_cgs),
             "exact; derived from k_B, c, and ħ",0,0,0,0,0,0,0},
            {{"Planck","h","plancks"},
             "kg*m^2/s",o2scl_const::o2scl_mks,
             o2scl_const::planck_f<fp_t>(o2scl_const::o2scl_mks),
             "exact",2,1,-1,0,0,0,0},
            {{"Planck","h","plancks"},
             "g*cm^2/s",o2scl_const::o2scl_cgs,
             o2scl_const::planck_f<fp_t>(o2scl_const::o2scl_cgs),
             "exact",0,0,0,0,0,0,0},
            {{"reduced Planck","hbar","ħ","reducedplancks"},
             "kg*m^2/s",o2scl_const::o2scl_mks,
             o2scl_const::hbar_f<fp_t>(o2scl_const::o2scl_mks),
             "exact; derived from the Planck constant",2,1,-1,0,0,0,0},
            {{"reduced Planck","hbar","ħ","reducedplancks"},
             "g*cm^2/s",o2scl_const::o2scl_cgs,
             o2scl_const::hbar_f<fp_t>(o2scl_const::o2scl_cgs),
             "exact; derived from the Planck constant",0,0,0,0,0,0,0},
            {{"Avogadro's number","na","avogadro"},
             "",fc_none,o2scl_const::avogadro_f<fp_t>(),
             "exact",0,0,0,0,0,0,0},
            {{"fine structure","alphaem","alpha","αem"},"",fc_none,
             o2scl_const::fine_structure_f<fp_t>(),
             "CODATA 2018",0,0,0,0,0,0,0},
            {{"pi","π"},"",fc_none,boost::math::constants::pi<fp_t>(),
             "exact",0,0,0,0,0,0,0},
            {{"zeta2","zeta(2)","ζ(2)"},"",fc_none,
             boost::math::constants::zeta_two<fp_t>(),
             "exact",0,0,0,0,0,0,0},
            {{"zeta3","zeta(3)","ζ(3)"},"",fc_none,
             boost::math::constants::zeta_three<fp_t>(),
             "exact",0,0,0,0,0,0,0},
            {{"pi2","pisquared","π²"},"",fc_none,
             boost::math::constants::pi_sqr<fp_t>(),
             "exact",0,0,0,0,0,0,0},
            {{"pi3","picubed","π³"},"",fc_none,
             boost::math::constants::pi<fp_t>()*
             boost::math::constants::pi_sqr<fp_t>(),
             "exact",0,0,0,0,0,0,0},
            {{"pi4","pifourth","π⁴"},"",fc_none,
             boost::math::constants::pi_sqr<fp_t>()*
             boost::math::constants::pi_sqr<fp_t>(),
             "exact",0,0,0,0,0,0,0},
            {{"rootpi","squarerootpi","√π"},"",fc_none,
             boost::math::constants::root_pi<fp_t>(),
             "exact",0,0,0,0,0,0,0},
            {{"Euler-Mascheroni","euler"},"",fc_none,
             boost::math::constants::euler<fp_t>(),
             "exact",0,0,0,0,0,0,0},
            {{"sin2thetaw","sin2θW","sin²θW"},"",fc_none,
             o2scl_const::sin2_theta_weak_f<fp_t>(),
             "PDG 2020 value",0,0,0,0,0,0,0},
            {{"gfermi","gf"},"s^4/m^4/kg^2",o2scl_const::o2scl_mks,
             o2scl_const::gfermi_f<fp_t>(),
             ((std::string)"derived from CODATA 2018 value for ")+
             "G_Fermi (identical to "+
             "PDG 2020 value) and CODATA 2018 value of electron volt",
             -4,-2,4,0,0,0,0},
            {{"gfermi","gf"},"s^4/cm^4/g^2",o2scl_const::o2scl_cgs,
             o2scl_const::gfermi_f<fp_t>(o2scl_const::o2scl_cgs),
             ((std::string)"derived from CODATA 2018 value for ")+
             "G_Fermi (identical to "+
             "PDG 2020 value) and CODATA 2018 value of electron volt",
             0,0,0,0,0,0,0},
            {{"gfermi","gf"},"1/GeV^2",0,
             o2scl_const::gfermi_gev2_f<fp_t>(),
             "CODATA 2018 (identical to PDG 2020 value)",0,0,0,0,0,0,0},
            {{"elementarycharge","electroncharge","e","chargeelectron",
               "qelectron"},"C",
               o2scl_const::o2scl_mks,o2scl_const::elem_charge_f<fp_t>(),
               "exact",0,0,1,0,1,0,0},
            {{"hbarc","ħc"},"MeV*fm",
             fc_other,o2scl_const::hc_mev_fm_f<fp_t>(),
             "derived from Plack constant",0,0,0,0,0,0,0},
            {{"hbarc","ħc"},"J*m",o2scl_const::o2scl_mks,
             o2scl_const::hbarc_f<fp_t>(o2scl_const::o2scl_mks),
             "derived from Plack constant",3,1,-2,0,0,0,0},
            {{"hbarc","ħc"},"erg*cm",o2scl_const::o2scl_cgs,
             o2scl_const::hbarc_f<fp_t>(o2scl_const::o2scl_cgs),
             "derived from Plack constant",0,0,0,0,0,0,0},
            {{"mass W","Wmass","mW","mW"},"GeV",
             fc_other,80.379,
             ((std::string)"https://pdg.lbl.gov/2021/tables/")+
             "contents_tables.html on 10/27/21",
             0,1,0,0,0,0,0},
            {{"mass Z","Zmass","mZ","mZ"},"GeV",
             fc_other,91.1876,
             ((std::string)"https://pdg.lbl.gov/2021/tables/contents")+
             "_tables.html on 10/27/21",
             0,1,0,0,0,0,0},
            {{"mass H","Hmass","mH","mH","mass higgs","higgs mass",
               "mH0","mH⁰"},"GeV",fc_other,125.25,
               ((std::string)"https://pdg.lbl.gov/2021/tables/contents")+
               "_tables.html on 10/27/21",
               0,1,0,0,0,0,0},
            {{"mass electron","electronmass","melectron","melec"},"kg",
             o2scl_const::o2scl_mks,o2scl_const::mass_electron_f<fp_t>(),
             "CODATA 2018",0,1,0,0,0,0,0},
            {{"mass electron","electronmass","melectron","melec"},"g",
             o2scl_const::o2scl_cgs,o2scl_const::mass_electron_f<fp_t>
             (o2scl_const::o2scl_cgs),"CODATA 2018",0,0,0,0,0,0,0},
            {{"mass muon","muonmass","mmuon"},"kg",o2scl_const::o2scl_mks,
             o2scl_const::mass_muon_f<fp_t>(),"CODATA 2018",0,1,0,0,0,0,0},
            {{"mass muon","muonmass","mmuon"},"g",o2scl_const::o2scl_cgs,
             o2scl_const::mass_muon_f<fp_t>(o2scl_const::o2scl_cgs),
             "CODATA 2018",0,0,0,0,0,0,0},
            {{"mass tau","taumass","mtau"},"kg",o2scl_const::o2scl_mks,
             o2scl_const::mass_tau_f<fp_t>(),"CODATA 2018",0,1,0,0,0,0,0},
            {{"mass tau","taumass","mtau"},"g",o2scl_const::o2scl_cgs,
             o2scl_const::mass_tau_f<fp_t>(o2scl_const::o2scl_cgs),
             "CODATA 2018",0,0,0,0,0,0,0},
            {{"mass neutron","neutronmass","mneutron","mneut"},"kg",
             o2scl_const::o2scl_mks,o2scl_const::mass_neutron_f<fp_t>(),
             "CODATA 2018",0,1,0,0,0,0,0},
            {{"mass neutron","neutronmass","mneutron","mneut"},"g",
             o2scl_const::o2scl_cgs,o2scl_const::mass_neutron_f<fp_t>
             (o2scl_const::o2scl_cgs),"CODATA 2018",0,0,0,0,0,0,0},
            {{"mass proton","protonmass","mproton","mprot"},"kg",
             o2scl_const::o2scl_mks,o2scl_const::mass_proton_f<fp_t>(),
             "CODATA 2018",0,1,0,0,0,0,0},
            {{"mass proton","protonmass","mproton","mprot"},"g",
             o2scl_const::o2scl_cgs,o2scl_const::mass_proton_f<fp_t>
             (o2scl_const::o2scl_cgs),"CODATA 2018",0,0,0,0,0,0,0},
            {{"mass deuteron","deuteronmass","mdeuteron","mdeut"},"kg",
             o2scl_const::o2scl_mks,o2scl_const::mass_deuteron_f<fp_t>(),
             "CODATA 2018",0,1,0,0,0,0,0},
            {{"mass deuteron","deuteronmass","mdeuteron","mdeut"},"g",
             o2scl_const::o2scl_cgs,o2scl_const::mass_deuteron_f<fp_t>
             (o2scl_const::o2scl_cgs),"CODATA 2018",0,0,0,0,0,0,0},
            {{"mass triton","tritonmass","mtriton"},"kg",
             o2scl_const::o2scl_mks,
             o2scl_const::mass_triton_f<fp_t>(),
             "CODATA 2018",0,1,0,0,0,0,0},
            {{"mass triton","tritonmass","mtriton"},"g",
             o2scl_const::o2scl_cgs,
             o2scl_const::mass_triton_f<fp_t>(o2scl_const::o2scl_cgs),
             "CODATA 2018",0,0,0,0,0,0,0},
            {{"mass helion","helionmass","mhelion"},"kg",
             o2scl_const::o2scl_mks,
             o2scl_const::mass_helion_f<fp_t>(),
             "CODATA 2018",0,1,0,0,0,0,0},
            {{"mass helion","helionmass","mhelion"},"g",
             o2scl_const::o2scl_cgs,
             o2scl_const::mass_helion_f<fp_t>(o2scl_const::o2scl_cgs),
             "CODATA 2018",0,0,0,0,0,0,0},
            {{"mass alpha","alphamass","malpha","mα"},"kg",
             o2scl_const::o2scl_mks,
             o2scl_const::mass_alpha_f<fp_t>(),
             "CODATA 2018",0,1,0,0,0,0,0},
            {{"mass alpha","alphamass","malpha","mα"},"g",
             o2scl_const::o2scl_cgs,
             o2scl_const::mass_alpha_f<fp_t>(o2scl_const::o2scl_cgs),
             "CODATA 2018",0,0,0,0,0,0,0},
            {{"mass lambda","lambdamass","mlambda","mΛ"},"MeV",0,
             o2scl_const::mass_lambda_MeV_f<fp_t>(),
             "\"OUR FIT\" value from PDG 2020",0,0,0,0,0,0,0},
            {{"mass sigma minus","sigmaminusmass",
               "msigma-","mΣ-","mΣ⁻"},"MeV",0,
             o2scl_const::mass_sigma_minus_MeV_f<fp_t>(),
             "\"OUR FIT\" value from PDG 2020",0,0,0,0,0,0,0},
            {{"mass sigma zero","sigmazeromass",
               "msigma0","mΣ0","mΣ⁰"},"MeV",0,
             o2scl_const::mass_sigma_zero_MeV_f<fp_t>(),
             "\"OUR FIT\" value from PDG 2020",0,0,0,0,0,0,0},
            {{"mass sigma plus","sigmaplusmass",
               "msigma+","mΣ+","mΣ⁺"},"MeV",0,
             o2scl_const::mass_sigma_plus_MeV_f<fp_t>(),
             "\"OUR FIT\" value from PDG 2020",0,0,0,0,0,0,0},
            {{"mass cascade zero","cascadezeromass","mcascade0","mxi0",
               "mΞ0","mΞ⁰"},
             "MeV",0,o2scl_const::mass_cascade_zero_MeV_f<fp_t>(),
             "\"OUR FIT\" value from PDG 2020",0,0,0,0,0,0,0},
            {{"mass cascade minus","cascademinusmass","mcascade-","mxi-",
               "mΞ-","mΞ⁻"},
             "MeV",0,o2scl_const::mass_cascade_minus_MeV_f<fp_t>(),
             "\"OUR FIT\" value from PDG 2020",0,0,0,0,0,0,0},
            {{"mass up","upmass","mup"},"MeV",0,
             o2scl_const::mass_up_MeV_f<fp_t>(),
             "\"OUR EVALUATION\" value from PDG 2020",0,0,0,0,0,0,0},
            {{"mass down","downmass","mdown"},"MeV",0,
             o2scl_const::mass_down_MeV_f<fp_t>(),
             "\"OUR EVALUATION\" value from PDG 2020",0,0,0,0,0,0,0},
            {{"mass strange","strangemass","mstrange"},"MeV",0,
             o2scl_const::mass_strange_MeV_f<fp_t>(),
             "\"OUR EVALUATION\" value from PDG 2020",0,0,0,0,0,0,0},
            {{"mass solar","solarmass","masssun",
               "sunmass","msun","modot","m☉"},
             "kg",o2scl_const::o2scl_mks,o2scl_const::solar_mass_f<fp_t>
             (o2scl_const::o2scl_mks),
             ((std::string)"derived from IAU's 2015 ")+
             "nominal value of the solar "+
             +"mass parameter divided by the CODATA 2018 value of the "+
             "gravitational constant",0,1,0,0,0,0,0},
            {{"mass solar","solarmass","masssun",
               "sunmass","msun","modot","m☉"},
             "g",o2scl_const::o2scl_cgs,o2scl_const::solar_mass_f<fp_t>
             (o2scl_const::o2scl_cgs),
             ((std::string)"derived from IAU's 2015 nominal ")+
             "value of the solar "+
             +"mass parameter divided by the CODATA 2018 value of the "+
             "gravitational constant",0,0,0,0,0,0,0},
            {{"mass mercury","mercurymass","mmercury","m☿"},
             "kg",o2scl_const::o2scl_mks,
             o2scl_const::mercury_mass_f<fp_t>(o2scl_const::o2scl_mks),
             "",0,1,0,0,0,0,0},
            {{"mass mercury","mercurymass","mmercury","m☿"},
             "g",o2scl_const::o2scl_cgs,
             o2scl_const::mercury_mass_f<fp_t>(o2scl_const::o2scl_cgs),
             "",0,0,0,0,0,0,0},
            {{"mass venus","venusmass","mvenus","m♀"},"kg",
             o2scl_const::o2scl_mks,
             o2scl_const::venus_mass_f<fp_t>(o2scl_const::o2scl_mks),
             "",0,1,0,0,0,0,0},
            {{"mass venus","venusmass","mvenus","m♀"},"g",
             o2scl_const::o2scl_cgs,
             o2scl_const::venus_mass_f<fp_t>(o2scl_const::o2scl_cgs),
             "",0,0,0,0,0,0,0},
            {{"mass earth","earthmass","mearth","m♁","m⊕","moplus"},
             "kg",o2scl_const::o2scl_mks,
             o2scl_const::earth_mass_f<fp_t>(o2scl_const::o2scl_mks),
             "IAU 2015 nominal value",0,1,0,0,0,0,0},
            {{"mass earth","earthmass","mearth","m♁","m⊕","moplus"},
             "g",o2scl_const::o2scl_cgs,
             o2scl_const::earth_mass_f<fp_t>(o2scl_const::o2scl_cgs),
             "IAU 2015 nominal value",0,0,0,0,0,0,0},
            {{"mass mars","marsmass","mmars","m♂"},"kg",
             o2scl_const::o2scl_mks,
             o2scl_const::mars_mass_f<fp_t>(o2scl_const::o2scl_mks),
             "",0,1,0,0,0,0,0},
            {{"mass mars","marsmass","mmars","m♂"},"g",
             o2scl_const::o2scl_cgs,
             o2scl_const::mars_mass_f<fp_t>(o2scl_const::o2scl_cgs),
             "",0,0,0,0,0,0,0},
            {{"mass jupiter","jupitermass","mjupiter","mjup","m♃"},
             "kg",o2scl_const::o2scl_mks,
             o2scl_const::jupiter_mass_f<fp_t>(o2scl_const::o2scl_mks),
             "IAU 2015 nominal value",0,1,0,0,0,0,0},
            {{"mass jupiter","jupitermass","mjupiter","mjup","m♃"},
             "g",o2scl_const::o2scl_cgs,
             o2scl_const::jupiter_mass_f<fp_t>(o2scl_const::o2scl_cgs),
             "IAU 2015 nominal value",0,0,0,0,0,0,0},
            {{"mass saturn","saturnmass","msaturn","m♄"},"kg",
             o2scl_const::o2scl_mks,
             o2scl_const::saturn_mass_f<fp_t>(o2scl_const::o2scl_mks),
             "",0,1,0,0,0,0,0},
            {{"mass saturn","saturnmass","msaturn","m♄"},"g",
             o2scl_const::o2scl_cgs,
             o2scl_const::saturn_mass_f<fp_t>(o2scl_const::o2scl_cgs),
             "",0,0,0,0,0,0,0},
            {{"mass uranus","uranusmass","muranus","m♅"},"kg",
             o2scl_const::o2scl_mks,
             o2scl_const::uranus_mass_f<fp_t>(o2scl_const::o2scl_mks),
             "",0,1,0,0,0,0,0},
            {{"mass uranus","uranusmass","muranus","m♅"},"g",
             o2scl_const::o2scl_cgs,
             o2scl_const::uranus_mass_f<fp_t>(o2scl_const::o2scl_cgs),
             "",0,0,0,0,0,0,0},
            {{"mass neptune","neptunemass","mneptune","m♆"},"kg",
             o2scl_const::o2scl_mks,
             o2scl_const::neptune_mass_f<fp_t>(o2scl_const::o2scl_mks),
             "",0,1,0,0,0,0,0},
            {{"mass neptune","neptunemass","mneptune","m♆"},"g",
             o2scl_const::o2scl_cgs,
             o2scl_const::neptune_mass_f<fp_t>(o2scl_const::o2scl_cgs),
             "",0,0,0,0,0,0,0},
            {{"mass pluto","plutomass","mpluto","m♇"},"kg",
             o2scl_const::o2scl_mks,
             o2scl_const::pluto_mass_f<fp_t>(o2scl_const::o2scl_mks),
             "",0,1,0,0,0,0,0},
            {{"mass pluto","plutomass","mpluto","m♇"},"g",
             o2scl_const::o2scl_cgs,
             o2scl_const::pluto_mass_f<fp_t>(o2scl_const::o2scl_cgs),
             "",0,0,0,0,0,0,0},
            {{"radius solar","solarradius","radiussun","sunradius",
               "rsun","r☉"},"m",o2scl_const::o2scl_mks,
               o2scl_const::solar_radius_f<fp_t>(o2scl_const::o2scl_mks),
               "",1,0,0,0,0,0,0},
            {{"radius solar","solarradius","radiussun","sunradius",
               "rsun","r☉"},"cm",o2scl_const::o2scl_cgs,
               o2scl_const::solar_radius_f<fp_t>(o2scl_const::o2scl_cgs),
               "",0,0,0,0,0,0,0},
            {{"radius mercury","mercuryradius","rmercury","r☿"},"m",
             o2scl_const::o2scl_mks,
             o2scl_const::mercury_radius_f<fp_t>(o2scl_const::o2scl_mks),"",
             1,0,0,0,0,0,0},
            {{"radius mercury","mercuryradius","rmercury","r☿"},"cm",
             o2scl_const::o2scl_cgs,
             o2scl_const::mercury_radius_f<fp_t>(o2scl_const::o2scl_cgs),"",
             0,0,0,0,0,0,0},
            {{"radius venus","venusradius","rvenus","r♀"},"m",
             o2scl_const::o2scl_mks,
             o2scl_const::venus_radius_f<fp_t>(o2scl_const::o2scl_mks),"",
             1,0,0,0,0,0,0},
            {{"radius venus","venusradius","rvenus","r♀"},"cm",
             o2scl_const::o2scl_cgs,
             o2scl_const::venus_radius_f<fp_t>(o2scl_const::o2scl_cgs),"",
             0,0,0,0,0,0,0},
            {{"radius earth equatorial","earthequatorialradius",
               "earthradiusequatorial","r♁eq","r⊕eq"},"m",
               o2scl_const::o2scl_mks,
               o2scl_const::earth_radius_eq_f<fp_t>(o2scl_const::o2scl_mks),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               1,0,0,0,0,0,0},
            {{"radius earth equatorial","earthequatorialradius",
               "earthradiusequatorial","r♁eq","r⊕eq"},"cm",
               o2scl_const::o2scl_cgs,
               o2scl_const::earth_radius_eq_f<fp_t>(o2scl_const::o2scl_cgs),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               0,0,0,0,0,0,0},
            {{"radius earth polar","earthpolarradius",
               "earthradiuspolar","r♁pol","r⊕pol"},"m",
               o2scl_const::o2scl_mks,
               o2scl_const::earth_radius_pol_f<fp_t>(o2scl_const::o2scl_mks),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               1,0,0,0,0,0,0},
            {{"radius earth polar","earthpolarradius",
               "earthradiuspolar","r♁pol","r⊕pol"},"cm",
               o2scl_const::o2scl_cgs,
               o2scl_const::earth_radius_pol_f<fp_t>(o2scl_const::o2scl_cgs),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               0,0,0,0,0,0,0},
            {{"radius mars equatorial","marsequatorialradius",
               "marsradiusequatorial","r♂eq"},"m",o2scl_const::o2scl_mks,
               o2scl_const::mars_radius_eq_f<fp_t>(o2scl_const::o2scl_mks),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               1,0,0,0,0,0,0},
            {{"radius mars equatorial","marsequatorialradius",
               "marsradiusequatorial","r♂eq"},"cm",o2scl_const::o2scl_cgs,
               o2scl_const::mars_radius_eq_f<fp_t>(o2scl_const::o2scl_cgs),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               0,0,0,0,0,0,0},
            {{"radius mars polar","marspolarradius",
               "marsradiuspolar","r♂pol"},"m",o2scl_const::o2scl_mks,
               o2scl_const::mars_radius_pol_f<fp_t>(o2scl_const::o2scl_mks),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               1,0,0,0,0,0,0},
            {{"radius mars polar","marspolarradius",
               "marsradiuspolar","r♂pol"},"cm",o2scl_const::o2scl_cgs,
               o2scl_const::mars_radius_pol_f<fp_t>(o2scl_const::o2scl_cgs),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               0,0,0,0,0,0,0},
            {{"radius jupiter equatorial","jupiterequatorialradius",
               "jupiterradiusequatorial","r♃eq"},"m",
               o2scl_const::o2scl_mks,
               o2scl_const::jupiter_radius_eq_f<fp_t>(o2scl_const::o2scl_mks),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               1,0,0,0,0,0,0},
            {{"radius jupiter equatorial","jupiterequatorialradius",
               "jupiterradiusequatorial","r♃eq"},"cm",
               o2scl_const::o2scl_cgs,
               o2scl_const::jupiter_radius_eq_f<fp_t>(o2scl_const::o2scl_cgs),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               0,0,0,0,0,0,0},
            {{"radius jupiter polar","jupiterpolarradius",
               "jupiterradiuspolar","r♃pol"},"m",o2scl_const::o2scl_mks,
               o2scl_const::jupiter_radius_pol_f<fp_t>(o2scl_const::o2scl_mks),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               1,0,0,0,0,0,0},
            {{"radius jupiter polar","jupiterpolarradius",
               "jupiterradiuspolar","r♃pol"},"cm",o2scl_const::o2scl_cgs,
               o2scl_const::jupiter_radius_pol_f<fp_t>(o2scl_const::o2scl_cgs),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               0,0,0,0,0,0,0},
            {{"radius saturn equatorial","saturnequatorialradius",
               "saturnradiusequatorial","r♄eq"},"m",
               o2scl_const::o2scl_mks,
               o2scl_const::saturn_radius_eq_f<fp_t>(o2scl_const::o2scl_mks),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               1,0,0,0,0,0,0},
            {{"radius saturn equatorial","saturnequatorialradius",
               "saturnradiusequatorial","r♄eq"},"cm",o2scl_const::o2scl_cgs,
               o2scl_const::saturn_radius_eq_f<fp_t>(o2scl_const::o2scl_cgs),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               0,0,0,0,0,0,0},
            {{"radius saturn polar","saturnpolarradius",
               "saturnradiuspolar","r♄pol"},"m",o2scl_const::o2scl_mks,
               o2scl_const::saturn_radius_pol_f<fp_t>(o2scl_const::o2scl_mks),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               1,0,0,0,0,0,0},
            {{"radius saturn polar","saturnpolarradius",
               "saturnradiuspolar","r♄pol"},"cm",o2scl_const::o2scl_cgs,
               o2scl_const::saturn_radius_pol_f<fp_t>(o2scl_const::o2scl_cgs),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               0,0,0,0,0,0,0},
            {{"radius uranus equatorial","uranusequatorialradius",
               "uranusradiusequatorial","r♅eq"},"m",o2scl_const::o2scl_mks,
               o2scl_const::uranus_radius_eq_f<fp_t>(o2scl_const::o2scl_mks),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               1,0,0,0,0,0,0},
            {{"radius uranus equatorial","uranusequatorialradius",
               "uranusradiusequatorial","r♅eq"},"cm",o2scl_const::o2scl_cgs,
               o2scl_const::uranus_radius_eq_f<fp_t>(o2scl_const::o2scl_cgs),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               0,0,0,0,0,0,0},
            {{"radius uranus polar","uranuspolarradius",
               "uranusradiuspolar","r♅pol"},"m",o2scl_const::o2scl_mks,
               o2scl_const::uranus_radius_pol_f<fp_t>(o2scl_const::o2scl_mks),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               1,0,0,0,0,0,0},
            {{"radius uranus polar","uranuspolarradius",
               "uranusradiuspolar","r♅pol"},"cm",o2scl_const::o2scl_cgs,
               o2scl_const::uranus_radius_pol_f<fp_t>(o2scl_const::o2scl_cgs),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               0,0,0,0,0,0,0},
            {{"radius neptune equatorial","neptuneequatorialradius",
               "neptuneradiusequatorial","r♆eq"},"m",o2scl_const::o2scl_mks,
               o2scl_const::neptune_radius_eq_f<fp_t>(o2scl_const::o2scl_mks),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               1,0,0,0,0,0,0},
            {{"radius neptune equatorial","neptuneequatorialradius",
               "neptuneradiusequatorial","r♆eq"},"cm",o2scl_const::o2scl_cgs,
               o2scl_const::neptune_radius_eq_f<fp_t>(o2scl_const::o2scl_cgs),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               0,0,0,0,0,0,0},
            {{"radius neptune polar","neptunepolarradius",
               "neptuneradiuspolar","r♆pol"},"m",o2scl_const::o2scl_mks,
               o2scl_const::neptune_radius_pol_f<fp_t>(o2scl_const::o2scl_mks),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               1,0,0,0,0,0,0},
            {{"radius neptune polar","neptunepolarradius",
               "neptuneradiuspolar","r♆pol"},"cm",o2scl_const::o2scl_cgs,
               o2scl_const::neptune_radius_pol_f<fp_t>(o2scl_const::o2scl_cgs),
               "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
               0,0,0,0,0,0,0},
            {{"radius pluto","plutoradius","rpluto","r♇"},"m",
             o2scl_const::o2scl_mks,
             o2scl_const::pluto_radius_f<fp_t>(o2scl_const::o2scl_mks),
             "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
             1,0,0,0,0,0,0},
            {{"radius pluto","plutoradius","rpluto","r♇"},"cm",
             o2scl_const::o2scl_cgs,
             o2scl_const::pluto_radius_f<fp_t>(o2scl_const::o2scl_cgs),
             "https://nssdc.gsfc.nasa.gov/planetary/factsheet/",
             0,0,0,0,0,0,0},
            {{"Rydberg"},"kg*m^2/s^2",o2scl_const::o2scl_mks,
             o2scl_const::rydberg_f<fp_t>(o2scl_const::o2scl_mks),
             "CODATA 2018",2,1,-2,0,0,0,0},
            {{"Rydberg"},"g*cm^2/s^2",o2scl_const::o2scl_cgs,
             o2scl_const::rydberg_f<fp_t>(o2scl_const::o2scl_cgs),
             "CODATA 2018",0,0,0,0,0,0,0},
            {{"tropical year","yeartropical"},"s",o2scl_const::o2scl_mks,
             o2scl_const::tropical_year_f<fp_t>(),
             ((std::string)"PDG 2021 (https://pdg.lbl.gov/2021/")+
             "reviews/contents_sports.html)",
             0,0,1,0,0,0,0},
            {{"sidereal year","yearsidereal"},"s",o2scl_const::o2scl_mks,
             o2scl_const::sidereal_year_f<fp_t>(),
             ((std::string)"PDG 2021 (https://pdg.lbl.gov/2021/")+
             "reviews/contents_sports.html)",
             0,0,1,0,0,0,0}
      };
      
    }

    /// Type for constant database (also used for list of matches)
    class const_entry {
      
    public:
      
      /// List of names for the constant, with the preferred name first
      std::vector<std::string> names;
      /// Unit
      std::string unit;
      /// Flag (currently in the range 0 to 4)
      int unit_flag;
      /// Value
      fp_t val;
      /// Source or reference for value
      std::string source;
      /// Power of length
      int m;
      /// Power of mass
      int k;
      /// Power of time
      int s;
      /// Power of temperature
      int K;
      /// Power of current
      int A;
      /// Power of moles
      int mol;
      /// Power of luminous intensity
      int cd;
      
    };

    /// \name Return values for find_nothrow()
    //@{
    static const int one_exact_match_unit_match=0;
    static const int one_exact_match_unit_diff=1;
    static const int exact_matches_no_unit=2;
    static const int exact_matches_unit_match=3;
    static const int exact_matches_unit_diff=4;
    static const int one_pattern_match_unit_match=5;
    static const int one_pattern_match_unit_diff=6;
    static const int pattern_matches_no_unit=7;
    static const int pattern_matches_unit_match=8;
    static const int pattern_matches_unit_diff=9;
    static const int no_matches=10;
    //@}

    /// \name List of constants and unit match function [protected]
    //@{
    /// Database of constant values
    std::vector<const_entry> list;

    /** \brief The function which decides if the requested unit matches
        the specified list entry

        Units match if 
        - the unit is unspecified (string of length zero) or 
        "none" and the flag is \ref fc_none
        - the unit is equal to "any" (case-insensitive comparison)
        - the unit is equal to the list unit (case-insensitive comparison)
        - the unit is "mks" (case-insensitive comparison) and the unit
        flag is either o2scl_mks or fc_none
        - the unit is "cgs" (case-insensitive comparison) and the unit
        flag is either o2scl_cgs or fc_none
    */
    bool unit_match_logic(std::string unit,
                          const const_entry &f) const {
      if (boost::iequals(unit,"any") ||
          ((boost::iequals(unit,"none") || unit.length()==0) &&
           f.unit_flag==fc_none) ||
          (boost::iequals(unit,"mks") &&
           (f.unit_flag==o2scl_const::o2scl_mks ||
            f.unit_flag==fc_none)) ||
          (boost::iequals(unit,"cgs") &&
           (f.unit_flag==o2scl_const::o2scl_cgs ||
            f.unit_flag==fc_none)) ||
          boost::iequals(unit,f.unit)) {
        return true;
      }
      return false;
    }
    //@}
    
    // FYI, from constants.h, we have:
    //
    // static const size_t o2scl_mks=1;
    // static const size_t o2scl_cgs=2;
    
    /// \name Other possible values of the unit flag
    //@{
    static const int fc_unknown=0;
    static const int fc_none=3;
    static const int fc_other=4;
    //@}

    /// \name Functions to show or modify the constant list
    //@{
    /** \brief Output the full list of constants to \c os 
     */
    void output_list(std::ostream &os) const  {
      for(size_t i=0;i<list.size();i++) {
        std::string s;
        s+=list[i].names[0];
        s+=" ";
        s+=o2scl::dtos(list[i].val);
        s+=" ";
        s+=list[i].unit;
        s+=" ";
        if (list[i].names.size()>1) {
          for(size_t j=1;j<list[i].names.size();j++) {
            s+='\''+list[i].names[j]+"\' ";
          }
        }
        std::vector<std::string> sv;
        rewrap(s,sv,75);
        os << sv[0];
        if (sv.size()>1) {
          os << "..." << std::endl;
        } else {
          os << std::endl;
        }
      }
      return;
    }

    /** \brief Output the full list of constants to \c os 
     */
    void output_list_full(std::ostream &os) const {
      os << "name unit flag value units (m,kg,s,K,A,mol,cd)" << std::endl;
      os << "  source" << std::endl;
      os << "  alternate names" << std::endl;
      os << "---------------------------------------" 
         << "---------------------------------------" << std::endl;
      for(size_t i=0;i<list.size();i++) {
        os << list[i].names[0] << " ";
        if (list[i].unit.length()==0) {
          os << "\"\" ";
        } else {
          os << list[i].unit << " ";
        }
        if (list[i].unit_flag==o2scl_const::o2scl_mks) {
          os << "MKS ";
        } else if (list[i].unit_flag==o2scl_const::o2scl_cgs) {
          os << "CGS ";
        } else if (list[i].unit_flag==fc_none) {
          os << "none ";
        } else if (list[i].unit_flag==fc_other) {
          os << "other ";
        } else {
          os << "unknown ";
        }
        os << list[i].val << " ";
        os << "(" << list[i].m << "," << list[i].k
           << "," << list[i].s << "," << list[i].K
           << "," << list[i].A << "," << list[i].mol
           << "," << list[i].cd << ")" << std::endl;
        std::vector<std::string> sv;
        rewrap(list[i].source,sv,77);
        for(size_t j=0;j<sv.size();j++) {
          os << "  " << sv[j] << std::endl;
        }
        if (list[i].names.size()>1) {
          os << "  ";
          for(size_t j=1;j<list[i].names.size();j++) {
            os << '\"' << list[i].names[j] << "\" ";
          }
          os << std::endl;
        } else {
          os << "  (no alternate names)" << std::endl;
        }
      }
      return;
    }


    /** \brief Output the full list of constants to 
        \c std::cout
    */
    void output_list_cout() const {
      output_list(std::cout);
      return;
    }

    /** \brief Output one entry from the constant database
        to \c os
    */
    void output(const find_constants::const_entry &c,
                std::ostream &os) const {
      os << "Name: " << c.names[0] << " unit: ";
      if (c.unit.length()==0) {
        os << "\"\" ";
      } else {
        os << c.unit << " ";
      }
      os << "flag: ";
      if (c.unit_flag==o2scl_const::o2scl_mks) {
        os << "MKS ";
      } else if (c.unit_flag==o2scl_const::o2scl_cgs) {
        os << "CGS ";
      } else if (c.unit_flag==fc_none) {
        os << "none ";
      } else if (c.unit_flag==fc_other) {
        os << "other ";
      } else {
        os << "unknown ";
      }
      os << "value: " << c.val << std::endl;
      os << "  (m:" << c.m << ",kg:" << c.k
         << ",s:" << c.s << ",K:" << c.K
         << ",A:" << c.A << ",mol:" << c.mol
         << ",cd:" << c.cd << ")" << std::endl;
      std::vector<std::string> sv;
      rewrap(c.source,sv,71);
      for(size_t j=0;j<sv.size();j++) {
        if (j==0) {
          os << "  Source: " << sv[0] << std::endl;
        } else {
          os << "  " << sv[j] << std::endl;
        }
      }
      if (c.names.size()>1) {
        os << "  Other names: ";
        for(size_t j=1;j<c.names.size();j++) {
          os << '\"' << c.names[j] << "\" ";
        }
        os << std::endl;
      } else {
        os << "  (no alternate names)" << std::endl;
      }
      return;
    }

   
    /** \brief Add a constant
     */
    void add_constant(const const_entry &f, int verbose=0) {

      if (verbose>1) {
        std::cout << "find_constants::add_constant() attempting "
                  << "to add constant " << f.names[0] 
                  << " with value " << f.val << std::endl;
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
        std::cout << "find_constants::add_constant() adding constant "
             << f.names[0] << " with value " << f.val << std::endl;
      }
      std::cout << "List.size(): " << list.size() << std::endl;
      list.push_back(f);
      std::cout << "List.size(): " << list.size() << std::endl;
  
      return;
    }

    
    /** \brief Remove a constant
     */
    void del_constant(std::string &name, int verbose=0) {

      if (verbose>1) {
        std::cout << "find_constants::add_constant() attempting to remove "
             << "constant named " << name << std::endl;
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
        typename std::vector<const_entry>::iterator it=list.begin();
        it+=i_match;
        if (verbose>1) {
          std::cout << "find_constants::add_constant() Removing "
               << "constant named " << name << " with value "
               << it->val << std::endl;
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

    //@}
    
  };

}

#endif
