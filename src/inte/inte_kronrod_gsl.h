/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Jerry Gagelman and Andrew W. Steiner
  
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
#ifndef O2SCL_INTE_GSL_KRONROD_H
#define O2SCL_INTE_GSL_KRONROD_H

/** \file inte_kronrod_gsl.h
    \brief File defining GSL-based integration coefficients, workspace, 
    and \ref o2scl::inte_kronrod_gsl
*/

#include <cmath>
#include <gsl/gsl_integration.h>
#include <o2scl/inte.h>
#include <o2scl/inte_gsl.h>
#include <o2scl/string_conv.h>
#include <o2scl/table_units.h>

/** \brief A namespace for the GSL adaptive Gauss-Kronrod integration
    coefficients

    The storage scheme follows that of QUADPACK: symmetry about the  
    origin permits storing only those abscissae in the interval 
    \f$ [0, 1). \f$ Similiarly for the weights. The odd-indexed
    abscissae \c xgk[1], \c xgk[3],... belong to the low-order Gauss 
    rule. For the full Gauss-Kronrod rule, the array \c wgk[] holds
    the weights corresponding to the abscissae \c xgk[]. For the
    low-order rule, \c wg[i] is the weight correpsonding to
    \c xgk[2*i+1].
 
    <b>Documentation from GSL</b>: \n    
    Gauss quadrature weights and kronrod quadrature abscissae and
    weights as evaluated with 80 decimal digit arithmetic by
    L. W. Fullerton, Bell Labs, Nov. 1981.
*/
namespace o2scl_inte_gk_coeffs {
  
  /** \brief Abscissae of the 15-point Kronrod rule */
  static const double qk15_xgk[8]={
                                   0.991455371120812639206854697526329,
                                   0.949107912342758524526189684047851,
                                   0.864864423359769072789712788640926,
                                   0.741531185599394439863864773280788,
                                   0.586087235467691130294144838258730,
                                   0.405845151377397166906606412076961,
                                   0.207784955007898467600689403773245,
                                   0.000000000000000000000000000000000
  };
  
  /** \brief Weights of the 7-point Gauss rule */
  static const double qk15_wg[4] =     
    {
     0.129484966168869693270611432679082,
     0.279705391489276667901467771423780,
     0.381830050505118944950369775488975,
     0.417959183673469387755102040816327
    };

  /** \brief Weights of the 15-point Kronrod rule */
  static const double qk15_wgk[8] =    
    {
     0.022935322010529224963732008058970,
     0.063092092629978553290700663189204,
     0.104790010322250183839876322541518,
     0.140653259715525918745189590510238,
     0.169004726639267902826583426598550,
     0.190350578064785409913256402421014,
     0.204432940075298892414161999234649,
     0.209482141084727828012999174891714
    };

  /** \brief Abscissae of the 21-point Kronrod rule */
  static const double qk21_xgk[11] =   
    {
     0.995657163025808080735527280689003,
     0.973906528517171720077964012084452,
     0.930157491355708226001207180059508,
     0.865063366688984510732096688423493,
     0.780817726586416897063717578345042,
     0.679409568299024406234327365114874,
     0.562757134668604683339000099272694,
     0.433395394129247190799265943165784,
     0.294392862701460198131126603103866,
     0.148874338981631210884826001129720,
     0.000000000000000000000000000000000
    };

  /** \brief Weights of the 10-point Gauss rule */
  static const double qk21_wg[5] =     
    {
     0.066671344308688137593568809893332,
     0.149451349150580593145776339657697,
     0.219086362515982043995534934228163,
     0.269266719309996355091226921569469,
     0.295524224714752870173892994651338
    };

  /** \brief Weights of the 21-point Kronrod rule */
  static const double qk21_wgk[11] =   
    {
     0.011694638867371874278064396062192,
     0.032558162307964727478818972459390,
     0.054755896574351996031381300244580,
     0.075039674810919952767043140916190,
     0.093125454583697605535065465083366,
     0.109387158802297641899210590325805,
     0.123491976262065851077958109831074,
     0.134709217311473325928054001771707,
     0.142775938577060080797094273138717,
     0.147739104901338491374841515972068,
     0.149445554002916905664936468389821
    };

  /** \brief Abscissae of the 31-point Kronrod rule */
  static const double qk31_xgk[16] =   
    {
     0.998002298693397060285172840152271,
     0.987992518020485428489565718586613,
     0.967739075679139134257347978784337,
     0.937273392400705904307758947710209,
     0.897264532344081900882509656454496,
     0.848206583410427216200648320774217,
     0.790418501442465932967649294817947,
     0.724417731360170047416186054613938,
     0.650996741297416970533735895313275,
     0.570972172608538847537226737253911,
     0.485081863640239680693655740232351,
     0.394151347077563369897207370981045,
     0.299180007153168812166780024266389,
     0.201194093997434522300628303394596,
     0.101142066918717499027074231447392,
     0.000000000000000000000000000000000
    };

  /** \brief Weights of the 15-point Gauss rule */
  static const double qk31_wg[8] =     
    {
     0.030753241996117268354628393577204,
     0.070366047488108124709267416450667,
     0.107159220467171935011869546685869,
     0.139570677926154314447804794511028,
     0.166269205816993933553200860481209,
     0.186161000015562211026800561866423,
     0.198431485327111576456118326443839,
     0.202578241925561272880620199967519
    };

  /** \brief Weights of the 31-point Kronrod rule */
  static const double qk31_wgk[16] =   
    {
     0.005377479872923348987792051430128,
     0.015007947329316122538374763075807,
     0.025460847326715320186874001019653,
     0.035346360791375846222037948478360,
     0.044589751324764876608227299373280,
     0.053481524690928087265343147239430,
     0.062009567800670640285139230960803,
     0.069854121318728258709520077099147,
     0.076849680757720378894432777482659,
     0.083080502823133021038289247286104,
     0.088564443056211770647275443693774,
     0.093126598170825321225486872747346,
     0.096642726983623678505179907627589,
     0.099173598721791959332393173484603,
     0.100769845523875595044946662617570,
     0.101330007014791549017374792767493
    };

  /** \brief Abscissae of the 41-point Kronrod rule */
  static const double qk41_xgk[21] =   
    {
     0.998859031588277663838315576545863,
     0.993128599185094924786122388471320,
     0.981507877450250259193342994720217,
     0.963971927277913791267666131197277,
     0.940822633831754753519982722212443,
     0.912234428251325905867752441203298,
     0.878276811252281976077442995113078,
     0.839116971822218823394529061701521,
     0.795041428837551198350638833272788,
     0.746331906460150792614305070355642,
     0.693237656334751384805490711845932,
     0.636053680726515025452836696226286,
     0.575140446819710315342946036586425,
     0.510867001950827098004364050955251,
     0.443593175238725103199992213492640,
     0.373706088715419560672548177024927,
     0.301627868114913004320555356858592,
     0.227785851141645078080496195368575,
     0.152605465240922675505220241022678,
     0.076526521133497333754640409398838,
     0.000000000000000000000000000000000
    };

  /** \brief Weights of the 20-point Gauss rule */
  static const double qk41_wg[11] =    
    {
     0.017614007139152118311861962351853,
     0.040601429800386941331039952274932,
     0.062672048334109063569506535187042,
     0.083276741576704748724758143222046,
     0.101930119817240435036750135480350,
     0.118194531961518417312377377711382,
     0.131688638449176626898494499748163,
     0.142096109318382051329298325067165,
     0.149172986472603746787828737001969,
     0.152753387130725850698084331955098
    };

  /** \brief Weights of the 41-point Kronrod rule */
  static const double qk41_wgk[21] =   
    {
     0.003073583718520531501218293246031,
     0.008600269855642942198661787950102,
     0.014626169256971252983787960308868,
     0.020388373461266523598010231432755,
     0.025882133604951158834505067096153,
     0.031287306777032798958543119323801,
     0.036600169758200798030557240707211,
     0.041668873327973686263788305936895,
     0.046434821867497674720231880926108,
     0.050944573923728691932707670050345,
     0.055195105348285994744832372419777,
     0.059111400880639572374967220648594,
     0.062653237554781168025870122174255,
     0.065834597133618422111563556969398,
     0.068648672928521619345623411885368,
     0.071054423553444068305790361723210,
     0.073030690332786667495189417658913,
     0.074582875400499188986581418362488,
     0.075704497684556674659542775376617,
     0.076377867672080736705502835038061,
     0.076600711917999656445049901530102
    };

  /** \brief Abscissae of the 51-point Kronrod rule */
  static const double qk51_xgk[26] =   
    {
     0.999262104992609834193457486540341,
     0.995556969790498097908784946893902,
     0.988035794534077247637331014577406,
     0.976663921459517511498315386479594,
     0.961614986425842512418130033660167,
     0.942974571228974339414011169658471,
     0.920747115281701561746346084546331,
     0.894991997878275368851042006782805,
     0.865847065293275595448996969588340,
     0.833442628760834001421021108693570,
     0.797873797998500059410410904994307,
     0.759259263037357630577282865204361,
     0.717766406813084388186654079773298,
     0.673566368473468364485120633247622,
     0.626810099010317412788122681624518,
     0.577662930241222967723689841612654,
     0.526325284334719182599623778158010,
     0.473002731445714960522182115009192,
     0.417885382193037748851814394594572,
     0.361172305809387837735821730127641,
     0.303089538931107830167478909980339,
     0.243866883720988432045190362797452,
     0.183718939421048892015969888759528,
     0.122864692610710396387359818808037,
     0.061544483005685078886546392366797,
     0.000000000000000000000000000000000
    };

  /** \brief Weights of the 25-point Gauss rule */
  static const double qk51_wg[13] =    
    {
     0.011393798501026287947902964113235,
     0.026354986615032137261901815295299,
     0.040939156701306312655623487711646,
     0.054904695975835191925936891540473,
     0.068038333812356917207187185656708,
     0.080140700335001018013234959669111,
     0.091028261982963649811497220702892,
     0.100535949067050644202206890392686,
     0.108519624474263653116093957050117,
     0.114858259145711648339325545869556,
     0.119455763535784772228178126512901,
     0.122242442990310041688959518945852,
     0.123176053726715451203902873079050
    };

  /** \brief Weights of the 51-point Kronrod rule */
  static const double qk51_wgk[26] =   
    {
     0.001987383892330315926507851882843,
     0.005561932135356713758040236901066,
     0.009473973386174151607207710523655,
     0.013236229195571674813656405846976,
     0.016847817709128298231516667536336,
     0.020435371145882835456568292235939,
     0.024009945606953216220092489164881,
     0.027475317587851737802948455517811,
     0.030792300167387488891109020215229,
     0.034002130274329337836748795229551,
     0.037116271483415543560330625367620,
     0.040083825504032382074839284467076,
     0.042872845020170049476895792439495,
     0.045502913049921788909870584752660,
     0.047982537138836713906392255756915,
     0.050277679080715671963325259433440,
     0.052362885806407475864366712137873,
     0.054251129888545490144543370459876,
     0.055950811220412317308240686382747,
     0.057437116361567832853582693939506,
     0.058689680022394207961974175856788,
     0.059720340324174059979099291932562,
     0.060539455376045862945360267517565,
     0.061128509717053048305859030416293,
     0.061471189871425316661544131965264,
     0.061580818067832935078759824240066
    };

  /** \brief Abscissae of the 61-point Kronrod rule */
  static const double qk61_xgk[31] =   
    {
     0.999484410050490637571325895705811,
     0.996893484074649540271630050918695,
     0.991630996870404594858628366109486,
     0.983668123279747209970032581605663,
     0.973116322501126268374693868423707,
     0.960021864968307512216871025581798,
     0.944374444748559979415831324037439,
     0.926200047429274325879324277080474,
     0.905573307699907798546522558925958,
     0.882560535792052681543116462530226,
     0.857205233546061098958658510658944,
     0.829565762382768397442898119732502,
     0.799727835821839083013668942322683,
     0.767777432104826194917977340974503,
     0.733790062453226804726171131369528,
     0.697850494793315796932292388026640,
     0.660061064126626961370053668149271,
     0.620526182989242861140477556431189,
     0.579345235826361691756024932172540,
     0.536624148142019899264169793311073,
     0.492480467861778574993693061207709,
     0.447033769538089176780609900322854,
     0.400401254830394392535476211542661,
     0.352704725530878113471037207089374,
     0.304073202273625077372677107199257,
     0.254636926167889846439805129817805,
     0.204525116682309891438957671002025,
     0.153869913608583546963794672743256,
     0.102806937966737030147096751318001,
     0.051471842555317695833025213166723,
     0.000000000000000000000000000000000
    };

  /** \brief Weights of the 30-point Gauss rule */
  static const double qk61_wg[15] =    
    {
     0.007968192496166605615465883474674,
     0.018466468311090959142302131912047,
     0.028784707883323369349719179611292,
     0.038799192569627049596801936446348,
     0.048402672830594052902938140422808,
     0.057493156217619066481721689402056,
     0.065974229882180495128128515115962,
     0.073755974737705206268243850022191,
     0.080755895229420215354694938460530,
     0.086899787201082979802387530715126,
     0.092122522237786128717632707087619,
     0.096368737174644259639468626351810,
     0.099593420586795267062780282103569,
     0.101762389748405504596428952168554,
     0.102852652893558840341285636705415
    };

  /** \brief Weights of the 61-point Kronrod rule */
  static const double qk61_wgk[31] =   
    {
     0.001389013698677007624551591226760,
     0.003890461127099884051267201844516,
     0.006630703915931292173319826369750,
     0.009273279659517763428441146892024,
     0.011823015253496341742232898853251,
     0.014369729507045804812451432443580,
     0.016920889189053272627572289420322,
     0.019414141193942381173408951050128,
     0.021828035821609192297167485738339,
     0.024191162078080601365686370725232,
     0.026509954882333101610601709335075,
     0.028754048765041292843978785354334,
     0.030907257562387762472884252943092,
     0.032981447057483726031814191016854,
     0.034979338028060024137499670731468,
     0.036882364651821229223911065617136,
     0.038678945624727592950348651532281,
     0.040374538951535959111995279752468,
     0.041969810215164246147147541285970,
     0.043452539701356069316831728117073,
     0.044814800133162663192355551616723,
     0.046059238271006988116271735559374,
     0.047185546569299153945261478181099,
     0.048185861757087129140779492298305,
     0.049055434555029778887528165367238,
     0.049795683427074206357811569379942,
     0.050405921402782346840893085653585,
     0.050881795898749606492297473049805,
     0.051221547849258772170656282604944,
     0.051426128537459025933862879215781,
     0.051494729429451567558340433647099
    };
  
}

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Integration workspace for the GSL integrators

      This is a simple rewrite of \c inte_integration_workspace
      into a class.

      QUADPACK workspace documentation:
      \verbatim
      c    parameters (meaning at output)
      c       limit  - integer
      c                maximum number of error estimates the list
      c                can contain
      c          last   - integer
      c                number of error estimates currently in the list
      c          maxerr - integer
      c                maxerr points to the nrmax-th largest error
      c                estimate currently in the list
      c          ermax  - double precision
      c                nrmax-th largest error estimate
      c                ermax = elist(maxerr)
      c          elist  - double precision
      c                vector of dimension last containing
      c                the error estimates
      c          iord   - integer
      c                vector of dimension last, the first k elements
      c                of which contain pointers to the error
      c                estimates, such that
      c                elist(iord(1)),...,  elist(iord(k))
      c                form a decreasing sequence, with
      c                k = last if last.le.(limit/2+2), and
      c                k = limit+1-last otherwise
      c          nrmax  - integer
      c      maxerr = iord(nrmax)
      \endverbatim
      \verbatim
      c     alist   - real
      c               vector of dimension at least limit, the first
      c                last  elements of which are the left
      c               end points of the subintervals in the partition
      c               of the given integration range (a,b)
      c            blist   - real
      c               vector of dimension at least limit, the first
      c                last  elements of which are the right
      c               end points of the subintervals in the partition
      c               of the given integration range (a,b)
      c            rlist   - real
      c               vector of dimension at least limit, the first
      c                last  elements of which are the
      c               integral approximations on the subintervals
      c            elist   - real
      c               vector of dimension at least limit, the first
      c                last  elements of which are the moduli of the
      c               absolute error estimates on the subintervals
      c            iord    - integer
      c               vector of dimension at least limit, the first k
      c               elements of which are pointers to the
      c               error estimates over the subintervals,
      c               such that elist(iord(1)), ...,
      c               elist(iord(k)) form a decreasing sequence,
      c               with k = last if last.le.(limit/2+2), and
      c               k = limit+1-last otherwise
      c            last    - integer
      c               number of subintervals actually produced in the
      c               subdivision process
      \endverbatim
      
  */
  class inte_workspace_gsl {
    
  public:
    
    inte_workspace_gsl();
    
    ~inte_workspace_gsl();

    /// Maximum number of subintervals allocated
    size_t limit;
    /// Current number of subintervals being used
    size_t size;		
    /// Counter for extrapolation routine
    size_t nrmax;
    /// Index of current subinterval
    size_t i;
    /// Depth of subdivisions reached
    size_t maximum_level;
    /// Left endpoints of subintervals
    double *alist;
    /// Right endpoints of subintervals
    double *blist;
    /// Integral approximations on subintervals
    double *rlist;
    /// Integral error estimates
    double *elist;
    /// Linear ordering vector for sort routine
    size_t *order;	
    /// Numbers of subdivisions made
    size_t *level;	
    
    /// Allocate a workspace
    int allocate(size_t sz);

    /// Free allocated workspace memory
    int free();
    
    /** \brief Initialize the workspace for an integration with limits
	\c a and \c b.
    */
    int initialise(double a, double b);

    /** \brief Update the workspace with the result and error
	from the first integration
    */
    int set_initial_result(double result, double error);

    /** \brief Retrieve the ith result from the 
	workspace stack
	
	The workspace variable \c i is used to specify which 
	interval is requested.
    */
    int retrieve(double *a, double *b, double *r, double *e) const;

    /** \brief Create a table from the current workspace
     */
    void make_table(o2scl::table_units<> &t);

    /** \brief Sort the workspace stack
	
	This routine maintains the descending ordering in the list of
	the local error estimated resulting from the interval
	subdivision process. at each call two error estimates are
	inserted using the sequential search method, top-down for the
	largest error estimate and bottom-up for the smallest error
	estimate.

	Originally written in QUADPACK by R. Piessens and E. de
	Doncker, translated into C for GSL by Brian Gough, and then
	rewritten for \o2.
    */
    int qpsrt();

    /** \brief Determine which new subinterval to add to the workspace
	stack and perform update
    */
    int update(double a1, double b1, double area1, double error1,
	       double a2, double b2, double area2, double error2);
    
    /// Add up all of the contributions to construct the final result
    double sum_results();

    /** \brief Test whether the proposed subdivision falls 
	before floating-point precision
    */
    int subinterval_too_small(double a1, double a2, double b2);

    /// Push a new interval to the workspace stack
    void append_interval(double a1, double b1, double area1, double error1);

  };

  /** \brief Basic Gauss-Kronrod integration class (GSL)
      
      This class provides the basic Gauss-Kronrod integration function
      and integration workspace for some of the GSL-based integration
      classes. 

      The main function of interest is \ref set_rule(), which
      sets the integration rule for the GSL integration classes
      which inherit from this base class. The argument to
      set rule should be selected from the following list:
      - 1: 15-point Gauss-Kronrod integration rule (default)
      - 2: 21-point rule
      - 3: 31-point rule
      - 4: 41-point rule
      - 5: 51-point rule
      - 6: 61-point rule
      Any value other than 1-6 forces the error handler to be
      called. All classes default to the 15-point rule,
      except for \ref o2scl::inte_qags_gsl which defaults to
      the 21-point rule for singularities. 

      The integration coefficients for use with this class and
      associated children are stored in the \ref o2scl_inte_gk_coeffs
      namespace.

      \future Convert 'double *' objects to ubvector
  */
  template<class func_t> class inte_kronrod_gsl : public inte_gsl,
                                                  public inte<func_t> {
      
  protected:
      
    /// The integration workspace
    inte_workspace_gsl *w;
      
    /// Size of Gauss-Kronrod arrays
    int n_gk;

    /// Gauss-Kronrod abscissae pointer
    const double *x_gk;
      
    /// Gauss weight pointer
    const double *w_g;
      
    /// Gauss-Kronrod weight pointer
    const double *w_gk;
      
    /// Scratch space
    double *f_v1;

    /// Scratch space
    double *f_v2;

  public:

    inte_kronrod_gsl() {
      w=new inte_workspace_gsl;
      w->allocate(1000);
      n_gk=0;
      set_rule(1);
    }
      
    ~inte_kronrod_gsl() {
      w->free();
      delete w;
      if (n_gk>0) {
        delete[] f_v1;
        delete[] f_v2;
      }
    }

    /** \brief Get a const reference to the workspace
     */
    inte_workspace_gsl &get_workspace() {
      return *w;
    }
    
    /** \brief Get the Gauss-Kronrod integration rule

        This returns the index of the GSL integration rule
        a number between 1 and 6 (inclusive)
    */
    int get_rule() {
      switch (n_gk) {
      case 8:
        return 1;
      case 11:
        return 2;
      case 16:
        return 3;
      case 21:
        return 4;
      case 26:
        return 5;
      case 31:
        return 6;
      }
      return 0;
    }
    
    /// Set the Gauss-Kronrod integration rule to be used
    void set_rule(int rule) {
	
      using namespace o2scl_inte_gk_coeffs;
	
      if (n_gk > 0) {
        delete[] f_v1;
        delete[] f_v2;
      }

      switch (rule) {
      case 1:
        n_gk=8;
        x_gk=qk15_xgk;
        w_g=qk15_wg;
        w_gk=qk15_wgk;
        break;
      case 2:
        n_gk=11;
        x_gk=qk21_xgk;
        w_g=qk21_wg;
        w_gk=qk21_wgk;
        break;
      case 3:
        n_gk=16;
        x_gk=qk31_xgk;
        w_g=qk31_wg;
        w_gk=qk31_wgk;
        break;
      case 4:
        n_gk=21;
        x_gk=qk41_xgk;
        w_g=qk41_wg;
        w_gk=qk41_wgk;
        break;
      case 5:
        n_gk=26;
        x_gk=qk51_xgk;
        w_g=qk51_wg;
        w_gk=qk51_wgk;
        break;
      case 6:
        n_gk=31;
        x_gk=qk61_xgk;
        w_g=qk61_wg;
        w_gk=qk61_wgk;
        break;
      default:
        O2SCL_ERR("Invalid rule.",exc_einval);
        break;
      }

      f_v1=new double[n_gk];
      f_v2=new double[n_gk];

      return;
    }

    /** \brief Set the limit for the number of subdivisions of 
        the integration region (default 1000)
	
        If the value of \c size is zero, the error handler will
        be called.
    */
    int set_limit(size_t lim) {
      if (lim==0) {
        O2SCL_ERR2("Requested zero limit in ",
                   "inte_kronrod_gsl::set_limit().",exc_einval);
      }
      w->free();
      w->allocate(lim);
      return 0;
    }

    /** \brief The base Gauss-Kronrod integration function template
	  
        Given abcissas and weights, this performs the integration of
        \c func between \c a and \c b, providing a result with
        uncertainties.
	
        The Gauss-Kronrod rule uses \f$ 2m+1 \f$ abscissae \f$ x_1,
        \ldots, x_{2m+1} \in (a, b) \f$ to estimate the integral of a
        function as a linear combination of values,
        \f[ 
        \int_a^b f~dx \; \approx \; 
        Q_m^{GK}f \; \equiv \;
        \sum_{k=1}^{2m+1} \beta_k f(x_k),
        \f]
        where the weights \f$ \beta_1,\ldots, \beta_{2m+1} \f$ are
        intrinsic to the abscissae. The data are designed so that the
        even-indexed abscissae yield a coarser estimate,
        \f[
        Q_m^{G}f \; \equiv \;
        \sum_{j=1}^{m} \alpha_j f(x_{2j}),
        \f]
        and their difference \f$ |Q_m^Gf - Q_m^{GK}f| \f$ is the "raw"
        error estimate. The various quantities that the function
        computes are
        \f{eqnarray*}{
        \mathtt{result} &=& Q_m^{GK}f, \\
        \mathtt{resabs} &=& Q_m^{GK}|f|, \\
        \mathtt{resasc} &=& Q_m^{GK}(|f| - \mu),
        \quad \mu \equiv \frac{Q_m^{GK}f}{b - a}.
        \f}
        The "absolute" error \c abserr is computed from the raw error
        value using the function \ref inte_gsl::rescale_error.
	
        This function is designed for use with the values given in the
        o2scl_inte_gk_coeffs namespace.
	
        This function never calls the error handler.

        \future This function, in principle, is an integration
        object in and of itself, and could be implemented separately
        as an object of type \ref o2scl::inte.
    */
    template<class func2_t> void gauss_kronrod_base
    (func2_t &func, double a, double b, double *result, 
     double *abserr, double *resabs, double *resasc) {
	
      const double center=0.5*(a+b);
      const double half_length=0.5*(b-a);
      const double abs_half_length=fabs (half_length);
	  
      double f_center;
      f_center=func(center);
	  
      double result_gauss=0.0;
      double result_kronrod=f_center*this->w_gk[this->n_gk-1];
 
      double result_abs=fabs (result_kronrod);
      double result_asc=0.0;
      double mean=0.0, err=0.0;
      
      if (this->n_gk % 2 == 0) {
        result_gauss=f_center*this->w_g[this->n_gk / 2-1];
      }
      
      // Sum up results from left half of interval
      for (int j=0; j < (this->n_gk-1) / 2; j++) {

        const int jtw=j*2+1;        
        const double abscissa=half_length*this->x_gk[jtw];
	      
        double fval1, fval2;
        fval1=func(center-abscissa);
        fval2=func(center+abscissa);
	      
        const double fsum=fval1+fval2;
        this->f_v1[jtw]=fval1;
        this->f_v2[jtw]=fval2;
        result_gauss+=this->w_g[j]*fsum;
        result_kronrod+=this->w_gk[jtw]*fsum;
        result_abs+=this->w_gk[jtw]*(fabs(fval1)+fabs(fval2));
      }
      
      // Sum up results from right half of interval
      for (int j=0; j < this->n_gk / 2; j++) {

        int jtwm1=j*2;
        const double abscissa=half_length*this->x_gk[jtwm1];

        double fval1, fval2;
        fval1=func(center-abscissa);
        fval2=func(center+abscissa);

        this->f_v1[jtwm1]=fval1;
        this->f_v2[jtwm1]=fval2;
        result_kronrod+=this->w_gk[jtwm1]*(fval1+fval2);
        result_abs+=this->w_gk[jtwm1]*(fabs(fval1)+fabs(fval2));
      }
      
      mean=result_kronrod*0.5;
      
      result_asc=this->w_gk[this->n_gk-1]*fabs(f_center-mean);
	
      for (int j=0;j<this->n_gk-1;j++) {
        result_asc+=this->w_gk[j]*(fabs(this->f_v1[j]-mean)+
                                   fabs(this->f_v2[j]-mean));
      }
      
      /* Scale by the width of the integration region */
      
      err=(result_kronrod-result_gauss)*half_length;
      
      result_kronrod*=half_length;
      result_abs*=abs_half_length;
      result_asc*=abs_half_length;
      
      *result=result_kronrod;
      *resabs=result_abs;
      *resasc=result_asc;
      *abserr=rescale_error(err,result_abs,result_asc);

      return;
    }

    /** \brief Integration wrapper for user-specified function type
     */
    virtual void gauss_kronrod
    (func_t &func, double a, double b,
     double *result, double *abserr, double *resabs, double *resasc) {
      return gauss_kronrod_base<func_t>
        (func,a,b,result,abserr,resabs,resasc);
					  
    }
    
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
