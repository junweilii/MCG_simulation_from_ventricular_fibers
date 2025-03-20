/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#define _pval pval
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if !NRNGPU
#undef exp
#define exp hoc_Exp
#endif
 
#define nrn_init _nrn_init__INa
#define _nrn_initial _nrn_initial__INa
#define nrn_cur _nrn_cur__INa
#define _nrn_current _nrn_current__INa
#define nrn_jacob _nrn_jacob__INa
#define nrn_state _nrn_state__INa
#define _net_receive _net_receive__INa 
#define rates rates__INa 
#define states states__INa 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *hoc_getarg(int);
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gnabar _p[0]
#define gnabar_columnindex 0
#define gnac _p[1]
#define gnac_columnindex 1
#define Tauact _p[2]
#define Tauact_columnindex 2
#define Tauinactf _p[3]
#define Tauinactf_columnindex 3
#define Tauinacts _p[4]
#define Tauinacts_columnindex 4
#define ina _p[5]
#define ina_columnindex 5
#define m _p[6]
#define m_columnindex 6
#define h _p[7]
#define h_columnindex 7
#define n _p[8]
#define n_columnindex 8
#define Dm _p[9]
#define Dm_columnindex 9
#define Dh _p[10]
#define Dh_columnindex 10
#define Dn _p[11]
#define Dn_columnindex 11
#define ena _p[12]
#define ena_columnindex 12
#define _g _p[13]
#define _g_columnindex 13
#define _ion_ena	*(_ppvar[0].get<double*>())
#define _ion_ina	*_ppvar[1].get<double*>()
#define _ion_dinadv	*_ppvar[2].get<double*>()
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_alp(void);
 static void _hoc_bet(void);
 static void _hoc_rates(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mechtype);
#endif
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 {"setdata_INa", _hoc_setdata},
 {"alp_INa", _hoc_alp},
 {"bet_INa", _hoc_bet},
 {"rates_INa", _hoc_rates},
 {0, 0}
};
#define alp alp_INa
#define bet bet_INa
 extern double alp( double , double );
 extern double bet( double , double );
 /* declare global and static user variables */
#define htau htau_INa
 double htau = 0;
#define hinf hinf_INa
 double hinf = 0;
#define mtau mtau_INa
 double mtau = 0;
#define minf minf_INa
 double minf = 0;
#define ntau ntau_INa
 double ntau = 0;
#define ninf ninf_INa
 double ninf = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 {"gnabar_INa", 0, 1e+09},
 {0, 0, 0}
};
 static HocParmUnits _hoc_parm_units[] = {
 {"mtau_INa", "ms"},
 {"htau_INa", "ms"},
 {"ntau_INa", "ms"},
 {"gnabar_INa", "S/cm2"},
 {"gnac_INa", "S/cm2"},
 {"Tauact_INa", "ms"},
 {"Tauinactf_INa", "ms"},
 {"Tauinacts_INa", "ms"},
 {"ina_INa", "mA/cm2"},
 {0, 0}
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 static double n0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 {"minf_INa", &minf_INa},
 {"hinf_INa", &hinf_INa},
 {"ninf_INa", &ninf_INa},
 {"mtau_INa", &mtau_INa},
 {"htau_INa", &htau_INa},
 {"ntau_INa", &ntau_INa},
 {0, 0}
};
 static DoubVec hoc_vdoub[] = {
 {0, 0, 0}
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, Memb_list*, int);
static void nrn_state(NrnThread*, Memb_list*, int);
 static void nrn_cur(NrnThread*, Memb_list*, int);
static void  nrn_jacob(NrnThread*, Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, Memb_list*, int);
static void _ode_matsol(NrnThread*, Memb_list*, int);
 
#define _cvode_ieq _ppvar[3].literal_value<int>()
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"INa",
 "gnabar_INa",
 "gnac_INa",
 "Tauact_INa",
 "Tauinactf_INa",
 "Tauinacts_INa",
 0,
 "ina_INa",
 0,
 "m_INa",
 "h_INa",
 "n_INa",
 0,
 0};
 static Symbol* _na_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 14, _prop);
 	/*initialize range parameters*/
 	gnabar = 0.004;
 	gnac = 3e-06;
 	Tauact = 1;
 	Tauinactf = 1;
 	Tauinacts = 1;
 	_prop->param = _p;
 	_prop->param_size = 14;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0] = &prop_ion->param[0]; /* ena */
 	_ppvar[1] = &prop_ion->param[3]; /* ina */
 	_ppvar[2] = &prop_ion->param[4]; /* _ion_dinadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 {0, 0}
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 extern "C" void _ina_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("na", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  register_nmodl_text_and_filename(_mechtype);
#endif
  hoc_register_prop_size(_mechtype, 14, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 INa ina.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double _zk ;
static int _reset;
static const char *modelname = "Cardiac sodium current";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(double);
 static int _deriv1_advance = 0;
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist2[3]; static double _dlist2[3];
 static double _savstate1[3], *_temp1 = _savstate1;
 static int _slist1[3], _dlist1[3];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / mtau ;
   Dh = ( hinf - h ) / htau ;
   Dn = ( ninf - n ) / ntau ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rates ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / htau )) ;
 Dn = Dn  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ntau )) ;
  return 0;
}
 /*END CVODE*/
 
static int states () {_reset=0;
 { static int _recurse = 0;
 int _counte = -1;
 if (!_recurse) {
 _recurse = 1;
 {int _id; for(_id=0; _id < 3; _id++) { _savstate1[_id] = _p[_slist1[_id]];}}
 error = newton<3>(_slist2, _p, states, _dlist2);
 _recurse = 0; if(error) {abort_run(error);}}
 {
   rates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / mtau ;
   Dh = ( hinf - h ) / htau ;
   Dn = ( ninf - n ) / ntau ;
   {int _id; for(_id=0; _id < 3; _id++) {
if (_deriv1_advance) {
 _dlist2[++_counte] = _p[_dlist1[_id]] - (_p[_slist1[_id]] - _savstate1[_id])/dt;
 }else{
_dlist2[++_counte] = _p[_slist1[_id]] - _savstate1[_id];}}}
 } }
 return _reset;}
 
double alp (  double _lv , double _li ) {
   double _lalp;
 if ( _li  == 0.0 ) {
     _lalp = ( - _lv - 47.0 ) / ( exp ( - 0.1 * ( _lv + 47.0 ) ) - 1.0 ) ;
     }
   else if ( _li  == 1.0 ) {
     _lalp = 0.126 * exp ( - 0.25 * ( _lv + 77.0 ) ) ;
     }
   else if ( _li  == 2.0 ) {
     _lalp = 0.055 * exp ( - 0.25 * ( _lv + 78.0 ) ) / ( exp ( - 0.2 * ( _lv + 78.0 ) ) + 1.0 ) ;
     }
   
return _lalp;
 }
 
static void _hoc_alp(void) {
  double _r;
   _r =  alp (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double bet (  double _lv , double _li ) {
   double _lbet;
 if ( _li  == 0.0 ) {
     _lbet = 40.0 * exp ( - 0.056 * ( _lv + 72.0 ) ) ;
     }
   else if ( _li  == 1.0 ) {
     _lbet = 1.7 / ( exp ( - 0.082 * ( _lv + 22.5 ) ) + 1.0 ) ;
     }
   else if ( _li  == 2.0 ) {
     _lbet = 0.3 / ( exp ( - 0.1 * ( _lv + 32.0 ) ) + 1.0 ) ;
     }
   
return _lbet;
 }
 
static void _hoc_bet(void) {
  double _r;
   _r =  bet (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static int  rates (  double _lv ) {
   double _la , _lb ;
 _la = alp ( _threadargscomma_ _lv , 0.0 ) ;
   _lb = bet ( _threadargscomma_ _lv , 0.0 ) ;
   mtau = 1.0 / ( _la + _lb ) ;
   minf = _la / ( _la + _lb ) ;
   _la = alp ( _threadargscomma_ _lv , 1.0 ) ;
   _lb = bet ( _threadargscomma_ _lv , 1.0 ) ;
   htau = 1.0 / ( _la + _lb ) ;
   hinf = _la / ( _la + _lb ) ;
   _la = alp ( _threadargscomma_ _lv , 2.0 ) ;
   _lb = bet ( _threadargscomma_ _lv , 2.0 ) ;
   ntau = 1.0 / ( _la + _lb ) ;
   ninf = _la / ( _la + _lb ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   _r = 1.;
 rates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 3;}
 
static void _ode_spec(NrnThread* _nt, Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
     _ode_spec1 ();
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 3; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(NrnThread* _nt, Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  m = m0;
  n = n0;
 {
   rates ( _threadargscomma_ v ) ;
   m = minf ;
   h = hinf ;
   n = ninf ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(NrnThread* _nt, Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ena = _ion_ena;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   ina = ( gnabar * m * m * m * h * n + gnac ) * ( v - ena ) ;
   }
 _current += ina;

} return _current;
}

static void nrn_cur(NrnThread* _nt, Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ena = _ion_ena;
 _g = _nrn_current(_v + .001);
 	{ double _dina;
  _dina = ina;
 _rhs = _nrn_current(_v);
  _ion_dinadv += (_dina - ina)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(NrnThread* _nt, Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(NrnThread* _nt, Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
double _dtsav = dt;
if (secondorder) { dt *= 0.5; }
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ena = _ion_ena;
 { error = _deriv1_advance = 1;
 derivimplicit(_ninits, 3, _slist1, _dlist1, _p, &t, dt, states, &_temp1);
_deriv1_advance = 0;
 if(error){fprintf(stderr,"at line 47 in file ina.mod:\nBREAKPOINT {\n"); nrn_complain(_p); abort_run(error);}
    if (secondorder) {
    int _i;
    for (_i = 0; _i < 3; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 } }}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = m_columnindex;  _dlist1[0] = Dm_columnindex;
 _slist1[1] = h_columnindex;  _dlist1[1] = Dh_columnindex;
 _slist1[2] = n_columnindex;  _dlist1[2] = Dn_columnindex;
 _slist2[0] = h_columnindex;
 _slist2[1] = m_columnindex;
 _slist2[2] = n_columnindex;
_first = 0;
}

#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mech_type) {
    const char* nmodl_filename = "ina.mod";
    const char* nmodl_file_text = 
  "TITLE Cardiac sodium current\n"
  ": from BEELER & REUTER, J.Physiol, 1977\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX INa\n"
  "	USEION na READ ena WRITE ina\n"
  "	RANGE  gnabar, gnac, ina, Tauact, Tauinactf, Tauinacts\n"
  "	GLOBAL minf, hinf, ninf, mtau, htau, ntau\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "	(mM) = (milli/liter)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	gnabar=0.004 (S/cm2) <0,1e9> \n"
  "	gnac = 0.000003 (S/cm2)\n"
  "	Tauact= 1 (ms) \n"
  "	Tauinactf=1 (ms) \n"
  "	Tauinacts=1 (ms)\n"
  "}\n"
  "\n"
  "STATE { : m h j\n"
  "	m h n \n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v (mV)\n"
  "	celsius (degC) : 37\n"
  "	ena (mV)\n"
  "	ina (mA/cm2)\n"
  "	minf hinf ninf\n"
  "	mtau (ms)\n"
  "	htau (ms)\n"
  "	ntau (ms)  \n"
  "}\n"
  "LOCAL k\n"
  "INITIAL { : m h j\n"
  "	rates(v)\n"
  "	m = minf\n"
  "	h = hinf\n"
  "	n = ninf\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states METHOD derivimplicit\n"
  "	ina = (gnabar*m*m*m*h*n + gnac)*(v - ena)\n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "	rates(v)\n"
  "	m' = (minf - m)/mtau\n"
  "	h' = (hinf - h)/htau\n"
  "	n' = (ninf - n)/ntau\n"
  "}\n"
  "UNITSOFF\n"
  "FUNCTION alp(v(mV),i) (/ms) { \n"
  ":LOCAL q10 :order m,h,j\n"
  "	if (i==0) {\n"
  "		alp = (- v - 47)/(exp(-0.1*(v + 47))-1)\n"
  "	}else if (i==1){\n"
  "		alp = 0.126*exp(-0.25*(v + 77))\n"
  "	}else if (i==2){\n"
  "		alp = 0.055*exp(-0.25*(v + 78))/(exp(-0.2*(v + 78))+1)\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION bet(v(mV),i)(/ms) { \n"
  ":LOCAL q10 :order m,h,n\n"
  "	if (i==0) {\n"
  "		bet = 40*exp(-0.056*(v + 72))\n"
  "	}else if (i==1){\n"
  "		bet = 1.7/(exp(-0.082*(v + 22.5))+1)\n"
  "	}else if (i==2){\n"
  "		bet = 0.3/(exp(-0.1*(v + 32))+1)\n"
  "	}\n"
  "}\n"
  "\n"
  "PROCEDURE rates(v) { : m h j\n"
  "LOCAL a, b\n"
  "	:TABLE minf, hinf, ninf, mtau, htau, ntau DEPEND celsius FROM -100 TO 100 WITH 200\n"
  "	a = alp(v,0)  b=bet(v,0)\n"
  "	mtau = 1/(a + b)\n"
  "	minf = a/(a + b)\n"
  "	a = alp(v,1)  b = bet(v,1)\n"
  "	htau = 1/(a + b)\n"
  "	hinf = a/(a + b)\n"
  "	a = alp(v,2)  b=bet(v,2)\n"
  "	ntau = 1/(a + b)\n"
  "	ninf = a/(a + b)\n"
  "}\n"
  "UNITSON\n"
  ;
    hoc_reg_nmodl_filename(mech_type, nmodl_filename);
    hoc_reg_nmodl_text(mech_type, nmodl_file_text);
}
#endif
