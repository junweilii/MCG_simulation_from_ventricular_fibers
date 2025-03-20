/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
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
 
#define nrn_init _nrn_init__HalfGap
#define _nrn_initial _nrn_initial__HalfGap
#define nrn_cur _nrn_cur__HalfGap
#define _nrn_current _nrn_current__HalfGap
#define nrn_jacob _nrn_jacob__HalfGap
#define nrn_state _nrn_state__HalfGap
#define _net_receive _net_receive__HalfGap 
#define setRandom setRandom__HalfGap 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *hoc_getarg(int);
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gmax _p[0]
#define gmax_columnindex 0
#define gmin _p[1]
#define gmin_columnindex 1
#define vhalf _p[2]
#define vhalf_columnindex 2
#define meang _p[3]
#define meang_columnindex 3
#define meant _p[4]
#define meant_columnindex 4
#define drift _p[5]
#define drift_columnindex 5
#define rg _p[6]
#define rg_columnindex 6
#define rt _p[7]
#define rt_columnindex 7
#define id _p[8]
#define id_columnindex 8
#define g _p[9]
#define g_columnindex 9
#define i _p[10]
#define i_columnindex 10
#define v _p[11]
#define v_columnindex 11
#define _g _p[12]
#define _g_columnindex 12
#define _tsav _p[13]
#define _tsav_columnindex 13
#define _nd_area  *_ppvar[0].get<double*>()
#define donotuse	*_ppvar[2].get<double*>()
#define _p_donotuse _ppvar[2].literal_value<void*>()
#define vgap	*_ppvar[3].get<double*>()
#define _p_vgap _ppvar[3].literal_value<void*>()
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 static int hoc_nrnpointerindex =  2;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static double _hoc_getpar(void*);
 static double _hoc_gv(void*);
 static double _hoc_mynormrand(void*);
 static double _hoc_setRandom(void*);
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
 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(Object* _ho) { void* create_point_process(int, Object*);
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt(void*);
 static double _hoc_loc_pnt(void* _vptr) {double loc_point_process(int, void*);
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(void* _vptr) {double has_loc_point(void*);
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(void* _vptr) {
 double get_loc_point_process(void*); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 {0, 0}
};
 static Member_func _member_func[] = {
 {"loc", _hoc_loc_pnt},
 {"has_loc", _hoc_has_loc},
 {"get_loc", _hoc_get_loc_pnt},
 {"getpar", _hoc_getpar},
 {"gv", _hoc_gv},
 {"mynormrand", _hoc_mynormrand},
 {"setRandom", _hoc_setRandom},
 {0, 0}
};
#define getpar getpar_HalfGap
#define gv gv_HalfGap
#define mynormrand mynormrand_HalfGap
 extern double getpar( _threadargsproto_ );
 extern double gv( _threadargsprotocomma_ double );
 extern double mynormrand( _threadargsprotocomma_ double , double );
 /* declare global and static user variables */
#define event event_HalfGap
 double event = 0;
#define slope4 slope4_HalfGap
 double slope4 = 10;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 {0, 0, 0}
};
 static HocParmUnits _hoc_parm_units[] = {
 {"slope4_HalfGap", "/millivolt"},
 {"event_HalfGap", "ms"},
 {"gmax", "nanosiemens"},
 {"gmin", "nanosiemens"},
 {"vhalf", "millivolt"},
 {"meang", "nanosiemens"},
 {"meant", "ms"},
 {"g", "nanosiemens"},
 {"i", "nanoamp"},
 {"vgap", "millivolt"},
 {0, 0}
};
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 {"slope4_HalfGap", &slope4_HalfGap},
 {"event_HalfGap", &event_HalfGap},
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
 static void _hoc_destroy_pnt(void* _vptr) {
   destroy_point_process(_vptr);
}
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"HalfGap",
 "gmax",
 "gmin",
 "vhalf",
 "meang",
 "meant",
 "drift",
 "rg",
 "rt",
 "id",
 0,
 "g",
 "i",
 0,
 0,
 "donotuse",
 "vgap",
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 14, _prop);
 	/*initialize range parameters*/
 	gmax = 1;
 	gmin = 1;
 	vhalf = 0;
 	meang = 30;
 	meant = 1e+06;
 	drift = 0;
 	rg = 0;
 	rt = 0;
 	id = 0;
  }
 	_prop->param = _p;
 	_prop->param_size = 14;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
 
#define _tqitem &(_ppvar[4])
 static void _net_receive(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 extern "C" void _halfgapm1_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 1,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  register_nmodl_text_and_filename(_mechtype);
#endif
  hoc_register_prop_size(_mechtype, 14, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "pointer");
  hoc_register_dparam_semantics(_mechtype, 3, "pointer");
  hoc_register_dparam_semantics(_mechtype, 4, "netsend");
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 HalfGap halfgapm1.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static const char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int setRandom(_threadargsproto_);
 
double gv ( _threadargsprotocomma_ double _lx ) {
   double _lgv;
 _lgv = ( gmax - gmin ) / ( 1.0 + exp ( slope4 * ( vhalf - _lx ) ) ) + gmin ;
   
return _lgv;
 }
 
static double _hoc_gv(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  gv ( _p, _ppvar, _thread, _nt, *getarg(1) );
 return(_r);
}
 
double getpar ( _threadargsproto_ ) {
   double _lgetpar;
 gmax = mynormrand ( _threadargscomma_ meang / 1.0 , rg ) * 1.0 ;
   if ( gmax < 0.0 ) {
     gmax = 0.0 ;
     }
   if ( gmin  != 0.0 ) {
     gmin = gmax ;
     }
   meang = meang + drift * meang ;
   rg = rg + drift * rg ;
   _lgetpar = mynormrand ( _threadargscomma_ meant / 1.0 , rt ) * 1.0 ;
   while ( _lgetpar <= 0.0 ) {
     _lgetpar = mynormrand ( _threadargscomma_ meant / 1.0 , rt ) * 1.0 ;
     }
   
return _lgetpar;
 }
 
static double _hoc_getpar(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  getpar ( _p, _ppvar, _thread, _nt );
 return(_r);
}
 
static void _net_receive (Point_process* _pnt, double* _args, double _lflag) 
{  double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _thread = (Datum*)0; _nt = (NrnThread*)_pnt->_vnt;   _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t;   if (_lflag == 1. ) {*(_tqitem) = nullptr;}
 {
   double _le ;
 if ( _lflag  == 1.0 ) {
     _le = getpar ( _threadargs_ ) ;
     net_send ( _tqitem, _args, _pnt, t +  _le , 1.0 ) ;
     }
   } }
 
/*VERBATIM*/
#ifndef NRN_VERSION_GTEQ_8_2_0
double nrn_random_pick(void* r); 
void* nrn_random_arg(int argpos);
#define RANDCAST
#else
#define RANDCAST (Rand*)
#endif
 
double mynormrand ( _threadargsprotocomma_ double _lmean , double _lvar ) {
   double _lmynormrand;
 
/*VERBATIM*/
	if (_p_donotuse) {
		double x = nrn_random_pick(RANDCAST _p_donotuse);
		_lmynormrand = x*_lvar + _lmean;
	}else{
		_lmynormrand = _lmean;
	}
 
return _lmynormrand;
 }
 
static double _hoc_mynormrand(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  mynormrand ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) );
 return(_r);
}
 
static int  setRandom ( _threadargsproto_ ) {
   
/*VERBATIM*/
 {
        void** pv = (void**)(&_p_donotuse);
        if (ifarg(1)) {
                *pv = nrn_random_arg(1);
        }else{
                *pv = (void*)0;
        }
 }
  return 0; }
 
static double _hoc_setRandom(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 setRandom ( _p, _ppvar, _thread, _nt );
 return(_r);
}

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
 {
   net_send ( _tqitem, nullptr, _ppvar[1].get<Point_process*>(), t +  event , 1.0 ) ;
   }

}
}

static void nrn_init(NrnThread* _nt, Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _tsav = -1e20;
#if EXTRACELLULAR
 _nd = _ml->_nodelist[_iml];
 if (_nd->_extnode) {
    _v = NODEV(_nd) +_nd->_extnode->_v[0];
 }else
#endif
 {
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 }
 v = _v;
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   double _lx ;
 if ( gmax  == gmin ) {
     g = gmax ;
     i = g * ( vgap - v ) * ( .001 ) ;
     }
   else {
     if ( id > 0.0 ) {
       _lx = v - vgap ;
       }
     else if ( id < 0.0 ) {
       _lx = vgap - v ;
       }
     else {
       
/*VERBATIM*/
			assert(0);
 }
     g = gv ( _threadargscomma_ _lx ) ;
     i = g * ( vgap - v ) * ( .001 ) ;
     }
   }
 _current += i;

} return _current;
}

static void nrn_cur(NrnThread* _nt, Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if EXTRACELLULAR
 _nd = _ml->_nodelist[_iml];
 if (_nd->_extnode) {
    _v = NODEV(_nd) +_nd->_extnode->_v[0];
 }else
#endif
 {
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 }
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) += _rhs;
  }else
#endif
  {
	NODERHS(_nd) += _rhs;
  }
  if (_nt->_nrn_fast_imem) { _nt->_nrn_fast_imem->_nrn_sav_rhs[_ni[_iml]] += _rhs; }
#if EXTRACELLULAR
 if (_nd->_extnode) {
   *_nd->_extnode->_rhs[0] += _rhs;
 }
#endif
 
}
 
}

static void nrn_jacob(NrnThread* _nt, Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) -= _g;
  }else
#endif
  {
	NODED(_nd) -= _g;
  }
  if (_nt->_nrn_fast_imem) { _nt->_nrn_fast_imem->_nrn_sav_d[_ni[_iml]] -= _g; }
#if EXTRACELLULAR
 if (_nd->_extnode) {
   *_nd->_extnode->_d[0] += _g;
 }
#endif
 
}
 
}

static void nrn_state(NrnThread* _nt, Memb_list* _ml, int _type) {

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mech_type) {
    const char* nmodl_filename = "halfgapm1.mod";
    const char* nmodl_file_text = 
  "NEURON {\n"
  "	POINT_PROCESS HalfGap\n"
  "	ELECTRODE_CURRENT i\n"
  "	RANGE g, i, meang, meant, rg, rt, drift\n"
  "	THREADSAFE : Only true if every instance has its own distinct Random\n"
  "	POINTER donotuse, vgap : A Normal Random generator with mean 1 and var 1.\n"
  "	RANGE id : For polarity of rectification and testing.\n"
  "		 : Should be equal and opposite for corresponding HalfGap\n"
  "		 : and otherwise distinct. For proper simulation results,\n"
  "		 : corresponding gaps should always have the same value\n"
  "		 : of g.\n"
  "	RANGE gmax, gmin, vhalf : Sigmoidal voltage sensitive conductance\n"
  "		: parameters. See gv(x) below. The sign of id defines\n"
  "		: the voltage polarity. If gmax == gmin, the gap is linear\n"
  "		: and id is not used.\n"
  "		: in pargap, gmin==gmax (linear) unless gmin is 0.\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	gmax = 1 (nanosiemens)\n"
  "	gmin = 1 (nanosiemens)\n"
  "	vhalf = 0 (millivolt)\n"
  "	slope4 = 10 (/millivolt)\n"
  "	meang = 30 (nanosiemens)\n"
  "	meant = 1000000 (ms)\n"
  "	drift = 0\n"
  "	rg=0\n"
  "	rt=0\n"
  "	event=0 (ms) : when gmax,gmin first assigned from meang,rg\n"
  "	id = 0\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	g (nanosiemens)\n"
  "	v (millivolt)\n"
  "	vgap (millivolt)\n"
  "	i (nanoamp)\n"
  "	donotuse\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	net_send(event,1)\n"
  "}\n"
  "\n"
  "\n"
  ": voltage sensitve gap conductance\n"
  ": for global variable time step, should be continuous to high order so\n"
  ": that performance does not suffer.\n"
  ": Argument is relative voltage at the positive polarity side.\n"
  "FUNCTION gv(x(millivolt))(nanosiemens) {\n"
  "	: sigmoid x >> vhalf means gv = gmax, x << vhalf means g = gmin\n"
  "	gv = (gmax - gmin)/(1 + exp(slope4*(vhalf - x))) + gmin\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	LOCAL x\n"
  "	if (gmax == gmin) { :linear gap junction\n"
  "		g = gmax\n"
  "		i = g * (vgap - v) * (.001)\n"
  "	}else{\n"
  "		: vgap > v means current is outward from this gap\n"
  "		if (id > 0 ) {\n"
  "			x = v - vgap :voltage relative to - side of gap\n"
  "		}else if (id < 0){\n"
  "			x = vgap - v : voltage relative to - side of gap\n"
  "		}else{\n"
  "VERBATIM\n"
  "			assert(0);\n"
  "ENDVERBATIM\n"
  "\n"
  "		}\n"
  "		g = gv(x)\n"
  "		i = g * (vgap - v) * (.001)\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION getpar() {\n"
  "	gmax=mynormrand(meang/1(nanosiemens),rg)*1(nanosiemens)\n"
  "	if (gmax<0) {gmax=0}\n"
  "	if (gmin != 0) {\n"
  "		gmin = gmax\n"
  "	}\n"
  "	meang=meang+drift*meang\n"
  "	rg=rg+drift*rg\n"
  "	getpar=mynormrand(meant/1(ms),rt)*1(ms)\n"
  "	WHILE(getpar <= 0) {\n"
  "		getpar = mynormrand(meant/1(ms), rt)*1(ms)\n"
  "	}\n"
  "}\n"
  "\n"
  "NET_RECEIVE (w) {\n"
  "	LOCAL e\n"
  "	if (flag == 1) { : from external\n"
  "		e = getpar()		:sets gmax,=gmin and next change\n"
  "		net_send(e, 1)\n"
  "	}\n"
  "}\n"
  "\n"
  ":Separate independent but reproducible streams for each instance.\n"
  ":For proper functioning, it is important that hoc Random distribution be\n"
  ": Random.Random123(id1, id2) <one could use MCellRan4 instead>\n"
  ": Random.normal(1,1)\n"
  ": and that corresponding HalfGap have the same id1, id2\n"
  ": A condition for correctness, that can be tested from hoc, is that\n"
  ": g (and also Random.seq()) for corresponding HalfGap have the same value.\n"
  ": If this is the case, then simulations with different numbers of processes\n"
  ": and different distibutions of gids should give quantitatively identical\n"
  ": results with the fixed step method and  (if cvode.use_long_double(1))\n"
  ": with the global variable time step method.\n"
  "\n"
  "VERBATIM\n"
  "#ifndef NRN_VERSION_GTEQ_8_2_0\n"
  "double nrn_random_pick(void* r); \n"
  "void* nrn_random_arg(int argpos);\n"
  "#define RANDCAST\n"
  "#else\n"
  "#define RANDCAST (Rand*)\n"
  "#endif\n"
  "ENDVERBATIM\n"
  "\n"
  "\n"
  "FUNCTION mynormrand(mean, var) {\n"
  "VERBATIM\n"
  "	if (_p_donotuse) {\n"
  "		double x = nrn_random_pick(RANDCAST _p_donotuse);\n"
  "		_lmynormrand = x*_lvar + _lmean;\n"
  "	}else{\n"
  "		_lmynormrand = _lmean;\n"
  "	}\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  "PROCEDURE setRandom() {\n"
  "VERBATIM\n"
  " {\n"
  "        void** pv = (void**)(&_p_donotuse);\n"
  "        if (ifarg(1)) {\n"
  "                *pv = nrn_random_arg(1);\n"
  "        }else{\n"
  "                *pv = (void*)0;\n"
  "        }\n"
  " }\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  ;
    hoc_reg_nmodl_filename(mech_type, nmodl_filename);
    hoc_reg_nmodl_text(mech_type, nmodl_file_text);
}
#endif
