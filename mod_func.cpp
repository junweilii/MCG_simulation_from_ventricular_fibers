#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern "C" void _cadynam_reg();
extern "C" void _gap_reg();
extern "C" void _halfgapm1_reg();
extern "C" void _ik1_reg();
extern "C" void _ina_reg();
extern "C" void _is_reg();
extern "C" void _ix1_reg();
extern "C" void _kcum_reg();
extern "C" void _nacum_reg();
extern "C" void _nipace_reg();

extern "C" void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," cadynam.mod");
fprintf(stderr," gap.mod");
fprintf(stderr," halfgapm1.mod");
fprintf(stderr," ik1.mod");
fprintf(stderr," ina.mod");
fprintf(stderr," is.mod");
fprintf(stderr," ix1.mod");
fprintf(stderr," kcum.mod");
fprintf(stderr," nacum.mod");
fprintf(stderr," nipace.mod");
fprintf(stderr, "\n");
    }
_cadynam_reg();
_gap_reg();
_halfgapm1_reg();
_ik1_reg();
_ina_reg();
_is_reg();
_ix1_reg();
_kcum_reg();
_nacum_reg();
_nipace_reg();
}
