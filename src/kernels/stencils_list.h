
#ifndef KA_STR
#define KA_STR(x)   #x
#define KA_SHOW(x)   printf("%s=%s\n", #x, KA_STR(x))
#endif
extern Parameters *gp; //@KADIR global parameter within a node
#include "stencils_spt_blk.ic"
#include "stencils_spt_blk_stat_sched.ic"
#include "stencils_1wf.ic"
#include "stencils_mwf.ic"
#include "stencils_femwf.ic"
#include "stencils_rsmwf.ic"
#include "stencils_rsfemwf.ic"
