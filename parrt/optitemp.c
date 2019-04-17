#include "petsc.h"
#include "nlopt.h"
#include "optitemp.h"

double optimize_generic(int DegFree, double *epsopt,
			double *lb, double *ub,
			nlopt_func obj, void *objdata,
			nlopt_func *constraint, void **constrdata, int nconstraints,
			alg_ alg,
			nlopt_result *result)
{

  nlopt_opt opt;
  nlopt_opt local_opt;
  int i;
  double maxf;

  opt = nlopt_create(alg.outer, DegFree);
  nlopt_set_lower_bounds(opt,lb);
  nlopt_set_upper_bounds(opt,ub);
  nlopt_set_maxeval(opt,alg.maxeval);
  nlopt_set_maxtime(opt,alg.maxtime);

  //if(alg.outer==11) nlopt_set_vector_storage(opt,4000);
  if(alg.inner){
    local_opt=nlopt_create(alg.inner, DegFree);
    nlopt_set_ftol_rel(local_opt, 1e-11);
    nlopt_set_maxeval(local_opt,10000);
    nlopt_set_local_optimizer(opt,local_opt);
  }

  if(nconstraints){
    for(i=0;i<nconstraints;i++){
      nlopt_add_inequality_constraint(opt,constraint[i],constrdata[i],1e-8);
    }
  }

  if(obj){
    if(alg.maxobj)
      nlopt_set_max_objective(opt,obj,objdata);
    else
      nlopt_set_min_objective(opt,obj,objdata);
    *result=nlopt_optimize(opt,epsopt,&maxf);
  }

  nlopt_destroy(opt);
  if(alg.inner) nlopt_destroy(local_opt);

  return maxf;
  
}
