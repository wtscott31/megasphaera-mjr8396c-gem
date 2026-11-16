"""Exact elimination of fixed variables before a native sparse QP solve."""
import numpy as np
from optlang.hybrid_interface import HybridProblem

def solve_reduced(self,solver):
    sp=self.build();names=list(self.variables);cnames=list(self.constraints)
    fixed={k:self.variable_lbs[k] for k in names if self.variable_lbs[k]==self.variable_ubs[k]}
    free=self.variables-set(fixed)
    p=HybridProblem();p.settings=dict(self.settings);p.direction=self.direction
    p.variables=free;p.variable_lbs={k:self.variable_lbs[k] for k in free};p.variable_ubs={k:self.variable_ubs[k] for k in free}
    shift={k:0.0 for k in cnames}
    for (cn,vn),value in self.constraint_coefs.items():
        if vn in fixed:shift[cn]+=value*fixed[vn]
    nonzero={cn for (cn,vn),coef in self.constraint_coefs.items() if vn in free and coef!=0}
    p.constraints=nonzero
    for k in set(cnames)-nonzero:
        if not self.constraint_lbs[k]-1e-8<=shift[k]<=self.constraint_ubs[k]+1e-8:raise ValueError('Inconsistent fixed-variable row '+k)
    p.constraint_lbs={k:self.constraint_lbs[k]-shift[k] for k in nonzero};p.constraint_ubs={k:self.constraint_ubs[k]-shift[k] for k in nonzero}
    p.constraint_coefs={(cn,vn):v for (cn,vn),v in self.constraint_coefs.items() if vn in free and cn in nonzero}
    p.obj_linear_coefs={k:v for k,v in self.obj_linear_coefs.items() if k in free}
    p.obj_quadratic_coefs={(a,b):v for (a,b),v in self.obj_quadratic_coefs.items() if a in free and b in free}
    for (a,b),v in self.obj_quadratic_coefs.items():
        if a in free and b in fixed:p.obj_linear_coefs[a]=p.obj_linear_coefs.get(a,0)+v*fixed[b]
        if b in free and a in fixed:p.obj_linear_coefs[b]=p.obj_linear_coefs.get(b,0)+v*fixed[a]
    solver(p)
    self.status=p.status;self.primals=fixed|p.primals
    x=np.array([self.primals[n] for n in names]);ax=sp.A@x;gradient=sp.P@x+sp.q;d=float(self.direction)
    rd=np.array([p.duals.get(n,0)*d for n in cnames]);vd=np.array([p.vduals.get(n,0)*d for n in names])
    fixed_indices=np.array([i for i,n in enumerate(names) if n in fixed],dtype=int)
    vd[fixed_indices]=(gradient-sp.A.T@rd)[fixed_indices]
    self.cprimals=dict(zip(cnames,ax));self.duals=dict(zip(cnames,rd*d));self.vduals=dict(zip(names,vd*d))
    obj=float(.5*x@(sp.P@x)+sp.q@x);self.obj_value=obj*d;self.info=p.info;self._solution=None
    cert=p.qp_certificates[-1];cert['objective']=obj;cert['variables']=dict(zip(names,map(float,x)))
    bounds=np.vstack([sp.bounds,sp.vbounds]);values=np.r_[ax,x];dual=np.r_[rd,vd]
    active=(bounds[:,0]!=bounds[:,1])&(np.abs(dual)>1e-12)
    slack=np.where(dual>=0,values-bounds[:,0],bounds[:,1]-values)
    cert.update(exact_fixed_variables_eliminated=len(fixed),empty_rows_eliminated=len(cnames)-len(nonzero),
        max_primal_violation=float(max(np.maximum(bounds[:,0]-values,0).max(),np.maximum(values-bounds[:,1],0).max())),
        max_stationarity_residual=float(np.abs(gradient-sp.A.T@rd-vd).max()),
        max_complementarity_residual=float(np.abs(dual[active]*slack[active]).max()))
    self.qp_certificates=getattr(self,'qp_certificates',[])+[cert]
