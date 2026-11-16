"""Explicit local Clarabel solver adapter for native MICOM/optlang hybrid.

Changes only QP backend dispatch, not MICOM's algorithms or model constraints.
The exact optlang quadratic matrix is transformed into equality/nonnegative
cones; full feasibility, KKT and complementarity diagnostics are retained.
"""
import numpy as np
import scipy.sparse as ss
import clarabel
from optlang.hybrid_interface import HybridProblem


def solve_clarabel(self):
    sp=self.build(); nv=len(self.variables); nc=len(self.constraints)
    a=ss.vstack([sp.A,ss.eye(nv,format='csc')],format='csc')
    bounds=np.vstack([sp.bounds,sp.vbounds]); lb,ub=bounds.T
    eq=np.flatnonzero(np.isfinite(lb)&np.isfinite(ub)&(lb==ub))
    hi=np.flatnonzero(np.isfinite(ub)&(lb!=ub))
    lo=np.flatnonzero(np.isfinite(lb)&(lb!=ub))
    aa=ss.vstack([a[eq],a[hi],-a[lo]],format='csc')
    b=np.r_[ub[eq],ub[hi],-lb[lo]]
    cones=[clarabel.ZeroConeT(len(eq)),clarabel.NonnegativeConeT(len(hi)+len(lo))]
    settings=clarabel.DefaultSettings()
    settings.verbose=bool(self.settings['verbose'])
    settings.max_iter=500
    settings.time_limit=float(self.settings['time_limit'] or 120)
    settings.tol_feas=1e-9
    settings.tol_gap_abs=1e-9
    settings.tol_gap_rel=1e-9
    settings.max_threads=1
    for name,value in getattr(self,'clarabel_options',{}).items():
        setattr(settings,name,value)
    solver=clarabel.DefaultSolver(ss.triu(sp.P,format='csc'),sp.q,aa,b,cones,settings)
    sol=solver.solve();x=np.asarray(sol.x);z=np.asarray(sol.z);slack=np.asarray(sol.s)
    status=str(sol.status)
    self.status={'Solved':'Optimal','PrimalInfeasible':'Infeasible','DualInfeasible':'Unbounded','MaxIterations':'Iteration limit reached','MaxTime':'Time limit reached'}.get(status,'Unknown')
    rd=np.zeros(nc+nv); n1=len(eq);n2=n1+len(hi)
    rd[eq]=-z[:n1];rd[hi]-=z[n1:n2];rd[lo]+=z[n2:]
    d=float(self.direction)
    self.primals=dict(zip(self.variables,x))
    self.cprimals=dict(zip(self.constraints,sp.A@x))
    self.duals=dict(zip(self.constraints,rd[:nc]*d))
    self.vduals=dict(zip(self.variables,rd[nc:]*d))
    self.obj_value=float(sol.obj_val)*d
    self.info=sol;self._solution=None
    ax=a@x
    cert={'status':status,'objective':float(sol.obj_val),
          'max_primal_violation':float(max(np.maximum(lb-ax,0).max(),np.maximum(ax-ub,0).max())),
          'max_stationarity_residual':float(np.abs(sp.P@x+sp.q+aa.T@z).max()),
          'max_complementarity_residual':float(np.abs(z[n1:]*slack[n1:]).max()),
          'dual_nonnegative_violation':float(np.maximum(-z[n1:],0).max()),
          'relative_primal_residual':sol.r_prim,'relative_dual_residual':sol.r_dual,
          'dual_objective':sol.obj_val_dual,'iterations':sol.iterations,'elapsed_seconds':sol.solve_time,
          'variables':dict(zip(self.variables,map(float,x)))}
    self.qp_certificates=getattr(self,'qp_certificates',[])+[cert]


def install():
    HybridProblem.solve_osqp=solve_clarabel


if __name__=='__main__':
    from optlang.hybrid_interface import Model,Variable,Objective,Constraint
    install()
    for direction,sign in [('min',1),('max',-1)]:
        m=Model();x=Variable('x',lb=0);y=Variable('y',lb=0)
        m.add([x,y,Constraint(x+y,lb=3)])
        m.objective=Objective(sign*(x*x+2*x*y+2*y*y-2*x-6*y),direction=direction)
        assert m.optimize()=='optimal'
        assert abs(x.primal-1)<1e-5 and abs(y.primal-2)<1e-5
        c=m.problem.qp_certificates[-1]
        assert c['max_primal_violation']<1e-8 and c['max_stationarity_residual']<1e-8
    print('Clarabel adapter analytic cross-term and objective-direction checks passed')
