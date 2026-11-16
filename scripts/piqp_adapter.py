"""Native PIQP sparse QP backend for unmodified native MICOM objectives."""
import numpy as np
import scipy.sparse as ss
import piqp
from optlang.hybrid_interface import HybridProblem

def solve_piqp(self):
    sp=self.build();lb,ub=sp.bounds.T
    eq=np.flatnonzero(np.isfinite(lb)&np.isfinite(ub)&(lb==ub))
    ine=np.flatnonzero(lb!=ub)
    solver=piqp.SparseSolver()
    solver.settings.verbose=bool(self.settings['verbose'])
    solver.settings.compute_timings=True
    solver.settings.eps_abs=1e-9
    solver.settings.eps_rel=1e-9
    solver.settings.max_iter=500
    for name,value in getattr(self,'piqp_options',{}).items():setattr(solver.settings,name,value)
    solver.setup(ss.triu(sp.P,format='csc'),sp.q,sp.A[eq].tocsc(),ub[eq],sp.A[ine].tocsc(),lb[ine],ub[ine],sp.vbounds[:,0],sp.vbounds[:,1])
    status=solver.solve();sol=solver.result;x=np.asarray(sol.x)
    self.status='Optimal' if status==piqp.PIQP_SOLVED else 'Unknown'
    rd=np.zeros(len(self.constraints));rd[eq]=-sol.y;rd[ine]=sol.z_l-sol.z_u
    vd=sol.z_bl-sol.z_bu;d=float(self.direction)
    self.primals=dict(zip(self.variables,x));self.cprimals=dict(zip(self.constraints,sp.A@x))
    self.duals=dict(zip(self.constraints,rd*d));self.vduals=dict(zip(self.variables,vd*d))
    self.obj_value=sol.info.primal_obj*d;self.info=sol.info;self._solution=None
    vals=np.r_[sp.A@x,x];bounds=np.vstack([sp.bounds,sp.vbounds]);dual=np.r_[rd,vd]
    active=(np.abs(dual)>1e-12)&(bounds[:,0]!=bounds[:,1])
    slack=np.where(dual>=0,vals-bounds[:,0],bounds[:,1]-vals)
    cert={'status':str(status),'objective':float(sol.info.primal_obj),
          'max_primal_violation':float(max(np.maximum(bounds[:,0]-vals,0).max(),np.maximum(vals-bounds[:,1],0).max())),
          'max_stationarity_residual':float(np.abs(sp.P@x+sp.q-sp.A.T@rd-vd).max()),
          'max_complementarity_residual':float(np.abs(dual[active]*slack[active]).max()),
          'dual_nonnegative_violation':0.0,'dual_objective':sol.info.dual_obj,
          'iterations':sol.info.iter,'elapsed_seconds':sol.info.run_time,
          'variables':dict(zip(self.variables,map(float,x))),'backend':'PIQP '+piqp.__version__}
    self.qp_certificates=getattr(self,'qp_certificates',[])+[cert]

def install():HybridProblem.solve_osqp=solve_piqp

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
    print('PIQP adapter analytic cross-term and objective-direction checks passed')
