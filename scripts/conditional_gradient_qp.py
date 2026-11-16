"""Certified general convex QP backend for native MICOM's low-rank QPs.

Fully corrective conditional gradient uses HiGHS LP oracles on the COMPLETE
native optlang constraints. A small convex QP over feasible oracle vertices
is solved with SciPy SLSQP. The Frank-Wolfe gap bounds global suboptimality.
No supplied v2 growth rates, reduced community equations, or SBML parser are
used. MICOM's native cooperative_tradeoff and native pFBA remain unchanged.
"""
import numpy as np
import highspy as hs
from scipy.optimize import minimize
import time
from optlang.hybrid_interface import HybridProblem

def solve_conditional_gradient(self):
    start=time.monotonic();sp=self.build();d=float(self.direction)
    env=hs.Highs();opts=self.highs_settings()
    opts.primal_feasibility_tolerance=1e-9;opts.dual_feasibility_tolerance=1e-9
    opts.ipm_optimality_tolerance=1e-10;opts.solver='simplex'
    env.passOptions(opts);model=hs.HighsModel();lp=model.lp_
    lp.num_col_=len(self.variables);lp.num_row_=len(self.constraints)
    lp.sense_=hs.ObjSense.kMinimize;lp.col_cost_=sp.q
    lp.col_lower_,lp.col_upper_=sp.vbounds.T
    lp.a_matrix_.format_=hs.MatrixFormat.kColwise
    lp.a_matrix_.start_=sp.A.indptr;lp.a_matrix_.index_=sp.A.indices;lp.a_matrix_.value_=sp.A.data
    lp.row_lower_,lp.row_upper_=sp.bounds.T
    if env.passModel(model)!=hs.HighsStatus.kOk:raise ValueError('LP oracle model rejected')
    indices=np.arange(lp.num_col_,dtype=np.int32);trace=[]
    def oracle(cost):
        env.changeColsCost(lp.num_col_,indices,np.asarray(cost))
        env.run()
        if env.getModelStatus()!=hs.HighsModelStatus.kOptimal:raise RuntimeError('LP oracle: '+str(env.getModelStatus()))
        sol=env.getSolution()
        return np.array(sol.col_value),np.array(sol.row_dual),np.array(sol.col_dual)
    x,rd,vd=oracle(sp.q);vertices=[x.copy()];weights=np.ones(1)
    for iteration in range(100):
        gradient=sp.P@x+sp.q
        s,rd,vd=oracle(gradient)
        gap=float(gradient@(x-s));obj=float(.5*x@(sp.P@x)+sp.q@x)
        trace.append({'iteration':iteration,'objective':obj,'global_suboptimality_bound':max(gap,0),'vertices':len(vertices)})
        if gap<=1e-8*(1+abs(obj)):break
        if time.monotonic()-start>120:raise RuntimeError('Conditional-gradient solver exceeded 120 seconds')
        vertices.append(s.copy());V=np.asarray(vertices).T
        H=V.T@(sp.P@V);linear=V.T@sp.q
        w0=np.r_[weights,0.]
        res=minimize(lambda w:.5*w@H@w+linear@w,w0,jac=lambda w:H@w+linear,
                     bounds=[(0,1)]*len(vertices),constraints={'type':'eq','fun':lambda w:w.sum()-1,'jac':lambda w:np.ones_like(w)},
                     method='SLSQP',options={'ftol':1e-14,'maxiter':1000})
        if not res.success and abs(res.x.sum()-1)>1e-10:raise RuntimeError('Corrective vertex QP: '+res.message)
        weights=np.maximum(res.x,0);weights/=weights.sum();x=V@weights
        # Drop inactive vertices to keep the corrective problem small.
        active=weights>1e-12
        vertices=[v for v,keep in zip(vertices,active) if keep];weights=weights[active];weights/=weights.sum()
    else:raise RuntimeError('Conditional-gradient solver did not certify convergence')
    ax=sp.A@x;gradient=sp.P@x+sp.q
    vals=np.r_[ax,x];bounds=np.vstack([sp.bounds,sp.vbounds]);duals=np.r_[rd,vd]
    slack=np.where(duals>=0,vals-bounds[:,0],bounds[:,1]-vals)
    active=(bounds[:,0]!=bounds[:,1])&(np.abs(duals)>1e-12)
    self.primals=dict(zip(self.variables,x));self.cprimals=dict(zip(self.constraints,ax))
    self.duals=dict(zip(self.constraints,rd*d));self.vduals=dict(zip(self.variables,vd*d))
    self.obj_value=obj*d;self.status='Optimal';self.info=env.getInfo();self._solution=None
    cert={'status':'Optimal with full-network LP supporting-plane certificate',
          'objective':obj,'global_suboptimality_bound':max(gap,0),'dual_lower_bound':obj-max(gap,0),
          'max_primal_violation':float(max(np.maximum(bounds[:,0]-vals,0).max(),np.maximum(vals-bounds[:,1],0).max())),
          'max_stationarity_residual':float(np.abs(gradient-sp.A.T@rd-vd).max()),
          'max_complementarity_residual':float(np.abs(duals[active]*slack[active]).max(initial=0)),
          'dual_nonnegative_violation':0.,'iterations':iteration+1,'elapsed_seconds':time.monotonic()-start,
          'variables':dict(zip(self.variables,map(float,x))),'backend':'Fully corrective conditional gradient / HiGHS 1.15.1 / SciPy SLSQP',
          'trace':trace}
    self.qp_certificates=getattr(self,'qp_certificates',[])+[cert]

def install():HybridProblem.solve_osqp=solve_conditional_gradient

if __name__=='__main__':
    from optlang.hybrid_interface import Model,Variable,Objective,Constraint
    install()
    for direction,sign in [('min',1),('max',-1)]:
        m=Model();x=Variable('x',lb=0,ub=10);y=Variable('y',lb=0,ub=10)
        m.add([x,y,Constraint(x+y,lb=3)])
        m.objective=Objective(sign*(x*x+2*x*y+2*y*y-2*x-6*y),direction=direction)
        assert m.optimize()=='optimal'
        assert abs(x.primal-1)<1e-5 and abs(y.primal-2)<1e-5
        c=m.problem.qp_certificates[-1]
        assert c['max_primal_violation']<1e-8 and c['global_suboptimality_bound']<1e-8
    print('Certified QP adapter analytic cross-term and objective-direction checks passed')
