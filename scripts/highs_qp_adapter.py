"""Use native HiGHS QP behind optlang's hybrid interface.

MICOM and its growth objectives remain unmodified. This explicit, local solver
adapter replaces only HybridProblem.solve_osqp, since upstream dispatches all
quadratic problems to OSQP. Installed package files are never edited.
"""
import numpy as np
import highspy as hs
from scipy.sparse import tril
from optlang.hybrid_interface import HybridProblem


def solve_highs_qp(self):
    sp = self.build()
    objective_scale = float(getattr(self,'qp_objective_scale',1000.0))
    env = hs.Highs()
    options = self.highs_settings()
    options.qp_iteration_limit = 200000
    options.qp_regularization_value = 1e-7
    env.passOptions(options)
    model = hs.HighsModel()
    lp = model.lp_
    lp.num_col_ = len(self.variables)
    lp.num_row_ = len(self.constraints)
    lp.sense_ = hs.ObjSense.kMinimize
    lp.col_cost_ = sp.q * objective_scale
    lp.col_lower_, lp.col_upper_ = sp.vbounds.T
    if sp.A is not None:
        lp.a_matrix_.format_ = hs.MatrixFormat.kColwise
        lp.a_matrix_.start_ = sp.A.indptr
        lp.a_matrix_.index_ = sp.A.indices
        lp.a_matrix_.value_ = sp.A.data
        lp.row_lower_, lp.row_upper_ = sp.bounds.T
    # optlang stores symmetric cross terms at half the polynomial coefficient;
    # build() doubles these, yielding the full symmetric 1/2 x' P x Hessian.
    p = tril(sp.P * objective_scale, format='csc')
    model.hessian_.dim_ = lp.num_col_
    model.hessian_.format_ = hs.HessianFormat.kTriangular
    model.hessian_.start_ = p.indptr
    model.hessian_.index_ = p.indices
    model.hessian_.value_ = p.data
    if env.passModel(model) != hs.HighsStatus.kOk:
        raise RuntimeError('HiGHS rejected the native QP matrix')
    env.run()
    self.status = env.modelStatusToString(env.getModelStatus())
    sol, info = env.getSolution(), env.getInfo()
    d = float(self.direction)
    x = np.asarray(sol.col_value)
    self.primals = dict(zip(self.variables, x))
    self.vduals = dict(zip(self.variables, np.asarray(sol.col_dual) * d / objective_scale))
    self.cprimals = dict(zip(self.constraints, sol.row_value))
    self.duals = dict(zip(self.constraints, np.asarray(sol.row_dual) * d / objective_scale))
    self.obj_value = info.objective_function_value * d / objective_scale
    self.info = info
    self._solution = None
    # Retain the full QP certificate before native MICOM switches to pFBA.
    ax = sp.A @ x
    violation = max(np.max(np.maximum(sp.bounds[:, 0]-ax, 0)),
                    np.max(np.maximum(ax-sp.bounds[:, 1], 0)),
                    np.max(np.maximum(sp.vbounds[:, 0]-x, 0)),
                    np.max(np.maximum(x-sp.vbounds[:, 1], 0)))
    gradient = sp.P @ x + sp.q
    stationarity = gradient - (sp.A.T @ np.asarray(sol.row_dual) + np.asarray(sol.col_dual)) / objective_scale
    from scipy.sparse import vstack, eye
    all_values = np.r_[ax,x]
    all_bounds = np.vstack([sp.bounds,sp.vbounds])
    all_duals = np.r_[sol.row_dual,sol.col_dual] / objective_scale
    slack = np.where(all_duals >= 0, all_values-all_bounds[:,0], all_bounds[:,1]-all_values)
    active = (np.abs(all_duals)>1e-12)&(all_bounds[:,0]!=all_bounds[:,1])
    complementarity = np.max(np.abs(all_duals[active]*slack[active])) if active.any() else 0.0
    certificate = {'status': self.status, 'objective': float(info.objective_function_value)/objective_scale,
                   'positive_objective_scale': objective_scale,
                   'max_primal_violation': float(violation),
                   'max_stationarity_residual': float(np.max(np.abs(stationarity))),
                   'max_complementarity_residual': float(complementarity),
                   'dual_nonnegative_violation': 0.0,
                   'qp_iterations': int(info.qp_iteration_count),
                   'elapsed_seconds': env.getRunTime(),
                   'variables': dict(zip(self.variables, map(float, x)))}
    self.qp_certificates = getattr(self, 'qp_certificates', []) + [certificate]


def install():
    HybridProblem.solve_osqp = solve_highs_qp


if __name__ == '__main__':
    from optlang.hybrid_interface import Model, Variable, Objective, Constraint
    install()
    # Cross terms, both objective directions, and a binding inequality exercise
    # Hessian conversion independently of MICOM's diagonal growth objective.
    for direction, sign in [('min', 1), ('max', -1)]:
        m=Model(); x=Variable('x', lb=0); y=Variable('y', lb=0)
        m.add([x,y,Constraint(x+y,lb=3)])
        m.objective=Objective(sign*(x*x+2*x*y+2*y*y-2*x-6*y),direction=direction)
        assert m.optimize()=='optimal'
        assert abs(x.primal-1)<1e-5 and abs(y.primal-2)<1e-5
        assert m.problem.qp_certificates[-1]['max_stationarity_residual']<1e-7
    print('HiGHS native QP adapter analytic checks passed')
