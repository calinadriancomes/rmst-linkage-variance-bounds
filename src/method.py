"""Attained empirical variance endpoints for exact marginal releases."""
import math
import numpy as np
from scipy.optimize import linear_sum_assignment
import observation_model as reference

def scores(records, tau):
    reference.validate(records, tau)
    r = np.asarray(records)
    t = r[:, 0]
    e = r[:, 1].astype(int)
    n = len(t)
    alltimes, codes = np.unique(t, return_inverse=True)
    count = np.bincount(codes, minlength=len(alltimes))
    deaths = np.bincount(codes, weights=e, minlength=len(alltimes))
    risk = np.cumsum(count[::-1])[::-1]
    keep = (deaths > 0) & (alltimes < tau)
    ts = alltimes[keep]
    y = risk[keep]
    d = deaths[keep]
    if not len(ts):
        return (float(tau), np.zeros(n), 0.0)
    S = np.cumprod(1 - d / y)
    pieces = np.diff(np.r_[ts, tau]) * S
    A = np.cumsum(pieces[::-1])[::-1]
    good = y > d
    positive = np.zeros(len(ts))
    negative = np.zeros(len(ts))
    positive[good] = A[good] * d[good] / (y[good] * (y[good] - d[good]))
    negative[good] = A[good] / (y[good] - d[good])
    j = np.searchsorted(ts, t, side='right') - 1
    u = np.zeros(n)
    ok = j >= 0
    u[ok] = np.cumsum(positive)[j[ok]]
    event = ok & (e == 1) & (t < tau)
    u[event] -= negative[j[event]]
    green = float(np.sum(A[good] ** 2 * d[good] / (y[good] * (y[good] - d[good]))))
    return (float(ts[0] + pieces.sum()), u, green)

def optimality_check(cost, graph, pi, sign):
    n = len(pi)
    c = np.where(graph, sign * cost, np.inf)
    matched = c[np.arange(n), pi]
    delta = c - matched[:, None]
    beta = np.zeros(n)
    for k in range(n + 1):
        updated = np.minimum(beta, np.min(beta[pi][:, None] + delta, axis=0))
        if np.max(np.abs(updated - beta)) < 1e-13 * max(1, float(np.max(cost))):
            beta = updated
            break
        beta = updated
    alpha = matched - beta[pi]
    value = float(cost[np.arange(n), pi].sum())
    violation = float(np.max((alpha[:, None] + beta[None, :] - c)[graph]))
    gap = float(abs(alpha.sum() + beta.sum() - sign * value))
    tol = 1e-09 * abs(value) if abs(value) > 1e-12 else 1e-12
    initial_violation = violation
    initial_gap = gap
    refinement_iterations = 0
    if violation > tol or gap > tol:
        for j in range(n + 1):
            updated = np.minimum(beta, np.min(beta[pi][:, None] + delta, axis=0))
            stable = np.array_equal(updated, beta)
            beta = updated
            refinement_iterations = j + 1
            if stable:
                break
        alpha = matched - beta[pi]
        violation = float(np.max((alpha[:, None] + beta[None, :] - c)[graph]))
        gap = float(abs(alpha.sum() + beta.sum() - sign * value))
    if violation > tol or gap > tol:
        raise AssertionError('optimality residual check failed')
    return dict(
        alpha=alpha.tolist(),
        beta=beta.tolist(),
        signed_objective=sign * value,
        gap=gap,
        max_dual_violation=violation,
        tolerance=tol,
        iterations=k + 1,
        refinement_iterations=refinement_iterations,
        initial_violation=initial_violation,
        initial_gap=initial_gap,
        residual_objective_allowance=gap + n * max(0.0, violation),
        adjusted_dual_objective=float(alpha.sum() + beta.sum() - n * max(0.0, violation)),
        numerical_status='residual_checked_after_refinement' if refinement_iterations else 'residual_checked')

def solve(up, ud, graph, check_optimality=False):
    cost = (up[:, None] - ud[None, :]) ** 2
    n = len(up)
    out = {}
    for label, sign in [('lower', 1), ('upper', -1)]:
        row, col = linear_sum_assignment(np.where(graph, sign * cost, np.inf))
        if not np.all(graph[row, col]):
            raise AssertionError('forbidden edge in assignment')
        out[label] = float(cost[row, col].sum())
        out[label + '_witness'] = col.tolist()
        if check_optimality:
            out[label + '_optimality_check'] = optimality_check(cost, graph, col, sign)
    out['edge_count'] = int(graph.sum())
    return out

def release(p, d, tau, model, level, check_optimality=False):
    if model not in ('separate', 'common'):
        raise ValueError("model must be 'separate' or 'common'")
    if len(p) != len(d):
        raise ValueError('unequal patient universe')
    if level == 'R1':
        pp, dd = (p, d)
        pc = dc = None
    elif level == 'R2':
        pc = reference.expand(reference.cells_from_table(
            len(p),
            reference.event_table(p),
            support=(tau, sum((t >= tau for t, e in p)))))
        dc = reference.expand(reference.cells_from_table(
            len(d),
            reference.event_table(d),
            support=(tau, sum((t >= tau for t, e in d)))))
        pp = [(l, a) for a, l, h in pc]
        dd = [(l, a) for a, l, h in dc]
    else:
        raise ValueError('unknown release')
    _, up, _ = scores(pp, tau)
    _, ud, _ = scores(dd, tau)
    graph = reference.dense_graph(
        p,
        d,
        model) if level == 'R1' else reference.dense_graph(pc,
        dc,
        model,
        'R2')
    out = solve(up, ud, graph, check_optimality)
    if check_optimality:
        out['P_scores'] = up.tolist()
        out['D_scores'] = ud.tolist()
        if level == 'R2':
            out['censor_cells_P'] = [(a, l, None if math.isinf(h) else h) for a, l, h in pc]
            out['censor_cells_D'] = [(a, l, None if math.isinf(h) else h) for a, l, h in dc]
            for label in ['lower', 'upper']:
                P, D = reference.instantiate(pc, dc, out[label + '_witness'], model)
                assert reference.event_table(P) == reference.event_table(p) and reference.event_table(D) == reference.event_table(d)
                assert sum((t >= tau for t, e in P)) == sum((t >= tau for t, e in p))
                assert sum((t >= tau for t, e in D)) == sum((t >= tau for t, e in d))
                assert np.max(np.abs(scores(P, tau)[1] - up)) < 1e-10 * max(1, tau)
                assert np.max(np.abs(scores(D, tau)[1] - ud)) < 1e-10 * max(1, tau)
                out[label + '_P_records'] = P
                out[label + '_D_records'] = D
    return out


def from_records(p, d, tau, model="separate", check_optimality=True):
    """R1: full exact marginal records on the same declared cohort."""
    result = release(p, d, tau, model, "R1", check_optimality)
    mp, up, vp = scores(p, tau)
    md, ud, vd = scores(d, tau)
    result.update(point=md-mp, rmst_P=mp, rmst_D=md, var_P=vp, var_D=vd)
    return result

def from_tables(n, table_P, table_D, tau, risk_P_tau, risk_D_tau, model="separate", check_optimality=True):
    """R2: full event tables (time, pre-event risk, multiplicity), including after tau.
    Horizon counts are pre-event risks. No cross-endpoint identifiers are inputs.
    """
    if model not in ("separate", "common"):
        raise ValueError("model must be 'separate' or 'common'")
    pc = reference.expand(reference.cells_from_table(n, table_P, support=(tau, risk_P_tau)))
    dc = reference.expand(reference.cells_from_table(n, table_D, support=(tau, risk_D_tau)))
    p = [(l,a) for a,l,h in pc]
    d = [(l,a) for a,l,h in dc]
    mp, up, vp = scores(p, tau)
    md, ud, vd = scores(d, tau)
    graph = reference.dense_graph(pc, dc, model, "R2")
    result = solve(up, ud, graph, check_optimality)
    result.update(point=md-mp, rmst_P=mp, rmst_D=md, var_P=vp, var_D=vd)
    return result

def upper_envelope_interval(point, variance_upper, normal_quantile=1.959963984540054):
    """Same-center envelope of compatible paired Wald intervals.
    Containment alone does not provide finite-sample nominal coverage.
    """
    if not math.isfinite(point) or not math.isfinite(variance_upper) or variance_upper < 0:
        raise ValueError("invalid interval inputs")
    if not math.isfinite(normal_quantile) or normal_quantile <= 0:
        raise ValueError("normal_quantile must be finite and positive")
    half = normal_quantile * math.sqrt(variance_upper)
    return point-half, point+half
