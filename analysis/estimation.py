"""Paired reference, unrestricted rearrangement and compatible-release comparisons."""
from pathlib import Path
import math
import sys
import time
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import observation_model as reference
from method import scores, release

def evaluate(p, d, tau, scale=1.0, common=False, check_optimality=False):
    """Evaluate the original marginal center and empirical variance comparators in study units."""
    if isinstance(p, np.ndarray):
        p = p.tolist()
    if isinstance(d, np.ndarray):
        d = d.tolist()
    if len(p) != len(d):
        raise ValueError('unequal patient universe')
    beg = time.perf_counter()
    mp, up, vp = scores(p, tau)
    md, ud, vd = scores(d, tau)
    if not np.all(reference.dense_graph(p, d, 'separate').diagonal()):
        raise AssertionError('actual separate pairing infeasible')
    if common and (not np.all(reference.dense_graph(p, d, 'common').diagonal())):
        raise AssertionError('actual common pairing infeasible')
    tol = 1e-10 * (tau * scale) ** 2
    for u, v in [(up, vp), (ud, vd)]:
        if abs(float(u @ u) - v) * scale ** 2 > max(tol, 1e-09 * abs(v * scale ** 2)):
            raise AssertionError('Greenwood/score identity')
    oracle = float(np.sum((up - ud) ** 2))
    su = np.sort(up)
    sd = np.sort(ud)
    unlo = float(np.sum((su - sd) ** 2))
    unhi = float(np.sum((su - sd[::-1]) ** 2))
    result = dict(
        rmst_P=mp * scale,
        rmst_D=md * scale,
        point=(md - mp) * scale,
        var_P=vp * scale ** 2,
        var_D=vd * scale ** 2,
        covariance=float(up @ ud) * scale ** 2,
        oracle=oracle * scale ** 2,
        independence=(vp + vd) * scale ** 2,
        Cauchy_lower=(math.sqrt(vp) - math.sqrt(vd)) ** 2 * scale ** 2,
        Cauchy_upper=(math.sqrt(vp) + math.sqrt(vd)) ** 2 * scale ** 2,
        unrestricted_lower=unlo * scale ** 2,
        unrestricted_upper=unhi * scale ** 2,
        levels={},
        seconds=None,
        ratio_oracle_defined=oracle > 0,
        ratio_unrestricted_defined=unhi > 0)
    levels = result['levels']
    for model in ['separate', 'common'] if common else ['separate']:
        for level in ['R1', 'R2']:
            b = release(p, d, tau, model, level, check_optimality)
            lo = b['lower'] * scale ** 2
            hi = b['upper'] * scale ** 2
            if lo > result['oracle'] + tol or hi < result['oracle'] - tol:
                raise AssertionError('oracle not included')
            if lo < result['unrestricted_lower'] - tol or hi > result['unrestricted_upper'] + tol:
                raise AssertionError('unrestricted nesting')
            b['lower_study_units'] = lo
            b['upper_study_units'] = hi
            b['covariance_lower_study_units'] = (result['var_P'] + result['var_D'] - hi) / 2
            b['covariance_upper_study_units'] = (result['var_P'] + result['var_D'] - lo) / 2
            levels[model + '_' + level] = b
        a = levels[model + '_R1']
        b = levels[model + '_R2']
        if a['lower_study_units'] < b['lower_study_units'] - tol or a['upper_study_units'] > b['upper_study_units'] + tol:
            raise AssertionError('R1/R2 nesting')
    if common:
        for level in ['R1', 'R2']:
            a = levels['common_' + level]
            b = levels['separate_' + level]
            if a['lower_study_units'] < b['lower_study_units'] - tol or a['upper_study_units'] > b['upper_study_units'] + tol:
                raise AssertionError('common/separate nesting')
    result['seconds'] = time.perf_counter() - beg
    return result
