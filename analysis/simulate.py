"""Study sampling definitions and exact represented-value lattice conversion."""
from fractions import Fraction
import math
import numpy as np

def ceil_ticks(values, den):
    """Ceiling of the represented value times the lattice denominator, without an epsilon."""
    if den not in (4, 10):
        raise ValueError('unsupported lattice spacing')
    if isinstance(values, Fraction):
        if values < 0:
            raise ValueError('negative latent time')
        z = values * den
        return -(-z.numerator // z.denominator)
    raw = np.asarray(values)
    if raw.dtype == object:
        return np.array([ceil_ticks(x, den) for x in raw.ravel()], dtype=np.int64).reshape(raw.shape)
    a = np.asarray(values, float)
    if np.any(~np.isfinite(a)) or np.any(a < 0):
        raise ValueError('invalid latent times')
    product = a * den
    if np.any(product >= 2 ** 62):
        raise ValueError('tick overflow')
    ticks = np.asarray(np.ceil(product), dtype=np.int64)
    af = a.ravel()
    pf = product.ravel()
    tf = ticks.ravel()
    for i in np.flatnonzero(pf == np.floor(pf)):
        num, d = af[i].as_integer_ratio()
        tf[i] = -(-num * den // d)
    return ticks


def observe(P, D, CP, CD):
    """Separate endpoint censoring with events preceding censors at equal times."""
    if np.any(np.asarray(P) > D):
        raise ValueError('latent P>D')
    return (np.c_[np.minimum(P, CP), np.asarray(P) <= CP], np.c_[np.minimum(D, CD), np.asarray(D) <= CD])


def generate(cell, rep, seed, rng=None):
    """Generate one prespecified cell/key. Lattice times are discretized before observation."""
    if rng is None:
        rng = np.random.default_rng(np.random.SeedSequence([seed, cell['id'], rep, 0]))
    n = cell['n']
    kind = cell['kind']
    scale = 1.0
    tau = cell['tau']
    if kind == 'FGM':
        P = rng.random(n)
        z = rng.random(n)
        k = cell['a'] * (1 - 2 * P)
        D = 2 + 2 * z / (1 + k + np.sqrt((1 + k) ** 2 - 4 * k * z))
        return (np.c_[P, np.ones(n)], np.c_[D, np.ones(n)], tau, scale)
    if kind in ('semi_markov_exponential', 'semi_markov_exponential_latent_lattice'):
        p0 = rng.exponential(1 / cell['lambda_P'], n)
        u = rng.random(n)
        w = rng.exponential(1 / cell['lambda_W'], n)
    elif kind == 'semi_markov_weibull':
        from scipy.special import gamma
        kp = cell['weibull_shape_P']
        kw = cell['weibull_shape_W']
        p0 = rng.weibull(kp, n) * cell['mean_P'] / gamma(1 + 1 / kp)
        u = rng.random(n)
        w = rng.weibull(kw, n) * cell['mean_W'] / gamma(1 + 1 / kw)
    else:
        raise ValueError('unrecognized sampling distribution')
    P = p0
    D = np.where(u < cell['death_first_probability'], p0, p0 + w)
    rate = cell.get('censor_rate', cell.get('base_censor_rate'))
    CP = rng.exponential(1 / rate, n)
    CD = CP.copy() if cell['censor_model'] == 'common' else rng.exponential(1 / rate, n)
    if cell['censor_model'] == 'informative_separate':
        f = np.exp(-cell['gamma'] * (p0 - 1))
        CP *= f
        CD *= f
    if kind == 'semi_markov_exponential_latent_lattice':
        lattice = cell['latent_lattice']
        den = lattice['spacing_denominator']
        P, D, CP, CD = [ceil_ticks(x, den) for x in (P, D, CP, CD)]
        tau = lattice['tau_ticks']
        scale = lattice['spacing_numerator'] / den
    if 'round_observed_times_to' in cell:
        raise ValueError('observed-time rounding is not part of this sampling definition')
    p, d = observe(P, D, CP, CD)
    return (p, d, tau, scale)
