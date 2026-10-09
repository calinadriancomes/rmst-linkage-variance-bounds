"""Marginal case-weight derivatives and exact-release compatibility."""
from fractions import Fraction as Q
from itertools import permutations
import math
import numpy as np

def validate(records, tau, allow_unidentified=False):
    if len(records) == 0 or not math.isfinite(float(tau)) or tau < 0:
        raise ValueError('invalid size/horizon')
    if any((not math.isfinite(float(t)) or t < 0 or d not in (0, 1) for t, d in records)):
        raise ValueError('invalid observations')
    maximum = max((t for t, d in records))
    if maximum < tau and (not all((d == 1 for t, d in records if t == maximum))):
        if not allow_unidentified:
            raise ValueError('unsupported positive-survival tail')

def exact_scores(records, tau, check_support=True):
    if check_support:
        validate(records, tau)
    n = len(records)
    tau = Q(tau)
    times = sorted({Q(t) for t, d in records if d and t < tau})
    risk = [[int(t >= s) for t, d in records] for s in times]
    events = [[int(d and Q(t) == s) for t, d in records] for s in times]
    ys = list(map(sum, risk))
    ds = list(map(sum, events))
    surv = []
    s = Q(1)
    for y, d in zip(ys, ds):
        s *= 1 - Q(d, y)
        surv.append(s)
    ends = times[1:] + [tau]
    pieces = [(h - l) * s for l, h, s in zip(times, ends, surv)]
    tails = [sum(pieces[j:], Q(0)) for j in range(len(times))]
    rmst = (times[0] if times else tau) + sum(pieces, Q(0))
    scores = [Q(0) for _ in records]
    green = Q(0)
    for j, (y, d, A) in enumerate(zip(ys, ds, tails)):
        if y == d:
            assert A == 0
            continue
        green += A * A * Q(d, y * (y - d))
        for i in range(n):
            scores[i] -= A * (events[j][i] - risk[j][i] * Q(d, y)) / (y - d)
    return (rmst, scores, green)

def float_scores(records, tau, weights=None, check_support=True):
    if check_support:
        validate(records, tau)
    x = np.asarray(records, float)
    t = x[:, 0]
    e = x[:, 1]
    n = len(t)
    w = np.ones(n) if weights is None else np.asarray(weights, float)
    ts = np.unique(t[(e == 1) & (t < tau)])
    if len(ts) == 0:
        return (float(tau), np.zeros(n), 0.0)
    Y = t[None, :] >= ts[:, None]
    E = (t[None, :] == ts[:, None]) & (e[None, :] == 1)
    y = Y @ w
    d = E @ w
    s = np.cumprod(1 - d / y)
    pieces = np.diff(np.r_[ts, tau]) * s
    A = np.cumsum(pieces[::-1])[::-1]
    good = y > d
    g = -(A[good] / (y[good] - d[good])) @ (E[good] - Y[good] * (d[good] / y[good])[:, None])
    v = np.sum(A[good] ** 2 * d[good] / (y[good] * (y[good] - d[good])))
    return (float(ts[0] + pieces.sum()), g, float(v))

def ledger_edge(p, d, model):
    x, a = p
    y, b = d
    if model == 'separate':
        return b == 0 or (x <= y if a else x < y)
    if model == 'common':
        return x <= y and (bool(a) or (b == 0 and x == y))
    raise ValueError('unknown censor model')

def event_table(records):
    return [(
        t,
        int(sum((x >= t for x, e in records))),
        sum((int(x == t and e == 1) for x, e in records))) for t in sorted({x for x, e in records if e})]

def cells_from_table(n, table, checkpoints=(), support=None):
    tab = sorted(table)
    cps = dict(checkpoints)
    if len(tab) != len({t for t, y, d in tab}):
        raise ValueError('duplicate event times')
    for j, (t, y, d) in enumerate(tab):
        if t < 0 or any((int(v) != v for v in (y, d))) or (not 0 < d <= y <= n):
            raise ValueError('bad table')
        if j and tab[j - 1][1] - tab[j - 1][2] < y:
            raise ValueError('negative censor count')
    tab = [(t, int(y), int(d)) for t, y, d in tab]
    boundaries = sorted(set([
        0] + [t for t,
        y,
        d in tab] + list(cps) + ([] if support is None else [support[0]])))
    risks = {t: y for t, y, d in tab}
    ev = {t: d for t, y, d in tab}
    if support is not None:
        if support[0] in cps and cps[support[0]] != support[1]:
            raise ValueError('conflicting support')
        cps[support[0]] = support[1]
    for t, r in cps.items():
        if t < 0 or int(r) != r or (not 0 <= r <= n):
            raise ValueError('bad checkpoint')
        if t in risks and risks[t] != r:
            raise ValueError('conflicting checkpoint')
        risks[t] = r
    if 0 not in risks:
        risks[0] = n
    elif risks[0] != n:
        raise ValueError('origin risk must be n')
    out = []
    for j, l in enumerate(boundaries):
        y = risks[l]
        d = ev.get(l, 0)
        if d:
            out.append((1, l, l, d))
        h = boundaries[j + 1] if j + 1 < len(boundaries) else float('inf')
        c = y - d - (risks[h] if np.isfinite(h) else 0)
        if c < 0:
            raise ValueError('negative censor count at checkpoint')
        if c:
            out.append((0, l, h, c))
    if sum((c[3] for c in out)) != n:
        raise ValueError('cell counts do not sum to n')
    return out

def expand(cells):
    return [(a, l, h) for a, l, h, k in cells for _ in range(k)]

def cell_edge(p, d, model):
    a, l, h = p
    b, L, H = d
    if model == 'separate':
        return b == 0 or (l <= L if a else l < L)
    if model != 'common':
        raise ValueError('unknown model')
    if a and b:
        return l <= L
    if a and (not b):
        return l < H
    if not a and b:
        return False
    return max(l, L) < min(h, H)

def instantiate(pc, dc, pi, model):
    p = [(l, a) for a, l, h in pc]
    d = [(l, a) for a, l, h in dc]
    for i, j in enumerate(pi):
        if not cell_edge(pc[i], dc[j], model):
            raise ValueError('forbidden cell match')
        if model == 'common' and dc[j][0] == 0:
            c = max(pc[i][1], dc[j][1])
            d[j] = (c, 0)
            if pc[i][0] == 0:
                p[i] = (c, 0)
    return (p, d)

def exhaustive(up, ud, graph):
    vs = []
    for pi in permutations(range(len(up))):
        if all((graph[i, j] for i, j in enumerate(pi))):
            vs.append((sum(((up[i] - ud[j]) ** 2 for i, j in enumerate(pi))), pi))
    if not vs:
        raise ValueError('infeasible')
    return (min(vs), max(vs), len(vs))

def exact_minflow(up, ud, graph, sign=1):
    n = len(up)
    source = 2 * n
    sink = source + 1
    vcount = sink + 1
    adj = [[] for _ in range(vcount)]

    def add(u, v, cost):
        f = [v, 1, cost, len(adj[v])]
        r = [u, 0, -cost, len(adj[u])]
        adj[u].append(f)
        adj[v].append(r)
    for i in range(n):
        add(source, i, Q(0))
        add(n + i, sink, Q(0))
    for i in range(n):
        for j in range(n):
            if graph[i, j]:
                add(i, n + j, sign * (up[i] - ud[j]) ** 2)
    value = Q(0)
    for k in range(n):
        dist = [None] * vcount
        dist[source] = Q(0)
        pred = [None] * vcount
        for it in range(vcount - 1):
            changed = False
            for u in range(vcount):
                if dist[u] is None:
                    continue
                for ix, (v, cap, cost, rev) in enumerate(adj[u]):
                    if cap and (dist[v] is None or dist[u] + cost < dist[v]):
                        dist[v] = dist[u] + cost
                        pred[v] = (u, ix)
                        changed = True
            if not changed:
                break
        if dist[sink] is None:
            raise ValueError('infeasible')
        value += dist[sink]
        v = sink
        while v != source:
            u, ix = pred[v]
            e = adj[u][ix]
            e[1] -= 1
            adj[v][e[3]][1] += 1
            v = u
    pi = []
    for i in range(n):
        pi.append(next((e[0] - n for e in adj[i] if n <= e[0] < 2 * n and e[1] == 0)))
    potentials = [Q(0)] * vcount
    for it in range(vcount):
        changed = False
        for u in range(vcount):
            for v, cap, cost, rev in adj[u]:
                if cap and potentials[v] > potentials[u] + cost:
                    potentials[v] = potentials[u] + cost
                    changed = True
        if not changed:
            break
    if changed:
        raise AssertionError('negative cycle')
    beta = potentials[n:2 * n]
    alpha = [sign * (up[i] - ud[j]) ** 2 - beta[j] for i, j in enumerate(pi)]
    assert all((alpha[i] + beta[j] <= sign * (up[i] - ud[j]) ** 2 for i in range(n) for j in range(n) if graph[i, j]))
    assert sum(alpha) + sum(beta) == value
    return (sign * value, pi, alpha, beta)

def dense_graph(p, d, model, level='R1'):
    if level == 'R1':
        p = np.asarray(p)
        d = np.asarray(d)
        x = p[:, 0, None]
        y = d[None, :, 0]
        a = p[:, 1, None]
        b = d[None, :, 1]
        return (
            b == 0) | np.where(a == 1,
            x <= y,
            x < y) if model == 'separate' else (x <= y) & ((a == 1) | (b == 0) & (x == y))
    p = np.asarray(p)
    d = np.asarray(d)
    a = p[:, 0, None]
    b = d[None, :, 0]
    l = p[:, 1, None]
    L = d[None, :, 1]
    h = p[:, 2, None]
    H = d[None, :, 2]
    if model == 'separate':
        return (b == 0) | np.where(a == 1, l <= L, l < L)
    return np.where(a == 1, np.where(b == 1, l <= L, l < H), (b == 0) & (np.maximum(l, L) < np.minimum(h, H)))
