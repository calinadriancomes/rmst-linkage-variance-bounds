"""Literal, support-conditioned and paired bootstrap definitions."""
import numpy as np

class Margin:
    """Weighted Kaplan–Meier calculation retaining support flags for every attempt."""

    def __init__(self, r, tau):
        r = np.asarray(r)
        self.n = len(r)
        self.tau = tau
        self.times, self.codes = np.unique(r[:, 0], return_inverse=True)
        self.e = r[:, 1]

    def curve(self, weights):
        total = np.bincount(self.codes, weights=weights, minlength=len(self.times))
        de = np.bincount(self.codes, weights=weights * self.e, minlength=len(self.times))
        y = np.cumsum(total[::-1])[::-1]
        keep = y > 0
        t = self.times[keep]
        d = de[keep]
        y = y[keep]
        s = np.cumprod(1 - d / y)
        before = t < self.tau
        tt = t[before]
        ss = s[before]
        mean = self.tau if len(tt) == 0 else tt[0] + np.sum(np.diff(np.r_[tt, self.tau]) * ss)
        supported = self.times[total > 0].max() >= self.tau or s[-1] == 0
        return (float(mean), t[d > 0], s[d > 0], bool(supported))


def daly_order(tp, sp, td, sd):
    """Literal event-time order filter for independently resampled endpoint records."""
    if not len(td):
        return False
    ix = np.searchsorted(tp, td, side='right') - 1
    return False if np.any(ix < 0) else bool(np.all(sd >= sp[ix]))


def bootstrap(p, d, tau, scale, rng, mode, target=499, cap=4990, injected=None):
    """Retain every attempted index pair, rejection and cap; no replacement of original samples."""
    if mode not in ['literal', 'support_conditioned', 'paired']:
        raise ValueError('unknown bootstrap mode')
    n = len(p)
    pm = Margin(p, tau)
    dm = Margin(d, tau)
    accepted = []
    attempts = []
    for attempt in range(1, cap + 1):
        if injected is None:
            ip = rng.integers(n, size=n)
            idd = ip if mode == 'paired' else rng.integers(n, size=n)
        else:
            try:
                ip, idd = next(injected)
            except StopIteration:
                break
            if mode == 'paired' and (not np.array_equal(ip, idd)):
                raise ValueError('paired index mismatch')
        wp = np.bincount(ip, minlength=n)
        wd = np.bincount(idd, minlength=n)
        mp, tp, sp, supP = pm.curve(wp)
        md, td, sd, supD = dm.curve(wd)
        if mode == 'support_conditioned' and (not (supP and supD)):
            status = 'SUPPORT_REJECT'
        elif mode != 'paired' and (not daly_order(tp, sp, td, sd)):
            status = 'ORDER_REJECT'
        else:
            status = 'ACCEPT'
        attempts.append(dict(
            attempt=attempt,
            status=status,
            support_P=supP,
            support_D=supD,
            rmst_P=mp * scale,
            rmst_D=md * scale))
        if status == 'ACCEPT':
            accepted.append((mp * scale, md * scale))
        if len(accepted) == target:
            break
    complete = len(accepted) == target
    out = dict(
        mode=mode,
        status='COMPLETE' if complete else 'CAPPED' if len(attempts) == cap else 'INJECTED_STREAM_ENDED',
        complete=complete,
        target=target,
        cap=cap,
        attempted=len(attempts),
        accepted=len(accepted),
        rejected_order=sum((x['status'] == 'ORDER_REJECT' for x in attempts)),
        rejected_support=sum((x['status'] == 'SUPPORT_REJECT' for x in attempts)),
        unsupported_attempts=sum((not (x['support_P'] and x['support_D']) for x in attempts)),
        unsupported_accepted=sum((x['status'] == 'ACCEPT' and (not (x['support_P'] and x['support_D'])) for x in attempts)),
        accepted_RMST=accepted,
        attempts=attempts)
    if complete:
        a = np.asarray(accepted)
        sd = a.std(axis=0, ddof=1)
        rho = float(np.corrcoef(a.T)[0, 1]) if np.all(sd > 0) else None
        out.update(
            sd_P=float(sd[0]),
            sd_D=float(sd[1]),
            correlation=rho,
            variance_difference=float(np.var(a[:, 1] - a[:, 0], ddof=1)))
    return out
