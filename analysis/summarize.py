"""Summaries of retained outcomes with separate support and method denominators."""
from pathlib import Path
import csv
import gzip
import json
import math
import sqlite3
import zlib
import numpy as np

def wilson(k, n, z=1.959963984540054):
    """Binomial Wilson limits and Monte Carlo standard error."""
    if not n:
        return (None, None, None, None)
    p = k / n
    den = 1 + z * z / n
    mid = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (p, mid - half, mid + half, math.sqrt(p * (1 - p) / n))


def distribution(values):
    """Finite-value summaries with an explicit denominator."""
    a = np.asarray(values, float)
    if not len(a):
        return dict(n=0, mean=None, MCSE=None, min=None, p05=None, median=None, p95=None, max=None)
    return dict(
        n=len(a),
        mean=float(a.mean()),
        MCSE=float(a.std(ddof=1) / math.sqrt(len(a))) if len(a) > 1 else None,
        min=float(a.min()),
        p05=float(np.quantile(a, 0.05)),
        median=float(np.median(a)),
        p95=float(np.quantile(a, 0.95)),
        max=float(a.max()))


def variances(result, boot=None):
    """Variance methods and bootstrap statuses for a retained original sample."""
    out = {k: dict(
        variance=result[k],
        status='COMPLETE',
        role='primary upper interval' if k.endswith('upper') else 'reference/diagnostic') for k in ['oracle',
        'independence',
        'Cauchy_lower',
        'Cauchy_upper',
        'unrestricted_lower',
        'unrestricted_upper']}
    for name, b in result['levels'].items():
        for end in ['lower', 'upper']:
            out[
                name + '_' + end] = dict(variance=b[end + '_study_units'],
                status='COMPLETE',
                role='primary upper interval' if end == 'upper' else 'lower-end diagnostic',
                covariance_range_width=b['covariance_upper_study_units'] - b['covariance_lower_study_units'],
                edge_count=b['edge_count'])
    if boot:
        for mode, b in boot.items():
            prefix = mode + '_'
            base = dict(
                status=b['status'],
                role='algorithmic bootstrap diagnostic',
                unsupported_accepted=b['unsupported_accepted'])
            out[prefix + 'bootstrap_SD'] = dict(base, variance=b.get('variance_difference'))
            if mode != 'paired':
                rho = b.get('correlation')
                v = None if rho is None else result['var_P'] + result['var_D'] - 2 * rho * math.sqrt(result['var_P'] * result['var_D'])
                out[
                    prefix + 'A3_original_SE'] = dict(base,
                    variance=v,
                    status=b['status'] if not b['complete'] else 'UNDEFINED_CORRELATION' if rho is None else 'COMPLETE')
    return out


def write_csv(path, rows, fields):
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fields, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def original_support(payload, status):
    """Return True/False for observed marginal support, or None before it is known.

    Re-evaluate stored observations rather than classifying every failed calculation
    as unsupported. This also reads older payloads without a support-stage flag.
    """
    if payload and all(k in payload for k in ("P_records", "D_records", "tau_internal")):
        import sys
        sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
        from observation_model import validate
        try:
            validate(payload["P_records"], payload["tau_internal"])
            validate(payload["D_records"], payload["tau_internal"])
        except (ValueError, TypeError, IndexError):
            return False
        return True
    if payload and isinstance(payload.get("support_valid"), bool):
        return payload["support_valid"]
    if status == "COMPLETED":
        return True
    if status == "SUPPORT_FAIL":
        return False
    return None


def export(root, design, truths, validation_run=False):
    """Read existing keyed outcomes and write summaries without drawing samples."""
    root = Path(root)
    out = root / 'reports'
    out.mkdir(exist_ok=True)
    db = sqlite3.connect(root / 'analysis_state/checkpoints.sqlite')
    cells = {x['id']: x for x in design['cells']}
    planned = {i: c['replicates'] for i, c in cells.items()}
    statuses = {i: {} for i in cells}
    support_counts = {i: {True: 0, False: 0, None: 0} for i in cells}
    bucket = {}
    bias = {i: [] for i in cells}
    bootcounts = {}
    keyrows = []
    bootstraprows = []
    methodfields = [
        'cell',
        'replicate',
        'method',
        'role',
        'status',
        'variance',
        'coverage',
        'lower',
        'upper',
        'width',
        'gain_unrestricted',
        'variance_oracle_ratio',
        'width_oracle_ratio',
        'covariance_range_width',
        'edge_count',
        'unsupported_accepted']
    mf = (out / 'method_results.csv').open('w', newline='', encoding='utf-8')
    mw = csv.DictWriter(mf, methodfields, extrasaction='ignore')
    mw.writeheader()
    af = gzip.open(out / 'bootstrap_attempts.csv.gz', 'wt', newline='', encoding='utf-8')
    aw = csv.DictWriter(af, ['cell', 'replicate', 'mode', 'attempt', 'status', 'support_P', 'support_D', 'rmst_P', 'rmst_D'])
    aw.writeheader()
    raw = gzip.open(out / 'lossless_replicates.jsonl.gz', 'wt', encoding='utf-8')
    for cid, rep, status, blob in db.execute('SELECT cell,rep,status,payload FROM jobs ORDER BY cell,rep'):
        p = json.loads(zlib.decompress(blob)) if blob is not None else None
        raw.write(json.dumps(
            dict(cell=cid, replicate=rep, status=status, payload=p),
            allow_nan=False,
            separators=(',', ':')) + '\n')
        statuses[cid][status] = statuses[cid].get(status, 0) + 1
        support_counts[cid][original_support(p, status)] += 1
        theta = truths[str(cid)]['theta']
        # A retained partial result is not a completed original calculation.
        r = p.get('result') if p and status == 'COMPLETED' else None
        keyrows.append(dict(
            cell=cid,
            replicate=rep,
            status=status,
            truth=theta,
            point=r['point'] if r else None,
            bias=r['point'] - theta if r else None,
            seconds=p.get('seconds') if p else None,
            error=p.get('error') if p else None,
            optimality_check=bool(p and p.get('optimality_check_reasons'))))
        selected = cid in design['bootstrap_comparator_block']['selected_cell_ids'] and rep < 500
        if selected:
            for mode in ['literal', 'support_conditioned', 'paired']:
                b = p.get('bootstrap', {}).get(mode) if p else None
                row = dict(
                    cell=cid,
                    replicate=rep,
                    mode=mode,
                    status=b['status'] if b else 'SKIPPED_' + status,
                    target=499,
                    cap=4990,
                    attempted=b['attempted'] if b else 0,
                    accepted=b['accepted'] if b else 0,
                    rejected_order=b['rejected_order'] if b else 0,
                    rejected_support=b['rejected_support'] if b else 0,
                    unsupported_attempts=b['unsupported_attempts'] if b else 0,
                    unsupported_accepted=b['unsupported_accepted'] if b else 0,
                    correlation=b.get('correlation') if b else None)
                bootstraprows.append(row)
                if b:
                    for a in b['attempts']:
                        aw.writerow(dict(cell=cid, replicate=rep, mode=mode, **a))
                counts = bootcounts.setdefault((cid, mode), dict(selected=0, attempted_cases=0, complete=0, capped=0, attempts=0, accepted=0, rejected_order=0, rejected_support=0, unsupported_attempts=0, unsupported_accepted=0, undefined_correlation=0))
                counts['selected'] += 1
                if b:
                    counts['attempted_cases'] += 1
                    counts['complete'] += b['complete']
                    counts['capped'] += b['status'] == 'CAPPED'
                    counts['undefined_correlation'] += b['complete'] and b.get('correlation') is None
                    for k in [
                        'accepted',
                        'rejected_order',
                        'rejected_support',
                        'unsupported_attempts',
                        'unsupported_accepted']:
                        counts[k] += b[k]
                    counts['attempts'] += b['attempted']
        if not r:
            continue
        bias[cid].append(r['point'] - theta)
        for name, x in variances(r, p.get('bootstrap')).items():
            v = x.get('variance')
            defined = x['status'] == 'COMPLETE' and v is not None
            row = dict(cell=cid, replicate=rep, method=name, **x)
            stats = bucket.setdefault((cid, name), dict(status={}, coverage=[], width=[], gain=[], variance_oracle_ratio=[], width_oracle_ratio=[], covariance_range_width=[], edges=[]))
            stats['status'][x['status']] = stats['status'].get(x['status'], 0) + 1
            if defined:
                if v < 0:
                    raise AssertionError('negative reported variance')
                half = design['normal_quantile'] * math.sqrt(v)
                width = 2 * half
                covered = abs(r['point'] - theta) <= half
                gain = 1 - math.sqrt(v / r['unrestricted_upper']) if r['unrestricted_upper'] > 0 else None
                ratio = v / r['oracle'] if r['oracle'] > 0 else None
                wr = math.sqrt(ratio) if ratio is not None else None
                row.update(
                    coverage=covered,
                    lower=r['point'] - half,
                    upper=r['point'] + half,
                    width=width,
                    gain_unrestricted=gain,
                    variance_oracle_ratio=ratio,
                    width_oracle_ratio=wr)
                stats['coverage'].append(covered)
                stats['width'].append(width)
                for key, value in [
                    ('gain', gain),
                    ('variance_oracle_ratio', ratio),
                    ('width_oracle_ratio', wr),
                    ('covariance_range_width', x.get('covariance_range_width')),
                    ('edges', x.get('edge_count'))]:
                    if value is not None:
                        stats[key].append(value)
            mw.writerow(row)
    mf.close()
    af.close()
    raw.close()
    db.close()
    summary = []
    denominators = []
    for cid in cells:
        st = statuses[cid]
        n = planned[cid]
        attempted = n - st.get('NOT_ATTEMPTED', 0)
        a = np.array(bias[cid])
        completed = st.get('COMPLETED', 0)
        support = support_counts[cid][True]
        denominators.append(dict(
            cell=cid,
            module=cells[cid]['module'],
            planned=n,
            attempted=attempted,
            support_valid=support,
            support_unknown_count=support_counts[cid][None],
            completed=completed,
            support_failed=st.get('SUPPORT_FAIL', 0),
            calculation_failed=st.get('CALCULATION_FAILED', 0),
            not_attempted=st.get('NOT_ATTEMPTED', 0),
            running=st.get('RUNNING', 0),
            completion_proportion=completed / attempted if attempted else None,
            support_proportion=support / attempted if attempted else None,
            support_includes_calculation_failures=True,
            bias=float(a.mean()) if len(a) else None,
            bias_MCSE=float(a.std(ddof=1) / math.sqrt(len(a))) if len(a) > 1 else None,
            RMSE=float(np.sqrt(np.mean(a * a))) if len(a) else None))
    for (cid, name), s in sorted(bucket.items()):
        cv = s['coverage']
        n = len(cv)
        k = sum(cv)
        rate, lo, hi, se = wilson(k, n)
        row = dict(
            cell=cid,
            module=cells[cid]['module'],
            method=name,
            planned=planned[cid],
            completed=statuses[cid].get('COMPLETED', 0),
            support_valid=support_counts[cid][True],
            method_completed=n,
            covered=k,
            coverage=rate,
            Wilson_lower=lo,
            Wilson_upper=hi,
            coverage_MCSE=se,
            covered_planned_operational_yield=k / planned[cid] if statuses[cid].get('NOT_ATTEMPTED', 0) == 0 and statuses[cid].get('RUNNING', 0) == 0 else None,
            status_counts=s['status'])
        for key in [
            'width',
            'gain',
            'variance_oracle_ratio',
            'width_oracle_ratio',
            'covariance_range_width',
            'edges']:
            row[key] = distribution(s[key])
        summary.append(row)
    result = dict(
        validation_run=validation_run,
        planned=sum(planned.values()),
        cell_denominators=denominators,
        method_summaries=summary,
        bootstrap_summaries=[dict(cell=c, mode=m, **v) for (c, m), v in sorted(bootcounts.items())],
        assumptions='Informative-censoring cells are assumption-violation diagnostics. Bootstrap outputs are algorithmic diagnostics.',
        zero_ratio_denominators='Undefined ratios are omitted from finite-value summaries and retained in their explicit n denominators.')
    (out / 'CELL_SUMMARIES.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    write_csv(
        out / 'replicate_status.csv',
        keyrows,
        ['cell', 'replicate', 'status', 'truth', 'point', 'bias', 'seconds', 'error', 'optimality_check'])
    write_csv(out / 'cell_denominators.csv', denominators, list(denominators[0]))
    write_csv(
        out / 'bootstrap_status.csv',
        bootstraprows,
        ['cell', 'replicate', 'mode', 'status', 'target', 'cap', 'attempted', 'accepted', 'rejected_order', 'rejected_support', 'unsupported_attempts', 'unsupported_accepted', 'correlation'])
    (out / 'output_format.md').write_text("# Analysis output format\n\n`analysis_state/checkpoints.sqlite` stores one outcome for every planned (cell, replicate) pair. Payloads are UTF-8 JSON compressed with zlib. `lossless_replicates.jsonl.gz` contains the observations, random-number branches, bound witnesses, optimality residual diagnostics and bootstrap attempts. Null denotes an undefined or unexecuted value.\n\n`replicate_status.csv` lists every planned sample. `cell_denominators.csv` separates valid marginal support from completed calculations and reports support failures, calculation failures, unattempted samples and running calculations. `completion_proportion` is completed/attempted; `support_proportion` is support_valid/attempted. Support can be unknown before observations are available.\n\n`method_results.csv` includes upper intervals, reference methods, lower-end diagnostics and bootstrap diagnostics. Intervals use the original point center and declared normal quantile. `CELL_SUMMARIES.json` reports support_valid, completed primary calculations and each method's completed intervals separately. Coverage, Wilson limits and MCSE use method-specific denominators. Covered/planned is operational yield, not unconditional coverage. Widths, gains, ratios and covariance-width summaries each carry a finite-value count. Bootstrap original-SE uses marginal Greenwood variances and accepted-draw correlation; filtered SD uses the same accepted draws. Undefined correlation remains undefined.\n\n`bootstrap_status.csv` includes skipped originals. `bootstrap_attempts.csv.gz` records executed attempts, rejection status, support flags and RMST pairs. Attempts are not replacement original samples. Capped cases retain partial accepted draws without completed-comparator coverage. Unsupported literal or paired attempts remain labeled as algorithmic output.\n", encoding='utf-8')
    return result
