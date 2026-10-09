"""Optional full analysis. The ordinary reproduction route uses stored results."""
from pathlib import Path
import argparse
import csv
import datetime
import hashlib
import json
import math
import sys
import time
import traceback
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import estimation as c
from simulate import generate
from bootstrap import bootstrap
from storage import Store
c.generate = generate

def memory():
    return {"measurement_available": False}


def bootstrap_job(p, d, tau, scale, entropy, mode, forced_indices=None):
    beg = time.process_time()
    rng = np.random.default_rng(np.random.SeedSequence(entropy)) if forced_indices is None else None
    result = bootstrap(
        p,
        d,
        tau,
        scale,
        rng,
        mode,
        499,
        4990,
        None if forced_indices is None else iter(forced_indices))
    return (result, dict(CPU_seconds=time.process_time() - beg, memory=memory()))


def clean(x):
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, np.generic):
        return x.item()
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    return x


def interval_metrics(result, z):
    variances = {key: result[
        key] for key in ['oracle',
        'independence',
        'Cauchy_lower',
        'Cauchy_upper',
        'unrestricted_lower',
        'unrestricted_upper']}
    for name, b in result['levels'].items():
        variances[name + '_lower'] = b['lower_study_units']
        variances[name + '_upper'] = b['upper_study_units']
    metrics = {}
    for name, v in variances.items():
        metrics[name + '_width'] = 2 * z * math.sqrt(v)
        if result['unrestricted_upper'] > 0:
            metrics[name + '_gain'] = 1 - math.sqrt(v / result['unrestricted_upper'])
        if result['oracle'] > 0:
            metrics[name + '_variance_oracle_ratio'] = v / result['oracle']
    return metrics


def execute(root, definition, injected=None):
    design = definition['design']
    dh = definition['design_SHA256']
    keys = [(cell['id'], rep) for cell in design['cells'] for rep in range(cell['replicates'])]
    store = Store(root / 'analysis_state/checkpoints.sqlite', dh, keys)
    pool = ProcessPoolExecutor(max_workers=3, mp_context=multiprocessing.get_context('spawn')) if injected is None or injected.get('bootstrap_indices') else None
    statefile = root / 'analysis_state/run_summary.json'
    start = time.perf_counter()
    runtimes = []
    stop = None

    def write_run_summary():
        counts = {s: n for s, n in store.db.execute('SELECT status,COUNT(*) FROM jobs GROUP BY status')}
        x = dict(
            design_SHA256=dh,
            planned=len(keys),
            status_counts=counts,
            attempted=len(keys) - counts.get('NOT_ATTEMPTED', 0),
            completed=counts.get('COMPLETED', 0),
            support_failed=counts.get('SUPPORT_FAIL', 0),
            status='CALCULATION_STOPPED' if stop else 'RUNNING' if counts.get('NOT_ATTEMPTED', 0) or counts.get('RUNNING', 0) else 'COMPLETED',
            exact_stop=stop,
            wall_seconds_this_invocation=time.perf_counter() - start,
            memory=memory(),
            UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
            simulation_seed=20261004 if injected is None else None,
            validation_run=injected is not None,
            no_replacement=True,
            no_optional_stopping=True)
        statefile.write_text(json.dumps(x, indent=2) + '\n')
        return x
    for cell in design['cells']:
        cid = cell['id']
        beg = time.perf_counter()
        cpu = time.process_time()
        supported = 0
        extrema = {}
        for rep in range(cell['replicates']):
            oldstate, payload = store.get(cid, rep)
            if oldstate not in ('NOT_ATTEMPTED', 'RUNNING'):
                if oldstate == 'COMPLETED':
                    supported += 1
                    oldresult = payload['result']
                    for name, value in interval_metrics(oldresult, design['normal_quantile']).items():
                        for label, better in [('min', lambda a, b: a < b), ('max', lambda a, b: a > b)]:
                            key = name + '_' + label
                            if key not in extrema or better(value, extrema[key][0]):
                                extrema[key] = (value, rep)
                continue
            resumed = store.start(cid, rep)
            diag = dict(
                cell=cid,
                replicate=rep,
                status='RUNNING',
                resumed_interrupted_key=resumed,
                seed_entropy=[design['base_seed'], cid, rep],
                branches={'generation': 0, 'literal': 1, 'support_conditioned': 2, 'paired': 3})
            t0 = time.perf_counter()
            try:
                p, d, tau, scale = c.generate(
                    cell,
                    rep,
                    design['base_seed']) if injected is None else injected['generate'](cell,
                    rep)
                diag.update(
                    P_records=p.tolist(),
                    D_records=d.tolist(),
                    tau_internal=tau,
                    scale=scale,
                    true_theta=definition['truth_values'][str(cid)]['theta'],
                    assumption_violation=cell['module'] == 'assumption_violation')
                if len(p) != cell['n'] or len(d) != cell['n']:
                    raise AssertionError('patient universe mismatch')
                try:
                    c.reference.validate(p, tau)
                    c.reference.validate(d, tau)
                except ValueError as exc:
                    if 'unsupported' not in str(exc):
                        raise
                    diag.update(
                        status='SUPPORT_FAIL',
                        support_valid=False,
                        error=str(exc),
                        bootstrap_status='SKIPPED_ORIGINAL_SUPPORT_FAIL')
                    store.finish(cid, rep, 'SUPPORT_FAIL', clean(diag))
                    write_run_summary()
                    continue
                diag['support_valid'] = True
                result = c.evaluate(p, d, tau, scale, cell.get('censor_model') == 'common', supported < 25)
                if cell['kind'] == 'FGM':
                    for b in result['levels'].values():
                        if abs(b['upper_study_units'] - result['unrestricted_upper']) > 1e-10 * cell['tau'] ** 2:
                            raise AssertionError('negative-control gain')
                metrics = interval_metrics(result, design['normal_quantile'])
                recheck_reasons = ['initial_sample_check'] if supported < 25 else []
                for name, value in metrics.items():
                    for label, better in [('min', lambda a, b: a < b), ('max', lambda a, b: a > b)]:
                        key = name + '_' + label
                        if key not in extrema or better(value, extrema[key][0]):
                            extrema[key] = (value, rep)
                            recheck_reasons.append('new_extreme:' + key)
                if recheck_reasons and supported >= 25:
                    result = c.evaluate(p, d, tau, scale, cell.get('censor_model') == 'common', True)
                diag.update(result=result, optimality_check_reasons=recheck_reasons)
                if cid in design['bootstrap_comparator_block']['selected_cell_ids'] and rep < 500:
                    diag['bootstrap'] = {}
                    if pool is not None:
                        futures = {}
                        for mode, branch in [('literal', 1), ('support_conditioned', 2), ('paired', 3)]:
                            forced = None if injected is None else injected['bootstrap_indices'](mode, 4990)
                            futures[
                                mode] = pool.submit(bootstrap_job,
                                p,
                                d,
                                tau,
                                scale,
                                [design['base_seed'], cid, rep, branch],
                                mode,
                                forced)
                        diag['bootstrap_resources'] = {}
                        for mode, future in futures.items():
                            diag['bootstrap'][mode], diag['bootstrap_resources'][mode] = future.result()
                    else:
                        for mode, branch in [('literal', 1), ('support_conditioned', 2), ('paired', 3)]:
                            diag[
                                'bootstrap'][mode] = injected['bootstrap'](p,
                                d,
                                tau,
                                scale,
                                None,
                                mode,
                                499,
                                4990)
                diag['status'] = 'COMPLETED'
                supported += 1
                diag['seconds'] = time.perf_counter() - t0
                store.finish(cid, rep, 'COMPLETED', clean(diag))
            except Exception as exc:
                stop = dict(cell=cid, replicate=rep, error=repr(exc), traceback=traceback.format_exc())
                diag.update(
                    status='CALCULATION_FAILED',
                    error=repr(exc),
                    traceback=stop['traceback'],
                    seconds=time.perf_counter() - t0)
                if store.status(cid, rep) == 'RUNNING':
                    store.finish(cid, rep, 'CALCULATION_FAILED', clean(diag))
                (root / 'analysis_state/calculation_failure.json').write_text(json.dumps(stop, indent=2) + '\n')
                write_run_summary()
                break
            if (rep + 1) % 25 == 0:
                x = write_run_summary()
                print(
                    'cell',
                    cid,
                    'index',
                    rep + 1,
                    'attempted',
                    x['attempted'],
                    'completed',
                    x['completed'],
                    'wall_s',
                    round(x['wall_seconds_this_invocation'], 1),
                    flush=True)
        runtimes.append(dict(
            cell=cid,
            wall_seconds=time.perf_counter() - beg,
            CPU_seconds=time.process_time() - cpu,
            supported=supported,
            memory=memory(),
            extrema=extrema))
        (root / 'analysis_state/cell_runtime.json').write_text(json.dumps(clean(runtimes), indent=2) + '\n')
        if stop:
            break
    final = write_run_summary()
    store.close()
    if pool is not None:
        pool.shutdown(wait=True)
    from summarize import export
    export(root, design, definition['truth_values'], validation_run=injected is not None)
    print(json.dumps(final, indent=2), flush=True)
    return 2 if stop else 0
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-full-study",
        action="store_true",
        help="Explicitly enable the optional 39,000-sample analysis")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not args.run_full_study:
        parser.error("Use code/check.py and code/reproduce.py for normal reproduction")
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    (root / "analysis_state").mkdir(exist_ok=True)
    if (root / "analysis_state/calculation_failure.json").exists():
        raise RuntimeError("An unresolved calculation failure prevents continuation")
    base = Path(__file__).resolve().parents[1]
    design = json.loads((base / "data/scenario_definitions.json").read_text())
    targets = json.loads((base / "data/population_targets.json").read_text())
    analysis_definition = json.loads((base / "data/analysis_definition.json").read_text())
    if len(design["cells"]) != 27 or sum(x["replicates"] for x in design["cells"]) != 39000 or design["base_seed"] != 20261004:
        raise ValueError("Study definition mismatch")
    if hashlib.sha256((base / "data/scenario_definitions.json").read_bytes()).hexdigest() != analysis_definition["design_file_SHA256"]:
        raise ValueError("Study definition checksum mismatch")
    if hashlib.sha256((base / "data/population_targets.json").read_bytes()).hexdigest() != analysis_definition["targets_file_SHA256"]:
        raise ValueError("Population target checksum mismatch")
    definition = {"design": design, "truth_values": targets, "design_SHA256": analysis_definition["design_SHA256"]}
    return execute(root, definition)

if __name__ == "__main__":
    sys.exit(main())
