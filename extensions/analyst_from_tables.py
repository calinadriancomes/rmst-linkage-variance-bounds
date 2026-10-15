"""Execute the Rotterdam terminal-table illustration from unlinked marginal tables.

The analyst's inputs are two CSV files and metadata.json; neither the
original Rotterdam patient records nor full post-horizon follow-up is read.
Run from the repository root:
    python extensions/analyst_from_tables.py
"""
import csv
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / 'src'))
from horizon_release import from_terminal_tables
from method import upper_envelope_interval


def load_inputs(folder):
    folder = Path(folder)
    meta = json.loads((folder / 'metadata.json').read_text(encoding='utf-8'))
    n = meta['n']; tau = meta['tau']
    assert isinstance(n, int) and n > 0
    assert isinstance(tau, (float, int)) and math.isfinite(tau) and tau > 0
    def table(name):
        with (folder / f'{name}_event_risk.csv').open(newline='', encoding='utf-8') as f:
            rdr = csv.DictReader(f)
            if rdr.fieldnames != ['time','n_risk','n_event']:
                raise ValueError('The event-risk CSV fields must be time,n_risk,n_event.')
            rows=[]
            for r in rdr:
                t = float(r['time']); y = int(r['n_risk']); d = int(r['n_event'])
                if not math.isfinite(t) or not 0 <= t < tau:
                    raise ValueError('Terminal release must include only finite event times strictly before tau.')
                if not 1 <= d <= y <= n:
                    raise ValueError('Invalid exact event or risk multiplicity.')
                rows.append((t,y,d))
        if rows != sorted(rows) or len({t for t,_,_ in rows})!=len(rows):
            raise ValueError('Event times must be strictly increasing and unique.')
        return rows
    if meta['censoring_model'] != 'separate':
        raise ValueError('This example requires separately observed endpoint follow-up.')
    return meta,table('P'),table('D')


def analyze(folder,check_optimality=True):
    meta,tp,td=load_inputs(folder)
    n=meta['n'];tau=meta['tau']
    yp=meta['P_risk_tau'];yd=meta['D_risk_tau']
    if not isinstance(yp,int) or not isinstance(yd,int) or not 0<=yp<=n or not 0<=yd<=n:
        raise ValueError('Exact integer pre-event horizon counts are required for both endpoints.')
    result=from_terminal_tables(n,tp,td,tau,yp,yd,model='separate',check_optimality=check_optimality)
    low,high=upper_envelope_interval(result['point'],result['upper'])
    output={k:float(result[k]) for k in ['point','rmst_P','rmst_D','var_P','var_D','lower','upper']}
    output.update(n=n,tau=tau,P_risk_tau=yp,D_risk_tau=yd,interval_lower=low,interval_upper=high,width=high-low,edge_count=result['edge_count'])
    if check_optimality:
        output['optimality']={b:{k:result[b+'_optimality_check'][k] for k in ('gap','max_dual_violation','tolerance','numerical_status')} for b in ['lower','upper']}
    return output


def main():
    folder = Path(sys.argv[1]) if len(sys.argv)>1 else HERE/'rotterdam_terminal_example'
    result=analyze(folder)
    print(json.dumps(result,indent=2))
    # Full-history comparator is used only AFTER calculating from input CSVs.
    if folder.resolve() == (HERE/'rotterdam_terminal_example').resolve():
        expected=json.loads((folder/'expected_results.json').read_text())
        for k in ['point','upper','interval_lower','interval_upper']:
            if not math.isclose(result[k], expected[k], rel_tol=1e-9, abs_tol=1e-8):
                raise AssertionError(f'Comparison to historical R2: {k} differs.')
        print('Terminal-input upper interval matches the stored full-history R2 result.')

if __name__=='__main__': main()
