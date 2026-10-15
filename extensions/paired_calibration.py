"""Additional paired-reference-only simulation diagnostic (not R1/R2 coverage).

The investigation used a fixed seed and 20,000 attempts for each of two
original parameter combinations, selected after inspection of the first study.
Run from the repository root:
    python extensions/paired_calibration.py
"""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.stats import norm
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from method import scores

SEED=2026101317
REPS=20000

def wilson(k,n):
    z=norm.ppf(.975);p=k/n;den=1+z*z/n
    h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    c=(p+z*z/(2*n))/den
    return [float(c-h),float(c+h)]

def run():
    outputs=[]
    for scen,(q,rate_w,c) in enumerate([(.25,.5,.25),(.5,2.,.75)]):
        rng=np.random.default_rng(np.random.SeedSequence([SEED,scen]));val=[];fail=0
        for rep in range(REPS):
            n=500;tau=2.
            p0=rng.exponential(1.,n)
            first=rng.uniform(size=n)<q
            dlatent=p0+np.where(first,0.,rng.exponential(1/rate_w,n))
            cp=rng.exponential(1/c,n);cd=rng.exponential(1/c,n)
            p=list(zip(np.minimum(p0,cp),p0<=cp))
            d=list(zip(np.minimum(dlatent,cd),dlatent<=cd))
            try:
                mp,up,_=scores(p,tau);md,ud,_=scores(d,tau)
            except ValueError:
                fail+=1;continue
            val.append((md-mp,float(np.sum((ud-up)**2))))
        vals=np.array(val);ests=vals[:,0];variances=vals[:,1]
        truth=(1-q)*((1-math.exp(-tau))/rate_w-math.exp(-rate_w*tau)/rate_w*(1-math.exp(-(1-rate_w)*tau))/(1-rate_w))
        k=int(np.sum(abs(ests-truth)<=norm.ppf(.975)*np.sqrt(variances)))
        outputs.append(dict(n=len(ests),covered=k,coverage=k/len(ests),wilson95=wilson(k,len(ests)),bias=float(ests.mean()-truth),rmse=float(np.sqrt(np.mean((ests-truth)**2))),variance_ratio=float(variances.mean()/ests.var(ddof=1)),mean_variance=float(variances.mean()),empirical_variance=float(ests.var(ddof=1)),q=q,lambda_W=rate_w,c=c,planned=REPS,support_failures=fail,n_patients=500,truth=truth))
    return dict(seed=SEED,replicates_per_scenario=REPS,scope='Subsequent paired-reference-only sensitivity simulation, not additional R1/R2 envelope evidence',results=outputs)

if __name__=='__main__':
    result=run()
    print(json.dumps(result,indent=2))
    reference=json.loads((Path(__file__).parent/'fresh_calibration_results.json').read_text())
    for got,expected in zip(result['results'],reference['results']):
        for key in ['n','covered','support_failures']:
            assert got[key]==expected[key],(key,got[key],expected[key])
        for key in ['coverage','variance_ratio','rmse','truth']:
            assert abs(got[key]-expected[key])<1e-12,(key,got[key],expected[key])
    print('The fixed-seed diagnostic matches the saved summary.')
