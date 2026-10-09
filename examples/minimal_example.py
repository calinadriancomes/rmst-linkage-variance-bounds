from pathlib import Path
import sys,json
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/"src"))
from method import from_records,from_tables,upper_envelope_interval
from observation_model import event_table
z=json.loads((root/"examples/tied_events.json").read_text())
p,d,tau=z["P"],z["D"],z["tau"]
r1=from_records(p,d,tau)
r2=from_tables(len(p),event_table(p),event_table(d),tau,sum(t>=tau for t,e in p),sum(t>=tau for t,e in d))
print(f"Restricted-mean contrast: {r1['point']:.6f}")
for label, result in [("R1",r1),("R2",r2)]:
    low,high=upper_envelope_interval(result["point"],result["upper"])
    print(f"{label}: variance endpoints [{result['lower']:.6f}, {result['upper']:.6f}]; "
          f"upper-variance interval [{low:.6f}, {high:.6f}]")
