from pathlib import Path
import sys,unittest
import numpy as np
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root/'src'))
from method import from_records,from_tables,scores
from observation_model import event_table
from horizon_release import administratively_censor,from_horizon_records,from_horizon_tables,from_terminal_tables,from_terminal_records

class HorizonChecks(unittest.TestCase):
    def test_weight_derivatives_and_horizon_events(self):
        r=[(0,1),(1,0),(2,1),(2,0),(4,1)]
        c=administratively_censor(r,2)
        self.assertEqual(c,[(0,1),(1,0),(2,1),(2,0),(2,0)])
        a=scores(r,2);b=scores(c,2)
        self.assertEqual(a[0],b[0]);np.testing.assert_array_equal(a[1],b[1]);self.assertEqual(a[2],b[2])
    def test_minimal_record_lower_counterexample(self):
        p=[(0,0),(1,1),(3,0)];d=[(0,0),(1,1),(3,1)]
        f=from_records(p,d,2);c=from_horizon_records(3,administratively_censor(p,2),administratively_censor(d,2),2)
        self.assertAlmostEqual(f['lower'],1/8);self.assertEqual(c['lower'],0);self.assertAlmostEqual(f['upper'],3/8);self.assertEqual(f['upper'],c['upper'])
    def test_minimal_table_lower_counterexample(self):
        p=[(0,0),(1,1),(4,1)];d=[(0,0),(1,1),(3,1)]
        f=from_tables(3,event_table(p),event_table(d),2,1,1)
        a=from_horizon_tables(3,[x for x in event_table(p) if x[0]<=2],[x for x in event_table(d) if x[0]<=2],2,1,1)
        self.assertAlmostEqual(f['lower'],1/8);self.assertEqual(a['lower'],0);self.assertAlmostEqual(f['upper'],3/8);self.assertEqual(a['upper'],f['upper'])
    def test_common_lower_and_upper(self):
        p=[(0,1),(1,1),(3,0),(4,1)];d=[(1,1),(3,0),(3,0),(4,1)]
        for level in ['records','tables']:
            if level=='records':
                f=from_records(p,d,2,'common');a=from_horizon_records(4,administratively_censor(p,2),administratively_censor(d,2),2,'common')
            else:
                f=from_tables(4,event_table(p),event_table(d),2,2,3,'common');a=from_horizon_tables(4,[x for x in event_table(p) if x[0]<=2],[x for x in event_table(d) if x[0]<=2],2,2,3,'common')
            self.assertEqual(f['lower'],a['lower']);self.assertEqual(f['upper'],a['upper'])
    def test_boundary_status_suppression_is_explicit(self):
        p=[(0,1),(2,1),(4,0)];d=[(1,1),(2,1),(4,0)]
        for model in ['separate','common']:
            f=from_tables(3,event_table(p),event_table(d),2,2,2,model)
            h=from_terminal_tables(3,[x for x in event_table(p) if x[0]<2],[x for x in event_table(d) if x[0]<2],2,2,2,model)
            self.assertAlmostEqual(f['upper'],h['upper'])
        with self.assertRaises(ValueError):from_terminal_tables(3,event_table(p),[(1,3,1)],2,2,2)
    def test_post_horizon_rows_and_risk_conflicts_rejected(self):
        with self.assertRaises(ValueError):from_horizon_tables(3,[(3,3,1)],[(1,3,1)],2,2,2)
        with self.assertRaises(ValueError):from_horizon_tables(3,[(2,3,1)],[(1,3,1)],2,2,2)
    def test_no_fabricated_support_or_dropped_persons(self):
        with self.assertRaises(ValueError):administratively_censor([(1,0)],2)
        with self.assertRaises(ValueError):from_horizon_records(2,[(0,1)],[(0,1)],2)
        a=from_horizon_records(2,[(0,1),(1,1)],[(0,1),(1,1)],2)
        self.assertEqual(a['upper'],0)
    def test_terminal_record_group_and_zero_horizon(self):
        a=from_terminal_records(3,[(0,0),(1,1)],[(0,0),(1,1)],2)
        self.assertAlmostEqual(a['upper'],3/8)
        z=from_terminal_tables(2,[],[],0,2,2,'common');self.assertEqual(z['upper'],0)

if __name__=='__main__':unittest.main(verbosity=2)
