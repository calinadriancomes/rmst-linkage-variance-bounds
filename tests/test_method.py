from pathlib import Path
import sys,itertools,unittest
from fractions import Fraction as F
import numpy as np
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root/"src"))
import observation_model as o
from method import scores,from_records,from_tables,upper_envelope_interval

class MethodChecks(unittest.TestCase):
    def test_strict_censor_boundary(self):
        self.assertFalse(o.ledger_edge((1,0),(1,1),"separate"))
        self.assertTrue(o.ledger_edge((1,1),(1,1),"separate"))
        self.assertTrue(o.ledger_edge((8,0),(1,0),"separate"))
        self.assertFalse(o.ledger_edge((1,0),(2,1),"common"))
    def test_exact_marginal_identity(self):
        for rec,tau in [([(1,1),(1,1),(3,0)],2),([(2,1),(3,0)],2),([(3,0),(3,0)],2),([(0,1),(3,0)],2),([(1,1),(1,1)],2)]:
            m,u,v=o.exact_scores(rec,tau);mm,uu,vv=scores(rec,tau)
            self.assertEqual(sum(x*x for x in u),v)
            self.assertAlmostEqual(float(m),mm,places=12)
            np.testing.assert_allclose([float(x) for x in u],uu,rtol=0,atol=1e-12)
            self.assertAlmostEqual(float(v),vv,places=12)
    def test_support_failures(self):
        for rec in [[(1,0)],[(1,1),(1,0)]]:
            with self.assertRaises(ValueError):scores(rec,2)
    def test_exact_endpoints_and_R2_api(self):
        fixtures=[([(1,1),(1,1),(3,0)],[(1,1),(2,1),(3,0)],2),
                  ([(1,0),(2,1),(3,0),(4,1)],[(2,1),(3,0),(4,1),(5,0)],3)]
        for p,d,tau in fixtures:
            _,up,_=o.exact_scores(p,tau);_,ud,_=o.exact_scores(d,tau)
            graph=o.dense_graph(p,d,"separate")
            lo,hi,count=o.exhaustive(up,ud,graph)
            r1=from_records(p,d,tau)
            self.assertAlmostEqual(r1["lower"],float(lo[0]),places=12)
            self.assertAlmostEqual(r1["upper"],float(hi[0]),places=12)
            for sign in [1,-1]:
                exact=o.exact_minflow(up,ud,graph,sign)
                self.assertEqual(exact[0],lo[0] if sign==1 else hi[0])
            pc=o.expand(o.cells_from_table(len(p),o.event_table(p),support=(tau,sum(t>=tau for t,e in p))))
            dc=o.expand(o.cells_from_table(len(d),o.event_table(d),support=(tau,sum(t>=tau for t,e in d))))
            _,u,_=o.exact_scores([(l,a) for a,l,h in pc],tau,False)
            _,v,_=o.exact_scores([(l,a) for a,l,h in dc],tau,False)
            low,high,_=o.exhaustive(u,v,o.dense_graph(pc,dc,"separate","R2"))
            r2=from_tables(len(p),o.event_table(p),o.event_table(d),tau,sum(t>=tau for t,e in p),sum(t>=tau for t,e in d))
            self.assertAlmostEqual(r2["lower"],float(low[0]),places=12)
            self.assertAlmostEqual(r2["upper"],float(high[0]),places=12)
            self.assertGreaterEqual(r1["lower"]+1e-12,r2["lower"])
            self.assertLessEqual(r1["upper"],r2["upper"]+1e-12)
    def test_infeasible_release(self):
        with self.assertRaises(ValueError):from_records([(2,0)],[(1,1)],1)
        with self.assertRaises(ValueError):from_tables(2,[(1,2,1)],[(1,2,1)],2,3,1)
    def test_interval_containment(self):
        a,b=upper_envelope_interval(1,4)
        self.assertLessEqual(a,1-1.959963984540054)
        self.assertGreaterEqual(b,1+1.959963984540054)
    def test_common_censoring_endpoints(self):
        p = [(1,1),(2,0),(1,1)]
        d = [(3,1),(2,0),(4,0)]
        tau = 2
        for level in ("R1", "R2"):
            if level == "R1":
                pp,dd=p,d
                result=from_records(p,d,tau,model="common")
            else:
                pc=o.expand(o.cells_from_table(3,o.event_table(p),support=(tau,1)))
                dc=o.expand(o.cells_from_table(3,o.event_table(d),support=(tau,3)))
                pp,dd=pc,dc
                result=from_tables(3,o.event_table(p),o.event_table(d),tau,1,3,model="common")
            rp=p if level=="R1" else [(l,a) for a,l,h in pp]
            rd=d if level=="R1" else [(l,a) for a,l,h in dd]
            _,u,_=o.exact_scores(rp,tau)
            _,v,_=o.exact_scores(rd,tau)
            lo,hi,n=o.exhaustive(u,v,o.dense_graph(pp,dd,"common",level))
            self.assertGreater(n,0)
            self.assertAlmostEqual(result["lower"],float(lo[0]),places=12)
            self.assertAlmostEqual(result["upper"],float(hi[0]),places=12)
    def test_invalid_public_arguments(self):
        p,d,tau=[(1,1),(3,0)],[(2,1),(3,0)],2
        with self.assertRaises(ValueError):
            from_records(p,d,tau,model="commmon")
        with self.assertRaises(ValueError):
            from_tables(2,o.event_table(p),o.event_table(d),tau,1,2,model="commmon")
        for z in (0,-1,float("inf"),float("nan")):
            with self.assertRaises(ValueError):
                upper_envelope_interval(1,1,z)
if __name__=="__main__":unittest.main()
