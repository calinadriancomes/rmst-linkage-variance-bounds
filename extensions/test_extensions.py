"""Tests of the data-only terminal example and its declared input limitations."""
import csv,json,sys,tempfile,unittest
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from analyst_from_tables import load_inputs,analyze

class TerminalExample(unittest.TestCase):
 def setUp(self):
  self.source=HERE/'rotterdam_terminal_example'
 def make(self,changes=None,field=None,replace=None):
  t=tempfile.TemporaryDirectory();self.addCleanup(t.cleanup);p=Path(t.name)
  for f in ('metadata.json','P_event_risk.csv','D_event_risk.csv'):(p/f).write_bytes((self.source/f).read_bytes())
  if changes:
   m=json.loads((p/'metadata.json').read_text());m.update(changes);(p/'metadata.json').write_text(json.dumps(m))
  if field:
   f=p/f'{field}_event_risk.csv';s=f.read_text();f.write_text(replace(s))
  return p
 def test_source_dimensions(self):
  m,tp,td=load_inputs(self.source)
  self.assertEqual((m['n'],m['P_risk_tau'],m['D_risk_tau']), (2982,1571,2084))
  self.assertEqual((len(tp),len(td)),(878,600))
 def test_numeric_reference(self):
  out=analyze(self.source,True)
  self.assertAlmostEqual(out['upper'],105.01687556404875,places=8)
  self.assertAlmostEqual(out['interval_lower'],214.20322028324568,places=8)
  self.assertAlmostEqual(out['interval_upper'],254.37375689690083,places=8)
  self.assertLessEqual(out['optimality']['upper']['gap'],out['optimality']['upper']['tolerance'])
 def test_missing_checkpoint(self):
  with self.assertRaises(ValueError):analyze(self.make({'P_risk_tau':None}))
 def test_incorrect_checkpoint(self):
  with self.assertRaises(ValueError):analyze(self.make({'P_risk_tau':2982}),False)
 def test_invalid_cohort_size(self):
  with self.assertRaises(ValueError):analyze(self.make({'n':100}),False)
 def test_posthorizon_row(self):
  with self.assertRaises(ValueError):load_inputs(self.make(field='D',replace=lambda s:s+'1827,1,1\n'))
 def test_horizon_boundary_row(self):
  with self.assertRaises(ValueError):load_inputs(self.make(field='D',replace=lambda s:s+'1826.25,1,1\n'))
 def test_noninteger_counts(self):
  with self.assertRaises(ValueError):load_inputs(self.make(field='P',replace=lambda s:s.replace('38.0,2981,1','38.0,2981,1.5',1)))
 def test_negative_censoring_cell(self):
  with self.assertRaises(ValueError):analyze(self.make(field='P',replace=lambda s:s.replace('38.0,2981,1','38.0,2,1',1)),False)
 def test_missing_matching(self):
  # A three-patient toy example with incompatible ordered endpoints.
  from horizon_release import from_terminal_tables
  with self.assertRaises(ValueError):from_terminal_tables(3,[(0.,3,1)],[(0.,3,3)],1.,0,0,model='separate',check_optimality=False)

if __name__=='__main__':unittest.main(verbosity=2)
