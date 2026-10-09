"""Small deterministic checks for reporting, not new simulation samples."""
from pathlib import Path
import json
import sys
import tempfile
import unittest
import numpy as np
root = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(root / 'src'), str(root / 'analysis')]
from method import optimality_check
from estimation import evaluate
from storage import Store
from summarize import export, original_support

class ReportingChecks(unittest.TestCase):
    def test_support_is_not_completion(self):
        p, d, tau = [(1, 1), (3, 0)], [(2, 1), (3, 0)], 2
        result = evaluate(p, d, tau)
        design = {'cells': [{'id': 0, 'replicates': 4, 'module': 'test'}],
                  'normal_quantile': 1.959963984540054,
                  'bootstrap_comparator_block': {'selected_cell_ids': []}}
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder); (path / 'analysis_state').mkdir()
            s = Store(path / 'analysis_state/checkpoints.sqlite', 'fixture', [(0, i) for i in range(4)])
            for i, status in enumerate(['COMPLETED', 'CALCULATION_FAILED', 'SUPPORT_FAIL']):
                s.start(0, i)
                payload = {'P_records': p, 'D_records': d, 'tau_internal': tau,
                           'result': result, 'support_valid': True}
                if status == 'SUPPORT_FAIL':
                    payload = {'P_records': [(1, 0), (1, 0)], 'D_records': d, 'tau_internal': tau}
                s.finish(0, i, status, payload)
            s.close()
            out = export(path, design, {'0': {'theta': .5}}, validation_run=True)
            den = out['cell_denominators'][0]
            self.assertEqual((den['planned'], den['attempted'], den['support_valid'], den['completed']), (4, 3, 2, 1))
            self.assertEqual((den['calculation_failed'], den['support_failed'], den['support_unknown_count']), (1, 1, 1))
            self.assertAlmostEqual(den['support_proportion'], 2/3)
            self.assertAlmostEqual(den['completion_proportion'], 1/3)
            for row in out['method_summaries']:
                self.assertEqual(row['method_completed'], 1)
                self.assertEqual(row['support_valid'], 2)
                self.assertIsNone(row['covered_planned_operational_yield'])
    def test_support_from_records(self):
        payload = {'P_records': [(1, 1), (3, 0)], 'D_records': [(2, 1), (3, 0)], 'tau_internal': 2}
        self.assertTrue(original_support(payload, 'CALCULATION_FAILED'))
        self.assertIsNone(original_support(None, 'NOT_ATTEMPTED'))
    def test_optimality_check_residual_allowance(self):
        cost = np.array([[1., 1.-5e-10], [1.-5e-10, 1.]])
        diag = optimality_check(cost, np.ones((2, 2), dtype=bool), np.array([0, 1]), 1)
        self.assertTrue(diag['numerical_status'].startswith('residual_checked'))
        self.assertGreaterEqual(diag['residual_objective_allowance'] + 1e-14, 2 - 2*(1-5e-10))
        self.assertLessEqual(diag['adjusted_dual_objective'], 2*(1-5e-10)+1e-14)

if __name__ == '__main__':
    unittest.main(verbosity=2)
