import unittest
from phiatomic import *
class TestAtomic(unittest.TestCase):
    def test_iron(self):
        self.assertEqual(orbit(complex(26/56)),(5,orbit(complex(26/56))[1],True))
    def test_real_boundary(self):
        self.assertFalse(orbit(.25+0j,1000)[2]); self.assertTrue(orbit(.26+0j,1000)[2])
    def test_collisions(self):
        self.assertEqual(parameters(6,6),parameters(4,4))
    def test_invalid(self):
        with self.assertRaises(ValueError): Nuclide(0,0,0)
        with self.assertRaises(ValueError): orbit(complex(float('nan')))
    def test_benchmark(self):
        rows=load(ROOT/'data/demo.csv'); a=benchmark(rows); self.assertEqual(a,benchmark(rows)); self.assertEqual(sum(f['test_n'] for f in a['folds']),len(rows)); self.assertTrue(all(math.isfinite(v) for v in a['pooled_mae_mev_per_nucleon'].values())); self.assertEqual(len(set.union(*(set(f['held_out_z']) for f in a['folds']))),len({r.z for r in rows}))
    def test_constant_features(self):
        self.assertAlmostEqual(fit_predict([[1],[1]],[2,4],[[1]])[0],3)
    def test_shell(self):
        self.assertEqual(shell_features(Nuclide(20,28,8)),[1,1])
if __name__=='__main__': unittest.main()
