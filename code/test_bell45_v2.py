"""Controles independants de Bell 45 v2 : matrices, frequences et journal."""
from pathlib import Path
import importlib.util
import unittest
import numpy as np

src=Path(__file__).resolve().parents[1]/'simu'/'Bell_test_45.py'
spec=importlib.util.spec_from_file_location('bell45',src)
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

class CoherenceTests(unittest.TestCase):
    def test_density_matrix_reference(self):
        x=np.array([[0.,1.],[1.,0.]])
        z=np.diag([1.,-1.]); eye=np.eye(2)
        for c in (0.,.3,1.):
            for phi in (-.43,0.,1.21):
                r=np.array([[np.cos(phi/2),-np.sin(phi/2)],[np.sin(phi/2),np.cos(phi/2)]])
                rho=np.diag([0.,.5,.5,0.]); rho[1,2]=rho[2,1]=-c/2
                u=np.kron(r,r); rho=u@rho@u.T
                for a,b in ((0.,0.),(.7,-.2),(np.pi/2,.33),(2.3,-1.9)):
                    p,q=m.correlation_components(a,b,phi,c)
                    for s,t in m.OUTCOMES:
                        pa=(eye+s*(np.cos(a)*z+np.sin(a)*x))/2
                        pb=(eye+t*(np.cos(b)*z+np.sin(b)*x))/2
                        expected=float(np.trace(rho@np.kron(pa,pb)))
                        self.assertAlmostEqual(expected,(1+s*t*(p+q))/4,places=13)

    def test_quadrature_and_positive_probabilities(self):
        phi=np.arange(128)*2*np.pi/128
        for c in (0.,.2,.7,1.):
            for a in np.linspace(-np.pi,np.pi,11):
                for b in np.linspace(-np.pi,np.pi,9):
                    p,q=m.correlation_components(a,b,phi,c)
                    self.assertLessEqual(np.max(np.abs(p+q)),1+1e-13)
                    self.assertAlmostEqual(float(np.mean(p+q)),m.theoretical_E(a,b,c),places=13)
                    probs=m.joint_probabilities(a,b,c)
                    self.assertGreaterEqual(probs.min(),-1e-13)
                    self.assertAlmostEqual(float(probs.sum()),1.)
                    self.assertAlmostEqual(float(probs[0]+probs[1]),.5)
                    self.assertAlmostEqual(float(probs[0]+probs[2]),.5)

    def test_recorded_pairs_and_frequencies(self):
        n=160000
        for c in (0.,.4,1.):
            for a,b in ((.0,.0),(.0,np.pi/2),(.27,1.12)):
                ao,bo,*_=m.sample_joint_pairs(a,b,n,m.make_rng(222),c)
                self.assertEqual(set(np.unique(ao)),{-1,1})
                self.assertEqual(set(np.unique(bo)),{-1,1})
                result=m.summarize_pairs(ao,bo)
                self.assertEqual(sum(result['counts_pp_pm_mp_mm']),n)
                expected=m.joint_probabilities(a,b,c)
                for actual,p in zip(result['joint_pp_pm_mp_mm'],expected):
                    self.assertLessEqual(abs(actual-p),6*np.sqrt(p*(1-p)/n)+1/n)
                self.assertLess(abs(result['marginal_A_plus']-.5),3/np.sqrt(n))
                self.assertLess(abs(result['marginal_B_plus']-.5),3/np.sqrt(n))
                pp,pm,mp,mm=result['counts_pp_pm_mp_mm']
                self.assertAlmostEqual(result['E'],(pp+mm-pm-mp)/n)
                if c==1 and a==b: self.assertTrue(np.all(ao==-bo))

    def test_reproducibility_and_invalid_inputs(self):
        first=m.sample_joint_pairs(.7,-.2,1000,m.make_rng(7),.6)
        second=m.sample_joint_pairs(.7,-.2,1000,m.make_rng(7),.6)
        for a,b in zip(first,second): np.testing.assert_array_equal(a,b)
        for c in (-.01,1.01,np.nan,np.inf):
            with self.assertRaises(ValueError): m.validate_coherence(c)
        for n in (0,-1,1.5):
            with self.assertRaises(ValueError): m.sample_joint_pairs(0,0,n,m.make_rng(0))

if __name__=='__main__': unittest.main(verbosity=2)
