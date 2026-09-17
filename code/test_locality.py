"""Controles des modeles paralleles, independants du cosinus cible conjoint."""
import unittest
import numpy as np
import audit_locality as audit


class ParallelTests(unittest.TestCase):
    def test_malus_from_jones_projections(self):
        phi=np.arange(256)*2*np.pi/256
        field_a=np.array([np.cos(phi/2),np.sin(phi/2)])
        field_b=np.array([-np.sin(phi/2),np.cos(phi/2)])
        for a,b in ((0.,0.),(.7,-.2),(1.5,.9)):
            analyzer_a=np.array([np.cos(a/2),np.sin(a/2)])
            analyzer_b=np.array([np.cos(b/2),np.sin(b/2)])
            mean_a=2*(analyzer_a@field_a)**2-1
            mean_b=2*(analyzer_b@field_b)**2-1
            self.assertAlmostEqual(float(np.mean(mean_a*mean_b)),-.5*np.cos(a-b),places=13)

    def test_parallel_responses_and_pointwise_chsh(self):
        rng=np.random.default_rng(13); n=150000
        phi=rng.uniform(-np.pi,np.pi,n); ua=rng.random(n); ub=rng.random(n)
        for alice,bob in ((lambda a:audit.malus_alice(a,phi,ua),lambda b:audit.malus_bob(b,phi,ub)),
                          (lambda a:audit.threshold_alice(a,phi),lambda b:audit.threshold_bob(b,phi))):
            a0,a1=alice(0.),alice(np.pi/2)
            b0,b1=bob(np.pi/4),bob(-np.pi/4)
            score=a0*b0+a0*b1+a1*b0-a1*b1
            np.testing.assert_array_equal(np.abs(score),2)
        # Frequencies from routing follow independent Jones projection result.
        a,b=.37,-.62
        measured=np.mean(audit.malus_alice(a,phi,ua)*audit.malus_bob(b,phi,ub))
        self.assertLess(abs(measured+.5*np.cos(a-b)),6/np.sqrt(n))


if __name__=='__main__': unittest.main(verbosity=2)
