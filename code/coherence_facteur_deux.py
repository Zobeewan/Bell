"""Decomposition du facteur deux par coherence de paire, sans modele causal local.

Calcul par matrices de densite et projecteurs, puis comparaison aux expressions
analytiques. Les angles sont des angles de Bloch, comme dans Bell_test_48.
Execution : python analyse/coherence_facteur_deux.py (dependance : numpy).
"""
from pathlib import Path
import json
import numpy as np

I = np.eye(2)
X = np.array([[0., 1.], [1., 0.]])
Y = np.array([[0., -1j], [1j, 0.]])
Z = np.diag([1., -1.])


def rotation(phi):
    c, s = np.cos(phi / 2), np.sin(phi / 2)
    return np.array([[c, -s], [s, c]])


def state(phi, coherence):
    rho = np.zeros((4, 4))
    rho[1, 1] = rho[2, 2] = .5
    rho[1, 2] = rho[2, 1] = -.5 * coherence
    u = np.kron(rotation(phi), rotation(phi))
    return u @ rho @ u.T


def obs(a):
    return np.cos(a) * Z + np.sin(a) * X


def joint(rho, a, b):
    return np.array([[np.trace(rho @ np.kron((I+s*obs(a))/2,
                                            (I+t*obs(b))/2)).real
                      for t in (1, -1)] for s in (1, -1)])


def correlation(p):
    return float(p[0, 0] - p[0, 1] - p[1, 0] + p[1, 1])


def reduced(rho):
    tensor = rho.reshape(2, 2, 2, 2)
    return np.einsum('abcb->ac', tensor), np.einsum('abad->bd', tensor)


def classical_gaussian_fields():
    """Classical shared two-component noise: full versus connected intensities.

    z is a circular complex Gaussian vector with covariance I.
    Alice receives z, Bob receives J z*. This is a source prescription, not a
    claim of a passive optical implementation or a binary detection model.
    Wick's identity gives <I_s I_t>=1+(1-st*cos(a-b))/2.
    """
    rng = np.random.default_rng(20260915)
    z = (rng.normal(size=(600000, 2)) + 1j*rng.normal(size=(600000, 2)))/np.sqrt(2)
    j = np.array([[0., 1.], [-1., 0.]])
    bob = z.conj() @ j.T
    a, b = .21, 1.03
    ia = np.abs(z @ rotation(a))**2
    ib = np.abs(bob @ rotation(b))**2
    g2 = ia.T @ ib / len(z)
    independent = np.outer(ia.mean(axis=0), ib.mean(axis=0))
    connected = g2 - independent
    signs = np.array([[1,-1],[-1,1]])
    exact_connected = (np.ones((2,2))-signs*np.cos(a-b))/2
    exact_g2 = np.ones((2,2))+exact_connected
    assert np.max(np.abs(g2-exact_g2)) < .025
    assert np.max(np.abs(connected-exact_connected)) < .025
    return dict(source='E_A=z; E_B=J conjugate(z); circular complex Gaussian z',
                samples=len(z), a=a, b=b,
                raw_intensity_matrix=g2.tolist(),
                connected_intensity_matrix=connected.tolist(),
                raw_normalized_correlation=correlation(g2/g2.sum()),
                raw_exact=-np.cos(a-b)/3,
                connected_normalized_correlation=correlation(connected/connected.sum()),
                connected_exact=-np.cos(a-b),
                note='The connected quantity subtracts the product of mean intensities; it is not a raw binary-event correlation.')


def main():
    phis = np.arange(64) * (2*np.pi/64)
    angles = np.linspace(-np.pi, np.pi, 15)
    worst = 0.
    output = []
    for c in (0., .25, .5, .75, 1.):
        rho_avg = np.mean([state(phi, c) for phi in phis], axis=0)
        expected_rho = (np.eye(4) - (1+c)/2 * (np.kron(X, X)+np.kron(Z, Z))
                        - c*np.kron(Y, Y))/4
        assert np.allclose(rho_avg, expected_rho)
        assert np.linalg.eigvalsh(rho_avg).min() >= -1e-12
        assert np.isclose(np.trace(rho_avg), 1)
        for phi in (.0, .37, 1.21):
            rho = state(phi, c)
            assert np.linalg.eigvalsh(rho).min() >= -1e-12
            for marginal in reduced(rho):
                assert np.allclose(marginal, I/2)
            for a in angles:
                for b in angles:
                    p = joint(rho, a, b)
                    e = -np.cos(a-phi)*np.cos(b-phi)-c*np.sin(a-phi)*np.sin(b-phi)
                    worst = max(worst, abs(correlation(p)-e))
                    assert p.min() >= -1e-12
                    assert np.isclose(p.sum(), 1)
                    assert np.allclose(p.sum(axis=0), .5)
                    assert np.allclose(p.sum(axis=1), .5)
        for a in angles:
            for b in angles:
                p = joint(rho_avg, a, b)
                e = -(1+c)/2*np.cos(a-b)
                formula = (np.ones((2, 2)) + np.array([[1,-1],[-1,1]])*e)/4
                assert np.allclose(p, formula)
                worst = max(worst, abs(correlation(p)-e))
        t = [float(np.trace(rho_avg @ np.kron(q,q)).real) for q in (X,Y,Z)]
        output.append(dict(coherence=c, visibility=(1+c)/2,
                           correlation_tensor_diagonal_xyz=t,
                           aligned_joint_probabilities=joint(rho_avg, 0, 0).tolist(),
                           rho_eigenvalues=np.linalg.eigvalsh(rho_avg).tolist()))

    # Independently verify that random relative phases give the same density matrix.
    c, phi = .6, .41
    phases = (np.arccos(c), -np.arccos(c))
    u = np.kron(rotation(phi), rotation(phi))
    density_from_amplitudes = np.zeros((4,4), dtype=complex)
    for chi in phases:
        psi = u @ (np.array([0,1,-np.exp(1j*chi),0])/np.sqrt(2))
        density_from_amplitudes += np.outer(psi, psi.conj())/2
    assert np.allclose(density_from_amplitudes, state(phi, c))

    # A signed cross contribution redistributes probability; it is not a second event.
    phi, a, b = .2, 1.1, .7
    p0 = joint(state(phi, 0), a, b)
    p1 = joint(state(phi, 1), a, b)
    cross = p1-p0
    expected_cross = -np.sin(a-phi)*np.sin(b-phi)/4*np.array([[1,-1],[-1,1]])
    assert np.allclose(cross, expected_cross)
    assert np.allclose(cross.sum(axis=0), 0)
    assert np.allclose(cross.sum(axis=1), 0)

    result = dict(scope='Quantum coherence benchmark; not a local event-generating mechanism.',
                  formula='E_c(a,b)=-(1+c)/2*cos(a-b), after uniform planar rotation',
                  maximum_identity_error=worst, states=output,
                  interference_example=dict(phi=phi,a=a,b=b,
                    incoherent=p0.tolist(),coherent=p1.tolist(),cross=cross.tolist()),
                  classical_gaussian_fields=classical_gaussian_fields(),
                  all_checks_passed=True)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(maximum_identity_error=worst, all_checks_passed=True,
                         states=output, classical_gaussian_fields=result['classical_gaussian_fields']), indent=2))


if __name__ == '__main__':
    main()
