"""Reference QUANTIQUE : couplage local B-C et coherence reduite de A-B.

U_BC = |0><0|_B x I_C + |1><1|_B x Ry(2*theta)_C.
Ce mecanisme agit sur un singulet deja prepare ; il n'explique pas sa creation
par une dynamique classique. q complexe, phase fixe et perte de module distingues.
"""
from pathlib import Path
import hashlib
import json
import csv
import numpy as np

I=np.eye(2); X=np.array([[0.,1.],[1.,0.]])
Y=np.array([[0.,-1j],[1j,0.]]); Z=np.diag([1.,-1.])


def rho_q(q):
    if abs(q)>1+1e-12: raise ValueError('Abs(q) must not exceed 1')
    rho=np.diag([0.,.5,.5,0.]).astype(complex)
    rho[1,2]=-np.conj(q)/2; rho[2,1]=-q/2
    return rho


def reduce_state(psi,keep):
    order=list(keep)+[i for i in range(3) if i not in keep]
    matrix=psi.reshape(2,2,2).transpose(order).reshape(2**len(keep),-1)
    return matrix@matrix.conj().T


def interact(theta):
    rotation=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    u_bc=np.zeros((4,4)); u_bc[:2,:2]=I; u_bc[2:,2:]=rotation
    singlet=np.array([0.,1.,-1.,0.])/np.sqrt(2)
    psi=np.kron(singlet,np.array([1.,0.]))
    unitary=np.kron(I,u_bc)
    return unitary@psi,unitary,psi


def tensor(rho):
    return np.array([[np.trace(rho@np.kron(a,b)).real for b in (X,Y,Z)] for a in (X,Y,Z)])


def concurrence(rho):
    yy=np.kron(Y,Y)
    eig=np.linalg.eigvals(rho@yy@rho.conj()@yy).real
    roots=np.sqrt(np.maximum(eig,0)); roots.sort()
    return max(0.,float(roots[-1]-roots[:-1].sum()))


def twirl(rho):
    average=np.zeros((4,4),complex)
    for phi in np.arange(32)*2*np.pi/32:
        r=np.array([[np.cos(phi/2),-np.sin(phi/2)],[np.sin(phi/2),np.cos(phi/2)]])
        u=np.kron(r,r); average+=u@rho@u.T/32
    return average


def smax(rho):
    values=np.linalg.svd(tensor(rho),compute_uv=False)
    return float(2*np.sqrt(np.sum(values[:2]**2)))


def main():
    rows=[]
    for theta in np.linspace(0,np.pi/2,61):
        psi,u,initial=interact(theta)
        ab=reduce_state(psi,(0,1)); ac=reduce_state(psi,(0,2)); bc=reduce_state(psi,(1,2))
        c=np.cos(theta)
        np.testing.assert_allclose(ab,rho_q(c),atol=1e-13)
        np.testing.assert_allclose(reduce_state(psi,(0,)),I/2,atol=1e-13)
        np.testing.assert_allclose(reduce_state(psi,(1,)),I/2,atol=1e-13)
        np.testing.assert_allclose(u.conj().T@psi,initial,atol=1e-13)
        cab,cac,cbc=map(concurrence,(ab,ac,bc))
        assert abs(cab-c)<1e-7 and max(cac,cbc)<1e-7
        tau=4*np.linalg.det(reduce_state(psi,(0,))).real-cab*cab-cac*cac
        assert abs(tau-np.sin(theta)**2)<1e-7
        sf,st=smax(ab),smax(twirl(ab))
        np.testing.assert_allclose([sf,st],[2*np.sqrt(1+c*c),np.sqrt(2)*(1+c)],atol=1e-13)
        np.testing.assert_allclose(tensor(twirl(ab)),np.diag([-(1+c)/2,-c,-(1+c)/2]),atol=1e-13)
        rows.append([theta,c,cab,cac,cbc,tau,sf,st,np.sin(theta)])
    phases=[]
    for delta in (0.,np.pi/4,np.pi/3,np.pi/2,np.pi):
        q=np.exp(1j*delta); rho=rho_q(q); t=tensor(rho)
        np.testing.assert_allclose(t,[[ -q.real,q.imag,0],[-q.imag,-q.real,0],[0,0,-1]],atol=1e-13)
        assert abs(concurrence(rho)-1)<1e-7
        assert abs(smax(rho)-2*np.sqrt(2))<1e-13
        phases.append(dict(phase_rad=delta,Re_q=float(q.real),R=float(abs(q)),
                           concurrence=concurrence(rho),Smax=smax(rho)))
    out=Path(__file__).parent/'resultats_third_system'; out.mkdir(exist_ok=True)
    with (out/'coupling.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['theta','c','C_AB','C_AC','C_BC','three_tangle','S_fixed','S_twirl','D_C']); w.writerows(rows)
    report={'scope':__doc__,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'phase_cases':phases,'angles_checked':len(rows),
            'c_half':{'theta':float(np.pi/3),'S_twirl':float(1.5*np.sqrt(2)),
                      'S_fixed':float(np.sqrt(5)),'three_tangle':.75,
                      'gaussian_phase_std_rad':float(np.sqrt(2*np.log(2))),
                      'gaussian_phase_std_deg':float(np.rad2deg(np.sqrt(2*np.log(2))))}}
    (out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    data=np.asarray(rows); degrees=np.rad2deg(data[:,0])
    fig,axs=plt.subplots(1,2,figsize=(10.8,4.2),layout='constrained')
    axs[0].plot(degrees,data[:,1],label='c = C_AB = cos(theta)',color='#236891')
    axs[0].plot(degrees,data[:,8],label='Distinguabilite par C',color='#bd652f')
    axs[0].plot(degrees,data[:,5],label='Intrication tripartite (tangle)',color='#8e5a87')
    axs[0].set(ylabel='Coherence / information / tangle',title='Couplage local a un troisieme systeme')
    axs[0].legend(fontsize=8)
    axs[1].plot(degrees,data[:,6],label='Preparation a axe fixe',color='#bd652f')
    axs[1].plot(degrees,data[:,7],label='Preparation moyennee du papier',color='#236891')
    axs[1].axhline(2,color='.5',linestyle=':')
    axs[1].set(ylabel='CHSH maximal',title='Les predictions dependent de la preparation')
    axs[1].legend(fontsize=8)
    for ax in axs:
        ax.set_xlabel('Parametre de couplage theta (degres)')
        ax.grid(alpha=.15); ax.spines[['top','right']].set_visible(False)
    fig.savefig(out/'third_system.png',dpi=200); plt.close(fig)
    print(json.dumps(report['c_half'],indent=2)); print('61 partial traces, local marginals, reversal and monogamy checks passed.')


if __name__=='__main__': main()
