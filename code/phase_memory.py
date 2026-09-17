"""Memoire classique de phase : preparation AB, stockage, bruit, contact BC.

Oscillateurs reduits en contact pendant les couplages. Pas de champ spatial,
pas de detecteur binaire, pas de conversion automatique de Q en correlation Bell.
Q_AB = moyenne d'ensemble de exp(i*(theta_B-theta_A)).
Unites de temps arbitraires ; frequences en rad/unite, D en rad^2/unite.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import numpy as np


def order(theta):
    return np.array([np.exp(1j*(theta[j]-theta[i])).mean() for i,j in ((0,1),(1,2),(0,2))])


def prepare(n,dt,duration,coupling,rng):
    theta=rng.uniform(-np.pi,np.pi,(3,n))
    for _ in range(round(duration/dt)):
        interaction=coupling*np.sin(theta[1]-theta[0])
        theta[0]+=dt*interaction
        theta[1]-=dt*interaction
    return theta


def evolve(initial,dt,duration,mode,*,noise=.05,detuning=.8,spread=.35,
           coupling=2.,contact_start=4.,contact_stop=8.,seed=17,save_every=10,snapshots=None):
    """All scenarios start from the same prepared population. Euler-Maruyama.

    In 'third', only B and C interact during their contact; A is isolated.
    'shared' retains AB coupling to show classical synchrony can be shared.
    'echo' reverses each stored detuning at duration/2 without inter-site coupling.
    """
    valid={'free','drift','diffusion','spread','echo','third','shared'}
    if mode not in valid: raise ValueError('Unknown scenario')
    theta=initial.copy(); rng=np.random.default_rng(seed)
    individual_detuning=rng.normal(0,spread,theta.shape[1])
    times=[0.]; values=[order(theta)]
    if snapshots is not None: snapshots.append(theta[:,:48].copy())
    steps=round(duration/dt)
    for k in range(steps):
        t=k*dt
        drift=np.zeros_like(theta)
        if mode=='drift': drift[1]=detuning
        if mode in ('spread','echo'):
            drift[1]=individual_detuning*(-1 if mode=='echo' and t>=duration/2 else 1)
        if mode=='shared':
            force=coupling*np.sin(theta[1]-theta[0])
            drift[0]+=force; drift[1]-=force
        if mode in ('third','shared') and contact_start<=t<contact_stop:
            force=coupling*np.sin(theta[2]-theta[1])
            drift[1]+=force; drift[2]-=force
        theta+=dt*drift
        if mode=='diffusion':
            theta[:2]+=np.sqrt(2*noise*dt)*rng.standard_normal(theta[:2].shape)
        if (k+1)%save_every==0 or k+1==steps:
            times.append((k+1)*dt); values.append(order(theta))
            if snapshots is not None: snapshots.append(theta[:,:48].copy())
    return np.asarray(times),np.asarray(values),theta


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replicas',type=int,default=4096)
    parser.add_argument('--dt',type=float,default=.01)
    parser.add_argument('--duration',type=float,default=12.)
    parser.add_argument('--prep-time',type=float,default=6.)
    parser.add_argument('--noise',type=float,default=.05)
    parser.add_argument('--detuning',type=float,default=.8)
    parser.add_argument('--spread',type=float,default=.35)
    parser.add_argument('--coupling',type=float,default=2.)
    parser.add_argument('--contact-start',type=float,default=4.)
    parser.add_argument('--contact-stop',type=float,default=8.)
    parser.add_argument('--seed',type=int,default=20260917)
    parser.add_argument('--output-dir',type=Path,default=Path(__file__).parent/'resultats_phase_memory')
    args=parser.parse_args(argv)
    numbers=[args.dt,args.duration,args.prep_time,args.noise,args.detuning,args.spread,args.coupling,args.contact_start,args.contact_stop]
    if not np.all(np.isfinite(numbers)): parser.error('Parametres finis requis')
    if args.replicas<100 or args.dt<=0 or args.duration<=0 or args.prep_time<=0 or min(args.noise,args.spread,args.coupling)<0:
        parser.error('Parametres hors domaine')
    if not 0<=args.contact_start<args.contact_stop<=args.duration:
        parser.error('Contact requis dans la duree de stockage')
    # Resolve all switches exactly on the integration grid.
    for t in (args.duration,args.prep_time,args.contact_start,args.contact_stop,args.duration/2):
        if not np.isclose(t/args.dt,round(t/args.dt)): parser.error('Les temps doivent etre des multiples de dt (demi-duree incluse)')
    if 4*args.coupling*args.dt>.2: parser.error('dt trop grand pour le couplage ; reduire dt')
    initial=prepare(args.replicas,args.dt,args.prep_time,args.coupling,np.random.default_rng(args.seed))
    configs={k:getattr(args,k) for k in ('noise','detuning','spread','coupling','contact_start','contact_stop')}
    modes=('free','drift','diffusion','spread','echo','third','shared')
    snapshots={m:[] for m in modes}
    runs={m:evolve(initial,args.dt,args.duration,m,**configs,seed=args.seed+1,snapshots=snapshots[m]) for m in modes}
    out=args.output_dir.resolve(); out.mkdir(parents=True,exist_ok=True)
    with (out/'trajectories.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['scenario','time','pair','Q_real','Q_imag','R','delta'])
        for mode,(times,qs,_) in runs.items():
            for t,row in zip(times,qs):
                for pair,q in zip(('AB','BC','AC'),row): w.writerow((mode,t,pair,q.real,q.imag,abs(q),np.angle(q)))
    # Retain the prepared phase population and each endpoint for audit.
    np.savez_compressed(out/'states.npz',prepared=initial,**{m:v[2] for m,v in runs.items()})
    visual_indices=np.unique(np.r_[np.arange(0,len(runs['free'][0]),2),len(runs['free'][0])-1])
    visual={'times':np.round(runs['free'][0][visual_indices],5).tolist(),'replicas':args.replicas,
            'contact':[args.contact_start,args.contact_stop],'modes':{}}
    for mode in modes:
        qs=runs[mode][1][visual_indices]
        samples=np.asarray(snapshots[mode])[visual_indices]
        differences=np.stack([samples[:,j]-samples[:,i] for i,j in ((0,1),(1,2),(0,2))],axis=1)
        differences=np.angle(np.exp(1j*differences))
        visual['modes'][mode]={'q':np.round(np.stack([qs.real,qs.imag],axis=-1),7).tolist(),
                               'phases':np.round(differences,4).tolist()}
    (out/'interactive_data.json').write_text(json.dumps(visual,separators=(',',':')),encoding='utf-8')
    q0=order(initial)[0]
    report={'scope':__doc__,'parameters':{k:v for k,v in vars(args).items() if k!='output_dir'},
            'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'numpy_version':np.__version__,'initial_R_AB':float(abs(q0)),
            'endpoints':{m:{pair:{'R':float(abs(q)),'real':float(q.real),'imag':float(q.imag)}
                            for pair,q in zip(('AB','BC','AC'),v[1][-1])} for m,v in runs.items()},
            'analytical_checks':{}}
    for mode in ('free','drift','diffusion','spread','echo'):
        times,qs,_=runs[mode]
        theory={'free':np.ones_like(times),'drift':np.exp(1j*args.detuning*times),
                'diffusion':np.exp(-2*args.noise*times),
                'spread':np.exp(-.5*(args.spread*times)**2),
                'echo':np.exp(-.5*(args.spread*np.minimum(times,args.duration-times))**2)}[mode]*q0
        error=float(np.max(np.abs(qs[:,0]-theory)))
        tolerance=1e-10 if mode in ('free','drift') else 6/np.sqrt(args.replicas)
        # Echo envelope is an ensemble expectation; exact refocusing is checked separately.
        assert error<tolerance,(mode,error,tolerance)
        report['analytical_checks'][mode]={'max_complex_error':error,'tolerance':tolerance}
    assert np.max(np.abs(runs['third'][2][0]-initial[0]))==0
    assert abs(runs['echo'][1][-1,0]-q0)<1e-10
    # Deterministic integration refinement for the contact model.
    fine_times,q_fine,_=evolve(initial,args.dt/2,args.duration,'third',**configs,seed=args.seed+1,save_every=20)
    np.testing.assert_allclose(fine_times,runs['third'][0],atol=1e-12)
    contact_error=float(np.max(np.abs(q_fine-runs['third'][1])))
    assert contact_error<.005,contact_error
    report['contact_dt_refinement_error']=contact_error
    (out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,2,figsize=(11,7.6),layout='constrained')
    times=runs['free'][0]
    colors={'free':'#497760','drift':'#9b4f7b','diffusion':'#bd652f','spread':'#7a768f','echo':'#47868b','third':'#236891','shared':'#497760'}
    for mode,label in (('free','Sans perturbation'),('drift','Decalage fixe de frequence'),('diffusion','Diffusion independante'),('spread','Dispersion des frequences')):
        axs[0,0].plot(times,np.abs(runs[mode][1][:,0]),color=colors[mode],label=label,linestyle='--' if mode=='drift' else '-')
    axs[0,0].set(title='A  Coherence de phase de la paire',ylabel='R_AB (moyenne d’ensemble)')
    axs[0,0].legend(fontsize=8)
    axs[0,1].plot(times,runs['drift'][1][:,0].real,color=colors['drift'],label='Partie reelle de Q_AB')
    axs[0,1].plot(times,np.abs(runs['drift'][1][:,0]),'--',color=colors['free'],label='Module R_AB')
    axs[0,1].set(title='B  Un dephasage tourne sans brouiller',ylabel='Q_AB reel / module')
    axs[0,1].legend(fontsize=8)
    for i,(pair,color) in enumerate(zip(('AB','BC','AC'),('#236891','#bd652f','#9b4f7b'))):
        axs[1,0].plot(times,np.abs(runs['third'][1][:,i]),color=color,label=pair)
    axs[1,0].axvspan(args.contact_start,args.contact_stop,color='.6',alpha=.12)
    axs[1,0].set(title='C  Contact B-C apres separation A-B',ylabel='Coherence R de chaque paire')
    axs[1,0].legend(fontsize=8)
    axs[1,1].plot(times,np.abs(runs['spread'][1][:,0]),color=colors['spread'],label='Dispersion sans echo')
    axs[1,1].plot(times,np.abs(runs['echo'][1][:,0]),color=colors['echo'],label='Dispersion puis inversion')
    axs[1,1].plot(times,np.abs(runs['shared'][1][:,0]),color=colors['shared'],linestyle='--',label='Contact BC, liaison AB maintenue')
    axs[1,1].axvline(args.duration/2,color='.6',linestyle=':')
    axs[1,1].set(title='D  Memoire recuperable et partage classique',ylabel='R_AB')
    axs[1,1].legend(fontsize=8)
    for ax in axs.flat:
        ax.set_xlabel('Temps apres preparation (unites arbitraires)')
        ax.grid(alpha=.15); ax.spines[['top','right']].set_visible(False)
    fig.savefig(out/'phase_memory.png',dpi=200); plt.close(fig)
    print(json.dumps({'initial_R_AB':report['initial_R_AB'],
          'final_R_AB':{m:v['AB']['R'] for m,v in report['endpoints'].items()},
          'contact_dt_refinement_error':contact_error},indent=2))


if __name__=='__main__': main()
