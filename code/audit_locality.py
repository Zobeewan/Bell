"""Source preparee avant choix : sampler conjoint et polariseurs PARALLELES.

Angles de Bloch a=2*alpha, b=2*beta. Sources d'orientations uniformes et
polarisations orthogonales. Pour Malus, les intensites normalisees servent
de probabilites de routage : cette regle de detection est explicite.
Le modele seuil est un controle local distinct, pas un polariseur de Malus.
"""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json
import numpy as np
import Bell_test_45 as bell


def malus_alice(a, phi, u):
    return np.where(u < np.cos((a-phi)/2)**2, 1, -1).astype(np.int8)


def malus_bob(b, phi, u):
    return np.where(u < np.sin((b-phi)/2)**2, 1, -1).astype(np.int8)


def threshold_alice(a, phi):
    return np.where(np.cos(a-phi) >= 0, 1, -1).astype(np.int8)


def threshold_bob(b, phi):
    return -threshold_alice(b, phi)


def model_responses(a, b, joint_source, optical_source):
    phi, ua, ub=optical_source
    return {
        'joint_c1': bell.respond_joint(a,b,joint_source,1.)[:2],
        'parallel_malus': (malus_alice(a,phi,ua),malus_bob(b,phi,ub)),
        'parallel_threshold': (threshold_alice(a,phi),threshold_bob(b,phi)),
    }


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pairs',type=int,default=1_000_000)
    parser.add_argument('--seed',type=int,default=20260917)
    parser.add_argument('--output-dir',type=Path,default=Path(__file__).parent/'resultats_locality')
    args=parser.parse_args(argv)
    if args.pairs < 1000: parser.error('--pairs doit etre au moins 1000')
    seeds=np.random.SeedSequence(args.seed).spawn(3)
    # Preparation complete BEFORE any setting is selected.
    joint_source=bell.prepare_source(args.pairs,np.random.default_rng(seeds[0]))
    optical_rng=np.random.default_rng(seeds[1])
    optical_source=(optical_rng.uniform(-np.pi,np.pi,args.pairs),
                    optical_rng.random(args.pairs),optical_rng.random(args.pairs))
    choices=np.random.default_rng(seeds[2])
    x=choices.integers(0,2,args.pairs); y=choices.integers(0,2,args.pairs)
    a=np.array([0.,np.pi/2])[x]; b=np.array([np.pi/4,-np.pi/4])[y]
    responses=model_responses(a,b,joint_source,optical_source)
    targets={'joint_c1':2*np.sqrt(2),'parallel_malus':np.sqrt(2),'parallel_threshold':2.}
    summary={}
    for name,(ao,bo) in responses.items():
        contexts=[]
        for i,j in ((0,0),(0,1),(1,0),(1,1)):
            mask=(x==i)&(y==j)
            context=bell.summarize_pairs(ao[mask],bo[mask]); context.update(x=i,y=j)
            contexts.append(context)
        s=abs(sum(k['E']*sign for k,sign in zip(contexts,(1,1,1,-1))))
        summary[name]={'S':s,'S_stderr':float(np.sqrt(sum(k['E_stderr']**2 for k in contexts))),
                       'expected_S':float(targets[name]),'contexts':contexts}
    # Interventions on EXACTLY the same source, including detector randomness.
    baseline=model_responses(0.,0.,joint_source,optical_source)
    change_a=model_responses(np.pi,0.,joint_source,optical_source)
    change_b=model_responses(0.,np.pi,joint_source,optical_source)
    for name,(ao,bo) in baseline.items():
        summary[name]['interventions']={
            'fraction_A_changed_by_beta':float(np.mean(ao!=change_b[name][0])),
            'fraction_B_changed_by_alpha':float(np.mean(bo!=change_a[name][1])),
            'fraction_product_changed_by_beta':float(np.mean(ao*bo!=change_b[name][0]*change_b[name][1]))}
    out=args.output_dir.resolve(); out.mkdir(parents=True,exist_ok=True)
    report={'seed':args.seed,'N_total':args.pairs,'source_prepared_before_choices':True,
            'angles':'Bloch; physical polarizers use half these angles',
            'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'bell45_sha256':hashlib.sha256(Path(bell.__file__).read_bytes()).hexdigest(),
            'numpy_version':np.__version__,'models':summary}
    (out/'audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    # Full run journal, no selection, all models evaluated on the chosen context.
    with gzip.open(out/'journal.csv.gz','wt',encoding='utf-8',newline='') as f:
        writer=csv.writer(f)
        writer.writerow(['trial','x','y','a','b']+[f'{name}_{side}' for name in responses for side in ('A','B')])
        writer.writerows(zip(range(args.pairs),x,y,a,b,*(v for ab in responses.values() for v in ab)))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    delta=np.linspace(-np.pi,np.pi,97)
    # Analytic curves and fresh evaluations using a fixed prepared population.
    selected_joint=bell.PreparedSource(*(v[:80000] for v in joint_source))
    selected_optical=tuple(v[:80000] for v in optical_source)
    curves={name:[] for name in responses}
    for d in delta:
        for name,(ao,bo) in model_responses(0.,-d,selected_joint,selected_optical).items():
            curves[name].append(float(np.mean(ao*bo)))
    fig,axs=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
    labels=['Conjoint c=1','Malus parallele','Seuil local']
    colors=['#236b8e','#bd6a30','#52947f']
    analytic=[-np.cos(delta),-.5*np.cos(delta),-1+2*np.abs(delta)/np.pi]
    for name,label,color,theory in zip(responses,labels,colors,analytic):
        axs[0].plot(delta,theory,color=color,label=label)
        axs[0].plot(delta[::4],np.array(curves[name])[::4],'.',color=color,markersize=4)
    axs[0].set(xlabel='a - b (angles de Bloch, rad)',ylabel='Moyenne AB',ylim=(-1.08,1.08))
    axs[0].legend(fontsize=8); axs[0].grid(alpha=.18)
    axs[1].bar(labels,[summary[k]['S'] for k in responses],color=colors,
               yerr=[1.96*summary[k]['S_stderr'] for k in responses],capsize=4)
    axs[1].axhline(2,color='#555',linestyle='--',label='Borne CHSH locale, MI')
    axs[1].axhline(2*np.sqrt(2),color='#236b8e',linestyle=':',label='2 sqrt(2)')
    axs[1].set(ylabel='|S|, choix apres preparation',ylim=(0,3.15))
    axs[1].tick_params(axis='x',labelsize=8); axs[1].legend(fontsize=8,loc='upper right')
    fig.savefig(out/'comparison.png',dpi=200); plt.close(fig)
    with (out/'curves.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.writer(f); w.writerow(['delta',*responses]); w.writerows(zip(delta,*curves.values()))
    print(json.dumps({name:{k:v for k,v in result.items() if k!='contexts'} for name,result in summary.items()},indent=2))


if __name__=='__main__': main()
