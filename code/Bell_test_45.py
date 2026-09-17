"""Bell 45 v2 - coherence de paire et frequences conjointes.
Loi : E_c(a,b)=-(1+c)/2*cos(a-b), c dans [0,1], c=1 par defaut.
Les paires sont echantillonnees CONJOINTEMENT avec les deux reglages.
Ce programme valide la loi statistique, pas une dynamique de detecteurs locaux.
Usage : python Bell_test_45.py --coherence 1 --no-show
Dependances : numpy, matplotlib. Original : Bell_test_45_original.py.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import platform
from typing import NamedTuple
import numpy as np

VERSION = "2.1"
OUTCOMES = ((1,1),(1,-1),(-1,1),(-1,-1))

def validate_coherence(c):
    c=float(c)
    if not np.isfinite(c) or not 0 <= c <= 1:
        raise ValueError("La coherence doit etre finie et comprise entre 0 et 1.")
    return c

def make_rng(seed):
    return np.random.default_rng(seed)

def hidden_vars(N_pairs,rng):
    """Reference originale : phase partagee et tirages independants."""
    phi=rng.uniform(-np.pi,np.pi,N_pairs)
    return phi,phi+np.pi,rng.random(N_pairs),rng.random(N_pairs)

def compute_branches(alpha,beta,phi_g,phi_d,rAs,rBs):
    ag=np.where(rAs<np.cos(.5*(alpha-phi_g))**2,1,-1)
    ad=np.where(rAs<np.cos(.5*(alpha-phi_d))**2,1,-1)
    bg=np.where(rBs<np.cos(.5*(beta-phi_g))**2,1,-1)
    bd=np.where(rBs<np.cos(.5*(beta-phi_d))**2,1,-1)
    return ag,ad,bg,bd

def correlation_components(alpha,beta,phi,coherence=1.):
    """Contributions a une esperance, pas deux resultats AB par paire."""
    c=validate_coherence(coherence)
    return -np.cos(alpha-phi)*np.cos(beta-phi),-c*np.sin(alpha-phi)*np.sin(beta-phi)

def theoretical_E(alpha,beta,coherence=1.):
    return -(1+validate_coherence(coherence))/2*np.cos(np.asarray(alpha)-beta)

def joint_probabilities(alpha,beta,coherence=1.):
    """Ordre ++,+-,-+,-- ; orientation phi moyennee."""
    e=theoretical_E(alpha,beta,coherence)
    return np.stack([(1+s*t*e)/4 for s,t in OUTCOMES],axis=-1)

class PreparedSource(NamedTuple):
    phi: np.ndarray
    alice_sign: np.ndarray
    agreement_uniform: np.ndarray


def prepare_source(N_pairs,rng):
    """Prepare tous les aleas sans recevoir aucun reglage.

    Etat statistique (phi, R, U), pas un modele de particule ou de champ.
    Les reglages peuvent etre tires APRES cet appel avec un autre generateur.
    """
    if not isinstance(N_pairs,(int,np.integer)) or N_pairs<1:
        raise ValueError("N_pairs doit etre un entier strictement positif.")
    phi=rng.uniform(-np.pi,np.pi,N_pairs)
    r=np.where(rng.random(N_pairs)<.5,1,-1).astype(np.int8)
    return PreparedSource(phi,r,rng.random(N_pairs))


def respond_joint(alpha,beta,source,coherence=1.):
    """Reponse deterministe a l'etat prepare ; acces EXPLICITE aux deux angles.

    A=R ; B=R si U < (1+E_phi(alpha,beta))/2, sinon B=-R.
    Preparer les aleas en avance assure leur independance des choix, mais
    ne retire pas la dependance de B au reglage alpha. Aucun champ n'evolue.
    Accepte aussi des tableaux de reglages, un choix par essai.
    """
    c=validate_coherence(coherence)
    phi,a,u=source
    pop,cross=correlation_components(alpha,beta,phi,c)
    e=pop+cross
    if np.any(~np.isfinite(e)) or np.any(np.abs(e)>1+1e-12):
        raise ValueError("Esperance conditionnelle invalide.")
    agree=u<np.clip((1+e)/2,0,1)
    b=np.where(agree,a,-a).astype(np.int8)
    return a,b,phi,pop,cross


def sample_joint_pairs(alpha,beta,N_pairs,rng,coherence=1.):
    """Echantillonneur conjoint historique, equivalent a preparation + reponse.

    Un couple AB par essai. L'ordre des nombres aleatoires est conserve depuis
    la v2.0 pour rendre les resultats precedents exactement reproductibles.
    L'interface respond_joint permet de choisir les angles apres preparation.
    """
    validate_coherence(coherence)
    return respond_joint(alpha,beta,prepare_source(N_pairs,rng),coherence)

def summarize_pairs(a,b):
    n=len(a); e=float(np.mean(a*b))
    counts=[int(np.count_nonzero((a==s)&(b==t))) for s,t in OUTCOMES]
    return dict(N=n,counts_pp_pm_mp_mm=counts,joint_pp_pm_mp_mm=[v/n for v in counts],
                E=e,E_stderr=float(np.sqrt(max(0,1-e*e)/n)),
                marginal_A_plus=float(np.mean(a==1)),marginal_B_plus=float(np.mean(b==1)))

def E_indistinguishable(alpha,beta,N_pairs=500_000,seed=None,coherence=1.):
    """Nom historique ; retourne maintenant la moyenne des produits AB."""
    a,b,*_=sample_joint_pairs(alpha,beta,N_pairs,make_rng(seed),coherence)
    return float(np.mean(a*b))

def E_components_detail(alpha,beta,N_pairs=100_000,seed=None,coherence=1.):
    phi=make_rng(seed).uniform(-np.pi,np.pi,N_pairs)
    p,q=correlation_components(alpha,beta,phi,coherence)
    return float(p.mean()),float(q.mean())

def simulate_E_of_delta(deltas,N_pairs=200_000,seed=2,coherence=1.):
    rng=make_rng(seed); values=[]
    for d in deltas:
        a,b,*_=sample_joint_pairs(0.,-float(d),N_pairs,rng,coherence)
        values.append(float(np.mean(a*b)))
    return np.array(values)

def rmse(a,b):
    return float(np.sqrt(np.mean((np.asarray(a)-b)**2)))

def convergence_test(deltas,N_list,seed=5,coherence=1.):
    target=theoretical_E(0.,-np.asarray(deltas),coherence)
    seeds=np.random.SeedSequence(seed).spawn(len(N_list))
    errors=[rmse(simulate_E_of_delta(deltas,int(n),s,coherence),target) for n,s in zip(N_list,seeds)]
    return np.asarray(N_list),np.asarray(errors)

def write_csv(path,header,rows):
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.writer(f); w.writerow(header); w.writerows(rows)

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coherence','--c',type=float,default=1.)
    parser.add_argument('--pairs',type=int,default=200_000,help='Paires par reglage')
    parser.add_argument('--seed',type=int,default=20260916)
    parser.add_argument('--no-show',action='store_true')
    parser.add_argument('--quick',action='store_true',help='Grilles reduites')
    parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent/'resultats_bell45_v2')
    args=parser.parse_args(argv); c=validate_coherence(args.coherence)
    if args.pairs<100: parser.error('--pairs doit etre au moins 100')
    import matplotlib
    if args.no_show: matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'savefig.dpi':200,'figure.dpi':110,'font.family':'DejaVu Sans'})
    out=args.output_dir.resolve(); out.mkdir(parents=True,exist_ok=True)
    rng=make_rng(args.seed); n=args.pairs
    deltas=np.linspace(-np.pi,np.pi,37 if args.quick else 73)
    target=theoretical_E(0.,-deltas,c); singlet=-np.cos(deltas)
    curve=[]; local=[]; population=[]; coherence=[]
    for d in deltas:
        a,b,phi,p,q=sample_joint_pairs(0.,-d,n,rng,c)
        curve.append(summarize_pairs(a,b)); population.append(float(p.mean())); coherence.append(float(q.mean()))
        ag,ad,bg,bd=compute_branches(0.,-d,*hidden_vars(n,rng))
        local.append(float(np.mean(ag*bd)))
    e_sim=np.array([r['E'] for r in curve]); e_se=np.array([r['E_stderr'] for r in curve])
    population=np.array(population); coherence=np.array(coherence)
    chsh=[]; journal=[]
    for i,(a,b) in enumerate(((0,np.pi/4),(0,-np.pi/4),(np.pi/2,np.pi/4),(np.pi/2,-np.pi/4))):
        ao,bo,phi,p,q=sample_joint_pairs(a,b,n,rng,c)
        r=summarize_pairs(ao,bo); r.update(alpha=a,beta=b,expected_E=float(theoretical_E(a,b,c))); chsh.append(r)
        for j in range(min(1000,n)):
            journal.append((i,j,a,b,phi[j],int(ao[j]),int(bo[j]),int(ao[j]*bo[j])))
    s_signed=sum(r['E']*sign for r,sign in zip(chsh,(1,1,1,-1)))
    s_abs=abs(s_signed); s_se=float(np.sqrt(sum(r['E_stderr']**2 for r in chsh)))
    s_exact=float(np.sqrt(2)*(1+c))
    ns=[625,1250,2500,5000,10000,20000,40000] if args.quick else [625,1250,2500,5000,10000,20000,40000,80000,160000,320000]
    ns,errors=convergence_test(deltas,ns,args.seed+1,c)
    fit,cov=np.polyfit(np.log(ns),np.log(errors),1,cov=True)
    slope,intercept=map(float,fit); slope_se=float(np.sqrt(cov[0,0]))
    expected_rmse=np.sqrt(np.mean(1-target**2)/ns)
    figures=[]
    def save(fig,name):
        fig.tight_layout(); fig.savefig(out/name); figures.append(name)

    fig,ax=plt.subplots(figsize=(8.5,4.8))
    ax.plot(deltas,singlet,color='#1f4b77',lw=2,label=r'Singulet : $-\cos\Delta$')
    if c!=1: ax.plot(deltas,target,'--',color='#cc6b32',label=f'Loi de paire, c = {c:g}')
    ax.errorbar(deltas,e_sim,yerr=1.96*e_se,fmt='o',ms=3,color='#b33845',label=f'Tirages conjoints, c = {c:g} (IC 95 %)')
    ax.plot(deltas,local,'--',color='#408573',lw=1.5,label=r'Branche indépendante : $-\cos\Delta/2$')
    ax.set(xlabel=r'$\Delta=a-b$ (rad)',ylabel=r'$\langle AB\rangle$',ylim=(-1.1,1.1),title='1  Corrélations des paires enregistrées')
    ax.grid(alpha=.2); ax.legend(fontsize=9); save(fig,'fig1_E_vs_delta.png')

    fig,ax=plt.subplots(figsize=(8.5,4.8))
    ax.plot(deltas,population,'-o',markevery=4,ms=4,color='#1f4b77',label=r'Populations : moyenne de $-\cos(a-\phi)\cos(b-\phi)$')
    ax.plot(deltas,coherence,'--x',markevery=(2,4),ms=5,color='#cc6b32',label=r'Cohérence : moyenne de $-c\sin(a-\phi)\sin(b-\phi)$')
    ax.plot(deltas,population+coherence,color='#34363a',lw=2,label='Somme des contributions à l’espérance')
    ax.plot(deltas,target,':',color='#408573',lw=2,label='Loi intégrée')
    ax.set(xlabel=r'$\Delta$ (rad)',ylabel='Contribution',title='2  Populations et cohérence croisée')
    ax.grid(alpha=.2); ax.legend(fontsize=8.5); save(fig,'fig2_branches.png')

    fig,ax=plt.subplots(figsize=(8.5,4.8))
    vals=[2.,s_abs,s_exact,2*np.sqrt(2)]
    bars=ax.bar(['Borne locale\nCHSH','Tirages\nconjoints',f'Loi c = {c:g}\n(théorie)','Singulet\n(théorie)'],vals,
                color=['#b1b7bd','#b33845','#cc6b32','#1f4b77'],width=.65)
    ax.errorbar(1,s_abs,yerr=1.96*s_se,color='black',capsize=5,fmt='none')
    for bar,val in zip(bars,vals): ax.text(bar.get_x()+bar.get_width()/2,val+.08,f'{val:.4f}',ha='center',fontsize=10)
    ax.axhline(2,color='#7e858c',ls='--',lw=1)
    ax.set(ylabel=r'$|S|$',ylim=(0,3.35),title='3  CHSH pour la loi conjointe échantillonnée')
    ax.grid(axis='y',alpha=.2); save(fig,'fig3_CHSH.png')

    fig,ax=plt.subplots(figsize=(8.5,4.8))
    ax.loglog(ns,errors,'o',color='#b33845',label='RMSE des tirages')
    ax.loglog(ns,np.exp(intercept)*ns**slope,'--',color='#1f4b77',label=f'Pente ajustée {slope:.3f} ± {slope_se:.3f}')
    ax.loglog(ns,expected_rmse,':',color='#408573',lw=2,label='Racine de l’erreur quadratique moyenne attendue')
    ax.set(xlabel='N (paires par réglage)',ylabel='RMSE par rapport à la loi pour c',title='4  Convergence statistique')
    ax.grid(which='both',alpha=.2); ax.legend(fontsize=8.5); save(fig,'fig4_convergence.png')

    fig,(ax,bx)=plt.subplots(1,2,figsize=(11,4.5))
    probs=np.array([r['joint_pp_pm_mp_mm'] for r in curve]); theory=joint_probabilities(0.,-deltas,c)
    for j,(label,marker,color) in enumerate(zip(('++','+-','-+','--'),('o','s','x','+'),('#1f4b77','#b33845','#cc6b32','#408573'))):
        ax.plot(deltas,theory[:,j],color=color,lw=1,alpha=.65)
        ax.plot(deltas[::3],probs[::3,j],marker,ms=4,color=color,label=label)
    ax.set(xlabel=r'$\Delta$ (rad)',ylabel='Fréquence conjointe',ylim=(-.03,.55),title='5a  Les quatre issues')
    ax.legend(ncol=2,fontsize=9)
    for key,label,color in (('marginal_A_plus','Alice +','#1f4b77'),('marginal_B_plus','Bob +','#b33845')):
        bx.plot(deltas,[r[key] for r in curve],'.-',ms=3,lw=.6,color=color,label=label)
    band=1.96*.5/np.sqrt(n)
    bx.axhspan(.5-band,.5+band,color='#b1b7bd',alpha=.3,label='Bande ponctuelle 95 % autour de 1/2')
    bx.axhline(.5,color='#34363a',ls='--',lw=1)
    spread=max(.008,4/np.sqrt(n))
    bx.set(xlabel=r'$\Delta$ (rad)',ylabel='Fréquence marginale',ylim=(.5-spread,.5+spread),title='5b  Marginales équilibrées')
    bx.legend(fontsize=8,loc='lower right'); ax.grid(alpha=.2); bx.grid(alpha=.2)
    save(fig,'fig5_probabilites.png')

    fig,ax=plt.subplots(figsize=(8.5,4.8))
    for cv,color in zip((0,.25,.5,.75,1),('#9ba4ad','#7597af','#408573','#cc6b32','#b33845')):
        ax.plot(deltas,theoretical_E(0,-deltas,cv),color=color,label=f'c = {cv:g}, visibilité {(1+cv)/2:g}')
    ax.set(xlabel=r'$\Delta$ (rad)',ylabel='Corrélation théorique',title='6  Cohérence partielle et visibilité')
    ax.grid(alpha=.2); ax.legend(fontsize=9); save(fig,'fig6_coherence.png')

    rows=[]
    for i,d in enumerate(deltas):
        r=curve[i]
        rows.append([d,r['E'],r['E_stderr'],target[i],singlet[i],local[i],population[i],coherence[i],
                     *r['joint_pp_pm_mp_mm'],r['marginal_A_plus'],r['marginal_B_plus']])
    write_csv(out/'courbes.csv',['delta','E_AB','stderr_E','E_c_theory','E_singlet','E_local_branch',
        'E_populations','E_coherence','P_pp','P_pm','P_mp','P_mm','P_A_plus','P_B_plus'],rows)
    write_csv(out/'convergence.csv',['N','RMSE_observed','sqrt_expected_MSE'],zip(ns,errors,expected_rmse))
    write_csv(out/'journal_extrait.csv',['context','trial_in_context','alpha','beta','phi','A','B','AB'],journal)
    result=dict(version=VERSION,coherence=c,pairs_per_setting=n,seed=args.seed,
        scope='Joint sampler using both settings; not an autonomous local detector dynamics.',
        formula='E_c=-(1+c)/2*cos(a-b)',outcome_order=['++','+-','-+','--'],
        chsh=dict(contexts=chsh,S_signed=s_signed,S_abs=s_abs,stderr=s_se,expected_abs=s_exact),
        curve_RMSE_vs_model=rmse(e_sim,target),curve_RMSE_vs_singlet=rmse(e_sim,singlet),
        maximum_marginal_deviation=max(abs(r[k]-.5) for r in curve for k in ('marginal_A_plus','marginal_B_plus')),
        convergence=dict(slope=slope,slope_stderr=slope_se,N=ns.tolist(),rmse=errors.tolist()),
        journal_note='First 1000 trials of each CHSH context; full counts in this JSON.',
        figures=figures,python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'resultats.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(f'Bell 45 v2 | c={c:g} | {n:,} paires/reglage | seed={args.seed}')
    print(f'CHSH |S|={s_abs:.6f} +/- {s_se:.6f} (1 sigma); theorie={s_exact:.6f}')
    print(f'RMSE courbe vs loi c : {result["curve_RMSE_vs_model"]:.6g}')
    print(f'Pente convergence : {slope:.4f} +/- {slope_se:.4f}')
    print(f'Ecart marginal maximal a 0.5 : {result["maximum_marginal_deviation"]:.6g}')
    print(f'Graphiques et donnees : {out}')
    print('Echantillonnage conjoint explicite ; pas une preuve de dynamique locale.')
    if args.no_show: plt.close('all')
    else: plt.show()
    return result

if __name__=='__main__':
    main()
