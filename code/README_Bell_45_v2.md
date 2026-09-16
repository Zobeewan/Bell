# Bell 45 v2

La version modifiée est `Bell_test_45.py`. L'original est conservé dans `Bell_test_45_original.py`.

La v2 génère un seul couple de signes par essai selon la loi conjointe

\[
P_c(A,B\mid a,b)=\frac14\left[1-AB\frac{1+c}{2}\cos(a-b)\right].
\]

Elle utilise la décomposition conditionnelle à l'orientation φ en populations et cohérence. Le générateur conjoint connaît les deux réglages. Ce programme n'est pas présenté comme une dynamique de deux détecteurs locaux autonomes. Les figures comparent la loi et les fréquences réellement tirées ; aucun doublage du produit AB n'est appliqué.

## Lancement

Depuis la racine du projet, avec votre Python habituel :

```powershell
python simu/Bell_test_45.py
```

Les figures s'affichent ensemble à la fin et sont aussi enregistrées. Pour enregistrer sans ouvrir de fenêtres :

```powershell
python simu/Bell_test_45.py --coherence 1 --no-show
python simu/Bell_test_45.py --coherence 0.5 --no-show --output-dir simu/resultats_c05
```

Paramètres : `--coherence` ou `--c` (défaut 1), `--pairs` (200000 par réglage), `--seed` (20260916), `--output-dir`, `--no-show`, `--quick` (grilles de calcul réduites). Le nombre de paires n'est pas réduit par `--quick` : utiliser aussi `--pairs 20000` pour un essai rapide.

Les dépendances sont NumPy et Matplotlib. SciPy n'est plus nécessaire. Si elles manquent à votre Python :

```powershell
python -m pip install numpy matplotlib
```

Un lanceur PowerShell utilise également le runtime Codex disponible sur cette machine et les dépendances déjà préparées dans le projet :

```powershell
powershell -ExecutionPolicy Bypass -File simu/Lancer_Bell_45_v2.ps1 -SansFenetre
```

Ce réglage de politique ne concerne que le processus lancé ; aucune politique Windows permanente n'est modifiée.

## Résultats

Le dossier par défaut `simu/resultats_bell45_v2` contient :

- `fig1_E_vs_delta.png` : courbe des produits AB et comparaison au singulet et à l'ancienne branche indépendante ;
- `fig2_branches.png` : contributions populations et cohérence à l'espérance ;
- `fig3_CHSH.png` : CHSH calculé à partir des fréquences conjointes ;
- `fig4_convergence.png` : convergence de l'erreur vers la loi au c choisi ;
- `fig5_probabilites.png` : quatre fréquences conjointes et deux marginales ;
- `fig6_coherence.png` : famille théorique pour plusieurs valeurs de c ;
- `resultats.json` : effectifs complets CHSH, incertitudes, paramètres, versions et empreinte du code ;
- `courbes.csv`, `convergence.csv` : données des tracés ;
- `journal_extrait.csv` : les 1000 premiers essais de chaque contexte CHSH. Cet extrait n'est pas le journal complet utilisé pour les résultats.

À c=1, exécution de référence : |S| = 2,82809 ± 0,00316 (un écart-type) ; cible 2,82843 ; RMSE angulaire 0,001365. Les marginales valent 1/2 en probabilité et fluctuent dans un échantillon fini.

La cohérence partielle contrôle une famille d'états, pas la fraction temporelle d'une interaction. Dans cette famille, c=0 est séparable, c>0 est intriqué et c=1 donne le singulet.

## Vérification et papier

```powershell
python analyse/test_bell45_v2.py
```

Ces contrôles comprennent une comparaison indépendante aux matrices de densité et projecteurs, les probabilités positives normalisées, les fréquences des événements, les marginales et la reproductibilité. Ils ne se limitent pas à vérifier une égalité avec la fonction qui produit les tirages.

Le papier est dans `output/pdf/Bell_Pair_Coherence_v2.pdf`. La source modifiable est le fichier Markdown de même nom, avec formules LaTeX et liens vers les références. Le générateur est `analyse/build_bell_v2_paper.py` (ReportLab, Matplotlib, Pillow et pypdf). Les chemins des figures sont relatifs au Markdown : conserver l'arborescence pour reconstruire le PDF.
