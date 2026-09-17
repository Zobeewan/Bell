# Mémoire de phase et couplage à un troisième système

Depuis la racine du dépôt, avec Python, NumPy et Matplotlib :

```sh
python code/phase_memory.py
python code/third_system_coherence.py
```

Les programmes enregistrent leurs figures, CSV et rapports JSON dans leurs dossiers `resultats_*`. Ils incluent des contrôles analytiques ; une assertion échoue si un contrôle n'est pas satisfait. La graine par défaut de la simulation classique est 20260917. `python code/phase_memory.py --help` décrit ses options.

## Ce qui est simulé

`phase_memory.py` prépare A et B par couplage sinusoïdal, les sépare, puis fait évoluer les phases. Il compare évolution libre, dérive fixe, diffusion, dispersion des fréquences, écho, contact BC seul et contact BC avec AB maintenu. Q_AB est la moyenne complexe de exp(i(theta_B-theta_A)) ; R_AB en est le module. Les temps et paramètres sont en unités arbitraires. Les termes de contact sont des oscillateurs réduits, sans géométrie de propagation ni détecteur. Une phase de champ et une cohérence de configurations quantiques ne sont pas identifiées automatiquement.

`third_system_coherence.py` applique au singulet AB une rotation contrôlée de C par B. La trace de C produit c = cos(theta). Ce calcul quantique prédit le changement de la paire sans changement des marginales ; il ne dérive pas le singulet d'un système classique. Il vérifie également le cas d'une phase complexe, les concurrences des trois paires, le résidu tripartite et le retour par la transformation inverse.

## Fichiers de résultats

- `resultats_phase_memory/trajectories.csv` : parties réelles, imaginaires et modules des cohérences AB, BC, AC.
- `states.npz` : phases initialement préparées et phases finales, avec noms de scénarios.
- `interactive_data.json` : 61 instants ; les vecteurs moyens utilisent les 4 096 réalisations, les nuages affichables en retiennent 48.
- `report.json` : paramètres, résultats, contrôles et empreinte du script.
- `resultats_third_system/coupling.csv` : cohérence, concurrences, tangle, distinguabilité et CHSH pour 61 couplages.

## Vérifications

Les scénarios libre et dérive sont comparés aux solutions complexes exactes. Diffusion et dispersion sont comparées à leurs moyennes analytiques avec une tolérance d'ensemble de 6/sqrt(N) ; cette tolérance est un contrôle numérique, pas un intervalle de confiance ajusté aux courbes. L'écho et l'absence de changement de A lors du contact BC sont vérifiés séparément. Le contrôle à demi-pas concerne le scénario déterministe de contact BC seul, pas la convergence forte du bruit stochastique.

Le calcul quantique vérifie les 61 traces partielles et leurs marginales, la concurrence et la relation tripartite, les CHSH via les valeurs singulières de T, ainsi que l'inversion de l'interaction. Les phases fixes sont testées avant la moyenne d'orientation.

Les nouvelles figures et leur interprétation sont dans les sections 11–15 du [papier v2.2](../paper/Bell_Pair_Coherence_v2.2.pdf). Pour reconstruire ce PDF, installer aussi ReportLab, Pillow et pypdf, puis exécuter `python code/build_bell_v2_paper.py`.
