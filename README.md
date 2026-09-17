# Bell
Bell inequality and pair consistency/synchronisation

[![DOI](https://zenodo.org/badge/DOI/10.5281/1372927596.svg)](https://doi.org/10.5281/zenodo.22795646)


## Version de travail 2.2

Le papier [Cohérence de paire et corrélations de Bell](paper/Bell_Pair_Coherence_v2.2.pdf) et sa [source modifiable](paper/Bell_Pair_Coherence_v2.2.md) présentent une décomposition des statistiques de paire. Le générateur conjoint est une réalisation de la loi cible, pas une dynamique locale démontrée.

Le [guide d'exécution](code/README_Bell_45_v2.md) décrit la simulation, les vérifications et la reconstruction du PDF. L'[audit des réponses locales](code/audit_locality.py) compare explicitement deux polariseurs en parallèle au générateur conjoint, avec préparation avant choix des réglages. Les résultats et journaux sont conservés avec le code.

Les modifications de cette version sont des fichiers de travail ; elles ne mettent pas automatiquement à jour une version déjà archivée sur Zenodo.

## Phase, mémoire et troisième système

La v2.2 ajoute la distinction entre phase fixe et dispersion, un canal quantique local B–C donnant c = cos(theta), une simulation classique de phases couplées par contact, un protocole de mesure hors du plan et une discussion de Laidlaw. La section 9 abrégée et la version 2.1 sont conservées.

- [Mémoire classique de phase](code/phase_memory.py) : sept scénarios, phases complexes, diffusion, écho et contact avec C. Aucun résultat binaire de Bell n'est généré par ce programme.
- [Canal quantique à trois systèmes](code/third_system_coherence.py) : trace partielle, concurrences, information tripartite et maxima CHSH ; le calcul part d'un singulet.
- [Guide des nouveaux calculs](code/README_Phase_Memory.md).

Ces calculs séparent une dynamique classique de synchronisation d'une réalisation quantique de la famille d'états. Leur correspondance physique reste à établir.
