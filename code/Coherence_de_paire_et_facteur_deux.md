# Cohérence de paire et facteur deux : une piste calculable

15 septembre 2026. Cette note développe une piste positive à partir de Bell 48, sans reprendre la démonstration de Bell. Elle sépare un calcul de référence quantique, un exemple de champs classiques et le mécanisme de détection qui reste à construire. Les calculs sont reproduits dans `coherence_facteur_deux.py` ; leurs résultats figurent dans le fichier JSON du même nom.

## Résultat principal

Dans le plan des mesures de Bell 48, le passage de −cos(a−b)/2 à −cos(a−b) admet une décomposition exacte :

\[
E_\phi(a,b)=-\cos(a-\phi)\cos(b-\phi)-\sin(a-\phi)\sin(b-\phi).
\]

Le premier terme est la contribution d'un mélange de deux configurations opposées. Le second est la contribution d'une cohérence entre ces configurations. Pour une orientation φ uniforme, chacune des deux contributions vaut −cos(a−b)/2. Le facteur deux a donc ici un sens précis : la contribution de cohérence égale, après moyenne, celle des populations.

Cette décomposition est obtenue dans le formalisme quantique. Elle fournit une cible pour une recherche de dynamique sous-jacente ; elle ne constitue pas une réalisation classique locale du processus de détection.

Les angles sont des angles de direction de spin/Bloch. Pour la polarisation linéaire des photons, les angles physiques des polariseurs correspondent à des angles de Bloch doubles.

## 1. Raccord au code existant

Dans `simu/Bell_test_48_heatmap.py`, les probabilités cos²((a−φ)/2) donnent les moyennes conditionnelles cos(a−φ) et −cos(b−φ) pour une branche opposée. Les tirages conditionnellement indépendants donnent donc :

\[
E_{\rm branche}(a,b|\phi)=-\cos(a-\phi)\cos(b-\phi).
\]

Sa moyenne uniforme est −cos(a−b)/2. Le mélange à poids égaux des configurations (+φ,−φ) et (−φ,+φ) possède exactement cette même corrélation, avec marginales 1/2 même à φ fixé. Le raccord au mélange symétrique porte sur la paire et sur les statistiques observées ; une branche orientée seule n'a pas ces marginales conditionnelles à φ fixé.

Le code calcule également la somme de deux scores de branche. Cette opération ne crée aucun terme d'interférence entre amplitudes. La piste étudiée ici consiste à rechercher une autre structure dans la paire, susceptible de modifier les probabilités de détection elles-mêmes.

## 2. Deux populations et une cohérence

Notons |+_φ −_φ⟩ et |−_φ +_φ⟩ les deux configurations de spin opposées sur un axe φ. Il s'agit de configurations physiques dans deux sous-systèmes, pas simplement d'un échange des noms de particules.

Le mélange incohérent est

\[
\rho_{\phi,0}=\tfrac12|+_\phi-_\phi\rangle\langle+_\phi-_\phi|
+\tfrac12|-_\phi+_\phi\rangle\langle-_\phi+_\phi|.
\]

Introduisons une cohérence réelle c, entre 0 et 1 :

\[
\rho_{\phi,c}=\rho_{\phi,0}-\frac c2\left(
|+_\phi-_\phi\rangle\langle-_\phi+_\phi|
+|-_\phi+_\phi\rangle\langle+_\phi-_\phi|\right).
\]

Ses valeurs propres sont (1+c)/2, (1−c)/2, 0, 0. Elle est positive et normalisée. Les deux états réduits sont I/2 pour tout c et tout φ. La cohérence est donc invisible dans toutes les mesures d'un seul membre de la paire.

Pour c=1, l'état est le singulet. Pour c=0, c'est le mélange incohérent. La comparaison entre ces deux états, incluant leur reconstruction tomographique, est présentée expérimentalement par [Virzí et collègues, 2019](https://arxiv.org/abs/1810.09331). Le calcul qui suit applique cette distinction à la moyenne planaire utilisée dans les simulations locales du dossier.

Avec l'observable σ(a)=cos(a)σ_z+sin(a)σ_x :

\[
\operatorname{Tr}[\rho_{\phi,c}\,\sigma(a)\otimes\sigma(b)]
=-\cos(a-\phi)\cos(b-\phi)-c\sin(a-\phi)\sin(b-\phi).
\]

Pour φ=0, le mélange est parfaitement anticorrélé sur z mais ne présente pas de corrélation sur x. Le singulet est parfaitement anticorrélé sur les deux axes. Cela correspond exactement à la différence entre mélange et singulet évoquée dans la conversation.

Après moyenne uniforme sur φ :

\[
\overline E_c(a,b)=-\frac{1+c}{2}\cos(a-b).
\]

Le facteur deux n'est pas une propriété universelle des corrélations classiques : il vient ici de la moyenne sur un cercle et de deux contributions de même poids. Par exemple une distribution uniforme de directions de spin sur la sphère donnerait une contribution de mélange −a·b/3, et non −a·b/2.

## 3. Ce qui change dans les quatre cases

Pour s,t∈{−1,+1}, les probabilités après moyenne planaire sont

\[
P_c(s,t|a,b)=\frac14\left[1-st\frac{1+c}{2}\cos(a-b)\right].
\]

À réglages alignés :

| c | P++ | P+− | P−+ | P−− | E |
|---|---|---|---|---|---|
| 0 | 1/8 | 3/8 | 3/8 | 1/8 | −1/2 |
| 1/2 | 1/16 | 7/16 | 7/16 | 1/16 | −3/4 |
| 1 | 0 | 1/2 | 1/2 | 0 | −1 |

Le terme ajouté à φ fixé vaut

\[
\Delta P_{st}=-\frac{st\,c}{4}\sin(a-\phi)\sin(b-\phi).
\]

Sa somme sur chaque ligne et chaque colonne est nulle. Ce terme est signé, contrairement à une probabilité d'événement supplémentaire. Il redistribue les poids entre les quatre cases. Les probabilités totales restent positives.

Avec deux résultats binaires et des marginales équilibrées, toute distribution conjointe est déjà entièrement fixée par un seul nombre E : P_st=(1+stE)/4. Une éventuelle complexité supplémentaire se trouve donc dans le processus qui produit cette table selon les réglages, dans le champ, dans sa mémoire ou dans les séries temporelles. Elle ne se cache pas dans un autre paramètre indépendant des quatre cases finales.

## 4. D'où vient le terme sinus × sinus ?

Pour une issue donnée, appelons u et v les amplitudes provenant des deux configurations de paire. L'intensité/probabilité du mélange est (|u|²+|v|²)/2. Une superposition cohérente antisymétrique donne

\[
\tfrac12|u-v|^2=\tfrac12(|u|^2+|v|^2)-\operatorname{Re}(uv^*).
\]

Le terme d'interférence est précisément celui qui manque au mélange. Pour les projecteurs de spin dans le plan, sa contribution à E est −sin(a−φ)sin(b−φ).

On peut produire ρ_{φ,c} en moyennant les états

\[
|\psi_{\phi,\chi}\rangle=\frac{|+_\phi-_\phi\rangle-e^{i\chi}|-_\phi+_\phi\rangle}{\sqrt2}
\]

avec ⟨e^{iχ}⟩=c réel. Une phase relative fixe χ=0 donne c=1 ; une phase entièrement brouillée donne c=0. Le code vérifie une réalisation intermédiaire avec les phases ±arccos(c).

On peut aussi associer à chaque configuration un état d'environnement |e_1⟩ ou |e_2⟩. Après élimination de cet environnement, le terme de cohérence est proportionnel à leur recouvrement. Cela donne un sens physique à l'indiscernabilité : ce qui importe est la persistance de la cohérence, et pas uniquement notre ignorance de la branche. Deux branches non identifiées par l'observateur peuvent être incohérentes.

À c=1, l'axe φ est une manière de décomposer le singulet invariant par rotation commune ; ce n'est pas la preuve qu'une paire possède secrètement deux spins classiques fixés sur cet axe.

## 5. Un exemple effectivement classique de corrélation de champ

Pour éviter de s'arrêter à une réécriture quantique, considérons un signal classique aléatoire à deux composantes complexes z, gaussien circulaire, avec ⟨zz†⟩=I. Une source fournit E_A=z et E_B=Jz*, où J est la matrice antisymétrique [[0,1],[−1,0]]. C'est une prescription de source corrélée ; on n'affirme pas qu'un simple montage optique passif la réalise.

Chaque côté projette son propre champ sur les deux canaux de son analyseur et calcule une intensité I_s^A ou I_t^B. Pour les angles de Bloch a,b dans le plan, le théorème des moments gaussiens donne :

\[
G_{st}=\langle I_s^A I_t^B\rangle
=1+\frac{1-st\cos(a-b)}2.
\]

Les intensités moyennes valent chacune 1. La partie connectée est donc

\[
C_{st}=G_{st}-\langle I_s^A\rangle\langle I_t^B\rangle
=\frac{1-st\cos(a-b)}2.
\]

Normalisée par sa somme, cette partie connectée donne exactement la table du singulet. Les intensités conjointes complètes, normalisées par leur propre somme, donnent en revanche E=−cos(a−b)/3. Le fond constant explique entièrement l'écart entre ces deux observables.

Les structures de quatrième ordre reliant champs et intensités sont le domaine de la relation de Siegert ; voir [Ferreira et collègues, 2020](https://arxiv.org/abs/2002.05425). La construction vectorielle ci-dessus est calculée explicitement ici. Des expériences de cohérence classique utilisant deux degrés de liberté d'un même faisceau sont également décrites par [Kagalwala et collègues, 2013](https://pubblicazioni.unicam.it/handle/11581/256381).

Cette expérience numérique confirme qu'une loi angulaire de type singulet peut se trouver dans une structure de cohérence classique. Elle ne transforme pas cette covariance d'intensité en un journal de résultats binaires. Il faudrait qu'un mécanisme de détection physique sélectionne cette contribution, avec un bilan complet des événements, pour aller plus loin.

Pour 600 000 réalisations et a−b=−0,82, le code donne : E complet=−0,22754, contre −0,22741 attendu ; E connecté=−0,68130, contre −0,68222 attendu. Ces nombres sont un contrôle Monte-Carlo des moments, pas un test de particules.

## 6. Cible pour une dynamique onde-particule

La question devient : existe-t-il, dans le champ et ses deux configurations de paire, un analogue physique de la cohérence croisée qui survive à la séparation et contribue aux événements détectés ?

Une recherche concrète peut suivre quatre étapes :

1. Définir deux configurations physiques de paire et les composantes du champ qui les accompagnent. Une permutation des étiquettes seule ne suffit pas à définir deux amplitudes physiques.
2. Mesurer leur cohérence croisée, éventuellement de quatrième ordre, avant de convertir les variables continues en deux signes. Une synchronisation des positions seule pourrait ne pas épuiser cette information.
3. Brouiller sélectivement la phase relative ou une variable de mémoire, en contrôlant que les distributions de chaque système restent identiques. Vérifier si la contribution conjointe diminue.
4. Construire explicitement la réponse des détecteurs au champ et calculer les quatre fréquences sans remplacer les événements par une covariance ou une intensité reconstruite.

Dans la famille de référence étudiée, une prédiction supplémentaire est disponible. Le tenseur de corrélation après moyenne planaire vaut diag(−(1+c)/2, −c, −(1+c)/2) dans les axes x,y,z. Si la visibilité dans le plan est V, la corrélation sur l'axe y vaut donc 1−2V. Une famille de Werner, par exemple, aurait −V sur les trois axes. Une même courbe planaire à demi-amplitude ne détermine pas toute la structure de paire ; mesurer la composante manquante permet de distinguer les modèles.

On ne dispose pas ici d'une dynamique classique locale complète produisant les clics du singulet. On dispose d'un terme manquant identifiable, d'un raccord exact au facteur deux de Bell 48, d'une famille positive conservant les marginales, d'un exemple classique de cohérence conjointe et de perturbations discriminantes à étudier.

## Vérification

Le script vérifie par produits de matrices, indépendamment des expressions trigonométriques : positivité, normalisation, marginales, probabilités conjointes et corrélations pour cinq valeurs de c, plusieurs axes de préparation et une grille de réglages. L'écart maximal des identités est 4,44×10⁻¹⁶. Il vérifie séparément la construction par phases aléatoires et la redistribution signée des quatre probabilités. Le modèle classique gaussien fait l'objet d'un contrôle Monte-Carlo indépendant.
