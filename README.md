# Projet Thermique - Simulation et Optimisation

Ce projet vise à résoudre l'équation de la chaleur en deux dimensions par la méthode des Éléments Finis. Nous résolvons :
 
1. Problème direct : Calculer le champ de température dans un four contenant différentes zones matérielles et des résistances chauffantes.
2. Problème inverse : Optimiser automatiquement la puissance de ces résistances pour atteindre une température précise sur la pièce à cuire.

---

## Architecture du Code

Le projet est écrit en C++ (calcul) et s'accompagne de scripts Python (visualisation). 

### Dossier `src/` et `include/` (Cœur de calcul)
* `maillage.cpp / .hpp` : Génération et gestion de la grille de calcul (nœuds, triangles). Permet l'identification spatiale des différentes zones comme l'air et la résine.
* `solver_direct.cpp / .hpp` : Assemblage des matrices élémentaires, application des conditions aux limites de Dirichlet et résolution du système linéaire pour obtenir le champ de température.
* `optimisation.cpp / .hpp` : Formulation et résolution du problème inverse via l'assemblage des matrices de masse locales. Contient également la logique de pénalisation pour lisser les puissances calculées.
* `export.cpp / .hpp` : Enregistrement des résultats et de la géométrie dans des fichiers texte lisibles par d'autres logiciels.
* `experiences.cpp / .hpp` : Comparer facilement différentes configurations.
* `main.cpp` : Point d'entrée du programme.

### Autres dossiers
* `scripts/` : Contient les utilitaires Python pour lire les données exportées et tracer. 
* `data/` : Dossier de destination des fichiers générés par l'exécutable C++.
* `img/` : Dossier de destination des graphiques générés par les scripts Python.

---

## Prérequis et Compilation

Le projet utilise CMake pour automatiser la compilation et s'appuie sur la bibliothèque Eigen pour l'algèbre linéaire. 

### Étapes de compilation :
Ouvrez un terminal à la racine du projet et tapez les commandes suivantes :

```bash
mkdir build
cd build
cmake ..
make
```

### Exécution :
Toujours depuis le dossier `build`, lancez l'exécutable :
```bash
./Four
```

Une fois le programme terminé, vous pouvez lancer les scripts d'affichage depuis le dossier `scripts/` :
```bash
python3 affichage.py visualisation_maillage.py
```

---

## Paramètres configurables

La plupart des paramètres physiques et géométriques sont modifiables dans les appels de fonctions dans le main. N'hésitez pas à les ajuster pour observer comment le modèle réagit.

### Paramètres géométriques et matériels
* `nx = 40`, `ny = 40` : Résolution du maillage. Augmenter ces valeurs affine les calculs mais augmente le temps de résolution et les risques d'instabilités du problème inverse.
* `k_air = 1.0` : Conductivité thermique de l'air.
* `k_resin = 10.0` : Conductivité thermique de l'objet à cuire.

### Paramètres thermiques (Problème direct)
* `T_haut = 50.0` : Température imposée sur la paroi supérieure (°C).
* `T_bas = 100.0` : Température imposée sur la paroi inférieure (°C).
* `alpha_manuel = 200000.0` : Puissance arbitraire (en Watts) attribuée aux résistances lors du test en mode direct.

### Positionnement des résistances
Vous pouvez déplacer ou ajouter des résistances en modifiant le vecteur de coordonnées.
* `pos_res = {{-0.7, 0.7}, {0.0, 0.7}, ...}` : Positions (x, y) des centres des résistances dans le four.

### Paramètres d'optimisation (Problème inverse)
* `T_cible = 250.0` : La température idéale que l'on souhaite atteindre dans la zone de la résine.
* `C_test = 1e-9` : Poids de la pénalité. Une valeur proche de 0 laisse le modèle libre au risque de valeurs extrêmes, une valeur plus grande lisse l'effort énergétique entre les résistances au détriment de la précision de cuisson.

---

## Documentation

Vous retrouverez à la racine du projet :
* `Recherche.pdf` : Document détaillé présentant les travaux de recherche, le détail des calculs et la formulation mathématique du problème.
* `Slides.pdf` : Présentation synthétique du projet pour l'oral.


