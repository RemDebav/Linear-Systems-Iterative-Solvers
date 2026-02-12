# Résolution de Systèmes Linéaires : Méthodes Itératives

Ce projet implémente et compare les principales méthodes itératives (Jacobi, Gauss-Seidel et SOR) sous MATLAB pour la résolution de systèmes de type $Ax = b$ de grande dimension.

## 🚀 Méthodes Implémentées
* **Jacobi** : Algorithme de base pour les matrices à diagonale dominante.
* **Gauss-Seidel** : Implémentation implicite pour une convergence accélérée.
* **SOR (Successive Over-Relaxation)** : Optimisation de la convergence via un facteur de relaxation $\omega$.

## 📊 Analyse de Performance
Le projet inclut des outils de mesure et de visualisation :
Compare le temps d'exécution des trois méthodes en fonction de la taille de la matrice.
Recherche expérimentale du meilleur facteur $\omega$ pour minimiser le nombre d'itérations du solveur SOR.
Utilitaires pour créer des matrices définies positives et symétriques afin de garantir la convergence.

## 📁 Structure du Projet
* `/src` : L'ensemble des scripts et des fonctions Matlab.
  - aux* : Fonction auxiliaire
  - method* : Fonction implémentant la méthode x.
  - analyse* : Script permettant d'analyser différentes méthodes et l'influence de certaines variables
* `/scripts` : Scripts de test et de comparaison.
* `/results` : Graphiques de performance et benchmarks.

## 🛠️ Utilisation
1. Clonez le dépôt.
2. Ajoutez le dossier `src` au chemin MATLAB.
3. Lancez `scripts/comparaison.m` pour visualiser les performances.

---
*Projet réalisé dans le cadre de mon cursus en Mathématiques Appliquées à Polytech Nice Sophia.*
