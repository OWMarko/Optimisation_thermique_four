#include "../include/solver_direct.hpp"
#include <vector>
#include <cmath>
#include <iostream>

// --- Assemblage de la matrice de rigidité (Exercice 11.1 - Q2 et Q3) ---
Eigen::SparseMatrix<double> assa(int Nbpt, int Nbtri, 
                                 const Eigen::MatrixXd& Coorneu, 
                                 const Eigen::VectorXi& Refneu, 
                                 const Eigen::MatrixXi& Numtri, 
                                 const Eigen::VectorXi& Reftri, 
                                 const Eigen::VectorXd& Conduc) {
    
    // En C++, on utilise une liste de triplets pour construire une matrice creuse efficacement
    std::vector<Eigen::Triplet<double>> triplets; // contient que trois valeurs, ligne Numtri(k,i), colonne Numtri(k,j), valeur calculée
    // on ne stocke pas les zéros, on stocke que les valeurs non nulles, c'est pour ça qu'on utilise une liste de triplets, c'est plus efficace que de construire une matrice dense et ensuite la convertir en sparse, car la conversion peut être coûteuse en temps et en mémoire si la matrice est grande et très creuse (O(N^2) en mémoire et O(N^3) en temps pour la conversion), alors que la liste de triplets permet de construire directement la matrice sparse sans passer par une matrice dense intermédiaire, c'est plus rapide et plus économe en mémoire, surtout pour les grands maillages où la matrice de rigidité est très creuse (beaucoup de zéros)
    
    // Optimisation R&D : on réserve la mémoire (9 termes par triangle)
    // 3 noeuds i, 3 noeuds j, matrice de taille 9, et un triangle a 3 sommets, donc 9 contributions non nulles par triangle, on réserve donc Nbtri * 9 pour éviter les reallocations dynamiques de la liste de triplets qui peuvent être coûteuses en temps
    triplets.reserve(Nbtri * 9); 
    
    // Boucle sur les triangles
    for (int k = 0; k < Nbtri; ++k) {
        
        // Numéros des 3 sommets du triangle k
        // On veut récuperer les numéros des sommets dans notre Numtri pour pouvoir ensuite récupérer les coordonnées de ces sommets dans Coorneu, et calculer les gradients et les intégrales pour assembler la matrice de rigidité, c'est la base de la méthode des éléments finis, on travaille avec les indices des sommets pour faire le lien entre la connectivité du maillage (Numtri) et les coordonnées géométriques (Coorneu)
        int n1 = Numtri(k, 0); 
        int n2 = Numtri(k, 1);
        int n3 = Numtri(k, 2);
        
        // Coordonnées X et Y des 3 sommets
        // C'est ce qu'on vient de faire, on récupère les coordonnées x et y de chaque sommet du triangle k en utilisant les indices n1, n2, n3 pour accéder à Coorneu, c'est essentiel pour calculer le déterminant Delta et la matrice du gradient Dp qui sont nécessaires pour assembler la matrice de rigidité élémentaire, et ensuite assembler la matrice globale Atotal
        double x1 = Coorneu(n1, 0), y1 = Coorneu(n1, 1);
        double x2 = Coorneu(n2, 0), y2 = Coorneu(n2, 1);
        double x3 = Coorneu(n3, 0), y3 = Coorneu(n3, 1);
        
        // Le déterminant Delta (qui vaut 2 fois l'aire du triangle)
        double Delta = std::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
        
        // La matrice du gradient Dp
        // On divise directement par Delta ici pour simplifier les calculs suivants
        // Déclaration propre et lisible
        Eigen::Matrix<double, 2, 3> Dp; // on fixe bien les dimensions de la matrice Dp, 2 lignes pour les composantes x et y du gradient, et 3 colonnes pour les 3 sommets du triangle, c'est plus clair que d'utiliser une matrice dynamique, on sait exactement la taille qu'on veut, et ça permet au compilateur d'optimiser les calculs, notamment en utilisant des instructions SIMD pour faire les opérations en parallèle sur les éléments de la matrice, c'est un avantage du C++ par rapport à d'autres langages plus dynamiques comme Python

        // Remplissage par blocs ou par colonnes
        Dp << (y2 - y3), (y3 - y1), (y1 - y2),
              (x3 - x2), (x1 - x3), (x2 - x1);

        Dp /= Delta; // Division scalaire sur toute la matrice d'un coup
        
        // Conductivité du matériau pour ce triangle
        double c = Conduc(Reftri(k));

        // Calcul de la matrice de rigidité locale 3x3
        // C'est l'équivalent vectorisé exact de la formule : c * Delta * (Dp(:,i)' * Dp(:,j)) / 2 
        Eigen::Matrix3d K_loc = (c * Delta / 2.0) * (Dp.transpose() * Dp); // comme c'est une matrice sym, Eigen comprendra et optimisera les calculs, on n'a pas besoin de faire les deux boucles pour i et j, on peut faire le produit matriciel Dp.transpose() * Dp qui va nous donner une matrice 3x3, et ensuite on multiplie par c * Delta / 2.0 pour obtenir la matrice de rigidité locale K_loc, c'est plus rapide que de faire les boucles imbriquées pour calculer chaque élément de K_loc un par un, et ça rend le code plus lisible et plus facile à maintenir, surtout si on veut faire des modifications ou des extensions dans le futur

        // c) Assemblage de la matrice globale via les triplets 
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) { 
                triplets.push_back(Eigen::Triplet<double>(Numtri(k, i), Numtri(k, j), K_loc(i, j))); // comme on a déja réserver en mémoire, rajouter à la fin du vecteur le triplet coute que O(1)
            }
        }
    }
    
    // Construction finale de la matrice creuse en un seul passage (O(N log N))
    Eigen::SparseMatrix<double> Atotal(Nbpt, Nbpt); 
    Atotal.setFromTriplets(triplets.begin(), triplets.end()); // au lieu de faire des insertions individuelles dans la matrice creuse qui peuvent être coûteuses en temps (O(N^2) dans le pire des cas), on utilise la méthode setFromTriplets qui construit la matrice de manière efficace à partir de la liste de triplets, c'est optimisé pour les grands maillages où la matrice de rigidité est très creuse, et ça rend le code plus propre et plus facile à comprendre, on voit clairement que la matrice est construite à partir des contributions locales de chaque triangle, et on évite les complexités liées à l'insertion individuelle dans une matrice creuse, de plus comme un sommet appartient souvent à plusieurs triangles, cette fonction fait le trie et additionne les doublons
    // c'est une fonctions en 3 étapes, la première est le tri en O(N log N), d'abord colonne puis lignes, ensuite il parcours de haut en bas pour voir les doublons (i,j), enfin il compresse en regardant le nombre exact de valeurs et alloue le strict nécessaire en mémoire pour la matrice creuse (Compressed Sparse Column - CSC) et mets les différentes valeurs
    return Atotal;
}

// --- Calcul du second membre pour une résistance (Exercice 11.1 - Q3) ---
// On utilise la formulation Delta/6.0 conforme au document du professeur
Eigen::VectorXd calculer_b_re(int Nbpt, int triangle_index, const Maillage& m) {
    // Initialisation d'un vecteur b rempli de zéros
    Eigen::VectorXd b = Eigen::VectorXd::Zero(Nbpt); // tous les triangles ont F = 0 sauf celui qui contient le résistor qui a F = 1
    
    // On extrait les identifiants puis les coordonnées cartésiennes des 3 sommets de l'unique triangle qui sert de résistance.
    // Indices des 3 sommets du triangle où se trouve la résistance
    int n1 = m.Numtri(triangle_index, 0);
    int n2 = m.Numtri(triangle_index, 1);
    int n3 = m.Numtri(triangle_index, 2);
    
    // Coordonnées pour calculer le Déterminant (Delta)
    double x1 = m.Coorneu(n1, 0), y1 = m.Coorneu(n1, 1);
    double x2 = m.Coorneu(n2, 0), y2 = m.Coorneu(n2, 1);
    double x3 = m.Coorneu(n3, 0), y3 = m.Coorneu(n3, 1);
    
    // Delta = |(x2-x1)(y3-y1) - (x3-x1)(y2-y1)|, Jacobienne
    double Delta = std::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
    
    // Contribution élémentaire selon le PDF : (F * Delta) / 6
    // On considère F = 1 pour calculer la signature thermique unitaire
    double contribution = Delta / 6.0;
    
    // On ajoute la contribution aux trois sommets du triangle
    // On injecte cette chaleur locale dans le grand vecteur global b aux adresses correspondantes aux sommets du triangle
    b(n1) += contribution;
    b(n2) += contribution;
    b(n3) += contribution;
    
    return b;
}

// --- Prise en compte des conditions de Dirichlet (Exercice 11.1 - Q4) ---
void elim(Eigen::SparseMatrix<double>& Atotal, Eigen::VectorXd& belim, 
          const Eigen::VectorXi& Refneu, int Refdir, double Valdir) {
    
    int n = Atotal.rows();
    
    // Étape 1 : Modification du second membre (b = b - A_2 * T_D)
    // On transfère l'influence de la température connue vers les nœuds inconnus
    // cette double boucle parcours les valeurs non nuls de la matrice creuse, dans un premier temps on parcours les colonnes k
    // ensuite nous parcourons les lignes it de ces colonnes
    for (int k = 0; k < Atotal.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Atotal, k); it; ++it) {
            int row = it.row(); // on récupère les indices 
            int col = it.col();
            
            // Si la colonne appartient au bord (elle représente un nœud à la température connue $T_D$) ET que la ligne n'est pas sur le bord (c'est un nœud interne inconnu), alors on est pile dans le bloc matriciel $A_2$ de la formule.
            if (Refneu(col) == Refdir && Refneu(row) != Refdir) {
                belim(row) -= it.value() * Valdir; // Application de -A_2 * T_D, c'est b - A_2 * T_D, on soustrait la contribution de la température connue du second membre pour les nœuds internes, c'est essentiel pour que le système linéaire modifié reflète correctement les conditions de Dirichlet imposées, et que la solution finale respecte ces conditions aux limites, c'est une étape clé dans la méthode des éléments finis pour traiter les conditions de Dirichlet de manière efficace et précise
            }
        }
    }
    
    // Étape 2 : Nettoyage de la matrice (Création des blocs 0 et I)
    for (int k = 0; k < Atotal.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Atotal, k); it; ++it) {
            int row = it.row();
            int col = it.col();
            
            // Si le terme touche une ligne ou une colonne du bord imposé
            //  si on touche un terme qui est sur une ligne du bord OU une colonne du bord (les blocs $A_2$, $A_2^T$ et $A_3$ ), on doit le modifier.
            if (Refneu(row) == Refdir || Refneu(col) == Refdir) {
                if (row == col) {
                    it.valueRef() = 1.0; // Bloc I : On met 1 sur la diagonale | Si le numéro de ligne est égal au numéro de colonne (c'est la diagonale du bloc $A_3$ ), on écrase la valeur thermique par 1.0. Cela crée la sous-matrice identité $I$.
                } else {
                    it.valueRef() = 0.0; // Blocs 0 : On efface les termes croisés |  : Si ce n'est pas la diagonale, on écrase la valeur par 0.0. On vient virtuellement d'effacer les blocs $A_2$ et $A_2^T$  pour les remplacer par des zéros, garantissant que le nœud du bord ne dépend plus de rien d'autre que de lui-même.
                }
            }
        }
    }
    
    // Étape 3 : Imposition de la température T_D sur le bloc droit
    // Enfin, on parcourt la liste de tous les nœuds de 0 à n-1.Si le nœud i appartient à la frontière à imposer, on remplace complètement sa ligne dans le second membre par la température désirée (Valdir).C'est le bloc du bas de la formule : $T_D$.
    for (int i = 0; i < n; ++i) {
        if (Refneu(i) == Refdir) {
            belim(i) = Valdir; // Application de T_D
        }
    }
}

// --- Calcul des signatures thermiques Tk pour chaque résistance ---
Eigen::MatrixXd computeResistanceSignatures(const Maillage& mesh,
                                            const Eigen::SparseMatrix<double>& stiffness,
                                            const std::vector<int>& resistor_indices) {
    int nr = resistor_indices.size();
    Eigen::MatrixXd Tres(mesh.Nbpt, nr);

    // 1. Préparation de la matrice (UNE SEULE FOIS)
    Eigen::SparseMatrix<double> A_signatures = stiffness;
    Eigen::VectorXd dummy_b = Eigen::VectorXd::Zero(mesh.Nbpt);
    
    // On modifie la matrice pour imposer 0 sur les bords (car Tk = T - T0, et T=T0 sur les bords)
    elim(A_signatures, dummy_b, mesh.Refneu, 1, 0.0);
    elim(A_signatures, dummy_b, mesh.Refneu, 2, 0.0);

    // 2. Factorisation (LE CALCUL LOURD, UNE SEULE FOIS)
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A_signatures);

    if(solver.info() != Eigen::Success) {
        std::cerr << "Erreur: La factorisation de la matrice des signatures a échoué." << std::endl;
    }

    // 3. Boucle ultra-rapide sur les résistances
    for (int k = 0; k < nr; ++k) {
        // Construction du second membre brut pour la résistance k
        Eigen::VectorXd b_k = calculer_b_re(mesh.Nbpt, resistor_indices[k], mesh);
        
        // Application manuelle de Dirichlet T=0 sur b_k pour rester cohérent avec la factorisation
        for (int i = 0; i < mesh.Nbpt; ++i) {
            if (mesh.Refneu(i) == 1 || mesh.Refneu(i) == 2) {
                b_k(i) = 0.0;
            }
        }

        // Résolution instantanée par substitution
        Tres.col(k) = solver.solve(b_k);
    }

    return Tres;
}

// --- Calcule la température de base T0 (sans résistances) ---
Eigen::VectorXd computeBaseTemperature(const Maillage& mesh, 
                                       const Eigen::SparseMatrix<double>& stiffness,
                                       double T_top, double T_bottom) {
    
    Eigen::SparseMatrix<double> A = stiffness;
    Eigen::VectorXd b = Eigen::VectorXd::Zero(mesh.Nbpt);

    // Conditions aux limites de Dirichlet
    elim(A, b, mesh.Refneu, 1, T_top); 
    elim(A, b, mesh.Refneu, 2, T_bottom); 

    // Résolution LDLT
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A); 
    
    if(solver.info() != Eigen::Success) {
        std::cerr << "Erreur: La factorisation de la matrice T0 a échoué." << std::endl;
    }

    return solver.solve(b);
}