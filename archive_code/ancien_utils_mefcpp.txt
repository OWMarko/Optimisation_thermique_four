#include "../include/utils_mef.hpp"
#include <fstream>
#include <cmath>
#include <iostream>

Maillage creer_maillage(int Nx, int Ny) {
    Maillage m;
    
    // Dimensions du four : carré de -1 à 1
    double x_min = -1.0, x_max = 1.0;
    double y_min = -1.0, y_max = 1.0;
    
    // Tailles totales
    m.Nbpt = (Nx + 1) * (Ny + 1);
    m.Nbtri = 2 * Nx * Ny;

    // Allocation mémoire optimisée avec Eigen
    // resize ne contient pas de données, contrairement à setZero qui initialise à 0
    // c'est mieux car on évite que le pc déplace les données en mémoire à chaque fois qu'on rajoute une ligne ou colonne
    // il faut aussi faire attention avec resize car il alloue une place en mémoire mais n'efface pas les données qu'il y avait dedans, pour faire la total on fait : m.Coorneu.setZero(m.Nbpt, 2), on remplace les ancienens valeurs par des zeros ensuite on remplit notre matrice avec une boucle for
    m.Coorneu.resize(m.Nbpt, 2); // on alloue en mémoire une matrice de Nbpt lignes et 2 colonnes pour X,Y car coorneu est une matrice
    m.Refneu.setZero(m.Nbpt); // Initialise tout à 0 (intérieur)
    
    m.Numtri.resize(m.Nbtri, 3); // la même, on a Nbtri lignes et 3 colonnes pour les 3 sommets de chaque triangle, il contient les indices des sommets de chaque triangle (point 0 : (0,0), point 1 : (1,0), point 2 : (0,1) pour le triangle 0, etc.) | on n'écrit pas les coordonnées dans Numtri, on écrit les indices des sommets, et ensuite avec Coorneu on peut récupérer les coordonnées exactes de ces sommets
    m.Reftri.setOnes(m.Nbtri); // vecteur Nbtri lignes et 1 colonne, initialisé à 1 (Air), tout le four est remplis d'air | on met tout à 1 par défaut, et ensuite on va mettre à 2 les triangles qui sont dans la zone de résine, ça nous permets de différencier les matériaux dans la suite du code, notamment pour l'assemblage de la matrice de rigidité et du système inverse, on va utiliser Reftri pour savoir si on est dans l'air ou dans la résine et appliquer les conductivités correspondantes

    // Pas d'espace
    double dx = (x_max - x_min) / Nx; // discrétisation uniforme en x, basique
    double dy = (y_max - y_min) / Ny; // la même discrétisation uniforme en y

    // 1. Création des Sommets et étiquetage des Bords
    // boucleO(Nx * Ny) = O(N), mieux que les truc de python car on peut vectoriser tout ça avec le compilateur C++ avec l'option -03 et utiliser les instructions SIMD du processeur pour faire les calculs en parallèle, contrairement à python qui est plus lent et ne peut pas faire ça, même avec des bibliothèques comme numpy, c'est pour ça que le C++ est souvent utilisé pour les calculs intensifs en ingénierie
    int index_noeud = 0; // compteur les boucles vont parcourir le maillage en 2D, mais on doit remplir les matrices en 1D (pour Coorneu et Refneu), donc on utilise un index qui va de 0 à Nbpt-1 pour remplir les données ligne par ligne
    for (int j = 0; j <= Ny; ++j) {
        for (int i = 0; i <= Nx; ++i) {
            double x = x_min + i * dx; // on part du min et on ajoute i fois le pas pour trouver la coordonnée x du nœud
            double y = y_min + j * dy; // la même pour y, on part du min et on ajoute j fois le pas pour trouver la coordonnée y du nœud
            
            // Accès avec parenthèses pour les matrices Eigen
            // on range x dans la première colonne de Coorneu, et y dans la deuxième colonne, pour chaque ligne index_noeud qui correspond à un nœud du maillage, on a les coordonnées x et y de ce nœud dans Coorneu, c'est plus facile pour la suite du code de faire comme ça plutôt que d'avoir deux vecteurs séparés pour x et y
            m.Coorneu(index_noeud, 0) = x; 
            m.Coorneu(index_noeud, 1) = y;
            
            // Étiquetage des frontières selon les données du projet
            if (y >= y_max - 1e-9) {
                m.Refneu(index_noeud) = 1; // on regarde si le sommet est tout en haut, si oui, c'est un Bord haut, Dirichlet TD = 50
            } else if (y <= y_min + 1e-9) {
                m.Refneu(index_noeud) = 2; // si le sommet est tout en bas, bord bas, Dirichlet TD = 100
            } else if (x <= x_min + 1e-9 || x >= x_max - 1e-9) {
                m.Refneu(index_noeud) = 3; // ni en haut ni en bas, on regarde s'il est à gauche ou à droite, Bords latéraux, Neumann
            }
            
            index_noeud++;
        }
    }

    // 2. Création des Triangles et étiquetage des Matériaux
    // Pour cela on découpe notre domaine en carré, puis on fait une diagonale
    int index_tri = 0;
    for (int j = 0; j < Ny; ++j) { // on utilise l'inégalité stricte car on ne fait pas sur les points mais sur les espaces entre les points (20 pts sur un axe : 19 points entre les extremes)
        for (int i = 0; i < Nx; ++i) {
            // Indices des 4 coins de la maille rectangulaire
            int n_bas_gauche  = j * (Nx + 1) + i; // on était en vision 2D mais on est en 1D, pour faire la converstion on prends le nombre totale de points sur une ligne qui est Nx + 1 et pour aller à la ligne j on doit sauter j fois les rangés de points, ça c'est pour une ranger, mais comme on a plusieurs rangé, on rajoute i pour le prochain
            int n_bas_droite  = n_bas_gauche + 1; // dans le vecteur le point bas droite est juste à côté du point bas gauche, donc on ajoute 1 pour avoir son indice
            int n_haut_gauche = (j + 1) * (Nx + 1) + i; // même principe, on doit juste bien avoir en tête la représentation du vecteur
            int n_haut_droite = n_haut_gauche + 1;
            
            // Calcul du centre de la maille pour déterminer le matériau
            double centre_x = x_min + (i + 0.5) * dx; 
            double centre_y = y_min + (j + 0.5) * dy;
            
            // L'objet en résine est au centre : [-0.5, 0.5] x [-0.2, 0.2]
            int ref_materiau = 1; // 1 = Air
            if (centre_x >= -0.5 && centre_x <= 0.5 && centre_y >= -0.2 && centre_y <= 0.2) {
                ref_materiau = 2; // 2 = Résine
            }

            // Triangle 1 (inférieur)
            m.Numtri(index_tri, 0) = n_bas_gauche;
            m.Numtri(index_tri, 1) = n_bas_droite;
            m.Numtri(index_tri, 2) = n_haut_gauche;
            m.Reftri(index_tri) = ref_materiau;
            index_tri++;

            // Triangle 2 (supérieur)
            m.Numtri(index_tri, 0) = n_bas_droite;
            m.Numtri(index_tri, 1) = n_haut_droite;
            m.Numtri(index_tri, 2) = n_haut_gauche;
            m.Reftri(index_tri) = ref_materiau;
            index_tri++;
        }
    }

    return m;
}

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
    // 3 noeuds i, 3 noeuds j, matrice de tailel 9, et un triangle a 3 sommets, donc 9 contributions non nulles par triangle, on réserve donc Nbtri * 9 pour éviter les reallocations dynamiques de la liste de triplets qui peuvent être coûteuses en temps
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


// Fonction pour trouver le triangle le plus proche d'une coordonnée (x, y)
int trouver_triangle_proche(const Maillage& m, double cible_x, double cible_y) {
    int meilleur_triangle = 0;
    double distance_min = 1e15; // Un chiffre très grand au départ

    for (int k = 0; k < m.Nbtri; ++k) {
        // Récupérer les sommets du triangle k
        int n1 = m.Numtri(k, 0);
        int n2 = m.Numtri(k, 1);
        int n3 = m.Numtri(k, 2);

        // Calculer le centre (barycentre) du triangle
        double centre_x = (m.Coorneu(n1, 0) + m.Coorneu(n2, 0) + m.Coorneu(n3, 0)) / 3.0;
        double centre_y = (m.Coorneu(n1, 1) + m.Coorneu(n2, 1) + m.Coorneu(n3, 1)) / 3.0;

        // Calculer la distance entre le centre du triangle et la cible
        double distance = std::sqrt(std::pow(centre_x - cible_x, 2) + std::pow(centre_y - cible_y, 2));

        // Si c'est le plus proche trouvé jusqu'à présent, on le garde
        if (distance < distance_min) {
            distance_min = distance;
            meilleur_triangle = k;
        }
    }

    return meilleur_triangle;
}



// --- Résolution du problème d'optimisation (Exercice 11.2) ---
// Version optimisée avec Eigen : plus compacte et lisible
void construire_systeme_opt(int nr, int Nbtri, const Maillage& m, 
                            const Eigen::VectorXd& T0, 
                            const Eigen::MatrixXd& Tres, 
                            double Tcui, 
                            Eigen::MatrixXd& Aopt, 
                            Eigen::VectorXd& bopt) {
    
    Aopt.setZero(nr, nr); 
    bopt.setZero(nr);     

    // 1. Matrice de masse élémentaire (standard P1) 
    Eigen::Matrix3d AK;
    AK << 2, 1, 1,
          1, 2, 1,
          1, 1, 2;
    AK /= 12.0;  // SANS LE DELTA QUI EST MESK : L'AIR

    for (int k = 0; k < Nbtri; ++k) {
        if (m.Reftri(k) == 2) { // Uniquement sur l'objet S | 'on calcule la somme uniquement sur les triangles internes à l'objet à cuire
            
            // Récupération des 3 indices du triangle k
            Eigen::Vector3i nodes = m.Numtri.row(k); 
            
            // Calcul du déterminant Delta (2 * Aire)
            double x1 = m.Coorneu(nodes(0), 0), y1 = m.Coorneu(nodes(0), 1);
            double x2 = m.Coorneu(nodes(1), 0), y2 = m.Coorneu(nodes(1), 1);
            double x3 = m.Coorneu(nodes(2), 0), y3 = m.Coorneu(nodes(2), 1);
            double Delta = std::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));

            // Températures locales T0 sur les 3 sommets du triangle
            Eigen::Vector3d t0_local;
            t0_local << T0(nodes(0)), T0(nodes(1)), T0(nodes(2));

            // Écart à la cible (Tcui - T0) localement
            Eigen::Vector3d target_diff = Eigen::Vector3d::Constant(Tcui) - t0_local;

            // 2. Assemblage Eigen-style
            for (int i = 0; i < nr; ++i) {
                // Signature de la résistance i sur les 3 sommets
                Eigen::Vector3d vi; // extrait les températures Tk(A1) etc produite par la résistance i sur les 3 sommets du triangle k, c'est la signature thermique de la résistance i sur ce triangle, elle nous dit comment la température change localement sur ce triangle lorsque la résistance i est activée, c'est essentiel pour calculer l'influence de cette résistance sur le champ de température global et pour construire le système d'optimisation qui nous permettra de trouver les puissances optimales des résistances pour atteindre la température cible dans la zone de résine
                vi << Tres(nodes(0), i), Tres(nodes(1), i), Tres(nodes(2), i);

                // Second membre b_opt (Produit scalaire via la matrice de masse)
                bopt(i) += Delta * vi.dot(AK * target_diff); // c'est Y^T * M * (Tcui - T0), on calcule d'abord le produit matriciel AK * target_diff qui nous donne un vecteur de taille 3, ensuite on fait le produit scalaire avec vi pour obtenir la contribution de la résistance i à l'écart à la cible sur ce triangle, et on multiplie par Delta pour prendre en compte la taille du triangle, c'est la contribution locale de ce triangle à b_opt(i), et en sommant sur tous les triangles de la résine, on obtient le second membre global b_opt(i) qui reflète l'influence de la résistance i sur l'écart à la cible dans toute la zone de résine

                for (int j = 0; j < nr; ++j) {
                    // Signature de la résistance j sur les 3 sommets
                    Eigen::Vector3d vj;
                    vj << Tres(nodes(0), j), Tres(nodes(1), j), Tres(nodes(2), j);

                    // Matrice A_opt (Interaction i et j) 
                    Aopt(i, j) += Delta * vi.dot(AK * vj); // c'est Y^T * M * Y
                }
            }
        }
    }
}

// ============================================================================
// FONCTIONS HELPER POUR STRUCTURER LE MAIN - Rendent le code lisible
// ============================================================================

Eigen::VectorXd computeBaseTemperature(const Maillage& mesh, 
                                       const Eigen::SparseMatrix<double>& stiffness,
                                       double T_top, double T_bottom) {
    
    // Assemblage de la matrice de rigidité globale
    Eigen::SparseMatrix<double> A = stiffness;
    Eigen::VectorXd b = Eigen::VectorXd::Zero(mesh.Nbpt);

    // Conditions aux limites (Dirichlet) : T_top au sommet, T_bottom en bas [cite: 210, 248]
    elim(A, b, mesh.Refneu, 1, T_top); // 50°C en haut
    elim(A, b, mesh.Refneu, 2, T_bottom); // 100°C en bas

    // Résolution du système linéaire, on cherche T0
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver; // A = L * D * L^T, L la lower traingulaire inf, D la diagonale, L^T matrice triangu supp, low transpose
    solver.compute(A); // décomposition de la matrice A 
    Eigen::VectorXd T0 = solver.solve(b); // $$L \times D \times L^T \times T_0 = b$$ | résoud  L*y = b , D * z = y, L^T * T0 = z
    // on sépare le solve et compute car on devra faire le calculs pour plusieurs resistance, donc le vecteur b va changer à chaque fois, mais la matrice A reste la même, donc on peut faire la décomposition une seule fois avec compute, et ensuite faire plusieurs solve pour différents b, c'est plus efficace que de refaire la décomposition à chaque fois

    return T0;
}

Eigen::MatrixXd computeResistanceSignatures(const Maillage& mesh,
                                            const Eigen::SparseMatrix<double>& stiffness,
                                            const std::vector<int>& resistor_indices) {
    int nr = resistor_indices.size();
    Eigen::MatrixXd Tres(mesh.Nbpt, nr);

    // 1. Préparation de la matrice (UNE SEULE FOIS)
    Eigen::SparseMatrix<double> A_signatures = stiffness;
    Eigen::VectorXd dummy_b = Eigen::VectorXd::Zero(mesh.Nbpt);
    
    // On modifie la matrice pour imposer 0 sur les bords 
    elim(A_signatures, dummy_b, mesh.Refneu, 1, 0.0);
    elim(A_signatures, dummy_b, mesh.Refneu, 2, 0.0);

    // 2. Factorisation (LE CALCUL LOURD, UNE SEULE FOIS)
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A_signatures);

    // 3. Boucle ultra-rapide sur les résistances
    for (int k = 0; k < nr; ++k) {
        
        // Construction du second membre brut pour la résistance k
        Eigen::VectorXd b_k = calculer_b_re(mesh.Nbpt, resistor_indices[k], mesh);
        
        // Au lieu d'appeler elim() qui referait toute la matrice, 
        // on applique manuellement l'effet de T=0 sur le vecteur b_k.
        // Mathématiquement, ça revient juste à forcer les cases des bords à 0.
        for (int i = 0; i < mesh.Nbpt; ++i) {
            if (mesh.Refneu(i) == 1 || mesh.Refneu(i) == 2) {
                b_k(i) = 0.0;
            }
        }

        // Résolution instantanée ! (Car le solver a déjà fait compute)
        Tres.col(k) = solver.solve(b_k);
    }

    return Tres;
}

void exportTemperatureField(const std::string& filename,
                            const Maillage& mesh,
                            const Eigen::VectorXd& temperature) {
    std::ofstream file(filename);
    file << "x,y,T\n";
    for (int i = 0; i < mesh.Nbpt; ++i) {
        file << mesh.Coorneu(i, 0) << "," << mesh.Coorneu(i, 1) << "," 
             << temperature(i) << "\n";
    }
    file.close();
}

void exportMeshGeometry(const std::string& filename, const Maillage& mesh) {
    std::ofstream file(filename);
    file << "n1,n2,n3\n";
    for (int i = 0; i < mesh.Nbtri; ++i) {
        file << mesh.Numtri(i, 0) << "," << mesh.Numtri(i, 1) << "," 
             << mesh.Numtri(i, 2) << "\n";
    }
    file.close();
}

Eigen::VectorXd solveInverseProblem(const Maillage& mesh,
                                    const Eigen::VectorXd& T0,
                                    const Eigen::MatrixXd& Tres,
                                    double target_temperature) {
    int nr = Tres.cols();

    // Construction du système d'optimisation
    Eigen::MatrixXd A_opt(nr, nr);
    Eigen::VectorXd b_opt(nr);
    construire_systeme_opt(nr, mesh.Nbtri, mesh, T0, Tres, target_temperature, A_opt, b_opt);

    // Résolution du petit système (nr x nr) pour trouver les puissances optimales
    Eigen::VectorXd alphas = A_opt.ldlt().solve(b_opt);

    return alphas;
}

// ============================================================================
// VÉRIFICATION DE LA TEMPÉRATURE DE CUISSON - Zone de Résine
// ============================================================================

double computeResinAverageTemperature(const Maillage& mesh,
                                       const Eigen::VectorXd& temperature,
                                       double target_temperature) {
    // Zone de la résine définie en coordonnées : [-0.5, 0.5] x [-0.2, 0.2]
    double x_min = -0.5, x_max = 0.5;
    double y_min = -0.2, y_max = 0.2;
    
    // Accumulation des températures et compte des nœuds dans la résine
    double sum_temperature = 0.0;
    int count_nodes = 0;

    // Parcourir tous les nœuds pour trouver ceux dans la zone de résine
    for (int i = 0; i < mesh.Nbpt; ++i) {
        double x = mesh.Coorneu(i, 0);
        double y = mesh.Coorneu(i, 1);

        // Vérifier si le nœud est dans la zone de résine
        if (x >= x_min && x <= x_max && y >= y_min && y <= y_max) {
            sum_temperature += temperature(i);
            count_nodes++;
        }
    }

    // Calcul de la température moyenne
    double avg_temp = (count_nodes > 0) ? sum_temperature / count_nodes : 0.0;

    // === AFFICHAGE DU DIAGNOSTIC ===
    std::cout << "\n========================================\n";
    std::cout << "  VERIFICATION CUISSON DE LA RESINE\n";
    std::cout << "========================================\n";
    std::cout << "Zone de résine : [-0.5, 0.5] x [-0.2, 0.2]\n";
    std::cout << "Nœuds détectés : " << count_nodes << "\n";
    std::cout << "Température moyenne : " << avg_temp << "°C\n";
    std::cout << "Température cible    : " << target_temperature << "°C\n";
    std::cout << "Écart                : " << (avg_temp - target_temperature) << "°C\n";

    // Diagnostic de qualité
    double tolerance = 5.0;  // Tolérance acceptable ±5°C
    double error_percent = std::abs(avg_temp - target_temperature) / target_temperature * 100.0;

    std::cout << "\n  Diagnostic : ";
    if (std::abs(avg_temp - target_temperature) < tolerance) {
        std::cout << "✓ EXCELLENTE - Cuisson optimale!\n";
    } else if (error_percent < 10.0) {
        std::cout << "✓ BON - Cuisson acceptable\n";
    } else if (avg_temp < target_temperature) {
        std::cout << "⚠ INSUFFISANT - Augmenter la puissance des résistances\n";
    } else {
        std::cout << "⚠ EXCESSIF - Diminuer la puissance des résistances\n";
    }
    std::cout << "========================================\n\n";

    return avg_temp;
}