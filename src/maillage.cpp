#include "../include/maillage.hpp"
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
            int n_bas_gauche  = j * (Nx + 1) + i; // on était en vision 2D mais on est en 1D, pour faire la converstion on prends le nombre totale de points sur une ligne qui est Nx + 1 et pour aller à la ligne j on doit sauter j fois les rangés de points, ça cest pour une ranger, mais comme on a plusieurs rangé, on rajoute i pour le prochain
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