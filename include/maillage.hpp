#pragma once
#include <Eigen/Dense>
#include <vector>

// Structure de données pour le maillage du four, 
// nous souhaitons stocker les informations sur les sommets, les triangles, et les références pour les conditions limites et les matériaux
struct Maillage {
    int Nbpt;     // Nombre de sommets
    int Nbtri;    // Nombre de triangles

    // Tableaux de données sous forme de matrices et vecteurs Eigen
    Eigen::MatrixXd Coorneu; // Matrice des coordonnées des sommets (Nbpt * 2) | Coorneu(0,0) = x du sommet 0, Coorneu(0,1) = y du sommet 0,
    Eigen::VectorXi Refneu; // Vecteur des références des sommets (pour les conditions limites, 0 = intérieur, 1 = bord haut, 2 = bord bas)
    Eigen::MatrixXi Numtri; // Matrice des indices des sommets pour chaque triangle (Nbtri * 3) | Numtri(0,0) = indice du sommet 1 du triangle 0, etc.
    Eigen::VectorXi Reftri; // Vecteur des références des milieux (Air/Résine)
};

// Fonction pour créer le maillage du four, avec les étiquetages des bords et des matériaux
Maillage creer_maillage(int Nx, int Ny); 

// Fonction qui nous servira à trouver l'indice du triangle le plus proche d'une coordonnée (x, y), pour placer les résistances dans le maillage
int trouver_triangle_proche(const Maillage& m, double cible_x, double cible_y);