#pragma once
#include <Eigen/Dense>
#include <vector>

struct Maillage {
    int Nbpt;     // Nombre de sommets
    int Nbtri;    // Nombre de triangles

    // Tableaux de données sous forme de matrices et vecteurs Eigen
    Eigen::MatrixXd Coorneu; // Matrice des coordonnées X,Y, stock tout à la suite | Contient pour chaque ligne i les coo x_0 et y_0 cartésiennes du fours
    Eigen::VectorXi Refneu;  // Vecteur des références des bords
    Eigen::MatrixXi Numtri;  // Matrice de connectivité des triangles | "triangle num 5 est formé par les sommets 0 1 2", et ensuite avec Coorneu on prends les coo de ces sommets et on a les coo exacte du triangle
    Eigen::VectorXi Reftri;  // Vecteur des références des milieux (Air/Résine)
};

// 1. Génération du maillage (Exercice 11.1)
Maillage creer_maillage(int Nx, int Ny); 

// Trouve l'indice du triangle le plus proche d'une coordonnée (x, y)
int trouver_triangle_proche(const Maillage& m, double cible_x, double cible_y);