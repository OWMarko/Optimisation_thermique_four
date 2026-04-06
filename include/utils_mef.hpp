#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>

// La structure basée sur les conventions du livre, 100% Eigen
struct Maillage {
    int Nbpt;     // Nombre de sommets
    int Nbtri;    // Nombre de triangles

    // Tableaux de données sous forme de matrices et vecteurs Eigen
    Eigen::MatrixXd Coorneu; // Matrice des coordonnées X,Y
    Eigen::VectorXi Refneu;  // Vecteur des références des bords
    Eigen::MatrixXi Numtri;  // Matrice de connectivité des triangles
    Eigen::VectorXi Reftri;  // Vecteur des références des milieux (Air/Résine)
};

// --- Déclarations calquées sur le PDF ---

// 1. Génération du maillage (Exercice 11.1)
Maillage creer_maillage(int Nx, int Ny); 

// 2. Assemblage de la matrice de rigidité
Eigen::SparseMatrix<double> assa(int Nbpt, int Nbtri, 
                                 const Eigen::MatrixXd& Coorneu, 
                                 const Eigen::VectorXi& Refneu, 
                                 const Eigen::MatrixXi& Numtri, 
                                 const Eigen::VectorXi& Reftri, 
                                 const Eigen::VectorXd& Conduc);

// 3. Prise en compte des conditions de Dirichlet
void elim(Eigen::SparseMatrix<double>& Atotal, Eigen::VectorXd& belim, 
          const Eigen::VectorXi& Refneu, int Refdir, double Valdir);

// 4. Assemblage du problème inverse
void asst(int Nbpt, int Nbtri, const Maillage& m, 
          Eigen::MatrixXd& Aopt, Eigen::VectorXd& bopt, 
          const Eigen::VectorXd& TO, const Eigen::MatrixXd& Tres, double Tcui);

// Dans include/utils_mef.hpp, ajoute cette ligne si elle n'y est pas :
Eigen::VectorXd calculer_b_re(int Nbpt, int triangle_index, const Maillage& m);

// --- Problème Inverse (Exercice 11.2) ---
// Calcule la matrice A_opt et le vecteur b_opt pour l'optimisation [cite: 233, 234]
void construire_systeme_opt(int nr, int Nbtri, const Maillage& m, 
                            const Eigen::VectorXd& T0, 
                            const Eigen::MatrixXd& Tres, 
                            double T_opt, 
                            Eigen::MatrixXd& A_opt, 
                            Eigen::VectorXd& b_opt);