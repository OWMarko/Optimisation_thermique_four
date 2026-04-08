#include "experiences.hpp"
#include "optimisation.hpp" // Pour appeler construire_systeme_opt
#include <fstream>
#include <iostream>
#include <cmath>

void generer_courbe_L(const Maillage& mesh,
                      const Eigen::VectorXd& T0,
                      const Eigen::MatrixXd& Tres,
                      double target_temp,
                      const std::string& filename) {
    
    int nr = Tres.cols(); // on grep le nombre de résistances à partir du nombre de colonnes de la matrice des signatures thermiques
    Eigen::MatrixXd A_opt(nr, nr); // on crée notre matrice A 
    Eigen::VectorXd b_opt(nr); 
    
    // On assemble le système de base une seule fois
    construire_systeme_opt(nr, mesh.Nbtri, mesh, T0, Tres, target_temp, A_opt, b_opt);

    std::ofstream file(filename);
    file << "C,residue,norm_alpha,error_resin\n";

    // On souhaite tester plusieurs valeurs de C, pour cela on crée une boucle basique
    for (double C = 1e-12; C <= 0.1; C *= 10.0)
    {
        Eigen::MatrixXd A_reg = A_opt + C * Eigen::MatrixXd::Identity(nr, nr); // notre système régularisé avec pénalité
        Eigen::VectorXd alpha_C = A_reg.ldlt().solve(b_opt); // on résoud le système 

        // On calcule le résidu et la norme alpha comme dans notre formule
        double residue = (A_opt * alpha_C - b_opt).norm();
        double norm_a = alpha_C.norm();
                
        file << C << "," << residue << "," << norm_a << "\n";
    }
    file.close();
}