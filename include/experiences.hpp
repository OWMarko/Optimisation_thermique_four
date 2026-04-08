#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>
#include "maillage.hpp"

// On souhaite créer une fonction qui balaye plusieurs valeurs de C pénalité
// pour analyser la stabilité de la solution optimisée en fonction de C et trouver le "coude" de la courbe en L 
// Pour ce faire on exporte les résultats dans un CSV contenant C, residue et norm_a, pour faire la courbe

void generer_courbe_L(const Maillage& mesh,
                      const Eigen::VectorXd& T0,
                      const Eigen::MatrixXd& Tres,
                      double target_temp,
                      const std::string& filename);

