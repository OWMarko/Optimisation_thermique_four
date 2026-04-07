#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>
#include "maillage.hpp"

// Calcule la courbe en L (L-curve) pour une plage de paramètres C
// Exporte les résultats dans un CSV pour analyse
void generer_courbe_L(const Maillage& mesh,
                      const Eigen::VectorXd& T0,
                      const Eigen::MatrixXd& Tres,
                      double target_temp,
                      const std::string& filename);

// Calcule l'erreur L2 relative entre la solution et la cible sur l'objet
double calculer_erreur_objet(const Maillage& mesh,
                             const Eigen::VectorXd& T_final,
                             double target_temp);