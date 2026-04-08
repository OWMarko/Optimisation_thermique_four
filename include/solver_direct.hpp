#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <string>
#include "maillage.hpp"

// Assemblage de la matrice de rigidité
Eigen::SparseMatrix<double> assa(int Nbpt,
                                 int Nbtri, 
                                 const Eigen::MatrixXd& Coorneu, 
                                 const Eigen::VectorXi& Refneu, 
                                 const Eigen::MatrixXi& Numtri, 
                                 const Eigen::VectorXi& Reftri, 
                                 const Eigen::VectorXd& Conduc);

// Prise en compte des conditions de Dirichlet
void elim(Eigen::SparseMatrix<double>& Atotal,
          Eigen::VectorXd& belim, 
          const Eigen::VectorXi& Refneu,
          int Refdir,
          double Valdir);

// Calcule le second membre pour une résistance
Eigen::VectorXd calculer_b_re(int Nbpt,
                              int triangle_index,
                              const Maillage& m);

// Calcule la température de base T0 sans résistances avec T_haut = 50 et T_bas = 100
Eigen::VectorXd computeBaseTemperature(const Maillage& mesh, 
                                       const Eigen::SparseMatrix<double>& matrice_elementaire,
                                       double T_top,
                                       double T_bottom);

// Calcule les signatures thermiques Tk pour chaque résistance
Eigen::MatrixXd computeResistanceSignatures(const Maillage& mesh,
                                            const Eigen::SparseMatrix<double>& matrice_elementaire,
                                            const std::vector<int>& resistor_indices);