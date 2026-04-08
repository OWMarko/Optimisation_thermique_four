#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include "maillage.hpp"


// Calcule la matrice A_opt et le vecteur b_opt
void construire_systeme_opt(int nr, int Nbtri, const Maillage& m, 
                            const Eigen::VectorXd& T0, 
                            const Eigen::MatrixXd& Tres, 
                            double T_opt, 
                            Eigen::MatrixXd& A_opt, 
                            Eigen::VectorXd& b_opt);

// Fonction qui résout le problème inverse et retourne les puissances optimales
Eigen::VectorXd solveInverseProblem(const Maillage& mesh,
                                    const Eigen::VectorXd& T0,
                                    const Eigen::MatrixXd& Tres,
                                    double target_temperature);

// Fonction qui résout le problème inverse avec régularisation et retourne les puissances optimales
Eigen::VectorXd solveInverseProblemPenalite(const Maillage& mesh,
                                            const Eigen::VectorXd& T0,
                                            const Eigen::MatrixXd& Tres,
                                            double target_temperature,
                                            double C);


// Fonction qui calcule et affiche la température moyenne de la résine
double computeResinAverageTemperature(const Maillage& mesh,
                                       const Eigen::VectorXd& temperature,
                                       double target_temperature = 250.0);