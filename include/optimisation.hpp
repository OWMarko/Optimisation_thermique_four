#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include "maillage.hpp"


// --- Problème Inverse (Exercice 11.2) ---
// Calcule la matrice A_opt et le vecteur b_opt pour l'optimisation [cite: 1582, 1583]
void construire_systeme_opt(int nr, int Nbtri, const Maillage& m, 
                            const Eigen::VectorXd& T0, 
                            const Eigen::MatrixXd& Tres, 
                            double T_opt, 
                            Eigen::MatrixXd& A_opt, 
                            Eigen::VectorXd& b_opt);

// Résout le problème inverse et retourne les puissances optimales
Eigen::VectorXd solveInverseProblem(const Maillage& mesh,
                                    const Eigen::VectorXd& T0,
                                    const Eigen::MatrixXd& Tres,
                                    double target_temperature);

// 2. CAS RÉGULARISÉ : Résout le problème inverse avec pénalité de Tikhonov (C)
// Permet d'éviter les résistances négatives ou trop extrêmes
Eigen::VectorXd solveInverseProblemPenalite(const Maillage& mesh,
                                            const Eigen::VectorXd& T0,
                                            const Eigen::MatrixXd& Tres,
                                            double target_temperature,
                                            double C);


// Calcule et affiche la température moyenne de la résine
// La résine se trouve dans la zone [-0.5, 0.5] x [-0.2, 0.2]
// Retourne la température moyenne et affiche un diagnostic vs température cible
double computeResinAverageTemperature(const Maillage& mesh,
                                       const Eigen::VectorXd& temperature,
                                       double target_temperature = 250.0);