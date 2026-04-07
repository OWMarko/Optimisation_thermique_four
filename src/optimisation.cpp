#include "../include/optimisation.hpp"
#include <iostream>
#include <cmath>

// --- Résolution du problème d'optimisation (Exercice 11.2) ---
// Version optimisée avec Eigen : plus compacte et lisible
void construire_systeme_opt(int nr, int Nbtri, const Maillage& m, 
                            const Eigen::VectorXd& T0, 
                            const Eigen::MatrixXd& Tres, 
                            double Tcui, 
                            Eigen::MatrixXd& Aopt, 
                            Eigen::VectorXd& bopt) {
    
    Aopt.setZero(nr, nr); 
    bopt.setZero(nr);     

    // 1. Matrice de masse élémentaire (standard P1) 
    Eigen::Matrix3d AK;
    AK << 2, 1, 1,
          1, 2, 1,
          1, 1, 2;
    AK /= 12.0;  // SANS LE DELTA QUI EST MESK : L'AIR

    for (int k = 0; k < Nbtri; ++k) {
        if (m.Reftri(k) == 2) { // Uniquement sur l'objet S | 'on calcule la somme uniquement sur les triangles internes à l'objet à cuire
            
            // Récupération des 3 indices du triangle k
            Eigen::Vector3i nodes = m.Numtri.row(k); 
            
            // Calcul du déterminant Delta (2 * Aire)
            double x1 = m.Coorneu(nodes(0), 0), y1 = m.Coorneu(nodes(0), 1);
            double x2 = m.Coorneu(nodes(1), 0), y2 = m.Coorneu(nodes(1), 1);
            double x3 = m.Coorneu(nodes(2), 0), y3 = m.Coorneu(nodes(2), 1);
            double Delta = std::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));

            // Températures locales T0 sur les 3 sommets du triangle
            Eigen::Vector3d t0_local;
            t0_local << T0(nodes(0)), T0(nodes(1)), T0(nodes(2));

            // Écart à la cible (Tcui - T0) localement
            Eigen::Vector3d target_diff = Eigen::Vector3d::Constant(Tcui) - t0_local;

            // 2. Assemblage Eigen-style
            for (int i = 0; i < nr; ++i) {
                // Signature de la résistance i sur les 3 sommets
                Eigen::Vector3d vi; 
                vi << Tres(nodes(0), i), Tres(nodes(1), i), Tres(nodes(2), i);

                // Second membre b_opt (Produit scalaire via la matrice de masse)
                bopt(i) += Delta * vi.dot(AK * target_diff); 

                for (int j = 0; j < nr; ++j) {
                    // Signature de la résistance j sur les 3 sommets
                    Eigen::Vector3d vj;
                    vj << Tres(nodes(0), j), Tres(nodes(1), j), Tres(nodes(2), j);

                    // Matrice A_opt (Interaction i et j) 
                    Aopt(i, j) += Delta * vi.dot(AK * vj); 
                }
            }
        }
    }
}

Eigen::VectorXd solveInverseProblem(const Maillage& mesh,
                                    const Eigen::VectorXd& T0,
                                    const Eigen::MatrixXd& Tres,
                                    double target_temperature) {
    int nr = Tres.cols();

    // Construction du système d'optimisation
    Eigen::MatrixXd A_opt(nr, nr);
    Eigen::VectorXd b_opt(nr);
    construire_systeme_opt(nr, mesh.Nbtri, mesh, T0, Tres, target_temperature, A_opt, b_opt);

    // Résolution du petit système (nr x nr) pour trouver les puissances optimales
    Eigen::VectorXd alphas = A_opt.ldlt().solve(b_opt);

    return alphas;
}

// --- Résolution du problème inverse AVEC régularisation de Tikhonov ---
Eigen::VectorXd solveInverseProblemPenalite(const Maillage& mesh,
                                            const Eigen::VectorXd& T0,
                                            const Eigen::MatrixXd& Tres,
                                            double target_temperature,
                                            double C) {
    int nr = Tres.cols();

    // 1. Construction du système d'optimisation classique
    Eigen::MatrixXd A_opt(nr, nr);
    Eigen::VectorXd b_opt(nr);
    construire_systeme_opt(nr, mesh.Nbtri, mesh, T0, Tres, target_temperature, A_opt, b_opt);

    // 2. Ajout de la pénalité de Tikhonov sur la diagonale (C * I)
    // Cela modifie la matrice A pour pénaliser les fortes variations des coefficients alpha
    A_opt += C * Eigen::MatrixXd::Identity(nr, nr);

    // 3. Résolution du système régularisé (nr x nr)
    Eigen::VectorXd alphas = A_opt.ldlt().solve(b_opt);

    return alphas;
}


double computeResinAverageTemperature(const Maillage& mesh,
                                       const Eigen::VectorXd& temperature,
                                       double target_temperature) {
    // Zone de la résine définie en coordonnées : [-0.5, 0.5] x [-0.2, 0.2]
    double x_min = -0.5, x_max = 0.5;
    double y_min = -0.2, y_max = 0.2;
    
    // Accumulation des températures et compte des nœuds dans la résine
    double sum_temperature = 0.0;
    int count_nodes = 0;

    // Parcourir tous les nœuds pour trouver ceux dans la zone de résine
    for (int i = 0; i < mesh.Nbpt; ++i) {
        double x = mesh.Coorneu(i, 0);
        double y = mesh.Coorneu(i, 1);

        // Vérifier si le nœud est dans la zone de résine
        if (x >= x_min && x <= x_max && y >= y_min && y <= y_max) {
            sum_temperature += temperature(i);
            count_nodes++;
        }
    }

    // Calcul de la température moyenne
    double avg_temp = (count_nodes > 0) ? sum_temperature / count_nodes : 0.0;

    // === AFFICHAGE DU DIAGNOSTIC ===
    std::cout << "\n========================================\n";
    std::cout << "  VERIFICATION CUISSON DE LA RESINE\n";
    std::cout << "========================================\n";
    std::cout << "Zone de résine : [-0.5, 0.5] x [-0.2, 0.2]\n";
    std::cout << "Nœuds détectés : " << count_nodes << "\n";
    std::cout << "Température moyenne : " << avg_temp << "°C\n";
    std::cout << "Température cible    : " << target_temperature << "°C\n";
    std::cout << "Écart                : " << (avg_temp - target_temperature) << "°C\n";

    // Diagnostic de qualité
    double tolerance = 5.0;  // Tolérance acceptable ±5°C
    double error_percent = std::abs(avg_temp - target_temperature) / target_temperature * 100.0;

    std::cout << "\n  Diagnostic : ";
    if (std::abs(avg_temp - target_temperature) < tolerance) {
        std::cout << "✓ EXCELLENTE - Cuisson optimale!\n";
    } else if (error_percent < 10.0) {
        std::cout << "✓ BON - Cuisson acceptable\n";
    } else if (avg_temp < target_temperature) {
        std::cout << "⚠ INSUFFISANT - Augmenter la puissance des résistances\n";
    } else {
        std::cout << "⚠ EXCESSIF - Diminuer la puissance des résistances\n";
    }
    std::cout << "========================================\n\n";

    return avg_temp;
}