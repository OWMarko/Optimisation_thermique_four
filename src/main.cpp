#include <iostream>
#include <vector>
#include "../include/utils_mef.hpp"

int main() {
    std::cout << "\n========================================\n";
    std::cout << "  SIMULATION THERMIQUE - FOUR (MEF)\n";
    std::cout << "========================================\n\n";

    // ========== ÉTAPE 1 : MAILLAGE ==========
    std::cout << "[1/4] GENERATION DU MAILLAGE...\n";
    Maillage four = creer_maillage(30, 30);
    
    // Assemblage de la matrice de rigidité (faite 1 seule fois)
    Eigen::VectorXd conductivities(3);
    conductivities << 0.0, 1.0, 10.0;
    Eigen::SparseMatrix<double> stiffness_matrix = assa(four.Nbpt, four.Nbtri,
                                                         four.Coorneu, four.Refneu, 
                                                         four.Numtri, four.Reftri, 
                                                         conductivities);

    // ========== ÉTAPE 2 : TEMPÉRATURE DE BASE (T0) ==========
    std::cout << "[2/4] CALCUL DE T0 (cas de base sans résistances)...\n";
    
    Eigen::VectorXd T0 = computeBaseTemperature(four, stiffness_matrix, 50.0, 100.0);
    
    exportMeshGeometry("../data/triangles.csv", four);
    exportTemperatureField("../data/noeuds_T0.csv", four, T0);

    // ========== ÉTAPE 3 : SIGNATURES DES RÉSISTANCES (Tk) ==========
    std::cout << "[3/4] CALCUL DES SIGNATURES THERMIQUES (Tk)...\n";
    
    // Définir les positions (x, y) souhaitées pour les 6 résistances
    // Exemple : 3 en haut de l'objet, 3 en bas
    std::vector<std::pair<double, double>> positions_souhaitees = {
        {-0.6,  0.4},  {0.0,  0.4},  {0.6,  0.4},  // Ligne du haut
        {-0.6, -0.4},  {0.0, -0.4},  {0.6, -0.4}   // Ligne du bas
    };

    // Trouver les indices des triangles correspondants
    std::vector<int> indices_resistors;
    for (const auto& pos : positions_souhaitees) {
        int tri_index = trouver_triangle_proche(four, pos.first, pos.second);
        indices_resistors.push_back(tri_index);
    }
    
    // Le reste du code ne change pas !
    Eigen::MatrixXd resistor_signatures = computeResistanceSignatures(four, stiffness_matrix, 
                                                                      indices_resistors);

    // ========== ÉTAPE 3.5 : PROBLÈME DIRECT (Test manuel) ==========
    std::cout << "[3.5/4] PROBLEME DIRECT (4 résistances manuelles à 25 000)...\n";
    
    // On fixe nous-mêmes la puissance (alpha) à 25 000 selon le PDF
    double alpha_manuel = 25000.0; 
    Eigen::VectorXd T_direct = T0;
    
    // Composition de la température par simple addition (superposition)
    for (int k = 0; k < (int)indices_resistors.size(); ++k) {
        T_direct += alpha_manuel * resistor_signatures.col(k);
    }
    
    std::cout << "      T_direct Max : " << T_direct.maxCoeff() << "°C\n";
    
    // On vérifie la température moyenne de la résine avec ce réglage arbitraire
    double resin_avg_direct = computeResinAverageTemperature(four, T_direct, 250.0);
    std::cout << "      Moyenne sur la résine : " << resin_avg_direct << "°C\n\n";

    // Export pour l'affichage (Tu pourras le comparer à la figure 11.5 du PDF !)
    exportTemperatureField("../data/noeuds_T_direct.csv", four, T_direct);

    // ========== ÉTAPE 4 : PROBLÈME INVERSE - OPTIMISATION ==========
    std::cout << "[4/4] RESOLUTION DU PROBLEME INVERSE...\n";
    
    double target_temp = 250.0;
    Eigen::VectorXd optimal_powers = solveInverseProblem(four, T0, resistor_signatures, target_temp);
    
    // Affichage des puissances optimales
    std::cout << "      Puissances optimales (W):\n";
    for (int k = 0; k < (int)indices_resistors.size(); ++k) {
        std::cout << "      Résistance " << (k+1) << " : " << optimal_powers(k) << " W\n";
    }
    
    // ASTUCE EIGEN : Composition finale instantanée sans boucle for !
    // T_final = T0 + sum(alpha_k * Tk)
    Eigen::VectorXd T_final = T0 + resistor_signatures * optimal_powers;

    // ========== VÉRIFICATION ET EXPORT ==========
    double resin_avg_temp = computeResinAverageTemperature(four, T_final, target_temp);
    exportTemperatureField("../data/noeuds_T_opt.csv", four, T_final);
    
    std::cout << "      Température moyenne de la résine : " << resin_avg_temp << "°C\n\n";

    std::cout << "========================================\n";
    std::cout << "  SIMULATION TERMINEE\n";
    std::cout << "========================================\n";

    return 0;
}