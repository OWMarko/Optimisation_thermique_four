#include <iostream>
#include <vector>
#include <iomanip>
#include "maillage.hpp"
#include "solver_direct.hpp"
#include "optimisation.hpp"
#include "export.hpp"
#include "experiences.hpp"

int main() {
    // --- Config Affichage ---
    std::cout << std::fixed << std::setprecision(2);

    std::cout << "\n========================================================\n";
    std::cout << "      PROJET R&D : OPTIMISATION THERMIQUE DU FOUR\n";
    std::cout << "========================================================\n";

    // ========== ÉTAPE 1 : MAILLAGE ==========
    int nx = 30, ny = 30; // Maillage fin
    std::cout << "\n GENERATION DU MAILLAGE\n";
    std::cout << "    -> Resolution : " << nx << "x" << ny << " elements\n";
    
    Maillage four = creer_maillage(nx, ny);
    
    std::cout << "    -> Noeuds     : " << four.Nbpt << "\n";
    std::cout << "    -> Triangles  : " << four.Nbtri << "\n";
    exportMeshGeometry("../data/triangles.csv", four);

    // ========== ÉTAPE 2 : PHYSIQUE DU DOMAINE ==========
    std::cout << "\n PROPRIETES PHYSIQUES ET CL\n";
    double k_air = 1.0, k_resin = 10.0;
    std::cout << "    -> Conductivite Air   : " << k_air << " W/m.K\n";
    std::cout << "    -> Conductivite Resine: " << k_resin << " W/m.K\n";
    
    Eigen::VectorXd conductivities(3);
    conductivities << 0.0, k_air, k_resin;
    auto stiffness_matrix = assa(four.Nbpt, four.Nbtri, four.Coorneu, four.Refneu, 
                                 four.Numtri, four.Reftri, conductivities);

    double T_haut = 50.0, T_bas = 100.0;
    std::cout << "    -> T_paroi_Haut       : " << T_haut << " C\n";
    std::cout << "    -> T_paroi_Bas        : " << T_bas << " C\n";
    
    Eigen::VectorXd T0 = computeBaseTemperature(four, stiffness_matrix, T_haut, T_bas);
    exportTemperatureField("../data/noeuds_T0.csv", four, T0); // Pour Fig 11.4

    // ========== ÉTAPE 3 : CONFIGURATION DES RESISTANCES ==========
    std::cout << "\n CONFIGURATION DES RESISTANCES\n";
    std::vector<std::pair<double, double>> pos_res = {
        {-0.9,  0.9}, {0.0,  0.9}, {0.9,  0.9}, // Haut : Coin gauche, Bord milieu, Coin droit
        {-0.9, -0.9}, {0.0, -0.9}, {0.9, -0.9}  // Bas  : Coin gauche, Bord milieu, Coin droit
    };
    
    std::vector<int> indices;
    for (size_t i = 0; i < pos_res.size(); ++i) {
        int idx = trouver_triangle_proche(four, pos_res[i].first, pos_res[i].second);
        indices.push_back(idx);
        std::cout << "    -> Resistance " << i+1 << " : position (" << pos_res[i].first 
                  << ", " << pos_res[i].second << ") | Triangle : " << idx << "\n";
    }

    Eigen::MatrixXd signatures = computeResistanceSignatures(four, stiffness_matrix, indices);
    exportTemperatureField("../data/noeuds_Tk_1.csv", four, signatures.col(0)); // Pour Fig 11.6

    // ========== ÉTAPE 4 : PROBLÈME DIRECT (Le test manuel !) ==========
    std::cout << "\n CALCUL DU PROBLEME DIRECT (Test Manuel)\n";
    double alpha_manuel = 25000.0;
    std::cout << "    -> Puissance uniforme : " << alpha_manuel << " W\n";
    Eigen::VectorXd T_direct = T0 + signatures * Eigen::VectorXd::Constant(indices.size(), alpha_manuel);
    exportTemperatureField("../data/noeuds_T_direct.csv", four, T_direct); // Pour Fig 11.5

    // ========== ÉTAPE 5 : ANALYSE DE RIGUEUR (COURBE EN L) ==========
    std::cout << "\n ANALYSE DE STABILITE (COURBE EN L)\n";
    double T_cible = 250.0;
    generer_courbe_L(four, T0, signatures, T_cible, "../data/l_curve.csv");
    std::cout << "    -> Etude de sensibilite exportee dans data/l_curve.csv\n"; // Pour Fig 6

    // ========== ÉTAPE 6 : OPTIMISATION ET COMPARAISON ==========
    std::cout << "\n RESOLUTION DU PROBLEME INVERSE\n";
    double C = 1e-9;
    std::cout << "    -> Temperature Cible : " << T_cible << " C\n";
    std::cout << "    -> Penalite Tikhonov : C = " << C << "\n";

    Eigen::VectorXd a_brut = solveInverseProblem(four, T0, signatures, T_cible);
    Eigen::VectorXd a_reg  = solveInverseProblemPenalite(four, T0, signatures, T_cible, C);

    // --- TABLEAU DES RÉSULTATS ---
    std::cout << "\n    RESULTATS DES PUISSANCES (W) :\n";
    std::cout << "    --------------------------------------------\n";
    std::cout << "    ID  |   SANS PENALITE   |   AVEC TIKHONOV   \n";
    std::cout << "    --------------------------------------------\n";
    for (size_t i = 0; i < indices.size(); ++i) {
        std::cout << "    " << i+1 << "   |   " << std::setw(13) << a_brut(i) 
                  << "   |   " << std::setw(13) << a_reg(i) << "\n";
    }
    std::cout << "    --------------------------------------------\n";

    // ========== ÉTAPE 7 : VALIDATION ET EXPORT ==========
    Eigen::VectorXd T_brut = T0 + signatures * a_brut;
    Eigen::VectorXd T_reg  = T0 + signatures * a_reg;

    std::cout << "\n DIAGNOSTIC FINAL DE CUISSON\n";
    
    // 1. Diagnostic du cas manuel (T_direct)
    std::cout << "    --- CAS DIRECT (PUISSANCE MANUELLE = 25 000 W) ---\n";
    computeResinAverageTemperature(four, T_direct, T_cible);

    // 2. Diagnostic du cas optimisé sans pénalité
    std::cout << "\n    --- CAS SANS PENALITE (OPTIMISE) ---\n";
    computeResinAverageTemperature(four, T_brut, T_cible);
    exportTemperatureField("../data/noeuds_T_opt_brut.csv", four, T_brut); // Pour Fig 4

    // 3. Diagnostic du cas optimisé avec Tikhonov
    std::cout << "\n    --- CAS AVEC TIKHONOV (REGULARISE) ---\n";
    computeResinAverageTemperature(four, T_reg, T_cible);
    exportTemperatureField("../data/noeuds_T_opt_reg.csv", four, T_reg); // Pour Fig 5

    std::cout << "\n========================================================\n";
    std::cout << "            SIMULATION TERMINEE AVEC SUCCES\n";
    std::cout << "========================================================\n\n";

    return 0;
}