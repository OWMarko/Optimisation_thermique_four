#include <iostream>
#include <vector>
#include <iomanip>
#include "maillage.hpp"
#include "solver_direct.hpp"
#include "optimisation.hpp"
#include "export.hpp"
#include "experiences.hpp"

int main() {
    std::cout << std::fixed << std::setprecision(2);

    std::cout << "\n========================================================\n";
    std::cout << "      OPTIMISATION DE LA TEMPERATURE D'UN FOUR \n";
    std::cout << "========================================================\n";

    // Maillage du domaine

    int nx = 40, ny = 40; 
    std::cout << "\n GENERATION DU MAILLAGE\n";
    std::cout << "    -> Resolution : " << nx << "x" << ny << " elements\n";
    
    Maillage four = creer_maillage(nx, ny);
    
    std::cout << "    -> Noeuds     : " << four.Nbpt << "\n";
    std::cout << "    -> Triangles  : " << four.Nbtri << "\n";

    exportMeshGeometry("../data/triangles.csv", four);

    // Informations sur les propriétés physiques et conditions limites

    std::cout << "\n PROPRIETES PHYSIQUES ET CONDITIONS LIMITES\n";
    double k_air = 1.0;
    double k_resin = 10.0;
    std::cout << "    -> Conductivite Air   : " << k_air << " W/m.K\n";
    std::cout << "    -> Conductivite Resine: " << k_resin << " W/m.K\n";
    
    Eigen::VectorXd conductivite(3); // vecteur contenant les différentes conductivités pour les éléments (0: air, 1: résine, 2: résistance)
    conductivite << 0.0, k_air, k_resin;

    auto matrice_elementaire = assa(four.Nbpt,
                                    four.Nbtri,
                                    four.Coorneu,
                                    four.Refneu, 
                                    four.Numtri,
                                    four.Reftri,
                                    conductivite); // matrice de rigidité globale du four

    double T_haut = 50.0;
    double T_bas = 100.0;
    std::cout << "    -> T_paroi_Haut       : " << T_haut << " C\n";
    std::cout << "    -> T_paroi_Bas        : " << T_bas << " C\n";
    
    Eigen::VectorXd T0 = computeBaseTemperature(four, matrice_elementaire, T_haut, T_bas); // vecteur température de base sans résistances
    exportTemperatureField("../data/noeuds_T0.csv", four, T0); 

    // Ajout des résistances et calcul de leurs signatures thermiques

    std::cout << "\n CONFIGURATION DES RESISTANCES\n";

    // Positions des résistances (6 résistances : 3 en haut, 3 en bas), arbitraire mais symétrique pour tester l'optimisation
    std::vector<std::pair<double, double>> pos_res =
    {
        {-0.7,  0.7}, {0.0,  0.7}, {0.7,  0.7}, // Haut : Coin gauche, Bord milieu, Coin droit
        {-0.5, -0.5}, {0.0, -0.5}, {0.5, -0.5}  // Bas  : Coin gauche, Bord milieu, Coin droit
    };
    
    // Pour une meilleure lisibilité, on affiche les positions et les indices des triangles où sont placées les résistances
    std::vector<int> indices;
    for (size_t i = 0; i < pos_res.size(); ++i)
    {
        int idx = trouver_triangle_proche(four, pos_res[i].first, pos_res[i].second);
        indices.push_back(idx);
        std::cout << "    -> Resistance " << i+1 << " : position (" << pos_res[i].first 
                  << ", " << pos_res[i].second << ") | Triangle : " << idx << "\n";
    }

    Eigen::MatrixXd signatures = computeResistanceSignatures(four, matrice_elementaire, indices); // matrice des signatures thermiques de chaque résistance (chaque colonne correspond à une résistance)
    exportTemperatureField("../data/noeuds_Tk_1.csv", four, signatures.col(0)); 

    std::cout << "\n CALCUL DU PROBLEME DIRECT | TEST MANUEL \n";
    double alpha_manuel = 25000.0;
    std::cout << "    -> Puissance uniforme : " << alpha_manuel << " W\n";
    Eigen::VectorXd T_direct = T0 + signatures * Eigen::VectorXd::Constant(indices.size(), alpha_manuel); // vecteur température

    exportTemperatureField("../data/noeuds_T_direct.csv", four, T_direct);

    // Résolution du problème inverse pour trouver les puissances optimales des résistances
    std::cout << "\n RESOLUTION DU PROBLEME INVERSE\n";
    double T_cible = 250.0;
    double C_test = 1e-9;
    std::cout << "    -> Temperature Cible : " << T_cible << " C\n";
    std::cout << "    -> Penalite : C = " << C_test << "\n";

    Eigen::VectorXd a_brut = solveInverseProblem(four, T0, signatures, T_cible);
    Eigen::VectorXd a_pen  = solveInverseProblemPenalite(four, T0, signatures, T_cible, C_test);

    // Génération de la courbe en L pour analyser la stabilité de la solution optimisée en fonction de la pénalité C 
    std::cout << "\n ANALYSE DE STABILITE (COURBE EN L)\n";
    generer_courbe_L(four, T0, signatures, T_cible, "../data/l_curve.csv");
    std::cout << "    -> Etude de sensibilite exportee dans data/l_curve.csv\n";

    // Résumé des résultats des puissances optimales pour chaque résistance avec et sans pénalité et notre cas manuel
    std::cout << "\n    RESULTATS DES PUISSANCES (W) :\n";
    std::cout << "    --------------------------------------------\n";
    std::cout << "    ID  |   SANS PENALITE   |   AVEC PENALITE  \n";
    std::cout << "    --------------------------------------------\n";
    for (size_t i = 0; i < indices.size(); ++i) {
        std::cout << "    " << i+1 << "   |   " << std::setw(13) << a_brut(i) 
                  << "   |   " << std::setw(13) << a_pen(i) << "\n";
    }
    std::cout << "    --------------------------------------------\n";

    Eigen::VectorXd T_brut = T0 + signatures * a_brut; // Vecteur température optimisée sans pénalité
    Eigen::VectorXd T_pen  = T0 + signatures * a_pen; // Vecteur température optimisée avec pénalité

    std::cout << "\n RESULTATS FINAUX \n";
    
    // 1) Notre cas manuel avec une puissance uniforme pour toutes les résistances
    std::cout << "    --- CAS DIRECT (PUISSANCE MANUELLE = 25 000 W) ---\n";
    computeResinAverageTemperature(four, T_direct, T_cible);

    // 2) Notre cas optimisé sans pénalité (solution brute du problème inverse)
    std::cout << "\n    --- CAS SANS PENALITE (OPTIMISE) ---\n";
    computeResinAverageTemperature(four, T_brut, T_cible);
    exportTemperatureField("../data/noeuds_T_opt_brut.csv", four, T_brut); 

    // 3) Notre cas optimisé avec pénalité (solution régularisée du problème inverse)
    std::cout << "\n    --- CAS AVEC PENALITE (REGULARISE) ---\n";
    computeResinAverageTemperature(four, T_pen, T_cible);
    exportTemperatureField("../data/noeuds_T_opt_pen.csv", four, T_pen);

    return 0;
}