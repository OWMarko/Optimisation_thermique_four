#include <iostream>
#include <fstream>
#include <vector>
#include "../include/utils_mef.hpp"

int main() {
    std::cout << "--- DEMARRAGE DU PROJET FOUR THERMIQUE ---" << std::endl;
    
    // 1. Génération du maillage (Ex 11.1 - Q1)
    // On utilise une résolution de 20x20 pour une meilleure précision visuelle
    Maillage four = creer_maillage(20, 20);
    std::cout << "-> Maillage cree : " << four.Nbpt << " sommets, " << four.Nbtri << " triangles." << std::endl;
    
    // 2. Définition des conductivités [cite: 248]
    // Indice 1: Air (c=1), Indice 2: Résine (c=10)
    Eigen::VectorXd Conduc(3);
    Conduc << 0.0, 1.0, 10.0;
    
    // 3. Assemblage de la matrice de rigidité globale (Ex 11.1 - Q2 et Q3)
    Eigen::SparseMatrix<double> Atotal = assa(four.Nbpt, four.Nbtri, 
                                              four.Coorneu, four.Refneu, 
                                              four.Numtri, four.Reftri, 
                                              Conduc);
    std::cout << "-> Matrice de rigidite assemblee." << std::endl;
                                              
    // 4. Calcul du champ de température T0 (sans résistances) [cite: 246, 247]
    Eigen::VectorXd b0 = Eigen::VectorXd::Zero(four.Nbpt);
    Eigen::SparseMatrix<double> A0 = Atotal;

    // Conditions aux limites pour T0 : 50 au sommet, 100 en bas [cite: 210, 248]
    elim(A0, b0, four.Refneu, 1, 50.0);
    elim(A0, b0, four.Refneu, 2, 100.0);

    std::cout << "-> Calcul de T0 (base thermique)..." << std::endl;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solveur;
    solveur.compute(A0);
    Eigen::VectorXd T0 = solveur.solve(b0);

    // 5. Export de la géométrie et de T0 pour Python
    std::cout << "-> Export des donnees (T0)..." << std::endl;
    std::ofstream f_noeuds("../data/noeuds_T0.csv");
    f_noeuds << "x,y,T\n";
    for(int i = 0; i < four.Nbpt; ++i) {
        f_noeuds << four.Coorneu(i, 0) << "," << four.Coorneu(i, 1) << "," << T0(i) << "\n";
    }
    f_noeuds.close();

    std::ofstream f_tri("../data/triangles.csv");
    f_tri << "n1,n2,n3\n";
    for(int i = 0; i < four.Nbtri; ++i) {
        f_tri << four.Numtri(i, 0) << "," << four.Numtri(i, 1) << "," << four.Numtri(i, 2) << "\n";
    }
    f_tri.close();

    // 6. CALCUL DES SIGNATURES DES RÉSISTANCES (Tk)
    // On place 4 résistances à des indices manuels répartis dans le four
    std::vector<int> indices_resistors = {0, 20, 38, 380, 400, 798};
    int nr = indices_resistors.size();
    
    // Matrice pour stocker chaque champ Tk en colonne [cite: 245]
    Eigen::MatrixXd Tres(four.Nbpt, nr);
    
    for (int k = 0; k < nr; ++k) {
        std::cout << "-> Calcul de la signature Tk pour la resistance " << k+1 << "..." << std::endl;
        
        Eigen::VectorXd b_k = calculer_b_re(four.Nbpt, indices_resistors[k], four);
        Eigen::SparseMatrix<double> A_k = Atotal;
        
        // Pour les signatures, la température aux bords est imposée à 0 [cite: 251]
        elim(A_k, b_k, four.Refneu, 1, 0.0);
        elim(A_k, b_k, four.Refneu, 2, 0.0);
        
        solveur.compute(A_k);
        Tres.col(k) = solveur.solve(b_k);
    }

    // 7. CALCUL DE LA FIGURE 11.5 (Combinaison linéaire)
    // On multiplie chaque signature par la valeur de la puissance alpha 
    double alpha = 25000.0;
    Eigen::VectorXd T_tot = T0; 

    for (int k = 0; k < nr; ++k) {
        T_tot += alpha * Tres.col(k);
    }

    // Export des données pour la Figure 11.5
    std::ofstream f_res("../data/noeuds_T_res.csv");
    f_res << "x,y,T\n";
    for(int i = 0; i < four.Nbpt; ++i) {
        f_res << four.Coorneu(i, 0) << "," << four.Coorneu(i, 1) << "," << T_tot(i) << "\n";
    }
    f_res.close();

    std::cout << "\n--- RESULTATS ---" << std::endl;
    std::cout << "T0 Min/Max : " << T0.minCoeff() << " / " << T0.maxCoeff() << std::endl;
    std::cout << "T_total Max (Figure 11.5) : " << T_tot.maxCoeff() << " degres." << std::endl;
    std::cout << "Fichiers CSV generes dans /data." << std::endl;

    // --- RÉSOLUTION DU PROBLÈME INVERSE ---
    std::cout << "-> Resolution du probleme inverse (Optimisation)..." << std::endl;
    
    Eigen::MatrixXd A_opt(nr, nr);
    Eigen::VectorXd b_opt(nr);
    double T_opt = 250.0; // Température cible [cite: 280]

    construire_systeme_opt(nr, four.Nbtri, four, T0, Tres, T_opt, A_opt, b_opt);

    // Résolution du petit système (nr x nr) pour trouver les alphas
    Eigen::VectorXd alphas = A_opt.ldlt().solve(b_opt);

    std::cout << "-> Valeurs optimales des resistances : " << std::endl;
    for(int k=0; k<nr; ++k) std::cout << "   Alpha " << k+1 << " = " << alphas(k) << std::endl;

    // Calcul du champ final par combinaison linéaire [cite: 224, 362]
    Eigen::VectorXd T_final = T0;
    for (int k = 0; k < nr; ++k) {
        T_final += alphas(k) * Tres.col(k);
    }

    // Export pour visualisation (Figure 11.7 / 11.8)
    std::ofstream f_opt("../data/noeuds_T_opt.csv");
    f_opt << "x,y,T\n";
    for(int i = 0; i < four.Nbpt; ++i) {
        f_opt << four.Coorneu(i, 0) << "," << four.Coorneu(i, 1) << "," << T_final(i) << "\n";
    }
    f_opt.close();

    // Export de la signature de la résistance 1 (Figure 11.6)
    std::ofstream f_sig("../data/noeuds_Tk_1.csv");
    f_sig << "x,y,T\n";
    for(int i = 0; i < four.Nbpt; ++i) {
    f_sig << four.Coorneu(i, 0) << "," << four.Coorneu(i, 1) << "," << Tres(i, 0) << "\n";
    }
    f_sig.close();

    return 0;
}