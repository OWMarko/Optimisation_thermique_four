#include "../include/utils_mef.hpp"

Maillage creer_maillage(int Nx, int Ny) {
    Maillage m;
    
    // Dimensions du four : carré de -1 à 1
    double x_min = -1.0, x_max = 1.0;
    double y_min = -1.0, y_max = 1.0;
    
    // Tailles totales
    m.Nbpt = (Nx + 1) * (Ny + 1);
    m.Nbtri = 2 * Nx * Ny;

    // Allocation mémoire optimisée avec Eigen
    m.Coorneu.resize(m.Nbpt, 2);
    m.Refneu.setZero(m.Nbpt); // Initialise tout à 0 (intérieur)
    
    m.Numtri.resize(m.Nbtri, 3);
    m.Reftri.setOnes(m.Nbtri); // Initialise tout à 1 (Air par défaut)

    // Pas d'espace
    double dx = (x_max - x_min) / Nx;
    double dy = (y_max - y_min) / Ny;

    // 1. Création des Sommets et étiquetage des Bords
    int index_noeud = 0;
    for (int j = 0; j <= Ny; ++j) {
        for (int i = 0; i <= Nx; ++i) {
            double x = x_min + i * dx;
            double y = y_min + j * dy;
            
            // Accès avec parenthèses pour les matrices Eigen
            m.Coorneu(index_noeud, 0) = x;
            m.Coorneu(index_noeud, 1) = y;
            
            // Étiquetage des frontières selon les données du projet
            if (y >= y_max - 1e-6) {
                m.Refneu(index_noeud) = 1; // Bord haut, Dirichlet TD = 50
            } else if (y <= y_min + 1e-6) {
                m.Refneu(index_noeud) = 2; // Bord bas, Dirichlet TD = 100
            } else if (x <= x_min + 1e-6 || x >= x_max - 1e-6) {
                m.Refneu(index_noeud) = 3; // Bords latéraux, Neumann
            }
            
            index_noeud++;
        }
    }

    // 2. Création des Triangles et étiquetage des Matériaux
    int index_tri = 0;
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            // Indices des 4 coins de la maille rectangulaire
            int n_bas_gauche  = j * (Nx + 1) + i;
            int n_bas_droite  = n_bas_gauche + 1;
            int n_haut_gauche = (j + 1) * (Nx + 1) + i;
            int n_haut_droite = n_haut_gauche + 1;
            
            // Calcul du centre de la maille pour déterminer le matériau
            double centre_x = x_min + (i + 0.5) * dx;
            double centre_y = y_min + (j + 0.5) * dy;
            
            // L'objet en résine est au centre : [-0.5, 0.5] x [-0.2, 0.2]
            int ref_materiau = 1; // 1 = Air
            if (centre_x >= -0.5 && centre_x <= 0.5 && centre_y >= -0.2 && centre_y <= 0.2) {
                ref_materiau = 2; // 2 = Résine
            }

            // Triangle 1 (inférieur)
            m.Numtri(index_tri, 0) = n_bas_gauche;
            m.Numtri(index_tri, 1) = n_bas_droite;
            m.Numtri(index_tri, 2) = n_haut_gauche;
            m.Reftri(index_tri) = ref_materiau;
            index_tri++;

            // Triangle 2 (supérieur)
            m.Numtri(index_tri, 0) = n_bas_droite;
            m.Numtri(index_tri, 1) = n_haut_droite;
            m.Numtri(index_tri, 2) = n_haut_gauche;
            m.Reftri(index_tri) = ref_materiau;
            index_tri++;
        }
    }

    return m;
}

// --- Assemblage de la matrice de rigidité (Exercice 11.1 - Q2 et Q3) ---
Eigen::SparseMatrix<double> assa(int Nbpt, int Nbtri, 
                                 const Eigen::MatrixXd& Coorneu, 
                                 const Eigen::VectorXi& Refneu, 
                                 const Eigen::MatrixXi& Numtri, 
                                 const Eigen::VectorXi& Reftri, 
                                 const Eigen::VectorXd& Conduc) {
    
    // En C++, on utilise une liste de triplets pour construire une matrice creuse efficacement
    std::vector<Eigen::Triplet<double>> triplets;
    
    // Optimisation R&D : on réserve la mémoire (9 termes par triangle)
    triplets.reserve(Nbtri * 9); 
    
    // Boucle sur les éléments (les triangles)
    for (int k = 0; k < Nbtri; ++k) {
        
        // Numéros des 3 sommets du triangle k
        int n1 = Numtri(k, 0);
        int n2 = Numtri(k, 1);
        int n3 = Numtri(k, 2);
        
        // Coordonnées X et Y des 3 sommets
        double x1 = Coorneu(n1, 0), y1 = Coorneu(n1, 1);
        double x2 = Coorneu(n2, 0), y2 = Coorneu(n2, 1);
        double x3 = Coorneu(n3, 0), y3 = Coorneu(n3, 1);
        
        // Le déterminant Delta (qui vaut 2 fois l'aire du triangle)
        double Delta = std::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
        
        // La matrice du gradient Dp
        // On divise directement par Delta ici pour simplifier les calculs suivants
        double Dp[2][3] = {
            {(y2 - y3)/Delta, (y3 - y1)/Delta, (y1 - y2)/Delta},
            {(x3 - x2)/Delta, (x1 - x3)/Delta, (x2 - x1)/Delta}
        };
        
        // Conductivité du matériau pour ce triangle
        double c = Conduc(Reftri(k)); 
        
        // Calcul des tableaux élémentaires et assemblage
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                // Formule de l'intégrale du gradient : c * Delta * (Dp' * Dp) / 2
                double val = c * Delta * (Dp[0][i] * Dp[0][j] + Dp[1][i] * Dp[1][j]) / 2.0;
                
                // On stocke la valeur (ligne, colonne, valeur)
                triplets.push_back(Eigen::Triplet<double>(Numtri(k, i), Numtri(k, j), val));
            }
        }
    }
    
    // Création de la matrice creuse globale d'ordre Nbpt
    Eigen::SparseMatrix<double> Atotal(Nbpt, Nbpt);
    Atotal.setFromTriplets(triplets.begin(), triplets.end());
    
    return Atotal;
}

// --- Prise en compte des conditions de Dirichlet (Exercice 11.1 - Q4) ---
void elim(Eigen::SparseMatrix<double>& Atotal, Eigen::VectorXd& belim, 
          const Eigen::VectorXi& Refneu, int Refdir, double Valdir) {
          
    // Méthode de pénalisation : on force l'équation avec un chiffre immense
    double tgv = 1e15; 
    
    // On parcourt tous les sommets
    for (int i = 0; i < Refneu.size(); ++i) {
        // Si le sommet appartient à la frontière que l'on veut fixer (Refdir)
        if (Refneu(i) == Refdir) {
            // On écrase la diagonale et le second membre
            Atotal.coeffRef(i, i) = tgv;
            belim(i) = Valdir * tgv;
        }
    }
}

// --- Calcul du second membre pour une résistance (Exercice 11.1 - Q3) ---
// On utilise la formulation Delta/6.0 conforme au document du professeur
Eigen::VectorXd calculer_b_re(int Nbpt, int triangle_index, const Maillage& m) {
    // Initialisation d'un vecteur b rempli de zéros
    Eigen::VectorXd b = Eigen::VectorXd::Zero(Nbpt);
    
    // Indices des 3 sommets du triangle où se trouve la résistance
    int n1 = m.Numtri(triangle_index, 0);
    int n2 = m.Numtri(triangle_index, 1);
    int n3 = m.Numtri(triangle_index, 2);
    
    // Coordonnées pour calculer le Déterminant (Delta)
    double x1 = m.Coorneu(n1, 0), y1 = m.Coorneu(n1, 1);
    double x2 = m.Coorneu(n2, 0), y2 = m.Coorneu(n2, 1);
    double x3 = m.Coorneu(n3, 0), y3 = m.Coorneu(n3, 1);
    
    // Calcul de Delta (2 * Aire du triangle)
    // Delta = |(x2-x1)(y3-y1) - (x3-x1)(y2-y1)|
    double Delta = std::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
    
    // Contribution élémentaire selon le PDF : (F * Delta) / 6
    // On considère F = 1 pour calculer la signature thermique unitaire
    double contribution = Delta / 6.0;
    
    // On ajoute la contribution aux trois sommets du triangle
    b(n1) += contribution;
    b(n2) += contribution;
    b(n3) += contribution;
    
    return b;
}

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
    AK /= 12.0; // On divise par 12 car le prof multiplie par Delta (2*Aire) [cite: 383]

    for (int k = 0; k < Nbtri; ++k) {
        if (m.Reftri(k) == 2) { // Uniquement sur l'objet S 
            
            // Récupération des 3 indices du triangle k
            Eigen::Vector3i nodes = m.Numtri.row(k); 
            
            // Calcul du déterminant Delta (2 * Aire) [cite: 383]
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

                // Second membre b_opt (Produit scalaire via la matrice de masse) [cite: 402]
                bopt(i) += Delta * vi.dot(AK * target_diff);

                for (int j = 0; j < nr; ++j) {
                    // Signature de la résistance j sur les 3 sommets
                    Eigen::Vector3d vj;
                    vj << Tres(nodes(0), j), Tres(nodes(1), j), Tres(nodes(2), j);

                    // Matrice A_opt (Interaction i et j) [cite: 392]
                    Aopt(i, j) += Delta * vi.dot(AK * vj);
                }
            }
        }
    }
}