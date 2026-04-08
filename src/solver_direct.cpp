#include "../include/solver_direct.hpp"
#include <vector>
#include <cmath>
#include <iostream>

// Assemblage de la matrice de rigidité
Eigen::SparseMatrix<double> assa(int Nbpt, int Nbtri, 
                                 const Eigen::MatrixXd& Coorneu, 
                                 const Eigen::VectorXi& Refneu, 
                                 const Eigen::MatrixXi& Numtri, 
                                 const Eigen::VectorXi& Reftri, 
                                 const Eigen::VectorXd& Conduc) {
    
    // On utilise la méthode des triplets pour construire la matrice creuse.
    // En même temps on réserve en mémoire 9 valeurs par triangle (3 noeuds i, 3 noeuds j)
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(Nbtri * 9); 
    
    for (int k = 0; k < Nbtri; ++k)
    {
        // On récupère les indices des 3 sommets du triangle k
        int n1 = Numtri(k, 0); 
        int n2 = Numtri(k, 1);
        int n3 = Numtri(k, 2);
        
        // On convertit les indices en coordonnées cartésiennes pour calculer le déterminant Delta et la matrice du gradient
        double x1 = Coorneu(n1, 0), y1 = Coorneu(n1, 1);
        double x2 = Coorneu(n2, 0), y2 = Coorneu(n2, 1);
        double x3 = Coorneu(n3, 0), y3 = Coorneu(n3, 1);
        
        // Le déterminant qui est égal à Air_Triangle * 2, peut être vu comme la jacobienne du changement de variable
        double Delta = std::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
        
        // La matrice du gradient
        Eigen::Matrix<double, 2, 3> Dp;

        // Remplissage par blocs ou par colonnes
        Dp << (y2 - y3), (y3 - y1), (y1 - y2),
              (x3 - x2), (x1 - x3), (x2 - x1);

        Dp /= Delta;
        
        // Conductivité du matériau pour le triangle k 
        double c = Conduc(Reftri(k));

        // Calcul de la matrice de rigidité locale 3x3
        //  c * Delta * (B(:,i)' * B(:,j)) / 2 
        Eigen::Matrix3d A_k = (c * Delta / 2.0) * (Dp.transpose() * Dp); 

        // Assemblage de la matrice globale via les triplets 
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            { 
                triplets.push_back(Eigen::Triplet<double>(Numtri(k, i), Numtri(k, j), A_k(i, j)));
            }
        }
    }
    
    // Construction finale de la matrice creuse à partir des triplets
    Eigen::SparseMatrix<double> Atotal(Nbpt, Nbpt); 
    Atotal.setFromTriplets(triplets.begin(), triplets.end());
    return Atotal;
}

// Calcul du second membre pour une résistance
Eigen::VectorXd calculer_b_re(int Nbpt,
                              int triangle_index,
                              const Maillage& m) {

    // Initialisation du second membre
    Eigen::VectorXd b = Eigen::VectorXd::Zero(Nbpt); 
    
    // Indices des 3 sommets du triangle où se trouve la résistance
    int n1 = m.Numtri(triangle_index, 0);
    int n2 = m.Numtri(triangle_index, 1);
    int n3 = m.Numtri(triangle_index, 2);
    
    // Coordonnées des 3 sommets du triangle
    double x1 = m.Coorneu(n1, 0), y1 = m.Coorneu(n1, 1);
    double x2 = m.Coorneu(n2, 0), y2 = m.Coorneu(n2, 1);
    double x3 = m.Coorneu(n3, 0), y3 = m.Coorneu(n3, 1);
    
    // Calcul du déterminant
    double Delta = std::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
    double contribution = Delta / 6.0; // notre constante pour chaque sommet du triangle
    
    // On ajoute la contribution aux trois sommets du triangle
    b(n1) += contribution;
    b(n2) += contribution;
    b(n3) += contribution;
    
    return b;
}

// Prise en compte des conditions de Dirichlet 
void elim(Eigen::SparseMatrix<double>& Atotal,
          Eigen::VectorXd& belim, 
          const Eigen::VectorXi& Refneu,
          int Refdir,
          double Valdir) {
    
    int n = Atotal.rows();
    
    // Etape 1 : Modification du second membre pour tenir compte de la contribution des nœuds de Dirichlet
    for (int k = 0; k < Atotal.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Atotal, k); it; ++it)
        {
            int row = it.row();
            int col = it.col();
            
            if (Refneu(col) == Refdir && Refneu(row) != Refdir)
            {
                belim(row) -= it.value() * Valdir; // b - A_2 * T_D
            }
        }
    }
    
    // Etape 2 : On parcourt la matrice pour savoir si le terme est dans A_2 ou A_3 ou A_2^T 
    for (int k = 0; k < Atotal.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Atotal, k); it; ++it)
        {
            int row = it.row();
            int col = it.col();
            
            // Si le terme touche une ligne ou une colonne du bord imposé
            if (Refneu(row) == Refdir || Refneu(col) == Refdir)
            {
                if (row == col)
                {
                    it.valueRef() = 1.0; // On connaît le comportement, condition Dirichlet : T_i = T_D = 50 ou 100
                } else
                {
                    it.valueRef() = 0.0; // On sait qu'il nous reste l'intéraction température va du libre vers le mur, condition Fourier, ici nulle, donc 0
                }
            }
        }
    }
    
    // Etape 3 : Imposition de la température T_D sur le bloc droit
    // La matrice est modifiée, il nous reste à appliquer la condition de Dirichlet 
    for (int i = 0; i < n; ++i)
    {
        if (Refneu(i) == Refdir)
        {
            belim(i) = Valdir; // Application de T_D
        }
    }
}

// Calcul des signatures thermiques Tk pour chaque résistance
Eigen::MatrixXd computeResistanceSignatures(const Maillage& mesh,
                                            const Eigen::SparseMatrix<double>& matrice_elementaire,
                                            const std::vector<int>& resistor_indices) {
    int nr = resistor_indices.size();
    Eigen::MatrixXd Tres(mesh.Nbpt, nr);

    // On crée notre matrice élementaire avec notre second membre avec influence des resistances
    Eigen::SparseMatrix<double> A_signatures = matrice_elementaire;
    Eigen::VectorXd b_simple = Eigen::VectorXd::Zero(mesh.Nbpt);
    
    // On modifie la matrice pour imposer 0 sur les bords
    elim(A_signatures, b_simple, mesh.Refneu, 1, 0.0);
    elim(A_signatures, b_simple, mesh.Refneu, 2, 0.0);

    // On factorise une seule fois la matrice des signatures, qui est la même pour toutes les résistances 
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A_signatures);

    if(solver.info() != Eigen::Success)
    {
        std::cerr << "Erreur: La factorisation de la matrice des signatures a échoué." << std::endl;
    }

    // On calcul les différentes signatures thermiques Tk pour chaque résistance k
    for (int k = 0; k < nr; ++k)
    {
        // Construction du second membre brut pour la résistance k
        Eigen::VectorXd b_k = calculer_b_re(mesh.Nbpt, resistor_indices[k], mesh);
        
        // Application manuelle de Dirichlet T=0 sur b_k pour rester cohérent avec la factorisation
        for (int i = 0; i < mesh.Nbpt; ++i)
        {
            if (mesh.Refneu(i) == 1 || mesh.Refneu(i) == 2)
            {
                b_k(i) = 0.0;
            }
        }

        // On place la signature thermique de la resistance k dans la colonne k de Tres
        Tres.col(k) = solver.solve(b_k);
    }

    return Tres;
}

// Calcule la température de base T0 (sans résistances) 
// Même logique que pour les signatures, mais avec un second membre nul et les conditions de Dirichlet imposées à 50 et 100
Eigen::VectorXd computeBaseTemperature(const Maillage& mesh, 
                                       const Eigen::SparseMatrix<double>& matrice_elementaire,
                                       double T_top,
                                       double T_bottom) {
    
    Eigen::SparseMatrix<double> A = matrice_elementaire;
    Eigen::VectorXd b = Eigen::VectorXd::Zero(mesh.Nbpt);

    // Conditions aux limites de Dirichlet
    elim(A, b, mesh.Refneu, 1, T_top); 
    elim(A, b, mesh.Refneu, 2, T_bottom); 

    // Résolution du système linéaire pour trouver la température de base T0
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A); 
    
    if(solver.info() != Eigen::Success) {
        std::cerr << "Erreur: La factorisation de la matrice T0 a échoué." << std::endl;
    }

    return solver.solve(b);
}