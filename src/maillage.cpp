#include "../include/maillage.hpp"
#include <fstream>
#include <cmath>
#include <iostream>

Maillage creer_maillage(int Nx,
                        int Ny) {
    Maillage m;
    
    // Dimensions du four : carré de -1 à 1
    double x_min = -1.0, x_max = 1.0;
    double y_min = -1.0, y_max = 1.0;
    
    // Tailles totales
    m.Nbpt = (Nx + 1) * (Ny + 1);
    m.Nbtri = 2 * Nx * Ny;

    // Allocation mémoire de nos objets
    m.Coorneu.resize(m.Nbpt, 2); 
    m.Refneu.setZero(m.Nbpt); 
    
    m.Numtri.resize(m.Nbtri, 3); 
    m.Reftri.setOnes(m.Nbtri); 

    // Pas d'espace
    double dx = (x_max - x_min) / Nx; // discrétisation uniforme en x
    double dy = (y_max - y_min) / Ny; 

    // 1) Création des sommets et étiquetage des bords
    int index_noeud = 0; 
    for (int j = 0; j <= Ny; ++j)
    {
        for (int i = 0; i <= Nx; ++i)
        {
            double x = x_min + i * dx; 
            double y = y_min + j * dy; 
            
            // On range les coordonnées dans la matrice Coorneu
            m.Coorneu(index_noeud, 0) = x; 
            m.Coorneu(index_noeud, 1) = y;
            
            // On regarde si le sommet est sur un bord 
            if (y >= y_max - 1e-9)
            {
                m.Refneu(index_noeud) = 1; // sommet en haut, bord haut, Dirichlet Td = 50
            } 
            else if (y <= y_min + 1e-9)
            {
                m.Refneu(index_noeud) = 2; // sommet en bas, bord bas, Dirichlet Td = 100
            }
            else if (x <= x_min + 1e-9 || x >= x_max - 1e-9)
            {
                m.Refneu(index_noeud) = 3; // sommet sur les bords gauche ou droit, bord isolant, Neumann Td = 0
            }
            
            index_noeud++;
        }
    }

    // 2) Création des triangles et étiquetage des matériaux
    int index_tri = 0;
    for (int j = 0; j < Ny; ++j)
    {
        for (int i = 0; i < Nx; ++i)
        {
            // Indices des 4 coins de la maille (i,j)
            int n_bas_gauche  = j * (Nx + 1) + i; 
            int n_bas_droite  = n_bas_gauche + 1; 
            int n_haut_gauche = (j + 1) * (Nx + 1) + i;
            int n_haut_droite = n_haut_gauche + 1;
            
            // Calcul du centre de la maille pour connaître le matériau dans le quel est le triangle
            double centre_x = x_min + (i + 0.5) * dx; 
            double centre_y = y_min + (j + 0.5) * dy;
            
            // L'objet en résine est en [-0.5, 0.5] x [-0.2, 0.2]
            int ref_materiau = 1; // air par défaut car majoritaire 
            if (centre_x >= -0.5 && centre_x <= 0.5 && centre_y >= -0.2 && centre_y <= 0.2)
            {
                ref_materiau = 2; // Résine
            }

            // Création du triangle 1 inférieur
            m.Numtri(index_tri, 0) = n_bas_gauche;
            m.Numtri(index_tri, 1) = n_bas_droite;
            m.Numtri(index_tri, 2) = n_haut_gauche;
            m.Reftri(index_tri) = ref_materiau;
            index_tri++;

            // Création du triangle 2 supérieur
            m.Numtri(index_tri, 0) = n_bas_droite;
            m.Numtri(index_tri, 1) = n_haut_droite;
            m.Numtri(index_tri, 2) = n_haut_gauche;
            m.Reftri(index_tri) = ref_materiau;
            index_tri++;
        }
    }

    return m;
}

// Fonction pour trouver le triangle le plus proche d'une coordonnée (x, y)
int trouver_triangle_proche(const Maillage& m,
                            double cible_x,
                            double cible_y) {
    int meilleur_triangle = 0;
    double distance_min = 1e15; 

    for (int k = 0; k < m.Nbtri; ++k)
    {
        // Récupérer les sommets du triangle k
        int n1 = m.Numtri(k, 0);
        int n2 = m.Numtri(k, 1);
        int n3 = m.Numtri(k, 2);

        // Calculer le centre (barycentre) du triangle
        double centre_x = (m.Coorneu(n1, 0) + m.Coorneu(n2, 0) + m.Coorneu(n3, 0)) / 3.0;
        double centre_y = (m.Coorneu(n1, 1) + m.Coorneu(n2, 1) + m.Coorneu(n3, 1)) / 3.0;

        // Calculer la distance entre le centre du triangle et la cible
        double distance = std::sqrt(std::pow(centre_x - cible_x, 2) + std::pow(centre_y - cible_y, 2));

        // Si c'est le plus proche trouvé jusqu'à présent, on le garde
        if (distance < distance_min)
        {
            distance_min = distance;
            meilleur_triangle = k;
        }
    }

    return meilleur_triangle;
}