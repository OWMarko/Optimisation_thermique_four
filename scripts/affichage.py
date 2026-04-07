import pandas as pd
import matplotlib
matplotlib.use('Qt5Agg') 
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import os
import numpy as np

if not os.path.exists('../img'): os.makedirs('../img')

# 1. Lecture des données (Assure-toi que le C++ exporte bien ces deux fichiers)
print("-> Lecture des données...")
df_T0      = pd.read_csv('../data/noeuds_T0.csv')
df_Tdir    = pd.read_csv('../data/noeuds_T_direct.csv')
df_Tk      = pd.read_csv('../data/noeuds_Tk_1.csv')
df_Tbrut   = pd.read_csv('../data/noeuds_T_opt_brut.csv') # Sans pénalité
df_Treg    = pd.read_csv('../data/noeuds_T_opt_reg.csv')  # Avec pénalité
df_l       = pd.read_csv('../data/l_curve.csv')
triangles  = pd.read_csv('../data/triangles.csv')

maillage = tri.Triangulation(df_T0['x'], df_T0['y'], triangles.values)

fig = plt.figure(figsize=(22, 11))
plt.suptitle("Analyse Comparative : Optimisation avec et sans Régularisation", fontsize=16, fontweight='bold')

# MODIFICATION ICI : zlim=None par défaut pour ne pas écraser les graphes
def plot_3d(pos, data, title, cmap='jet', zlim=None):
    ax = fig.add_subplot(pos, projection='3d')
    ax.plot_trisurf(maillage, data, cmap=cmap, edgecolor='none')
    ax.set_title(title, fontweight='bold')
    if zlim: ax.set_zlim(zlim)
    ax.view_init(elev=25, azim=-45)

# --- GRILLE DE RÉSULTATS ---
# Échelles adaptées pour bien voir les pentes et les bosses
plot_3d(231, df_T0['T'], '1. T0 (Base)', zlim=(0, 150))
plot_3d(232, df_Tdir['T'], '2. T_direct (Manuel)', zlim=(0, 200))
plot_3d(233, df_Tk['T'], '3. Tk (Signature)', cmap='viridis', zlim=None)

# --- COMPARAISON DES DEUX OPTIMISATIONS ---
# On fixe à 300 pour bien mettre en évidence les plateaux à 250°C
plot_3d(234, df_Tbrut['T'], '4. T_opt (SANS pénalité)', zlim=(0, 300))
plot_3d(235, df_Treg['T'], '5. T_opt (AVEC pénalité)', zlim=(0, 300))

# --- COURBE EN L POUR JUSTIFIER LE CHOIX ---
ax6 = fig.add_subplot(236)
ax6.loglog(df_l['residue'], df_l['norm_alpha'], 'o-b', label='Courbe en L')
ax6.set_title('6. Sensibilité : Courbe en L', fontweight='bold')
ax6.set_xlabel('Résidu (Erreur)')
ax6.set_ylabel('Norme Alpha (Puissance)')
ax6.grid(True, which="both", alpha=0.3)
ax6.legend()

plt.tight_layout()
plt.savefig('../img/Resultats_Complets_Four.png')
print("-> Rapport comparatif généré.")
plt.show()