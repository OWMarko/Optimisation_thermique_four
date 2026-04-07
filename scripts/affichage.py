import pandas as pd
import matplotlib
matplotlib.use('Qt5Agg') 
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import os

# 1. Vérification du dossier de destination
dossier_img = '../img'
if not os.path.exists(dossier_img):
    os.makedirs(dossier_img)

# 2. Lecture des données exportées par le C++
print("-> Lecture de l'ensemble des données...")
df_T0   = pd.read_csv('../data/noeuds_T0.csv')      # Fig 11.4 
df_Tdir = pd.read_csv('../data/noeuds_T_direct.csv')# Fig 11.5 (Problème direct)
df_Tk   = pd.read_csv('../data/noeuds_Tk_1.csv')    # Fig 11.6 
df_Topt = pd.read_csv('../data/noeuds_T_opt.csv')   # Fig 11.7 / 11.8 
triangles = pd.read_csv('../data/triangles.csv')    # Connectivité 

# Vérification cohérence des données
print(f"   Nœuds T0: {len(df_T0)}, Triangles: {len(triangles)}, Nœuds T_opt: {len(df_Topt)}")

# 3. Création de la triangulation commune
maillage = tri.Triangulation(df_T0['x'], df_T0['y'], triangles.values)

# 4. Création d'une figure 2x3 (pour accueillir 5 graphiques)
fig = plt.figure(figsize=(22, 12)) # Légèrement élargie pour la 3ème colonne
plt.suptitle("Synthèse du Projet : Optimisation Thermique du Four", fontsize=16, fontweight='bold')

# --- FIG 11.4 : T0 (Base thermique) ---
ax1 = fig.add_subplot(231, projection='3d')
surf1 = ax1.plot_trisurf(maillage, df_T0['T'], cmap='jet', edgecolor='none')
ax1.set_title('Figure 11.4 : Sans résistance (T0)')
ax1.set_zlabel('T (°C)')
fig.colorbar(surf1, ax=ax1, shrink=0.5)

# --- FIG 11.5 : T_direct (Test Manuel avec 4 résistances) ---
ax2 = fig.add_subplot(232, projection='3d')
surf2 = ax2.plot_trisurf(maillage, df_Tdir['T'], cmap='jet', edgecolor='none')
ax2.set_title('Figure 11.5 : Problème Direct (α = 25 000)')
ax2.set_zlabel('T (°C)')
fig.colorbar(surf2, ax=ax2, shrink=0.5)

# --- FIG 11.6 : Tk (Signature unitaire) ---
ax3 = fig.add_subplot(233, projection='3d')
surf3 = ax3.plot_trisurf(maillage, df_Tk['T'], cmap='viridis', edgecolor='none')
ax3.set_title('Figure 11.6 : Signature d\'une résistance (Tk)')
ax3.set_zlabel('T (Unit)')
fig.colorbar(surf3, ax=ax3, shrink=0.5)

# --- FIG 11.8 : T_optimisée (Plateau 250°C avec 4 résistances) ---
ax4 = fig.add_subplot(234, projection='3d')
surf4 = ax4.plot_trisurf(maillage, df_Topt['T'], cmap='jet', edgecolor='none')
ax4.set_title('Figure 11.8 : Température Optimisée (Cible 250°C)')
ax4.set_zlabel('T (°C)')
ax4.set_zlim(0, 350)  # Pour bien voir le plateau
fig.colorbar(surf4, ax=ax4, shrink=0.5)

# --- FIG Comparaison : T0 vs T_optimisée ---
ax5 = fig.add_subplot(235, projection='3d')
diff = df_Topt['T'] - df_T0['T']
surf5 = ax5.plot_trisurf(maillage, diff, cmap='RdYlBu_r', edgecolor='none')
ax5.set_title('Différence : T_optimisée - T0')
ax5.set_zlabel('ΔT (°C)')
fig.colorbar(surf5, ax=ax5, shrink=0.5)

# Ajustement de la vue pour tous les subplots
for ax in [ax1, ax2, ax3, ax4, ax5]:
    ax.view_init(elev=25, azim=-45)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])

# 5. Sauvegarde et Affichage
chemin_save = os.path.join(dossier_img, 'Resultats_Complets_Four.png')
plt.savefig(chemin_save, dpi=300)
print(f"-> Rapport visuel complet sauvegardé dans '{chemin_save}'")

print("-> Ouverture de la fenêtre interactive...")
plt.show()