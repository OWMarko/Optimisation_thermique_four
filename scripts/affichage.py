import pandas as pd
import matplotlib
# Utilisation du moteur interactif Qt5 pour WSL
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
df_T0   = pd.read_csv('../data/noeuds_T0.csv')      # Fig 11.4 [cite: 247]
df_Tres = pd.read_csv('../data/noeuds_T_res.csv')   # Fig 11.5 
df_Tk   = pd.read_csv('../data/noeuds_Tk_1.csv')    # Fig 11.6 [cite: 257]
df_Topt = pd.read_csv('../data/noeuds_T_opt.csv')   # Fig 11.7 [cite: 281]
triangles = pd.read_csv('../data/triangles.csv')    # Connectivité [cite: 181]

# 3. Création de la triangulation commune
maillage = tri.Triangulation(df_T0['x'], df_T0['y'], triangles.values)

# 4. Création d'une figure 2x2
fig = plt.figure(figsize=(18, 12))
plt.suptitle("Synthèse du Projet : Optimisation Thermique du Four", fontsize=16)

# --- FIG 11.4 : T0 (Base thermique) ---
ax1 = fig.add_subplot(221, projection='3d')
surf1 = ax1.plot_trisurf(maillage, df_T0['T'], cmap='jet', edgecolor='none')
ax1.set_title('Figure 11.4 : Sans résistance (T0)')
ax1.set_zlabel('T (°C)')
fig.colorbar(surf1, ax=ax1, shrink=0.5)

# --- FIG 11.5 : T_total (Alpha = 25 000) ---
ax2 = fig.add_subplot(222, projection='3d')
surf2 = ax2.plot_trisurf(maillage, df_Tres['T'], cmap='jet', edgecolor='none')
ax2.set_title('Figure 11.5 : Réglage manuel (Alpha=25000)')
ax2.set_zlabel('T (°C)')
fig.colorbar(surf2, ax=ax2, shrink=0.5)

# --- FIG 11.6 : Tk (Signature unitaire) ---
ax3 = fig.add_subplot(223, projection='3d')
surf3 = ax3.plot_trisurf(maillage, df_Tk['T'], cmap='viridis', edgecolor='none')
ax3.set_title('Figure 11.6 : Signature d\'une résistance (Tk)')
ax3.set_zlabel('T (Unit)')
fig.colorbar(surf3, ax=ax3, shrink=0.5)

# --- FIG 11.7 : T_optimisée (Plateau 250°C) ---
ax4 = fig.add_subplot(224, projection='3d')
surf4 = ax4.plot_trisurf(maillage, df_Topt['T'], cmap='jet', edgecolor='none')
ax4.set_title('Figure 11.7 : Température Optimisée (Cible 250°C)')
ax4.set_zlabel('T (°C)')
ax4.set_zlim(0, 350) # Pour bien voir le plateau
fig.colorbar(surf4, ax=ax4, shrink=0.5)

# Ajustement de la vue pour tous les subplots
for ax in [ax1, ax2, ax3, ax4]:
    ax.view_init(elev=25, azim=-45)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])

# 5. Sauvegarde et Affichage
chemin_save = os.path.join(dossier_img, 'Resultats_Complets_Four.png')
plt.savefig(chemin_save, dpi=300)
print(f"-> Rapport visuel complet sauvegardé dans '{chemin_save}'")

print("-> Ouverture de la fenêtre interactive...")
plt.show()