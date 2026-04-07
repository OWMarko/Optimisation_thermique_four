import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.patches as patches
import os

# 1. Lecture des données exportées par le C++
print("-> Lecture du maillage...")
df_T0 = pd.read_csv('../data/noeuds_T0.csv')
triangles = pd.read_csv('../data/triangles.csv')

print(f"   Nœuds: {len(df_T0)}, Triangles: {len(triangles)}")

# 2. Création de la triangulation
maillage = tri.Triangulation(df_T0['x'], df_T0['y'], triangles.values)

# 3. Création de la figure 2D
fig, ax = plt.subplots(figsize=(12, 12))

# === Affichage du maillage ===
# Tracer avec des lignes fines pour voir la structure
ax.triplot(maillage, linewidth=0.5, color='blue', alpha=0.6)

# === Affichage du contour du four ===
four_rect = patches.Rectangle((-1, -1), 2, 2, linewidth=3, 
                               edgecolor='black', facecolor='none', 
                               label='Four ([-1,1]²)')
ax.add_patch(four_rect)

# === Affichage du rectangle de résine ===
resine_rect = patches.Rectangle((-0.5, -0.2), 1, 0.4, linewidth=2.5, 
                                edgecolor='red', facecolor='red', 
                                alpha=0.15, label='Résine ([-0.5,0.5]×[-0.2,0.2])')
ax.add_patch(resine_rect)

# === Affichage des nœuds ===
ax.scatter(df_T0['x'], df_T0['y'], s=1, c='darkblue', alpha=0.3)

# === Configuration de la figure ===
ax.set_aspect('equal')
ax.set_xlim(-1.2, 1.2)
ax.set_ylim(-1.2, 1.2)
ax.set_xlabel('x', fontsize=12)
ax.set_ylabel('y', fontsize=12)
ax.set_title('Visualisation du Maillage - Four Thermique', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.2, linestyle='--')
ax.legend(loc='upper right', fontsize=11)

# === Ajout d'informations ===
info_text = f"Maillage: {len(df_T0)} nœuds, {len(triangles)} triangles (P1)"
ax.text(0.02, 0.98, info_text, transform=ax.transAxes, 
        fontsize=10, verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

plt.tight_layout()

# 4. Sauvegarde
dossier_img = '../img'
if not os.path.exists(dossier_img):
    os.makedirs(dossier_img)

chemin_save = os.path.join(dossier_img, 'Maillage_Four.png')
plt.savefig(chemin_save, dpi=300, bbox_inches='tight')
print(f"-> Maillage sauvegardé dans '{chemin_save}'")

print("-> Ouverture de la fenêtre interactive...")
plt.show()
