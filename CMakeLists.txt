cmake_minimum_required(VERSION 3.10)
project(Four_Thermique)

# 1. Standard C++ et Robustesse
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 2. Drapeaux de compilation (Flags)
# On garde tes -Wall -Wextra pour la propreté du code
# On ajoute -march=native : Crucial en HPC pour que Eigen utilise 
# toutes les capacités SIMD (AVX, FMA) de ton processeur actuel.
if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /O2 /arch:AVX2")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -march=native")
endif()

# 3. Gestion d'Eigen
# Note : Sur un système propre, on préfère find_package, 
# mais je garde ton chemin explicite pour garantir que ça compile chez toi immédiatement.
set(EIGEN_INCLUDE_DIR "/usr/include/eigen3")
include_directories(${EIGEN_INCLUDE_DIR})

# 4. Inclusion du dossier des headers
# Cela permet de faire #include "maillage.hpp" au lieu de #include "../include/maillage.hpp"
include_directories(include)

# 5. Collecte automatique des sources
# Cette commande cherche tous les .cpp dans src/ (maillage, solver, optimisation, export, main)
file(GLOB_RECURSE SOURCES "src/*.cpp")

# 6. Création de l'exécutable
add_executable(${PROJECT_NAME} ${SOURCES})

# Message de log pour vérifier au moment du cmake ..
message(STATUS "Fichiers sources détectés : ${SOURCES}")