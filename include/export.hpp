#pragma once
#include <Eigen/Dense>
#include <string>
#include "maillage.hpp"

// Exporte un champ de température dans un fichier CSV
void exportTemperatureField(const std::string& filename,
                            const Maillage& mesh,
                            const Eigen::VectorXd& temperature);

// Exporte la géométrie du maillage
void exportMeshGeometry(const std::string& filename, const Maillage& mesh);