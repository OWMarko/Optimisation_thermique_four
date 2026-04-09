#include "../include/export.hpp"
#include <fstream>

void exportTemperatureField(const std::string& filename,
                            const Maillage& mesh,
                            const Eigen::VectorXd& temperature){
    std::ofstream file(filename);
    file << "x,y,T\n";
    for (int i = 0; i < mesh.Nbpt; ++i)
    {
        file << mesh.Coorneu(i, 0) << "," << mesh.Coorneu(i, 1) << "," 
             << temperature(i) << "\n";
    }
    file.close();
}

void exportMeshGeometry(const std::string& filename,
                        const Maillage& mesh) {
    std::ofstream file(filename);
    file << "n1,n2,n3\n";
    for (int i = 0; i < mesh.Nbtri; ++i)
    {
        file << mesh.Numtri(i, 0) << "," << mesh.Numtri(i, 1) << "," 
             << mesh.Numtri(i, 2) << "\n";
    }
    file.close();
}