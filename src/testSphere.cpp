#include <iostream>
#include <vector>
#include "polyscope/polyscope.h"
#include "geometrycentral/surface/geometry.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/utilities/vector3.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"


using namespace geometrycentral;
using namespace geometrycentral::surface;

int main() {
    // Initialiser polyscope
    polyscope::init();

    // Créer une nouvelle géométrie de position de sommet pour stocker les données de la sphère
    std::unique_ptr<VertexPositionGeometry> geo;
    std::unique_ptr<ManifoldSurfaceMesh> mesh; 
    // Créer les données pour la sphère
    int numDivisions = 50;
    std::vector<Vector3> positions;
    std::vector<std::vector<size_t>> indices;

    for (int i = 0; i <= numDivisions; i++) {
        for (int j = 0; j <= numDivisions; j++) {
            double theta = i * (2 * M_PI) / numDivisions;
            double phi = j * M_PI / numDivisions;
            double x = sin(phi) * cos(theta);
            double y = sin(phi) * sin(theta);
            double z = cos(phi);
            positions.push_back({x, y, z});
        }
    }

    /*
    // Ajouter les sommets à la géométrie
    for (int i = 0; i < numDivisions; i++) {
        for (int j = 0; j < numDivisions; j++) {
            unsigned long int p1 = i * (numDivisions + 1) + j;
            unsigned long int p2 = p1 + numDivisions + 1;
            unsigned long int p3 = p1 + 1;
            unsigned long int p4 = p2 + 1;
            indices.push_back({p1, p2, p4});
            indices.push_back({p1, p4, p3});
        }
    }

    std::tie(mesh, geo) = makeManifoldSurfaceMeshAndGeometry(positions, indices);


    // Ajouter la géométrie à polyscope
    polyscope::registerSurfaceMesh("sphere", geo->inputVertexPositions, mesh->getFaceVertexList());
*/
      // Define the vertex positions
    std::unique_ptr<VertexPositionGeometry> geometry;

  /*std::vector<Vector3> vertexPositions = {
      {0., 0., 0.}, {1., 0., 0.}, {1., 1., 0.}, {0., 1., 0.},
      {0., 0., 1.}, {1., 0., 1.}, {1., 1., 1.}, {0., 1., 1.}};
  */
  
  std::vector<std::vector<size_t>> faceIndices ;
  for (int i = 0; i < numDivisions; i++) {
        for (int j = 0; j < numDivisions; j++) {
            unsigned long int p1 = i * (numDivisions + 1) + j;
            unsigned long int p2 = p1 + numDivisions + 1;
            unsigned long int p3 = p1 + 1;
            unsigned long int p4 = p2 + 1;
            faceIndices.push_back({p1, p2, p4, p3});
        }
    }

  std::tie(mesh, geometry) = makeManifoldSurfaceMeshAndGeometry(faceIndices, positions);
  polyscope::registerSurfaceMesh("cube", geometry->inputVertexPositions, mesh->getFaceVertexList());
  


    // Visualiser la sphère
    polyscope::show();

    // Create another sphere with another center and a different radius
    std::vector<Vector3> positions2;
    std::vector<std::vector<size_t>> indices2;
    float radius = 0.5;
    for (int i = 0; i <= numDivisions; i++) {
        for (int j = 0; j <= numDivisions; j++) {
            double theta = i * (2 * M_PI) / numDivisions;
            double phi = j * M_PI / numDivisions;
            double x = radius * sin(phi) * cos(theta) + 1;
            double y = radius * sin(phi) * sin(theta) + 1;
            double z = radius * cos(phi) + 1;
            positions2.push_back({x, y, z});
        }
    }

    return 0;
}
