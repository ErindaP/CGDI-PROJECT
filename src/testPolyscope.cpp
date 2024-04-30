#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/utilities/vector3.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"


#include<ctime>

#include<new>

using namespace geometrycentral;
using namespace geometrycentral::surface;


// Create a structure for a cube
struct Cube
{
  int minIndex;
  int maxIndex;
};


struct Sphere
{
  Vector3 center;
  int radius;
};





// Create a structure to store an object, its mesh, its geometry, and its forces
struct Object
{
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  std::vector<float> forces;
  std::string name; // Name of the object
  int type; // Type of the object (0: cube, 1: sphere, 2: cylinder)
  Cube cube; // Cube structure
  Sphere sphere; // Sphere structure


};


// Create a vector to store the pointers to the objects
std::vector<Object*> objects;




// Create a manifold surface mesh
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// Define the force vector
std::vector<float> forces = {0., 0., 0.};


// Define the number of objects
int numObjects = 0;
bool isForceApplied = false;
bool isAllForcesApplied = false;


// Function to apply a force on a surface mesh
void applyForce(ManifoldSurfaceMesh& mesh, std::vector<float> forces)
{
  // Get the vertex positions
  VertexData<Vector3> vertexPositions = geometry->inputVertexPositions;

  // Apply the force to each vertex
  for (Vertex v : mesh.vertices())
  {
    Vector3 force = {forces[0], forces[1], forces[2]};
    vertexPositions[v] += force;
  }

  // Update the vertex positions
  geometry->vertexPositions = vertexPositions;
  

}

void applyForceObject(Object& object)
{
  // Get the vertex positions
  std::cout <<" Debug1" << std::endl;
  
  VertexData<Vector3> vertexPositions = object.geometry->inputVertexPositions;
  std::cout <<" Debug" << std::endl;
  // Apply the force to each vertex
  Vector3 force = {object.forces[0], object.forces[1], object.forces[2]};

  for (Vertex v : object.mesh->vertices())
  {
    vertexPositions[v] += force;
  }

  // Update the vertex positions
  object.geometry->vertexPositions = vertexPositions;
  if (object.type == 1)
  {

    object.sphere.center += force;
  }

}

// Function that apply each object's force to its mesh
void applyAllForces()
{
  std::cout << "Applying all forces" << std::endl;
  std::cout << objects[0] << std::endl;
  std::cout << objects[0]->name << std::endl;
  
  for (int i = 0; i < objects.size(); i++)
  {
    std::cout << "Applying force to object " << (objects[i])->name << std::endl;
    applyForceObject(*objects[i]);
    polyscope::getSurfaceMesh(objects[i]->name)->updateVertexPositions(objects[i]->geometry->inputVertexPositions);

  }
  
}


// Function to generate a random sphere
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> generateRandomSphere(Vector3 center, int radius)
{
  // Define the vertex positions
  std::vector<Vector3> vertexPositions;
  std::vector<std::vector<size_t>> faceIndices;

  // Create the sphere
  for (int i = 0; i <= 50; i++) {
    for (int j = 0; j <= 50; j++) {
      double theta = i * (2 * M_PI) / 50;
      double phi = j * M_PI / 50;
      double x = center[0] + radius * sin(phi) * cos(theta);
      double y = center[1] + radius * sin(phi) * sin(theta);
      double z = center[2] + radius * cos(phi);
      vertexPositions.push_back({x, y, z});
    }
  }

  for (int i = 0; i < 50; i++) {
    for (int j = 0; j < 50; j++) {
      unsigned long int p1 = i * 51 + j;
      unsigned long int p2 = p1 + 51;
      unsigned long int p3 = p1 + 1;
      unsigned long int p4 = p2 + 1;
      faceIndices.push_back({p1, p2, p4, p3});
    }
  }

  return makeManifoldSurfaceMeshAndGeometry(faceIndices, vertexPositions);
}


// Function to generate a random cube



std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> generateRandomCube(Vector3 p1, Vector3 p2){


  // Define the vertex positions

  std::vector<Vector3> vertexPositions = {
      {p1[0], p1[1], p1[2]}, {p2[0], p1[1], p1[2]}, {p2[0], p2[1], p1[2]}, {p1[0], p2[1], p1[2]},
      {p1[0], p1[1], p2[2]}, {p2[0], p1[1], p2[2]}, {p2[0], p2[1], p2[2]}, {p1[0], p2[1], p2[2]}
  };
  std::vector<std::vector<size_t>> faceIndices = {{0, 1, 2, 3}, {7, 6, 5, 4}, {0, 4, 5, 1},
                                                  {1, 5, 6, 2}, {2, 6, 7, 3}, {3, 7, 4, 0}};


  // Register the vertices as a point cloud 
  polyscope::registerPointCloud("cubePoints", vertexPositions);
  // Add a value to each point corresponding to the vertex index
  std::vector<double> pointValues(vertexPositions.size());
  for (size_t i = 0; i < vertexPositions.size(); i++) {
    pointValues[i] = i;
  }
  polyscope::getPointCloud("cubePoints")->addScalarQuantity("vertex index", pointValues);



  return makeManifoldSurfaceMeshAndGeometry(faceIndices, vertexPositions);


}



// Function to generate random objects
void generateRandomObjects(int numObjects)
{
  // Generate random objects
  for (int i = 0; i < numObjects; i++)
  {
    // Create a new object
    Object* object = new Object;

    int choice = rand() % 2 ;
    std::cout << "Choice: " << choice << std::endl;
    if (choice == 0)
      {
        // Generate a random sphere
        // Select a random point and generate a sphere around it
        Vector3 center = {rand() % 10, rand() % 10, rand() % 10};
        int radius = (rand() % 2)+1;
        Sphere* sphere = new Sphere;
        sphere->center = center;
        sphere->radius = radius;
        object->type = 1;
        object->name = "sphere" + std::to_string(i);
        object->sphere = *sphere;
        std::tie(object->mesh, object->geometry) = generateRandomSphere(center, radius);
        std::cout << "Sphere center: " << center << " radius: " << radius << std::endl;
      }
    else
      {




      // Generate a random cube
        // Select two random points and generate a cube from them
      Vector3 p1 = {rand() % 10, rand() % 10, rand() % 10};
      //Vector3 p2 = {rand() % 10, rand() % 10, rand() % 10};
      Vector3 p2 = {p1[0] + 0.5+ (rand()%2), p1[1] + 0.5 + (rand()%3), p1[2] + 0.5+(rand()%3)};

      Cube* cube = new Cube;

      object->cube = *cube;
      object->type = 0;



      std::tie(object->mesh, object->geometry) = generateRandomCube(p1, p2);
      object->name = "cube" + std::to_string(i);

      // Find the points that have min x y z and max x y z among the 8 points
      Vector3 min = object->geometry->inputVertexPositions[object->mesh->vertex(0)];
      Vector3 max = object->geometry->inputVertexPositions[object->mesh->vertex(0)];
      int minIndex = 0;
      int maxIndex = 0;
      for (int j = 1 ; j < 8; j++)
      {
        Vector3 p = object->geometry->inputVertexPositions[object->mesh->vertex(j)];
        if (p[0] <= min[0] && p[1] <= min[1] && p[2] <= min[2])
        {
          min = p;
          minIndex = j;
        }
        if (p[0] >= max[0] && p[1] >= max[1] && p[2] >= max[2])
        {
          max = p;
          maxIndex = j;
        }
      }
      object->cube.minIndex = minIndex;
      object->cube.maxIndex = maxIndex;
      std::cout << "minIndex: " << minIndex << " maxIndex: " << maxIndex << std::endl;
      }
      






    

    // Register the mesh with polyscope
    polyscope::registerSurfaceMesh(object->name, object->geometry->inputVertexPositions, object->mesh->getFaceVertexList());
    // Generate a random small force vector
    std::vector<float> forces = {(rand() % 10) /  100., (rand() % 10 ) / 100., (rand() % 10 )/ 100.};
    std::cout << object->name << std::endl;
    std::cout << "Force: " << forces[0] << " " << forces[1] << " " << forces[2] << std::endl;
    // Store the object
    object->forces = forces;

    objects.push_back(object);



  }
}


// Function to check collisions between objects (AABB) if there is a collision invert the force vector
void checkCollisions()
{
  // Check collisions between objects
  for (int i = 0; i < objects.size(); i++)
  {


    if (objects[i]->type == 1)
    {
      // Check if the object hit the walls
      Vector3 rad = {objects[i]->sphere.radius, objects[i]->sphere.radius, objects[i]->sphere.radius};
      Vector3 min = objects[i]->sphere.center - rad;
      Vector3 max = objects[i]->sphere.center + rad ;

      if (min[0] < 0 || max[0] > 15 || min[1] < 0 || max[1] > 15 || min[2] < 0 || max[2] > 15)
      {
        /*
        objects[i]->forces[0] = -objects[i]->forces[0];
        objects[i]->forces[1] = -objects[i]->forces[1];
        objects[i]->forces[2] = -objects[i]->forces[2];
        */
       // Make the symetric of the force vector with respect to the wall 
        if (min[0] < 0)
        {
          objects[i]->forces[0] = -objects[i]->forces[0];
        }
        if (max[0] > 15)
        {
          objects[i]->forces[0] = -objects[i]->forces[0];
        }
        if (min[1] < 0)
        {
          objects[i]->forces[1] = -objects[i]->forces[1];
        }
        if (max[1] > 15)
        {
          objects[i]->forces[1] = -objects[i]->forces[1];
        }
        if (min[2] < 0)
        {
          objects[i]->forces[2] = -objects[i]->forces[2];
        }
        if (max[2] > 15)
        {
          objects[i]->forces[2] = -objects[i]->forces[2];
        }
      }
    }
    else
    {
      // Check if the object hit the walls
      Vector3 min = objects[i]->geometry->inputVertexPositions[objects[i]->mesh->vertex(objects[i]->cube.minIndex)];
      Vector3 max = objects[i]->geometry->inputVertexPositions[objects[i]->mesh->vertex(objects[i]->cube.maxIndex)];

      if (min[0] < 0 || max[0] > 15 || min[1] < 0 || max[1] > 15 || min[2] < 0 || max[2] > 15)
      {
        /*
        objects[i]->forces[0] = -objects[i]->forces[0];
        objects[i]->forces[1] = -objects[i]->forces[1];
        objects[i]->forces[2] = -objects[i]->forces[2];
        */
        // Make the symetric of the force vector with respect to the wall
        if (min[0] < 0)
        {
          objects[i]->forces[0] = -objects[i]->forces[0];
        }
        if (max[0] > 15)
        {
          objects[i]->forces[0] = -objects[i]->forces[0];
        }
        if (min[1] < 0)
        {
          objects[i]->forces[1] = -objects[i]->forces[1];
        }
        if (max[1] > 15)
        {
          objects[i]->forces[1] = -objects[i]->forces[1];
        }
        if (min[2] < 0)
        {
          objects[i]->forces[2] = -objects[i]->forces[2];
        }
        if (max[2] > 15)
        {
          objects[i]->forces[2] = -objects[i]->forces[2];
        }


      }
    }

    /*
    // if object hit the walls invert the force vector
    Vector3 min = objects[i]->geometry->inputVertexPositions[objects[i]->mesh->vertex(0)];
    Vector3 max = objects[i]->geometry->inputVertexPositions[objects[i]->mesh->vertex(6)];

    if (min[0] < 0 || max[0] > 15 || min[1] < 0 || max[1] > 15 || min[2] < 0 || max[2] > 15)
    {
      objects[i]->forces[0] = -objects[i]->forces[0];
      objects[i]->forces[1] = -objects[i]->forces[1];
      objects[i]->forces[2] = -objects[i]->forces[2];
    }

    */


    for (int j = i + 1; j < objects.size(); j++)
    {
      // Check if there is a collision
      // If there is a collision invert the force vector
      // Bounding box test : Vertex 0 and Vertex 6 are the min and max of the bounding box

      if (objects[i]->type == 1 && objects[j]->type == 1)
      {
        // Check if the two spheres collide
        Vector3 center1 = objects[i]->sphere.center;
        Vector3 center2 = objects[j]->sphere.center;
        int radius1 = objects[i]->sphere.radius;
        int radius2 = objects[j]->sphere.radius;

        if (sqrt(pow(center1[0] - center2[0], 2) + pow(center1[1] - center2[1], 2) + pow(center1[2] - center2[2], 2)) < radius1 + radius2)
        {
          // Invert the force vector
          objects[i]->forces[0] = -objects[i]->forces[0];
          objects[i]->forces[1] = -objects[i]->forces[1];
          objects[i]->forces[2] = -objects[i]->forces[2];

          objects[j]->forces[0] = -objects[j]->forces[0];
          objects[j]->forces[1] = -objects[j]->forces[1];
          objects[j]->forces[2] = -objects[j]->forces[2];
        }
      }
      else if (objects[i]->type == 0 && objects[j]->type == 0)
      {
        // Check if the two cubes collide

        Vector3 min1 = objects[i]->geometry->inputVertexPositions[objects[i]->mesh->vertex(objects[i]->cube.minIndex)];
        Vector3 max1 = objects[i]->geometry->inputVertexPositions[objects[i]->mesh->vertex(objects[i]->cube.maxIndex)];

        Vector3 min2 = objects[j]->geometry->inputVertexPositions[objects[j]->mesh->vertex(objects[j]->cube.minIndex)];
        Vector3 max2 = objects[j]->geometry->inputVertexPositions[objects[j]->mesh->vertex(objects[j]->cube.maxIndex)];

        std::cout << "min1: " << min1 << " max1: " << max1 << std::endl;
        std::cout << "min2: " << min2 << " max2: " << max2 << std::endl;

        if (min1[0] < max2[0] && max1[0] > min2[0] &&
            min1[1] < max2[1] && max1[1] > min2[1] &&
            min1[2] < max2[2] && max1[2] > min2[2])
        {
          // Invert the force vector
          objects[i]->forces[0] = -objects[i]->forces[0];
          objects[i]->forces[1] = -objects[i]->forces[1];
          objects[i]->forces[2] = -objects[i]->forces[2];

          objects[j]->forces[0] = -objects[j]->forces[0];
          objects[j]->forces[1] = -objects[j]->forces[1];
          objects[j]->forces[2] = -objects[j]->forces[2];
        }
      }
      else
      {
        // Check if the cube and the sphere collide
        if (objects[i]->type == 0)
        {
          // Check if the cube and the sphere collide
          Vector3 min = objects[i]->geometry->inputVertexPositions[objects[i]->mesh->vertex(objects[i]->cube.minIndex)];
          Vector3 max = objects[i]->geometry->inputVertexPositions[objects[i]->mesh->vertex(objects[i]->cube.maxIndex)];

          Vector3 center = objects[j]->sphere.center;
          int radius = objects[j]->sphere.radius;

          if (center[0] + radius > min[0] && center[0] - radius < max[0] &&
              center[1] + radius > min[1] && center[1] - radius < max[1] &&
              center[2] + radius > min[2] && center[2] - radius < max[2])
          {
            // Invert the force vector
            objects[i]->forces[0] = -objects[i]->forces[0];
            objects[i]->forces[1] = -objects[i]->forces[1];
            objects[i]->forces[2] = -objects[i]->forces[2];

            objects[j]->forces[0] = -objects[j]->forces[0];
            objects[j]->forces[1] = -objects[j]->forces[1];
            objects[j]->forces[2] = -objects[j]->forces[2];
          }
        }
        else
        {
          // Check if the cube and the sphere collide
          Vector3 min = objects[j]->geometry->inputVertexPositions[objects[j]->mesh->vertex(objects[j]->cube.minIndex)];
          Vector3 max = objects[j]->geometry->inputVertexPositions[objects[j]->mesh->vertex(objects[j]->cube.maxIndex)];

          Vector3 center = objects[i]->sphere.center;
          int radius = objects[i]->sphere.radius;

          if (center[0] + radius > min[0] && center[0] - radius < max[0] &&
              center[1] + radius > min[1] && center[1] - radius < max[1] &&
              center[2] + radius > min[2] && center[2] - radius < max[2])
          {
            // Invert the force vector
            objects[i]->forces[0] = -objects[i]->forces[0];
            objects[i]->forces[1] = -objects[i]->forces[1];
            objects[i]->forces[2] = -objects[i]->forces[2];

            objects[j]->forces[0] = -objects[j]->forces[0];
            objects[j]->forces[1] = -objects[j]->forces[1];
            objects[j]->forces[2] = -objects[j]->forces[2];




    }
  }
      }
    }
  }
}




// Callback function to be called when the user changes the force vector
void myCallback()
{
  ImGui::PushItemWidth(100);
  // Field to store the force vector
  ImGui::InputFloat("x", &forces[0]);
  ImGui::InputFloat("y", &forces[1]);
  ImGui::InputFloat("z", &forces[2]);



  if (ImGui::Button("Apply force")) {
    // Apply the force
    //applyForce(*mesh, forces);
    isForceApplied = !isForceApplied;

  }

  if (isForceApplied) {
    applyForce(*mesh, forces);
    polyscope::getSurfaceMesh("cube")->updateVertexPositions(geometry->inputVertexPositions);

}
  ImGui::PopItemWidth();

  ImGui::InputInt("Number of objects", &numObjects);


  if (ImGui::Button("Generate random objects")){

    generateRandomObjects(numObjects);



  }

  if (ImGui::Button("Apply all forces")) {
    isAllForcesApplied = !isAllForcesApplied;
  }

  if (isAllForcesApplied) {
    applyAllForces();
    checkCollisions();
  }
    


  
}




int main(int argc, char **argv) {
  // Seed the random number generator
  srand(time(NULL));


  // Initialize polyscope
  polyscope::init();


  // Make a cube mesh
  // Define the vertex positions
  /*
  std::vector<Vector3> vertexPositions = {
      {0., 0., 0.}, {1., 0., 0.}, {1., 1., 0.}, {0., 1., 0.},
      {0., 0., 1.}, {1., 0., 1.}, {1., 1., 1.}, {0., 1., 1.}};
  std::vector<std::vector<size_t>> faceIndices = {{0, 1, 2, 3}, {7, 6, 5, 4}, {0, 4, 5, 1},
                                                  {1, 5, 6, 2}, {2, 6, 7, 3}, {3, 7, 4, 0}};

  auto cubeMesh = polyscope::registerSurfaceMesh("cube", vertexPositions, faceIndices);
  cubeMesh->setSurfaceColor({0.1, 0.1, 0.8});
  */ 


  // Define a manifold surface mesh (a cube)



  /*

  // Define the vertex positions
  std::vector<Vector3> vertexPositions = {
      {0., 0., 0.}, {1., 0., 0.}, {1., 1., 0.}, {0., 1., 0.},
      {0., 0., 1.}, {1., 0., 1.}, {1., 1., 1.}, {0., 1., 1.}};
  std::vector<std::vector<size_t>> faceIndices = {{0, 1, 2, 3}, {7, 6, 5, 4}, {0, 4, 5, 1},
                                                  {1, 5, 6, 2}, {2, 6, 7, 3}, {3, 7, 4, 0}};
  std::tie(mesh, geometry) = makeManifoldSurfaceMeshAndGeometry(faceIndices, vertexPositions);

  // Register the mesh with polyscope
  polyscope::registerSurfaceMesh("cube", geometry->inputVertexPositions, mesh->getFaceVertexList());

  */










  // Specify the callback function
  polyscope::state::userCallback = myCallback;

  // Give control to the polyscope gui
  polyscope::show();


// J'ai pu traiter jusqu'à la question 7  que je n'ai pas vraiment fini




  return EXIT_SUCCESS;
}
