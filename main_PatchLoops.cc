#include <deque>
#include <cassert>
#include <csignal>
#include <unistd.h>

#include "common.h"

#include "happly.h"
#include "AABBTree.h"
#include "cOctree.h"
#include "cSimVolume.h"
#include "MeshSampling.h"
#include "fibonacci.hpp" 
#include "cGrid3d.h"
#include "SpatialHashTable.h"

// ----------------------------------------------------------------------------
// Import mesh with holes, patch them with facets
// ----------------------------------------------------------------------------

int main(int argc, char **argv)
{
  MyMesh mesh;
  mesh.request_face_normals();
  mesh.request_vertex_normals();
  mesh.request_vertex_colors();
  // IO::Options ropt;
 
  // -------------------- read mesh
 
  // TODO : extract into OpenMeshMP

<<<<<<< HEAD
  double step = argc > 2? atof(argv[2]) : 0.125;
  int nIter = argc > 3? atoi(argv[3]) : 100;

=======
>>>>>>> 2c02bf3f30c7ffab9d1ea0db43f7aa75558a0a73
  printf("read mesh %s\n", argv[1]);

  OpenMesh::IO::Options ropt;
  ropt += OpenMesh::IO::Options::Binary; 
  ropt += OpenMesh::IO::Options::VertexColor;

  if ( ! OpenMesh::IO::read_mesh(mesh, argv[1], ropt))
  {
    std::cerr << "Error loading mesh from file "<< std::endl;
    return 1;
  }
 
  mesh.update_normals();

<<<<<<< HEAD
  PatchLoops(mesh, step);
=======
  PatchLoops(mesh, 0.125);
>>>>>>> 2c02bf3f30c7ffab9d1ea0db43f7aa75558a0a73

  return 0;
}
