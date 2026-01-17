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

//./build/Build/bin/PatchLoops ~/Documents/Models/topology-optimized-bottle-opener-1.snapshot.2/Bottle_Opener_remesh_hole_pair.ply 0.07
//./build/Build/bin/PatchLoops ~/Documents/Models/topology-optimized-bottle-opener-1.snapshot.2/Bottle_Opener_remesh_hole_pair.ply 0.125

int main(int argc, char **argv)
{
  MyMesh mesh;
  mesh.request_face_normals();
  mesh.request_vertex_normals();
  mesh.request_vertex_colors();
  // IO::Options ropt;
 
  // -------------------- read mesh
 
  // TODO : extract into OpenMeshMP

  double step = argc > 2? atof(argv[2]) : 0.125;
  int nIter = argc > 3? atoi(argv[3]) : 100;

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

  PatchLoops(mesh, step);

  return 0;
}
