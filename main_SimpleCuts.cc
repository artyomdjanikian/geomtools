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
// Import mesh of arbitrary topology, compute simple cuts and remesh it with quads
// making edge lines follow the simple cuts and sharp edges

int main(int argc, char **argv)
{
  int nIter = argc > 2? atoi(argv[2]) : 100;
  int nCycles = argc > 3? atoi(argv[3]) : 3;
  double step = argc > 4? atof(argv[4]) : 5.0;
  double dT = argc > 5? atof(argv[5]) : 0.005;
  double ff = argc > 6? atof(argv[6]) : 1.0;


  MyMesh mesh;
  mesh.request_face_normals();
  mesh.request_vertex_normals();
  mesh.request_vertex_colors();
  // IO::Options ropt;
 
  // -------------------- read mesh
 
  // TODO : extract into OpenMeshMP

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

  //  SDFTest(nIter, mesh, tree);

  // OctreeTest(nIter, mesh, tree);

  //  __SimVolumeTest(nIter, nCycles, bounds, step, dT, ff);

  //  SimulateCloth(mesh);

  // PatchLoopTest(mesh, 0.125);

  // Bottle_Opener_remeshplys.ply
  QuadRemeshTest(mesh, step, nIter);

  // SimulateSandCastle(mesh, step, nIter);



  return 0;
}
