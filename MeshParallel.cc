#include "common.h"

// minimalistic representation of triangle mesh, supports parallel Laplacian smoothing
void MeshParallel()
{
  int nVertices = 1000000;

  std::vector<int> adjVec;
  std::vector<int> firstAdj;
  std::vector<int> numAdj;

  std::vector<MyMesh::Point> points;

  printf("setup\n");
  for(int iVertex = 0; iVertex < nVertices; iVertex++) {

    points.push_back({1.0*rand(), 1.0*rand(), 1.0*rand()});

    firstAdj.push_back(adjVec.size());

    int nAdj = 6;

    numAdj.push_back(6);

    for(int iAdj = 0; iAdj < nAdj; iAdj++)
      adjVec.push_back(rand() % nVertices);
  }

  std::vector<MyMesh::Point> nextPoints;
  nextPoints.resize(points.size());

  for(int iIter = 0; iIter < 1000; iIter++) {

    if(!(iIter % 10))
      printf("%d.\n", iIter);

#pragma omp parallel for
    for(size_t iV = 0; iV < numAdj.size(); iV++) {

      nextPoints[iV] = {0.0, 0.0, 0.0};

      for(int iN = 0; iN < numAdj[iV]; iN++)
        nextPoints[iV] += points[adjVec[firstAdj[iV] + iN]];

      nextPoints[iV] *= 1.0/numAdj[iV];
    }

    std::swap(points, nextPoints);
  }

  printf("done\n");
}
