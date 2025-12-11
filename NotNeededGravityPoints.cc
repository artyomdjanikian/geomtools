  auto GravitySimulation = [](MyMesh &mesh, int nIter) 
  {
    // gravity attraction of mass points constrained to the mesh surface

    AABBTree tree;
    tree.Build(mesh, 20);
  
    // gravity simulation
    std::vector<OpenMesh::Vec3f> normals;
    std::vector<OpenMesh::Vec3f> points;

    normals.resize(mesh.n_vertices());
    points.resize(mesh.n_vertices());

    for(auto vertexHandle : mesh.vertices()) {
        auto pos = mesh.point(vertexHandle);
        auto normal = mesh.normal(vertexHandle);

        points[vertexHandle.idx()] = pos;
        normals[vertexHandle.idx()] = normal;
    }

    for(int iTry = 0; iTry < nIter; iTry++) {
      std::vector<double> tangs;
      std::vector<OpenMesh::Vec3f> forces;

      tangs.resize(mesh.n_vertices());
      forces.resize(mesh.n_vertices());

      for(auto vertexHandle : mesh.vertices()) {

        auto pos = points[vertexHandle.idx()];
        auto normal = normals[vertexHandle.idx()];
        OpenMesh::Vec3f force(0.0, 0.0, 0.0);
        for(auto gravityHandle : mesh.vertices()) {
          if(vertexHandle != gravityHandle) {

            auto gPos = points[gravityHandle.idx()];

            OpenMesh::Vec3f radVec = gPos-pos;
            double dist = radVec.norm();

            if(dist > 1e-5) {

              force += radVec * 1.0/(dist*dist*dist);

            }
          }
        }

        force -= normal*normal.dot(force);
        forces[vertexHandle.idx()] = force;

        tangs[vertexHandle.idx()] = force.norm();
      }

      double maxTang = 0.0;
      double averTang = 0.0;
  
      for(auto tang : tangs) {
        maxTang = std::max(tang, maxTang);
        averTang += tang;
      }

      averTang /= tangs.size();

      printf("%.6f maxTang %.6f averTang\n", maxTang, averTang);

      double factor = averTang/0.1;

      happly::PLYExport eply;

      double averDist = 0.0;

      for(auto vertexHandle : mesh.vertices()) {
        int ii = vertexHandle.idx();

        auto pos = mesh.point(vertexHandle);

        double len = tangs[ii]/averTang;

        double darkFactor = 1.0;

        if (len > 1.0) {

          darkFactor = 1.0/len;

          len = 1.0;
        }

        double red = darkFactor*len;
        double green = darkFactor*4.0*len*(1.0-len);
        double blue = darkFactor*1.0-len;


        pos += forces[ii]/factor;

        auto np = tree.FindNearestPoint(toVec(pos));

        if(np.facetIndex != -1) {
          points[ii] = OpenMesh::Vec3f(np.pnt[0], np.pnt[1], np.pnt[2]);

          normals[ii] = mesh.normal(MyMesh::FaceHandle(np.facetIndex));

          averDist += (mesh.point(vertexHandle)-points[ii]).norm();
        }

        // TODO : project onto the mesh, update points, normals

        eply.AddVertex(points[ii], red, green, blue);
      }

      printf("  aver shift %.6f\n", averDist/mesh.n_vertices());

      eply.WritePLY("gravpoints.ply");
    }
  };

  //GravitySimulation(mesh, 50);