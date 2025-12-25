include (VCICommon)

find_package(GTest CONFIG REQUIRED)
enable_testing()

# TODO : add_compile_options(-Wall -Wextra -pedantic)

vci_add_executable (SimpleCuts main_SimpleCuts.cc cSimVolume.cc Dijxtra.cc AABBTree.cc MeshParallel.cc MeshFloodFill.cc SimulateCloth.cc MeshTools.cc cOctree.cc MeshSampling.cc MeshCrossField.cc VolumePoints.cc MeshTopoCuts.cc SDF.cc  QuadRemesh.cc tools3d.cc raytri.cc tribox3.cc tritri_isectline.cc)
vci_add_executable (PatchLoops main_PatchLoops.cc MeshPatchLoops.cc Dijxtra.cc PointSet.cc AABBTree.cc MeshParallel.cc MeshFloodFill.cc SimulateCloth.cc MeshTools.cc cOctree.cc MeshSampling.cc MeshCrossField.cc VolumePoints.cc MeshTopoCuts.cc SDF.cc tools3d.cc raytri.cc tribox3.cc tritri_isectline.cc)
vci_add_executable (Tests test_main.cc tools3d.cc cSimVolume.cc Dijxtra.cc AABBTree.cc MeshParallel.cc MeshFloodFill.cc SimulateCloth.cc MeshTools.cc cOctree.cc MeshSampling.cc MeshCrossField.cc VolumePoints.cc MeshTopoCuts.cc SDF.cc  QuadRemesh.cc tools3d.cc raytri.cc tribox3.cc tritri_isectline.cc)

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
     target_link_libraries(SimpleCuts OpenMeshCore OpenMeshTools OpenMP::OpenMP_CXX)
     target_link_libraries(PatchLoops OpenMeshCore OpenMeshTools OpenMP::OpenMP_CXX)

else()
     target_link_libraries(SimpleCuts OpenMeshCore OpenMeshTools)
     target_link_libraries(PatchLoops OpenMeshCore OpenMeshTools OpenMP::OpenMP_CXX)

endif()

target_link_libraries(Tests PRIVATE gtest OpenMeshCore OpenMeshTools gtest_main pthread)
