#include "vcInternalModels.h"
#include "vcInternalModelsData.h"

#include "gl/vcMesh.h"
#include "vcPolygonModel.h"

static int gRefCount = 0;
vcMesh *gInternalMeshes[vcInternalMeshType_Count] = {};
vcPolygonModel *gInternalModels[vcInternalModelType_Count] = {};

udResult vcInternalModels_Init()
{
  udResult result;

  ++gRefCount;
  UD_ERROR_IF(gRefCount != 1, udR_Success);

  UD_ERROR_CHECK(vcMesh_Create(&gInternalMeshes[vcInternalMeshType_ScreenQuad], vcP3UV2VertexLayout, int(udLengthOf(vcP3UV2VertexLayout)), screenQuadVertices, 4, screenQuadIndices, 6));
  UD_ERROR_CHECK(vcMesh_Create(&gInternalMeshes[vcInternalMeshType_FlippedScreenQuad], vcP3UV2VertexLayout, int(udLengthOf(vcP3UV2VertexLayout)), flippedScreenQuadVertices, 4, flippedScreenQuadIndices, 6));
  UD_ERROR_CHECK(vcMesh_Create(&gInternalMeshes[vcInternalMeshType_ImGuiQuad], vcImGuiVertexLayout, int(udLengthOf(vcImGuiVertexLayout)), imGuiQuadVertices, 4, imGuiQuadIndices, 6));

  UD_ERROR_CHECK(vcMesh_Create(&gInternalMeshes[vcInternalMeshType_WorldQuad], vcP3N3UV2VertexLayout, int(udLengthOf(vcP3N3UV2VertexLayout)), worldQuadVertices, 4, worldQuadIndices, 6, vcMF_IndexShort));
  UD_ERROR_CHECK(vcMesh_Create(&gInternalMeshes[vcInternalMeshType_Billboard], vcP3UV2VertexLayout, (int)udLengthOf(vcP3UV2VertexLayout), billboardVertices, (uint32_t)udLengthOf(billboardVertices), billboardIndices, (uint32_t)udLengthOf(billboardIndices)));
  UD_ERROR_CHECK(vcMesh_Create(&gInternalMeshes[vcInternalMeshType_Sphere], vcP3N3UV2VertexLayout, (int)udLengthOf(vcP3N3UV2VertexLayout), pSphereVertices, (uint32_t)udLengthOf(sphereVerticesFltArray), sphereIndices, (uint32_t)udLengthOf(sphereIndices)));
  UD_ERROR_CHECK(vcMesh_Create(&gInternalMeshes[vcInternalMeshType_Tube], vcP3N3UV2VertexLayout, (int)udLengthOf(vcP3N3UV2VertexLayout), pTubeVertices, (uint32_t)udLengthOf(tubeVerticesFltArray), tubeIndices, (uint32_t)udLengthOf(tubeIndices)));

  UD_ERROR_CHECK(vcMesh_Create(&gInternalMeshes[vcInternalMeshType_Orbit], vcP3N3VertexLayout, (int)udLengthOf(vcP3N3VertexLayout), pOrbitVertices, (uint32_t)udLengthOf(orbitVerticesFltArray), orbitIndices, (uint32_t)udLengthOf(orbitIndices)));
  UD_ERROR_CHECK(vcMesh_Create(&gInternalMeshes[vcInternalMeshType_Compass], vcP3N3VertexLayout, int(udLengthOf(vcP3N3VertexLayout)), pCompassVerts, (uint32_t)udLengthOf(compassVertsFltArray), compassIndices, (uint32_t)udLengthOf(compassIndices)));

  UD_ERROR_CHECK(vcPolygonModel_CreateFromRawVertexData(&gInternalModels[vcInternalModelType_Cube], (void *)cubeVerticesFltArray, (uint32_t)udLengthOf(cubeVerticesFltArray), vcP3N3UV2VertexLayout, (int)(udLengthOf(vcP3N3UV2VertexLayout)), cubeIndices, (uint32_t)udLengthOf(cubeIndices)));
  UD_ERROR_CHECK(vcPolygonModel_CreateFromRawVertexData(&gInternalModels[vcInternalModelType_Sphere], (void *)sphereVerticesFltArray, (uint32_t)udLengthOf(sphereVerticesFltArray), vcP3N3UV2VertexLayout, (int)(udLengthOf(vcP3N3UV2VertexLayout)), sphereIndices, (uint32_t)udLengthOf(sphereIndices)));
  UD_ERROR_CHECK(vcPolygonModel_CreateFromRawVertexData(&gInternalModels[vcInternalModelType_Cylinder], (void *)cylinderVerticesFltArray, (uint32_t)udLengthOf(cylinderVerticesFltArray), vcP3N3UV2VertexLayout, (int)(udLengthOf(vcP3N3UV2VertexLayout)), cylinderIndices, (uint32_t)udLengthOf(cylinderIndices)));
  UD_ERROR_CHECK(vcPolygonModel_CreateFromRawVertexData(&gInternalModels[vcInternalModelType_Tube], (void *)tubeVerticesFltArray, (uint32_t)udLengthOf(tubeVerticesFltArray), vcP3N3UV2VertexLayout, (int)udLengthOf(vcP3N3UV2VertexLayout), tubeIndices, (uint32_t)udLengthOf(tubeIndices)));
  UD_ERROR_CHECK(vcPolygonModel_CreateFromRawVertexData(&gInternalModels[vcInternalModelType_Quad], (void *)worldQuadVertices, 4, vcP3N3UV2VertexLayout, (int)udLengthOf(vcP3N3UV2VertexLayout), worldQuadIndices, 6));

  result = udR_Success;
epilogue:
  if (result != udR_Success)
    vcInternalModels_Deinit();

  return result;
}

udResult vcInternalModels_Deinit()
{
  udResult result;
  --gRefCount;
  UD_ERROR_IF(gRefCount != 0, udR_Success);

  for (int i = 0; i < vcInternalMeshType_Count; ++i)
    vcMesh_Destroy(&gInternalMeshes[i]);

  for (int i = 0; i < vcInternalModelType_Count; ++i)
    vcPolygonModel_Destroy(&gInternalModels[i]);

  result = udR_Success;
epilogue:
  return result;
}
