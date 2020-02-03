#include "vcTileRenderer.h"
#include "vcQuadTree.h"
#include "vcGIS.h"
#include "vcSettings.h"

#include "gl/vcGLState.h"
#include "gl/vcRenderShaders.h"
#include "gl/vcShader.h"
#include "gl/vcMesh.h"

#include "udThread.h"
#include "udFile.h"
#include "udPlatformUtil.h"
#include "udChunkedArray.h"
#include "udStringUtil.h"

#include "stb_image.h"

#include <vector>

// Debug tiles with colour information
#define VISUALIZE_DEBUG_TILES 0

vcTexture *pDEMTexture[2] = {};

enum
{
  TileVertexResolution = 21, // This should match GPU struct size
  TileIndexResolution = (TileVertexResolution - 1),

  MaxTextureUploadsPerFrame = 3,
};

static const float sTileFadeSpeed = 2.15f;

struct vcTileRenderer
{
  float frameDeltaTime;
  float totalTime;

  vcSettings *pSettings;
  vcQuadTree quadTree;

  vcMesh *pTileMeshes[16];
  vcTexture *pEmptyTileTexture;

  udDouble3 cameraPosition;

  std::vector<std::vector<vcQuadTreeNode*>> *pRenderQueue;
  std::vector<vcQuadTreeNode*> *pTransparentTiles;

  // cache textures
  struct vcTileCache
  {
    volatile bool keepLoading;
    udThread *pThreads[8];
    udSemaphore *pSemaphore;
    udMutex *pMutex;
    udChunkedArray<vcQuadTreeNode*> tileLoadList;
    udChunkedArray<vcQuadTreeNode*> tileTimeoutList;
  } cache;

  struct
  {
    vcShader *pProgram;
    vcShaderConstantBuffer *pConstantBuffer;
    vcShaderSampler *uniform_texture;
    vcShaderSampler *uniform_dem0;
    vcShaderSampler *uniform_dem1;

    struct
    {
      udFloat4x4 projectionMatrix;
      udFloat4x4 viewMatrix;
      udFloat4 eyePositions[4];
      udFloat4 colour;
      udFloat4 demUVs[2 * 4];
      udFloat4 colourUV;
    } everyObject;
  } presentShader;
};

// This functionality here for now until the cache module is implemented
bool vcTileRenderer_TryWriteTile(const char *filename, void *pFileData, size_t fileLen)
{
  udFile *pFile = nullptr;
  if (udFile_Open(&pFile, filename, udFOF_Create | udFOF_Write) == udR_Success)
  {
    udFile_Write(pFile, pFileData, fileLen);
    udFile_Close(&pFile);
    return true;
  }

  return false;
}

// This functionality here for now. In the future will be migrated to udPlatformUtils.
udResult vcTileRenderer_CreateDirRecursive(const char *pFolderPath)
{
  udResult result = udR_Success;
  char *pMutableDirectoryPath = nullptr;

  UD_ERROR_NULL(pFolderPath, udR_InvalidParameter_);

  pMutableDirectoryPath = udStrdup(pFolderPath);
  UD_ERROR_NULL(pMutableDirectoryPath, udR_MemoryAllocationFailure);

  for (uint32_t i = 0;; ++i)
  {
    if (pMutableDirectoryPath[i] == '\0' || pMutableDirectoryPath[i] == '/' || pMutableDirectoryPath[i] == '\\')
    {
      pMutableDirectoryPath[i] = '\0';
      result = udCreateDir(pMutableDirectoryPath);

      // TODO: handle directories already existing
      // TODO: handle path not found
      //if (result != udR_Success)
      //  UD_ERROR_HANDLE();

      pMutableDirectoryPath[i] = pFolderPath[i];

      if (pMutableDirectoryPath[i] == '\0')
        break;
    }
  }

epilogue:
  udFree(pMutableDirectoryPath);

  return result;
}

bool vcTileRenderer_ShouldLoadNode(vcQuadTreeNode *pNode)
{
  return pNode->renderInfo.tryLoad && pNode->touched && (pNode->renderInfo.loadStatus == vcNodeRenderInfo::vcTLS_InQueue);
}

uint32_t vcTileRenderer_LoadThread(void *pThreadData)
{
  vcTileRenderer *pRenderer = (vcTileRenderer*)pThreadData;
  vcTileRenderer::vcTileCache *pCache = &pRenderer->cache;

  while (pCache->keepLoading)
  {
    int loadStatus = udWaitSemaphore(pCache->pSemaphore, 1000);

    if (loadStatus != 0 && pCache->tileLoadList.length == 0)
      continue;

    while (pCache->tileLoadList.length > 0 && pCache->keepLoading)
    {
      udLockMutex(pCache->pMutex);

      // TODO: Store in priority order and recalculate on insert/delete
      int best = -1;
      vcQuadTreeNode *pNode = nullptr;
      udDouble3 tileCenter = udDouble3::zero();
      double bestDistancePrioritySqr = FLT_MAX;

      for (int i = 0; i < (int)pCache->tileLoadList.length; ++i)
      {
        pNode = pCache->tileLoadList[i];

        if (!vcTileRenderer_ShouldLoadNode(pNode))
        {
          pNode->renderInfo.loadStatus = vcNodeRenderInfo::vcTLS_None;
          pCache->tileLoadList.RemoveSwapLast(i);
          --i;
          continue;
        }

        tileCenter = udDouble3::create(pNode->tileCenter, pRenderer->pSettings->maptiles.mapHeight);
        double distanceToCameraSqr = udMagSq3(tileCenter - pRenderer->cameraPosition);

        // root (special case)
        if (pNode == &pRenderer->quadTree.nodes.pPool[pRenderer->quadTree.rootIndex])
        {
          best = i;
          break;
        }

        bool betterNode = true;
        if (best != -1)
        {
          vcQuadTreeNode *pBestNode = pCache->tileLoadList[best];

          // priorities: visibility > failed to render visible area > distance
          betterNode = pNode->visible && !pBestNode->visible;
          if (pNode->visible == pBestNode->visible)
          {
            betterNode = !pNode->rendered && pBestNode->rendered;
            if (pNode->rendered == pBestNode->rendered)
            {
              betterNode = distanceToCameraSqr < bestDistancePrioritySqr;
            }
          }
        }

        if (betterNode)
        {
          bestDistancePrioritySqr = distanceToCameraSqr;
          best = int(i);
        }
      }

      if (best == -1)
      {
        udReleaseMutex(pCache->pMutex);
        break;
      }

      vcQuadTreeNode *pBestNode = pCache->tileLoadList[best];
      pBestNode->renderInfo.loadStatus = vcNodeRenderInfo::vcTLS_Downloading;
      pCache->tileLoadList.RemoveSwapLast(best);

      char localFileName[vcMaxPathLength];
      char serverAddress[vcMaxPathLength];

      udSprintf(localFileName, "%s/%s/%d/%d/%d.%s", pRenderer->pSettings->cacheAssetPath, udUUID_GetAsString(pRenderer->pSettings->maptiles.tileServerAddressUUID), pBestNode->slippyPosition.z, pBestNode->slippyPosition.x, pBestNode->slippyPosition.y, pRenderer->pSettings->maptiles.tileServerExtension);
      udSprintf(serverAddress, "%s/%d/%d/%d.%s", pRenderer->pSettings->maptiles.tileServerAddress, pBestNode->slippyPosition.z, pBestNode->slippyPosition.x, pBestNode->slippyPosition.y, pRenderer->pSettings->maptiles.tileServerExtension);
      udReleaseMutex(pCache->pMutex);

      bool downloadingFromServer = true;
      char *pTileURL = serverAddress;

      if (udFileExists(localFileName) == udR_Success)
      {
        pTileURL = localFileName;
        downloadingFromServer = false;
      }

      udResult result = udR_Failure_;
      void *pFileData = nullptr;
      int64_t fileLen = -1;
      int width = 0;
      int height = 0;
      int channelCount = 0;
      uint8_t *pData = nullptr;

      UD_ERROR_CHECK(udFile_Load(pTileURL, &pFileData, &fileLen));
      UD_ERROR_IF(fileLen == 0, udR_InternalError);

      // Node has been invalidated since download started
      if (!pBestNode->touched)
      {
        // TODO: Put into LRU texture cache (but for now just throw it out)
        pBestNode->renderInfo.loadStatus = vcNodeRenderInfo::vcTLS_None;

        // Even though the node is now invalid - since we the data, it may be put into local disk cache
        UD_ERROR_SET(udR_Success);
      }

      pData = stbi_load_from_memory((stbi_uc*)pFileData, (int)fileLen, (int*)&width, (int*)&height, (int*)&channelCount, 4);
      UD_ERROR_NULL(pData, udR_InternalError);

      //pBestNode->renderInfo.transparency = 0.0f;
      pBestNode->renderInfo.width = width;
      pBestNode->renderInfo.height = height;
      pBestNode->renderInfo.pData = udMemDup(pData, sizeof(uint32_t)*width*height, 0, udAF_None);

      pBestNode->renderInfo.loadStatus = vcNodeRenderInfo::vcTLS_Downloaded;

epilogue:

      if (result != udR_Success)
      {
        pBestNode->renderInfo.loadStatus = vcNodeRenderInfo::vcTLS_Failed;
        if (result == udR_Pending)
        {
          pBestNode->renderInfo.loadStatus = vcNodeRenderInfo::vcTLS_InQueue;

          udLockMutex(pCache->pMutex);
          if (pBestNode->slippyPosition.z <= 10)
          {
            // TODO: server prioritizes these tiles, so will be available much sooner. Requeue immediately
            pCache->tileLoadList.PushBack(pBestNode);
          }
          else
          {
            pBestNode->renderInfo.timeoutTime = pRenderer->totalTime + 15.0f; // 15 seconds
            pCache->tileTimeoutList.PushBack(pBestNode); // timeout it, inserted last
          }
          udReleaseMutex(pCache->pMutex);
        }
      }

      // This functionality here for now until the cache module is implemented
      if (pFileData && fileLen > 0 && downloadingFromServer && pCache->keepLoading)
      {
        if (!vcTileRenderer_TryWriteTile(localFileName, pFileData, fileLen))
        {
          size_t index = 0;
          char localFolderPath[vcMaxPathLength];
          udSprintf(localFolderPath, "%s", localFileName);
          if (udStrrchr(localFileName, "\\/", &index) != nullptr)
            localFolderPath[index] = '\0';

          if (vcTileRenderer_CreateDirRecursive(localFolderPath) == udR_Success)
            vcTileRenderer_TryWriteTile(localFileName, pFileData, fileLen);
        }
      }

      udFree(pFileData);
      stbi_image_free(pData);
    }
  }

  return 0;
}

int up = 1 << 0;
int right = 1 << 1;
int down = 1 << 2;
int left = 1 << 3;

// build meshes
int meshConfigurations[] =
{
  0,

  up,                      //  ^
  right,                   //   >
  down,                    //  .
  left,                    // <

  up | right,              //  ^>
  up | down,               //  ^.
  up | left,               // <^

  right | down,            //  .>
  right | left,            // < >

  down | left,             // <.

  up | right | down,       //  ^.>
  up | left | right,       // <^>
  up | left | down,        // <^.

  down | left | right,     // <.>

  up | left | right | down // <^.>
};


void vcTileRenderer_BuildMeshVertices(vcP3Vertex *pVerts, int *pIndicies, udFloat2 minUV, udFloat2 maxUV, int collapseEdgeMask)
{
  for (int y = 0; y < TileIndexResolution; ++y)
  {
    for (int x = 0; x < TileIndexResolution; ++x)
    {
      int index = y * TileIndexResolution + x;
      int vertIndex = y * TileVertexResolution + x;

      // TODO: once figured out, remove commented code from all the conditionals statements below
      pIndicies[index * 6 + 0] = vertIndex + TileVertexResolution;
      pIndicies[index * 6 + 1] = vertIndex + 1;
      pIndicies[index * 6 + 2] = vertIndex;

      pIndicies[index * 6 + 3] = vertIndex + TileVertexResolution;
      pIndicies[index * 6 + 4] = vertIndex + TileVertexResolution + 1;
      pIndicies[index * 6 + 5] = vertIndex + 1;


      // corner cases
      if ((collapseEdgeMask & down) && (collapseEdgeMask & right) && x >= (TileIndexResolution - 2) && y >= (TileIndexResolution - 2))
      {
        if (x == TileIndexResolution - 2 && y == TileIndexResolution - 2)
        {
          //pIndicies[index * 6 + 0] = vertIndex + TileVertexResolution;
          //pIndicies[index * 6 + 1] = vertIndex + 1;
          //pIndicies[index * 6 + 2] = vertIndex;

          //pIndicies[index * 6 + 3] = vertIndex + TileVertexResolution;
          //pIndicies[index * 6 + 4] = vertIndex + TileVertexResolution + 1;
          //pIndicies[index * 6 + 5] = vertIndex + 1;
        }
        else if (x == TileIndexResolution - 1 && y == TileIndexResolution - 2)
        {
          //pIndicies[index * 6 + 0] = vertIndex + TileVertexResolution;
          //pIndicies[index * 6 + 1] = vertIndex + 1;
          //pIndicies[index * 6 + 2] = vertIndex;

          // collapse
          pIndicies[index * 6 + 3] = vertIndex;
          pIndicies[index * 6 + 4] = vertIndex;
          pIndicies[index * 6 + 5] = vertIndex;
        }
        else if (x == TileIndexResolution - 2 && y == TileIndexResolution - 1)
        {
          //pIndicies[index * 6 + 0] = vertIndex + TileVertexResolution;
          //pIndicies[index * 6 + 1] = vertIndex + 1;
          //pIndicies[index * 6 + 2] = vertIndex;

          // collapse
          pIndicies[index * 6 + 3] = vertIndex;
          pIndicies[index * 6 + 4] = vertIndex;
          pIndicies[index * 6 + 5] = vertIndex;
        }
        else // x == 1 && y == 1
        {
          pIndicies[index * 6 + 0] = vertIndex + TileVertexResolution - 1;
          pIndicies[index * 6 + 1] = vertIndex + TileVertexResolution + 1;
          pIndicies[index * 6 + 2] = vertIndex;

          pIndicies[index * 6 + 3] = vertIndex + TileVertexResolution + 1;
          pIndicies[index * 6 + 4] = vertIndex - TileVertexResolution + 1;
          pIndicies[index * 6 + 5] = vertIndex;
        }
      }
      else if ((collapseEdgeMask & down) && (collapseEdgeMask & left) && x <= 1 && y >= (TileIndexResolution - 2))
      {
        if (x == 0 && y == TileIndexResolution - 2)
        {
          // collapse
          pIndicies[index * 6 + 0] = vertIndex + 1;
          //pIndicies[index * 6 + 1] = vertIndex + 1;
          pIndicies[index * 6 + 2] = vertIndex + 1;

          // re-orient
          pIndicies[index * 6 + 3] = vertIndex + TileVertexResolution + 1;
          pIndicies[index * 6 + 4] = vertIndex + 1;
          pIndicies[index * 6 + 5] = vertIndex;
        }
        else if (x == 1 && y == TileIndexResolution - 2)
        {
          //pIndicies[index * 6 + 0] = vertIndex + TileVertexResolution;
          //pIndicies[index * 6 + 1] = vertIndex + 1;
          //pIndicies[index * 6 + 2] = vertIndex;
          //
          //pIndicies[index * 6 + 3] = vertIndex + TileVertexResolution;
          //pIndicies[index * 6 + 4] = vertIndex + TileVertexResolution + 1;
          //pIndicies[index * 6 + 5] = vertIndex + 1;
        }
        else if (x == 0 && y == TileIndexResolution - 1)
        {
          // re-orient
          pIndicies[index * 6 + 0] = vertIndex + TileVertexResolution;
          pIndicies[index * 6 + 1] = vertIndex + 1;
          pIndicies[index * 6 + 2] = vertIndex - TileVertexResolution;

          // collapse
          pIndicies[index * 6 + 3] = vertIndex + TileVertexResolution;
          pIndicies[index * 6 + 4] = vertIndex + TileVertexResolution + 2;
          pIndicies[index * 6 + 5] = vertIndex + 1;
        }
        else // x == 1 && y == 1
        {
          // collapse
          //pIndicies[index * 6 + 0] = vertIndex + TileVertexResolution;
          pIndicies[index * 6 + 1] = vertIndex + TileVertexResolution;
          pIndicies[index * 6 + 2] = vertIndex + TileVertexResolution;

          // re-orient
          pIndicies[index * 6 + 3] = vertIndex + TileVertexResolution + 1;
          pIndicies[index * 6 + 4] = vertIndex + 1;
          pIndicies[index * 6 + 5] = vertIndex;
        }
      }
      else if ((collapseEdgeMask & up) && (collapseEdgeMask & right) && x >= (TileIndexResolution - 2) && y <= 1)
      {
        if (x == TileIndexResolution - 2 && y == 0)
        {
          // collapse
          pIndicies[index * 6 + 0] = vertIndex;
          pIndicies[index * 6 + 1] = vertIndex;
          pIndicies[index * 6 + 2] = vertIndex;

          // re-orient triangle
          //pIndicies[index * 6 + 3] = vertIndex + TileVertexResolution;
          //pIndicies[index * 6 + 4] = vertIndex + TileVertexResolution + 1;
          pIndicies[index * 6 + 5] = vertIndex;
        }
        else if (x == TileIndexResolution - 1 && y == 0)
        {
          pIndicies[index * 6 + 0] = vertIndex + TileVertexResolution;
          pIndicies[index * 6 + 1] = vertIndex + 1;
          pIndicies[index * 6 + 2] = vertIndex - 1;

          pIndicies[index * 6 + 3] = vertIndex + TileVertexResolution;
          pIndicies[index * 6 + 4] = vertIndex + TileVertexResolution + TileVertexResolution + 1;
          pIndicies[index * 6 + 5] = vertIndex + 1;
        }
        else if (x == TileIndexResolution - 2 && y == 1)
        {
          //pIndicies[index * 6 + 0] = vertIndex + TileVertexResolution;
          //pIndicies[index * 6 + 1] = vertIndex + 1;
          //pIndicies[index * 6 + 2] = vertIndex;
          //
          //pIndicies[index * 6 + 3] = vertIndex + TileVertexResolution;
          //pIndicies[index * 6 + 4] = vertIndex + TileVertexResolution + 1;
          //pIndicies[index * 6 + 5] = vertIndex + 1;
        }
        else // x == TileIndexResolution -1 && y == 1
        {
          // collapse
          pIndicies[index * 6 + 0] = vertIndex;
          pIndicies[index * 6 + 1] = vertIndex;
          //pIndicies[index * 6 + 2] = vertIndex;

          // re-orient triangle
          //pIndicies[index * 6 + 3] = vertIndex + TileVertexResolution;
          //pIndicies[index * 6 + 4] = vertIndex + TileVertexResolution + 1;
          pIndicies[index * 6 + 5] = vertIndex;
        }
      }
      else if ((collapseEdgeMask & up) && (collapseEdgeMask & left) && x <= 1 && y <= 1)
      {
        if (x == 0 && y == 0)
        {
          pIndicies[index * 6 + 0] = vertIndex + TileVertexResolution + 1;
          pIndicies[index * 6 + 1] = vertIndex + 2;
          pIndicies[index * 6 + 2] = vertIndex;

          pIndicies[index * 6 + 3] = vertIndex + TileVertexResolution + TileVertexResolution;
          pIndicies[index * 6 + 4] = vertIndex + TileVertexResolution + 1;
          pIndicies[index * 6 + 5] = vertIndex;
        }
        else if (x == 1 && y == 0)
        {
          // collapse
          pIndicies[index * 6 + 0] = vertIndex + 1;
          //pIndicies[index * 6 + 1] = vertIndex + 1;
          pIndicies[index * 6 + 2] = vertIndex + 1;

          //pIndicies[index * 6 + 3] = vertIndex + TileVertexResolution;
          //pIndicies[index * 6 + 4] = vertIndex + TileVertexResolution + 1;
          //pIndicies[index * 6 + 5] = vertIndex + 1;
        }
        else if (x == 0 && y == 1)
        {
          // collapse
          pIndicies[index * 6 + 0] = vertIndex + 1;
          //pIndicies[index * 6 + 1] = vertIndex + 1;
          pIndicies[index * 6 + 2] = vertIndex + 1;

          //pIndicies[index * 6 + 3] = vertIndex + TileVertexResolution;
          //pIndicies[index * 6 + 4] = vertIndex + TileVertexResolution + 1;
          //pIndicies[index * 6 + 5] = vertIndex + 1;
        }
        else // x==1, y == 1
        {
          // do nothing
        }
      }
      else if (y == 0 && (collapseEdgeMask & up))
      {
        if ((x & 1) == 0)
        {
          //pIndicies[index * 6 + 0] = vertIndex + TileVertexResolution;
          pIndicies[index * 6 + 1] = vertIndex + 2;
          //pIndicies[index * 6 + 2] = vertIndex;

          //pIndicies[index * 6 + 3] = vertIndex + TileVertexResolution;
          //pIndicies[index * 6 + 4] = vertIndex + TileVertexResolution + 1;
          pIndicies[index * 6 + 5] = vertIndex + 2;
        }
        else
        {
          // collapse
          pIndicies[index * 6 + 0] = vertIndex + 1;
          //pIndicies[index * 6 + 1] = vertIndex + 1;
          pIndicies[index * 6 + 2] = vertIndex + 1;

          //pIndicies[index * 6 + 3] = vertIndex + TileVertexResolution;
          //pIndicies[index * 6 + 4] = vertIndex + TileVertexResolution + 1;
          //pIndicies[index * 6 + 5] = vertIndex + 1;
        }
      }
      else if (y == TileIndexResolution - 1 && (collapseEdgeMask & down))
      {
        if ((x & 1) == 0)
        {
          //pIndicies[index * 6 + 0] = vertIndex + TileVertexResolution;
          //pIndicies[index * 6 + 1] = vertIndex + 1;
          //pIndicies[index * 6 + 2] = vertIndex;

          // collapse
          //pIndicies[index * 6 + 3] = vertIndex + TileVertexResolution;
          pIndicies[index * 6 + 4] = vertIndex + TileVertexResolution;
          pIndicies[index * 6 + 5] = vertIndex + TileVertexResolution;
        }
        else
        {
          pIndicies[index * 6 + 0] = vertIndex + TileVertexResolution - 1;
          //pIndicies[index * 6 + 1] = vertIndex + 1;
          //pIndicies[index * 6 + 2] = vertIndex;

          pIndicies[index * 6 + 3] = vertIndex + TileVertexResolution - 1;
          //pIndicies[index * 6 + 4] = vertIndex + TileVertexResolution + 1;
          //pIndicies[index * 6 + 5] = vertIndex + 1;
        }
      }
      else if (x == TileIndexResolution - 1 && (collapseEdgeMask & right))
      {
        if ((y & 1) == 0)
        {
          // pIndicies[index * 6 + 0] = vertIndex + TileVertexResolution;
           //pIndicies[index * 6 + 1] = vertIndex + 1;
           //pIndicies[index * 6 + 2] = vertIndex;

           // collapse
          pIndicies[index * 6 + 3] = vertIndex + 1;
          pIndicies[index * 6 + 4] = vertIndex + 1;
          //pIndicies[index * 6 + 5] = vertIndex + 1;
        }
        else
        {
          //pIndicies[index * 6 + 0] = vertIndex + TileVertexResolution;
          pIndicies[index * 6 + 1] = vertIndex - TileVertexResolution + 1;
          //pIndicies[index * 6 + 2] = vertIndex;

          //pIndicies[index * 6 + 3] = vertIndex + TileVertexResolution;
          //pIndicies[index * 6 + 4] = vertIndex + TileVertexResolution + 1;
          pIndicies[index * 6 + 5] = vertIndex - TileVertexResolution + 1;
        }
      }
      else if (x == 0 && (collapseEdgeMask & left))
      {
        if ((y & 1) == 0)
        {
          pIndicies[index * 6 + 0] = vertIndex + TileVertexResolution + TileVertexResolution;
          //pIndicies[index * 6 + 1] = vertIndex + 1;
          //pIndicies[index * 6 + 2] = vertIndex;

          pIndicies[index * 6 + 3] = vertIndex + TileVertexResolution + TileVertexResolution;
          //pIndicies[index * 6 + 4] = vertIndex + TileVertexResolution + 1;
          //pIndicies[index * 6 + 5] = vertIndex + 1;
        }
        else
        {
          // collapse
          //pIndicies[index * 6 + 0] = vertIndex + TileVertexResolution;
          pIndicies[index * 6 + 1] = vertIndex + TileVertexResolution;
          pIndicies[index * 6 + 2] = vertIndex + TileVertexResolution;

          //pIndicies[index * 6 + 3] = vertIndex + TileVertexResolution;
          //pIndicies[index * 6 + 4] = vertIndex + TileVertexResolution + 1;
          //pIndicies[index * 6 + 5] = vertIndex + 1;
        }
      }
    }
  }

  float normalizeVertexPositionScale = float(TileVertexResolution) / (TileVertexResolution - 1); // ensure verts are [0, 1]
  for (int y = 0; y < TileVertexResolution; ++y)
  {
    for (int x = 0; x < TileVertexResolution; ++x)
    {
      uint32_t index = y * TileVertexResolution + x;
      float normX = ((float)(x) / TileVertexResolution) * normalizeVertexPositionScale;
      float normY = ((float)(y) / TileVertexResolution) * normalizeVertexPositionScale;
      pVerts[index].position.x = minUV.x + normX * (maxUV.x - minUV.x);
      pVerts[index].position.y = minUV.y + normY * (maxUV.y - minUV.y);
      pVerts[index].position.z = (float)index;
      //
      //float normX = ((float)(x) / TileVertexResolution) * normalizeVertexPositionScale;
      //float normY = ((float)(y) / TileVertexResolution) * normalizeVertexPositionScale;
      //pVerts[index].position.x = minUV.x + normX * (maxUV.x - minUV.x);
      //pVerts[index].position.y = minUV.y + normY * (maxUV.y - minUV.y);
      
    }
  }
}

udResult vcTileRenderer_Create(vcTileRenderer **ppTileRenderer, vcSettings *pSettings)
{
  udResult result;
  vcTileRenderer *pTileRenderer = nullptr;
  vcP3Vertex verts[TileVertexResolution * TileVertexResolution] = {};
  int indicies[TileIndexResolution * TileIndexResolution * 6] = {};
  uint32_t greyPixel = 0xf3f3f3ff;
  const char *tiles[] = { "D:\\git\\vaultclient\\builds\\assets\\S28E152.hgt", "D:\\git\\vaultclient\\builds\\assets\\S28E153.hgt" };
  UD_ERROR_NULL(ppTileRenderer, udR_InvalidParameter_);

  pTileRenderer = udAllocType(vcTileRenderer, 1, udAF_Zero);
  UD_ERROR_NULL(pTileRenderer, udR_MemoryAllocationFailure);

  vcQuadTree_Create(&pTileRenderer->quadTree, pSettings);

  pTileRenderer->pSettings = pSettings;

  pTileRenderer->cache.pSemaphore = udCreateSemaphore();
  pTileRenderer->cache.pMutex = udCreateMutex();
  pTileRenderer->cache.keepLoading = true;
  pTileRenderer->cache.tileLoadList.Init(128);
  pTileRenderer->cache.tileTimeoutList.Init(128);

  for (size_t i = 0; i < udLengthOf(pTileRenderer->cache.pThreads); ++i)
    UD_ERROR_CHECK(udThread_Create(&pTileRenderer->cache.pThreads[i], vcTileRenderer_LoadThread, pTileRenderer));

  UD_ERROR_IF(!vcShader_CreateFromText(&pTileRenderer->presentShader.pProgram, g_tileVertexShader, g_tileFragmentShader, vcP3VertexLayout), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetConstantBuffer(&pTileRenderer->presentShader.pConstantBuffer, pTileRenderer->presentShader.pProgram, "u_EveryObject", sizeof(pTileRenderer->presentShader.everyObject)), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetSamplerIndex(&pTileRenderer->presentShader.uniform_texture, pTileRenderer->presentShader.pProgram, "u_texture"), udR_InternalError);
  vcShader_GetSamplerIndex(&pTileRenderer->presentShader.uniform_dem0, pTileRenderer->presentShader.pProgram, "u_dem0");
  vcShader_GetSamplerIndex(&pTileRenderer->presentShader.uniform_dem1, pTileRenderer->presentShader.pProgram, "u_dem1");

  // build mesh variants
  for (int i = 0; i < 16; ++i)
  {
    vcTileRenderer_BuildMeshVertices(verts, indicies, udFloat2::create(0.0f, 0.0f), udFloat2::create(1.0f, 1.0f), meshConfigurations[i]);
    vcMesh_Create(&pTileRenderer->pTileMeshes[i], vcP3VertexLayout, (int)udLengthOf(vcP3VertexLayout), verts, TileVertexResolution * TileVertexResolution, indicies, TileIndexResolution * TileIndexResolution * 6);
  }

  UD_ERROR_CHECK(vcTexture_Create(&pTileRenderer->pEmptyTileTexture, 1, 1, &greyPixel));

  for (int i = 0; i < 2; ++i)
  {
    void *pFileData;
    int64_t fileLen;
    udFile_Load(tiles[i], &pFileData, &fileLen);
    int outputSize = 3601;
    int inputSize = 3601;
    uint16_t *pRealignedPixels = udAllocType(uint16_t, outputSize * outputSize, udAF_Zero);
    uint16_t lastValidHeight = 0;
    for (int y = 0; y < outputSize; ++y)
    {
      for (int x = 0; x < outputSize; ++x)
      {
        uint16_t p = ((uint16_t *)pFileData)[y * inputSize + x];
        p = ((p & 0xff00) >> 8) | ((p & 0x00ff) << 8);
        if (p > 40000)
          p = lastValidHeight;
        lastValidHeight = p;
        pRealignedPixels[y * outputSize + x] = p;
      }
    }
    vcTexture_Create(&pDEMTexture[i], outputSize, outputSize, pRealignedPixels, vcTextureFormat_R16, vcTFM_Linear, false, vcTWM_Clamp);//, vcTCF_None, 16);
  }
  pTileRenderer->pTransparentTiles = new std::vector<vcQuadTreeNode*>();
  pTileRenderer->pRenderQueue = new std::vector<std::vector<vcQuadTreeNode*>>();
  for (int i = 0; i < MaxVisibleTileLevel; ++i)
    pTileRenderer->pRenderQueue->push_back(std::vector<vcQuadTreeNode*>());

  *ppTileRenderer = pTileRenderer;
  pTileRenderer = nullptr;
  result = udR_Success;

epilogue:
  if (pTileRenderer)
    vcTileRenderer_Destroy(&pTileRenderer);

  return result;
}

udResult vcTileRenderer_Destroy(vcTileRenderer **ppTileRenderer)
{
  if (ppTileRenderer == nullptr || *ppTileRenderer == nullptr)
    return udR_InvalidParameter_;

  vcTileRenderer *pTileRenderer = *ppTileRenderer;

  pTileRenderer->cache.keepLoading = false;

  for (size_t i = 0; i < udLengthOf(pTileRenderer->cache.pThreads); ++i)
    udIncrementSemaphore(pTileRenderer->cache.pSemaphore);

  for (size_t i = 0; i < udLengthOf(pTileRenderer->cache.pThreads); ++i)
  {
    udThread_Join(pTileRenderer->cache.pThreads[i]);
    udThread_Destroy(&pTileRenderer->cache.pThreads[i]);
  }

  udDestroyMutex(&pTileRenderer->cache.pMutex);
  udDestroySemaphore(&pTileRenderer->cache.pSemaphore);

  pTileRenderer->cache.tileLoadList.Deinit();
  pTileRenderer->cache.tileTimeoutList.Deinit();

  vcShader_ReleaseConstantBuffer(pTileRenderer->presentShader.pProgram, pTileRenderer->presentShader.pConstantBuffer);
  vcShader_DestroyShader(&(pTileRenderer->presentShader.pProgram));
  for (int i = 0; i < 16; ++i)
    vcMesh_Destroy(&pTileRenderer->pTileMeshes[i]);
  vcTexture_Destroy(&pTileRenderer->pEmptyTileTexture);

  delete pTileRenderer->pTransparentTiles;
  delete pTileRenderer->pRenderQueue;
  pTileRenderer->pTransparentTiles = nullptr;
  pTileRenderer->pRenderQueue = nullptr;

  vcQuadTree_Destroy(&(*ppTileRenderer)->quadTree);
  udFree(*ppTileRenderer);
  *ppTileRenderer = nullptr;

  vcTexture_Destroy(&pDEMTexture[0]);
  vcTexture_Destroy(&pDEMTexture[1]);
  return udR_Success;
}

bool vcTileRenderer_UpdateTileTexture(vcTileRenderer *pTileRenderer, vcQuadTreeNode *pNode)
{
  vcTileRenderer::vcTileCache *pTileCache = &pTileRenderer->cache;
  if (pNode->renderInfo.loadStatus == vcNodeRenderInfo::vcTLS_None)
  {
    pNode->renderInfo.loadStatus = vcNodeRenderInfo::vcTLS_InQueue;

    pNode->renderInfo.pData = nullptr;
    pNode->renderInfo.pTexture = nullptr;
    pNode->renderInfo.timeoutTime = pTileRenderer->totalTime;

    pTileCache->tileLoadList.PushBack(pNode);
    udIncrementSemaphore(pTileCache->pSemaphore);
  }

  pNode->renderInfo.tryLoad = true;

  if (pNode->renderInfo.loadStatus == vcNodeRenderInfo::vcTLS_Downloaded)
  {
    pNode->renderInfo.loadStatus = vcNodeRenderInfo::vcTLS_Loaded;
    pNode->renderInfo.tryLoad = false;

    vcTexture_Create(&pNode->renderInfo.pTexture, pNode->renderInfo.width, pNode->renderInfo.height, pNode->renderInfo.pData, vcTextureFormat_RGBA8, vcTFM_Linear, true, vcTWM_Clamp, vcTCF_None, 16);
    udFree(pNode->renderInfo.pData);

    return true;
  }

  return false;
}

void vcTileRenderer_UpdateTextureQueuesRecursive(vcTileRenderer *pTileRenderer, vcQuadTreeNode *pNode, int &tileUploadCount)
{
  if (tileUploadCount >= MaxTextureUploadsPerFrame)
    return;

  if (!vcQuadTree_IsLeafNode(pNode))
  {
    for (int c = 0; c < 4; ++c)
      vcTileRenderer_UpdateTextureQueuesRecursive(pTileRenderer, &pTileRenderer->quadTree.nodes.pPool[pNode->childBlockIndex + c], tileUploadCount);
  }

  if (pNode->renderInfo.loadStatus != vcNodeRenderInfo::vcTLS_Loaded && vcQuadTree_IsVisibleLeafNode(&pTileRenderer->quadTree, pNode))
  {
    if (vcTileRenderer_UpdateTileTexture(pTileRenderer, pNode))
      ++tileUploadCount;
  }

  //if (pNode->renderInfo.loadStatus == vcNodeRenderInfo::vcTLS_Loaded && !pNode->renderInfo.fadingIn && pNode->renderInfo.transparency == 0.0f)
  //{
  //  pNode->renderInfo.fadingIn = true;
  //  pTileRenderer->pTransparentTiles->push_back(pNode);
  //}
}

void vcTileRenderer_UpdateTextureQueues(vcTileRenderer *pTileRenderer)
{
  // invalidate current queue
  for (size_t i = 0; i < pTileRenderer->cache.tileLoadList.length; ++i)
    pTileRenderer->cache.tileLoadList[i]->renderInfo.tryLoad = false;

  // Limit the max number of tiles uploaded per frame
  // TODO: use timings instead
  int tileUploadCount = 0;

  // update visible tiles textures
  vcTileRenderer_UpdateTextureQueuesRecursive(pTileRenderer, &pTileRenderer->quadTree.nodes.pPool[pTileRenderer->quadTree.rootIndex], tileUploadCount);

  // always request root
  //vcTileRenderer_UpdateTileTexture(pTileRenderer, &pTileRenderer->quadTree.nodes.pPool[pTileRenderer->quadTree.rootIndex]);

  // remove from the queue any tiles that are invalid
  for (int i = 0; i < (int)pTileRenderer->cache.tileLoadList.length; ++i)
  {
    if (!pTileRenderer->cache.tileLoadList[i]->renderInfo.tryLoad)
    {
      pTileRenderer->cache.tileLoadList[i]->renderInfo.loadStatus = vcNodeRenderInfo::vcTLS_None;
      pTileRenderer->cache.tileLoadList.RemoveSwapLast(i);
      --i;
    }
  }

  // Note: this is a FIFO queue, so only need to check the head
  while (pTileRenderer->cache.tileTimeoutList.length > 0)
  {
    vcQuadTreeNode *pNode = pTileRenderer->cache.tileTimeoutList[0];
    if (vcTileRenderer_ShouldLoadNode(pNode) && (pNode->renderInfo.timeoutTime - pTileRenderer->totalTime) > 0.0f)
      break;

    pTileRenderer->cache.tileLoadList.PushBack(pNode);
    pTileRenderer->cache.tileTimeoutList.RemoveSwapLast(0);
  }

  // TODO: For each tile in cache, LRU destroy
}

void vcTileRenderer_Update(vcTileRenderer *pTileRenderer, const double deltaTime, vcGISSpace *pSpace, const udDouble3 worldCorners[4], const udInt3 &slippyCoords, const udDouble3 &cameraWorldPos, const udDouble4x4 &viewProjectionMatrix)
{
  pTileRenderer->frameDeltaTime = (float)deltaTime;
  pTileRenderer->totalTime += pTileRenderer->frameDeltaTime;
  pTileRenderer->cameraPosition = cameraWorldPos;

  double slippyCornersViewSize = udMag3(worldCorners[1] - worldCorners[2]) * 0.5;
  vcQuadTreeViewInfo viewInfo =
  {
    pSpace,
    slippyCoords,
    cameraWorldPos,
    slippyCornersViewSize,
    pTileRenderer->pSettings->maptiles.mapHeight,
    viewProjectionMatrix,
    MaxVisibleTileLevel
  };

  uint64_t startTime = udPerfCounterStart();

  vcQuadTree_Update(&pTileRenderer->quadTree, viewInfo);

  udLockMutex(pTileRenderer->cache.pMutex);
  vcTileRenderer_UpdateTextureQueues(pTileRenderer);
  udReleaseMutex(pTileRenderer->cache.pMutex);

  pTileRenderer->quadTree.metaData.generateTimeMs = udPerfCounterMilliseconds(startTime);
}

bool vcTileRenderer_NodeHasValidBounds(vcQuadTreeNode *pNode)
{
  return !((pNode->tileExtents.x <= UD_EPSILON || pNode->tileExtents.y <= UD_EPSILON) ||
    pNode->worldBounds[0].x > pNode->worldBounds[1].x || pNode->worldBounds[2].x > pNode->worldBounds[3].x ||
    pNode->worldBounds[2].y > pNode->worldBounds[0].y || pNode->worldBounds[3].y > pNode->worldBounds[1].y);
}

bool vcTileRenderer_IsRootNode(vcTileRenderer *pTileRenderer, vcQuadTreeNode *pNode)
{
  return (pNode == &pTileRenderer->quadTree.nodes.pPool[pTileRenderer->quadTree.rootIndex]);
}

bool vcTileRenderer_CanNodeDraw(vcQuadTreeNode *pNode)
{
  if (!pNode->renderInfo.pTexture)// || pNode->renderInfo.fadingIn)
    return false;

  return vcTileRenderer_NodeHasValidBounds(pNode);
}

#include "udGeoZone.h"
#include "vcGIS.h"

bool vcTileRenderer_DrawNode(vcTileRenderer *pTileRenderer, vcQuadTreeNode *pNode, vcMesh *pMesh, const udDouble4x4 &view, bool parentCanDraw)
{
  vcTexture *pTexture = pNode->renderInfo.pDrawTexture;
//  float tileTransparency = pNode->renderInfo.transparency * pTileRenderer->pSettings->maptiles.transparency;
  if (pTexture == nullptr)
  {
//#if !VISUALIZE_DEBUG_TILES
//    if (!vcTileRenderer_IsRootNode(pTileRenderer, pNode) && parentCanDraw)
//      return false;
//#endif
//
    pTexture = pTileRenderer->pEmptyTileTexture;
//    tileTransparency = pTileRenderer->pSettings->maptiles.transparency;
  }

  //int slippyLayerDescendAmount = 1;//udMin((MAX_SLIPPY_LEVEL - slippyTileCoord.z), gSlippyLayerDescendAmount[3]);
  //udDouble3 tileBounds[9];
  //for (int t = 0; t < 9; ++t)
  //{
  //  udInt2 slippySampleCoord = udInt2::create((pNode->slippyPosition.x * (1 << slippyLayerDescendAmount)) + (t % 3),
  //    (pNode->slippyPosition.y * (1 << slippyLayerDescendAmount)) + (t / 3));
  //  vcGIS_SlippyToLocal(&pTileRenderer->quadTree.gisSpace, &tileBounds[t], slippySampleCoord, pNode->slippyPosition.z + slippyLayerDescendAmount);
  //  //tileBounds[t] = udDouble2::create(localCorners[t].x, localCorners[t].y);
  //}

  for (int t = 0; t < 4; ++t)//TileVertexResolution * TileVertexResolution; ++t)
  {
    udFloat4 eyeSpaceVertexPosition = udFloat4::create(view * udDouble4::create(pNode->worldBounds[t], 0.0, 1.0));
    pTileRenderer->presentShader.everyObject.eyePositions[t] = eyeSpaceVertexPosition;
  }

  //for (int t = 0; t < TileVertexResolution * TileVertexResolution; ++t)
  //{
  //  udFloat4 eyeSpaceVertexPosition = udFloat4::create(view * udDouble4::create(pNode->worldBounds[t], 0.0, 1.0));
  //  pTileRenderer->presentShader.everyObject.eyePositions[t] = eyeSpaceVertexPosition;
  //}

  //pTileRenderer->presentShader.everyObject.colour = udFloat4::create(1.f, 1.f, 1.f, tileTransparency);
//#if VISUALIZE_DEBUG_TILES
//  pTileRenderer->presentShader.everyObject.colour.w = 1.0f;
//  if (!pNode->renderInfo.pTexture)
//  {
//    pTileRenderer->presentShader.everyObject.colour.x = pNode->level / 21.0f;
//    if (!pNode->visible)
//      pTileRenderer->presentShader.everyObject.colour.z = pNode->level / 21.0f;
//  }
//#endif

  //S28E153
  udDouble3 r0 = udGeoZone_LatLongToCartesian(pTileRenderer->quadTree.gisSpace.zone, udDouble3::create(-28.0, 152, 0));
  udDouble3 r1 = udGeoZone_LatLongToCartesian(pTileRenderer->quadTree.gisSpace.zone, udDouble3::create(-27.0, 153, 0));

  //S28E152
  udDouble3 r2 = udGeoZone_LatLongToCartesian(pTileRenderer->quadTree.gisSpace.zone, udDouble3::create(-28.0, 153.0, 0));
  udDouble3 r3 = udGeoZone_LatLongToCartesian(pTileRenderer->quadTree.gisSpace.zone, udDouble3::create(-27.0, 154.0, 0));
  // left brisbane (S28E152.hgt), right brisbane (S28E153.hgt)
  //udDouble2 mins[] = { udDouble2::create(400781.82958118513, 6902797.6293904129), udDouble2::create(500000.00000000000, 6902394.7726541311) };
  //udDouble2 maxs[] = { udDouble2::create(500000.00000000000, 7013171.6474111192), udDouble2::create(598325.33504640602, 7013564.7575185951) };


  //{x = 400781.82958145085 y = 12986828.352577262 z = 0.00000000000000000 }

  // left brisbane (S28E152.hgt), right brisbane (S28E153.hgt)
  //udDouble2 mins[] = { udDouble2::create(400781.82958118513, 6902797.6293904129), udDouble2::create(500000.00000000000, 6902394.7726541311) };
  //udDouble2 maxs[] = { udDouble2::create(500000.00000000000, 7013171.6474111192), udDouble2::create(598325.33504640602, 7013564.7575185951) };
  //udDouble2 mins[] = { udDouble2::create(401674.66495384468, 6902394.7726660399), udDouble2::create(500000.00000000000, 6902797.6294023134) };
  //udDouble2 maxs[] = { udDouble2::create(500000.00000000000, 6902797.6294023134), udDouble2::create(598325.33504615526, 6902394.7726660399) };
  udDouble2 mins[] = { udDouble2::create(r0.x, r0.y), udDouble2::create(r2.x, r2.y) };
  udDouble2 maxs[] = { udDouble2::create(r1.x, r1.y), udDouble2::create(r3.x, r3.y) };

  mins[0].x += 100.0; // manual correction (because its busted)
  maxs[0].y += 365.0; // manual correction (because its busted)
  udDouble2 ranges[] = { maxs[0] - mins[0], maxs[1] - mins[1] };

  bool in0Bounds = !(pNode->worldBounds[1].x < mins[0].x || pNode->worldBounds[0].x > maxs[0].x || pNode->worldBounds[1].y < mins[0].y || pNode->worldBounds[3].y > maxs[0].y);
  bool in1Bounds = !(pNode->worldBounds[1].x < mins[1].x || pNode->worldBounds[0].x > maxs[1].x || pNode->worldBounds[1].y < mins[1].y || pNode->worldBounds[3].y > maxs[1].y);
  if (!in0Bounds && !in1Bounds)
    return true;

  for (int d = 0; d < 2; ++d)
  {
    for (int t = 0; t < 4; ++t)
    {
      double u2 = (pNode->worldBounds[t].x - mins[d].x) / ranges[d].x;
      double v2 = (pNode->worldBounds[t].y - mins[d].y) / ranges[d].y;

      pTileRenderer->presentShader.everyObject.demUVs[d * 4 + t].x = float(u2);
      pTileRenderer->presentShader.everyObject.demUVs[d * 4 + t].y = float(1.0 - v2);
    }
  }

  pTileRenderer->presentShader.everyObject.colourUV = udFloat4::create(
    pNode->renderInfo.uvStart.x,
    pNode->renderInfo.uvStart.y,
    pNode->renderInfo.uvEnd.x - pNode->renderInfo.uvStart.x,
    pNode->renderInfo.uvEnd.y - pNode->renderInfo.uvStart.y);

  vcShader_BindTexture(pTileRenderer->presentShader.pProgram, pTexture, 0, pTileRenderer->presentShader.uniform_texture);
  vcShader_BindTexture(pTileRenderer->presentShader.pProgram, pDEMTexture[0], 1, pTileRenderer->presentShader.uniform_dem0);
  vcShader_BindTexture(pTileRenderer->presentShader.pProgram, pDEMTexture[1], 2, pTileRenderer->presentShader.uniform_dem1);


  vcShader_BindConstantBuffer(pTileRenderer->presentShader.pProgram, pTileRenderer->presentShader.pConstantBuffer, &pTileRenderer->presentShader.everyObject, sizeof(pTileRenderer->presentShader.everyObject));
  vcMesh_Render(pMesh, TileIndexResolution * TileIndexResolution * 2); // 2 tris per quad

  pNode->rendered = true;
  ++pTileRenderer->quadTree.metaData.nodeRenderCount;

  return true;
}

void vcTileRenderer_RecursiveSetRendered(vcTileRenderer *pTileRenderer, vcQuadTreeNode *pNode, bool rendered)
{
  pNode->rendered = pNode->rendered || rendered;
  if (!vcQuadTree_IsLeafNode(pNode))
  {
    for (int c = 0; c < 4; ++c)
      vcTileRenderer_RecursiveSetRendered(pTileRenderer, &pTileRenderer->quadTree.nodes.pPool[pNode->childBlockIndex + c], pNode->rendered);
  }
}

// 'true' indicates the node was able to render itself (or it didn't want to render itself).
// 'false' indicates that the nodes ancestor needs to be rendered.
bool vcTileRenderer_RecursiveBuildRenderQueue(vcTileRenderer *pTileRenderer, vcQuadTreeNode *pNode, vcQuadTreeNode *pBestTexturedAncestor)
{
  if (!pNode->touched)
  {
    vcQuadTreeNode *pParentNode = &pTileRenderer->quadTree.nodes.pPool[pNode->parentIndex];

    // parent can render itself (opaque), so this tile is not needed
    //if (pParentNode->renderInfo.pTexture && !pParentNode->renderInfo.fadingIn)
    //  return false;

    // re-test visibility
    pNode->visible = pParentNode->visible && vcQuadTree_IsNodeVisible(&pTileRenderer->quadTree, pNode);
  }

  if (!pNode->visible)
    return false;

  pNode->renderInfo.pDrawTexture = nullptr;
  pNode->renderInfo.uvStart = udFloat2::zero();
  pNode->renderInfo.uvEnd = udFloat2::one();
  if (vcTileRenderer_CanNodeDraw(pNode))
  {
    pNode->renderInfo.pDrawTexture = pNode->renderInfo.pTexture;
    pBestTexturedAncestor = pNode;
  }
  else if (pBestTexturedAncestor != nullptr)
  {
    // will be using best ancestor
    // TODO use morten indicies
    pNode->renderInfo.pDrawTexture = pBestTexturedAncestor->renderInfo.pDrawTexture;
    int depthDiff = pNode->level - pBestTexturedAncestor->level;

    // TODO: probably a better way of doing this...
    udDouble2 boundsRange = pBestTexturedAncestor->worldBounds[3] - pBestTexturedAncestor->worldBounds[0];
    pNode->renderInfo.uvStart = udFloat2::create((pNode->worldBounds[0] - pBestTexturedAncestor->worldBounds[0]) / boundsRange);
    pNode->renderInfo.uvEnd = udFloat2::create(udDouble2::one() - (pBestTexturedAncestor->worldBounds[3] - pNode->worldBounds[3]) / boundsRange);
  }

  bool childrenNeedThisTileRendered = vcQuadTree_IsLeafNode(pNode);
  if (!childrenNeedThisTileRendered)
  {
    for (int c = 0; c < 4; ++c)
    {
      vcQuadTreeNode *pChildNode = &pTileRenderer->quadTree.nodes.pPool[pNode->childBlockIndex + c];
      childrenNeedThisTileRendered = !vcTileRenderer_RecursiveBuildRenderQueue(pTileRenderer, pChildNode, pBestTexturedAncestor) || childrenNeedThisTileRendered;
    }

    return false; // atm only render leaf children
  }

  //if (pNode->renderInfo.fadingIn)
  //  return false;

  if (childrenNeedThisTileRendered)
  {
    //if (!pNode->renderInfo.pTexture && canParentDraw)//(!vcTileRenderer_IsRootNode(pTileRenderer, pNode) && canParentDraw))
    //  return false;

    //if (!vcTileRenderer_NodeHasValidBounds(pNode))
    //  return false;

    pTileRenderer->pRenderQueue->at(pNode->level).push_back(pNode);
  }

  // This child doesn't need parent to draw itself
  return true;
}

// Depth first rendering, using stencil to ensure no overdraw
void vcTileRenderer_DrawRenderQueue(vcTileRenderer *pTileRenderer, const udDouble4x4 &view)
{
  for (int i = MaxVisibleTileLevel - 1; i >= 0; --i)
  {
    for (size_t t = 0; t < pTileRenderer->pRenderQueue->at(i).size(); ++t)
    {
      vcQuadTreeNode *pNode = pTileRenderer->pRenderQueue->at(i).at(t);
      int meshIndex = 0;
      for (int mc = 0; mc < 16; ++mc)
      {
        if (meshConfigurations[mc] == pNode->neighbours)
        {
          meshIndex = mc;
          break;
        }
      }

      vcTileRenderer_DrawNode(pTileRenderer, pNode, pTileRenderer->pTileMeshes[meshIndex], view, false);
    }
  }
}

#include "gl/opengl/vcOpenGL.h"

void vcTileRenderer_Render(vcTileRenderer *pTileRenderer, const udDouble4x4 &view, const udDouble4x4 &proj)
{
  vcQuadTreeNode *pRootNode = &pTileRenderer->quadTree.nodes.pPool[pTileRenderer->quadTree.rootIndex];
  if (!pRootNode->touched) // can occur on failed re-roots
    return;

  for (int i = 0; i < MaxVisibleTileLevel; ++i)
    pTileRenderer->pRenderQueue->at(i).clear();

  udDouble4x4 viewWithMapTranslation = view * udDouble4x4::translation(0, 0, pTileRenderer->pSettings->maptiles.mapHeight);

  vcGLStencilSettings stencil = {};
  stencil.writeMask = 0xFF;
  stencil.compareFunc = vcGLSSF_Equal;
  stencil.compareValue = 0;
  stencil.compareMask = 0xFF;
  stencil.onStencilFail = vcGLSSOP_Keep;
  stencil.onDepthFail = vcGLSSOP_Keep;
  stencil.onStencilAndDepthPass = vcGLSSOP_Increment;

  vcGLState_SetFaceMode(vcGLSFM_Solid, vcGLSCM_None);
  vcGLState_SetDepthStencilMode(vcGLSDM_LessOrEqual, true, nullptr);//&stencil);

  if (pTileRenderer->pSettings->maptiles.transparency >= 1.0f)
    vcGLState_SetBlendMode(vcGLSBM_None);
  else
    vcGLState_SetBlendMode(vcGLSBM_Interpolative);

  if (pTileRenderer->pSettings->maptiles.blendMode == vcMTBM_Overlay)
  {
    vcGLState_SetViewportDepthRange(0.0f, 0.0f);
    vcGLState_SetDepthStencilMode(vcGLSDM_Always, false, &stencil);
  }
  else if (pTileRenderer->pSettings->maptiles.blendMode == vcMTBM_Underlay)
  {
    vcGLState_SetViewportDepthRange(1.0f, 1.0f);
  }

  for (int i = 0; i < 1; ++i)
  {
    vcShader_Bind(pTileRenderer->presentShader.pProgram);
    pTileRenderer->presentShader.everyObject.projectionMatrix = udFloat4x4::create(proj);
    pTileRenderer->presentShader.everyObject.viewMatrix = udFloat4x4::create(view);
    pTileRenderer->presentShader.everyObject.colour = udFloat4::create(0.f, 0.f, 0.f, 0.0f);//pTileRenderer->pSettings->maptiles.transparency);

    if (i == 1)
    {
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      glPolygonOffset(1.0f, -0.1f);
      pTileRenderer->presentShader.everyObject.colour = udFloat4::create(0.f, 1.f, 0.f, 1.f);//pTileRenderer->pSettings->maptiles.transparency);
      vcGLState_SetDepthStencilMode(vcGLSDM_LessOrEqual, true, nullptr);//&stencil);
    }

    //pRootNode->renderInfo.transparency = 1.0f;
    vcTileRenderer_RecursiveBuildRenderQueue(pTileRenderer, pRootNode, nullptr);
    vcTileRenderer_DrawRenderQueue(pTileRenderer, viewWithMapTranslation);

    // Render the root tile again (if it hasn't already been rendered normally) to cover up gaps between tiles
    //if (!pRootNode->rendered && pRootNode->renderInfo.pTexture && vcTileRenderer_NodeHasValidBounds(pRootNode))
    //  vcTileRenderer_DrawNode(pTileRenderer, pRootNode, pTileRenderer->pFullTileMesh, viewWithMapTranslation, false);

    // Draw transparent tiles
    /*if (pTileRenderer->pTransparentTiles->size() > 0)
    {
      // We know there will always be a stenciled opaque tile behind every transparent tile, so draw
      // with no depth testing, but stencil testing for map tiles
      stencil.writeMask = 0xFF;
      stencil.compareFunc = vcGLSSF_NotEqual;
      stencil.compareValue = 0;
      stencil.compareMask = 0xFF;
      stencil.onStencilFail = vcGLSSOP_Keep;
      stencil.onDepthFail = vcGLSSOP_Keep;
      stencil.onStencilAndDepthPass = vcGLSSOP_Keep;

      vcGLState_SetDepthStencilMode(vcGLSDM_Always, false, &stencil);
      vcGLState_SetBlendMode(vcGLSBM_Interpolative);
      for (auto tile : (*pTileRenderer->pTransparentTiles))
      {
        tile->renderInfo.transparency = udMin(1.0f, tile->renderInfo.transparency + pTileRenderer->frameDeltaTime * sTileFadeSpeed);
        if (tile->visible && vcTileRenderer_NodeHasValidBounds(tile))
          vcTileRenderer_DrawNode(pTileRenderer, tile, pTileRenderer->pFullTileMesh, viewWithMapTranslation, false);
      }

      for (int i = 0; i < int(pTileRenderer->pTransparentTiles->size()); ++i)
      {
        if (pTileRenderer->pTransparentTiles->at(i)->renderInfo.transparency >= 1.0f)
        {
          pTileRenderer->pTransparentTiles->at(i)->renderInfo.fadingIn = false;
          pTileRenderer->pTransparentTiles->erase(pTileRenderer->pTransparentTiles->begin() + i);
          --i;
        }
      }
    }
    */

    vcTileRenderer_RecursiveSetRendered(pTileRenderer, pRootNode, pRootNode->rendered);

    vcGLState_SetViewportDepthRange(0.0f, 1.0f);
    vcGLState_SetDepthStencilMode(vcGLSDM_LessOrEqual, true, nullptr);
    vcGLState_SetBlendMode(vcGLSBM_None);
    vcShader_Bind(nullptr);
  }

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glPolygonOffset(1.0f, 0.0f);

#if VISUALIZE_DEBUG_TILES
  printf("touched=%d, visible=%d, rendered=%d, leaves=%d\n", pTileRenderer->quadTree.metaData.nodeTouchedCount, pTileRenderer->quadTree.metaData.visibleNodeCount, pTileRenderer->quadTree.metaData.nodeRenderCount, pTileRenderer->quadTree.metaData.leafNodeCount);
#endif
}

void vcTileRenderer_ClearTiles(vcTileRenderer *pTileRenderer)
{
  udLockMutex(pTileRenderer->cache.pMutex);

  pTileRenderer->pTransparentTiles->clear();
  pTileRenderer->cache.tileLoadList.Clear();
  vcQuadTree_Reset(&pTileRenderer->quadTree);

  udReleaseMutex(pTileRenderer->cache.pMutex);
}
