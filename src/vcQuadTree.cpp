#include "vcQuadTree.h"
#include "vcGIS.h"
#include "gl/vcTexture.h"
#include "udPlatformUtil.h"

#define INVALID_NODE_INDEX 0xffffffff

enum
{
  NodeChildCount = 4,
};

// Returns -1=outside, 0=inside, >0=partial (bits of planes crossed)
static int vcQuadTree_FrustumTest(const udDouble4 frustumPlanes[6], const udDouble3 &boundCenter, const udDouble3 &boundExtents)
{
  int partial = 0;

  for (int i = 0; i < 6; ++i)
  {
    double distToCenter = udDot4(udDouble4::create(boundCenter, 1.0), frustumPlanes[i]);
    //optimized for case where boxExtents are all same: udFloat radiusBoxAtPlane = udDot3(boxExtents, udAbs(udVec3(curPlane)));
    double radiusBoxAtPlane = udDot3(boundExtents, udAbs(frustumPlanes[i].toVector3()));
    if (distToCenter < -radiusBoxAtPlane)
      return -1; // Box is entirely behind at least one plane
    else if (distToCenter <= radiusBoxAtPlane) // If spanned (not entirely infront)
      partial |= (1 << i);
  }

  return partial;
}

udResult vcQuadTree_ExpandCapacity(vcQuadTree *pQuadTree)
{
  udResult result = udR_Success;
  vcQuadTreeNode *pNewNodeMemory = nullptr;
  uint32_t newCapacity = pQuadTree->nodes.capacity * 2;

  pNewNodeMemory = udAllocType(vcQuadTreeNode, newCapacity, udAF_Zero);
  UD_ERROR_NULL(pNewNodeMemory, udR_MemoryAllocationFailure);

  if (pQuadTree->nodes.pPool)
  {
    memcpy(pNewNodeMemory, pQuadTree->nodes.pPool, pQuadTree->nodes.capacity * sizeof(vcQuadTreeNode));
    udFree(pQuadTree->nodes.pPool);
  }

  pQuadTree->nodes.pPool = pNewNodeMemory;
  pQuadTree->nodes.capacity = newCapacity;
  pQuadTree->metaData.memUsageMB = (sizeof(vcQuadTreeNode) * newCapacity) >> 20;

epilogue:
  return result;
}

uint32_t vcQuadTree_GetNodeIndex(vcQuadTree *pQuadTree, const udInt3 &slippy)
{
  vcQuadTreeNode *pNode = nullptr;
  for (uint32_t i = 0; i < pQuadTree->nodes.used; ++i)
  {
    pNode = &pQuadTree->nodes.pPool[i];
    if (pNode->isUsed && pNode->slippyPosition == slippy)
      return i;
  }

  return INVALID_NODE_INDEX;
}

uint32_t vcQuadTree_FindFreeChildBlock(vcQuadTree *pQuadTree)
{
  for (uint32_t blockIndex = 0; blockIndex < pQuadTree->nodes.used; blockIndex += NodeChildCount)
  {
    if (!pQuadTree->nodes.pPool[blockIndex].isUsed) // whole blocks are marked as unused
      return blockIndex;
  }

  if (pQuadTree->nodes.used >= pQuadTree->nodes.capacity)
    return INVALID_NODE_INDEX;

  pQuadTree->nodes.used += NodeChildCount;
  return pQuadTree->nodes.used - NodeChildCount;
}

double vcQuadTree_PointToRectDistance(udDouble3 edges[9], const udDouble3 &point)
{
  static const udInt2 edgePairs[] =
  {
    udInt2::create(0, 2), // top
    udInt2::create(0, 6), // left
    udInt2::create(6, 8), // bottom
    udInt2::create(2, 8), // right
  };
  
  double closestEdgeDistance = 0.0;
  
  // test each edge to find minimum distance to quadrant shape (3d)
  for (int e = 0; e < 4; ++e)
  {
    udDouble3 p1 = edges[edgePairs[e].x];
    udDouble3 p2 = edges[edgePairs[e].y];
  
    udDouble3 edge = p2 - p1;
    double r = udDot3(edge, (point - p1)) / udMagSq3(edge);
  
    // TODO: tile heights (DEM)
    udDouble3 closestPointOnEdge = p1 + udClamp(r, 0.0, 1.0) * edge;
  
    double distToEdge = udMag3(closestPointOnEdge - point);
    closestEdgeDistance = (e == 0) ? distToEdge : udMin(closestEdgeDistance, distToEdge);
  }
  
  return closestEdgeDistance;
}

void vcQuadTree_CleanupNode(vcQuadTreeNode *pNode)
{
  vcTexture_Destroy(&pNode->renderInfo.pTexture);
  udFree(pNode->renderInfo.pData);

  memset(pNode, 0, sizeof(vcQuadTreeNode));
}

#include "udGeoZone.h"
#include "vcGIS.h"

void vcQuadTree_CalculateNodeBounds(vcQuadTree *pQuadTree, vcQuadTreeNode *pNode)
{
  udDouble3 boundsMin = udDouble3::create(FLT_MAX, FLT_MAX, FLT_MAX);
  udDouble3 boundsMax = udDouble3::create(-FLT_MAX, -FLT_MAX, -FLT_MAX);
  for (int edge = 0; edge < 9; ++edge)
  {
    udInt2 slippySampleCoord = udInt2::create((pNode->slippyPosition.x * 2) + (edge % 3),
      (pNode->slippyPosition.y * 2) + (edge / 3));

    vcGIS_SlippyToLocal(&pQuadTree->gisSpace, &pNode->worldBounds[edge], slippySampleCoord, pNode->slippyPosition.z + 1);

    boundsMin = udMin(pNode->worldBounds[edge], boundsMin);
    boundsMax = udMax(pNode->worldBounds[edge], boundsMax);
  }

  pNode->tileCenter = (boundsMax + boundsMin) * 0.5;
  pNode->tileExtents = (boundsMax - boundsMin) * 0.5;
}

void vcQuadTree_InitNode(vcQuadTree *pQuadTree, uint32_t slotIndex, const udInt3 &childSlippy)
{
  vcQuadTreeNode *pNode = &pQuadTree->nodes.pPool[slotIndex];

  pNode->isUsed = true;
  pNode->childBlockIndex = INVALID_NODE_INDEX;
  pNode->parentIndex = INVALID_NODE_INDEX;
  pNode->slippyPosition = childSlippy;

  vcQuadTree_CalculateNodeBounds(pQuadTree, pNode);
}

bool vcQuadTree_IsNodeVisible(const vcQuadTree *pQuadTree, const vcQuadTreeNode *pNode)
{
  return -1 < vcQuadTree_FrustumTest(pQuadTree->frustumPlanes, pNode->tileCenter, pNode->tileExtents);
}

inline bool vcQuadTree_ShouldSubdivide(double distance, int depth)
{
  // trial and error'd this heuristic
  const int RootRegionSize = 16000000;
  return distance < (RootRegionSize >> depth);
}

void vcQuadTree_RecurseGenerateTree(vcQuadTree *pQuadTree, uint32_t currentNodeIndex, int currentDepth)
{
  vcQuadTreeNode *pCurrentNode = &pQuadTree->nodes.pPool[currentNodeIndex];
  pCurrentNode->childMask = 0;

  if (currentDepth >= pQuadTree->metaData.maxTreeDepth)
  {
    ++pQuadTree->metaData.leafNodeCount;
    return;
  }

  if (pCurrentNode->childBlockIndex == INVALID_NODE_INDEX)
  {
    uint32_t freeBlockIndex = vcQuadTree_FindFreeChildBlock(pQuadTree);
    if (freeBlockIndex == INVALID_NODE_INDEX)
      return;

    udInt3 childSlippy = udInt3::create(pCurrentNode->slippyPosition.x << 1, pCurrentNode->slippyPosition.y << 1, pCurrentNode->slippyPosition.z + 1);
    for (uint32_t i = 0; i < NodeChildCount; ++i)
      vcQuadTree_InitNode(pQuadTree, freeBlockIndex + i, childSlippy + udInt3::create(i & 1, i >> 1, 0));

    pCurrentNode->childBlockIndex = freeBlockIndex;
  }

  udInt2 pViewSlippyCoords;
  vcGIS_LocalToSlippy(&pQuadTree->gisSpace, &pViewSlippyCoords, pQuadTree->cameraWorldPosition, pQuadTree->slippyCoords.z + currentDepth + 1);

  udInt2 mortenIndices[] = { {0, 0}, {1, 0}, {0, 1}, {1, 1} };

  //subdivide
  // 0 == bottom left
  // 1 == bottom right
  // 2 == top left
  // 3 == top right
  for (uint32_t childQuadrant = 0; childQuadrant < NodeChildCount; ++childQuadrant)
  {
    pCurrentNode->childMask |= 1 << childQuadrant;

    uint32_t childIndex = pCurrentNode->childBlockIndex + childQuadrant;
    vcQuadTreeNode *pChildNode = &pQuadTree->nodes.pPool[childIndex];
    pChildNode->touched = true;

    // leave this here as we could be fixing up a re-root
    pChildNode->parentIndex = currentNodeIndex;
    pChildNode->visible = pCurrentNode->visible && vcQuadTree_IsNodeVisible(pQuadTree, pChildNode);
    pChildNode->morten.x = pCurrentNode->morten.x | (mortenIndices[childQuadrant].x << (31 - pChildNode->level));
    pChildNode->morten.y = pCurrentNode->morten.y | (mortenIndices[childQuadrant].y << (31 - pChildNode->level));

    // TODO: tile heights (DEM)
    double distanceToQuadrant = pQuadTree->cameraDistanceZeroAltitude;

    int32_t slippyManhattanDist = udAbs(pViewSlippyCoords.x - pChildNode->slippyPosition.x) + udAbs(pViewSlippyCoords.y - pChildNode->slippyPosition.y);
    if (slippyManhattanDist != 0)
      distanceToQuadrant = vcQuadTree_PointToRectDistance(pChildNode->worldBounds, pQuadTree->cameraWorldPosition);

    ++pQuadTree->metaData.nodeTouchedCount;
    if (pChildNode->visible)
      ++pQuadTree->metaData.visibleNodeCount;
    else if ((pQuadTree->pSettings->maptiles.mapOptions & vcTRF_OnlyRequestVisibleTiles) != 0)
      continue;

    int totalDepth = pQuadTree->slippyCoords.z + currentDepth;
    bool alwaysSubdivide = pChildNode->visible && totalDepth < 3;
    if (alwaysSubdivide || vcQuadTree_ShouldSubdivide(distanceToQuadrant, totalDepth))
      vcQuadTree_RecurseGenerateTree(pQuadTree, childIndex, currentDepth + 1);
    else
      ++pQuadTree->metaData.leafNodeCount;
  }
}


udInt2 decodeMorten(udInt2 &m, int d)
{
  int mask = ~(0xffffffff - ((2 << d) - 1));
  int shift = 31 - d;
  return { (m.x >> shift) &mask, (m.y >> shift) &mask };
}

vcQuadTreeNode *vcQuadTree_FindNodeFromMorten(vcQuadTree *pQuadTree, vcQuadTreeNode *pNode, const udInt2 &targetMorten, int targetDepth)
{
  if (vcQuadTree_IsLeafNode(pNode) || pNode->level >= targetDepth)
    return pNode;

  // TODO: handle outside bounds (e.g. morten.x < 0 || morten.y < 0 || morten.x >= ?? || morten.y >= ??)

  int shift = targetDepth - (pNode->level + 1);
  udInt2 thisLevel = { (targetMorten.x >> shift) & 1, (targetMorten.y >> shift) & 1 };
  int childIndex = thisLevel.x + thisLevel.y * 2;

  vcQuadTreeNode *pChildNode = &pQuadTree->nodes.pPool[pNode->childBlockIndex + childIndex];
  return vcQuadTree_FindNodeFromMorten(pQuadTree, pChildNode, targetMorten, targetDepth);
}

void vcQuadTree_CalculateNeighbours(vcQuadTree *pQuadTree, vcQuadTreeNode *pNode)
{
  if (vcQuadTree_IsLeafNode(pNode))
  {
    udInt2 morten = decodeMorten(pNode->morten, pNode->level);
    vcQuadTreeNode *pRootNode = &pQuadTree->nodes.pPool[pQuadTree->rootIndex];

    vcQuadTreeNode *pNeighbourUp = vcQuadTree_FindNodeFromMorten(pQuadTree, pRootNode, morten + udInt2::create(0, -1), pNode->level);
    vcQuadTreeNode *pNeighbourRight = vcQuadTree_FindNodeFromMorten(pQuadTree, pRootNode, morten + udInt2::create(1, 0), pNode->level);
    vcQuadTreeNode *pNeighbourDown = vcQuadTree_FindNodeFromMorten(pQuadTree, pRootNode, morten + udInt2::create(0, 1), pNode->level);
    vcQuadTreeNode *pNeighbourLeft = vcQuadTree_FindNodeFromMorten(pQuadTree, pRootNode, morten + udInt2::create(-1, 0), pNode->level);

    pNode->neighbours = 0;
    pNode->neighbours |= 0x1 * int(pNeighbourUp->level < pNode->level);
    pNode->neighbours |= 0x2 * int(pNeighbourRight->level < pNode->level);
    pNode->neighbours |= 0x4 * int(pNeighbourDown->level < pNode->level);
    pNode->neighbours |= 0x8 * int(pNeighbourLeft->level < pNode->level);
  }
  else
  {
    for (int c = 0; c < NodeChildCount; ++c)
    {
      vcQuadTreeNode *pChildNode = &pQuadTree->nodes.pPool[pNode->childBlockIndex + c];
      vcQuadTree_CalculateNeighbours(pQuadTree, pChildNode);
    }
  }
}

void vcQuadTree_CleanupNodes(vcQuadTree *pQuadTree)
{
  for (uint32_t i = 0; i < pQuadTree->nodes.used; ++i)
    vcQuadTree_CleanupNode(&pQuadTree->nodes.pPool[i]);

  pQuadTree->nodes.used = 0;
}

void vcQuadTree_Create(vcQuadTree *pQuadTree, vcSettings *pSettings)
{
  pQuadTree->nodes.capacity = 1024; // best guess, should hold enough nodes for usage
  vcQuadTree_ExpandCapacity(pQuadTree);

  vcQuadTree_CleanupNodes(pQuadTree);
  vcQuadTree_Reset(pQuadTree);

  pQuadTree->pSettings = pSettings;
}

void vcQuadTree_Destroy(vcQuadTree *pQuadTree)
{
  vcQuadTree_CleanupNodes(pQuadTree);

  udFree(pQuadTree->nodes.pPool);
  pQuadTree->nodes.capacity = 0;
}

void vcQuadTree_Reset(vcQuadTree *pQuadTree)
{
  memset(&pQuadTree->metaData, 0, sizeof(vcQuadTreeMetaData));
  pQuadTree->completeRerootRequired = true;
}

uint32_t vcQuadTree_NodeIndexToBlockIndex(uint32_t nodeIndex)
{
  return nodeIndex & ~3;
}

void vcQuadTree_InitRootBlock(vcQuadTree *pQuadTree)
{
  uint32_t rootBlockIndex = vcQuadTree_NodeIndexToBlockIndex(pQuadTree->rootIndex);
  for (uint32_t c = 0; c < NodeChildCount; ++c)
  {
    pQuadTree->nodes.pPool[rootBlockIndex + c].parentIndex = INVALID_NODE_INDEX;
    pQuadTree->nodes.pPool[rootBlockIndex + c].level = 0;
  }

  pQuadTree->nodes.pPool[rootBlockIndex].morten = udInt2::zero();
}

void vcQuadTree_Reroot(vcQuadTree *pQuadTree, const udInt3 &slippyCoords)
{
  // First determine if that node already exists
  uint32_t newRootIndex = vcQuadTree_GetNodeIndex(pQuadTree, slippyCoords);
  if (newRootIndex == INVALID_NODE_INDEX)
  {
    // Could not find node, add a new root one layer higher
    vcQuadTreeNode *pOldRoot = &pQuadTree->nodes.pPool[pQuadTree->rootIndex];
    uint32_t newRootBlockIndex = vcQuadTree_FindFreeChildBlock(pQuadTree);
    if (newRootBlockIndex == INVALID_NODE_INDEX || pOldRoot->slippyPosition.z < slippyCoords.z)
    {
      // signal a complete re-root is required
      pQuadTree->completeRerootRequired = true;
      return;
    }

    udInt3 newRootSlippy = udInt3::create(pOldRoot->slippyPosition.x >> 1, pOldRoot->slippyPosition.y >> 1, pOldRoot->slippyPosition.z - 1);
    udInt3 newRootBlockSlippy = udInt3::create((newRootSlippy.x >> 1) << 1, (newRootSlippy.y >> 1) << 1, newRootSlippy.z);
    for (uint32_t c = 0; c < NodeChildCount; ++c)
    {
      vcQuadTree_InitNode(pQuadTree, newRootBlockIndex + c, newRootBlockSlippy + udInt3::create(c & 1, c >> 1, 0));

      // preserve spatial position of root in its root block
      if (pQuadTree->nodes.pPool[newRootBlockIndex + c].slippyPosition == newRootSlippy)
      {
        uint32_t oldRootBlockIndex = vcQuadTree_NodeIndexToBlockIndex(pQuadTree->rootIndex);
        newRootIndex = newRootBlockIndex + c;

        // hook up old root and new root
        for (uint32_t c2 = 0; c2 < NodeChildCount; ++c2)
          pQuadTree->nodes.pPool[oldRootBlockIndex + c2].parentIndex = newRootIndex;

        pQuadTree->nodes.pPool[newRootIndex].childBlockIndex = oldRootBlockIndex;
      }
    }
  }
  else
  {
    // New root is existing node
    uint32_t rootBlockIndex = vcQuadTree_NodeIndexToBlockIndex(newRootIndex);
    for (uint32_t c = 0; c < NodeChildCount; ++c)
      pQuadTree->nodes.pPool[rootBlockIndex + c].parentIndex = INVALID_NODE_INDEX;
  }

  pQuadTree->rootIndex = newRootIndex;
  pQuadTree->completeRerootRequired = false;

  vcQuadTree_InitRootBlock(pQuadTree);
}

void vcQuadTree_CompleteReroot(vcQuadTree *pQuadTree, const udInt3 &slippyCoords)
{
  uint32_t freeBlockIndex = vcQuadTree_FindFreeChildBlock(pQuadTree);
  if (freeBlockIndex == INVALID_NODE_INDEX)
  {
    // complete re-root has failed, must wait for nodes to be freed at a later time
    return;
  }

  udInt3 rootBlockCoords = udInt3::create((slippyCoords.x >> 1) << 1, (slippyCoords.y >> 1) << 1, slippyCoords.z);
  for (uint32_t c = 0; c < NodeChildCount; ++c)
  {
    vcQuadTree_InitNode(pQuadTree, freeBlockIndex + c, rootBlockCoords + udInt3::create(c & 1, c >> 1, 0));

    // preserve spatial position of root in its block
    pQuadTree->rootIndex = (pQuadTree->nodes.pPool[freeBlockIndex + c].slippyPosition == slippyCoords) ? (freeBlockIndex + c) : pQuadTree->rootIndex;
  }

  pQuadTree->completeRerootRequired = false;
  vcQuadTree_InitRootBlock(pQuadTree);
}

void vcQuadTree_ConditionalReroot(vcQuadTree *pQuadTree, const udInt3 &slippyCoords)
{
  while (!pQuadTree->completeRerootRequired && pQuadTree->nodes.pPool[pQuadTree->rootIndex].slippyPosition != slippyCoords)
    vcQuadTree_Reroot(pQuadTree, slippyCoords);

  if (pQuadTree->completeRerootRequired)
    vcQuadTree_CompleteReroot(pQuadTree, slippyCoords);
}

void vcQuadTree_Update(vcQuadTree *pQuadTree, const vcQuadTreeViewInfo &viewInfo)
{
  bool zoneChangeOccurred = pQuadTree->gisSpace.SRID != viewInfo.pSpace->SRID;
  pQuadTree->gisSpace = *viewInfo.pSpace;

  // invalidate so we can detect nodes that need pruning
  for (uint32_t i = 0; i < pQuadTree->nodes.used; ++i)
  {
    vcQuadTreeNode *pNode = &pQuadTree->nodes.pPool[i];
    pNode->rendered = false;
    pNode->touched = false;
    pNode->visible = false;

    if (zoneChangeOccurred)
      vcQuadTree_CalculateNodeBounds(pQuadTree, pNode);
  }

  pQuadTree->slippyCoords = viewInfo.slippyCoords;
  pQuadTree->cameraWorldPosition = viewInfo.cameraPosition;
  pQuadTree->cameraDistanceZeroAltitude = udMag3(pQuadTree->cameraWorldPosition - viewInfo.cameraPositionZeroAltitude);

  pQuadTree->metaData.nodeTouchedCount = 0;
  pQuadTree->metaData.leafNodeCount = 0;
  pQuadTree->metaData.visibleNodeCount = 0;
  pQuadTree->metaData.nodeRenderCount = 0;
  pQuadTree->metaData.maxTreeDepth = udMax(0, (viewInfo.maxVisibleTileLevel - 1) - viewInfo.slippyCoords.z);

  // extract frustum planes
  udDouble4x4 transposedViewProjection = udTranspose(viewInfo.viewProjectionMatrix);
  pQuadTree->frustumPlanes[0] = transposedViewProjection.c[3] + transposedViewProjection.c[0]; // Left
  pQuadTree->frustumPlanes[1] = transposedViewProjection.c[3] - transposedViewProjection.c[0]; // Right
  pQuadTree->frustumPlanes[2] = transposedViewProjection.c[3] + transposedViewProjection.c[1]; // Bottom
  pQuadTree->frustumPlanes[3] = transposedViewProjection.c[3] - transposedViewProjection.c[1]; // Top
  pQuadTree->frustumPlanes[4] = transposedViewProjection.c[3] + transposedViewProjection.c[2]; // Near
  pQuadTree->frustumPlanes[5] = transposedViewProjection.c[3] - transposedViewProjection.c[2]; // Far
  // Normalize the planes
  for (int j = 0; j < 6; ++j)
    pQuadTree->frustumPlanes[j] /= udMag3(pQuadTree->frustumPlanes[j]);

  vcQuadTree_ConditionalReroot(pQuadTree, viewInfo.slippyCoords);

  // TODO: should this go in InitRootBlock()?
  //pQuadTree->nodes.pPool[pQuadTree->rootIndex].morten = udInt2::zero();

  // Must re-check `completeRerootRequired` condition here, because complete re-rooting can fail (finite amount of nodes)
  if (!pQuadTree->completeRerootRequired)
  {
    // validate the entire root block
    uint32_t rootBlockIndex = vcQuadTree_NodeIndexToBlockIndex(pQuadTree->rootIndex);
    for (uint32_t c = 0; c < NodeChildCount; ++c)
      pQuadTree->nodes.pPool[rootBlockIndex + c].touched = true;
    pQuadTree->nodes.pPool[pQuadTree->rootIndex].visible = true;

    vcQuadTree_RecurseGenerateTree(pQuadTree, pQuadTree->rootIndex, 0);
    vcQuadTree_Prune(pQuadTree);
    vcQuadTree_CalculateNeighbours(pQuadTree, &pQuadTree->nodes.pPool[pQuadTree->rootIndex]);
  }
}

bool vcQuadTree_IsBlockUsed(vcQuadTree *pQuadTree, uint32_t blockIndex)
{
  bool isUsed = false;
  for (uint32_t c = 0; c < NodeChildCount && !isUsed; ++c)
    isUsed = pQuadTree->nodes.pPool[blockIndex + c].isUsed;
  return isUsed;
}

bool vcQuadTree_ShouldFreeBlock(vcQuadTree *pQuadTree, uint32_t blockIndex)
{
  for (uint32_t c = 0; c < NodeChildCount; ++c)
  {
    vcQuadTreeNode *pChildNode = &pQuadTree->nodes.pPool[blockIndex + c];
    if (pChildNode->touched || pChildNode->renderInfo.loadStatus == vcNodeRenderInfo::vcTLS_Downloading)
      return false;
  }

  // determine if this block is being used for rendering
  for (uint32_t c = 0; c < NodeChildCount; ++c)
  {
    vcQuadTreeNode *pNode = &pQuadTree->nodes.pPool[blockIndex + c];

    // case #1: its a leaf node that could be being used for rendering by an ancestor
    if (pNode->renderInfo.pTexture)
    {
      uint32_t parentIndex = pNode->parentIndex;
      while (parentIndex != INVALID_NODE_INDEX)
      {
        vcQuadTreeNode *pParentNode = &pQuadTree->nodes.pPool[parentIndex];
        if (pParentNode->touched)
        {
          // We have an ancestor that has no texture
          if (!pParentNode->renderInfo.pTexture)
            return false;

          break;
        }
        else if (pParentNode->renderInfo.pTexture)
          return true;

        parentIndex = pParentNode->parentIndex;
      }
    }
  }

  return true;
}

void vcQuadTree_Prune(vcQuadTree *pQuadTree)
{
  for (uint32_t blockIndex = 0; blockIndex < pQuadTree->nodes.used; blockIndex += NodeChildCount)
  {
    if (pQuadTree->nodes.pPool[blockIndex].isUsed && vcQuadTree_ShouldFreeBlock(pQuadTree, blockIndex))
    {
      vcQuadTreeNode *pNode = &pQuadTree->nodes.pPool[blockIndex];

      // inform parent its children are dead
      if (pNode->parentIndex != INVALID_NODE_INDEX)
      {
        vcQuadTreeNode *pParentNode = &pQuadTree->nodes.pPool[pNode->parentIndex];
        pParentNode->childMask = 0;
        pParentNode->childBlockIndex = INVALID_NODE_INDEX;
      }

      // cleanup block
      for (uint32_t c = 0; c < NodeChildCount; ++c)
      {
        pNode = &pQuadTree->nodes.pPool[blockIndex + c];

        // inform children its parent is dead (this can still happen if e.g. child mid-download)
        uint32_t childBlockIndex = pNode->childBlockIndex;
        if (childBlockIndex != INVALID_NODE_INDEX)
        {
          for (uint32_t c2 = 0; c2 < NodeChildCount; ++c2)
            pQuadTree->nodes.pPool[childBlockIndex + c2].parentIndex = INVALID_NODE_INDEX;
        }

        vcQuadTree_CleanupNode(pNode);
      }
    }
  }

  // trim
  while (pQuadTree->nodes.used > 0 && !vcQuadTree_IsBlockUsed(pQuadTree, pQuadTree->nodes.used - NodeChildCount))
    pQuadTree->nodes.used -= NodeChildCount;
}
