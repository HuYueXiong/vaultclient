#include "vcWaterNode.h"

#include "vcStrings.h"
#include "vcState.h"

#include "vcWaterRenderer.h"
#include "vcRender.h"

#include "imgui.h"
#include "imgui_ex/vcImGuiSimpleWidgets.h"

vcWater::vcWater(vdkProject *pProject, vdkProjectNode *pNode, vcState *pProgramState) :
  vcSceneItem(pProject, pNode, pProgramState)
{
  m_pWaterRenderer = nullptr;

  vcWaterRenderer_Create(&m_pWaterRenderer);
  m_loadStatus = vcSLS_Loaded;

  OnNodeUpdate(pProgramState);
}

void vcWater::OnNodeUpdate(vcState *pProgramState)
{
  ChangeProjection(pProgramState->gis.zone);
}

void vcWater::AddToScene(vcState * /*pProgramState*/, vcRenderData *pRenderData)
{
  pRenderData->waterVolumes.PushBack(m_pWaterRenderer);
}

void vcWater::ApplyDelta(vcState * /*pProgramState*/, const udDouble4x4 & /*delta*/)
{

}

void vcWater::HandleImGui(vcState * /*pProgramState*/, size_t * /*pItemID*/)
{
  //TODO: Water Settings
}

void vcWater::Cleanup(vcState * /*pProgramState*/)
{
  vcWaterRenderer_Destroy(&m_pWaterRenderer);
}

void vcWater::ChangeProjection(const udGeoZone &newZone)
{
  if (m_pWaterRenderer == nullptr)
    return;

  vcWaterRenderer_ClearAllVolumes(m_pWaterRenderer);

  m_pivot = udGeoZone_LatLongToCartesian(newZone, udDouble3::create(m_pNode->pCoordinates[0], m_pNode->pCoordinates[1], m_pNode->pCoordinates[2]));

  udDouble2 *pPoints = udAllocType(udDouble2, m_pNode->geomCount, udAF_None);

  for (int i = 0; i < m_pNode->geomCount; ++i)
  {
    pPoints[i] = udGeoZone_LatLongToCartesian(newZone, udDouble3::create(m_pNode->pCoordinates[i * 3 + 0], m_pNode->pCoordinates[i * 3 + 1], m_pNode->pCoordinates[i * 3 + 2])).toVector2();
  }

  if (m_pNode->pFirstChild != nullptr)
  {
    vdkProjectNode *pIslandNode = m_pNode->pFirstChild;

    do
    {
      udDouble2 *pIslandPoints = udAllocType(udDouble2, pIslandNode->geomCount, udAF_None);

      for (int i = 0; i < pIslandNode->geomCount; ++i)
      {
        pIslandPoints[i] = udGeoZone_LatLongToCartesian(newZone, udDouble3::create(pIslandNode->pCoordinates[i * 3 + 0], pIslandNode->pCoordinates[i * 3 + 1], pIslandNode->pCoordinates[i * 3 + 2])).toVector2();
      }

      udFree(pIslandPoints);
      pIslandNode = pIslandNode->pNextSibling;
    } while (pIslandNode != nullptr);
  }

  vcWaterRenderer_AddVolume(m_pWaterRenderer, pPoints, m_pNode->geomCount);
  udFree(pPoints);
}

udDouble3 vcWater::GetLocalSpacePivot()
{
  return m_pivot;
}
