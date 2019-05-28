#include "vcI3S.h"

#include "vcState.h"
#include "vcStrings.h"
#include "vcRender.h"

#include "udPlatform.h"

#include "imgui.h"
#include "imgui_ex/vcImGuiSimpleWidgets.h"

vcI3S::vcI3S(vdkProjectNode *pNode, vcState *pProgramState) :
  vcSceneItem(pNode, pProgramState),
  m_pSceneRenderer(nullptr)
{
  if (vcSceneLayerRenderer_Create(&m_pSceneRenderer, pProgramState->pWorkerPool, pNode->pURI) == udR_Success)
  {
    m_loadStatus = vcSLS_Loaded;

    udGeoZone *pInternalZone = vcSceneLayer_GetPreferredZone(m_pSceneRenderer->pSceneLayer);
    if (pInternalZone != nullptr)
    {
      m_pPreferredProjection = udAllocType(udGeoZone, 1, udAF_Zero);
      m_pCurrentProjection = udAllocType(udGeoZone, 1, udAF_Zero);
      memcpy(m_pPreferredProjection, pInternalZone, sizeof(udGeoZone));
      memcpy(m_pCurrentProjection, pInternalZone, sizeof(udGeoZone));
    }
  }
  else
  {
    m_loadStatus = vcSLS_Failed;
  }
}

vcI3S::~vcI3S()
{
  vcSceneLayerRenderer_Destroy(&m_pSceneRenderer);
}

void vcI3S::AddToScene(vcState * /*pProgramState*/, vcRenderData *pRenderData)
{
  if (!m_visible || m_pSceneRenderer == nullptr) // Nothing to render
    return;

  pRenderData->sceneLayers.PushBack(m_pSceneRenderer);
}

void vcI3S::ApplyDelta(vcState * /*pProgramState*/, const udDouble4x4 & /*delta*/)
{
  //I3S doesn't support being moved
}

void vcI3S::HandleImGui(vcState * /*pProgramState*/, size_t * /*pItemID*/)
{
  ImGui::TextUnformatted(vcString::Get("sceneExplorerI3S"));
}

void vcI3S::Cleanup(vcState * /*pProgramState*/)
{
  // Do stuff
}

void vcI3S::ChangeProjection(const udGeoZone &newZone)
{
  //TODO: Support this

  // Call the base version
  vcSceneItem::ChangeProjection(newZone);
}

udDouble3 vcI3S::GetLocalSpacePivot()
{
  if (m_pSceneRenderer == nullptr)
    return udDouble3::zero();

  return vcSceneLayer_GetCenter(m_pSceneRenderer->pSceneLayer);
}