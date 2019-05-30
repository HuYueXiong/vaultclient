#include "vcPOI.h"

#include "vcState.h"
#include "vcRender.h"
#include "vcStrings.h"

#include "vcFenceRenderer.h"

#include "udMath.h"
#include "udFile.h"
#include "udStringUtil.h"

#include "imgui.h"
#include "imgui_ex/vcImGuiSimpleWidgets.h"

vcPOI::vcPOI(vdkProjectNode *pNode, vcState *pProgramState) :
  vcSceneItem(pNode, pProgramState)
{
  m_nameColour = 0xFFFFFFFF;
  m_backColour = 0x7F000000;
  m_namePt = vcLFS_Medium;

  m_showArea = false;
  m_showLength = false;

  memset(&m_line, 0, sizeof(m_line));

  m_line.selectedPoint = -1; // Sentinel for no point selected
  m_line.numPoints = m_pNode->geomCount;

  m_line.colourPrimary = 0xFFFFFFFF;
  m_line.colourSecondary = 0xFFFFFFFF;
  m_line.lineWidth = 1.0;
  m_line.closed = (m_pNode->geomtype == vdkPGT_Polygon);
  m_line.lineStyle = vcRRIM_Arrow;
  m_line.fenceMode = vcRRVM_Fence;

  if (m_line.numPoints > 0)
  {
    m_line.pPoints = udAllocType(udDouble3, m_line.numPoints, udAF_Zero);
    memcpy(m_line.pPoints, m_pNode->pCoordinates, sizeof(udDouble3)*m_line.numPoints);
  }

  if (pProgramState->gis.isProjected)
  {
    m_pPreferredProjection = udAllocType(udGeoZone, 1, udAF_Zero);
    memcpy(m_pPreferredProjection, &pProgramState->gis.zone, sizeof(udGeoZone));
    m_pCurrentProjection = udAllocType(udGeoZone, 1, udAF_Zero);
    memcpy(m_pCurrentProjection, &pProgramState->gis.zone, sizeof(udGeoZone));
    m_line.pPoints[0] = udGeoZone_ToCartesian(*m_pCurrentProjection, m_line.pPoints[0], true);
  }

  m_pLabelText = nullptr;

  m_pLabelInfo = udAllocType(vcLabelInfo, 1, udAF_Zero);
  m_pLabelInfo->pText = m_pNode->pName;
  if (m_line.numPoints > 0)
    m_pLabelInfo->worldPosition = m_line.pPoints[0];
  m_pLabelInfo->textSize = m_namePt;
  m_pLabelInfo->textColourRGBA = vcIGSW_BGRAToRGBAUInt32(m_nameColour);
  m_pLabelInfo->backColourRGBA = vcIGSW_BGRAToRGBAUInt32(m_backColour);

  m_pFence = nullptr;

  OnNodeUpdate();

  m_loadStatus = vcSLS_Loaded;
}

void vcPOI::OnNodeUpdate()
{
  double tempDouble;
  const char *pTemp;

  vdkProjectNode_GetMetadataUint(m_pNode, "nameColour", &m_nameColour, 0xFFFFFFFF);
  vdkProjectNode_GetMetadataUint(m_pNode, "backColour", &m_backColour, 0x7F000000);

  vdkProjectNode_GetMetadataString(m_pNode, "textSize", &pTemp, "Medium");
  if (udStrEquali(pTemp, "x-small") || udStrEquali(pTemp, "small"))
    m_pLabelInfo->textSize = vcLFS_Small;
  else if (udStrEquali(pTemp, "large") || udStrEquali(pTemp, "x-large"))
    m_pLabelInfo->textSize = vcLFS_Large;
  else
    m_pLabelInfo->textSize = vcLFS_Medium;

  vdkProjectNode_GetMetadataBool(m_pNode, "showLength", &m_showLength, true);
  vdkProjectNode_GetMetadataBool(m_pNode, "showArea", &m_showArea, true);
  vdkProjectNode_GetMetadataUint(m_pNode, "lineColourPrimary", &m_line.colourPrimary, 0xFFFFFFFF);
  vdkProjectNode_GetMetadataUint(m_pNode, "lineColourSecondary", &m_line.colourSecondary, 0xFFFFFFFF);
  vdkProjectNode_GetMetadataDouble(m_pNode, "lineWidth", (double*)&tempDouble, 1.0);
  m_line.lineWidth = tempDouble;

  vdkProjectNode_GetMetadataString(m_pNode, "lineStyle", &pTemp, vcFRIMStrings[0]);
  int i = 0;
  for (; i < vcRRIM_Count; ++i)
    if (udStrEquali(pTemp, vcFRIMStrings[i]))
      break;
  m_line.lineStyle = (vcFenceRendererImageMode)i;

  vdkProjectNode_GetMetadataString(m_pNode, "lineMode", &pTemp, vcFRVMStrings[0]);
  for (i = 0; i < vcRRVM_Count; ++i)
    if (udStrEquali(pTemp, vcFRVMStrings[i]))
      break;
  m_line.fenceMode = (vcFenceRendererVisualMode)i;

  UpdatePoints();
}

void vcPOI::AddToScene(vcState *pProgramState, vcRenderData *pRenderData)
{
  // if POI is invisible or if it exceeds maximum visible POI distance
  if (!m_visible || (pProgramState->settings.camera.cameraMode != vcCM_OrthoMap && udMag3(m_pLabelInfo->worldPosition - pProgramState->pCamera->position) > pProgramState->settings.presentation.POIFadeDistance))
    return;

  //TODO: Super hacky- remove later
  m_pProject = pProgramState->sceneExplorer.pProject;

  if (m_pFence != nullptr)
    pRenderData->fences.PushBack(m_pFence);

  if (m_pLabelInfo != nullptr)
  {
    if (m_showLength || m_showArea)
      m_pLabelInfo->pText = m_pLabelText;
    else
      m_pLabelInfo->pText = m_pNode->pName;

    pRenderData->labels.PushBack(m_pLabelInfo);
  }
}

void vcPOI::ApplyDelta(vcState * /*pProgramState*/, const udDouble4x4 &delta)
{
  //bool isGeolocated = m_pZone != nullptr;
  if (m_line.selectedPoint == -1) // We need to update all the points
  {
    for (int i = 0; i < m_line.numPoints; ++i)
      m_line.pPoints[i] = (delta * udDouble4x4::translation(m_line.pPoints[i])).axis.t.toVector3();
  }
  else
  {
    m_line.pPoints[m_line.selectedPoint] = (delta * udDouble4x4::translation(m_line.pPoints[m_line.selectedPoint])).axis.t.toVector3();
  }

  UpdatePoints();
  UpdateProjectGeometry();
}

void vcPOI::UpdatePoints()
{
  // Calculate length, area and label position
  m_calculatedLength = 0;
  m_calculatedArea = 0;
  udDouble3 averagePosition = udDouble3::zero();

  int j = ((m_line.numPoints == 0) ? 0 : m_line.numPoints - 1);

  for (int i = 0; i < m_line.numPoints; i++)
  {
    if (m_showArea && m_line.closed && m_line.numPoints > 2) // Area requires at least 2 points
      m_calculatedArea = m_calculatedArea + (m_line.pPoints[j].x + m_line.pPoints[i].x) * (m_line.pPoints[j].y - m_line.pPoints[i].y);

    if (m_line.closed || i > 0) // Calculate length
      m_calculatedLength += udMag3(m_line.pPoints[j] - m_line.pPoints[i]);

    averagePosition += m_line.pPoints[i];

    j = i;
  }

  m_calculatedArea = udAbs(m_calculatedArea) / 2;
  m_pLabelInfo->worldPosition = averagePosition / m_line.numPoints;

  if (m_showArea && m_showLength)
    udSprintf(&m_pLabelText, "%s\n%s: %.3f\n%s: %.3f", m_pNode->pName, vcString::Get("scenePOILineLength"), m_calculatedLength, vcString::Get("scenePOIArea"), m_calculatedArea);
  else if (m_showLength)
    udSprintf(&m_pLabelText, "%s\n%s: %.3f", m_pNode->pName, vcString::Get("scenePOILineLength"), m_calculatedLength);
  else if (m_showArea)
    udSprintf(&m_pLabelText, "%s\n%s: %.3f", m_pNode->pName, vcString::Get("scenePOIArea"), m_calculatedArea);

  // update the fence renderer as well
  if (m_line.numPoints > 1)
  {
    if (m_pFence == nullptr)
      vcFenceRenderer_Create(&m_pFence);

    vcFenceRendererConfig config;
    config.visualMode = m_line.fenceMode;
    config.imageMode = m_line.lineStyle;
    config.bottomColour = vcIGSW_BGRAToImGui(m_line.colourSecondary);
    config.topColour = vcIGSW_BGRAToImGui(m_line.colourPrimary);
    config.ribbonWidth = m_line.lineWidth;
    config.textureScrollSpeed = 1.f;
    config.textureRepeatScale = 1.f;

    vcFenceRenderer_SetConfig(m_pFence, config);

    vcFenceRenderer_ClearPoints(m_pFence);
    vcFenceRenderer_AddPoints(m_pFence, m_line.pPoints, m_line.numPoints, m_line.closed);
  }
}

void vcPOI::UpdateProjectGeometry()
{
  if (m_pCurrentProjection == nullptr)
    return; // Can't update if we aren't sure what zone we're in currently

  udDouble3 *pGeom = udAllocType(udDouble3, m_line.numPoints, udAF_Zero);

  // Change all points in the POI to the new projection
  for (int i = 0; i < m_line.numPoints; ++i)
    pGeom[i] = udGeoZone_ToLatLong(*m_pCurrentProjection, m_line.pPoints[i], true);

  if (m_line.closed)
    vdkProjectNode_SetGeometry(m_pProject, m_pNode, vdkPGT_Polygon, m_line.numPoints, (double*)pGeom);
  else
    vdkProjectNode_SetGeometry(m_pProject, m_pNode, vdkPGT_LineString, m_line.numPoints, (double*)pGeom);
}

void vcPOI::HandleImGui(vcState *pProgramState, size_t *pItemID)
{
  if (vcIGSW_ColorPickerU32(udTempStr("%s##POIColour%zu", vcString::Get("scenePOILabelColour"), *pItemID), &m_nameColour, ImGuiColorEditFlags_None))
  {
    m_pLabelInfo->textColourRGBA = vcIGSW_BGRAToRGBAUInt32(m_nameColour);
    vdkProjectNode_SetMetadataUint(m_pNode, "nameColour", m_nameColour);
  }

  if (vcIGSW_ColorPickerU32(udTempStr("%s##POIBackColour%zu", vcString::Get("scenePOILabelBackgroundColour"), *pItemID), &m_backColour, ImGuiColorEditFlags_None))
  {
    m_pLabelInfo->backColourRGBA = vcIGSW_BGRAToRGBAUInt32(m_backColour);
    vdkProjectNode_SetMetadataUint(m_pNode, "backColour", m_backColour);
  }

  const char *labelSizeOptions[] = { vcString::Get("scenePOILabelSizeNormal"), vcString::Get("scenePOILabelSizeSmall"), vcString::Get("scenePOILabelSizeLarge") };
  if (ImGui::Combo(udTempStr("%s##POILabelSize%zu", vcString::Get("scenePOILabelSize"), *pItemID), (int*)&m_pLabelInfo->textSize, labelSizeOptions, (int)udLengthOf(labelSizeOptions)))
  {
    UpdatePoints();
    const char *pTemp;
    switch (m_pLabelInfo->textSize)
    {
    case vcLFS_Small:
      pTemp = "small";
      break;
    case vcLFS_Medium:
      pTemp = "medium";
      break;
    case vcLFS_Large:
      pTemp = "large";
      break;
    }
    vdkProjectNode_SetMetadataString(m_pNode, "textSize", pTemp);
  }

  if (m_line.numPoints > 1)
  {
    if (!pProgramState->cameraInput.flyThroughActive)
    {
      if (ImGui::Button(vcString::Get("scenePOIPerformFlyThrough")))
      {
        pProgramState->cameraInput.inputState = vcCIS_FlyingThrough;
        pProgramState->cameraInput.flyThroughActive = false; // set false to activate MoveTo, will be set true in vcCamera before next frame
        pProgramState->cameraInput.pObjectInfo = &m_line;
      }
    }
    else
    {
      if (ImGui::Button(vcString::Get("scenePOICancelFlyThrough")))
      {
        pProgramState->cameraInput.inputState = vcCIS_None;
        pProgramState->cameraInput.flyThroughActive = false;
        pProgramState->cameraInput.flyThroughPoint = 0;
        pProgramState->cameraInput.pObjectInfo = nullptr;
      }
    }

    if (ImGui::SliderInt(vcString::Get("scenePOISelectedPoint"), &m_line.selectedPoint, -1, m_line.numPoints - 1))
      m_line.selectedPoint = udClamp(m_line.selectedPoint, -1, m_line.numPoints - 1);

    if (m_line.selectedPoint != -1)
    {
      if (ImGui::InputScalarN(udTempStr("%s##POIPointPos%zu", vcString::Get("scenePOIPointPosition"), *pItemID), ImGuiDataType_Double, &m_line.pPoints[m_line.selectedPoint].x, 3))
        UpdatePoints();

      if (ImGui::Button(vcString::Get("scenePOIRemovePoint")))
        RemovePoint(m_line.selectedPoint);
    }

    if (ImGui::TreeNode("%s##POILineSettings%zu", vcString::Get("scenePOILineSettings"), *pItemID))
    {
      if (ImGui::Checkbox(udTempStr("%s##POIShowLength%zu", vcString::Get("scenePOILineShowLength"), *pItemID), &m_showLength))
      {
        UpdatePoints();
        vdkProjectNode_SetMetadataBool(m_pNode, "showLength", m_showLength);
      }

      if (ImGui::Checkbox(udTempStr("%s##POIShowArea%zu", vcString::Get("scenePOILineShowArea"), *pItemID), &m_showArea))
      {
        UpdatePoints();
        vdkProjectNode_SetMetadataBool(m_pNode, "showArea", m_showArea);
      }

      if (ImGui::Checkbox(udTempStr("%s##POILineClosed%zu", vcString::Get("scenePOILineClosed"), *pItemID), &m_line.closed))
        UpdatePoints();

      if (vcIGSW_ColorPickerU32(udTempStr("%s##POILineColourPrimary%zu", vcString::Get("scenePOILineColour1"), *pItemID), &m_line.colourPrimary, ImGuiColorEditFlags_None))
      {
        UpdatePoints();
        vdkProjectNode_SetMetadataUint(m_pNode, "lineColourPrimary", m_line.colourPrimary);
      }

      if (vcIGSW_ColorPickerU32(udTempStr("%s##POILineColourSecondary%zu", vcString::Get("scenePOILineColour2"), *pItemID), &m_line.colourSecondary, ImGuiColorEditFlags_None))
      {
        UpdatePoints();
        vdkProjectNode_SetMetadataUint(m_pNode, "lineColourSecondary", m_line.colourSecondary);
      }

      if (ImGui::SliderFloat(udTempStr("%s##POILineWidth%zu", vcString::Get("scenePOILineWidth"), *pItemID), &m_line.lineWidth, 0.01f, 1000.f, "%.2f", 3.f))
      {
        UpdatePoints();
        vdkProjectNode_SetMetadataDouble(m_pNode, "lineWidth", (double)m_line.lineWidth);
      }

      const char *lineOptions[] = { vcString::Get("scenePOILineStyleArrow"), vcString::Get("scenePOILineStyleGlow"), vcString::Get("scenePOILineStyleSolid") };
      if (ImGui::Combo(udTempStr("%s##POILineColourSecondary%zu", vcString::Get("scenePOILineStyle"), *pItemID), (int *)&m_line.lineStyle, lineOptions, (int)udLengthOf(lineOptions)))
      {
        UpdatePoints();
        vdkProjectNode_SetMetadataString(m_pNode, "lineStyle", vcFRIMStrings[m_line.lineStyle]);
      }

      const char *fenceOptions[] = { vcString::Get("scenePOILineOrientationVert"), vcString::Get("scenePOILineOrientationHorz") };
      if (ImGui::Combo(udTempStr("%s##POIFenceStyle%zu", vcString::Get("scenePOILineOrientation"), *pItemID), (int *)&m_line.fenceMode, fenceOptions, (int)udLengthOf(fenceOptions)))
      {
        UpdatePoints();
        vdkProjectNode_SetMetadataString(m_pNode, "lineMode", vcFRVMStrings[m_line.fenceMode]);
      }

      ImGui::TreePop();
    }
  }

  // Handle hyperlinks
  const char *pHyperlink = m_metadata.Get("hyperlink").AsString();
  if (pHyperlink != nullptr)
  {
    ImGui::TextWrapped("%s: %s", vcString::Get("scenePOILabelHyperlink"), pHyperlink);
    if (udStrEndsWithi(pHyperlink, ".png") || udStrEndsWithi(pHyperlink, ".jpg"))
    {
      ImGui::SameLine();
      if (ImGui::Button(vcString::Get("scenePOILabelOpenHyperlink")))
        pProgramState->pLoadImage = udStrdup(pHyperlink);
    }
  }
}

void vcPOI::AddPoint(const udDouble3 &position)
{
  udDouble3 *pNewPoints = udAllocType(udDouble3, m_line.numPoints + 1, udAF_Zero);

  memcpy(pNewPoints, m_line.pPoints, sizeof(udDouble3) * m_line.numPoints);
  pNewPoints[m_line.numPoints] = (m_pCurrentProjection == nullptr) ? position : udGeoZone_ToCartesian(*m_pCurrentProjection, position, true);
  udFree(m_line.pPoints);
  m_line.pPoints = pNewPoints;

  ++m_line.numPoints;
  UpdatePoints();
  UpdateProjectGeometry();
}

void vcPOI::RemovePoint(int index)
{
  if (index < (m_line.numPoints - 1))
    memmove(m_line.pPoints + index, m_line.pPoints + index + 1, sizeof(udDouble3) * (m_line.numPoints - index - 1));

  --m_line.numPoints;
  UpdatePoints();
  UpdateProjectGeometry();
}

void vcPOI::ChangeProjection(const udGeoZone &newZone)
{
  // If new zone is different and points fit within this zone's bounds, assign it to this new zone
  if (m_pCurrentProjection == nullptr || m_pCurrentProjection->srid != newZone.srid)
  {
    bool withinBounds = true;
    for (int i = 0; i < m_line.numPoints; ++i)
    {
      udDouble3 *pCoord = (udDouble3*)m_pNode->pCoordinates + (3 * i);
      if (pCoord->y < newZone.latLongBoundMin.x || pCoord->y > newZone.latLongBoundMax.x || pCoord->x < newZone.latLongBoundMin.y || pCoord->x > newZone.latLongBoundMax.y)
      {
        withinBounds = false;
        break;
      }
    }
    if (withinBounds)
    {
      m_pPreferredProjection = udAllocType(udGeoZone, 1, udAF_Zero);
      memcpy(m_pPreferredProjection, &newZone, sizeof(udGeoZone));
      m_pCurrentProjection = udAllocType(udGeoZone, 1, udAF_Zero);
      memcpy(m_pCurrentProjection, &newZone, sizeof(udGeoZone));

      // Change coords from longlat to cartesian within the new zone
      for (int i = 0; i < m_line.numPoints; ++i)
        ((udDouble3*)m_line.pPoints)[i] = udGeoZone_ToCartesian(newZone, *(udDouble3*)m_pNode->pCoordinates, true);
    }
  }

  UpdatePoints();

  // Call the parent version - m_pCurrentProjection is updated in here
  vcSceneItem::ChangeProjection(newZone);
}

void vcPOI::Cleanup(vcState * /*pProgramState*/)
{
  udFree(m_line.pPoints);
  udFree(m_pLabelText);

  vcFenceRenderer_Destroy(&m_pFence);
  udFree(m_pLabelInfo);
}

void vcPOI::SetCameraPosition(vcState *pProgramState)
{
  pProgramState->pCamera->position = m_pLabelInfo->worldPosition;
}

udDouble4x4 vcPOI::GetWorldSpaceMatrix()
{
  if (m_line.selectedPoint == -1)
    return udDouble4x4::translation(m_pLabelInfo->worldPosition);
  else
    return udDouble4x4::translation(m_line.pPoints[m_line.selectedPoint]);
}
