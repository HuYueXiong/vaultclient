#ifndef vcViewpoint_h__
#define vcViewpoint_h__

#include "vcSceneItem.h"
#include "vcCamera.h"
#include "vdkRenderContext.h"

struct vdkPointCloud;
struct vcState;

class vcViewpoint : public vcSceneItem
{
private:
  udDouble3 m_CameraPosition;
  udDouble3 m_CameraRotation;

public:
  vcViewpoint(vdkProject *pProject, vdkProjectNode *pNode, vcState *pProgramState);
  ~vcViewpoint() {};

  void OnNodeUpdate(vcState *pProgramState);

  void AddToScene(vcState *pProgramState, vcRenderData *pRenderData);
  void ApplyDelta(vcState *pProgramState, const udDouble4x4 &delta);
  void HandleImGui(vcState *pProgramState, size_t *pItemID);
  void Cleanup(vcState *pProgramState);
  void ChangeProjection(const udGeoZone &newZone);

  void SetCameraPosition(vcState *pProgramState);
  udDouble4x4 GetWorldSpaceMatrix();
};

#endif //vcViewpoint_h__
