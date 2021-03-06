#ifndef vcSettingsUI_h__
#define vcSettingsUI_h__

#include "imgui.h"

struct vcState;
struct vcVisualizationSettings;

void vcSettingsUI_Show(vcState *pProgramState);
bool vcSettingsUI_LangCombo(vcState *pProgramState);
void vcSettingsUI_VisualizationSettings(vcState *pProgramState, vcVisualizationSettings *pVisualizationSettings);

#endif
