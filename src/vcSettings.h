#ifndef vcSettings_h__
#define vcSettings_h__

#include "SDL2/SDL.h"

#include "udPlatform.h"
#include "udPlatformUtil.h"
#include "udMath.h"
#include "udChunkedArray.h"
#include "udUUID.h"

#include "vcCamera.h"
#include "gl/vcGLState.h"
#include "gl/vcTexture.h"

#include "imgui.h"

enum vcMapTileBlendMode
{
  vcMTBM_Hybrid,
  vcMTBM_Overlay,
  vcMTBM_Underlay,

  vcMTBM_Count
};

enum vcVisualizatationMode
{
  vcVM_Default,
  vcVM_Colour,
  vcVM_Intensity,
  vcVM_Classification,
  vcVM_Displacement
};

enum vcAnchorStyle
{
  vcAS_None,
  vcAS_Orbit,
  vcAS_Compass,

  vcAS_Count
};

enum
{
  vcMaxPathLength = 512,
  vcMetadataMaxLength = 256,
};

enum vcPresentationMode
{
  vcPM_Hide,
  vcPM_Show,
  vcPM_Responsive
};

enum vcWindowLayout
{
  vcWL_SceneLeft,
  vcWL_SceneRight
};

enum vcSettingsUIRegions
{
  vcSR_Appearance,
  vcSR_Inputs,
  vcSR_Maps,
  vcSR_Visualisations,
  vcSR_KeyBindings,
  vcSR_ConvertDefaults,
  vcSR_Screenshot,
  vcSR_Connection,

  vcSR_ReleaseNotes,
  vcSR_Update,
  vcSR_About,
  vcSR_Count
};

enum vcSettingCategory
{
  vcSC_Appearance,
  vcSC_InputControls,
  vcSC_MapsElevation,
  vcSC_Visualization,
  vcSC_Convert,
  vcSC_Languages,
  vcSC_Bindings,
  vcSC_Screenshot,
  vcSC_Connection,

  vcSC_All
};

enum vcTileRendererMapQuality
{
  vcTRMQ_Low,
  vcTRMQ_Medium,
  vcTRMQ_High,

  vcTRMQ_Total,
};

enum vcTileRendererFlags
{
  vcTRF_None = 0,
  vcTRF_OnlyRequestVisibleTiles = 0x1,
};

struct vcLanguageOption
{
  // Arbitrary limits
  char languageName[32];
  char filename[16];
};

struct vcVisualizationSettings
{
  vcVisualizatationMode mode;

  int minIntensity;
  int maxIntensity;

  udFloat2 displacement;

  bool useCustomClassificationColours;
  bool customClassificationToggles[256];
  uint32_t customClassificationColors[256];
  const char *customClassificationColorLabels[256];
};

struct vcSettings
{
  bool noLocalStorage; //If set to true; cannot save or load from local storage
  const char *pSaveFilePath;
  char cacheAssetPath[vcMaxPathLength];

  bool onScreenControls;

  struct
  {
    int styleIndex;
    bool showCameraInfo;
    bool showProjectionInfo;
    bool showDiagnosticInfo;
    bool showEuclideonLogo;
    bool showAdvancedGIS;
    bool loginRenderLicense;

    bool showSkybox;
    udFloat4 skyboxColour;

    float saturation;

    vcAnchorStyle mouseAnchor;
    bool showCompass;
    float POIFadeDistance;
    float imageRescaleDistance;
    bool limitFPSInBackground;

    int pointMode;

    vcWindowLayout layout;
    int sceneExplorerSize;
    float convertLeftPanelPercentage;
    bool columnSizeCorrect; // Not saved, this is updated when the columns have been set to the correct size
  } presentation;

  struct
  {
    int xpos;
    int ypos;
    int width;
    int height;
    bool maximized;
    bool touchscreenFriendly;

    bool isFullscreen;
    bool useNativeUI;

    char languageCode[16];
  } window;

  struct
  {
    bool rememberServer;
    char serverURL[vcMaxPathLength];

    bool rememberUsername;
    char username[vcMaxPathLength];

    char proxy[vcMaxPathLength];
    char proxyTestURL[vcMaxPathLength];
    char autoDetectProxyURL[vcMaxPathLength];
    char userAgent[vcMaxPathLength];
    bool ignoreCertificateVerification;

    bool requiresProxyAuth;

    bool rememberProxyUsername;
    char proxyUsername[vcMaxPathLength];
    char proxyPassword[vcMaxPathLength];

    // These are not saved but are required testing information
    udInterlockedBool testing;
    bool tested;
    int testStatus;
  } loginInfo;

  vcVisualizationSettings visualization;

  struct
  {
    struct
    {
      bool enable;
      int width;
      float threshold;
      udFloat4 colour;
    } edgeOutlines;

    struct
    {
      bool enable;
      udFloat4 minColour;
      udFloat4 maxColour;
      float startHeight;
      float endHeight;
    } colourByHeight;

    struct
    {
      bool enable;
      udFloat4 colour;
      float startDepth;
      float endDepth;
    } colourByDepth;

    struct
    {
      bool enable;
      udFloat4 colour;
      float distances;
      float bandHeight;
      float rainbowRepeat;
      float rainbowIntensity;
    } contours;
  } postVisualization;

  vcCameraSettings camera;

  struct
  {
    bool mapEnabled;

    vcTileRendererMapQuality mapQuality;
    vcTileRendererFlags mapOptions;

    float mapHeight;
    char tileServerAddress[vcMaxPathLength];
    udUUID tileServerAddressUUID;
    char tileServerExtension[4];

    vcMapTileBlendMode blendMode;
    float transparency;

    bool mouseInteracts;
  } maptiles;

  struct
  {
    bool enable;
    udFloat4 colour;
    float thickness;
  } objectHighlighting;

  struct
  {
    char tempDirectory[vcMaxPathLength];
    struct
    {
      bool isDirty;
      char filename[vcMaxPathLength];
      vcTexture *pTexture;
      int width;
      int height;
    } watermark;

    char author[vcMetadataMaxLength];
    char comment[vcMetadataMaxLength];
    char copyright[vcMetadataMaxLength];
    char license[vcMetadataMaxLength];
  } convertdefaults;

  struct
  {
    bool taking;
    uint8_t* pPixels;

    bool viewShot;

    int resolutionIndex; // Not saved
    udUInt2 resolution;

    char outputPath[vcMaxPathLength];
  } screenshot;

  // These are experimental features that will eventually be removed or moved to another setting.
  // They will mostly be exposed via the System->Experiments menu to hide them away from most users
  struct
  {
    //No current experiments
  } experimental;

  udChunkedArray<vcLanguageOption> languageOptions;

  vcPresentationMode responsiveUI;
  int hideIntervalSeconds;

  struct
  {
    bool enable;
    int range;
  } mouseSnap;
};

// Settings Limits (vcSL prefix)
const double vcSL_GlobalLimit = 40000000.0; // these limits can apply to any sliders/inputs related to position or distance
const double vcSL_GlobalLimitSmall = 1000000.0;
const float vcSL_GlobalLimitf = (float)vcSL_GlobalLimit;
const float vcSL_GlobalLimitSmallf = (float)vcSL_GlobalLimitSmall;

const float vcSL_CameraMinMoveSpeed = 0.5f;
const float vcSL_CameraMaxMoveSpeed = 500.f;
const float vcSL_CameraFieldOfViewMin = 5;
const float vcSL_CameraFieldOfViewMax = 100;

const float vcSL_OSCPixelRatio = 100.f;

const float vcSL_POIFaderMin = 0.f;
const float vcSL_POIFaderMax = 1000000.f;

const float vcSL_ImageRescaleMin = 0.f;
const float vcSL_ImageRescaleMax = 1000000.f;

const float vcSL_MapHeightMin = -1000.f;
const float vcSL_MapHeightMax = 1000.f;

const float vcSL_OpacityMin = 0.f;
const float vcSL_OpacityMax = 1.f;

const float vcSL_IntensityMin = 0.f;
const float vcSL_IntensityMax = 65535.f;

const int vcSL_EdgeHighlightMin = 1;
const int vcSL_EdgeHighlightMax = 10;
const float vcSL_EdgeHighlightThresholdMin = 0.001f;
const float vcSL_EdgeHighlightThresholdMax = 10.f;

const float vcSL_ColourByHeightMin = 0.f;
const float vcSL_ColourByHeightMax = 1000.f;
const float vcSL_ColourByDepthMin = 0.f;
const float vcSL_ColourByDepthMax = 1000.f;
const float vcSL_ContourDistanceMin = 0.f;
const float vcSL_ContourDistanceMax = 100.f;
const float vcSL_ContourBandHeightMin = 0.f;
const float vcSL_ContourBandHeightMax = 10.f;

const int vcSL_MouseSnapRangeMax = 50;

// Settings Functions
bool vcSettings_Load(vcSettings *pSettings, bool forceReset = false, vcSettingCategory group = vcSC_All);
bool vcSettings_Save(vcSettings *pSettings);

void vcSettings_Cleanup(vcSettings *pSettings);

// Uses udTempStr internally.
const char *vcSettings_GetAssetPath(const char *pFilename);

// Provides a handler for "asset://" files
udResult vcSettings_RegisterAssetFileHandler();

// Various settings helpers
udResult vcSettings_UpdateLanguageOptions(vcSettings *pSettings);

#endif // !vcSettings_h__
