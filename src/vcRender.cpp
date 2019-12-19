#include "vcRender.h"

#include "gl/vcRenderShaders.h"
#include "gl/vcFramebuffer.h"
#include "gl/vcShader.h"
#include "gl/vcGLState.h"
#include "vcFenceRenderer.h"
#include "vcWaterRenderer.h"
#include "vcTileRenderer.h"
#include "vcCompass.h"
#include "vcState.h"
#include "vcVoxelShaders.h"

#include "vcInternalModels.h"
#include "vcSceneLayerRenderer.h"
#include "vcCamera.h"

#include "stb_image.h"
#include <vector>

#include "atmosphere/model.h"

atmosphere::Model *pModel = nullptr;

enum
{
  // directX framebuffer can only be certain increments
  vcRender_SceneSizeIncrement = 32,

  // certain effects don't need to be at 100% resolution (e.g. outline). 0 is highest quality
  vcRender_OutlineEffectDownscale = 1,

  // number of buffers for primary rendering passes
  vcRender_RenderBufferCount = 2,
};

// Temp hard-coded view shed properties
static const int ViewShedMapCount = 3;
static const udUInt2 ViewShedMapRes = udUInt2::create(640 * ViewShedMapCount, 1920);

struct vcViewShedRenderContext
{
  // re-use this memory
  float *pDepthBuffer;
  vcTexture *pUDDepthTexture;

  vcTexture *pDepthTex;
  vcTexture *pDummyColour;
  vcFramebuffer *pFramebuffer;

  vdkRenderView *pRenderView;
};

struct vcUDRenderContext
{
  vdkRenderContext *pRenderer;
  vdkRenderView *pRenderView;
  uint32_t *pColorBuffer;
  float *pDepthBuffer;

  vcFramebuffer *pFramebuffer;
  vcTexture *pColourTex;
  vcTexture *pDepthTex;

  struct
  {
    vcShader *pProgram;
    vcShaderSampler *uniform_texture;
    vcShaderSampler *uniform_depth;
  } presentShader;

  struct
  {
    vcShader *pProgram;
    vcShaderConstantBuffer *uniform_params;
    vcShaderSampler *uniform_texture;
    vcShaderSampler *uniform_depth;

    struct
    {
      udFloat4 id;
    } params;

  } splatIdShader;
};

struct vcRenderContext
{
  vcViewShedRenderContext viewShedRenderingContext;

  udUInt2 sceneResolution;
  udUInt2 originalSceneResolution;

  uint32_t activeRenderTarget;

  vcFramebuffer *pFramebuffer[vcRender_RenderBufferCount];
  vcTexture *pTexture[vcRender_RenderBufferCount];
  vcTexture *pDepthTexture[vcRender_RenderBufferCount];

  vcFramebuffer *pAuxiliaryFramebuffers[2];
  vcTexture *pAuxiliaryTextures[2];
  udUInt2 effectResolution;

  struct
  {
    vcShader *pProgram;
    vcShaderSampler *uniform_texture;
    vcShaderSampler *uniform_depth;
    vcShaderConstantBuffer *uniform_vertParams;
    vcShaderConstantBuffer *uniform_fragParams;

    struct
    {
      udFloat4 screenParams;  // sampleStepX, sampleStepSizeY, near plane, far plane
      udFloat4x4 inverseViewProjection;
      udFloat4x4 inverseProjection;

      // outlining
      udFloat4 outlineColour;
      udFloat4 outlineParams;   // outlineWidth, threshold, (unused), (unused)

      // colour by height
      udFloat4 colourizeHeightColourMin;
      udFloat4 colourizeHeightColourMax;
      udFloat4 colourizeHeightParams; // min world height, max world height, (unused), (unused)

      // colour by depth
      udFloat4 colourizeDepthColour;
      udFloat4 colourizeDepthParams; // min depth, max depth, (unused), (unused)

      // contours
      udFloat4 contourColour;
      udFloat4 contourParams; // contour distance, contour band height, contour rainbow repeat rate, contour rainbow factoring
    } fragParams;

    struct
    {
      udFloat4 outlineStepSize; // outlineStepSize.xy (in uv space), (unused), (unused)
    } vertParams;

  } visualizationShader;

  struct
  {
    vcShader *pProgram;
    vcShaderSampler *uniform_texture;
    vcShaderSampler *uniform_depth;
    vcShaderConstantBuffer *uniform_params;

    struct
    {
      udFloat4 screenParams;  // sampleStepX, sampleStepSizeY, near plane, far plane
    } params;

  } fxaaShader;

  struct
  {
    vcShader *pProgram;
    vcShaderSampler *uniform_depth;
    vcShaderSampler *uniform_shadowMapAtlas;
    vcShaderConstantBuffer *uniform_params;

    struct
    {
      udFloat4x4 shadowMapVP[ViewShedMapCount];
      udFloat4x4 inverseProjection;
      udFloat4 visibleColour;
      udFloat4 notVisibleColour;
      udFloat4 nearFarPlane; // .zw unused
    } params;

  } shadowShader;

  vcUDRenderContext udRenderContext;
  vcFenceRenderer *pDiagnosticFences;

  vcTileRenderer *pTileRenderer;
  vcAnchor *pCompass;

  float previousFrameDepth;
  udFloat2 currentMouseUV;

  struct
  {
    vcShader *pProgram;
    vcShaderSampler *uniform_texture;
    vcShaderConstantBuffer *uniform_MatrixBlock;

    vcTexture *pSkyboxTexture;
  } skyboxShaderPanorama;

  struct
  {
    vcShader *pProgram;
    vcShaderConstantBuffer *uniform_params;
    vcShaderSampler *uniform_texture;

    vcTexture *pLogoTexture;

    struct
    {
      udFloat4 colour;
      udFloat4 imageSize;
    } params;
  } skyboxShaderTintImage;

  struct
  {
    udUInt2 location;

    vcFramebuffer *pFramebuffer;
    vcTexture *pTexture;
    vcTexture *pDepth;
  } picking;

  struct
  {
    vcShader *pProgram;
    vcShaderSampler *uniform_texture;
    vcShaderConstantBuffer *uniform_params;

    struct
    {
      udFloat4 stepSize;
    } params;
  } blurShader;

  struct
  {
    vcShader *pProgram;
    vcShaderSampler *uniform_texture;
    vcShaderConstantBuffer *uniform_params;

    struct
    {
      udFloat4 stepSizeThickness;  // (stepSize.xy, outline thickness, inner overlay strength)
      udFloat4 colour;    // rgb, (unused)
    } params;
  } selectionShader;
};

udResult vcRender_RecreateUDView(vcState *pProgramState, vcRenderContext *pRenderContext);
udResult vcRender_RenderUD(vcState *pProgramState, vcRenderContext *pRenderContext, vdkRenderView *pRenderView, vcCamera *pCamera, vcRenderData &renderData, bool doPick);

constexpr double kPi = 3.1415926;
constexpr double kSunAngularRadius = 0.00935 / 2.0;
constexpr double kSunSolidAngle = kPi * kSunAngularRadius * kSunAngularRadius;
constexpr double kLengthUnitInMeters = 1.0;

const char kVertexShader[] = R"(
    #version 330
    uniform mat4 model_from_view;
    uniform mat4 view_from_clip;
    layout(location = 0) in vec4 vertex;
    out vec3 view_ray;
    out vec2 v_uv;
    void main() {
      view_ray =
          (model_from_view * vec4((view_from_clip * vertex).xyz, 0.0)).xyz;
      gl_Position = vertex;
      v_uv = vertex.xy * vec2(0.5, 0.5) + vec2(0.5, 0.5);
    })";

//#include "atmosphere/demo/demo.glsl.inc"

enum Luminance {
  // Render the spectral radiance at kLambdaR, kLambdaG, kLambdaB.
  NONE,
  // Render the sRGB luminance, using an approximate (on the fly) conversion
  // from 3 spectral radiance values only (see section 14.3 in <a href=
  // "https://arxiv.org/pdf/1612.04336.pdf">A Qualitative and Quantitative
  //  Evaluation of 8 Clear Sky Models</a>).
  APPROXIMATE,
  // Render the sRGB luminance, precomputed from 15 spectral radiance values
  // (see section 4.4 in <a href=
  // "http://www.oskee.wz.cz/stranka/uploads/SCCG10ElekKmoch.pdf">Real-time
  //  Spectral Scattering in Large-scale Natural Participating Media</a>).
  PRECOMPUTED
};

#include "gl/opengl/vcOpenGL.h"

GLuint vertex_shader_;
GLuint fragment_shader_;
GLuint program_;

bool use_constant_solar_spectrum_ = false;
bool use_ozone_ = true;
Luminance use_luminance_ = NONE;
bool use_half_precision_ = true;
bool use_combined_textures_ = true;
bool do_white_balance_ = false;

double view_distance_meters_ = 9000.0;
double view_zenith_angle_radians_ = 1.47;
double view_azimuth_angle_radians_ = -0.1;
double sun_zenith_angle_radians_ = 1.3;
double sun_azimuth_angle_radians_ = 2.9;
double exposure_ = 10.0;

udUInt2 sceneRes = udUInt2::zero();

void InitModel()
{

  // Values from "Reference Solar Spectral Irradiance: ASTM G-173", ETR column
  // (see http://rredc.nrel.gov/solar/spectra/am1.5/ASTMG173/ASTMG173.html),
  // summed and averaged in each bin (e.g. the value for 360nm is the average
  // of the ASTM G-173 values for all wavelengths between 360 and 370nm).
  // Values in W.m^-2.
  constexpr int kLambdaMin = 360;
  constexpr int kLambdaMax = 830;
  constexpr double kSolarIrradiance[48] = {
    1.11776, 1.14259, 1.01249, 1.14716, 1.72765, 1.73054, 1.6887, 1.61253,
    1.91198, 2.03474, 2.02042, 2.02212, 1.93377, 1.95809, 1.91686, 1.8298,
    1.8685, 1.8931, 1.85149, 1.8504, 1.8341, 1.8345, 1.8147, 1.78158, 1.7533,
    1.6965, 1.68194, 1.64654, 1.6048, 1.52143, 1.55622, 1.5113, 1.474, 1.4482,
    1.41018, 1.36775, 1.34188, 1.31429, 1.28303, 1.26758, 1.2367, 1.2082,
    1.18737, 1.14683, 1.12362, 1.1058, 1.07124, 1.04992
  };
  // Values from http://www.iup.uni-bremen.de/gruppen/molspec/databases/
  // referencespectra/o3spectra2011/index.html for 233K, summed and averaged in
  // each bin (e.g. the value for 360nm is the average of the original values
  // for all wavelengths between 360 and 370nm). Values in m^2.
  constexpr double kOzoneCrossSection[48] = {
    1.18e-27, 2.182e-28, 2.818e-28, 6.636e-28, 1.527e-27, 2.763e-27, 5.52e-27,
    8.451e-27, 1.582e-26, 2.316e-26, 3.669e-26, 4.924e-26, 7.752e-26, 9.016e-26,
    1.48e-25, 1.602e-25, 2.139e-25, 2.755e-25, 3.091e-25, 3.5e-25, 4.266e-25,
    4.672e-25, 4.398e-25, 4.701e-25, 5.019e-25, 4.305e-25, 3.74e-25, 3.215e-25,
    2.662e-25, 2.238e-25, 1.852e-25, 1.473e-25, 1.209e-25, 9.423e-26, 7.455e-26,
    6.566e-26, 5.105e-26, 4.15e-26, 4.228e-26, 3.237e-26, 2.451e-26, 2.801e-26,
    2.534e-26, 1.624e-26, 1.465e-26, 2.078e-26, 1.383e-26, 7.105e-27
  };
  // From https://en.wikipedia.org/wiki/Dobson_unit, in molecules.m^-2.
  constexpr double kDobsonUnit = 2.687e20;
  // Maximum number density of ozone molecules, in m^-3 (computed so at to get
  // 300 Dobson units of ozone - for this we divide 300 DU by the integral of
  // the ozone density profile defined below, which is equal to 15km).
  constexpr double kMaxOzoneNumberDensity = 300.0 * kDobsonUnit / 15000.0;
  // Wavelength independent solar irradiance "spectrum" (not physically
  // realistic, but was used in the original implementation).
  constexpr double kConstantSolarIrradiance = 1.5;
  constexpr double kBottomRadius = 6378137.0;//6360000.0;
  constexpr double kTopRadius = 6420000.0 + 18137.0;
  constexpr double kRayleigh = 1.24062e-6;
  constexpr double kRayleighScaleHeight = 8000.0;
  constexpr double kMieScaleHeight = 1200.0;
  constexpr double kMieAngstromAlpha = 0.0;
  constexpr double kMieAngstromBeta = 5.328e-3;
  constexpr double kMieSingleScatteringAlbedo = 0.9;
  constexpr double kMiePhaseFunctionG = 0.8;
  constexpr double kGroundAlbedo = 0.1;
  const double max_sun_zenith_angle =
    (use_half_precision_ ? 102.0 : 120.0) / 180.0 * kPi;

  atmosphere::DensityProfileLayer
    rayleigh_layer(0.0, 1.0, -1.0 / kRayleighScaleHeight, 0.0, 0.0);
  atmosphere::DensityProfileLayer mie_layer(0.0, 1.0, -1.0 / kMieScaleHeight, 0.0, 0.0);
  // Density profile increasing linearly from 0 to 1 between 10 and 25km, and
  // decreasing linearly from 1 to 0 between 25 and 40km. This is an approximate
  // profile from http://www.kln.ac.lk/science/Chemistry/Teaching_Resources/
  // Documents/Introduction%20to%20atmospheric%20chemistry.pdf (page 10).
  std::vector<atmosphere::DensityProfileLayer> ozone_density;
  ozone_density.push_back(
    atmosphere::DensityProfileLayer(25000.0, 0.0, 0.0, 1.0 / 15000.0, -2.0 / 3.0));
  ozone_density.push_back(
    atmosphere::DensityProfileLayer(0.0, 0.0, 0.0, -1.0 / 15000.0, 8.0 / 3.0));

  std::vector<double> wavelengths;
  std::vector<double> solar_irradiance;
  std::vector<double> rayleigh_scattering;
  std::vector<double> mie_scattering;
  std::vector<double> mie_extinction;
  std::vector<double> absorption_extinction;
  std::vector<double> ground_albedo;
  for (int l = kLambdaMin; l <= kLambdaMax; l += 10) {
    double lambda = static_cast<double>(l) * 1e-3;  // micro-meters
    double mie =
      kMieAngstromBeta / kMieScaleHeight * pow(lambda, -kMieAngstromAlpha);
    wavelengths.push_back(l);
    if (use_constant_solar_spectrum_) {
      solar_irradiance.push_back(kConstantSolarIrradiance);
    }
    else {
      solar_irradiance.push_back(kSolarIrradiance[(l - kLambdaMin) / 10]);
    }
    rayleigh_scattering.push_back(kRayleigh * pow(lambda, -4));
    mie_scattering.push_back(mie * kMieSingleScatteringAlbedo);
    mie_extinction.push_back(mie);
    absorption_extinction.push_back(use_ozone_ ?
      kMaxOzoneNumberDensity * kOzoneCrossSection[(l - kLambdaMin) / 10] :
      0.0);
    ground_albedo.push_back(kGroundAlbedo);
  }

  pModel = new atmosphere::Model(wavelengths, solar_irradiance, kSunAngularRadius,
    kBottomRadius, kTopRadius, { rayleigh_layer }, rayleigh_scattering,
    { mie_layer }, mie_scattering, mie_extinction, kMiePhaseFunctionG,
    ozone_density, absorption_extinction, ground_albedo, max_sun_zenith_angle,
    kLengthUnitInMeters, use_luminance_ == PRECOMPUTED ? 15 : 3,
    use_combined_textures_, use_half_precision_);
    pModel->Init();


    vertex_shader_ = glCreateShader(GL_VERTEX_SHADER);
    const char *const vertex_shader_source = kVertexShader;
    glShaderSource(vertex_shader_, 1, &vertex_shader_source, NULL);
    glCompileShader(vertex_shader_);

    const std::string fragment_shader_str =
      "#version 330\n" +
      std::string(use_luminance_ != NONE ? "#define USE_LUMINANCE\n" : "") +
      "const float kLengthUnitInMeters = " +
      std::to_string(kLengthUnitInMeters) + ";\n" +
      R"demo(
/**
 * Copyright (c) 2017 Eric Bruneton
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holders nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/*<h2>atmosphere/demo/demo.glsl</h2>
<p>This GLSL fragment shader is used to render our demo scene, which consists of
a sphere S on a purely spherical planet P. It is rendered by "ray tracing", i.e.
the vertex shader outputs the view ray direction, and the fragment shader
computes the intersection of this ray with the spheres S and P to produce the
final pixels. The fragment shader also computes the intersection of the light
rays with the sphere S, to compute shadows, as well as the intersections of the
view ray with the shadow volume of S, in order to compute light shafts.
<p>Our fragment shader has the following inputs and outputs:
*/

uniform vec3 camera;
uniform float exposure;
uniform vec3 white_point;
uniform vec3 earth_center;
uniform vec3 sun_direction;
uniform vec2 sun_size;
uniform sampler2D u_colour;
uniform sampler2D u_depth;
uniform mat4 u_inverseViewProjection;
uniform sampler2D u_stars;

in vec3 view_ray;
in vec2 v_uv;
layout(location = 0) out vec4 color;

/*
<p>It uses the following constants, as well as the following atmosphere
rendering functions, defined externally (by the <code>Model</code>'s
<code>GetShader()</code> shader). The <code>USE_LUMINANCE</code> option is used
to select either the functions returning radiance values, or those returning
luminance values (see <a href="../model.h.html">model.h</a>).
*/

const float PI = 3.14159265;
//const vec3 kSphereCenter = vec3(0.0, 0.0, 0.0) / kLengthUnitInMeters;
//const float kSphereRadius = 1.0 / kLengthUnitInMeters;
//const vec3 kSphereAlbedo = vec3(0.8);
const vec3 kGroundAlbedo = vec3(1.0, 0.0, 1.04);

#ifdef USE_LUMINANCE
#define GetSolarRadiance GetSolarLuminance
#define GetSkyRadiance GetSkyLuminance
#define GetSkyRadianceToPoint GetSkyLuminanceToPoint
#define GetSunAndSkyIrradiance GetSunAndSkyIlluminance
#endif

vec3 GetSolarRadiance();
vec3 GetSkyRadiance(vec3 camera, vec3 view_ray, float shadow_length,
    vec3 sun_direction, out vec3 transmittance);
vec3 GetSkyRadianceToPoint(vec3 camera, vec3 point, float shadow_length,
    vec3 sun_direction, out vec3 transmittance);
vec3 GetSunAndSkyIrradiance(
    vec3 p, vec3 normal, vec3 sun_direction, out vec3 sky_irradiance);

/*<h3>Shadows and light shafts</h3>
<p>The functions to compute shadows and light shafts must be defined before we
can use them in the main shader function, so we define them first. Testing if
a point is in the shadow of the sphere S is equivalent to test if the
corresponding light ray intersects the sphere, which is very simple to do.
However, this is only valid for a punctual light source, which is not the case
of the Sun. In the following function we compute an approximate (and biased)
soft shadow by taking the angular size of the Sun into account:
*/

float GetSunVisibility(vec3 point, vec3 sun_direction, float sceneDepth) {
  return float(sceneDepth == 1.0);
  //vec3 p = point - kSphereCenter;
  //float p_dot_v = dot(p, sun_direction);
  //float p_dot_p = dot(p, p);
  //float ray_sphere_center_squared_distance = p_dot_p - p_dot_v * p_dot_v;
  //float distance_to_intersection = -p_dot_v - sqrt(
  //    kSphereRadius * kSphereRadius - ray_sphere_center_squared_distance);
  //if (distance_to_intersection > 0.0) {
  //  // Compute the distance between the view ray and the sphere, and the
  //  // corresponding (tangent of the) subtended angle. Finally, use this to
  //  // compute an approximate sun visibility.
  //  float ray_sphere_distance =
  //      kSphereRadius - sqrt(ray_sphere_center_squared_distance);
  //  float ray_sphere_angular_distance = -ray_sphere_distance / p_dot_v;
  //  return smoothstep(1.0, 0.0, ray_sphere_angular_distance / sun_size.x);
  //}
  //return 1.0;
}

/*
<p>The sphere also partially occludes the sky light, and we approximate this
effect with an ambient occlusion factor. The ambient occlusion factor due to a
sphere is given in <a href=
"http://webserver.dmt.upm.es/~isidoro/tc3/Radiation%20View%20factors.pdf"
>Radiation View Factors</a> (Isidoro Martinez, 1995). In the simple case where
the sphere is fully visible, it is given by the following function:
*/

float GetSkyVisibility(vec3 point, float sceneDepth) {
  return float(sceneDepth == 1.0);
  //vec3 p = point - kSphereCenter;
  //float p_dot_p = dot(p, p);
  //return
  //    1.0 + p.z / sqrt(p_dot_p) * kSphereRadius * kSphereRadius / p_dot_p;
} 

/*
<p>To compute light shafts we need the intersections of the view ray with the
shadow volume of the sphere S. Since the Sun is not a punctual light source this
shadow volume is not a cylinder but a cone (for the umbra, plus another cone for
the penumbra, but we ignore it here):
<svg width="505px" height="200px">
  <style type="text/css"><![CDATA[
    circle { fill: #000000; stroke: none; }
    path { fill: none; stroke: #000000; }
    text { font-size: 16px; font-style: normal; font-family: Sans; }
    .vector { font-weight: bold; }
  ]]></style>
  <path d="m 10,75 455,120"/>
  <path d="m 10,125 455,-120"/>
  <path d="m 120,50 160,130"/>
  <path d="m 138,70 7,0 0,-7"/>
  <path d="m 410,65 40,0 m -5,-5 5,5 -5,5"/>
  <path d="m 20,100 430,0" style="stroke-dasharray:8,4,2,4;"/>
  <path d="m 255,25 0,155" style="stroke-dasharray:2,2;"/>
  <path d="m 280,160 -25,0" style="stroke-dasharray:2,2;"/>
  <path d="m 255,140 60,0" style="stroke-dasharray:2,2;"/>
  <path d="m 300,105 5,-5 5,5 m -5,-5 0,40 m -5,-5 5,5 5,-5"/>
  <path d="m 265,105 5,-5 5,5 m -5,-5 0,60 m -5,-5 5,5 5,-5"/>
  <path d="m 260,80 -5,5 5,5 m -5,-5 85,0 m -5,5 5,-5 -5,-5"/>
  <path d="m 335,95 5,5 5,-5 m -5,5 0,-60 m -5,5 5,-5 5,5"/>
  <path d="m 50,100 a 50,50 0 0 1 2,-14" style="stroke-dasharray:2,1;"/>
  <circle cx="340" cy="100" r="60" style="fill: none; stroke: #000000;"/>
  <circle cx="340" cy="100" r="2.5"/>
  <circle cx="255" cy="160" r="2.5"/>
  <circle cx="120" cy="50" r="2.5"/>
  <text x="105" y="45" class="vector">p</text>
  <text x="240" y="170" class="vector">q</text>
  <text x="425" y="55" class="vector">s</text>
  <text x="135" y="55" class="vector">v</text>
  <text x="345" y="75">R</text>
  <text x="275" y="135">r</text>
  <text x="310" y="125">ρ</text>
  <text x="215" y="120">d</text>
  <text x="290" y="80">δ</text>
  <text x="30" y="95">α</text>
</svg>
<p>Noting, as in the above figure, $\bp$ the camera position, $\bv$ and $\bs$
the unit view ray and sun direction vectors and $R$ the sphere radius (supposed
to be centered on the origin), the point at distance $d$ from the camera is
$\bq=\bp+d\bv$. This point is at a distance $\delta=-\bq\cdot\bs$ from the
sphere center along the umbra cone axis, and at a distance $r$ from this axis
given by $r^2=\bq\cdot\bq-\delta^2$. Finally, at distance $\delta$ along the
axis the umbra cone has radius $\rho=R-\delta\tan\alpha$, where $\alpha$ is
the Sun's angular radius. The point at distance $d$ from the camera is on the
shadow cone only if $r^2=\rho^2$, i.e. only if
\begin{equation}
(\bp+d\bv)\cdot(\bp+d\bv)-((\bp+d\bv)\cdot\bs)^2=
(R+((\bp+d\bv)\cdot\bs)\tan\alpha)^2
\end{equation}
Developping this gives a quadratic equation for $d$:
\begin{equation}
ad^2+2bd+c=0
\end{equation}
where
<ul>
<li>$a=1-l(\bv\cdot\bs)^2$,</li>
<li>$b=\bp\cdot\bv-l(\bp\cdot\bs)(\bv\cdot\bs)-\tan(\alpha)R(\bv\cdot\bs)$,</li>
<li>$c=\bp\cdot\bp-l(\bp\cdot\bs)^2-2\tan(\alpha)R(\bp\cdot\bs)-R^2$,</li>
<li>$l=1+\tan^2\alpha$</li>
</ul>
From this we deduce the two possible solutions for $d$, which must be clamped to
the actual shadow part of the mathematical cone (i.e. the slab between the
sphere center and the cone apex or, in other words, the points for which
$\delta$ is between $0$ and $R/\tan\alpha$). The following function implements
these equations:
*/

void GetSphereShadowInOut(vec3 view_direction, vec3 sun_direction,
    out float d_in, out float d_out) {
  //vec3 pos = camera - kSphereCenter;
  //float pos_dot_sun = dot(pos, sun_direction);
  //float view_dot_sun = dot(view_direction, sun_direction);
  //float k = sun_size.x;
  //float l = 1.0 + k * k;
  //float a = 1.0 - l * view_dot_sun * view_dot_sun;
  //float b = dot(pos, view_direction) - l * pos_dot_sun * view_dot_sun -
  //    k * kSphereRadius * view_dot_sun;
  //float c = dot(pos, pos) - l * pos_dot_sun * pos_dot_sun -
  //    2.0 * k * kSphereRadius * pos_dot_sun - kSphereRadius * kSphereRadius;
  //float discriminant = b * b - a * c;
  //if (discriminant > 0.0) {
  //  d_in = max(0.0, (-b - sqrt(discriminant)) / a);
  //  d_out = (-b + sqrt(discriminant)) / a;
  //  // The values of d for which delta is equal to 0 and kSphereRadius / k.
  //  float d_base = -pos_dot_sun / view_dot_sun;
  //  float d_apex = -(pos_dot_sun + kSphereRadius / k) / view_dot_sun;
  //  if (view_dot_sun > 0.0) {
  //    d_in = max(d_in, d_apex);
  //    d_out = a > 0.0 ? min(d_out, d_base) : d_base;
  //  } else {
  //    d_in = a > 0.0 ? max(d_in, d_base) : d_base;
  //    d_out = min(d_out, d_apex);
  //  }
  //} else {
  //  d_in = 0.0;
  //  d_out = 0.0;
  //}
    d_in = 0.0;
    d_out = 0.0;
}

/*<h3>Main shading function</h3>
<p>Using these functions we can now implement the main shader function, which
computes the radiance from the scene for a given view ray. This function first
tests if the view ray intersects the sphere S. If so it computes the sun and
sky light received by the sphere at the intersection point, combines this with
the sphere BRDF and the aerial perspective between the camera and the sphere.
It then does the same with the ground, i.e. with the planet sphere P, and then
computes the sky radiance and transmittance. Finally, all these terms are
composited together (an opacity is also computed for each object, using an
approximate view cone - sphere intersection factor) to get the final radiance.
<p>We start with the computation of the intersections of the view ray with the
shadow volume of the sphere, because they are needed to get the aerial
perspective for the sphere and the planet:
*/

void main() {
  // Normalized view direction vector.
  vec3 view_direction = normalize(view_ray);
  // Tangent of the angle subtended by this fragment.
  float fragment_angular_size =
      length(dFdx(view_ray) + dFdy(view_ray)) / length(view_ray);

  vec4 starsTexture = texture(u_stars, view_direction.xy * 0.5 + 0.5);
  float sceneDepth = texture(u_depth, v_uv).x;
  vec4 sceneColour = texture(u_colour, v_uv);
  sceneColour.xyz = pow(sceneColour.xyz, vec3(2.2));

  gl_FragDepth = sceneDepth;

  float shadow_in;
  float shadow_out;
  GetSphereShadowInOut(view_direction, sun_direction, shadow_in, shadow_out);

  // Hack to fade out light shafts when the Sun is very close to the horizon.
  float lightshaft_fadein_hack = smoothstep(
      0.02, 0.04, dot(normalize(camera - earth_center), sun_direction));

)demo" +
R"demo(
/*
<p>We then test whether the view ray intersects the sphere S or not. If it does,
we compute an approximate (and biased) opacity value, using the same
approximation as in <code>GetSunVisibility</code>:
*/

  // Compute the distance between the view ray line and the sphere center,
  // and the distance between the camera and the intersection of the view
  // ray with the sphere (or NaN if there is no intersection).
  vec3 p = vec3(0.0);//camera - kSphereCenter;
  float p_dot_v = dot(p, view_direction);
  float p_dot_p = dot(p, p);
  //float ray_sphere_center_squared_distance = p_dot_p - p_dot_v * p_dot_v;
  //float distance_to_intersection = -p_dot_v - sqrt(
  //    kSphereRadius * kSphereRadius - ray_sphere_center_squared_distance);

  float distance_to_intersection = sqrt(-1.0);
  vec4 pp = u_inverseViewProjection * vec4(v_uv * 2.0 - vec2(1.0), sceneDepth * 2.0 - 1.0, 1.0);
  pp /= pp.w;
  if (sceneDepth < 1.0)
    distance_to_intersection = length(camera - pp.xyz);

  // Compute the radiance reflected by the sphere, if the ray intersects it.
  float sphere_alpha = 0.0;
  vec3 geometry_radiance = vec3(0.0);
  if (distance_to_intersection > 0.0) {
    // Compute the distance between the view ray and the sphere, and the
    // corresponding (tangent of the) subtended angle. Finally, use this to
    // compute the approximate analytic antialiasing factor sphere_alpha.
    //float ray_sphere_distance =
    //    kSphereRadius - sqrt(ray_sphere_center_squared_distance);
    //float ray_sphere_angular_distance = -ray_sphere_distance / p_dot_v;
    sphere_alpha = 1.0;//min(1.0, 1.0 - distance_to_intersection * 0.001);
      //  min(ray_sphere_angular_distance / fragment_angular_size, 1.0);

/*
<p>We can then compute the intersection point and its normal, and use them to
get the sun and sky irradiance received at this point. The reflected radiance
follows, by multiplying the irradiance with the sphere BRDF:
*/
    vec3 point = camera + view_direction * distance_to_intersection;
    vec3 normal = vec3(0,0,1);//normalize(point - kSphereCenter);

    // Compute the radiance reflected by the sphere.
    vec3 sky_irradiance;
    vec3 sun_irradiance = GetSunAndSkyIrradiance(
        point - earth_center, normal, sun_direction, sky_irradiance);
    geometry_radiance =
        sceneColour.xyz * (1.0 / PI) * (sun_irradiance + sky_irradiance);

/*
<p>Finally, we take into account the aerial perspective between the camera and
the sphere, which depends on the length of this segment which is in shadow:
*/
    float shadow_length = 0;//
        //max(0.0, min(shadow_out, distance_to_intersection) - shadow_in) *
        //lightshaft_fadein_hack;
    vec3 transmittance;
    vec3 in_scatter = GetSkyRadianceToPoint(camera - earth_center,
        point - earth_center, shadow_length, sun_direction, transmittance);
    geometry_radiance = geometry_radiance * transmittance + in_scatter;
  }

/*
<p>In the following we repeat the same steps as above, but for the planet sphere
P instead of the sphere S (a smooth opacity is not really needed here, so we
don't compute it. Note also how we modulate the sun and sky irradiance received
on the ground by the sun and sky visibility factors):
*/

  // Compute the distance between the view ray line and the Earth center,
  // and the distance between the camera and the intersection of the view
  // ray with the ground (or NaN if there is no intersection).
  p = camera - earth_center;
  p_dot_v = dot(p, view_direction);
  p_dot_p = dot(p, p);
  float ray_earth_center_squared_distance = p_dot_p - p_dot_v * p_dot_v;
  distance_to_intersection = -p_dot_v - sqrt(
      earth_center.z * earth_center.z - ray_earth_center_squared_distance);

  // Compute the radiance reflected by the ground, if the ray intersects it.
  float ground_alpha = 0.0;
  vec3 ground_radiance = vec3(0.0);
  if (distance_to_intersection > 0.0) {
    vec3 point = camera + view_direction * distance_to_intersection;
    vec3 normal = normalize(point - earth_center);

    // Compute the radiance reflected by the ground.
    vec3 sky_irradiance;
    vec3 sun_irradiance = GetSunAndSkyIrradiance(
        point - earth_center, normal, sun_direction, sky_irradiance);
    ground_radiance = kGroundAlbedo * (1.0 / PI) * (
        sun_irradiance * GetSunVisibility(point, sun_direction, sceneDepth) +
        sky_irradiance * GetSkyVisibility(point, sceneDepth));

    float shadow_length =
        max(0.0, min(shadow_out, distance_to_intersection) - shadow_in) *
        lightshaft_fadein_hack;
    vec3 transmittance;
    vec3 in_scatter = GetSkyRadianceToPoint(camera - earth_center,
        point - earth_center, shadow_length, sun_direction, transmittance);
    ground_radiance = ground_radiance * transmittance + in_scatter;
    ground_alpha = 1.0;
  }

/*
<p>Finally, we compute the radiance and transmittance of the sky, and composite
together, from back to front, the radiance and opacities of all the objects of
the scene:
*/

  // Compute the radiance of the sky.
  float shadow_length = max(0.0, shadow_out - shadow_in) *
      lightshaft_fadein_hack;
  vec3 transmittance;
  vec3 radiance = GetSkyRadiance(
      camera - earth_center, view_direction, shadow_length, sun_direction,
      transmittance);

  // If the view ray intersects the Sun, add the Sun radiance.
  if (dot(view_direction, sun_direction) > sun_size.y) {
    radiance = radiance + transmittance * GetSolarRadiance();
  }

  radiance = mix(radiance, ground_radiance, ground_alpha);
  radiance = mix(radiance, geometry_radiance, sphere_alpha);

  color.rgb = 
      pow(vec3(1.0) - exp(-radiance / white_point * exposure), vec3(1.0 / 2.2));
  color.a = 1.0;

  float stars_fadein_hack = clamp(smoothstep(
      -0.2, 0.05, dot(normalize(camera - earth_center), sun_direction)), 0.0, 1.0);

  float height_stars_fadein_hack = clamp(smoothstep(10000, 60000, camera.z), 0.0, 1.0);
  vec3 t = pow(transmittance, vec3(1.0 / 2.2));
  color.rgb = mix(color.rgb, starsTexture.xyz, ((1.0 - stars_fadein_hack) + height_stars_fadein_hack)*min(1.0, dot(t,t)));
}
)demo";

    const char *fragment_shader_source = fragment_shader_str.c_str();
    fragment_shader_ = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader_, 1, &fragment_shader_source, NULL);
    glCompileShader(fragment_shader_);

    if (program_ != 0) {
      glDeleteProgram(program_);
    }
    program_ = glCreateProgram();
    glAttachShader(program_, vertex_shader_);
    glAttachShader(program_, fragment_shader_);
    glAttachShader(program_, pModel->shader());
    glLinkProgram(program_);
    glDetachShader(program_, vertex_shader_);
    glDetachShader(program_, fragment_shader_);
    glDetachShader(program_, pModel->shader());

    /*
    <p>Finally, it sets the uniforms of this program that can be set once and for
    all (in our case this includes the <code>Model</code>'s texture uniforms,
    because our demo app does not have any texture of its own):
    */

    glUseProgram(program_);
    pModel->SetProgramUniforms(program_, 0, 1, 2, 3);
    double white_point_r = 1.0;
    double white_point_g = 1.0;
    double white_point_b = 1.0;
    if (do_white_balance_) {
      atmosphere::Model::ConvertSpectrumToLinearSrgb(wavelengths, solar_irradiance,
        &white_point_r, &white_point_g, &white_point_b);
      double white_point = (white_point_r + white_point_g + white_point_b) / 3.0;
      white_point_r /= white_point;
      white_point_g /= white_point;
      white_point_b /= white_point;
    }
    glUniform3f(glGetUniformLocation(program_, "white_point"),
      white_point_r, white_point_g, white_point_b);

    glUniform2f(glGetUniformLocation(program_, "sun_size"),
      tan(kSunAngularRadius),
      cos(kSunAngularRadius));
}

void DrawAtmosphere(vcState *pProgramState, vcRenderContext *pRenderContext, vcTexture *pSceneColour, vcTexture *pSceneDepth)
{
  glUseProgram(program_);
  pModel->SetProgramUniforms(program_, 0, 1, 2, 3);

  if (!pProgramState->gis.isProjected || pProgramState->gis.zone.projection >= udGZPT_TransverseMercator) //TODO: Fix this list
  {
    udDouble3 earthCenterMaybe = pProgramState->camera.position;

    if (pProgramState->gis.isProjected)
      earthCenterMaybe.z = -pProgramState->gis.zone.semiMajorAxis;
    else
      earthCenterMaybe.z = -6378137.000;

    //4978

    glUniform3f(glGetUniformLocation(program_, "earth_center"),
      earthCenterMaybe.x, earthCenterMaybe.y, earthCenterMaybe.z);//-kBottomRadius / kLengthUnitInMeters);
  }
  else
  {
    udGeoZone destZone = {};
    udGeoZone_SetFromSRID(&destZone, 4978);
    udDouble3 earthCenterMaybe = udGeoZone_TransformPoint(pProgramState->camera.position, pProgramState->gis.zone, destZone);

    //4978

    glUniform3f(glGetUniformLocation(program_, "earth_center"),
      earthCenterMaybe.x, earthCenterMaybe.y, earthCenterMaybe.z);//-kBottomRadius / kLengthUnitInMeters);
  }

  const float kFovY = 50.0 / 180.0 * kPi;
  const float kTanFovY = tan(kFovY / 2.0);
  float aspect_ratio = static_cast<float>(sceneRes.x) / sceneRes.y;

  //// Transform matrix from clip space to camera space (i.e. the inverse of a
  //// GL_PROJECTION matrix).
  //float view_from_clip[16] = {
  //  kTanFovY * aspect_ratio, 0.0, 0.0, 0.0,
  //  0.0, kTanFovY, 0.0, 0.0,
  //  0.0, 0.0, 0.0, -1.0,
  //  0.0, 0.0, 1.0, 1.0
  //};
  //glUniformMatrix4fv(glGetUniformLocation(program_, "view_from_clip"), 1, true,
  //  view_from_clip);

  udFloat4x4 inverseProjection = udFloat4x4::create(udInverse(pProgramState->camera.matrices.projection));
  udFloat4x4 inverseView = udFloat4x4::create(udInverse(pProgramState->camera.matrices.view));
  udFloat4x4 inverseViewProjection = udFloat4x4::create(pProgramState->camera.matrices.inverseViewProjection);
  glUniformMatrix4fv(glGetUniformLocation(program_, "view_from_clip"), 1, false, inverseProjection.a);
  glUniformMatrix4fv(glGetUniformLocation(program_, "u_inverseViewProjection"), 1, false, inverseViewProjection.a);


  // Unit vectors of the camera frame, expressed in world space.
  float cos_z = cos(view_zenith_angle_radians_);
  float sin_z = sin(view_zenith_angle_radians_);
  float cos_a = cos(view_azimuth_angle_radians_);
  float sin_a = sin(view_azimuth_angle_radians_);
  float ux[3] = { -sin_a, cos_a, 0.0 };
  float uy[3] = { -cos_z * cos_a, -cos_z * sin_a, sin_z };
  float uz[3] = { sin_z * cos_a, sin_z * sin_a, cos_z };
  float l = view_distance_meters_ / kLengthUnitInMeters;

  // Transform matrix from camera frame to world space (i.e. the inverse of a
  // GL_MODELVIEW matrix).
  float model_from_view[16] = {
    ux[0], uy[0], uz[0], uz[0] * l,
    ux[1], uy[1], uz[1], uz[1] * l,
    ux[2], uy[2], uz[2], uz[2] * l,
    0.0, 0.0, 0.0, 1.0
  };
  VERIFY_GL();
  GLuint camLoc = glGetUniformLocation(program_, "camera");
  GLuint exposureLoc = glGetUniformLocation(program_, "exposure");
  VERIFY_GL();
  glUniform3f(camLoc,
    pProgramState->camera.position.x,
    pProgramState->camera.position.y,
    pProgramState->camera.position.z);
  VERIFY_GL();
  glUniform1f(exposureLoc,
    use_luminance_ != NONE ? exposure_ * 1e-5 : exposure_);
  VERIFY_GL();
  glUniformMatrix4fv(glGetUniformLocation(program_, "model_from_view"),
    1, false, inverseView.a);
  VERIFY_GL();
  glUniform3f(glGetUniformLocation(program_, "sun_direction"),
    cos(sun_azimuth_angle_radians_) * sin(sun_zenith_angle_radians_),
    sin(sun_azimuth_angle_radians_) * sin(sun_zenith_angle_radians_),
    cos(sun_zenith_angle_radians_));
  VERIFY_GL();

  VERIFY_GL();
  glActiveTexture(GL_TEXTURE4);
  glBindTexture(GL_TEXTURE_2D, pSceneColour->id);
  glUniform1i(glGetUniformLocation(program_, "u_colour"), 4);

  VERIFY_GL();
  glActiveTexture(GL_TEXTURE5);
  glBindTexture(GL_TEXTURE_2D, pSceneDepth->id);
  glUniform1i(glGetUniformLocation(program_, "u_depth"), 5);

  glActiveTexture(GL_TEXTURE6);
  glBindTexture(GL_TEXTURE_2D, pRenderContext->skyboxShaderPanorama.pSkyboxTexture->id);
  glUniform1i(glGetUniformLocation(program_, "u_stars"), 6);

  //glBindVertexArray(full_screen_quad_vao_);
  //glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
  vcMesh_Render(gInternalMeshes[vcInternalMeshType_ScreenQuad]);
  glBindVertexArray(0);
}

int previous_mouse_x_ = 0;
int previous_mouse_y_ = 0;

void HandleFakeInput()
{
  int mouse_x = ImGui::GetIO().MousePos.x;
  int mouse_y = ImGui::GetIO().MousePos.y;

  if (ImGui::GetIO().MouseDown[0])
  {
    constexpr double kScale = 500.0;
    if (ImGui::GetIO().KeyCtrl) {
      sun_zenith_angle_radians_ -= (previous_mouse_y_ - mouse_y) / kScale;
      sun_zenith_angle_radians_ =
        udMax(0.0, udMin(kPi, sun_zenith_angle_radians_));
      sun_azimuth_angle_radians_ += (previous_mouse_x_ - mouse_x) / kScale;
    }
    else {
      view_zenith_angle_radians_ += (previous_mouse_y_ - mouse_y) / kScale;
      view_zenith_angle_radians_ =
        udMax(0.0, udMin(kPi / 2.0, view_zenith_angle_radians_));
      view_azimuth_angle_radians_ += (previous_mouse_x_ - mouse_x) / kScale;
    }
  }
  previous_mouse_x_ = mouse_x;
  previous_mouse_y_ = mouse_y;
}

udResult vcRender_Init(vcState *pProgramState, vcRenderContext **ppRenderContext, udWorkerPool *pWorkerPool, const udUInt2 &sceneResolution)
{
  sceneRes = sceneResolution;
  udResult result;
  vcRenderContext *pRenderContext = nullptr;

  InitModel();

  UD_ERROR_NULL(ppRenderContext, udR_InvalidParameter_);

  UD_ERROR_CHECK(vcInternalModels_Init());

  pRenderContext = udAllocType(vcRenderContext, 1, udAF_Zero);
  UD_ERROR_NULL(pRenderContext, udR_MemoryAllocationFailure);

  pRenderContext->viewShedRenderingContext.pDepthBuffer = udAllocType(float, ViewShedMapRes.x * ViewShedMapRes.y, udAF_Zero);
  UD_ERROR_CHECK(vcTexture_Create(&pRenderContext->viewShedRenderingContext.pUDDepthTexture, ViewShedMapRes.x, ViewShedMapRes.y, nullptr, vcTextureFormat_D32F, vcTFM_Nearest, false, vcTWM_Repeat, vcTCF_Dynamic));
  UD_ERROR_CHECK(vcTexture_Create(&pRenderContext->viewShedRenderingContext.pDepthTex, ViewShedMapRes.x, ViewShedMapRes.y, nullptr, vcTextureFormat_D24S8, vcTFM_Nearest, false, vcTWM_Repeat, vcTCF_RenderTarget));
  UD_ERROR_CHECK(vcTexture_Create(&pRenderContext->viewShedRenderingContext.pDummyColour, ViewShedMapRes.x, ViewShedMapRes.y, nullptr, vcTextureFormat_BGRA8, vcTFM_Nearest, false, vcTWM_Repeat, vcTCF_RenderTarget));
  UD_ERROR_IF(!vcFramebuffer_Create(&pRenderContext->viewShedRenderingContext.pFramebuffer, pRenderContext->viewShedRenderingContext.pDummyColour, pRenderContext->viewShedRenderingContext.pDepthTex), udR_InternalError);

  UD_ERROR_IF(!vcShader_CreateFromText(&pRenderContext->udRenderContext.presentShader.pProgram, g_udVertexShader, g_udFragmentShader, vcP3UV2VertexLayout), udR_InternalError);
  UD_ERROR_IF(!vcShader_CreateFromText(&pRenderContext->visualizationShader.pProgram, g_VisualizationVertexShader, g_VisualizationFragmentShader, vcP3UV2VertexLayout), udR_InternalError);
  UD_ERROR_IF(!vcShader_CreateFromText(&pRenderContext->shadowShader.pProgram, g_ViewShedVertexShader, g_ViewShedFragmentShader, vcP3UV2VertexLayout), udR_InternalError);
  UD_ERROR_IF(!vcShader_CreateFromText(&pRenderContext->skyboxShaderPanorama.pProgram, g_vcSkyboxVertexShaderPanorama, g_vcSkyboxFragmentShaderPanorama, vcP3UV2VertexLayout), udR_InternalError);
  UD_ERROR_IF(!vcShader_CreateFromText(&pRenderContext->skyboxShaderTintImage.pProgram, g_vcSkyboxVertexShaderImageColour, g_vcSkyboxFragmentShaderImageColour, vcP3UV2VertexLayout), udR_InternalError);
  UD_ERROR_IF(!vcShader_CreateFromText(&pRenderContext->udRenderContext.splatIdShader.pProgram, g_udVertexShader, g_udSplatIdFragmentShader, vcP3UV2VertexLayout), udR_InternalError);
  UD_ERROR_IF(!vcShader_CreateFromText(&pRenderContext->fxaaShader.pProgram, g_FXAAVertexShader, g_FXAAFragmentShader, vcP3UV2VertexLayout), udR_InternalError);

  UD_ERROR_CHECK(vcTexture_AsyncCreateFromFilename(&pRenderContext->skyboxShaderPanorama.pSkyboxTexture, pWorkerPool, "asset://assets/skyboxes/stars.jpg", vcTFM_Linear));
  UD_ERROR_CHECK(vcCompass_Create(&pRenderContext->pCompass));

  UD_ERROR_IF(!vcShader_Bind(pRenderContext->visualizationShader.pProgram), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetSamplerIndex(&pRenderContext->visualizationShader.uniform_texture, pRenderContext->visualizationShader.pProgram, "u_texture"), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetSamplerIndex(&pRenderContext->visualizationShader.uniform_depth, pRenderContext->visualizationShader.pProgram, "u_depth"), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetConstantBuffer(&pRenderContext->visualizationShader.uniform_vertParams, pRenderContext->visualizationShader.pProgram, "u_vertParams", sizeof(pRenderContext->visualizationShader.vertParams)), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetConstantBuffer(&pRenderContext->visualizationShader.uniform_fragParams, pRenderContext->visualizationShader.pProgram, "u_fragParams", sizeof(pRenderContext->visualizationShader.fragParams)), udR_InternalError);

  UD_ERROR_IF(!vcShader_Bind(pRenderContext->shadowShader.pProgram), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetSamplerIndex(&pRenderContext->shadowShader.uniform_depth, pRenderContext->shadowShader.pProgram, "u_depth"), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetSamplerIndex(&pRenderContext->shadowShader.uniform_shadowMapAtlas, pRenderContext->shadowShader.pProgram, "u_shadowMapAtlas"), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetConstantBuffer(&pRenderContext->shadowShader.uniform_params, pRenderContext->shadowShader.pProgram, "u_params", sizeof(pRenderContext->shadowShader.params)), udR_InternalError);

  UD_ERROR_IF(!vcShader_Bind(pRenderContext->skyboxShaderPanorama.pProgram), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetSamplerIndex(&pRenderContext->skyboxShaderPanorama.uniform_texture, pRenderContext->skyboxShaderPanorama.pProgram, "u_texture"), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetConstantBuffer(&pRenderContext->skyboxShaderPanorama.uniform_MatrixBlock, pRenderContext->skyboxShaderPanorama.pProgram, "u_EveryFrame", sizeof(udFloat4x4)), udR_InternalError);

  UD_ERROR_IF(!vcShader_Bind(pRenderContext->skyboxShaderTintImage.pProgram), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetSamplerIndex(&pRenderContext->skyboxShaderTintImage.uniform_texture, pRenderContext->skyboxShaderTintImage.pProgram, "u_texture"), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetConstantBuffer(&pRenderContext->skyboxShaderTintImage.uniform_params, pRenderContext->skyboxShaderTintImage.pProgram, "u_EveryFrame", sizeof(pRenderContext->skyboxShaderTintImage.params)), udR_InternalError);

  UD_ERROR_IF(!vcShader_Bind(pRenderContext->udRenderContext.presentShader.pProgram), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetSamplerIndex(&pRenderContext->udRenderContext.presentShader.uniform_texture, pRenderContext->udRenderContext.presentShader.pProgram, "u_texture"), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetSamplerIndex(&pRenderContext->udRenderContext.presentShader.uniform_depth, pRenderContext->udRenderContext.presentShader.pProgram, "u_depth"), udR_InternalError);

  UD_ERROR_IF(!vcShader_Bind(pRenderContext->udRenderContext.splatIdShader.pProgram), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetConstantBuffer(&pRenderContext->udRenderContext.splatIdShader.uniform_params, pRenderContext->udRenderContext.splatIdShader.pProgram, "u_params", sizeof(pRenderContext->udRenderContext.splatIdShader.params)), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetSamplerIndex(&pRenderContext->udRenderContext.splatIdShader.uniform_depth, pRenderContext->udRenderContext.splatIdShader.pProgram, "u_depth"), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetSamplerIndex(&pRenderContext->udRenderContext.splatIdShader.uniform_texture, pRenderContext->udRenderContext.splatIdShader.pProgram, "u_texture"), udR_InternalError);

  UD_ERROR_IF(!vcShader_Bind(pRenderContext->fxaaShader.pProgram), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetSamplerIndex(&pRenderContext->fxaaShader.uniform_texture, pRenderContext->fxaaShader.pProgram, "u_texture"), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetSamplerIndex(&pRenderContext->fxaaShader.uniform_depth, pRenderContext->fxaaShader.pProgram, "u_depth"), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetConstantBuffer(&pRenderContext->fxaaShader.uniform_params, pRenderContext->fxaaShader.pProgram, "u_params", sizeof(pRenderContext->fxaaShader.params)), udR_InternalError);

  UD_ERROR_IF(!vcShader_CreateFromText(&pRenderContext->blurShader.pProgram, g_BlurVertexShader, g_BlurFragmentShader, vcP3UV2VertexLayout), udR_InternalError);
  UD_ERROR_IF(!vcShader_Bind(pRenderContext->blurShader.pProgram), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetSamplerIndex(&pRenderContext->blurShader.uniform_texture, pRenderContext->blurShader.pProgram, "u_texture"), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetConstantBuffer(&pRenderContext->blurShader.uniform_params, pRenderContext->blurShader.pProgram, "u_EveryFrame", sizeof(pRenderContext->blurShader.params)), udR_InternalError);

  UD_ERROR_IF(!vcShader_CreateFromText(&pRenderContext->selectionShader.pProgram, g_HighlightVertexShader, g_HighlightFragmentShader, vcP3UV2VertexLayout), udR_InternalError);
  UD_ERROR_IF(!vcShader_Bind(pRenderContext->selectionShader.pProgram), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetSamplerIndex(&pRenderContext->selectionShader.uniform_texture, pRenderContext->selectionShader.pProgram, "u_texture"), udR_InternalError);
  UD_ERROR_IF(!vcShader_GetConstantBuffer(&pRenderContext->selectionShader.uniform_params, pRenderContext->selectionShader.pProgram, "u_EveryFrame", sizeof(pRenderContext->selectionShader.params)), udR_InternalError);

  UD_ERROR_CHECK(vcPolygonModel_CreateShaders());
  UD_ERROR_CHECK(vcImageRenderer_Init());

  UD_ERROR_IF(!vcShader_Bind(nullptr), udR_InternalError);

  UD_ERROR_CHECK(vcTileRenderer_Create(&pRenderContext->pTileRenderer, &pProgramState->settings));
  UD_ERROR_CHECK(vcFenceRenderer_Create(&pRenderContext->pDiagnosticFences));

  UD_ERROR_CHECK(vcRender_ResizeScene(pProgramState, pRenderContext, sceneResolution.x, sceneResolution.y));

  *ppRenderContext = pRenderContext;
  pRenderContext = nullptr;
  result = udR_Success;
epilogue:

  if (pRenderContext != nullptr)
    vcRender_Destroy(pProgramState, &pRenderContext);

  return result;
}

udResult vcRender_Destroy(vcState *pProgramState, vcRenderContext **ppRenderContext)
{
  if (ppRenderContext == nullptr || *ppRenderContext == nullptr)
    return udR_Success;

  udResult result;
  vcRenderContext *pRenderContext = nullptr;

  UD_ERROR_NULL(ppRenderContext, udR_InvalidParameter_);

  pRenderContext = *ppRenderContext;
  *ppRenderContext = nullptr;

  if (pProgramState->pVDKContext != nullptr)
  {
    if (pRenderContext->viewShedRenderingContext.pRenderView != nullptr && vdkRenderView_Destroy(&pRenderContext->viewShedRenderingContext.pRenderView) != vE_Success)
      UD_ERROR_SET(udR_InternalError);

    if (pRenderContext->udRenderContext.pRenderView != nullptr && vdkRenderView_Destroy(&pRenderContext->udRenderContext.pRenderView) != vE_Success)
      UD_ERROR_SET(udR_InternalError);

    if (vdkRenderContext_Destroy(&pRenderContext->udRenderContext.pRenderer) != vE_Success)
      UD_ERROR_SET(udR_InternalError);
  }

  vcShader_DestroyShader(&pRenderContext->udRenderContext.presentShader.pProgram);
  vcShader_DestroyShader(&pRenderContext->visualizationShader.pProgram);
  vcShader_DestroyShader(&pRenderContext->fxaaShader.pProgram);
  vcShader_DestroyShader(&pRenderContext->shadowShader.pProgram);
  vcShader_DestroyShader(&pRenderContext->skyboxShaderPanorama.pProgram);
  vcShader_DestroyShader(&pRenderContext->skyboxShaderTintImage.pProgram);
  vcShader_DestroyShader(&pRenderContext->udRenderContext.splatIdShader.pProgram);
  vcShader_DestroyShader(&pRenderContext->blurShader.pProgram);
  vcShader_DestroyShader(&pRenderContext->selectionShader.pProgram);

  vcTexture_Destroy(&pRenderContext->skyboxShaderPanorama.pSkyboxTexture);
  UD_ERROR_CHECK(vcCompass_Destroy(&pRenderContext->pCompass));

  UD_ERROR_CHECK(vcPolygonModel_DestroyShaders());
  UD_ERROR_CHECK(vcImageRenderer_Destroy());

  udFree(pRenderContext->udRenderContext.pColorBuffer);
  udFree(pRenderContext->udRenderContext.pDepthBuffer);
  udFree(pRenderContext->viewShedRenderingContext.pDepthBuffer);

  UD_ERROR_CHECK(vcTileRenderer_Destroy(&pRenderContext->pTileRenderer));
  UD_ERROR_CHECK(vcFenceRenderer_Destroy(&pRenderContext->pDiagnosticFences));

  UD_ERROR_CHECK(vcInternalModels_Deinit());
  result = udR_Success;

epilogue:
  vcTexture_Destroy(&pRenderContext->viewShedRenderingContext.pUDDepthTexture);
  vcTexture_Destroy(&pRenderContext->viewShedRenderingContext.pDepthTex);
  vcTexture_Destroy(&pRenderContext->viewShedRenderingContext.pDummyColour);
  vcFramebuffer_Destroy(&pRenderContext->viewShedRenderingContext.pFramebuffer);

  vcTexture_Destroy(&pRenderContext->udRenderContext.pColourTex);
  vcTexture_Destroy(&pRenderContext->udRenderContext.pDepthTex);
  vcFramebuffer_Destroy(&pRenderContext->udRenderContext.pFramebuffer);
  for (int i = 0; i < vcRender_RenderBufferCount; ++i)
  {
    vcTexture_Destroy(&pRenderContext->pTexture[i]);
    vcTexture_Destroy(&pRenderContext->pDepthTexture[i]);
    vcFramebuffer_Destroy(&pRenderContext->pFramebuffer[i]);
  }

  vcTexture_Destroy(&pRenderContext->picking.pTexture);
  vcTexture_Destroy(&pRenderContext->picking.pDepth);
  vcFramebuffer_Destroy(&pRenderContext->picking.pFramebuffer);

  for (int i = 0; i < 2; ++i)
  {
    vcTexture_Destroy(&pRenderContext->pAuxiliaryTextures[i]);
    vcFramebuffer_Destroy(&pRenderContext->pAuxiliaryFramebuffers[i]);
  }

  udFree(pRenderContext);
  return result;
}

udResult vcRender_SetVaultContext(vcState *pProgramState, vcRenderContext *pRenderContext)
{
  udResult result = udR_Success;

  UD_ERROR_NULL(pRenderContext, udR_InvalidParameter_);

  if (vdkRenderContext_Create(pProgramState->pVDKContext, &pRenderContext->udRenderContext.pRenderer) != vE_Success)
    UD_ERROR_SET(udR_InternalError);

epilogue:
  return result;
}

udResult vcRender_ResizeScene(vcState *pProgramState, vcRenderContext *pRenderContext, const uint32_t width, const uint32_t height)
{
  sceneRes.x = width;
  sceneRes.y = height;
  udResult result = udR_Success;

  uint32_t widthIncr = width + (width % vcRender_SceneSizeIncrement != 0 ? vcRender_SceneSizeIncrement - width % vcRender_SceneSizeIncrement : 0);
  uint32_t heightIncr = height + (height % vcRender_SceneSizeIncrement != 0 ? vcRender_SceneSizeIncrement - height % vcRender_SceneSizeIncrement : 0);

  UD_ERROR_NULL(pRenderContext, udR_InvalidParameter_);
  UD_ERROR_IF(width == 0, udR_InvalidParameter_);
  UD_ERROR_IF(height == 0, udR_InvalidParameter_);

  pRenderContext->sceneResolution.x = widthIncr;
  pRenderContext->sceneResolution.y = heightIncr;
  pRenderContext->originalSceneResolution.x = width;
  pRenderContext->originalSceneResolution.y = height;

  //Resize CPU Targets
  udFree(pRenderContext->udRenderContext.pColorBuffer);
  udFree(pRenderContext->udRenderContext.pDepthBuffer);

  pRenderContext->udRenderContext.pColorBuffer = udAllocType(uint32_t, pRenderContext->sceneResolution.x * pRenderContext->sceneResolution.y, udAF_Zero);
  UD_ERROR_NULL(pRenderContext->udRenderContext.pColorBuffer, udR_MemoryAllocationFailure);

  pRenderContext->udRenderContext.pDepthBuffer = udAllocType(float, pRenderContext->sceneResolution.x * pRenderContext->sceneResolution.y, udAF_Zero);
  UD_ERROR_NULL(pRenderContext->udRenderContext.pDepthBuffer, udR_MemoryAllocationFailure);

  //Resize GPU Targets
  vcTexture_Destroy(&pRenderContext->udRenderContext.pColourTex);
  vcTexture_Destroy(&pRenderContext->udRenderContext.pDepthTex);
  vcFramebuffer_Destroy(&pRenderContext->udRenderContext.pFramebuffer);

  UD_ERROR_CHECK(vcTexture_Create(&pRenderContext->udRenderContext.pColourTex, pRenderContext->sceneResolution.x, pRenderContext->sceneResolution.y, pRenderContext->udRenderContext.pColorBuffer, vcTextureFormat_BGRA8, vcTFM_Nearest, false, vcTWM_Clamp, vcTCF_Dynamic));
  UD_ERROR_CHECK(vcTexture_Create(&pRenderContext->udRenderContext.pDepthTex, pRenderContext->sceneResolution.x, pRenderContext->sceneResolution.y, pRenderContext->udRenderContext.pDepthBuffer, vcTextureFormat_D32F, vcTFM_Nearest, false, vcTWM_Clamp, vcTCF_Dynamic));

  for (int i = 0; i < vcRender_RenderBufferCount; ++i)
  {
    vcTexture_Destroy(&pRenderContext->pTexture[i]);
    vcTexture_Destroy(&pRenderContext->pDepthTexture[i]);
    vcFramebuffer_Destroy(&pRenderContext->pFramebuffer[i]);
  }

  vcTexture_Destroy(&pRenderContext->picking.pTexture);
  vcTexture_Destroy(&pRenderContext->picking.pDepth);
  vcFramebuffer_Destroy(&pRenderContext->picking.pFramebuffer);

  for (int i = 0; i < 2; ++i)
  {
    vcTexture_Destroy(&pRenderContext->pAuxiliaryTextures[i]);
    vcFramebuffer_Destroy(&pRenderContext->pAuxiliaryFramebuffers[i]);
  }

  for (int i = 0; i < vcRender_RenderBufferCount; ++i)
  {
    UD_ERROR_CHECK(vcTexture_Create(&pRenderContext->pTexture[i], widthIncr, heightIncr, nullptr, vcTextureFormat_RGBA16F, vcTFM_Linear, false, vcTWM_Clamp, vcTCF_RenderTarget));
    UD_ERROR_CHECK(vcTexture_Create(&pRenderContext->pDepthTexture[i], widthIncr, heightIncr, nullptr, vcTextureFormat_D24S8, vcTFM_Linear, false, vcTWM_Clamp, vcTCF_RenderTarget | vcTCF_AsynchronousRead));
    UD_ERROR_IF(!vcFramebuffer_Create(&pRenderContext->pFramebuffer[i], pRenderContext->pTexture[i], pRenderContext->pDepthTexture[i]), udR_InternalError);
  }

  pRenderContext->effectResolution.x = widthIncr >> vcRender_OutlineEffectDownscale;
  pRenderContext->effectResolution.y = heightIncr >> vcRender_OutlineEffectDownscale;
  pRenderContext->effectResolution.x = pRenderContext->effectResolution.x + (pRenderContext->effectResolution.x % vcRender_SceneSizeIncrement != 0 ? vcRender_SceneSizeIncrement - pRenderContext->effectResolution.x % vcRender_SceneSizeIncrement : 0);
  pRenderContext->effectResolution.y = pRenderContext->effectResolution.y + (pRenderContext->effectResolution.y % vcRender_SceneSizeIncrement != 0 ? vcRender_SceneSizeIncrement - pRenderContext->effectResolution.y % vcRender_SceneSizeIncrement : 0);

  for (int i = 0; i < 2; ++i)
  {
    UD_ERROR_CHECK(vcTexture_Create(&pRenderContext->pAuxiliaryTextures[i], pRenderContext->effectResolution.x, pRenderContext->effectResolution.y, nullptr, vcTextureFormat_RGBA8, vcTFM_Linear, false, vcTWM_Clamp, vcTCF_RenderTarget));
    UD_ERROR_IF(!vcFramebuffer_Create(&pRenderContext->pAuxiliaryFramebuffers[i], pRenderContext->pAuxiliaryTextures[i]), udR_InternalError);
  }

  UD_ERROR_CHECK(vcTexture_Create(&pRenderContext->picking.pTexture, pRenderContext->effectResolution.x, pRenderContext->effectResolution.y, nullptr, vcTextureFormat_BGRA8, vcTFM_Nearest, false, vcTWM_Clamp, vcTCF_RenderTarget));
  UD_ERROR_CHECK(vcTexture_Create(&pRenderContext->picking.pDepth, pRenderContext->effectResolution.x, pRenderContext->effectResolution.y, nullptr, vcTextureFormat_D24S8, vcTFM_Nearest, false, vcTWM_Clamp, vcTCF_RenderTarget));
  UD_ERROR_IF(!vcFramebuffer_Create(&pRenderContext->picking.pFramebuffer, pRenderContext->picking.pTexture, pRenderContext->picking.pDepth), udR_InternalError);

  if (pProgramState->pVDKContext)
    UD_ERROR_CHECK(vcRender_RecreateUDView(pProgramState, pRenderContext));

epilogue:
  return result;
}

// Asychronously read a 1x1 region of last frames depth buffer 
udResult vcRender_AsyncReadFrameDepth(vcRenderContext *pRenderContext)
{
  udResult result = udR_Success;

  if (pRenderContext->currentMouseUV.x < 0 || pRenderContext->currentMouseUV.x > 1 || pRenderContext->currentMouseUV.y < 0 || pRenderContext->currentMouseUV.y > 1)
    return result;

  uint8_t depthBytes[4] = {};
  udUInt2 pickLocation = { (uint32_t)(pRenderContext->currentMouseUV.x * pRenderContext->sceneResolution.x), (uint32_t)(pRenderContext->currentMouseUV.y * pRenderContext->sceneResolution.y) };
#if GRAPHICS_API_OPENGL
  pickLocation.y = pRenderContext->sceneResolution.y - pickLocation.y - 1; // upside-down
#endif

  static const int readBufferIndex = 0;
  UD_ERROR_IF(!vcTexture_EndReadPixels(pRenderContext->pDepthTexture[readBufferIndex], pickLocation.x, pickLocation.y, 1, 1, depthBytes), udR_InternalError); // read previous copy
  UD_ERROR_IF(!vcTexture_BeginReadPixels(pRenderContext->pDepthTexture[readBufferIndex], pickLocation.x, pickLocation.y, 1, 1, depthBytes, pRenderContext->pFramebuffer[readBufferIndex]), udR_InternalError); // begin copy for next frame read

  // 24 bit unsigned int -> float
#if GRAPHICS_API_OPENGL || GRAPHICS_API_METAL
  pRenderContext->previousFrameDepth = uint32_t((depthBytes[3] << 16) | (depthBytes[2] << 8) | (depthBytes[1] << 0)) / ((1 << 24) - 1.0f);
  //uint8_t stencil = depthBytes[0];
#else
  pRenderContext->previousFrameDepth = uint32_t((depthBytes[2] << 16) | (depthBytes[1] << 8) | (depthBytes[0] << 0)) / ((1 << 24) - 1.0f);
  //uint8_t stencil = depthBytes[3];
#endif

  // fbo state may not be valid (e.g. first read back will be '0')
  if (pRenderContext->previousFrameDepth == 0.0f)
    pRenderContext->previousFrameDepth = 1.0f;

epilogue:
  return result;
}

void vcRenderSkybox(vcState *pProgramState, vcRenderContext *pRenderContext)
{
  // Draw the skybox only at the far plane, where there is no geometry.
  vcGLState_SetDepthStencilMode(vcGLSDM_LessOrEqual, false);

  if (pProgramState->settings.presentation.showSkybox)
  {
    udFloat4x4 viewMatrixF = udFloat4x4::create(pProgramState->camera.matrices.view);
    udFloat4x4 projectionMatrixF = udFloat4x4::create(pProgramState->camera.matrices.projectionNear);
    udFloat4x4 inverseViewProjMatrixF = projectionMatrixF * viewMatrixF;
    inverseViewProjMatrixF.axis.t = udFloat4::create(0, 0, 0, 1);
    inverseViewProjMatrixF.inverse();

    vcShader_Bind(pRenderContext->skyboxShaderPanorama.pProgram);
    vcShader_BindTexture(pRenderContext->skyboxShaderPanorama.pProgram, pRenderContext->skyboxShaderPanorama.pSkyboxTexture, 0, pRenderContext->skyboxShaderPanorama.uniform_texture);
    vcShader_BindConstantBuffer(pRenderContext->skyboxShaderPanorama.pProgram, pRenderContext->skyboxShaderPanorama.uniform_MatrixBlock, &inverseViewProjMatrixF, sizeof(inverseViewProjMatrixF));
  }
  else
  {
    pRenderContext->skyboxShaderTintImage.params.colour = pProgramState->settings.presentation.skyboxColour;

    int x = 0;
    int y = 0;

    vcShader_Bind(pRenderContext->skyboxShaderTintImage.pProgram);

    if (vcTexture_GetSize(pRenderContext->skyboxShaderTintImage.pLogoTexture, &x, &y) == udR_Success)
    {
      pRenderContext->skyboxShaderTintImage.params.imageSize.x = (float)x / pRenderContext->sceneResolution.x;
      pRenderContext->skyboxShaderTintImage.params.imageSize.y = -(float)y / pRenderContext->sceneResolution.y;

      vcShader_BindTexture(pRenderContext->skyboxShaderTintImage.pProgram, pRenderContext->skyboxShaderTintImage.pLogoTexture, 0, pRenderContext->skyboxShaderTintImage.uniform_texture);
    }
    else
    {
      pRenderContext->skyboxShaderTintImage.params.colour.w = 0.0;
    }

    vcShader_BindConstantBuffer(pRenderContext->skyboxShaderTintImage.pProgram, pRenderContext->skyboxShaderTintImage.uniform_params, &pRenderContext->skyboxShaderTintImage.params, sizeof(pRenderContext->skyboxShaderTintImage.params));
  }

  vcGLState_SetViewportDepthRange(1.0f, 1.0f);

  vcMesh_Render(gInternalMeshes[vcInternalMeshType_ScreenQuad]);

  vcGLState_SetViewportDepthRange(0.0f, 1.0f);

  vcShader_Bind(nullptr);
}

void vcRender_SplatUDWithId(vcState *pProgramState, vcRenderContext *pRenderContext, float id)
{
  udUnused(pProgramState); // Some configurations no longer use this parameter

  vcShader_Bind(pRenderContext->udRenderContext.splatIdShader.pProgram);

  vcShader_BindTexture(pRenderContext->udRenderContext.splatIdShader.pProgram, pRenderContext->udRenderContext.pColourTex, 0, pRenderContext->udRenderContext.splatIdShader.uniform_texture);
  vcShader_BindTexture(pRenderContext->udRenderContext.splatIdShader.pProgram, pRenderContext->udRenderContext.pDepthTex, 1, pRenderContext->udRenderContext.splatIdShader.uniform_depth);

  pRenderContext->udRenderContext.splatIdShader.params.id = udFloat4::create(0.0f, 0.0f, 0.0f, id);
  vcShader_BindConstantBuffer(pRenderContext->udRenderContext.splatIdShader.pProgram, pRenderContext->udRenderContext.splatIdShader.uniform_params, &pRenderContext->udRenderContext.splatIdShader.params, sizeof(pRenderContext->udRenderContext.splatIdShader.params));

  vcMesh_Render(gInternalMeshes[vcInternalMeshType_ScreenQuad]);
}

void vcRender_SplatUD(vcState *pProgramState, vcRenderContext *pRenderContext, vcTexture *pColour, vcTexture *pDepth)
{
  udUnused(pProgramState);

  vcShader_Bind(pRenderContext->udRenderContext.presentShader.pProgram);

  vcShader_BindTexture(pRenderContext->udRenderContext.presentShader.pProgram, pColour, 0, pRenderContext->udRenderContext.presentShader.uniform_texture);
  vcShader_BindTexture(pRenderContext->udRenderContext.presentShader.pProgram, pDepth, 1, pRenderContext->udRenderContext.presentShader.uniform_depth);

  vcMesh_Render(gInternalMeshes[vcInternalMeshType_ScreenQuad]);
}

void vcRenderTerrain(vcState *pProgramState, vcRenderContext *pRenderContext)
{
  if (pProgramState->gis.isProjected && pProgramState->settings.maptiles.mapEnabled)
  {
    udDouble4x4 cameraMatrix = pProgramState->camera.matrices.camera;
    udDouble4x4 viewProjection = pProgramState->camera.matrices.viewProjection;

#ifndef GIT_BUILD
    static bool debugDetachCamera = false;
    static udDouble4x4 gRealCameraMatrix = udDouble4x4::identity();
    if (!debugDetachCamera)
      gRealCameraMatrix = pProgramState->camera.matrices.camera;

    cameraMatrix = gRealCameraMatrix;
    viewProjection = pProgramState->camera.matrices.projection * udInverse(cameraMatrix);
#endif
    udDouble3 localCamPos = cameraMatrix.axis.t.toVector3();

    // Corners [nw, ne, sw, se]
    udDouble3 localCorners[4];
    udInt2 slippyCorners[4];

    int currentZoom = 21;

    double farPlane = pProgramState->settings.camera.farPlane;

    // Cardinal Limits
    localCorners[0] = localCamPos + udDouble3::create(-farPlane, +farPlane, 0);
    localCorners[1] = localCamPos + udDouble3::create(+farPlane, +farPlane, 0);
    localCorners[2] = localCamPos + udDouble3::create(-farPlane, -farPlane, 0);
    localCorners[3] = localCamPos + udDouble3::create(+farPlane, -farPlane, 0);

    for (int i = 0; i < 4; ++i)
      vcGIS_LocalToSlippy(&pProgramState->gis, &slippyCorners[i], localCorners[i], currentZoom);

    while (currentZoom > 0 && (slippyCorners[0] != slippyCorners[1] || slippyCorners[1] != slippyCorners[2] || slippyCorners[2] != slippyCorners[3]))
    {
      --currentZoom;

      for (int i = 0; i < 4; ++i)
        slippyCorners[i] /= 2;
    }

    for (int i = 0; i < 4; ++i)
      vcGIS_SlippyToLocal(&pProgramState->gis, &localCorners[i], slippyCorners[0] + udInt2::create(i & 1, i / 2), currentZoom);

    vcTileRenderer_Update(pRenderContext->pTileRenderer, pProgramState->deltaTime, &pProgramState->gis, localCorners, udInt3::create(slippyCorners[0], currentZoom), localCamPos, viewProjection);
    vcTileRenderer_Render(pRenderContext->pTileRenderer, pProgramState->camera.matrices.view, pProgramState->camera.matrices.projection);
  }
}

void vcRender_FXAAPass(vcState *pProgramState, vcRenderContext *pRenderContext)
{
  udUnused(pProgramState);

  if (!pProgramState->settings.presentation.antiAliasingOn)
    return;

  vcGLState_SetBlendMode(vcGLSBM_None);
  vcGLState_SetFaceMode(vcGLSFM_Solid, vcGLSCM_Back);
  vcGLState_SetDepthStencilMode(vcGLSDM_Always, false);

  pRenderContext->activeRenderTarget = 1 - pRenderContext->activeRenderTarget;
  vcFramebuffer_Bind(pRenderContext->pFramebuffer[pRenderContext->activeRenderTarget], vcFramebufferClearOperation_All, 0x00FF8080);

  vcShader_Bind(pRenderContext->fxaaShader.pProgram);
  vcShader_BindTexture(pRenderContext->fxaaShader.pProgram, pRenderContext->pTexture[1 - pRenderContext->activeRenderTarget], 0, pRenderContext->fxaaShader.uniform_texture);
  vcShader_BindTexture(pRenderContext->fxaaShader.pProgram, pRenderContext->pDepthTexture[1 - pRenderContext->activeRenderTarget], 1, pRenderContext->fxaaShader.uniform_depth);

  pRenderContext->fxaaShader.params.screenParams.x = (1.0f / pRenderContext->sceneResolution.x);
  pRenderContext->fxaaShader.params.screenParams.y = (1.0f / pRenderContext->sceneResolution.y);

  vcShader_BindConstantBuffer(pRenderContext->fxaaShader.pProgram, pRenderContext->fxaaShader.uniform_params, &pRenderContext->fxaaShader.params, sizeof(pRenderContext->fxaaShader.params));

  vcMesh_Render(gInternalMeshes[vcInternalMeshType_ScreenQuad]);
}

void vcRender_VisualizationPass(vcState *pProgramState, vcRenderContext *pRenderContext)
{
  vcGLState_SetDepthStencilMode(vcGLSDM_Always, true);

  pRenderContext->activeRenderTarget = 1 - pRenderContext->activeRenderTarget;
  vcFramebuffer_Bind(pRenderContext->pFramebuffer[pRenderContext->activeRenderTarget], vcFramebufferClearOperation_All, 0x00FF8080);

  vcShader_Bind(pRenderContext->visualizationShader.pProgram);
  vcShader_BindTexture(pRenderContext->visualizationShader.pProgram, pRenderContext->pTexture[1 - pRenderContext->activeRenderTarget], 0, pRenderContext->visualizationShader.uniform_texture);
  vcShader_BindTexture(pRenderContext->visualizationShader.pProgram, pRenderContext->pDepthTexture[1 - pRenderContext->activeRenderTarget], 1, pRenderContext->visualizationShader.uniform_depth);

  float nearPlane = pProgramState->settings.camera.nearPlane;
  float farPlane = pProgramState->settings.camera.farPlane;

  // edge outlines
  int outlineWidth = pProgramState->settings.postVisualization.edgeOutlines.width;
  float outlineEdgeThreshold = pProgramState->settings.postVisualization.edgeOutlines.threshold;
  udFloat4 outlineColour = pProgramState->settings.postVisualization.edgeOutlines.colour;
  if (!pProgramState->settings.postVisualization.edgeOutlines.enable)
    outlineColour.w = 0.0f;

  // colour by height
  udFloat4 colourByHeightMinColour = pProgramState->settings.postVisualization.colourByHeight.minColour;
  if (!pProgramState->settings.postVisualization.colourByHeight.enable)
    colourByHeightMinColour.w = 0.f;
  udFloat4 colourByHeightMaxColour = pProgramState->settings.postVisualization.colourByHeight.maxColour;
  if (!pProgramState->settings.postVisualization.colourByHeight.enable)
    colourByHeightMaxColour.w = 0.f;
  float colourByHeightStartHeight = pProgramState->settings.postVisualization.colourByHeight.startHeight;
  float colourByHeightEndHeight = pProgramState->settings.postVisualization.colourByHeight.endHeight;

  // colour by depth
  udFloat4 colourByDepthColour = pProgramState->settings.postVisualization.colourByDepth.colour;
  if (!pProgramState->settings.postVisualization.colourByDepth.enable)
    colourByDepthColour.w = 0.f;
  float colourByDepthStart = pProgramState->settings.postVisualization.colourByDepth.startDepth;
  float colourByDepthEnd = pProgramState->settings.postVisualization.colourByDepth.endDepth;

  // contours
  udFloat4 contourColour = pProgramState->settings.postVisualization.contours.colour;
  float contourDistances = pProgramState->settings.postVisualization.contours.distances;
  float contourBandHeight = pProgramState->settings.postVisualization.contours.bandHeight;
  float contourRainboxRepeatRate = pProgramState->settings.postVisualization.contours.rainbowRepeat;
  float contourRainboxIntensity = pProgramState->settings.postVisualization.contours.rainbowIntensity;

  if (!pProgramState->settings.postVisualization.contours.enable)
  {
    contourColour.w = 0.f;
    contourRainboxIntensity = 0.f;
  }

  pRenderContext->visualizationShader.fragParams.inverseViewProjection = udFloat4x4::create(pProgramState->camera.matrices.inverseViewProjection);
  pRenderContext->visualizationShader.fragParams.inverseProjection = udFloat4x4::create(udInverse(pProgramState->camera.matrices.projection));
  pRenderContext->visualizationShader.fragParams.screenParams.x = (1.0f / pRenderContext->sceneResolution.x);
  pRenderContext->visualizationShader.fragParams.screenParams.y = (1.0f / pRenderContext->sceneResolution.y);
  pRenderContext->visualizationShader.fragParams.screenParams.z = nearPlane;
  pRenderContext->visualizationShader.fragParams.screenParams.w = farPlane;
  pRenderContext->visualizationShader.fragParams.outlineColour = outlineColour;
  pRenderContext->visualizationShader.fragParams.outlineParams.x = (float)outlineWidth;
  pRenderContext->visualizationShader.fragParams.outlineParams.y = outlineEdgeThreshold;
  pRenderContext->visualizationShader.fragParams.colourizeHeightColourMin = colourByHeightMinColour;
  pRenderContext->visualizationShader.fragParams.colourizeHeightColourMax = colourByHeightMaxColour;
  pRenderContext->visualizationShader.fragParams.colourizeHeightParams.x = colourByHeightStartHeight;
  pRenderContext->visualizationShader.fragParams.colourizeHeightParams.y = colourByHeightEndHeight;
  pRenderContext->visualizationShader.fragParams.colourizeDepthColour = colourByDepthColour;
  pRenderContext->visualizationShader.fragParams.colourizeDepthParams.x = colourByDepthStart;
  pRenderContext->visualizationShader.fragParams.colourizeDepthParams.y = colourByDepthEnd;
  pRenderContext->visualizationShader.fragParams.contourColour = contourColour;
  pRenderContext->visualizationShader.fragParams.contourParams.x = contourDistances;
  pRenderContext->visualizationShader.fragParams.contourParams.y = contourBandHeight;
  pRenderContext->visualizationShader.fragParams.contourParams.z = contourRainboxRepeatRate;
  pRenderContext->visualizationShader.fragParams.contourParams.w = contourRainboxIntensity;

  pRenderContext->visualizationShader.vertParams.outlineStepSize.x = outlineWidth * (1.0f / pRenderContext->sceneResolution.x);
  pRenderContext->visualizationShader.vertParams.outlineStepSize.y = outlineWidth * (1.0f / pRenderContext->sceneResolution.y);

  vcShader_BindConstantBuffer(pRenderContext->visualizationShader.pProgram, pRenderContext->visualizationShader.uniform_vertParams, &pRenderContext->visualizationShader.vertParams, sizeof(pRenderContext->visualizationShader.vertParams));
  vcShader_BindConstantBuffer(pRenderContext->visualizationShader.pProgram, pRenderContext->visualizationShader.uniform_fragParams, &pRenderContext->visualizationShader.fragParams, sizeof(pRenderContext->visualizationShader.fragParams));

  vcMesh_Render(gInternalMeshes[vcInternalMeshType_ScreenQuad]);
}

void vcRender_ApplyViewShed(vcRenderContext *pRenderContext)
{
  vcShader_Bind(pRenderContext->shadowShader.pProgram);
  vcShader_BindTexture(pRenderContext->shadowShader.pProgram, pRenderContext->pDepthTexture[0], 0, pRenderContext->shadowShader.uniform_depth);
  vcShader_BindTexture(pRenderContext->shadowShader.pProgram, pRenderContext->viewShedRenderingContext.pDepthTex, 1, pRenderContext->shadowShader.uniform_shadowMapAtlas);

  vcShader_BindConstantBuffer(pRenderContext->udRenderContext.presentShader.pProgram, pRenderContext->shadowShader.uniform_params, &pRenderContext->shadowShader.params, sizeof(pRenderContext->shadowShader.params));

  vcMesh_Render(gInternalMeshes[vcInternalMeshType_ScreenQuad]);
}

void vcRender_RenderAndApplyViewSheds(vcState *pProgramState, vcRenderContext *pRenderContext, vcRenderData &renderData)
{
  if (renderData.viewSheds.length == 0)
    return;

  udUInt2 singleRenderSize = udUInt2::create(ViewShedMapRes.x / ViewShedMapCount, ViewShedMapRes.y);
  udFloat2 atlasSize = udFloat2::create((float)ViewShedMapRes.x, (float)ViewShedMapRes.y);

  // TODO: aabb vs. point test to render only affecting geometry
  bool doUDRender = renderData.models.length > 0;
  bool doPolygonRender = false;
  for (size_t p = 0; p < renderData.polyModels.length && !doPolygonRender; ++p)
    doPolygonRender = !renderData.polyModels[p].HasFlag(vcRenderPolyInstance::RenderFlags_Transparent);

  if (!doUDRender && !doPolygonRender)
    return;

  if (pRenderContext->viewShedRenderingContext.pRenderView == nullptr)
    vdkRenderView_Create(pProgramState->pVDKContext, &pRenderContext->viewShedRenderingContext.pRenderView, pRenderContext->udRenderContext.pRenderer, singleRenderSize.x, singleRenderSize.y);

  for (size_t v = 0; v < renderData.viewSheds.length; ++v)
  {
    vcViewShedData *pViewShedData = &renderData.viewSheds[v];

    vcCameraSettings cameraSettings = {};
    cameraSettings.nearPlane = pViewShedData->nearFarPlane.x;
    cameraSettings.farPlane = pViewShedData->nearFarPlane.y;
    cameraSettings.fieldOfView = pViewShedData->fieldOfView;

    // set up cameras for these renders
    vcCamera shadowRenderCameras[ViewShedMapCount] = {};
    for (int r = 0; r < ViewShedMapCount; ++r)
    {
      shadowRenderCameras[r].position = pViewShedData->position;

      double rot = (UD_DEG2RAD(360.0) / ViewShedMapCount) * r;
      shadowRenderCameras[r].eulerRotation = udDouble3::create(-rot, 0, 0);
      vcCamera_UpdateMatrices(&shadowRenderCameras[r], cameraSettings, atlasSize, nullptr);

      pRenderContext->shadowShader.params.shadowMapVP[r] = udFloat4x4::create(shadowRenderCameras[r].matrices.projectionUD * (shadowRenderCameras[r].matrices.view * udInverse(pProgramState->camera.matrices.view)));
    }

    // Texture uploads first (Unlimited Detail)
    if (doUDRender)
    {
      for (int r = 0; r < ViewShedMapCount; ++r)
      {
        // configure UD render to only render into portion of buffer
        vdkRenderView_SetTargetsWithPitch(pRenderContext->viewShedRenderingContext.pRenderView, nullptr, 0, pRenderContext->viewShedRenderingContext.pDepthBuffer + r * singleRenderSize.x, 0, ViewShedMapRes.x * 4);
        vdkRenderView_SetMatrix(pRenderContext->viewShedRenderingContext.pRenderView, vdkRVM_Projection, shadowRenderCameras[r].matrices.projectionUD.a);

        // render UD
        vcRender_RenderUD(pProgramState, pRenderContext, pRenderContext->viewShedRenderingContext.pRenderView, &shadowRenderCameras[r], renderData, false);
      }

      vcTexture_UploadPixels(pRenderContext->viewShedRenderingContext.pUDDepthTexture, pRenderContext->viewShedRenderingContext.pDepthBuffer, ViewShedMapRes.x, ViewShedMapRes.y);
    }

    vcGLState_SetBlendMode(vcGLSBM_None);
    vcGLState_SetFaceMode(vcGLSFM_Solid, vcGLSCM_None);
    vcGLState_SetDepthStencilMode(vcGLSDM_LessOrEqual, true);

    vcGLState_SetViewport(0, 0, ViewShedMapRes.x, ViewShedMapRes.y);
    vcFramebuffer_Bind(pRenderContext->viewShedRenderingContext.pFramebuffer, vcFramebufferClearOperation_All);

    if (doUDRender)
      vcRender_SplatUD(pProgramState, pRenderContext, pRenderContext->viewShedRenderingContext.pDummyColour, pRenderContext->viewShedRenderingContext.pUDDepthTexture);

    if (doPolygonRender)
    {
      for (int r = 0; r < ViewShedMapCount; ++r)
      {
        udDouble4x4 viewProjection = shadowRenderCameras[r].matrices.projection * shadowRenderCameras[r].matrices.view;
        vcGLState_SetViewport(r * singleRenderSize.x, 0, singleRenderSize.x, ViewShedMapRes.y);

        for (size_t p = 0; p < renderData.polyModels.length; ++p)
        {
          vcRenderPolyInstance *pInstance = &renderData.polyModels[p];
          if (pInstance->HasFlag(vcRenderPolyInstance::RenderFlags_Transparent))
            continue;

          if (pInstance->renderType == vcRenderPolyInstance::RenderType_Polygon)
            vcPolygonModel_Render(pInstance->pModel, pInstance->worldMat, viewProjection, vcPMP_Shadows);
          else if (pInstance->renderType == vcRenderPolyInstance::RenderType_SceneLayer)
            vcSceneLayerRenderer_Render(pInstance->pSceneLayer, pInstance->worldMat, viewProjection, shadowRenderCameras[r].position, ViewShedMapRes, nullptr, true);
        }
      }
    }

    pRenderContext->shadowShader.params.inverseProjection = udFloat4x4::create(udInverse(pProgramState->camera.matrices.projection));
    pRenderContext->shadowShader.params.nearFarPlane = udFloat4::create(pViewShedData->nearFarPlane.x, pViewShedData->nearFarPlane.y, 0.0f, 0.0f);
    pRenderContext->shadowShader.params.visibleColour = pViewShedData->visibleColour;
    pRenderContext->shadowShader.params.notVisibleColour = pViewShedData->notVisibleColour;

    vcGLState_SetDepthStencilMode(vcGLSDM_Always, false);
    vcGLState_SetFaceMode(vcGLSFM_Solid, vcGLSCM_Back);
    vcGLState_SetBlendMode(vcGLSBM_Additive);

    vcFramebuffer_Bind(pRenderContext->pFramebuffer[1]); // assumed this is the working target
    vcGLState_SetViewport(0, 0, pRenderContext->sceneResolution.x, pRenderContext->sceneResolution.y);
    vcRender_ApplyViewShed(pRenderContext);
  }

  vcGLState_ResetState();
}

void vcRender_OpaquePass(vcState *pProgramState, vcRenderContext *pRenderContext, vcRenderData &renderData)
{
  vcFramebuffer_Bind(pRenderContext->pFramebuffer[pRenderContext->activeRenderTarget], vcFramebufferClearOperation_All, 0xFFFF8080);

  vcGLState_ResetState();

  // UD
  vcRender_SplatUD(pProgramState, pRenderContext, pRenderContext->udRenderContext.pColourTex, pRenderContext->udRenderContext.pDepthTex);

  // Polygon Models
  {
    vcGLState_SetBlendMode(vcGLSBM_None);
    vcGLState_SetFaceMode(vcGLSFM_Solid, vcGLSCM_Back);
    vcGLState_SetDepthStencilMode(vcGLSDM_LessOrEqual, true);

    vcSceneLayer_BeginFrame();

    udFloat4 whiteColour = udFloat4::one();
    for (size_t i = 0; i < renderData.polyModels.length; ++i)
    {
      vcRenderPolyInstance *pInstance = &renderData.polyModels[i];
      if (pInstance->HasFlag(vcRenderPolyInstance::RenderFlags_Transparent))
        continue;

      udFloat4 *pTintOverride = nullptr;
      if (pInstance->HasFlag(vcRenderPolyInstance::RenderFlags_IgnoreTint))
        pTintOverride = &whiteColour;

      vcGLState_SetFaceMode(vcGLSFM_Solid, pInstance->cullFace);

      if (pInstance->renderType == vcRenderPolyInstance::RenderType_Polygon)
        vcPolygonModel_Render(pInstance->pModel, pInstance->worldMat, pProgramState->camera.matrices.viewProjection, vcPMP_Standard, pInstance->pDiffuseOverride, pTintOverride);
      else if (pInstance->renderType == vcRenderPolyInstance::RenderType_SceneLayer)
        vcSceneLayerRenderer_Render(pInstance->pSceneLayer, pInstance->worldMat, pProgramState->camera.matrices.viewProjection, pProgramState->camera.position, pRenderContext->sceneResolution);
    }

    vcGLState_SetFaceMode(vcGLSFM_Solid, vcGLSCM_Back);
    vcSceneLayer_EndFrame();

    for (size_t i = 0; i < renderData.waterVolumes.length; ++i)
      vcWaterRenderer_Render(renderData.waterVolumes[i], pProgramState->camera.matrices.view, pProgramState->camera.matrices.viewProjection, pRenderContext->skyboxShaderPanorama.pSkyboxTexture, pProgramState->deltaTime);
  }

  vcRender_AsyncReadFrameDepth(pRenderContext); // note: one frame behind
}

void vcRender_TransparentPass(vcState *pProgramState, vcRenderContext *pRenderContext, vcRenderData &renderData)
{
  vcGLState_SetBlendMode(vcGLSBM_Interpolative);
  vcGLState_SetDepthStencilMode(vcGLSDM_LessOrEqual, false);

  // Images
  {
    vcGLState_SetFaceMode(vcGLSFM_Solid, vcGLSCM_Front);

    for (size_t i = 0; i < renderData.images.length; ++i)
    {
      static const double distScalar = 1000.0; // Param

      double zScale = 1.0;
      zScale -= udMag3(pProgramState->camera.position - renderData.images[i]->position) / distScalar;

      if (zScale < 0) // too far
        continue;

      vcImageRenderer_Render(renderData.images[i], pProgramState->camera.matrices.viewProjection, pRenderContext->sceneResolution, zScale);
    }
  }

  // Fences
  {
    vcGLState_SetFaceMode(vcGLSFM_Solid, vcGLSCM_None);
    if (pProgramState->settings.presentation.showDiagnosticInfo)
      vcFenceRenderer_Render(pRenderContext->pDiagnosticFences, pProgramState->camera.matrices.viewProjection, pProgramState->deltaTime);

    for (size_t i = 0; i < renderData.fences.length; ++i)
      vcFenceRenderer_Render(renderData.fences[i], pProgramState->camera.matrices.viewProjection, pProgramState->deltaTime);
  }

  udFloat4 transparentColour = udFloat4::create(1, 1, 1, 0.65f);
  for (size_t i = 0; i < renderData.polyModels.length; ++i)
  {
    vcRenderPolyInstance *pInstance = &renderData.polyModels[i];
    if (!pInstance->HasFlag(vcRenderPolyInstance::RenderFlags_Transparent))
      continue;

    vcGLState_SetFaceMode(vcGLSFM_Solid, pInstance->cullFace);

    if (pInstance->renderType == vcRenderPolyInstance::RenderType_Polygon)
      vcPolygonModel_Render(pInstance->pModel, pInstance->worldMat, pProgramState->camera.matrices.viewProjection, vcPMP_Standard, pInstance->pDiffuseOverride, &transparentColour);
    else if (pInstance->renderType == vcRenderPolyInstance::RenderType_SceneLayer)
      vcSceneLayerRenderer_Render(pInstance->pSceneLayer, pInstance->worldMat, pProgramState->camera.matrices.viewProjection, pProgramState->camera.position, pRenderContext->sceneResolution);
  }

  vcGLState_SetFaceMode(vcGLSFM_Solid, vcGLSCM_Back);
}

void vcRender_BeginFrame(vcState *pProgramState, vcRenderContext *pRenderContext, vcRenderData &renderData)
{
  udUnused(pProgramState);
  udUnused(pRenderContext);

  // Would be nice to use 'pRenderContext->activeRenderTarget' here, but this causes
  // a single frame 'flicker' if options are changed at run time.
  renderData.pSceneTexture = pRenderContext->pTexture[0];//pProgramState->settings.presentation.antiAliasingOn ? 0 : 1];
  renderData.sceneScaling = udFloat2::one();

  pRenderContext->activeRenderTarget = 0;

  // TODO (EVC-835): fix scene scaling
  // udFloat2::create(float(pRenderContext->originalSceneResolution.x) / pRenderContext->sceneResolution.x, float(pRenderContext->originalSceneResolution.y) / pRenderContext->sceneResolution.y);
}

void vcRender_ApplySelectionBuffer(vcState *pProgramState, vcRenderContext *pRenderContext)
{
  udFloat2 sampleStepSize = udFloat2::create(1.0f / pRenderContext->effectResolution.x, 1.0f / pRenderContext->effectResolution.y);

  vcGLState_SetBlendMode(vcGLSBM_Interpolative);
  vcGLState_SetDepthStencilMode(vcGLSDM_Always, false);

  pRenderContext->selectionShader.params.stepSizeThickness.x = sampleStepSize.x;
  pRenderContext->selectionShader.params.stepSizeThickness.y = sampleStepSize.y;
  pRenderContext->selectionShader.params.stepSizeThickness.z = (float)(2 << vcRender_OutlineEffectDownscale) * (pProgramState->settings.objectHighlighting.thickness - 1.0f); // roughly
  pRenderContext->selectionShader.params.stepSizeThickness.w = 0.2f;
  pRenderContext->selectionShader.params.colour = pProgramState->settings.objectHighlighting.colour;

  vcShader_Bind(pRenderContext->selectionShader.pProgram);

  vcShader_BindTexture(pRenderContext->selectionShader.pProgram, pRenderContext->pAuxiliaryTextures[0], 0, pRenderContext->selectionShader.uniform_texture);
  vcShader_BindConstantBuffer(pRenderContext->selectionShader.pProgram, pRenderContext->selectionShader.uniform_params, &pRenderContext->selectionShader.params, sizeof(pRenderContext->selectionShader.params));

  vcMesh_Render(gInternalMeshes[vcInternalMeshType_ScreenQuad]);
}

udFloat4 vcRender_EncodeIdAsColour(uint32_t id)
{
  return udFloat4::create(0.0f, ((id & 0xff) / 255.0f), ((id & 0xff00) >> 8) / 255.0f, 1.0f);// ((id & 0xff0000) >> 16) / 255.0f);// ((id & 0xff000000) >> 24) / 255.0f);
}

bool vcRender_DrawSelectedGeometry(vcState *pProgramState, vcRenderContext *pRenderContext, vcRenderData &renderData)
{
  bool active = false;

  // check UD first
  uint32_t modelIndex = 0; // index is based on certain models
  for (size_t i = 0; i < renderData.models.length; ++i)
  {
    if (renderData.models[i]->m_visible && renderData.models[i]->m_loadStatus == vcSLS_Loaded)
    {
      if (renderData.models[i]->IsSubitemSelected(0))
      {
        float splatId = 1.0f / 255.0f;

        vcRender_SplatUDWithId(pProgramState, pRenderContext, splatId);
        active = true;
      }
      ++modelIndex;
    }
  }

  udFloat4 selectionMask = udFloat4::create(1.0f); // mask selected object
  for (size_t i = 0; i < renderData.polyModels.length; ++i)
  {
    vcRenderPolyInstance *pInstance = &renderData.polyModels[i];
    if ((pInstance->sceneItemInternalId == 0 && pInstance->pSceneItem->m_selected) || (pInstance->sceneItemInternalId != 0 && pInstance->pSceneItem->IsSubitemSelected(pInstance->sceneItemInternalId)))
    {
      vcGLState_SetFaceMode(vcGLSFM_Solid, pInstance->cullFace);

      if (pInstance->renderType == vcRenderPolyInstance::RenderType_Polygon)
        vcPolygonModel_Render(pInstance->pModel, pInstance->worldMat, pProgramState->camera.matrices.viewProjection, vcPMP_ColourOnly, nullptr, &selectionMask);
      else if (pInstance->renderType == vcRenderPolyInstance::RenderType_SceneLayer)
        vcSceneLayerRenderer_Render(pInstance->pSceneLayer, pInstance->worldMat, pProgramState->camera.matrices.viewProjection, pProgramState->camera.position, pRenderContext->sceneResolution, &selectionMask);

      active = true;
    }

  }

  vcGLState_SetFaceMode(vcGLSFM_Solid, vcGLSCM_Back);

  return active;
}

bool vcRender_CreateSelectionBuffer(vcState *pProgramState, vcRenderContext *pRenderContext, vcRenderData &renderData)
{
  if (pProgramState->settings.objectHighlighting.colour.w == 0.0f) // disabled
    return false;

  vcGLState_SetDepthStencilMode(vcGLSDM_Always, false);
  vcGLState_SetFaceMode(vcGLSFM_Solid, vcGLSCM_Back);
  vcGLState_SetBlendMode(vcGLSBM_None);
  vcGLState_SetViewport(0, 0, pRenderContext->effectResolution.x, pRenderContext->effectResolution.y);

  vcFramebuffer_Bind(pRenderContext->pAuxiliaryFramebuffers[0], vcFramebufferClearOperation_All);

  if (!vcRender_DrawSelectedGeometry(pProgramState, pRenderContext, renderData))
    return false;

  // blur disabled at thickness value of 1.0
  if (pProgramState->settings.objectHighlighting.thickness > 1.0f)
  {
    udFloat2 sampleStepSize = udFloat2::create(1.0f / pRenderContext->effectResolution.x, 1.0f / pRenderContext->effectResolution.y);

    vcGLState_SetBlendMode(vcGLSBM_None);
    for (int i = 0; i < 2; ++i)
    {
      vcFramebuffer_Bind(pRenderContext->pAuxiliaryFramebuffers[1 - i], vcFramebufferClearOperation_All);

      pRenderContext->blurShader.params.stepSize.x = i == 0 ? sampleStepSize.x : 0.0f;
      pRenderContext->blurShader.params.stepSize.y = i == 0 ? 0.0f : sampleStepSize.y;

      vcShader_Bind(pRenderContext->blurShader.pProgram);

      vcShader_BindTexture(pRenderContext->blurShader.pProgram, pRenderContext->pAuxiliaryTextures[i], 0, pRenderContext->blurShader.uniform_texture);
      vcShader_BindConstantBuffer(pRenderContext->blurShader.pProgram, pRenderContext->blurShader.uniform_params, &pRenderContext->blurShader.params, sizeof(pRenderContext->blurShader.params));

      vcMesh_Render(gInternalMeshes[vcInternalMeshType_ScreenQuad]);

      vcShader_BindTexture(pRenderContext->blurShader.pProgram, nullptr, 0, pRenderContext->blurShader.uniform_texture);
    }
  }
  return true;
}

void vcRender_RenderScene(vcState *pProgramState, vcRenderContext *pRenderContext, vcRenderData &renderData, vcFramebuffer *pDefaultFramebuffer)
{
  udUnused(pDefaultFramebuffer);

  udUnused(pDefaultFramebuffer);

  float aspect = pRenderContext->sceneResolution.x / (float)pRenderContext->sceneResolution.y;

  // Render and upload UD buffers
  {
    vcRender_RenderUD(pProgramState, pRenderContext, pRenderContext->udRenderContext.pRenderView, &pProgramState->camera, renderData, true);
    vcTexture_UploadPixels(pRenderContext->udRenderContext.pColourTex, pRenderContext->udRenderContext.pColorBuffer, pRenderContext->sceneResolution.x, pRenderContext->sceneResolution.y);
    vcTexture_UploadPixels(pRenderContext->udRenderContext.pDepthTex, pRenderContext->udRenderContext.pDepthBuffer, pRenderContext->sceneResolution.x, pRenderContext->sceneResolution.y);
  }
  
  bool selectionBufferActive = vcRender_CreateSelectionBuffer(pProgramState, pRenderContext, renderData);
  
  vcGLState_SetDepthStencilMode(vcGLSDM_LessOrEqual, true);
  vcGLState_SetViewport(0, 0, pRenderContext->sceneResolution.x, pRenderContext->sceneResolution.y);
  
  vcRender_OpaquePass(pProgramState, pRenderContext, renderData); // first pass
  vcRender_VisualizationPass(pProgramState, pRenderContext);

  vcGLState_SetViewport(0, 0, pRenderContext->sceneResolution.x, pRenderContext->sceneResolution.y);

  vcFramebuffer_Bind(pRenderContext->pFramebuffer[pRenderContext->activeRenderTarget]);

  vcRender_RenderAndApplyViewSheds(pProgramState, pRenderContext, renderData);

  vcRenderSkybox(pProgramState, pRenderContext); // Drawing skybox after opaque geometry saves a bit on fill rate.
  vcRenderTerrain(pProgramState, pRenderContext);
  vcRender_TransparentPass(pProgramState, pRenderContext, renderData);


  pRenderContext->activeRenderTarget = 1 - pRenderContext->activeRenderTarget;
  vcFramebuffer_Bind(pRenderContext->pFramebuffer[pRenderContext->activeRenderTarget], vcFramebufferClearOperation_All, 0x000000ff);
  HandleFakeInput();

  DrawAtmosphere(pProgramState, pRenderContext, pRenderContext->pTexture[1 - pRenderContext->activeRenderTarget], pRenderContext->pDepthTexture[1 - pRenderContext->activeRenderTarget]);

  vcRender_FXAAPass(pProgramState, pRenderContext);

  //if (selectionBufferActive)
   // vcRender_ApplySelectionBuffer(pProgramState, pRenderContext);

  if (pProgramState->settings.presentation.mouseAnchor != vcAS_None && (pProgramState->pickingSuccess || pProgramState->isUsingAnchorPoint))
  {
    udDouble4x4 mvp = pProgramState->camera.matrices.viewProjection * udDouble4x4::translation(pProgramState->isUsingAnchorPoint ? pProgramState->worldAnchorPoint : pProgramState->worldMousePosCartesian);
    vcGLState_SetFaceMode(vcGLSFM_Solid, vcGLSCM_Back);

    // Render highlighting any occlusion
    vcGLState_SetBlendMode(vcGLSBM_Additive);
    vcGLState_SetDepthStencilMode(vcGLSDM_Greater, false);
    vcCompass_Render(pRenderContext->pCompass, pProgramState->settings.presentation.mouseAnchor, mvp, udDouble4::create(0.0, 0.15, 1.0, 0.5));

    // Render non-occluded
    vcGLState_SetBlendMode(vcGLSBM_Interpolative);
    vcGLState_SetDepthStencilMode(vcGLSDM_Less, true);
    vcCompass_Render(pRenderContext->pCompass, pProgramState->settings.presentation.mouseAnchor, mvp);

    vcGLState_ResetState();
  }

  if (pProgramState->settings.presentation.showCompass)
  {
    udDouble4x4 cameraRotation = udDouble4x4::rotationYPR(pProgramState->camera.matrices.camera.extractYPR());
    vcGLState_SetFaceMode(vcGLSFM_Solid, vcGLSCM_Back);
    vcGLState_SetDepthStencilMode(vcGLSDM_Always, false);

    if (!pProgramState->gis.isProjected)
    {
      vcCompass_Render(pRenderContext->pCompass, vcAS_Compass, udDouble4x4::perspectiveZO(vcLens30mm, aspect, 0.01, 2.0) * udDouble4x4::translation(vcLens30mm * 0.45 * aspect, 1.0, -vcLens30mm * 0.45) * udDouble4x4::scaleUniform(vcLens30mm / 20.0) * udInverse(cameraRotation));
    }
    else
    {
      udDouble3 currentLatLong = udGeoZone_CartesianToLatLong(pProgramState->gis.zone, pProgramState->camera.position);
      currentLatLong.x = udClamp(currentLatLong.x, -90.0, 89.0);
      udDouble3 norther = udGeoZone_LatLongToCartesian(pProgramState->gis.zone, udDouble3::create(currentLatLong.x + 1.0, currentLatLong.y, currentLatLong.z));
      udDouble4x4 north = udDouble4x4::lookAt(pProgramState->camera.position, norther);
      vcCompass_Render(pRenderContext->pCompass, vcAS_Compass, udDouble4x4::perspectiveZO(vcLens30mm, aspect, 0.01, 2.0) * udDouble4x4::translation(vcLens30mm * 0.45 * aspect, 1.0, -vcLens30mm * 0.45) * udDouble4x4::scaleUniform(vcLens30mm / 20.0) * udDouble4x4::rotationYPR(north.extractYPR()) * udInverse(cameraRotation));
    }

    vcGLState_ResetState();
  }

  vcShader_Bind(nullptr);
  vcGLState_SetViewport(0, 0, pRenderContext->sceneResolution.x, pRenderContext->sceneResolution.y);
}

void vcRender_SceneImGui(vcState *pProgramState, vcRenderContext *pRenderContext, const vcRenderData &renderData)
{
  // Labels
  for (size_t i = 0; i < renderData.labels.length; ++i)
    vcLabelRenderer_Render(renderData.labels[i], pProgramState->camera.matrices.viewProjection, pRenderContext->sceneResolution);
}

udResult vcRender_RecreateUDView(vcState *pProgramState, vcRenderContext *pRenderContext)
{
  udResult result = udR_Success;

  UD_ERROR_NULL(pRenderContext, udR_InvalidParameter_);

  if (pRenderContext->udRenderContext.pRenderView && vdkRenderView_Destroy(&pRenderContext->udRenderContext.pRenderView) != vE_Success)
    UD_ERROR_SET(udR_InternalError);

  if (vdkRenderView_Create(pProgramState->pVDKContext, &pRenderContext->udRenderContext.pRenderView, pRenderContext->udRenderContext.pRenderer, pRenderContext->sceneResolution.x, pRenderContext->sceneResolution.y) != vE_Success)
    UD_ERROR_SET(udR_InternalError);

  if (vdkRenderView_SetTargets(pRenderContext->udRenderContext.pRenderView, pRenderContext->udRenderContext.pColorBuffer, 0, pRenderContext->udRenderContext.pDepthBuffer) != vE_Success)
    UD_ERROR_SET(udR_InternalError);

  if (vdkRenderView_SetMatrix(pRenderContext->udRenderContext.pRenderView, vdkRVM_Projection, pProgramState->camera.matrices.projectionUD.a) != vE_Success)
    UD_ERROR_SET(udR_InternalError);

epilogue:
  return result;
}

udResult vcRender_RenderUD(vcState *pProgramState, vcRenderContext *pRenderContext, vdkRenderView *pRenderView, vcCamera *pCamera, vcRenderData &renderData, bool doPick)
{
  if (pRenderContext == nullptr)
    return udR_InvalidParameter_;

  vdkRenderInstance *pModels = nullptr;
  vcUDRSData *pVoxelShaderData = nullptr;

  int numVisibleModels = 0;

  vdkRenderView_SetMatrix(pRenderView, vdkRVM_Projection, pCamera->matrices.projectionUD.a);
  vdkRenderView_SetMatrix(pRenderView, vdkRVM_View, pCamera->matrices.view.a);

  if (renderData.models.length > 0)
  {
    pModels = udAllocType(vdkRenderInstance, renderData.models.length, udAF_None);
    pVoxelShaderData = udAllocType(vcUDRSData, renderData.models.length, udAF_None);
  }

  double maxDistSqr = pProgramState->settings.camera.farPlane * pProgramState->settings.camera.farPlane;
  pProgramState->pSceneWatermark = nullptr;

  for (size_t i = 0; i < renderData.models.length; ++i)
  {
    if (renderData.models[i]->m_visible && renderData.models[i]->m_loadStatus == vcSLS_Loaded)
    {
      // Copy to the contiguous array
      pModels[numVisibleModels].pPointCloud = renderData.models[i]->m_pPointCloud;
      memcpy(&pModels[numVisibleModels].matrix, renderData.models[i]->m_sceneMatrix.a, sizeof(pModels[numVisibleModels].matrix));
      pModels[numVisibleModels].modelFlags = vdkRMF_None;

      pModels[numVisibleModels].pVoxelShader = vcVoxelShader_Black;
      pModels[numVisibleModels].pVoxelUserData = &pVoxelShaderData[numVisibleModels];

      pVoxelShaderData[numVisibleModels].pModel = renderData.models[i];
      pVoxelShaderData[numVisibleModels].pProgramData = pProgramState;

      switch (pProgramState->settings.visualization.mode)
      {
      case vcVM_Intensity:
        if (vdkAttributeSet_GetOffsetOfStandardAttribute(&renderData.models[i]->m_pointCloudHeader.attributes, vdkSA_Intensity, &pVoxelShaderData[numVisibleModels].attributeOffset) == vE_Success)
        {
          pModels[numVisibleModels].pVoxelShader = vcVoxelShader_Intensity;

          pVoxelShaderData[numVisibleModels].data.intensity.maxIntensity = (uint16_t)pProgramState->settings.visualization.maxIntensity;
          pVoxelShaderData[numVisibleModels].data.intensity.minIntensity = (uint16_t)pProgramState->settings.visualization.minIntensity;
          pVoxelShaderData[numVisibleModels].data.intensity.intensityRange = (float)(pProgramState->settings.visualization.maxIntensity - pProgramState->settings.visualization.minIntensity);
        }

        break;
      case vcVM_Classification:
        if (vdkAttributeSet_GetOffsetOfStandardAttribute(&renderData.models[i]->m_pointCloudHeader.attributes, vdkSA_Classification, &pVoxelShaderData[numVisibleModels].attributeOffset) == vE_Success)
        {
          pModels[numVisibleModels].pVoxelShader = vcVoxelShader_Classification;
        }

        break;
      case vcVM_Displacement:
        if (vdkAttributeSet_GetOffsetOfNamedAttribute(&renderData.models[i]->m_pointCloudHeader.attributes, "udDisplacement", &pVoxelShaderData[numVisibleModels].attributeOffset) == vE_Success)
        {
          pModels[numVisibleModels].pVoxelShader = vcVoxelShader_Displacement;

          pVoxelShaderData[numVisibleModels].data.displacement.minThreshold = pProgramState->settings.visualization.displacement.x;
          pVoxelShaderData[numVisibleModels].data.displacement.maxThreshold = pProgramState->settings.visualization.displacement.y;
        }

        break;
      default: //Includes vcVM_Colour
        if (vdkAttributeSet_GetOffsetOfStandardAttribute(&renderData.models[i]->m_pointCloudHeader.attributes, vdkSA_ARGB, &pVoxelShaderData[numVisibleModels].attributeOffset) == vE_Success)
        {
          pModels[numVisibleModels].pVoxelShader = vcVoxelShader_Colour;
        }
        break;
      }

      ++numVisibleModels;

      if (renderData.models[i]->m_hasWatermark)
      {
        udDouble3 distVector = pCamera->position - renderData.models[i]->GetWorldSpacePivot();

        double cameraDistSqr = udMagSq(distVector);
        if (cameraDistSqr < maxDistSqr)
        {
          maxDistSqr = cameraDistSqr;

          if (renderData.models[i]->m_pWatermark == nullptr) // Load the watermark
          {
            const char *pWatermarkStr = renderData.models[i]->m_metadata.Get("Watermark").AsString();
            if (pWatermarkStr)
            {
              uint8_t *pImage = nullptr;
              size_t imageLen = 0;
              if (udBase64Decode(&pImage, &imageLen, pWatermarkStr) == udR_Success)
              {
                int imageWidth, imageHeight, imageChannels;
                unsigned char *pImageData = stbi_load_from_memory(pImage, (int)imageLen, &imageWidth, &imageHeight, &imageChannels, 4);
                vcTexture_Create(&renderData.models[i]->m_pWatermark, imageWidth, imageHeight, pImageData, vcTextureFormat_RGBA8, vcTFM_Nearest, false);
                free(pImageData);
              }

              udFree(pImage);
            }
          }

          pProgramState->pSceneWatermark = renderData.models[i]->m_pWatermark;
        }
      }
    }
  }

  if (pProgramState->settings.presentation.showDiagnosticInfo && pProgramState->gis.isProjected)
  {
    vcFenceRenderer_ClearPoints(pRenderContext->pDiagnosticFences);

    float z = 0;
    if (pProgramState->settings.maptiles.mapEnabled)
      z = pProgramState->settings.maptiles.mapHeight;

    udInt2 slippyCurrent;
    udDouble3 localCurrent;

    double xRange = pProgramState->gis.zone.latLongBoundMax.x - pProgramState->gis.zone.latLongBoundMin.x;
    double yRange = pProgramState->gis.zone.latLongBoundMax.y - pProgramState->gis.zone.latLongBoundMin.y;

    if (xRange > 0 || yRange > 0)
    {
      std::vector<udDouble3> corners;

      for (int i = 0; i < yRange; ++i)
      {
        vcGIS_LatLongToSlippy(&slippyCurrent, udDouble3::create(pProgramState->gis.zone.latLongBoundMin.x, pProgramState->gis.zone.latLongBoundMin.y + i, 0), 21);
        vcGIS_SlippyToLocal(&pProgramState->gis, &localCurrent, slippyCurrent, 21);
        corners.push_back(udDouble3::create(localCurrent.x, localCurrent.y, z));
      }
      for (int i = 0; i < xRange; ++i)
      {
        vcGIS_LatLongToSlippy(&slippyCurrent, udDouble3::create(pProgramState->gis.zone.latLongBoundMin.x + i, pProgramState->gis.zone.latLongBoundMax.y, 0), 21);
        vcGIS_SlippyToLocal(&pProgramState->gis, &localCurrent, slippyCurrent, 21);
        corners.push_back(udDouble3::create(localCurrent.x, localCurrent.y, z));
      }
      for (int i = 0; i < yRange; ++i)
      {
        vcGIS_LatLongToSlippy(&slippyCurrent, udDouble3::create(pProgramState->gis.zone.latLongBoundMax.x, pProgramState->gis.zone.latLongBoundMax.y - i, 0), 21);
        vcGIS_SlippyToLocal(&pProgramState->gis, &localCurrent, slippyCurrent, 21);
        corners.push_back(udDouble3::create(localCurrent.x, localCurrent.y, z));
      }
      for (int i = 0; i < xRange; ++i)
      {
        vcGIS_LatLongToSlippy(&slippyCurrent, udDouble3::create(pProgramState->gis.zone.latLongBoundMax.x - i, pProgramState->gis.zone.latLongBoundMin.y, 0), 21);
        vcGIS_SlippyToLocal(&pProgramState->gis, &localCurrent, slippyCurrent, 21);
        corners.push_back(udDouble3::create(localCurrent.x, localCurrent.y, z));
      }
      corners.push_back(udDouble3::create(corners[0]));

      vcFenceRenderer_AddPoints(pRenderContext->pDiagnosticFences, &corners[0], corners.size());
    }
  }

  vdkRenderPicking picking = {};
  picking.x = (uint32_t)((float)renderData.mouse.position.x / (float)pRenderContext->originalSceneResolution.x * (float)pRenderContext->sceneResolution.x);
  picking.y = (uint32_t)((float)renderData.mouse.position.y / (float)pRenderContext->originalSceneResolution.y * (float)pRenderContext->sceneResolution.y);

  vdkRenderOptions renderOptions;
  memset(&renderOptions, 0, sizeof(vdkRenderOptions));

  if (doPick)
  {
    pProgramState->udModelPickedIndex = -1;
    renderOptions.pPick = &picking;
  }

  renderOptions.pFilter = renderData.pQueryFilter;
  renderOptions.pointMode = (vdkRenderContextPointMode)pProgramState->settings.presentation.pointMode;

  vdkError result = vdkRenderContext_Render(pRenderContext->udRenderContext.pRenderer, pRenderView, pModels, numVisibleModels, &renderOptions);

  if (result == vE_Success)
  {
    if (doPick && picking.hit)
    {
      // More to be done here
      pProgramState->pickingSuccess = true;
      pProgramState->worldMousePosCartesian = udDouble3::create(picking.pointCenter[0], picking.pointCenter[1], picking.pointCenter[2]);

      uint32_t j = 0;
      for (size_t i = 0; i < renderData.models.length; ++i)
      {
        if (renderData.models[i]->m_visible && renderData.models[i]->m_loadStatus == vcSLS_Loaded)
        {
          if (j == picking.modelIndex)
          {
            pProgramState->udModelPickedIndex = (int)i;
            break;
          }
          ++j;
        }
      }
    }
  }
  else
  {
    //TODO: Clear the buffers
  }

  udFree(pModels);
  udFree(pVoxelShaderData);
  return udR_Success;
}

void vcRender_ClearTiles(vcRenderContext *pRenderContext)
{
  if (pRenderContext == nullptr)
    return;

  vcTileRenderer_ClearTiles(pRenderContext->pTileRenderer);
}

void vcRender_ClearPoints(vcRenderContext *pRenderContext)
{
  if (pRenderContext == nullptr)
    return;

  vcFenceRenderer_ClearPoints(pRenderContext->pDiagnosticFences);
}

vcRenderPickResult vcRender_PolygonPick(vcState *pProgramState, vcRenderContext *pRenderContext, vcRenderData &renderData, bool doSelectRender)
{
  vcRenderPickResult result = {};

  pRenderContext->currentMouseUV = udFloat2::create((float)renderData.mouse.position.x / (float)pRenderContext->originalSceneResolution.x, (float)renderData.mouse.position.y / (float)pRenderContext->originalSceneResolution.y);
  if (pRenderContext->currentMouseUV.x < 0 || pRenderContext->currentMouseUV.x > 1 || pRenderContext->currentMouseUV.y < 0 || pRenderContext->currentMouseUV.y > 1)
    return result;

  pRenderContext->picking.location.x = (uint32_t)(pRenderContext->currentMouseUV.x * pRenderContext->effectResolution.x);
  pRenderContext->picking.location.y = (uint32_t)(pRenderContext->currentMouseUV.y * pRenderContext->effectResolution.y);

  double currentDist = pProgramState->settings.camera.farPlane;
  float pickDepth = 1.0f;

  if (doSelectRender && (renderData.models.length > 0 || renderData.polyModels.length > 0))
  {
    // render pickable geometry with id encoded in colour
    vcGLState_SetBlendMode(vcGLSBM_None);
    vcGLState_SetDepthStencilMode(vcGLSDM_LessOrEqual, true);
    vcGLState_SetFaceMode(vcGLSFM_Solid, vcGLSCM_Back);

    vcFramebuffer_Bind(pRenderContext->picking.pFramebuffer, vcFramebufferClearOperation_All);

    vcGLState_SetViewport(0, 0, pRenderContext->effectResolution.x, pRenderContext->effectResolution.y);
    vcGLState_Scissor(pRenderContext->picking.location.x, pRenderContext->picking.location.y, pRenderContext->picking.location.x + 1, pRenderContext->picking.location.y + 1);

    {
      uint32_t modelId = 1; // note: start at 1, because 0 is 'null'

      // Polygon Models
      for (size_t i = 0; i < renderData.polyModels.length; ++i)
      {
        vcRenderPolyInstance *pInstance = &renderData.polyModels[i];
        udFloat4 idAsColour = vcRender_EncodeIdAsColour((uint32_t)(modelId++));

        vcGLState_SetFaceMode(vcGLSFM_Solid, pInstance->cullFace);

        if (pInstance->renderType == vcRenderPolyInstance::RenderType_Polygon)
          vcPolygonModel_Render(pInstance->pModel, pInstance->worldMat, pProgramState->camera.matrices.viewProjection, vcPMP_ColourOnly, nullptr, &idAsColour);
        else if (pInstance->renderType == vcRenderPolyInstance::RenderType_SceneLayer)
          vcSceneLayerRenderer_Render(pInstance->pSceneLayer, pInstance->worldMat, pProgramState->camera.matrices.viewProjection, pProgramState->camera.position, pRenderContext->sceneResolution, &idAsColour);

      }

      vcGLState_SetFaceMode(vcGLSFM_Solid, vcGLSCM_Back);
    }

    udUInt2 readLocation = { pRenderContext->picking.location.x, pRenderContext->picking.location.y };
    uint8_t colourBytes[4] = {};
    uint8_t depthBytes[4] = {};

#if GRAPHICS_API_OPENGL
    // read upside down
    readLocation.y = pRenderContext->effectResolution.y - pRenderContext->picking.location.y - 1;
#endif

    // Synchronously read back data
    vcTexture_BeginReadPixels(pRenderContext->picking.pTexture, readLocation.x, readLocation.y, 1, 1, colourBytes, pRenderContext->picking.pFramebuffer);
    vcTexture_BeginReadPixels(pRenderContext->picking.pDepth, readLocation.x, readLocation.y, 1, 1, depthBytes, pRenderContext->picking.pFramebuffer);

    vcGLState_SetViewport(0, 0, pRenderContext->sceneResolution.x, pRenderContext->sceneResolution.y);

    // 24 bit unsigned int -> float
#if GRAPHICS_API_OPENGL || GRAPHICS_API_METAL
    pickDepth = uint32_t((depthBytes[3] << 16) | (depthBytes[2] << 8) | (depthBytes[1] << 0)) / ((1 << 24) - 1.0f);
    //uint8_t stencil = depthBytes[0];
#else
    pickDepth = uint32_t((depthBytes[2] << 16) | (depthBytes[1] << 8) | (depthBytes[0] << 0)) / ((1 << 24) - 1.0f);
    //uint8_t stencil = depthBytes[3];
#endif

    // note `-1`, and BGRA format
    int pickedPolygonId = (int)((colourBytes[1] << 0) | (colourBytes[0] << 8)) - 1;
    if (pickedPolygonId != -1)
    {
      result.success = true;
      result.pPolygon = &renderData.polyModels[pickedPolygonId];
    }
  }
  else
  {
    result.success = true;
    pickDepth = pRenderContext->previousFrameDepth;
  }

  if (pickDepth == 0.0)
    pickDepth = 1.0;

  if (result.success)
  {
    // note: upside down (1.0 - uv.y)
    udDouble4 clipPos = udDouble4::create(pRenderContext->currentMouseUV.x * 2.0 - 1.0, (1.0 - pRenderContext->currentMouseUV.y) * 2.0 - 1.0, pickDepth, 1.0);
#if GRAPHICS_API_OPENGL
    clipPos.z = clipPos.z * 2.0 - 1.0;
#endif
    udDouble4 pickPosition = pProgramState->camera.matrices.inverseViewProjection * clipPos;
    pickPosition = pickPosition / pickPosition.w;
    result.position = pickPosition.toVector3();

    currentDist = udMag3(result.position - pProgramState->camera.position);
  }

  if (pProgramState->settings.maptiles.mapEnabled && pProgramState->settings.maptiles.mouseInteracts)// check map tiles
  {
    udPlane<double> mapPlane = udPlane<double>::create({ 0, 0, pProgramState->settings.maptiles.mapHeight }, { 0, 0, 1 });

    double hitDistance = 0.0;
    udDouble3 hitPoint = {};

    if (mapPlane.intersects(pProgramState->camera.worldMouseRay, &hitPoint, &hitDistance))
    {
      if (hitDistance < (currentDist - pProgramState->settings.camera.nearPlane))
      {
        result.success = true;
        result.position = hitPoint;
        result.pPolygon = nullptr;
        result.pModel = nullptr;
      }
    }
  }

  return result;
}
