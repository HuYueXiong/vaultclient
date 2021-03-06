cbuffer u_cameraPlaneParams
{
  float s_CameraNearPlane;
  float s_CameraFarPlane;
  float u_unused1;
  float u_unused2;
};

struct PS_INPUT
{
  float4 pos : SV_POSITION;
  float2 uv : TEXCOORD0;
  float2 edgeSampleUV0 : TEXCOORD1;
  float2 edgeSampleUV1 : TEXCOORD2;
  float2 edgeSampleUV2 : TEXCOORD3;
  float2 edgeSampleUV3 : TEXCOORD4;
};

struct PS_OUTPUT
{
  float4 Color0 : SV_Target;
  float Depth0 : SV_Depth;
};

sampler sampler0;
Texture2D texture0;

sampler sampler1;
Texture2D texture1;

cbuffer u_fragParams : register(b0)
{
  float4 u_screenParams;  // sampleStepSizeX, sampleStepSizeY, (unused), (unused)
  float4x4 u_inverseViewProjection;
  float4x4 u_inverseProjection;

  // outlining
  float4 u_outlineColour;
  float4 u_outlineParams;   // outlineWidth, threshold, (unused), (unused)

  // colour by height
  float4 u_colourizeHeightColourMin;
  float4 u_colourizeHeightColourMax;
  float4 u_colourizeHeightParams; // min world height, max world height, (unused), (unused)

  // colour by depth
  float4 u_colourizeDepthColour;
  float4 u_colourizeDepthParams; // min distance, max distance, (unused), (unused)

  // contours
  float4 u_contourColour;
  float4 u_contourParams; // contour distance, contour band height, contour rainbow repeat rate, contour rainbow factoring
};

float logToLinearDepth(float logDepth)
{
  float a = s_CameraFarPlane / (s_CameraFarPlane - s_CameraNearPlane);
  float b = s_CameraFarPlane * s_CameraNearPlane / (s_CameraNearPlane - s_CameraFarPlane);
  float worldDepth = pow(2.0, logDepth * log2(s_CameraFarPlane + 1.0)) - 1.0;
  return a + b / worldDepth;
}

float getNormalizedPosition(float v, float min, float max)
{
  return clamp((v - min) / (max - min), 0.0, 1.0);
}

// note: an adjusted depth is packed into the returned .w component
// this is to show the edge highlights against the skybox
float4 edgeHighlight(PS_INPUT input, float3 col, float depth, float logDepth, float4 outlineColour, float edgeOutlineThreshold)
{
  float4 eyePosition = mul(u_inverseProjection, float4(input.uv.x * 2.0 - 1.0, (1.0 - input.uv.y) * 2.0 - 1.0, depth, 1.0));
  eyePosition /= eyePosition.w;

  float sampleDepth0 = texture1.Sample(sampler1, input.edgeSampleUV0).x;
  float sampleDepth1 = texture1.Sample(sampler1, input.edgeSampleUV1).x;
  float sampleDepth2 = texture1.Sample(sampler1, input.edgeSampleUV2).x;
  float sampleDepth3 = texture1.Sample(sampler1, input.edgeSampleUV3).x;

  float4 eyePosition0 = mul(u_inverseProjection, float4(input.edgeSampleUV0.x * 2.0 - 1.0, (1.0 - input.edgeSampleUV0.y) * 2.0 - 1.0, logToLinearDepth(sampleDepth0), 1.0));
  float4 eyePosition1 = mul(u_inverseProjection, float4(input.edgeSampleUV1.x * 2.0 - 1.0, (1.0 - input.edgeSampleUV1.y) * 2.0 - 1.0, logToLinearDepth(sampleDepth1), 1.0));
  float4 eyePosition2 = mul(u_inverseProjection, float4(input.edgeSampleUV2.x * 2.0 - 1.0, (1.0 - input.edgeSampleUV2.y) * 2.0 - 1.0, logToLinearDepth(sampleDepth2), 1.0));
  float4 eyePosition3 = mul(u_inverseProjection, float4(input.edgeSampleUV3.x * 2.0 - 1.0, (1.0 - input.edgeSampleUV3.y) * 2.0 - 1.0, logToLinearDepth(sampleDepth3), 1.0));

  eyePosition0 /= eyePosition0.w;
  eyePosition1 /= eyePosition1.w;
  eyePosition2 /= eyePosition2.w;
  eyePosition3 /= eyePosition3.w;

  float3 diff0 = eyePosition.xyz - eyePosition0.xyz;
  float3 diff1 = eyePosition.xyz - eyePosition1.xyz;
  float3 diff2 = eyePosition.xyz - eyePosition2.xyz;
  float3 diff3 = eyePosition.xyz - eyePosition3.xyz;

  float isEdge = 1.0 - step(length(diff0), edgeOutlineThreshold) * step(length(diff1), edgeOutlineThreshold) * step(length(diff2), edgeOutlineThreshold) * step(length(diff3), edgeOutlineThreshold);

  float3 edgeColour = lerp(col.xyz, u_outlineColour.xyz, u_outlineColour.w);
  float edgeLogDepth = min(min(min(sampleDepth0, sampleDepth1), sampleDepth2), sampleDepth3);
  return lerp(float4(col.xyz, logDepth), float4(edgeColour, edgeLogDepth), isEdge);
}

float3 hsv2rgb(float3 c)
{
  float4 K = float4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
  float3 p = abs(frac(c.xxx + K.xyz) * 6.0 - K.www);
  return c.z * lerp(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

float3 contourColour(float3 col, float3 fragWorldPosition)
{
  float contourDistance = u_contourParams.x;
  float contourBandHeight = u_contourParams.y;
  float contourRainboxRepeat = u_contourParams.z;
  float contourRainboxIntensity = u_contourParams.w;

  float3 rainbowColour = hsv2rgb(float3(fragWorldPosition.z * (1.0 / contourRainboxRepeat), 1.0, 1.0));
  float3 baseColour = lerp(col.xyz, rainbowColour, contourRainboxIntensity);

  float isContour = 1.0 - step(contourBandHeight, fmod(abs(fragWorldPosition.z), contourDistance));
  return lerp(baseColour, u_contourColour.xyz, isContour * u_contourColour.w);
}

float3 colourizeByHeight(float3 col, float3 fragWorldPosition)
{
  float2 worldColourMinMax = u_colourizeHeightParams.xy;

  float minMaxColourStrength = getNormalizedPosition(fragWorldPosition.z, worldColourMinMax.x, worldColourMinMax.y);

  float3 minColour = lerp(col.xyz, u_colourizeHeightColourMin.xyz, u_colourizeHeightColourMin.w);
  float3 maxColour = lerp(col.xyz, u_colourizeHeightColourMax.xyz, u_colourizeHeightColourMax.w);
  return lerp(minColour, maxColour, minMaxColourStrength);
}

float3 colourizeByEyeDistance(float3 col, float3 fragEyePos)
{
  float2 depthColourMinMax = u_colourizeDepthParams.xy;

  float depthColourStrength = getNormalizedPosition(length(fragEyePos), depthColourMinMax.x, depthColourMinMax.y);
  return lerp(col.xyz, u_colourizeDepthColour.xyz, depthColourStrength * u_colourizeDepthColour.w);
}

PS_OUTPUT main(PS_INPUT input)
{
  PS_OUTPUT output;

  float4 col = texture0.Sample(sampler0, input.uv);
  float logDepth = texture1.Sample(sampler1, input.uv).x;
  float depth = logToLinearDepth(logDepth);

  // TODO: I'm fairly certain this is actually wrong (world space calculations), and will have precision issues
  float4 fragWorldPosition = mul(u_inverseViewProjection, float4(input.uv.x * 2.0 - 1.0, (1.0 - input.uv.y) * 2.0 - 1.0, depth, 1.0));
  fragWorldPosition /= fragWorldPosition.w;

  float4 fragEyePosition = mul(u_inverseProjection, float4(input.uv.x * 2.0 - 1.0, (1.0 - input.uv.y) * 2.0 - 1.0, depth, 1.0));
  fragEyePosition /= fragEyePosition.w;

  col.xyz = colourizeByHeight(col.xyz, fragWorldPosition.xyz);
  col.xyz = colourizeByEyeDistance(col.xyz, fragEyePosition.xyz);

  col.xyz = contourColour(col.xyz, fragWorldPosition.xyz);

  float edgeOutlineWidth = u_outlineParams.x;
  float edgeOutlineThreshold = u_outlineParams.y;
  float4 outlineColour = u_outlineColour;
  if (edgeOutlineWidth > 0.0 && u_outlineColour.w > 0.0)
  {
    float4 edgeResult = edgeHighlight(input, col.xyz, depth, logDepth, outlineColour, edgeOutlineThreshold);
    col.xyz = edgeResult.xyz;
    logDepth = edgeResult.w; // to preserve outlines, depth written may be adjusted
  }

  output.Color0 = float4(col.xyz, 1.0);// UD always opaque
  output.Depth0 = logDepth;
  return output;
}
