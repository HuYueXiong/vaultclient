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
  float3 normal : COLOR0;
  float4 colour : COLOR1;
  float3 sunDirection : COLOR2;
  float4 fragClipPosition : COLOR3;
  float2 fLogDepth : TEXCOORD0;
};

struct PS_OUTPUT
{
  float4 Color0 : SV_Target;
  float Depth0 : SV_Depth;
};

PS_OUTPUT main(PS_INPUT input)
{
  PS_OUTPUT output;

  float3 fakeEyeVector = normalize(input.fragClipPosition.xyz / input.fragClipPosition.w);
  float3 worldNormal = input.normal * float3(2.0, 2.0, 2.0) - float3(1.0, 1.0, 1.0);
  float ndotl = 0.5 + 0.5 * -dot(input.sunDirection, worldNormal);
  float edotr = max(0.0, -dot(-fakeEyeVector, worldNormal));
  edotr = pow(edotr, 60.0);
  float3 sheenColour = float3(1.0, 1.0, 0.9);
  output.Color0 = float4(input.colour.a * (ndotl * input.colour.xyz + edotr * sheenColour), 1.0);

  float halfFcoef = 1.0 / log2(s_CameraFarPlane + 1.0);
  output.Depth0 = log2(input.fLogDepth.x) * halfFcoef;

  return output;
}
