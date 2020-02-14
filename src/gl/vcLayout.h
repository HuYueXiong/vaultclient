#ifndef vcLayout_h__
#define vcLayout_h__

#include "udMath.h"

enum vcVertexLayoutTypes
{
  vcVLT_Position2, //Vec2
  vcVLT_Position3, //Vec3
  vcVLT_Position4, //Vec4
  vcVLT_TextureCoords2, //Vec2
  vcVLT_ColourBGRA, //uint32_t
  vcVLT_Normal3, //Vec3
  vcVLT_RibbonInfo4, // Vec4
  vcVLT_QuadCorner, // Vec2

  vcVLT_Unsupported,

  vcVLT_TotalTypes
};

// Vertex Formats
//vertex structure naming = "vc" N*(layoutTypeAcronym) "Vertex"
//layoutTypeAcronym = acronym for the `vcVertexLayoutTypes` type data definition
struct vcP3Vertex
{
  udFloat3 position;
};

struct vcUV2Vertex
{
  udFloat2 uv;
};

struct vcP3N3Vertex
{
  udFloat3 position;
  udFloat3 normal;
};

struct vcP3UV2Vertex
{
  udFloat3 position;
  udFloat2 uv;
};

struct vcP3N3UV2Vertex
{
  udFloat3 position;
  udFloat3 normal;
  udFloat2 uv;
};

struct vcP3UV2RI4Vertex
{
  udFloat3 position;
  udFloat2 uv;
  udFloat4 ribbonInfo; // xyz: expand vector; w: pair id (0 or 1)
};

struct vcP4C1Vertex
{
  udFloat4 pos;
  uint32_t color;
};

struct vcP4C1QC2Vertex
{
  udFloat4 pos;
  uint32_t color;
  float corner[2];
};

const vcVertexLayoutTypes vcP3VertexLayout[] = { vcVLT_Position3 }; // tiles
const vcVertexLayoutTypes vcP2VertexLayout[] = { vcVLT_Position2 }; // water
const vcVertexLayoutTypes vcP3N3VertexLayout[] = { vcVLT_Position3, vcVLT_Normal3 }; // compass / anchor
const vcVertexLayoutTypes vcP3UV2VertexLayout[] = { vcVLT_Position3, vcVLT_TextureCoords2 }; // screen quad
const vcVertexLayoutTypes vcP3N3UV2VertexLayout[] = { vcVLT_Position3, vcVLT_Normal3, vcVLT_TextureCoords2 }; // polygon model
const vcVertexLayoutTypes vcP3UV2RI4VertexLayout[] = { vcVLT_Position3, vcVLT_TextureCoords2, vcVLT_RibbonInfo4 }; // ribbon
const vcVertexLayoutTypes vcP4C1VertexLayout[] = { vcVLT_Position4, vcVLT_ColourBGRA }; // gpu renderer (point)
const vcVertexLayoutTypes vcP4C1QC2VertexLayout[] = { vcVLT_Position4, vcVLT_ColourBGRA, vcVLT_QuadCorner }; // gpu renderer (quad)

// ImGui
const vcVertexLayoutTypes vcImGuiVertexLayout[] = { vcVLT_Position2, vcVLT_TextureCoords2, vcVLT_ColourBGRA };

uint32_t vcLayout_GetSize(const vcVertexLayoutTypes layoutType);
uint32_t vcLayout_GetSize(const vcVertexLayoutTypes *pLayout, int numTypes);
uint32_t vcLayout_GetSize(const vcVertexLayoutTypes layoutType);

#endif//vcLayout_h__
