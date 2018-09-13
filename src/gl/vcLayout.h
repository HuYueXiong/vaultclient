#ifndef vcLayout_h__
#define vcLayout_h__

#include "udPlatform/udMath.h"

enum vcVertexLayoutTypes
{
  vcVLT_Position2, //Vec2
  vcVLT_Position3, //Vec3
  vcVLT_TextureCoords2, //Vec2
  vcVLT_ColourBGRA, //uint32_t
  vcVLT_Normal3, //Vec3

  vcVLT_TotalTypes
};

struct vcSimpleVertex
{
  udFloat3 position;
  udFloat2 uv;
};
const vcVertexLayoutTypes vcSimpleVertexLayout[] = { vcVLT_Position3, vcVLT_TextureCoords2 };

uint32_t vcLayout_GetSize(const vcVertexLayoutTypes *pLayout, int numTypes);

#endif//vcLayout_h__
