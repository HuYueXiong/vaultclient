#ifndef vcAtmosphereRenderer_h__
#define vcAtmosphereRenderer_h__

#include "udResult.h"

struct vcAtmosphereRenderer;

udResult vcAtmosphereRenderer_Create(vcAtmosphereRenderer **ppAtmosphereRenderer);
udResult vcAtmosphereRenderer_Destroy(vcAtmosphereRenderer **ppAtmosphereRenderer);

bool vcAtmosphereRenderer_Render(vcAtmosphereRenderer *pAtmosphereRenderer);

#endif//vcAtmosphereRenderer_h__
