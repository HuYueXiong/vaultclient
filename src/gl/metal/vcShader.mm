#import "vcMetal.h"
#import "gl/vcShader.h"
#import <Metal/Metal.h>

#import "udPlatformUtil.h"

// Takes shader function names instead of shader description string
bool vcShader_CreateFromText(vcShader **ppShader, const char *pVertexShader, const char *pFragmentShader, const vcVertexLayoutTypes * pVertLayout, uint32_t totalTypes)
{
  if (ppShader == nullptr || pVertexShader == nullptr || pFragmentShader == nullptr)
    return false;

  NSError *err = nil;
  
  vcShader *pShader = udAllocType(vcShader, 1, udAF_Zero);

  MTLVertexDescriptor *vertexDesc = [[MTLVertexDescriptor alloc] init];
  
  ptrdiff_t accumulatedOffset = 0;
  for (uint32_t i = 0; i < totalTypes; ++i)
  {
    vertexDesc.attributes[i].bufferIndex = 0;
    vertexDesc.attributes[i].offset = accumulatedOffset;
    
    switch (pVertLayout[i])
    {
      case vcVLT_Position2:
        vertexDesc.attributes[i].format = MTLVertexFormatFloat2;
        accumulatedOffset += 2 * sizeof(float);
        break;
      case vcVLT_Position3:
        vertexDesc.attributes[i].format = MTLVertexFormatFloat3;
        accumulatedOffset += 3 * sizeof(float);
        break;
      case vcVLT_TextureCoords2:
        vertexDesc.attributes[i].format = MTLVertexFormatFloat2;
        accumulatedOffset += 2 * sizeof(float);
        break;
      case vcVLT_RibbonInfo4:
        vertexDesc.attributes[i].format = MTLVertexFormatFloat4;
        accumulatedOffset += 4 * sizeof(float);
        break;
      case vcVLT_ColourBGRA:
        vertexDesc.attributes[i].format = MTLVertexFormatUInt;
        accumulatedOffset += 1 * sizeof(uint32_t);
        break;
      case vcVLT_Normal3:
        vertexDesc.attributes[i].format = MTLVertexFormatFloat3;
        accumulatedOffset += 3 * sizeof(float);
        break;
      case vcVLT_TotalTypes:
        break;
    }
    
    vertexDesc.layouts[0].stride = accumulatedOffset;
    vertexDesc.layouts[0].stepFunction = MTLVertexStepFunctionPerVertex;
    vertexDesc.layouts[0].stepRate = 1;
  }
  
  id<MTLFunction> vFunc = [_library newFunctionWithName:[NSString stringWithUTF8String:pVertexShader]];
  id<MTLFunction> fFunc = [_library newFunctionWithName:[NSString stringWithUTF8String:pFragmentShader]];
  
  MTLRenderPipelineDescriptor *pDesc = [[MTLRenderPipelineDescriptor alloc] init];
  pDesc.vertexDescriptor = vertexDesc;
  pDesc.vertexFunction = vFunc;
  pDesc.fragmentFunction = fFunc;
  
  pDesc.colorAttachments[0].pixelFormat = MTLPixelFormatBGRA8Unorm;
#if UDPLATFORM_IOS || UDPLATFORM_IOS_SIMULATOR
  pDesc.depthAttachmentPixelFormat = MTLPixelFormatDepth32Float;
  pDesc.stencilAttachmentPixelFormat = MTLPixelFormatDepth32Float;
#elif UDPLATFORM_OSX
  pDesc.depthAttachmentPixelFormat = MTLPixelFormatDepth24Unorm_Stencil8;
  pDesc.stencilAttachmentPixelFormat = MTLPixelFormatDepth24Unorm_Stencil8;
#else
# error "Unknown platform!"
#endif
    
  id<MTLRenderPipelineState> state = [_device newRenderPipelineStateWithDescriptor:pDesc error:&err];
#ifdef METAL_DEBUG
  if (!state)
  {
    NSLog(@"Build pipeline state failed: %@", err);
    return false;
  }
#endif
  pShader->ID = (uint32_t)_viewCon.renderer.pipelines.count;
  
  [_viewCon.renderer.pipelines addObject:state];
  [_viewCon.renderer.pipeDescs addObject:pDesc];
  
  // Bind here ensures depth/stencil state is constructed before first use
  [_viewCon.renderer bindPipeline:pShader];
  
  *ppShader = pShader;
  pShader = nullptr;

  return (*ppShader != nullptr);
}

void vcShader_DestroyShader(vcShader **ppShader)
{
  if (ppShader == nullptr || *ppShader == nullptr)
    return;

  udFree(*ppShader);
}

bool vcShader_Bind(vcShader *pShader)
{
  if (pShader != nullptr)
  {
    [_viewCon.renderer bindPipeline:pShader];
  }
  return true;
}

bool vcShader_BindTexture(vcShader *pShader, vcTexture *pTexture, uint16_t samplerIndex, vcShaderSampler *pSampler/* = nullptr*/)
{
  udUnused(pShader);
  if (pTexture == nullptr)
    return false;

  [_viewCon.renderer bindTexture:pTexture index:samplerIndex];

  if (pSampler)
  {
    [_viewCon.renderer bindSampler:pSampler index:samplerIndex];
  }
  return true;
}

bool vcShader_GetConstantBuffer(vcShaderConstantBuffer **ppBuffer, vcShader *pShader, const char *pBufferName, const size_t bufferSize)
{
  if (ppBuffer == nullptr || pShader == nullptr || bufferSize == 0)
    return false;
  
  *ppBuffer = nullptr;
  
  for (int i = 0; i < pShader->numBufferObjects; ++i)
  {
    if (udStrEquali(pShader->bufferObjects[i].name, pBufferName) && (bufferSize == pShader->bufferObjects[i].expectedSize))
    {
      *ppBuffer = &pShader->bufferObjects[i];
      return true;
    }
  }
  
  vcShaderConstantBuffer *temp = udAllocType(vcShaderConstantBuffer, 1, udAF_Zero);
  temp->expectedSize = bufferSize;
  temp->ID = (uint32_t)_viewCon.renderer.constantBuffers.count;
  udStrcpy(temp->name, 32, pBufferName);
  
  [_viewCon.renderer.constantBuffers addObject:[_device newBufferWithLength:bufferSize options:MTLStorageModeShared]];
  
  pShader->bufferObjects[pShader->numBufferObjects] = *temp;
  ++pShader->numBufferObjects;
  
  *ppBuffer = temp;
  temp = nullptr;
  
  return true;
}

bool vcShader_BindConstantBuffer(vcShader *pShader, vcShaderConstantBuffer *pBuffer, const void *pData, const size_t bufferSize)
{
  if (pShader == nullptr || pBuffer == nullptr || pData == nullptr || bufferSize == 0)
    return false;

  int found = -1;
  for (int i = 0; i < pShader->numBufferObjects; ++i)
    if (udStrEquali(pShader->bufferObjects[i].name, pBuffer->name))
      found = i;
  
  if (found < 0)
  {
    pShader->bufferObjects[pShader->numBufferObjects] = *pBuffer;
    ++pShader->numBufferObjects;
  }
  
  if (pBuffer->expectedSize >= bufferSize)
  {
    // !!!!!!!!!!!
    // Differences in packing of structs in our code vs metal shader source makes this not work out, may help to move structs to seperate header
    //memcpy(_viewCon.renderer.constantBuffers[pBuffer->ID].contents, pData, bufferSize);
    //[_viewCon.renderer.constantBuffers[pBuffer->ID] didModifyRange:NSMakeRange(0, bufferSize)];
    
    [_viewCon.renderer.constantBuffers replaceObjectAtIndex:pBuffer->ID withObject:[_device newBufferWithBytes:pData length:bufferSize options:MTLStorageModeShared]];
  }
  else
  {
    [_viewCon.renderer.constantBuffers replaceObjectAtIndex:pBuffer->ID withObject:[_device newBufferWithBytes:pData length:bufferSize options:MTLStorageModeShared]];
  }
  
  return true;
}

bool vcShader_ReleaseConstantBuffer(vcShader *pShader, vcShaderConstantBuffer *pBuffer)
{
  if (pShader == nullptr || pBuffer == nullptr)
    return false;
  
  // TODO
    
  return true;
}

bool vcShader_GetSamplerIndex(vcShaderSampler **ppSampler, vcShader *pShader, const char *pSamplerName)
{
  if (pShader == nullptr)
    return false;
  
  for (int i = 0; i < pShader->numBufferObjects; ++i)
  {
    if (udStrEquali(pShader->samplerIndexes[i].name, pSamplerName))
    {
      *ppSampler = &pShader->samplerIndexes[i];
      return true;
    }
  }

  return false;
}