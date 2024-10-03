#pragma once
#ifdef FRC_EXPORTS
#define FRC_API __declspec(dllexport)
#else
#define FRC_API __declspec(dllimport)
#endif

extern "C" FRC_API void uninitialize();
extern "C" FRC_API int initialize(__int32 nRowFull, __int32 nColFull, __int32 subImgSizeIn, __int16 intensityOffsetIn, __int32 I2eFactorIn, unsigned __int16 maxI);
extern "C" FRC_API void frc1(unsigned __int16* imgIn, float* frc, __int32 roiLeft, __int32 roiTop);
extern "C" FRC_API void downscale2(unsigned __int16* imgIn, unsigned __int32 nRow, unsigned __int32 nCol, unsigned __int16* imgOut);
extern "C" FRC_API void downscale16(unsigned __int16* imgIn, unsigned __int32 nRow, unsigned __int32 nCol, unsigned __int16* imgOut);
extern "C" FRC_API void downscale8(unsigned __int16* imgIn, unsigned __int32 nRow, unsigned __int32 nCol, unsigned __int16* imgOut);