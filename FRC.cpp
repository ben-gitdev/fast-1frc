//
#pragma warning(disable:4996) //disable one warning for fopen 

#include "pch.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Windows.h>
//#include "mkl_dfti.h"
#include "mkl.h"
#include<time.h>
#include <omp.h>
#include"FRC.h"


//function definations
long double logFactorial(unsigned __int16 n);
void binomSplit(unsigned __int16* imgIn, float* imgOutA, float* imgOutB, __int32 roiLeft, __int32 roiTop, __int32 fftImgHeight, __int32 fftImgWidth);
void unpackCCE(float* imgIn, MKL_Complex8* imgOut, int height, int width);
void createFourierRing(unsigned __int16* fourierRing, int fourierRingWidth, int fourierRingHeight);
void frcFromFFTs(float* frc, float* imgOutA, float* imgOutB, unsigned __int16* fourierRing, int fftImgHeight, int fftImgWidth, int fourierRingWidth, int frcLength);
void writeFloat(wchar_t* fileName, float* data, int dataLength);
void generateRandomNumbers();//unsigned __int16* randomNumbers, int nRow, int nCol, int maxElectronCount
int generateBinomTable();
unsigned __int16 maxIntensity(unsigned __int16* imgIn);
//void frc1(unsigned __int16* imgIn, float* frc);
//int initialize();
//void uninitialize();


// global variables
DFTI_DESCRIPTOR_HANDLE fft = NULL;
unsigned __int16 nRow;
unsigned __int16 nCol;
unsigned __int16 maxElectronCount;
//unsigned __int16 bitDepth;

int fftImgWidth;
int fftImgHeight;
float* imgOutA;
float* imgOutB;
unsigned __int16* randomNumbers;
unsigned __int16* binomNumbers;

int fourierRingWidth;
int fourierRingHeight;
unsigned __int16* fourierRing;


int roiNColFull;
int roiNRowFull;


DFTI_DESCRIPTOR_HANDLE fftSub = NULL;
__int32 subImgSize;
int fftSubImgWidth;
int fourierRingSubWidth;
int fourierRingSubHeight;
int fftSubImgHeight;
unsigned __int16* fourierRingSub;
float* imgSubOutA;
float* imgSubOutB;

__int16 intensityOffset = 0;
__int32 I2eFactor = 4;



/*
 ---------------------------------------------
|                                             |
|             HELPER FUNCTIONS                |
|                                             |
 ---------------------------------------------
*/


//==================binomial split====================
void binomSplit(unsigned __int16* imgIn, float* imgOutA, float* imgOutB, __int32 roiLeft, __int32 roiTop, __int32 fftImgHeightIn, __int32 fftImgWidthIn) {//int nRow, int nCol,, unsigned __int16* randomNumbers, unsigned __int16* binomNumbers
//#pragma omp parallel for shared(imgIn, imgOutA, imgOutB, roiTop, roiLeft, fftImgHeightIn, fftImgWidthIn)//nCol,nRow,,randomNumbers,binomNumbers
	for (int i = 0; i < fftImgHeightIn; i++) {
		int iPlusRoiTopxnCol = (i + roiTop) * nCol;
		for (int j = 0; j < fftImgWidthIn; j++) {
			int index = iPlusRoiTopxnCol + j + roiLeft;
			unsigned __int16 pixelI = imgIn[index];
			unsigned __int16 eCount;
			if (pixelI > intensityOffset) {
				if ( I2eFactor == 1 ) { eCount = pixelI - intensityOffset; }
				else if (I2eFactor == 2) {
					eCount = (pixelI - intensityOffset) / 2;
				}
				else if (I2eFactor == 4) {
					eCount = (pixelI - intensityOffset) / 4;
				}
				else {
					eCount = (pixelI - intensityOffset) / I2eFactor;
				}	
			}
			else
				eCount = 0;
			unsigned __int16 randomNumber = randomNumbers[index];
			unsigned __int16 pixelA = binomNumbers[eCount * maxElectronCount + randomNumber];
			int imgOutInd = i * fftImgWidthIn + j;
			imgOutA[imgOutInd] = float(pixelA);
			imgOutB[imgOutInd] = float(pixelI - pixelA);
		}
	}
}



//=================unpack CCE complex FFT result===============================
//CCE complex FFT result saves FFT result in real float/double array
//it only stores FFT results in qudrants 2 and 3, 
//data in qudrants 1 and 4 are rotationally symatric to data in 2 and 3
//data arrangement:
// real00, imag00, real01, imag01, real02, imag02......
// real10, imag10, real11, imag11, real12, imag12......
// real20, imag20, real21, imag21, real22, imag22......
// ......
//realxy is real part and imagxy is imaginary part of the result
void unpackCCE(float* imgIn, MKL_Complex8* imgOut, int height, int width) {
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			MKL_Complex8 val;
			if (j < width / 2 + 1)
			{
				val.real = imgIn[i * (width / 2 + 1) * 2 + 2 * j];
				val.imag = imgIn[i * (width / 2 + 1) * 2 + 2 * j + 1];
			}
			else // unpack CCE format result
			{
				if (i == 0)
				{
					val.real = imgIn[(width - j) * 2];
					val.imag = imgIn[(width - j) * 2 + 1] * (-1);
				}
				else
				{
					val.real = imgIn[(height - i) * (width / 2 + 1) * 2 + (width - j) * 2];
					val.imag = imgIn[(height - i) * (width / 2 + 1) * 2 + (width - j) * 2 + 1] * (-1);
				}
			}
			imgOut[i * width + j] = val;
		}
	}
}


//===================== create fourier rings===================================
//due to symetry of fft result, only need to create rings for quadrants 2 and 3.
void createFourierRing(unsigned __int16* fourierRing, int fourierRingWidth, int fourierRingHeight) {//
	float widthHeightFactor = float(fourierRingWidth - 1) / float(fourierRingHeight / 2);
	for (int i = 0; i < fourierRingHeight; i++) {
		for (int j = 0; j < fourierRingWidth; j++) {
			if (i < fourierRingHeight / 2 + 1)
				fourierRing[i * fourierRingWidth + j] = (unsigned __int16) (round(sqrt(i * i + float(j * j) / widthHeightFactor / widthHeightFactor)));//
			if (i >= fourierRingHeight / 2 + 1)
				fourierRing[i * fourierRingWidth + j] = (unsigned __int16) (round(sqrt((i - fourierRingHeight) * (i - fourierRingHeight) + float(j * j) / widthHeightFactor / widthHeightFactor)));//unsigned __int16
		}
	}
}


//=======================calculate FRC from FFT results===================================
void frcFromFFTs(float* frc, float* fftA, float* fftB, unsigned __int16* fourierRing, int fftImgHeight, int fftImgWidth, int fourierRingWidth, int frcLength) {
	//int fftImgWidth = (roiNCol / 2 + 1) * 2;
	//buffers for the three sums in fourier ring equations, use 1024 FRC points
	// max frc length in code set to 2048 (more than enough for most cases), if user requsts longer frc, data after 2048 will be 0;
	float sum1[2048] = { 0 };
	float sum2[2048] = { 0 };
	float sum3[2048] = { 0 };

	// go through fourier result, place results in corresponding rings and calculate the sums
//#pragma omp parallel for //shared(fftImgWidth,fftImgHeight, fourierRing,fourierRingWidth,fftA,fftB, frcLength,sum1,sum2,sum3)
//#pragma omp  simd 
	for (int i = 1; i < fftImgHeight; i++) {
		int ixfourierRingWidth = i * fourierRingWidth;
		int ixfftImgWidth = i * fftImgWidth;
		for (int j = 1; j < fftImgWidth / 2 - 1; j++) {
			// i and j start from 1 to leave the center cross of fft out to reduce fft edge effect when not using a window, this reduces noise of the frc curve
			int frcIndex = fourierRing[ixfourierRingWidth + j];
			int limit;
			if (frcLength < 2048) limit = frcLength;
			else limit = 2048;
			if (frcIndex < limit) {
				// due to symmetry of fft result, frc equation could be simplified to the following
				int fftIndex = ixfftImgWidth + j * 2;
				int fftIndexPlus1 = fftIndex + 1;
				sum1[frcIndex] += (fftA[fftIndex] * fftB[fftIndex] + fftA[fftIndexPlus1] * fftB[fftIndexPlus1]);//realA*realB + imagA*imagB
				sum2[frcIndex] += (fftA[fftIndex] * fftA[fftIndex] + fftA[fftIndexPlus1] * fftA[fftIndexPlus1]);//realA*realA + imagA*imagA
				sum3[frcIndex] += (fftB[fftIndex] * fftB[fftIndex] + fftB[fftIndexPlus1] * fftB[fftIndexPlus1]);//realB*realB + imagB*imagB
			}
		}
	}
	//calculate final frc curve
	if (frcLength <= 2048) {
		for (int i = 1; i < frcLength; i++) {
				frc[i] = sum1[i] / sqrt(sum2[i] * sum3[i]);
			}
	}
	else {
		for (int i = 1; i < 2048; i++) {
			frc[i] = sum1[i] / sqrt(sum2[i] * sum3[i]);
		}
		for (int i = 2048; i < frcLength; i++) {
			frc[i] = 0;
		}
	}
	frc[0] = 1;
}


//=================write float data to disk in raw format=====================
void writeFloat(wchar_t* fileName, float* data, int dataLength) {
	HANDLE hFile = CreateFile(
		fileName,     // Filename
		GENERIC_WRITE,          // Desired access
		FILE_SHARE_READ,        // Share mode
		NULL,                   // Security attributes
		CREATE_ALWAYS,             // Creates a new file, only if it doesn't already exist
		FILE_ATTRIBUTE_NORMAL,  // Flags and attributes
		NULL);                  // Template file handle
	DWORD bytesWritten;
	WriteFile(
		hFile,            // Handle to the file
		data,  // Buffer to write
		dataLength * sizeof(float),   // Buffer size
		&bytesWritten,    // Bytes written
		nullptr);         // Overlapped
	CloseHandle(hFile);// Close the handle once we don't need it.
}


//=======================generate random numbers=============================
// random number array size same as input image
void generateRandomNumbers() {//unsigned __int16* randomNumbers, int nRow, int nCol, int maxElectronCount
	srand(12345);
	for (int i = 0; i < nRow; i++) {
		for (int j = 0; j < nCol; j++) {
			randomNumbers[i * nCol + j] = (unsigned __int16) (double(rand()) / double(RAND_MAX) * (maxElectronCount - 1));//
		}
	}
}


//=====================generate table of binomial numbers===========================
//binomial numbers of different N and k combination
int generateBinomTable() {
	unsigned __int16* binomNOccur;
	binomNOccur = (unsigned __int16*)malloc(maxElectronCount * sizeof(unsigned __int16));
	if (binomNOccur == NULL)
		return 1;
	long double* logFactorials;
	logFactorials = (long double*)malloc(maxElectronCount * sizeof(long double));
	if (logFactorials == NULL)
		return 1;

	for (int i = 0; i < maxElectronCount; i++) {
		logFactorials[i] = logFactorial(i);
	}

	// calculate the binom numbers table for the binomial distribution: each number should appear x times for different binomial paramters
	// each row are all the numbers that forms binomial distribution of N = rowIndex, k=colIndex, p = 0.5
	for (int i = 0; i < maxElectronCount; i++) {
		double binomPdf;
		long double logBinomPdf;
		for (int j = 0; j < maxElectronCount; j++) {
			if (j > i) {
				binomPdf = 0;
			}
			else if (j == i) {
				binomPdf = pow(0.5, j);
			}
			else { //factorials in the binom pdf function can be too huge to hold in standard c data types, calculate log(binom pdf) instead
				logBinomPdf = logFactorials[i] - logFactorials[j] - logFactorials[i - j] + i * logl(0.5);
				binomPdf = pow(exp(1), logBinomPdf);
			}
			binomNOccur[j] = round(binomPdf * maxElectronCount);
		}
		// create actual binomial number table from nOccur
		int index = 0;
		for (int j = 0; j < maxElectronCount; j++) {
			if (binomNOccur[j] > 0) {
				for (int k = index; k < index + binomNOccur[j]; k++) {
					if (k < maxElectronCount) binomNumbers[i * maxElectronCount + k] = j;
				}
				index += binomNOccur[j];
			}
		}
		if (index < maxElectronCount) { //due to roundings in calculating binom probability and rounding, we might have spaces not filled, then fill it with expectation of the distrubution.
			for (int k = index; k < maxElectronCount; k++) {
				binomNumbers[i * maxElectronCount + k] = i / 2;
			}
		}
	}
	free(binomNOccur);
	free(logFactorials);
	return 0;
}


//===================calculate log of factorial of a number=====================
long double logFactorial(unsigned __int16 n) {
	long double result = 0;
	for (long double i = 1; i <= n; i++) {
		result += logl(i);
	}
	return result;
}

// find max intensity of an image, images size is global nRow x nCol;
unsigned __int16 maxIntensity(unsigned __int16* imgIn) {
	unsigned __int16 maxI = 0;
	unsigned __int16 pixI = 0;
	for (int i = 10; i < nRow - 10; i += 2) {
		for (int j = 10; j < nCol - 10; j += 2) {
			pixI = imgIn[i * nCol + j];
			if (pixI > maxI)
				maxI = pixI;
		}
	}
	return maxI;
}

/*
 ---------------------------------------------
|                                             |
|                DLL FUNCTIONS                |
|                                             |
 ---------------------------------------------
*/
//=====================CALCULATE FRC FROM SINGLE INPUT IMAGE=======================

void frc1(unsigned __int16* imgIn, float* frc, __int32 roiLeftIn, __int32 roiTopIn) {	// input buffer frc length is always 1024 for full image, and subImgSize/2 for sub image
	MKL_LONG status;
	if (roiLeftIn < 0 || roiTopIn < 0) { // do 1frc on full image
		int roiLeft = 0;
		int roiTop = 0;
		//int roiNRow = roiNRowFull;
		//int roiNCol = roiNColFull;
		int frcLength = roiNRowFull / 2;
		binomSplit(imgIn, imgOutA, imgOutB, roiLeft, roiTop, fftImgHeight, fftImgWidth);
		status = DftiComputeForward(fft, imgOutA);
		status = DftiComputeForward(fft, imgOutB);
		frcFromFFTs(frc, imgOutA, imgOutB, fourierRing, fftImgHeight, fftImgWidth, fourierRingWidth, frcLength);
	}
	else { // do 1frc on roi subimage
		int roiLeft = roiLeftIn;
		int roiTop = roiTopIn;
		if (roiLeft + subImgSize > nCol) {
			roiLeft = nCol - subImgSize;
		}
		if (roiTop + subImgSize > nRow) {
			roiTop = nRow - subImgSize;
		}
		//int roiNRow = subImgSize;
		//int roiNCol = subImgSize;
		int frcLength = subImgSize / 2;
		binomSplit(imgIn, imgSubOutA, imgSubOutB, roiLeft, roiTop, fftSubImgHeight, fftSubImgWidth);
		status = DftiComputeForward(fftSub, imgSubOutA);
		status = DftiComputeForward(fftSub, imgSubOutB);
		frcFromFFTs(frc, imgSubOutA, imgSubOutB, fourierRingSub, fftSubImgHeight, fftSubImgWidth, fourierRingSubWidth, frcLength);
	}
}




int initialize(__int32 nRowFull, __int32 nColFull, __int32 subImgSizeIn, __int16 intensityOffsetIn, __int32 I2eFactorIn, unsigned __int16 maxI) {
	intensityOffset = intensityOffsetIn;
	I2eFactor = I2eFactorIn;
	nRow = nRowFull;
	nCol = nColFull;
	//bitDepth = 12; // according to lightning manual, intensity to electron conversion factor is 0.27, ~=0.25, so 12 bit img becomes 10 bit electron. Accordingly, in function binomsplit, the input image is divided by 4
	maxElectronCount = (maxI - intensityOffset) / I2eFactor;
	
	roiNColFull = nCol;
	roiNRowFull = nRow;
	fftImgWidth = (roiNColFull / 2 + 1) * 2;
	fftImgHeight = roiNRowFull;
	fourierRingWidth = roiNColFull / 2 + 1;
	fourierRingHeight = fftImgHeight;


	// initialize values for sub image frc calculation.
	subImgSize = subImgSizeIn; //sub image size recommended to be multiples of 2, e.g. 256 for best speed
	fftSubImgWidth = (subImgSize / 2 + 1) * 2;
	fftSubImgHeight = subImgSize;
	fourierRingSubWidth = subImgSize / 2 + 1;
	fourierRingSubHeight = fftSubImgHeight;



	randomNumbers = (unsigned __int16*)malloc(nRow * nCol * sizeof(unsigned __int16));
	if (randomNumbers == NULL)
		return 1;
	generateRandomNumbers();//randomNumbers, nRow, nCol, maxElectronCount

	binomNumbers = (unsigned __int16*)malloc(maxElectronCount * maxElectronCount * sizeof(unsigned __int16));
	if (binomNumbers == NULL)
		return 1;
	int status1 = generateBinomTable();
	if (status1 > 0)
		return status1;

	fourierRing = (unsigned __int16*)malloc(fourierRingHeight * fourierRingWidth * sizeof(unsigned __int16));
	if (fourierRing == NULL)
		return 1;
	createFourierRing(fourierRing, fourierRingWidth, fourierRingHeight);

	fourierRingSub = (unsigned __int16*)malloc(fourierRingSubHeight * fourierRingSubWidth * sizeof(unsigned __int16));
	if (fourierRingSub == NULL)
		return 1;
	createFourierRing(fourierRingSub, fourierRingSubWidth, fourierRingSubHeight);


	imgOutA = (float*)mkl_malloc(fftImgHeight * fftImgWidth * sizeof(float), 64);
	if (imgOutA == NULL)
		return 1;
	imgOutB = (float*)mkl_malloc(fftImgHeight * fftImgWidth * sizeof(float), 64);
	if (imgOutB == NULL)
		return 1;

	imgSubOutA = (float*)mkl_malloc(fftSubImgHeight * fftSubImgWidth * sizeof(float), 64);
	if (imgSubOutA == NULL)
		return 1;
	imgSubOutB = (float*)mkl_malloc(fftSubImgHeight * fftSubImgWidth * sizeof(float), 64);
	if (imgSubOutB == NULL)
		return 1;

	//setup fft
	MKL_LONG status;
	MKL_LONG dim_sizes[2] = { fftImgHeight, roiNColFull };
	MKL_LONG inputStride[3] = { 0, (dim_sizes[1] / 2 + 1) * 2, 1 };
	MKL_LONG outputStride[3] = { 0, (dim_sizes[1] / 2 + 1), 1 };
	status = DftiCreateDescriptor(&fft, DFTI_SINGLE, DFTI_REAL, 2, dim_sizes);
	status = DftiSetValue(fft, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
	status = DftiSetValue(fft, DFTI_INPUT_STRIDES, inputStride);
	status = DftiSetValue(fft, DFTI_OUTPUT_STRIDES, outputStride);
	status = DftiCommitDescriptor(fft);
	if (status > 0)
		return status;

	MKL_LONG dim_sizesSub[2] = { fftSubImgHeight, subImgSize };
	MKL_LONG inputStrideSub[3] = { 0, (dim_sizesSub[1] / 2 + 1) * 2, 1 };
	MKL_LONG outputStrideSub[3] = { 0, (dim_sizesSub[1] / 2 + 1), 1 };
	status = DftiCreateDescriptor(&fftSub, DFTI_SINGLE, DFTI_REAL, 2, dim_sizesSub);
	status = DftiSetValue(fftSub, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
	status = DftiSetValue(fftSub, DFTI_INPUT_STRIDES, inputStrideSub);
	status = DftiSetValue(fftSub, DFTI_OUTPUT_STRIDES, outputStrideSub);
	status = DftiCommitDescriptor(fftSub);

	return status;
}


void uninitialize() {
	free(randomNumbers);
	free(binomNumbers);
	free(fourierRing);
	free(fourierRingSub);
	mkl_free(imgOutA);
	mkl_free(imgOutB);
	DftiFreeDescriptor(&fft);
	mkl_free(imgSubOutA);
	mkl_free(imgSubOutB);
	DftiFreeDescriptor(&fftSub);
}

void downscale2(unsigned __int16* imgIn, unsigned __int32 nRow, unsigned __int32 nCol, unsigned __int16* imgOut) {
	int nRowHalf = nRow / 2;
	int nColHalf = nCol / 2;
//#pragma omp parallel for shared(imgIn, nRow, nCol, imgOut, nRowHalf, nColHalf)
//#pragma omp parallel for collapse(2)
	for (int i = 0; i < nRowHalf; i++) {
		int ix2xnCol = i * 2 * nCol;
		int ixnColHalf = i * nColHalf;
		for (int j = 0; j < nColHalf; j++) {
			int sum = 0;
			int ind1 = ix2xnCol + j * 2;
			sum += imgIn[ind1];
			sum += imgIn[ind1 + 1];
			int ind3 = ind1 + nCol;
			sum += imgIn[ind3];
			sum += imgIn[ind3 + 1];
			imgOut[ixnColHalf + j] = (unsigned __int16) ((sum / 4 ) + ((sum >>1)&1));
		}
	}
}

void downscale16(unsigned __int16* imgIn, unsigned __int32 nRow, unsigned __int32 nCol, unsigned __int16* imgOut) {
	int nRowDown = nRow / 16;
	int nColDown = nCol / 16;
#pragma omp parallel for collapse(2) num_threads(4)
	for (int i = 0; i < nRowDown; i++) {
		int ix16 = i * 16;
		for (int j = 0; j < nColDown; j++) {
			int jx16 = j * 16;
			unsigned __int32 sum = 0;
			for (int k = 0; k < 16; k++) {
				int ix16PluskxnCol = (ix16 + k)*nCol;
				for (int l = 0; l < 16; l++)
					sum += imgIn[ix16PluskxnCol + jx16 + l];
			}
			imgOut[i * nColDown + j] = (unsigned __int16)(sum / 256) + ((sum >> 7) & 1);
		}
	}
}

void downscale8(unsigned __int16* imgIn, unsigned __int32 nRow, unsigned __int32 nCol, unsigned __int16* imgOut) {
	int nRowDown = nRow / 8;
	int nColDown = nCol / 8;
#pragma omp parallel for collapse(2) num_threads(4)
	for (int i = 0; i < nRowDown; i++) {
		int ix8 = i * 8;
		for (int j = 0; j < nColDown; j++) {
			int jx8 = j * 8;
			unsigned __int32 sum = 0;
			for (int k = 0; k < 8; k++) {
				int ix8PluskxnCol = (ix8 + k) * nCol;
				for (int l = 0; l < 8; l++)
					sum += imgIn[ix8PluskxnCol + jx8 + l];
			}
			imgOut[i * nColDown + j] = (unsigned __int16)(sum / 64) + ((sum >> 5) & 1);
		}
	}
}