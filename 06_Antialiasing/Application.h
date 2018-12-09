// Application.h: interface for the Application class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_APPLICATION_H__3387B79A_B69F_491D_B782_81D9CAFAAB0F__INCLUDED_)
#define AFX_APPLICATION_H__3387B79A_B69F_491D_B782_81D9CAFAAB0F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#define AAKERNEL_SIZE 6

#include "Gz.h"
#include "rend.h"

class Application  
{
public:
	Application();
	virtual ~Application();
	
public:
	GzRender*  m_pRender[AAKERNEL_SIZE];		// the renderer
	GzInput*   m_pUserInput;
	char* m_pFrameBuffer;	// Frame Buffer
	GzPixel* m_pixelbuffer;
	int   m_nWidth;			// width of Frame Buffer
	int   m_nHeight;		// height of Frame Buffer

	virtual int Render() = 0; // Pass user input data and call renderer
	virtual int NormalizePixelBuffer() = 0;	// Normalize all pixel buffers in renders using AAFilter weight values.
	virtual int GzFlushDisplay2File(FILE* outfile) = 0; /* write pixels to ppm file based on display class -- "P6 %d %d 255\r" */
	virtual int GzFlushDisplay2FrameBuffer() = 0;
};

#endif // !defined(AFX_APPLICATION_H__3387B79A_B69F_491D_B782_81D9CAFAAB0F__INCLUDED_)
