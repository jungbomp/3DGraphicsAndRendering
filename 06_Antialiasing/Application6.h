// Application4.h: interface for the Application4 class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_APPLICATION6_H__43A7FA9C_6CD6_4A79_9567_2354BFEFAFFB__INCLUDED_)
#define AFX_APPLICATION6_H__43A7FA9C_6CD6_4A79_9567_2354BFEFAFFB__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Application.h"

class Application6 : public Application  
{
public:
	Application6();
	virtual ~Application6();
	
	int	Initialize();
	virtual int Render(); 
	virtual int NormalizePixelBuffer(); // Normalize all pixel buffers in renders using AAFilter weight values.
	virtual int GzFlushDisplay2File(FILE* outfile); /* write pixels to ppm file based on display class -- "P6 %d %d 255\r" */
	virtual int GzFlushDisplay2FrameBuffer(); /* write pixels to framebuffer */
	int Clean();
};

#endif // !defined(AFX_APPLICATION6_H__43A7FA9C_6CD6_4A79_9567_2354BFEFAFFB__INCLUDED_)
