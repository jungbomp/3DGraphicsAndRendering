#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
/*   CS580 HW   */
#include    "stdafx.h"  
#include	"Gz.h"


GzRender::GzRender(int xRes, int yRes)
{
/* HW1.1 create a framebuffer for MS Windows display:
 -- set display resolution
 -- allocate memory for framebuffer : 3 bytes(b, g, r) x width x height
 -- allocate memory for pixel buffer
 */
	framebuffer = new char[xRes*yRes*3];
	pixelbuffer = new GzPixel[xRes*yRes];
	xres = xRes;
	yres = yRes;
}

GzRender::~GzRender()
{
/* HW1.2 clean up, free buffer memory */
	delete[] framebuffer;
	delete[] pixelbuffer;

}

int GzRender::GzDefault()
{
/* HW1.3 set pixel buffer to some default values - start a new frame */
	memset(framebuffer, '\0', sizeof(char)*xres*yres * 3);
	memset(pixelbuffer, '\0', sizeof(GzPixel)*xres*yres);

	for (int i = 0; i < xres*yres; i++) {
		pixelbuffer[i].red = 0x0800;
		pixelbuffer[i].green = 0x0700;
		pixelbuffer[i].blue = 0x0600;
	}

	return GZ_SUCCESS;
}


int GzRender::GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
/* HW1.4 write pixel values into the buffer */
	if (i < 0 || j < 0 || xres <= i || yres <= j) {
		return GZ_FAILURE;
	}

	if (i == 245 && j == 12) {
		int a = 0;
	}

	GzPixel& pixel = pixelbuffer[i + (j * xres)];
	pixel.red = r;
	pixel.green = g;
	pixel.blue = b;
	pixel.alpha = a;
	pixel.z = z;

	/*framebuffer[3 * (j*xres + i) + 0] = (b >> 4);
	framebuffer[3 * (j*xres + i) + 1] = (g >> 4);
	framebuffer[3 * (j*xres + i) + 2] = (r >> 4);*/

	return GZ_SUCCESS;
}


int GzRender::GzGet(int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
/* HW1.5 retrieve a pixel information from the pixel buffer */
	if (i < 0 || j < 0 || xres <= i || yres <= j) {
		return GZ_FAILURE;
	}

	GzPixel pixel = pixelbuffer[i + (j * xres)];
	*r = pixel.red;
	*g = pixel.green;
	*b = pixel.blue;
	*a = pixel.alpha;
	*z = pixel.z;

	return GZ_SUCCESS;
}


int GzRender::GzFlushDisplay2File(FILE* outfile)
{
/* HW1.6 write image to ppm file -- "P6 %d %d 255\n" */
	char* buf = new char[xres * yres * 3];
	memset(buf, '\0', sizeof(char) * xres * yres * 3);

	for (int i = 0, j = 0; i < xres*yres; i++) {
		buf[j++] = (0x0FFF < pixelbuffer[i].red ? 0xFF : (pixelbuffer[i].red >> 4));
		buf[j++] = (0x0FFF < pixelbuffer[i].green ? 0xFF : (pixelbuffer[i].green >> 4));
		buf[j++] = (0x0FFF < pixelbuffer[i].blue ? 0xFF : (pixelbuffer[i].blue >> 4));
	}
	
	fprintf(outfile, "P6 %d %d %d\n", xres, yres, 255);
	fwrite(buf, sizeof(char), xres*yres * 3, outfile);
	delete[] buf;
	
	return GZ_SUCCESS;
}

int GzRender::GzFlushDisplay2FrameBuffer()
{
/* HW1.7 write pixels to framebuffer: 
	- put the pixels into the frame buffer
	- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red 
	- NOT red, green, and blue !!!
*/
	for (int i = 0, j = 0; i < xres*yres; i++) {
		framebuffer[j++] = (0x0FFF < pixelbuffer[i].blue ? 0xFF : (pixelbuffer[i].blue >> 4));
		framebuffer[j++] = (0x0FFF < pixelbuffer[i].green ? 0xFF : (pixelbuffer[i].green >> 4));
		framebuffer[j++] = (0x0FFF < pixelbuffer[i].red ? 0xFF : (pixelbuffer[i].red >> 4));
	}

	return GZ_SUCCESS;
}