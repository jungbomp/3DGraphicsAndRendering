#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"

/***********************************************/
/* HW1 methods: copy here the methods from HW1 */

GzRender::GzRender(int xRes, int yRes)
{
/* HW1.1 create a framebuffer for MS Windows display:
 -- set display resolution
 -- allocate memory for framebuffer : 3 bytes(b, g, r) x width x height
 -- allocate memory for pixel buffer
 */
	framebuffer = new char[xRes * yRes * 3];
	pixelbuffer = new GzPixel[xRes * yRes];
	xres = xRes;
	yres = yRes;
}

GzRender::~GzRender()
{
/* HW1.2 clean up, free buffer memory */
	if (framebuffer) {
		delete[] framebuffer;
		framebuffer = NULL;
	}

	if (pixelbuffer) {
		delete[] pixelbuffer;
		pixelbuffer = NULL;
	}
}

int GzRender::GzDefault()
{
/* HW1.3 set pixel buffer to some default values - start a new frame */
	memset(framebuffer, '\0', sizeof(char) * xres * yres * 3);
	
	for (int i = 0; i < xres * yres; i++) {
		pixelbuffer[i].red = 0x0800;
		pixelbuffer[i].green = 0x0700;
		pixelbuffer[i].blue = 0x0600;
		pixelbuffer[i].alpha = 0;
		pixelbuffer[i].z = INT_MAX;
	}

	return GZ_SUCCESS;
}


int GzRender::GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
/* HW1.4 write pixel values into the buffer */
	if (i < 0 || j < 0 || xres <= i || yres <= j) {
		return GZ_FAILURE;
	}

	GzPixel& pixel = pixelbuffer[ARRAY(i, j)];
	pixel.red = r;
	pixel.green = g;
	pixel.blue = b;
	pixel.alpha = a;
	pixel.z = z;

	return GZ_SUCCESS;
}


int GzRender::GzGet(int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
/* HW1.5 retrieve a pixel information from the pixel buffer */
	if (i < 0 || j < 0 || xres <= i || yres <= j) {
		return GZ_FAILURE;
	}

	GzPixel pixel = pixelbuffer[ARRAY(i, j)];
	*r = pixel.red;
	*g = pixel.green;
	*b = pixel.blue;
	*a = pixel.alpha;
	*z = pixel.z;

	return GZ_SUCCESS;
}


int GzRender::GzFlushDisplay2File(FILE* outfile)
{
/* HW1.6 write image to ppm file -- "P6 %d %d 255\r" */
	char* buf = new char[xres * yres * 3];
	memset(buf, '\0', sizeof(char) * xres * yres * 3);

	for (int i = 0, j = 0; i < xres * yres; i++) {
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


/***********************************************/
/* HW2 methods: implement from here */

int GzRender::GzPutAttribute(int numAttributes, GzToken	*nameList, GzPointer *valueList)
{
	/* HW 2.1
	-- Set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
	-- In later homeworks set shaders, interpolaters, texture maps, and lights
	*/
	for (int i = 0; i < numAttributes; i++) {
		if (GZ_RGB_COLOR == nameList[i]) {
			flatcolor[RED] = ((float*)valueList[i])[RED];
			flatcolor[GREEN] = ((float*)valueList[i])[GREEN];
			flatcolor[BLUE] = ((float*)valueList[i])[BLUE];
		}
	}

	return GZ_SUCCESS;
}

int GzRender::GzPutTriangle(int	numParts, GzToken *nameList, GzPointer *valueList)
/* numParts - how many names and values */
{
	/* HW 2.2
	-- Pass in a triangle description with tokens and values corresponding to
		  GZ_NULL_TOKEN:		do nothing - no values
		  GZ_POSITION:		3 vert positions in model space
	-- Invoke the rastrizer/scanline framework
	-- Return error code
	*/
	GzCoord vertices[3];
	for (int i = 0; i < numParts; i++) {
		if (GZ_POSITION == nameList[i]) {
			memcpy(vertices, valueList[i], sizeof(GzCoord) * 3);
			
			// sort vertices
			if (vertices[0][Y] > vertices[2][Y] || (vertices[0][Y] == vertices[2][Y] && vertices[0][X] > vertices[2][Y])) {
				GzCoord coord;
				memcpy(coord, vertices[0], sizeof(GzCoord));
				memcpy(vertices[0], vertices[2], sizeof(GzCoord));
				memcpy(vertices[2], coord, sizeof(GzCoord));
			}

			if (vertices[1][Y] > vertices[2][Y] || (vertices[1][Y] == vertices[2][Y] && vertices[1][X] > vertices[2][Y])) {
				GzCoord coord;
				memcpy(coord, vertices[1], sizeof(GzCoord));
				memcpy(vertices[1], vertices[2], sizeof(GzCoord));
				memcpy(vertices[2], coord, sizeof(GzCoord));
			}

			if (vertices[0][Y] > vertices[1][Y] || (vertices[0][Y] == vertices[1][Y] && vertices[0][X] > vertices[1][Y])) {
				GzCoord coord;
				memcpy(coord, vertices[0], sizeof(GzCoord));
				memcpy(vertices[0], vertices[1], sizeof(GzCoord));
				memcpy(vertices[1], coord, sizeof(GzCoord));
			}

			// set check area - lefttop - bottomright
			int left = 0;
			int top = 0;
			int right = 0;
			int bottom = 0;

			// get edges
			const int A = 0; // A = dY
			const int B = 1; // B = -dX
			const int C = 2; // C = dX*Y1-dY*X1
			const int D = 3; // D = dZ
			const int SLOPE_X = 4; // dX/dY (X2-X1) / (Y2-Y1)
			const int SLOPE_Z = 5; // dZ/dY (Z2-Z1) / (Y2-Y1)
			int head[3] = { 1, 2, 2 };
			int tail[3] = { 0, 0, 1 };
			float edges[3][6] = { 0., };
			for (int i = 0; i < sizeof(edges) / sizeof(edges[i]); i++) {
				edges[i][A] = vertices[head[i]][Y] - vertices[tail[i]][Y];
				edges[i][B] = -(vertices[head[i]][X] - vertices[tail[i]][X]);
				edges[i][C] = edges[i][B] * vertices[tail[i]][Y] - edges[i][A] * vertices[tail[i]][X];
				edges[i][D] = vertices[head[i]][Z] - vertices[tail[i]][Z];
				edges[i][SLOPE_X] = (vertices[head[i]][X] - vertices[tail[i]][X]) / edges[i][A];
				edges[i][SLOPE_Z] = (vertices[head[i]][Z] - vertices[tail[i]][Z]) / (vertices[head[i]][Y] - vertices[tail[i]][Y]);
			}

			int leftStartVert = 0;
			int rightStartVert = 0;
			int leftEdge = 0; // dX/dy < dX/dY
			int rightEdge = 1; // dX/dY > dX/dY

			if (vertices[0][Y] == vertices[1][Y]) {
				rightStartVert = 1;
				leftEdge = 1;
				rightEdge = 2;
			} else {
				if (edges[1][SLOPE_X] < edges[0][SLOPE_X]) {
					leftEdge = 1;
					rightEdge = 0;
				}
			}
			
			for (int i = 0; i < 2; i++) {
				// top/left rule -> Pixeling when the pixel is on the top, left, or in the triangle.
				float topY = ceil(vertices[0][Y]);
				if (vertices[0][Y] != vertices[1][Y] && (vertices[0][Y] == topY)) {
					topY += 1.0;
				}
				
				for (int row = (int)topY; row < (int)ceil(vertices[i+1][Y]); row++) {
					float dY = (float)row - vertices[i][Y];
					float leftX = vertices[leftStartVert][X] + (edges[leftEdge][SLOPE_X] * dY);
					float rightX = vertices[rightStartVert][X] + (edges[rightEdge][SLOPE_X] * dY);
					float leftZ = vertices[leftStartVert][Z] + (edges[leftEdge][SLOPE_Z] * dY);
					float rightZ = vertices[rightStartVert][Z] + (edges[rightEdge][SLOPE_Z] * dY);
					float slopeZ = (rightZ - leftZ) / (rightX - leftX);

					for (int col = (int)ceil(leftX); col < (int)ceil(rightX); col++) {
						float dX = (float)col - leftX;
						GzDepth zVal = (GzDepth)(leftZ + slopeZ * dX);
						GzPixel pixel;

						if (GZ_FAILURE == GzGet(col, row, &pixel.red, &pixel.green, &pixel.blue, &pixel.alpha, &pixel.z)) {
							continue;
						}

						if (zVal < pixel.z) {
							GzPut(col, row, ctoi(flatcolor[RED]), ctoi(flatcolor[GREEN]), ctoi(flatcolor[BLUE]), 1, zVal);
						}
					}
				}

				if (i < 1 && vertices[0][Y] != vertices[1][Y] && vertices[1][Y] != vertices[2][Y]) {
					GzCoord midVertex;
					midVertex[Y] = vertices[1][Y];
					float dY = vertices[1][Y] - vertices[0][Y];
					midVertex[X] = vertices[0][X] + (edges[1][SLOPE_X] * dY);
					midVertex[Z] = vertices[0][Z] + (edges[1][SLOPE_Z] * dY);
					
					if (midVertex[X] < vertices[1][X]) {
						memcpy(vertices[0], midVertex, sizeof(GzCoord));
						leftEdge = 1;
						rightEdge = 2;
					} else {
						memcpy(vertices[0], vertices[1], sizeof(GzCoord));
						memcpy(vertices[1], midVertex, sizeof(GzCoord));

						leftEdge = 2;
						rightEdge = 1;
					}

					rightStartVert = 1;
				} else {
					break;
				}
			}
		}
	}

	return GZ_SUCCESS;
}

