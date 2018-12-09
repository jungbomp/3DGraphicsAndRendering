/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"

#define PI (float) 3.14159265358979323846

int GzRender::GzRotXMat(float degree, GzMatrix mat)
{
/* HW 3.1
// Create rotate matrix : rotate along x axis
// Pass back the matrix using mat value
*/
	if (!mat) return GZ_FAILURE;

	float r = degree * PI / 180.0;
	float sin = sinf(r);
	float cos = cosf(r);

	GzMatrix rotXMat = { 1,		0,		0,		0,
						 0,		cos,	-sin,	0,
						 0,		sin,	cos,	0,
						 0,		0,		0,		1 };
	
	memcpy(mat, rotXMat, sizeof(rotXMat));

	return GZ_SUCCESS;
}

int GzRender::GzRotYMat(float degree, GzMatrix mat)
{
/* HW 3.2
// Create rotate matrix : rotate along y axis
// Pass back the matrix using mat value
*/
	if (!mat) return GZ_FAILURE;

	float r = degree * PI / 180.0;
	float sin = sinf(r);
	float cos = cosf(r);

	GzMatrix rotYMat = { cos,	0,		sin,	0,
						 0,		1,		0,		0,
						 -sin,	0,		cos,	0,
						 0,		0,		0,		1 };

	memcpy(mat, rotYMat, sizeof(rotYMat));

	return GZ_SUCCESS;
}

int GzRender::GzRotZMat(float degree, GzMatrix mat)
{
/* HW 3.3
// Create rotate matrix : rotate along z axis
// Pass back the matrix using mat value
*/
	if (!mat) return GZ_FAILURE;

	float r = degree * PI / 180.0;
	float sin = sinf(r);
	float cos = cosf(r);

	GzMatrix rotZMat = { cos,	-sin,	0,		0,
						 sin,	cos,	0,		0,
						 0,		0,		1,		0,
						 0,		0,		0,		1 };

	memcpy(mat, rotZMat, sizeof(rotZMat));

	return GZ_SUCCESS;
}

int GzRender::GzTrxMat(GzCoord translate, GzMatrix mat)
{
/* HW 3.4
// Create translation matrix
// Pass back the matrix using mat value
*/
	if (!mat) return GZ_FAILURE;

	GzMatrix trxMat = { 1, 0, 0, translate[X],
						0, 1, 0, translate[Y],
						0, 0, 1, translate[Z],
						0, 0, 0, 1 };

	memcpy(mat, trxMat, sizeof(trxMat));
	
	return GZ_SUCCESS;
}


int GzRender::GzScaleMat(GzCoord scale, GzMatrix mat)
{
/* HW 3.5
// Create scaling matrix
// Pass back the matrix using mat value
*/
	if (!mat) return GZ_FAILURE;

	GzMatrix scaleMat = {	scale[X],	0,			0,			0,
							0,			scale[Y],	0,			0,
							0,			0,			scale[Z],	0,
							0,			0,			0,			1 };

	memcpy(mat, scaleMat, sizeof(scaleMat));

	return GZ_SUCCESS;
}


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

/* HW 3.6
- setup Xsp and anything only done once 
- init default camera 
*/ 
	matlevel = -1;

	m_camera.position[X] = DEFAULT_IM_X;
	m_camera.position[Y] = DEFAULT_IM_Y;
	m_camera.position[Z] = DEFAULT_IM_Z;

	m_camera.lookat[X] = 0.0;
	m_camera.lookat[Y] = 0.0;
	m_camera.lookat[Z] = 0.0;

	m_camera.worldup[X] = 0.0;
	m_camera.worldup[Y] = 1.0;
	m_camera.worldup[Z] = 0.0;

	m_camera.FOV = DEFAULT_FOV;

	for (int i = 0; i < sizeof(GzMatrix) / sizeof(Xsp[i]); i++) {
		for (int j = 0; j < sizeof(Xsp[i]) / sizeof(float); j++) {
			Xsp[i][j] = 0.0;
			m_camera.Xiw[i][j] = 0.0;
			m_camera.Xpi[i][j] = 0.0;
		}
	}
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
		pixelbuffer[i].alpha = 1;
		pixelbuffer[i].z = INT_MAX;
	}

	return GZ_SUCCESS;
}

int GzRender::GzBeginRender()
{
	/* HW 3.7
	- setup for start of each frame - init frame buffer color,alpha,z
	- compute Xiw and projection xform Xpi from camera definition
	- init Ximage - put Xsp at base of stack, push on Xpi and Xiw
	- now stack contains Xsw and app can push model Xforms when needed
	*/
	float oneOverD = tanf(m_camera.FOV * PI / 180.0 / 2.0);

	Xsp[0][0] = (float)xres / 2.0;
	Xsp[0][3] = (float)xres / 2.0;
	Xsp[1][1] = -((float)yres / 2.0);
	Xsp[1][3] = (float)yres / 2.0;
	Xsp[2][2] = INT_MAX * oneOverD;
	Xsp[3][3] = 1.0;

	m_camera.Xpi[0][0] = 1.0;
	m_camera.Xpi[1][1] = 1.0;
	m_camera.Xpi[2][2] = 1.0;
	m_camera.Xpi[3][3] = 1.0;
	m_camera.Xpi[3][2] = oneOverD;

	// compute Xiw
	// Compute Z axis
	GzCoord zAxis;
	zAxis[X] = m_camera.lookat[X] - m_camera.position[X];
	zAxis[Y] = m_camera.lookat[Y] - m_camera.position[Y];
	zAxis[Z] = m_camera.lookat[Z] - m_camera.position[Z];

	float length = sqrtf(zAxis[X] * zAxis[X] + zAxis[Y] * zAxis[Y] + zAxis[Z] * zAxis[Z]);
	if (0.0 == length) return GZ_FAILURE;

	zAxis[X] = zAxis[X] / length;
	zAxis[Y] = zAxis[Y] / length;
	zAxis[Z] = zAxis[Z] / length;

	// Compute Y axis
	GzCoord yAxis;
	float upZ = m_camera.worldup[X] * zAxis[X] + m_camera.worldup[Y] * zAxis[Y] + m_camera.worldup[Z] * zAxis[Z];
	yAxis[X] = m_camera.worldup[X] - upZ * zAxis[X];
	yAxis[Y] = m_camera.worldup[Y] - upZ * zAxis[Y];
	yAxis[Z] = m_camera.worldup[Z] - upZ * zAxis[Z];

	length = sqrtf(yAxis[X] * yAxis[X] + yAxis[Y] * yAxis[Y] + yAxis[Z] * yAxis[Z]);
	if (0.0 == length) return GZ_FAILURE;	

	yAxis[X] = yAxis[X] / length;
	yAxis[Y] = yAxis[Y] / length;
	yAxis[Z] = yAxis[Z] / length;

	// Compute X axis
	GzCoord xAxis;
	xAxis[X] = yAxis[Y] * zAxis[Z] - yAxis[Z] * zAxis[Y];
	xAxis[Y] = yAxis[Z] * zAxis[X] - yAxis[X] * zAxis[Z];
	xAxis[Z] = yAxis[X] * zAxis[Y] - yAxis[Y] * zAxis[X];

	// Xiw
	GzCoord* pCoords[3] = { &xAxis, &yAxis, &zAxis };
	for (int i = 0; i < sizeof(pCoords) / sizeof(GzCoord*); i++) {
		m_camera.Xiw[i][0] = (*pCoords[i])[X];
		m_camera.Xiw[i][1] = (*pCoords[i])[Y];
		m_camera.Xiw[i][2] = (*pCoords[i])[Z];
		m_camera.Xiw[i][3] = -((*pCoords[i])[X] * m_camera.position[X] + (*pCoords[i])[Y] * m_camera.position[Y] + (*pCoords[i])[Z] * m_camera.position[Z]);
	}

	m_camera.Xiw[3][3] = 1.0;

	// init Ximage - put Xsp at base of stack
	matlevel = 0;
	memcpy(Ximage[matlevel], Xsp, sizeof(Xsp));
	
	// push on Xpi
	GzPushMatrix(m_camera.Xpi);

	// push Xiw
	GzPushMatrix(m_camera.Xiw);

	return GZ_SUCCESS;
}

int GzRender::GzPutCamera(GzCamera camera)
{
/* HW 3.8 
/*- overwrite renderer camera structure with new camera definition
*/
	memcpy(m_camera.position, camera.position, sizeof(camera.position));
	memcpy(m_camera.lookat, camera.lookat, sizeof(camera.lookat));
	memcpy(m_camera.worldup, camera.worldup, sizeof(camera.worldup));

	m_camera.FOV = camera.FOV;

	return GZ_SUCCESS;	
}

int GzRender::GzPushMatrix(GzMatrix	matrix)
{
/* HW 3.9 
- push a matrix onto the Ximage stack
- check for stack overflow
*/
	if (!matrix) {
		return GZ_FAILURE;
	}

	for (int i = 0; i < 4; i++)	{
		for (int j = 0; j < 4; j++) {
			Ximage[matlevel + 1][i][j] = Ximage[matlevel][i][0] * matrix[0][j]
				+ Ximage[matlevel][i][1] * matrix[1][j]
				+ Ximage[matlevel][i][2] * matrix[2][j]
				+ Ximage[matlevel][i][3] * matrix[3][j];
		}
	}

	matlevel++;
	
	return GZ_SUCCESS;
}

int GzRender::GzPopMatrix()
{
/* HW 3.10
- pop a matrix off the Ximage stack
- check for stack underflow
*/
	if (matlevel < 0) {
		return GZ_FAILURE;
	}

	matlevel--;

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

int GzRender::GzPutTriangle(int numParts, GzToken *nameList, GzPointer *valueList)
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

			for (int j = 0; j < 3; j++) {
				GzCoord homoVert;
				float W;

				// Model to screen
				memcpy(homoVert, vertices[j], sizeof(GzCoord));
				W = Ximage[matlevel][3][0] * homoVert[X] + Ximage[matlevel][3][1] * homoVert[Y] + Ximage[matlevel][3][2] * homoVert[Z] + Ximage[matlevel][3][3] * 1.0;
				if (W == 0.0) return GZ_FAILURE;
				
				vertices[j][X] = (Ximage[matlevel][0][0] * homoVert[X] + Ximage[matlevel][0][1] * homoVert[Y] + Ximage[matlevel][0][2] * homoVert[Z] + Ximage[matlevel][0][3] * 1.0) / W;
				vertices[j][Y] = (Ximage[matlevel][1][0] * homoVert[X] + Ximage[matlevel][1][1] * homoVert[Y] + Ximage[matlevel][1][2] * homoVert[Z] + Ximage[matlevel][1][3] * 1.0) / W;
				vertices[j][Z] = (Ximage[matlevel][2][0] * homoVert[X] + Ximage[matlevel][2][1] * homoVert[Y] + Ximage[matlevel][2][2] * homoVert[Z] + Ximage[matlevel][2][3] * 1.0) / W;

				if (vertices[j][Z] <= 0.0) return GZ_SUCCESS;
			}

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
			}
			else {
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

				for (int row = (int)topY; row < (int)ceil(vertices[i + 1][Y]); row++) {
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
					}
					else {
						memcpy(vertices[0], vertices[1], sizeof(GzCoord));
						memcpy(vertices[1], midVertex, sizeof(GzCoord));

						leftEdge = 2;
						rightEdge = 1;
					}

					rightStartVert = 1;
				}
				else {
					break;
				}
			}
		}
	}

	return GZ_SUCCESS;
}

