/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"

#define PI (float) 3.14159265358979323846

// Utility functions
bool SwapMemory(GzPointer ptr1, GzPointer ptr2, int size) {
	if (NULL == ptr1 || NULL == ptr2) return false;

	byte buf[1024];
	byte* pBuf = buf;

	if (1024 < size) {
		pBuf = new byte[size];
	}

	memcpy(pBuf, ptr1, size);
	memcpy(ptr1, ptr2, size);
	memcpy(ptr2, pBuf, size);

	if (1024 < size) {
		delete[] pBuf;
	}

	return true;
}

float DotVectors(const GzCoord vec1, const GzCoord vec2)
{
	return vec1[X] * vec2[X] + vec1[Y] * vec2[Y] + vec1[Z] * vec2[Z];
}

inline bool NormalizeVector(GzCoord vec) {
	float length = sqrtf(vec[X] * vec[X] + vec[Y] * vec[Y] + vec[Z] * vec[Z]);
	if (0 == length) {
		return false;
	}

	vec[X] = vec[X] / length;
	vec[Y] = vec[Y] / length;
	vec[Z] = vec[Z] / length;

	return true;
}

bool CrossVector(const GzCoord vec1, const GzCoord vec2, GzCoord ret) {
	if (NULL == ret || NULL == vec1 || NULL == vec2) return false;

	ret[X] = vec1[Y] * vec2[Z] - vec1[Z] * vec2[Y];
	ret[Y] = vec1[Z] * vec2[X] - vec1[X] * vec2[Z];
	ret[Z] = vec1[X] * vec2[Y] - vec1[Y] * vec2[X];

	return true;
}

bool homogeneousMultiplyMatrixToVector(const GzMatrix mat, GzCoord vec) {
	GzCoord homoVec;
	float W;

	// Model to screen
	memcpy(homoVec, vec, sizeof(GzCoord));
	W = mat[Z + 1][X] * homoVec[X] + mat[Z + 1][Y] * homoVec[Y] + mat[Z + 1][Z] * homoVec[Z] + mat[Z + 1][Z + 1] * 1.0;
	if (W == 0.0) return false;

	vec[X] = (mat[X][X] * homoVec[X] + mat[X][Y] * homoVec[Y] + mat[X][Z] * homoVec[Z] + mat[X][Z + 1] * 1.0) / W;
	vec[Y] = (mat[Y][X] * homoVec[X] + mat[Y][Y] * homoVec[Y] + mat[Y][Z] * homoVec[Z] + mat[Y][Z + 1] * 1.0) / W;
	vec[Z] = (mat[Z][X] * homoVec[X] + mat[Z][Y] * homoVec[Y] + mat[Z][Z] * homoVec[Z] + mat[Z][Z + 1] * 1.0) / W;

	return true;
}

// Shading Color = Sum_lights(specular + diffuse + ambient components)
//               = (Ks * sum[Ie(dot(R,E)^s]) + (Kd * sum[Ie(dot(N,L)]) + (Ka * Ia) 
// Color, Ie, Ks, Kd, Ka, Ia - are all RGB color vectors([0,1])
// R = reflected ray direction vector (normalized), R = 2(dot(N,L))N - L
// E = eye ray direction vector (normalized)
// N = surface normal vector (normalized)
// L = light ray direction vector (normalized)
// s = exponential factor of specular
bool ShadingColor(const GzCoord imageNormal,
	const int numLights, const GzLight lights[],
	const GzLight ambientLight,
	const float spec,
	const GzColor Ks, GzColor Kd, GzColor Ka, 
	const GzTextureIndex uv,
	const GzTexture tex_fun,
	const int interp_mode,
	GzColor colorResult)
{
	// Intialize color to black just in case
	colorResult[RED] = colorResult[GREEN] = colorResult[BLUE] = 0;

	// So we must sum over all the lights.
	GzColor specular, diffuse;
	specular[RED] = specular[GREEN] = specular[BLUE] = diffuse[RED] = diffuse[GREEN] = diffuse[BLUE] = 0;

	// Surface normal vector
	GzCoord N;
	memcpy(N, imageNormal, sizeof(GzCoord));

	// In image space, the direction to the eye (e.g. camera) is simply (0, 0, -1), normalized vector is also (0, 0, -1)
	GzCoord E = { 0, 0, -1 };

	GzColor textureColor;

	for (int i = 0; i < numLights; i++) {
		// Compute dot(N,L)
		float NdotL = DotVectors(N, lights[i].direction);

		// Test the illumination and the viewing direction relative to the normal of a surface
		//	- Sign of dot(N,L) and dot(N,E)
		//		- Both positive : compute lighting model
		//		- Both negative : flop normal and compute lighting model on backside of surface
		//		- Each has different sigh : light and eye are on opposite sides of the surface so the light contributes zero - skip it

		// Compute dot(N,E)
		float NdotE = DotVectors(N, E);
		if (NdotL * NdotE < 0) {
			continue;
		}

		if (NdotL < 0 && NdotE < 0) {
			// Flip normal
			N[X] *= -1;
			N[Y] *= -1;
			N[Z] *= -1;

			// Flip NdotL & NdotE
			NdotL *= -1;
			NdotE *= -1;
		}

		// L = light ray direction vector
		const GzCoord& L = lights[i].direction;

		// Compute reflected ray. R = 2(N dot L)N - L   
		GzCoord R;
		R[X] = 2 * NdotL *  N[X] - L[X];
		R[Y] = 2 * NdotL *  N[Y] - L[Y];
		R[Z] = 2 * NdotL *  N[Z] - L[Z];

		// make sure reflected ray is normalized
		NormalizeVector(R);

		// Compute dot(R,E)
		float RdotE = DotVectors(R, E);
		// Dot(R,E) calculations must be clamped to zero to maintain in [0,1] bounded range
		//	- Dot(R,E) may be negative for front or back surface illuminations
		if (RdotE < 0)
			RdotE = 0;

		// Ie = RGB color vectors of light
		const GzColor& Ie = lights[i].color;
		// now add in the Kd and Ks contributions from this light

		specular[RED] += Ie[RED] * pow(RdotE, spec);
		specular[GREEN] += Ie[GREEN] * pow(RdotE, spec);
		specular[BLUE] += Ie[BLUE] * pow(RdotE, spec);

		diffuse[RED] += Ie[RED] * NdotL;
		diffuse[GREEN] += Ie[GREEN] * NdotL;
		diffuse[BLUE] += Ie[BLUE] * NdotL;
	}

	if (tex_fun) {
		if (GZ_COLOR == interp_mode) {
			/*
			* For simplicity, we'll set KT = Kd = Ks = Ka for Gouraud textures
			* KT is the texture color (RGB triple)
			*
			* Gouraud or color interpolation is done as light-intensity interpolation from modified vertex shading calculation
			* Delay the multiplication of cumulative light intensities by Ka, Kd, Ks
			* We'll multiply by pixel-by-pixel texture color at pixel rasterization time
			*/
			// C = (KT) * (sumOverLights[ lightIntensity ( R dot E )^s]  + sumOverLights[ lightIntensity ( N dot L )]) + Ia)
			const GzColor& Ia = ambientLight.color;
			colorResult[RED] = specular[RED] + diffuse[RED] + Ia[RED];
			colorResult[GREEN] = specular[GREEN] + diffuse[GREEN] + Ia[GREEN];
			colorResult[BLUE] = specular[BLUE] + diffuse[BLUE] + Ia[BLUE];

			return true;
		}

		if (GZ_NORMALS == interp_mode) {
			
			tex_fun(uv[U], uv[V], textureColor);

			Ka = textureColor;
			Kd = textureColor;
		}
	}

	// Compute specular component 
	specular[RED] *= Ks[RED];
	specular[GREEN] *= Ks[GREEN];
	specular[BLUE] *= Ks[BLUE];

	// Compute diffuse component 
	diffuse[RED] *= Kd[RED];
	diffuse[GREEN] *= Kd[GREEN];
	diffuse[BLUE] *= Kd[BLUE];

	// Ia =  RGB color vectors of ambient light
	const GzColor& Ia = ambientLight.color;

	// Compute ambient component Ka * Ia
	GzCoord ambient;
	ambient[RED] = Ka[RED] * Ia[RED];
	ambient[GREEN] = Ka[GREEN] * Ia[GREEN];
	ambient[BLUE] = Ka[BLUE] * Ia[BLUE];

	//	Color = Sum_lights(specular + diffuse + ambient components) 
	// add all components together
	colorResult[RED] = specular[RED] + diffuse[RED] + ambient[RED];
	colorResult[GREEN] = specular[GREEN] + diffuse[GREEN] + ambient[GREEN];
	colorResult[BLUE] = specular[BLUE] + diffuse[BLUE] + ambient[BLUE];

	return true;
}

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

	memcpy(mat, rotXMat, sizeof(GzMatrix));

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

	memcpy(mat, rotYMat, sizeof(GzMatrix));

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

	memcpy(mat, rotZMat, sizeof(GzMatrix));

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

	memcpy(mat, trxMat, sizeof(GzMatrix));

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

	memcpy(mat, scaleMat, sizeof(GzMatrix));

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

	// Xsp = [	xs/2,	0,		0,		xs/2
	//			0,		-ys/2,	0,		ys/2
	//			0,		0,		MAXINT,	0
	//			0,		0,		0,		1		]
	Xsp[0][0] = static_cast<float>(xres / 2.0);
	Xsp[0][3] = static_cast<float>(xres / 2.0);
	Xsp[1][1] = -(static_cast<float>(yres / 2.0));
	Xsp[1][3] = static_cast<float>(yres / 2.0);
	Xsp[2][2] = INT_MAX;
	Xsp[3][3] = 1.0;

	numlights = 0;

	tex_fun = 0;
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
	// Compute Xpi from camera difinition
	// 1/d = tan(FOV/2)
	// Xpi = [	1,		0,		0,		0
	//			0,		1,		0,		0
	//			0,		0,		1/d,	0
	//			0,		0,		1/d,	1	]
	float oneOverD = tanf(m_camera.FOV * PI / 180.0 / 2.0);
	m_camera.Xpi[0][0] = 1.0;
	m_camera.Xpi[1][1] = 1.0;
	m_camera.Xpi[2][2] = oneOverD;
	m_camera.Xpi[3][3] = 1.0;
	m_camera.Xpi[3][2] = oneOverD;

	// compute Xiw
	// Compute Z axis (Z = Vec_cl / || Vec_cl ||, c = camera's position, l = look at point, Vec_cl = l - c)
	GzCoord zAxis;
	zAxis[X] = m_camera.lookat[X] - m_camera.position[X];
	zAxis[Y] = m_camera.lookat[Y] - m_camera.position[Y];
	zAxis[Z] = m_camera.lookat[Z] - m_camera.position[Z];

	// || Vec_cl ||
	if (!NormalizeVector(zAxis)) return GZ_FAILURE;

	// Compute Y axis (Y = up' / || up' ||, up' = the up vector in image space, up' = up - dot(up,Z)Z, up = camera's up)
	GzCoord yAxis;
	float upZ = DotVectors(m_camera.worldup, zAxis);
	yAxis[X] = m_camera.worldup[X] - upZ * zAxis[X];
	yAxis[Y] = m_camera.worldup[Y] - upZ * zAxis[Y];
	yAxis[Z] = m_camera.worldup[Z] - upZ * zAxis[Z];

	// Y = up' / || up' ||
	if (!NormalizeVector(yAxis)) return GZ_FAILURE;

	// Compute X axis (X = Y x Z = i(YyZz - YzZy) + j(YzZx - YxZz) + k(YxZy - YyZx), in a left-hand coordnate)
	GzCoord xAxis;
	xAxis[X] = yAxis[Y] * zAxis[Z] - yAxis[Z] * zAxis[Y];
	xAxis[Y] = yAxis[Z] * zAxis[X] - yAxis[X] * zAxis[Z];
	xAxis[Z] = yAxis[X] * zAxis[Y] - yAxis[Y] * zAxis[X];

	// Xwi = [	Xx,		Yx,		Zx,		Cx
	//			Xy,		Yy,		Zy,		Cy
	//			Xz,		Yz,		Zz,		Cz
	//			0,		0,		0,		1	]
	// This matrix is T * R, where T is [Cx, Cy, Cz] and R is [X, Y, Z]
	// So, Xwi = T * R
	// Xiw is the inverse of Xwi, inverse(R) * inverse(T)
	// The translation column elements become dot-products of two vectors during concatenation
	// Xiw = [	Xx,		Xy,		Xz,		-dot(X,C)
	//			Yx,		Yy,		Yz,		-dot(Y,C)
	//			Zx,		Zy,		Zz,		-dot(Z,C)
	//			0,		0,		0,		1			]

	// Xiw
	GzCoord* pCoords[3] = { &xAxis, &yAxis, &zAxis };
	for (int i = 0; i < sizeof(pCoords) / sizeof(GzCoord*); i++) {
		m_camera.Xiw[i][0] = (*pCoords[i])[X];
		m_camera.Xiw[i][1] = (*pCoords[i])[Y];
		m_camera.Xiw[i][2] = (*pCoords[i])[Z];
		m_camera.Xiw[i][3] = -DotVectors((*pCoords[i]), m_camera.position);
	}

	m_camera.Xiw[3][3] = 1.0;

	// init Ximage - put Xsp at base of stack
	matlevel = 0;
	memcpy(Ximage[matlevel], Xsp, sizeof(GzMatrix));

	// push identity matrix into Xnorm.
	GzMatrix iMat = { 1, 0, 0, 0,
					  0, 1, 0, 0,
					  0, 0, 1, 0,
					  0, 0, 0, 1 };
	memcpy(Xnorm[matlevel], iMat, sizeof(GzMatrix));

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
	if (!matrix || (MATLEVELS - 2 < matlevel)) {
		return GZ_FAILURE;
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			Ximage[matlevel + 1][i][j] = Ximage[matlevel][i][0] * matrix[0][j]
				+ Ximage[matlevel][i][1] * matrix[1][j]
				+ Ximage[matlevel][i][2] * matrix[2][j]
				+ Ximage[matlevel][i][3] * matrix[3][j];
		}
	}

	if (matlevel < 1) {
		// Push Xsp and Xpi as identity [I]
		// When matlevel 0, it means that the stack has one element which is Xsp. So next matrix is Xpi
		GzMatrix iMat = { 1, 0, 0, 0,
						  0, 1, 0, 0,
						  0, 0, 1, 0,
						  0, 0, 0, 1 };

		memcpy(Xnorm[matlevel + 1], iMat, sizeof(GzMatrix));
	} else {
		// Transform norms to image space with separate stack
		//	- Xn used only for N
		//	- A No scaling or translation allowed in Xn
		//		- Make sure Xn is a unitary rotation matrix so the resulting N remains normalized

		// Xn must only contain pure rotations (unitary). This means that 3 upper elements of right column should be 0
		GzMatrix mat;
		memcpy(mat, matrix, sizeof(GzMatrix));
		mat[X][Z + 1] = 0;
		mat[Y][Z + 1] = 0;
		mat[Z][Z + 1] = 0;
		mat[Z + 1][Z + 1] = 1;

		// To ensure Xn is a unitary rotation, pre-process each Xn
		//	- Use knowledge of unitary R matrix properties
		//	- Length of any row/col vector must = 1
		// Compute a scale factor and apply to all elements in the matrix
		//	- use any row/col and compute scale factor (K = (a^2 + b^2 + c^2)^1/2)
		//	- divide all elements of 3x3 R : R' = R / K
		//	- R' is normalized (unitary) rotation matrix

		// Choose arbitrary row or column, in here, choose first row
		float K = sqrtf(mat[X][X] * mat[X][X] + mat[Y][X] * mat[Y][X] + mat[Z][X] * mat[Z][X]);

		// Divide all element by K
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				mat[i][j] *= (1 / K);
			}
		}

		//  Push onto the stack
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				Xnorm[matlevel + 1][i][j] = Xnorm[matlevel][i][0] * mat[0][j]
					+ Xnorm[matlevel][i][1] * mat[1][j]
					+ Xnorm[matlevel][i][2] * mat[2][j]
					+ Xnorm[matlevel][i][3] * mat[3][j];
			}
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

/*
- GzPutAttribute() must accept the following tokens/values:

- GZ_RGB_COLOR					GzColor		default flat-shade color
- GZ_INTERPOLATE				int			shader interpolation mode
- GZ_DIRECTIONAL_LIGHT			GzLight
- GZ_AMBIENT_LIGHT            	GzLight		(ignore direction)
- GZ_AMBIENT_COEFFICIENT		GzColor		Ka reflectance
- GZ_DIFFUSE_COEFFICIENT		GzColor		Kd reflectance
- GZ_SPECULAR_COEFFICIENT       GzColor		Ks coef's
- GZ_DISTRIBUTION_COEFFICIENT   float		spec power
*/
	for (int i = 0; i < numAttributes; i++) {
		if (GZ_RGB_COLOR == nameList[i]) {
			GzColor* pColor = static_cast<GzColor*>(valueList[i]);
			flatcolor[RED] = pColor[0][RED];
			flatcolor[GREEN] = pColor[0][GREEN];
			flatcolor[BLUE] = pColor[0][BLUE];
		} else if (GZ_INTERPOLATE == nameList[i]) {
			interp_mode = (static_cast<int*>(valueList[i]))[0];
		} else if (GZ_DIRECTIONAL_LIGHT == nameList[i]) {
			if (numlights < MAX_LIGHTS) {
				memcpy(&lights[numlights], static_cast<GzLight*>(valueList[i]), sizeof(GzLight));
				NormalizeVector(lights[i].direction);
				numlights++;
			}
		} else if (GZ_AMBIENT_LIGHT == nameList[i]) {
			memcpy(ambientlight.color, (static_cast<GzLight*>(valueList[i]))->color, sizeof(GzColor));
		} else if (GZ_AMBIENT_COEFFICIENT == nameList[i]) {
			GzColor* pColor = static_cast<GzColor*>(valueList[i]);
			Ka[RED] = pColor[0][RED];
			Ka[GREEN] = pColor[0][GREEN];
			Ka[BLUE] = pColor[0][BLUE];
		} else if (GZ_DIFFUSE_COEFFICIENT == nameList[i]) {
			GzColor* pColor = static_cast<GzColor*>(valueList[i]);
			Kd[RED] = pColor[0][RED];
			Kd[GREEN] = pColor[0][GREEN];
			Kd[BLUE] = pColor[0][BLUE];
		} else if (GZ_SPECULAR_COEFFICIENT == nameList[i]) {
			GzColor* pColor = static_cast<GzColor*>(valueList[i]);
			Ks[RED] = pColor[0][RED];
			Ks[GREEN] = pColor[0][GREEN];
			Ks[BLUE] = pColor[0][BLUE];
		} else if (GZ_DISTRIBUTION_COEFFICIENT == nameList[i]) {
			spec = (static_cast<float*>(valueList[i]))[0];
		} else if (GZ_TEXTURE_MAP == nameList[i]) {
			tex_fun = static_cast<GzTexture>(valueList[i]);
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
-- Return error code
*/
/*
-- Xform positions of verts using matrix on top of stack 
-- Clip - just discard any triangle with any vert(s) behind view plane 
		- optional: test for triangles with all three verts off-screen (trivial frustum cull)
-- invoke triangle rasterizer  
*/
	GzCoord screenVertices[3];
	GzCoord imageVertices[3];
	GzCoord imageNormals[3];

	GzTextureIndex uvList[3];
	for (int i = 0; i < numParts; i++) {
		if (GZ_POSITION == nameList[i]) {
			// Copy model space vertices into matrix variable.
			memcpy(screenVertices, static_cast<GzCoord*>(valueList[i]), sizeof(GzCoord) * 3);
			memcpy(imageVertices, static_cast<GzCoord*>(valueList[i]), sizeof(GzCoord) * 3);
		} else if (GZ_NORMAL == nameList[i]) {
			// Copy model space normals into matrix variable.
			memcpy(imageNormals, static_cast<GzCoord*>(valueList[i]), sizeof(GzCoord) * 3);
		} else if (GZ_TEXTURE_INDEX == nameList[i]) {
			memcpy(uvList, static_cast<GzTextureIndex*>(valueList[i]), sizeof(GzTextureIndex) * 3);
		}
	}

	for (int j = 0; j < 3; j++) {
		// Transform model space vertices into screen space.
		if (!homogeneousMultiplyMatrixToVector(Ximage[matlevel], screenVertices[j])) {
			return GZ_FAILURE;
		}

		// Do not draw a trangle when the one of them is behind of the view plane
		if (screenVertices[j][Z] <= 0.0) return GZ_SUCCESS;

		// Transform model space vertices into image space
		if (!homogeneousMultiplyMatrixToVector(Xnorm[matlevel], imageVertices[j])) {
			return GZ_FAILURE;
		}

		// Transform model space normals into image space
		if (!homogeneousMultiplyMatrixToVector(Xnorm[matlevel], imageNormals[j])) {
			return GZ_FAILURE;
		}
	}

	// sort screenVertices
	if (screenVertices[0][Y] > screenVertices[2][Y] || (screenVertices[0][Y] == screenVertices[2][Y] && screenVertices[0][X] > screenVertices[2][Y])) {
		SwapMemory(screenVertices[0], screenVertices[2], sizeof(GzCoord));

		// Swap other vertices corresponding to screenVertices
		SwapMemory(imageVertices[0], imageVertices[2], sizeof(GzCoord));
		SwapMemory(imageNormals[0], imageNormals[2], sizeof(GzCoord));

		SwapMemory(uvList[0], uvList[2], sizeof(GzTextureIndex));
	}

	if (screenVertices[1][Y] > screenVertices[2][Y] || (screenVertices[1][Y] == screenVertices[2][Y] && screenVertices[1][X] > screenVertices[2][Y])) {
		SwapMemory(screenVertices[1], screenVertices[2], sizeof(GzCoord));

		// Swap other vertices corresponding to screenVertices
		SwapMemory(imageVertices[1], imageVertices[2], sizeof(GzCoord));
		SwapMemory(imageNormals[1], imageNormals[2], sizeof(GzCoord));

		SwapMemory(uvList[1], uvList[2], sizeof(GzTextureIndex));
	}

	if (screenVertices[0][Y] > screenVertices[1][Y] || (screenVertices[0][Y] == screenVertices[1][Y] && screenVertices[0][X] > screenVertices[1][Y])) {
		SwapMemory(screenVertices[0], screenVertices[1], sizeof(GzCoord));

		// Swap other verticescorresponding to screenVertices
		SwapMemory(imageVertices[0], imageVertices[1], sizeof(GzCoord));
		SwapMemory(imageNormals[0], imageNormals[1], sizeof(GzCoord));

		SwapMemory(uvList[0], uvList[1], sizeof(GzTextureIndex));
	}

	// Define color variables for each vertex used in Gouraud shading
	GzColor vertColors[3];

	// Define color plane in order to interpolate each color for Gouraud shading
	float colorPlane[3][4] = { 0., };

	// Define normal plane in order to interpolate each normal for Phong shading
	float normalPlane[3][4] = { 0., };

	// Define a normal for GZ_FLAT mode, which is in image space
	GzCoord imageTrangeFaceNormal;

	// A general 3D plane equation has four terms: Ax + By + Cz + D = 0
	// The Cross-product of two tri edges produces (A,B,C) vector
	// Create edge vectors (X,Y,Z)0 and (X,Y,Z)1 from any two pairs of vert
	// (X,Y,Z)0 x (X,Y,Z)1 = (A,B,C) = norm to plane of edges (and tri)
	// Use (A,B,C) and plug any vertex coord into a general 3D plane equation and solve for D

	int interpolate_head[2] = { 1, 2 };
	int interpolate_tail[2] = { 0, 1 };

	// Define vectors for interpolation
	GzCoord interpolateVectors[2];
	for (int i = 0; i < sizeof(interpolateVectors) / sizeof(GzCoord); i++) {
		interpolateVectors[i][X] = screenVertices[interpolate_head[i]][X] - screenVertices[interpolate_tail[i]][X];
		interpolateVectors[i][Y] = screenVertices[interpolate_head[i]][Y] - screenVertices[interpolate_tail[i]][Y];
		interpolateVectors[i][Z] = screenVertices[interpolate_head[i]][Z] - screenVertices[interpolate_tail[i]][Z];
	}

	const int A = 0;
	const int B = 1;
	const int C = 2;
	const int D = 3;
	float plane[4] = { 0., };

	// Compute cross product of two vector
	CrossVector(interpolateVectors[0], interpolateVectors[1], plane);
	plane[D] = -(plane[A] * screenVertices[0][X] + plane[B] * screenVertices[0][Y] + plane[C] * screenVertices[0][Z]);

	// Define texture plane in order to interpolate each texture coordination
	float texturePlane[2][4];
	if (true) {
		GzTextureIndex perspTextureCoords[3];

		for (int i = 0; i < 3; i++) {
			// Ax + By + Cz + D = 0 => z = -(Ax + By + D) / C
			float interpolatedZ = -(plane[A] * screenVertices[i][X] + plane[B] * screenVertices[i][Y] + plane[D]) / plane[C];
			perspTextureCoords[i][U] = uvList[i][U] / ((interpolatedZ / (INT_MAX - interpolatedZ)) + 1);
			perspTextureCoords[i][V] = uvList[i][V] / ((interpolatedZ / (INT_MAX - interpolatedZ)) + 1);
		}

		for (int i = 0; i < 2; i++) {
			interpolateVectors[0][Z] = perspTextureCoords[interpolate_head[0]][i] - perspTextureCoords[interpolate_tail[0]][i];
			interpolateVectors[1][Z] = perspTextureCoords[interpolate_head[1]][i] - perspTextureCoords[interpolate_tail[1]][i];

			CrossVector(interpolateVectors[0], interpolateVectors[1], texturePlane[i]);
			texturePlane[i][D] = -(texturePlane[i][A] * screenVertices[0][X] + texturePlane[i][B] * screenVertices[0][Y] + texturePlane[i][C] * perspTextureCoords[0][i]);
		}
	}

	switch (interp_mode) {
		case GZ_FLAT: {
			// Define vectors for interpolation
			GzCoord vec[2];
			for (int i = 0; i < sizeof(interpolateVectors) / sizeof(GzCoord); i++) {
				vec[i][X] = imageVertices[interpolate_head[i]][X] - imageVertices[interpolate_tail[i]][X];
				vec[i][Y] = imageVertices[interpolate_head[i]][Y] - imageVertices[interpolate_tail[i]][Y];
				vec[i][Z] = imageVertices[interpolate_head[i]][Z] - imageVertices[interpolate_tail[i]][Z];
			}

			CrossVector(vec[0], vec[1], imageTrangeFaceNormal);
		} break;

		case GZ_COLOR:
			for (int i = 0; i < 3; i++) {
				// Compute the color at each vertex
				ShadingColor(imageNormals[i], numlights, lights, ambientlight, spec, Ks, Kd, Ka, uvList[i], tex_fun, interp_mode, vertColors[i]);
			}

			for (int i = 0; i < 3; i++) {
				// Set up color interpolation 
				interpolateVectors[0][Z] = vertColors[interpolate_head[0]][i] - vertColors[interpolate_tail[0]][i];
				interpolateVectors[1][Z] = vertColors[interpolate_head[1]][i] - vertColors[interpolate_tail[1]][i];

				CrossVector(interpolateVectors[0], interpolateVectors[1], colorPlane[i]);
				colorPlane[i][D] = -(colorPlane[i][A] * screenVertices[0][X] + colorPlane[i][B] * screenVertices[0][Y] + colorPlane[i][C] * vertColors[0][i]);
			}
			break;

		case GZ_NORMALS:
			for (int i = 0; i < 3; i++) {
				interpolateVectors[0][Z] = imageNormals[interpolate_head[0]][i] - imageNormals[interpolate_tail[0]][i];
				interpolateVectors[1][Z] = imageNormals[interpolate_head[1]][i] - imageNormals[interpolate_tail[1]][i];

				CrossVector(interpolateVectors[0], interpolateVectors[1], normalPlane[i]);
				normalPlane[i][D] = -(normalPlane[i][A] * screenVertices[0][X] + normalPlane[i][B] * screenVertices[0][Y] + normalPlane[i][C] * imageNormals[0][i]);
			}
			break;

		default:
			break;
	}

	// Setup edges for scan-line rasterization 
	const int SLOPE_X = 0; // dX/dY (X2-X1) / (Y2-Y1)
	const int SLOPE_Z = 1; // dZ/dY (Z2-Z1) / (Y2-Y1)
	const int SLOPE_NORM_X = 2; // norm_dX/dy
	const int SLOPE_NORM_Y = 3; // norm_dY/dy
	const int SLOPE_NORM_Z = 4; // norm_dZ/dy
	const int SLOPE_R = 5; // color_dR/dy
	const int SLOPE_G = 6; // color_dG/dy
	const int SLOPE_B = 7; // color_dB/dy

	int head[3] = { 1, 2, 2 };
	int tail[3] = { 0, 0, 1 };
	float edges[3][8] = { 0., };
	for (int i = 0; i < sizeof(edges) / sizeof(edges[i]); i++) {
		edges[i][SLOPE_X] = (screenVertices[head[i]][X] - screenVertices[tail[i]][X]) / (screenVertices[head[i]][Y] - screenVertices[tail[i]][Y]);
		edges[i][SLOPE_Z] = (screenVertices[head[i]][Z] - screenVertices[tail[i]][Z]) / (screenVertices[head[i]][Y] - screenVertices[tail[i]][Y]);
		edges[i][SLOPE_NORM_X] = (imageNormals[head[i]][X] - imageNormals[tail[i]][X]) / (screenVertices[head[i]][Y] - screenVertices[tail[i]][Y]);
		edges[i][SLOPE_NORM_Y] = (imageNormals[head[i]][Y] - imageNormals[tail[i]][Y]) / (screenVertices[head[i]][Y] - screenVertices[tail[i]][Y]);
		edges[i][SLOPE_NORM_Z] = (imageNormals[head[i]][Z] - imageNormals[tail[i]][Z]) / (screenVertices[head[i]][Y] - screenVertices[tail[i]][Y]);
		edges[i][SLOPE_R] = (vertColors[head[i]][RED] - vertColors[tail[i]][RED]) / (screenVertices[head[i]][Y] - screenVertices[tail[i]][Y]);
		edges[i][SLOPE_G] = (vertColors[head[i]][GREEN] - vertColors[tail[i]][GREEN]) / (screenVertices[head[i]][Y] - screenVertices[tail[i]][Y]);
		edges[i][SLOPE_B] = (vertColors[head[i]][BLUE] - vertColors[tail[i]][BLUE]) / (screenVertices[head[i]][Y] - screenVertices[tail[i]][Y]);
	}

	// Decide left edge and right edge.
	int leftStartVert = 0;
	int rightStartVert = 0;
	int leftEdge = 0; // dX/dy < dX/dY
	int rightEdge = 1; // dX/dY > dX/dY

	if (screenVertices[0][Y] == screenVertices[1][Y]) {
		rightStartVert = 1;
		leftEdge = 1;
		rightEdge = 2;
	} else {
		if (edges[1][SLOPE_X] < edges[0][SLOPE_X]) {
			leftEdge = 1;
			rightEdge = 0;
		}
	}

	// Rasterize with scan-line
	for (int i = 0; i < 2; i++) {
		// top/left rule -> Pixeling when the pixel is on the top, left, or in the triangle.
		float topY = ceil(screenVertices[0][Y]);
		if (screenVertices[0][Y] != screenVertices[1][Y] && (screenVertices[0][Y] == topY)) {
			topY += 1.0;
		}

		for (int row = static_cast<int>(topY); row < static_cast<int>(ceil(screenVertices[i + 1][Y])); row++) {
			float dY = static_cast<float>(row) - screenVertices[i][Y];
			float leftX = screenVertices[leftStartVert][X] + (edges[leftEdge][SLOPE_X] * dY);
			float rightX = screenVertices[rightStartVert][X] + (edges[rightEdge][SLOPE_X] * dY);
			float leftZ = screenVertices[leftStartVert][Z] + (edges[leftEdge][SLOPE_Z] * dY);
			float rightZ = screenVertices[rightStartVert][Z] + (edges[rightEdge][SLOPE_Z] * dY);
			float slopeZ = (rightZ - leftZ) / (rightX - leftX);

			for (int col = static_cast<int>(ceil(leftX)); col < static_cast<int>(ceil(rightX)); col++) {
				float dX = static_cast<float>(col) - leftX;
				GzDepth zVal = static_cast<GzDepth>(leftZ + slopeZ * dX);

				// Ax + By + Cz + D = 0 => z = -(Ax + By + D) / C
				float interpolatedZ = -(plane[A] * col + plane[B] * row + plane[D]) / plane[C];

				GzPixel pixel;
				if (GZ_FAILURE == GzGet(col, row, &pixel.red, &pixel.green, &pixel.blue, &pixel.alpha, &pixel.z)) {
					continue;
				}

				if (!(zVal < pixel.z)) continue;

				GzColor color;
				switch (interp_mode) {
					case GZ_FLAT: // Flat Shading
								  // Nomalize normal
						NormalizeVector(imageTrangeFaceNormal);

						ShadingColor(imageTrangeFaceNormal, numlights, lights, ambientlight, spec, Ks, Kd, Ka, NULL, NULL, interp_mode, color);
						break;

					case GZ_COLOR: // Gourad Shading - interpolate color 
						color[RED] = -(colorPlane[RED][A] * col + colorPlane[RED][B] * row + colorPlane[RED][D]) / colorPlane[RED][C];
						color[GREEN] = -(colorPlane[GREEN][A] * col + colorPlane[GREEN][B] * row + colorPlane[GREEN][D]) / colorPlane[GREEN][C];
						color[BLUE] = -(colorPlane[BLUE][A] * col + colorPlane[BLUE][B] * row + colorPlane[BLUE][D]) / colorPlane[BLUE][C];

						if (tex_fun) {
							GzTextureIndex interpolatedTextureCoord;
							interpolatedTextureCoord[U] = -(texturePlane[U][A] * col + texturePlane[U][B] * row + texturePlane[U][D]) / texturePlane[U][C];
							interpolatedTextureCoord[V] = -(texturePlane[V][A] * col + texturePlane[V][B] * row + texturePlane[V][D]) / texturePlane[V][C];

							interpolatedTextureCoord[U] = interpolatedTextureCoord[U] * ((interpolatedZ / (INT_MAX - interpolatedZ)) + 1);
							interpolatedTextureCoord[V] = interpolatedTextureCoord[V] * ((interpolatedZ / (INT_MAX - interpolatedZ)) + 1);

							GzColor textureColor;
							tex_fun(interpolatedTextureCoord[U], interpolatedTextureCoord[V], textureColor);
							color[RED] *= textureColor[RED];
							color[GREEN] *= textureColor[GREEN];
							color[BLUE] *= textureColor[BLUE];
						}
						break;

					case GZ_NORMAL: { // Phong Shading - interpolate normal and compute color
						GzCoord interpolatedNormal;
						interpolatedNormal[X] = -(normalPlane[X][A] * col + normalPlane[X][B] * row + normalPlane[X][D]) / normalPlane[X][C];
						interpolatedNormal[Y] = -(normalPlane[Y][A] * col + normalPlane[Y][B] * row + normalPlane[Y][D]) / normalPlane[Y][C];
						interpolatedNormal[Z] = -(normalPlane[Z][A] * col + normalPlane[Z][B] * row + normalPlane[Z][D]) / normalPlane[Z][C];

						// Nomalize normal
						NormalizeVector(interpolatedNormal);

						if (!tex_fun) {
							// Compute color
							ShadingColor(interpolatedNormal, numlights, lights, ambientlight, spec, Ks, Kd, Ka, NULL, NULL, interp_mode, color);
							break;
						}

						GzTextureIndex interpolatedTextureCoord;
						interpolatedTextureCoord[U] = -(texturePlane[U][A] * col + texturePlane[U][B] * row + texturePlane[U][D]) / texturePlane[U][C];
						interpolatedTextureCoord[V] = -(texturePlane[V][A] * col + texturePlane[V][B] * row + texturePlane[V][D]) / texturePlane[V][C];

						interpolatedTextureCoord[U] = interpolatedTextureCoord[U] * ((interpolatedZ / (INT_MAX - interpolatedZ)) + 1);
						interpolatedTextureCoord[V] = interpolatedTextureCoord[V] * ((interpolatedZ / (INT_MAX - interpolatedZ)) + 1);
						
						// Compute color
						ShadingColor(interpolatedNormal, numlights, lights, ambientlight, spec, Ks, Kd, Ka, interpolatedTextureCoord, tex_fun, interp_mode, color);
					} break;

					default:
						break;
				}

				GzPut(col, row, ctoi(color[RED]), ctoi(color[GREEN]), ctoi(color[BLUE]), 1, zVal);
			}
		}

		if (i < 1 && screenVertices[0][Y] != screenVertices[1][Y] && screenVertices[1][Y] != screenVertices[2][Y]) {
			GzCoord midVertex;
			midVertex[Y] = screenVertices[1][Y];
			float dY = screenVertices[1][Y] - screenVertices[0][Y];
			midVertex[X] = screenVertices[0][X] + (edges[1][SLOPE_X] * dY);
			midVertex[Z] = screenVertices[0][Z] + (edges[1][SLOPE_Z] * dY);

			if (midVertex[X] < screenVertices[1][X]) {
				memcpy(screenVertices[0], midVertex, sizeof(GzCoord));
				leftEdge = 1;
				rightEdge = 2;
			} else {
				memcpy(screenVertices[0], screenVertices[1], sizeof(GzCoord));
				memcpy(screenVertices[1], midVertex, sizeof(GzCoord));

				leftEdge = 2;
				rightEdge = 1;
			}

			rightStartVert = 1;
		} else {
			break;
		}
	}

	return GZ_SUCCESS;
}

