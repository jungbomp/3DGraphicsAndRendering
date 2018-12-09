/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"
#include	"math.h"

GzColor	*image = NULL;
int xs, ys;
int reset = 1;

#define STONE_THRESHOLD 0.02f

const int SAMPLE_SIZE = 1024;
bool gIsInit = false;
int p[SAMPLE_SIZE + SAMPLE_SIZE + 2];
float g3[SAMPLE_SIZE + SAMPLE_SIZE + 2][3];
float g2[SAMPLE_SIZE + SAMPLE_SIZE + 2][2];
float g1[SAMPLE_SIZE + SAMPLE_SIZE + 2];

void addSamples(long xi, long yi, long zi, long maxOrder, float at[3], float *F, unsigned long *ID)
{
	float dx, dy, dz, fx, fy, fz, d2;
	long count, i, j, index;
	unsigned long seed, this_id;

	char poissonCount[256] =	{
		4,3,1,1,1,2,4,2,2,2,5,1,0,2,1,2,2,0,4,3,2,1,2,1,3,2,2,4,2,2,5,1,2,3,2,2,2,2,2,3,
		2,4,2,5,3,2,2,2,5,3,3,5,2,1,3,3,4,4,2,3,0,4,2,2,2,1,3,2,2,2,3,3,3,1,2,0,2,1,1,2,
		2,2,2,5,3,2,3,2,3,2,2,1,0,2,1,1,2,1,2,2,1,3,4,2,2,2,5,4,2,4,2,2,5,4,3,2,2,5,4,3,
		3,3,5,2,2,2,2,2,3,1,1,4,2,1,3,3,4,3,2,4,3,3,3,4,5,1,4,2,4,3,1,2,3,5,3,2,1,3,1,3,
		3,3,2,3,1,5,5,4,2,2,4,1,3,4,1,5,3,3,5,3,4,3,2,2,1,1,1,1,1,2,4,5,4,5,4,2,1,5,1,1,
		2,3,3,3,2,5,2,3,3,2,0,2,1,1,4,2,1,3,2,1,2,2,3,2,5,5,3,4,5,5,2,4,4,5,3,2,2,2,1,4,
		2,3,3,4,2,5,4,2,4,2,2,2,4,5,3,2
	};

	seed = 702395077 * xi + 915488749 * yi + 2120969693 * zi;
	count = poissonCount[seed >> 24];
	seed = 1402024253 * seed + 586950981;

	float temp = 1.0 / 4294967296.0;
	this_id = seed;
	for (j = 0; j < count; j++)  {
		seed = 1402024253 * seed + 586950981;
		fx = (seed + 0.5) * temp;
		seed = 1402024253 * seed + 586950981;
		fy = (seed + 0.5) * temp;
		seed = 1402024253 * seed + 586950981;
		fz = (seed + 0.5) * temp;
		seed = 1402024253 * seed + 586950981;

		dx = xi + fx - at[0];
		dy = yi + fy - at[1];
		dz = zi + fz - at[2];
		d2 = dx * dx + dy * dy + dz * dz;

		if (d2 < F[maxOrder - 1]) {
			index = maxOrder;
			while (index > 0 && d2 < F[index - 1]) index--;

			for (i = maxOrder - 2; i >= index; i--) {
				F[i + 1] = F[i];
				ID[i + 1] = ID[i];
			}

			F[index] = d2;
			ID[index] = this_id;
		}
	}
}

void WorleyNoise3D(float at[3], long maxOrder, float *F, unsigned long *ID)
{
	float x2, y2, z2, mx2, my2, mz2;
	float newAt[3];
	long intAt[3], i;

	const float DENSITY = 0.398150;

	for (i = 0; i < maxOrder; i++) F[i] = 999999.9;

	newAt[0] = DENSITY * at[0];
	newAt[1] = DENSITY * at[1];
	newAt[2] = DENSITY * at[2];

	intAt[0] = int(floor(newAt[0]));
	intAt[1] = int(floor(newAt[1]));
	intAt[2] = int(floor(newAt[2]));

	addSamples(intAt[0], intAt[1], intAt[2], maxOrder, newAt, F, ID);

	x2 = newAt[0] - intAt[0];
	y2 = newAt[1] - intAt[1];
	z2 = newAt[2] - intAt[2];
	mx2 = (1.0 - x2) * (1.0 - x2);
	my2 = (1.0 - y2) * (1.0 - y2);
	mz2 = (1.0 - z2) * (1.0 - z2);
	x2 *= x2;
	y2 *= y2;
	z2 *= z2;

	if (x2 < F[maxOrder - 1])  addSamples(intAt[0] - 1, intAt[1], intAt[2], maxOrder, newAt, F, ID);
	if (y2 < F[maxOrder - 1])  addSamples(intAt[0], intAt[1] - 1, intAt[2], maxOrder, newAt, F, ID);
	if (z2 < F[maxOrder - 1])  addSamples(intAt[0], intAt[1], intAt[2] - 1, maxOrder, newAt, F, ID);

	if (mx2 < F[maxOrder - 1]) addSamples(intAt[0] + 1, intAt[1], intAt[2], maxOrder, newAt, F, ID);
	if (my2 < F[maxOrder - 1]) addSamples(intAt[0], intAt[1] + 1, intAt[2], maxOrder, newAt, F, ID);
	if (mz2 < F[maxOrder - 1]) addSamples(intAt[0], intAt[1], intAt[2] + 1, maxOrder, newAt, F, ID);

	if (x2 + y2 < F[maxOrder - 1]) addSamples(intAt[0] - 1, intAt[1] - 1, intAt[2], maxOrder, newAt, F, ID);
	if (x2 + z2 < F[maxOrder - 1]) addSamples(intAt[0] - 1, intAt[1], intAt[2] - 1, maxOrder, newAt, F, ID);
	if (y2 + z2 < F[maxOrder - 1]) addSamples(intAt[0], intAt[1] - 1, intAt[2] - 1, maxOrder, newAt, F, ID);
	if (mx2 + my2 < F[maxOrder - 1]) addSamples(intAt[0] + 1, intAt[1] + 1, intAt[2], maxOrder, newAt, F, ID);
	if (mx2 + mz2 < F[maxOrder - 1]) addSamples(intAt[0] + 1, intAt[1], intAt[2] + 1, maxOrder, newAt, F, ID);
	if (my2 + mz2 < F[maxOrder - 1]) addSamples(intAt[0], intAt[1] + 1, intAt[2] + 1, maxOrder, newAt, F, ID);
	if (x2 + my2 < F[maxOrder - 1]) addSamples(intAt[0] - 1, intAt[1] + 1, intAt[2], maxOrder, newAt, F, ID);
	if (x2 + mz2 < F[maxOrder - 1]) addSamples(intAt[0] - 1, intAt[1], intAt[2] + 1, maxOrder, newAt, F, ID);
	if (y2 + mz2 < F[maxOrder - 1]) addSamples(intAt[0], intAt[1] - 1, intAt[2] + 1, maxOrder, newAt, F, ID);
	if (mx2 + y2 < F[maxOrder - 1]) addSamples(intAt[0] + 1, intAt[1] - 1, intAt[2], maxOrder, newAt, F, ID);
	if (mx2 + z2 < F[maxOrder - 1]) addSamples(intAt[0] + 1, intAt[1], intAt[2] - 1, maxOrder, newAt, F, ID);
	if (my2 + z2 < F[maxOrder - 1]) addSamples(intAt[0], intAt[1] + 1, intAt[2] - 1, maxOrder, newAt, F, ID);

	if (x2 + y2 + z2 < F[maxOrder - 1]) addSamples(intAt[0] - 1, intAt[1] - 1, intAt[2] - 1, maxOrder, newAt, F, ID);
	if (x2 + y2 + mz2 < F[maxOrder - 1]) addSamples(intAt[0] - 1, intAt[1] - 1, intAt[2] + 1, maxOrder, newAt, F, ID);
	if (x2 + my2 + z2 < F[maxOrder - 1]) addSamples(intAt[0] - 1, intAt[1] + 1, intAt[2] - 1, maxOrder, newAt, F, ID);
	if (x2 + my2 + mz2 < F[maxOrder - 1]) addSamples(intAt[0] - 1, intAt[1] + 1, intAt[2] + 1, maxOrder, newAt, F, ID);
	if (mx2 + y2 + z2 < F[maxOrder - 1]) addSamples(intAt[0] + 1, intAt[1] - 1, intAt[2] - 1, maxOrder, newAt, F, ID);
	if (mx2 + y2 + mz2 < F[maxOrder - 1]) addSamples(intAt[0] + 1, intAt[1] - 1, intAt[2] + 1, maxOrder, newAt, F, ID);
	if (mx2 + my2 + z2 < F[maxOrder - 1]) addSamples(intAt[0] + 1, intAt[1] + 1, intAt[2] - 1, maxOrder, newAt, F, ID);
	if (mx2 + my2 + mz2 < F[maxOrder - 1]) addSamples(intAt[0] + 1, intAt[1] + 1, intAt[2] + 1, maxOrder, newAt, F, ID);

	for (i = 0; i < maxOrder; i++) {
		F[i] = sqrt(F[i]) * (1.0 / DENSITY);
	}
}

float GetPerlinNoise(float vec[2])
{
	int bx0, bx1, by0, by1, b00, b10, b01, b11;
	float rx0, rx1, ry0, ry1, *q, sx, sy, a, b, t, u, v;
	int i, j;

	const int SEED = 94;

	if (!gIsInit) {
		srand(SEED);

		int i, j, k;

		for (i = 0; i < SAMPLE_SIZE; i++) {
			p[i] = i;
			g1[i] = (float)((rand() % (SAMPLE_SIZE + SAMPLE_SIZE)) - SAMPLE_SIZE) / SAMPLE_SIZE;
			for (j = 0; j < 2; j++)
				g2[i][j] = (float)((rand() % (SAMPLE_SIZE + SAMPLE_SIZE)) - SAMPLE_SIZE) / SAMPLE_SIZE;

			float s = (float)sqrt(g2[i][0] * g2[i][0] + g2[i][1] * g2[i][1]);
			g2[i][0] = g2[i][0] * (1.0 / s);
			g2[i][1] = g2[i][1] * (1.0 / s);

			for (j = 0; j < 3; j++)
				g3[i][j] = (float)((rand() % (SAMPLE_SIZE + SAMPLE_SIZE)) - SAMPLE_SIZE) / SAMPLE_SIZE;

			s = (float)sqrt(g3[i][0] * g3[i][0] + g3[i][1] * g3[i][1] + g3[i][2] * g3[i][2]);
			g3[i][0] = g3[i][0] * (1.0 / s);
			g3[i][1] = g3[i][1] * (1.0 / s);
			g3[i][2] = g3[i][2] * (1.0 / s);
		}

		while (--i) {
			k = p[i];
			p[i] = p[j = rand() % SAMPLE_SIZE];
			p[j] = k;
		}

		for (i = 0; i < SAMPLE_SIZE + 2; i++) {
			p[SAMPLE_SIZE + i] = p[i];
			g1[SAMPLE_SIZE + i] = g1[i];
			for (j = 0; j < 2; j++)
				g2[SAMPLE_SIZE + i][j] = g2[i][j];
			for (j = 0; j < 3; j++)
				g3[SAMPLE_SIZE + i][j] = g3[i][j];
		}

		gIsInit = true;
	}
	
	t = vec[0] + 0x1000;
	bx0 = ((int)t) & (SAMPLE_SIZE-1);
	bx1 = (bx0 + 1) & (SAMPLE_SIZE - 1);
	rx0 = t - (int)t;
	rx1 = rx0 - 1.0f;

	t = vec[1] + 0x1000;
	by0 = ((int)t) & (SAMPLE_SIZE - 1);
	by1 = (by0 + 1) & (SAMPLE_SIZE-1);
	ry0 = t - (int)t;
	ry1 = ry0 - 1.0f;

	i = p[bx0];
	j = p[bx1];

	b00 = p[i + by0];
	b10 = p[j + by0];
	b01 = p[i + by1];
	b11 = p[j + by1];

	sx = rx0 * rx0 * (3.0f - 2.0f * rx0);
	sy = ry0 * ry0 * (3.0f - 2.0f * ry0);

	q = g2[b00];
	u = (rx0 * q[0] + ry0 * q[1]);
	q = g2[b10];
	v = (rx1 * q[0] + ry0 * q[1]);
	a = u + sx * (v - u);

	q = g2[b01];
	u = (rx0 * q[0] + ry1 * q[1]);
	q = g2[b11];
	v = (rx1 * q[0] + ry1 * q[1]);
	b = u + sx * (v - u);

	return (a + sy * (b - a));
}

float PerlinNoise2D(float vec[2], int octaves, float frequency, float amplitude) {
	int terms = octaves;
	float freq = frequency;
	float result = 0.0f;
	float amp = amplitude;

	vec[0] *= frequency;
	vec[1] *= frequency;

	for (int i = 0; i<terms; i++)
	{
		result += GetPerlinNoise(vec)*amp;
		vec[0] *= 2.0f;
		vec[1] *= 2.0f;
		amp *= 0.5f;
	}

	return result;
}

/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
	unsigned char pixel[3];
	unsigned char dummy;
	char foo[8];
	int i, j;
	FILE *fd;

	if (reset) {          /* open and load texture file */
		fd = fopen("texture", "rb");
		if (fd == NULL) {
			fprintf(stderr, "texture file not found\n");
			exit(-1);
		}

		fscanf(fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
		image = (GzColor*)malloc(sizeof(GzColor)*(xs + 1)*(ys + 1));
		if (image == NULL) {
			fprintf(stderr, "malloc for texture image failed\n");
			exit(-1);
		}

		for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */
			fread(pixel, sizeof(pixel), 1, fd);
			image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
			image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
			image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
		}

		reset = 0;          /* init is done */
		fclose(fd);
	}

	/* bounds-test u,v to make sure nothing will overflow image array bounds */
	/* determine texture cell corner values and perform bilinear interpolation */
	/* set color to interpolated GzColor value and return */

	// Cut the u and v value when they exceed 1 or below 0
	u = max(0, min(1, u));
	v = max(0, min(1, v));
	
	// Do texture look up at a screen pixel: u,v -> Color = f(u,v)
	// Scale u,v range to x, y texture image size (x-size, y-size)
	// Scaling is done with multiply of pixel u,v by the size of the texture image 
	// (x-size - 1) and (y-size - 1) to get a float coord (x, y), which is
	// a sample point in the texture image
	float texelX = u * (xs - 1);
	float texelY = v * (ys - 1);

	// Color for the triangle pixel is found in the texture map image,
	// but the sample point lies between texels (pixels in the texture map)
	float left = floor(texelX);
	float right = ceil(texelX);
	float top = floor(texelY);
	float bottom = ceil(texelY);

	/*
	* The variable names below correspond to this layout:
	*	A (left, top)		B (right, top)
	*
	*   D (left, bottom)	C (right, bottom)
	*/
	GzColor A, B, C, D;
	memcpy(A, image[(int)(top * xs + left)], sizeof(GzColor));
	memcpy(B, image[(int)(top * xs + right)], sizeof(GzColor));
	memcpy(C, image[(int)(bottom * xs + right)], sizeof(GzColor));
	memcpy(D, image[(int)(bottom * xs + left)], sizeof(GzColor));

	// Apply simple biliner interpolation filter (a.k.a triangle or pyramid reconstruction filter
	// Color(p) = stC + (1-s)tD + s(1-t)B + (1-s)(1-t)A
	// s, t are fractional distances [0, 1]
	// A, B, C, D are pixel RGB colors at neighboring integer-coord texels
	// E.G., Find Color (2.4, 3.7) where: s = .4   t = .7
	GzColor colorTopRow, colorBottomRow;
	float s = texelX - left;
	float t = texelY - top;

	color[RED] = s * t * C[RED] + (1 - s) * t * D[RED] + s * (1 - t) * B[RED] + (1 - s) * (1 - t) * A[RED];
	color[GREEN] = s * t * C[GREEN] + (1 - s) * t * D[GREEN] + s * (1 - t) * B[GREEN] + (1 - s) * (1 - t) * A[GREEN];
	color[BLUE] = s * t * C[BLUE] + (1 - s) * t * D[BLUE] + s * (1 - t) * B[BLUE] + (1 - s) * (1 - t) * A[BLUE];

	return GZ_SUCCESS;
}

/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color)
{
	if (u < .005 || u > .995 || v < .005 || v > .995) {
		color[RED] = color[BLUE] = color[GREEN] = 0;
		return GZ_SUCCESS;
	}

	float F[2];
	unsigned long ID[2];

	float at[3];
	at[0] = u * 3;
	at[1] = v * 3;
	at[2] = u * v;

	WorleyNoise3D(at, 2, F, ID);
	float distance = F[1] - F[0];
	float colorID = ID[0];

	if (distance > STONE_THRESHOLD) {
		float perlinNoise;
		float vec[2] = { u * 3, v * 3 };
		perlinNoise = (1 + PerlinNoise2D(vec, 4, 4, 1)) / 2;
		color[RED] = color[BLUE] = color[GREEN] = perlinNoise * 0.5f;
	} else {
		color[RED] = color[GREEN] = color[BLUE] = 0;
	}

	return GZ_SUCCESS;
}

/* Free texture memory */
int GzFreeTexture()
{
	if (image != NULL)
		free(image);
	return GZ_SUCCESS;
}

