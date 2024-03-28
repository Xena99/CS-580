/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"

GzColor	*image;
int xs;
int ys;
int reset = 1;

/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
	if (reset)
	{
		FILE *file = fopen("tex.ppm", "rb");
		if (!file) {
			exit(EXIT_FAILURE);
		}

		char header[3];
		fscanf(file, "%2s\n", header);
		if (strcmp(header, "P6") != 0) {
			fprintf(stderr, "Unsupported texture format (only binary PPM supported)\n");
			exit(EXIT_FAILURE);
		}
		OutputDebugString(header);
		int maxval;
		fscanf(file, "\n%d %d\n%d\n", &xs, &ys, &maxval);
		if (maxval != 255) {
			fprintf(stderr, "Unsupported max value in texture file (must be 255)\n");
			exit(EXIT_FAILURE);
		}
		
		image = (GzColor *)malloc(xs * ys * sizeof(GzColor));
		if (!*image) {
			fprintf(stderr, "Failed to allocate memory for texture image\n");
			exit(EXIT_FAILURE);
		}

		for (int i = 0; i < xs * ys; i++) {
			unsigned char pixel[3];
			fread(pixel, sizeof(pixel), 1, file);
			(image)[i][0] = (float)pixel[0] / 255.0f;
			(image)[i][1] = (float)pixel[1] / 255.0f;
			(image)[i][2] = (float)pixel[2] / 255.0f;
		}

		reset = 0;          /* init is done */
		fclose(file);
	}
	// Wrap texture coordinates
	u = fmod(u, 1.0f);
	v = fmod(v, 1.0f);
	if (u < 0) u += 1.0f; // Ensure positive values
	if (v < 0) v += 1.0f;

	// Compute the location on the texture map
	float xLocation = u * (xs - 1);
	float yLocation = v * (ys - 1);

	// Compute the four surrounding pixels
	int p00x = static_cast<int>(floor(xLocation));
	int p00y = static_cast<int>(floor(yLocation));
	int p10x = p00x + 1 < xs ? p00x + 1 : p00x; // Clamp to edge of texture
	int p10y = p00y;
	int p01x = p00x;
	int p01y = p00y + 1 < ys ? p00y + 1 : p00y; // Clamp to edge of texture
	int p11x = p10x;
	int p11y = p01y;

	// Fractional part of the location
	float xFrac = xLocation - p00x;
	float yFrac = yLocation - p00y;

	// Get the colors of the four pixels
	GzColor p00Color, p10Color, p01Color, p11Color;
	memcpy(p00Color, image[p00y * xs + p00x], sizeof(GzColor));
	memcpy(p10Color, image[p10y * xs + p10x], sizeof(GzColor));
	memcpy(p01Color, image[p01y * xs + p01x], sizeof(GzColor));
	memcpy(p11Color, image[p11y * xs + p11x], sizeof(GzColor));

	// Linear interpolations
	GzColor p0010RGB, p0111RGB, pOutputRGB;
	for (int i = 0; i < 3; i++) {
		p0010RGB[i] = (1 - xFrac) * p00Color[i] + xFrac * p10Color[i];
		p0111RGB[i] = (1 - xFrac) * p01Color[i] + xFrac * p11Color[i];
		pOutputRGB[i] = (1 - yFrac) * p0010RGB[i] + yFrac * p0111RGB[i];
	}

	// Assign the interpolated color to output
	memcpy(color, pOutputRGB, sizeof(GzColor));

	return GZ_SUCCESS;

}



/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color)
{
	return 0;
}

