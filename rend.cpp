#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include "disp.h"
#include <stdio.h>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>


int GzNewRender(GzRender **render, GzRenderClass renderClass, GzDisplay *display)
{
	/*
	- malloc a renderer struct
	- keep closed until BeginRender inits are done
	- span interpolator needs pointer to display for pixel writes
	- check for legal class GZ_Z_BUFFER_RENDER
	*/
	if (renderClass != GZ_Z_BUFFER_RENDER)
	{
		return GZ_FAILURE;
	}
	GzRender *render1 = (GzRender *)malloc(sizeof(GzRender));
	*render = render1;
	(*render)->display = display;
	return GZ_SUCCESS;
}


int GzFreeRender(GzRender *render)
{
	/*
	-free all renderer resources
	*/
	if (render == NULL) {
		return GZ_FAILURE;
	}
	free(render);
	return GZ_SUCCESS;
}


int GzBeginRender(GzRender	*render)
{
	/*
	- set up for start of each frame - init frame buffer
	*/
	if (render == NULL) {
		return GZ_FAILURE;
	}
	render->open = 1;
	render->flatcolor[RED] = 0;
	render->flatcolor[GREEN] = 0;
	render->flatcolor[BLUE] = 0;
	return GZ_SUCCESS;
}


int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList,
	GzPointer *valueList) /* void** valuelist */
{
	/*
	- set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
	- later set shaders, interpolaters, texture maps, and lights
	*/
	if (render == NULL) {
		return GZ_FAILURE;
	}

	for (int i = 0; i < numAttributes; i++) {
		if (nameList[i] == GZ_RGB_COLOR) {
			GzColor *color = (GzColor *)valueList[0];
			if ((*color)[0] < 0) {
				(*color)[0] = 0;
			}
			if ((*color)[1] < 0) {
				(*color)[1] = 0;
			}
			if ((*color)[2] < 0) {
				(*color)[2] = 0;
			}
			if ((*color)[0] > 4095) {
				(*color)[0] = 4095;
			}
			if ((*color)[1] > 4095) {
				(*color)[1] = 4095;
			}
			if ((*color)[2] > 4095) {
				(*color)[2] = 4095;
			}
			render->flatcolor[RED] = (*color)[0];
			render->flatcolor[GREEN] = (*color)[1];
			render->flatcolor[BLUE] = (*color)[2];

		}
	}

	return GZ_SUCCESS;
}


void computeTriangleColor(const GzCoord& normal, GzColor& triColor) {
	// Define the light vector
	
}

short ctoi(float color)		/* convert float color to GzIntensity short */
{
  return(short)((int)(color * ((1 << 12) - 1)));
}

float EdgeFunctionCW(const GzCoord& a, const GzCoord& b, int x, int y) {
	return (b[1] - a[1]) * (x - a[0]) - (b[0] - a[0]) * (y - a[1]);
}

void sortVertices(GzCoord* verticesList) {
	// Sort vertices by Y-coordinate in ascending order
	for (int i = 0; i < 2; ++i) {
		for (int j = i + 1; j < 3; ++j) {
			if (verticesList[i][1] > verticesList[j][1]) {
				// Swap the entire GzCoord
				for (int k = 0; k < 3; ++k) {
					std::swap(verticesList[i][k], verticesList[j][k]);
				}
			}
		}
	}

	// Handle cases where Y-coordinates are approximately equal
	if (std::abs(verticesList[0][1] - verticesList[1][1]) < 0.001f) {
		if (verticesList[0][0] > verticesList[1][0]) {
			for (int k = 0; k < 3; ++k) {
				std::swap(verticesList[0][k], verticesList[1][k]);
			}
		}
	}

	if (std::abs(verticesList[1][1] - verticesList[2][1]) < 0.001f) {
		if (verticesList[1][0] > verticesList[2][0]) {
			for (int k = 0; k < 3; ++k) {
				std::swap(verticesList[1][k], verticesList[2][k]);
			}
		}
	}
}


// Assume GzRender, GzToken, GzPointer, etc. are defined elsewhere
int GzPutTriangle(GzRender* render, int numParts, GzToken* nameList, GzPointer* valueList) {

	GzCoord* verticesList = (GzCoord*)valueList[0];

	sortVertices(verticesList);
	
	//Get Bounding box
	int minx = (int)(min(min(verticesList[0][0], verticesList[1][0]), verticesList[2][0]) + 0.5);
	int maxx = (int)(max(max(verticesList[0][0], verticesList[1][0]), verticesList[2][0]) + 0.5);
	int miny = (int)(min(min(verticesList[0][1], verticesList[1][1]), verticesList[2][1]) + 0.5);
	int maxy = (int)(max(max(verticesList[0][1], verticesList[1][1]), verticesList[2][1]) + 0.5);

	// Iterate over each pixel within the bounding box
	for (int y = miny; y <= maxy; ++y) {
		for (int x = minx; x <= maxx; ++x) {
			// Barycentric coordinates for clockwise order
			float f12 = EdgeFunctionCW(verticesList[1], verticesList[2], x, y);
			float f20 = EdgeFunctionCW(verticesList[2], verticesList[0], x, y);
			float f01 = EdgeFunctionCW(verticesList[0], verticesList[1], x, y);

			// Calculate barycentric weights
			float denom = f12 + f20 + f01;
			float alpha = f12 / denom;
			float beta = f20 / denom;
			float gamma = f01 / denom;

			if ((alpha >= 0 && beta >= 0 && gamma >= 0)) {
				// Calculate depth (z) using barycentric coordinates
				int z_inter = alpha * verticesList[0][2] + beta * verticesList[1][2] + gamma * verticesList[2][2];

				GzIntensity a, b, g, r;
				GzDepth z;
				GzGetDisplay(render->display, x, y, &r, &g, &b, &a, &z);
				if (z_inter <  z || z_inter == 0) {
					// Update Z-buffer with the new depth value
					GzPutDisplay(render->display, x, y, ctoi(render->flatcolor[0]),
						ctoi(render->flatcolor[1]), ctoi(render->flatcolor[2]), 1, z_inter);
				}
			}
		}
	}
	return GZ_SUCCESS;
}
