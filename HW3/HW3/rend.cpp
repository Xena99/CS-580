/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include "disp.h"
#include <algorithm>
#include <iostream>

#define Pi 3.1415926

int GzRotXMat(float degree, GzMatrix mat) {
	if (mat == nullptr) return GZ_FAILURE;

	float Rad = degree * Pi / 180;
	// Initialize to identity matrix
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			mat[i][j] = (i == j) ? 1.0f : 0.0f;
		}
	}
	// Set rotation values
	mat[1][1] = cos(Rad);
	mat[1][2] = -1.0f * sin(Rad);
	mat[2][1] = sin(Rad);
	mat[2][2] = cos(Rad);

	return GZ_SUCCESS;
}


int GzRotYMat(float degree, GzMatrix mat)
{
	if (mat == nullptr)
		return GZ_FAILURE;

	// Initialize to identity matrix
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			mat[i][j] = (i == j) ? 1.0f : 0.0f;
		}
	}
	float Rad = degree * Pi / 180;

	mat[0][0] = cos(Rad);
	mat[0][2] = sin(Rad);
	mat[2][0] = -1.0f * sin(Rad);
	mat[2][2] = cos(Rad);

	return GZ_SUCCESS;
}


int GzRotZMat(float degree, GzMatrix mat)
{
	if (mat == nullptr) return GZ_FAILURE;

	// Initialize to identity matrix
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			mat[i][j] = (i == j) ? 1.0f : 0.0f;
		}
	}

	float Rad = degree * Pi / 180;

	mat[0][0] = cos(Rad);
	mat[0][1] = -1.0f * sin(Rad);
	mat[1][0] = sin(Rad);
	mat[1][1] = cos(Rad);

	return GZ_SUCCESS;
}


int GzTrxMat(GzCoord translate, GzMatrix mat)
{
	if (mat == nullptr) return GZ_FAILURE;

	// Initialize to identity matrix
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			mat[i][j] = (i == j) ? 1.0f : 0.0f;
		}
	}

	mat[0][3] = translate[0];
	mat[1][3] = translate[1];
	mat[2][3] = translate[2];
	return GZ_SUCCESS;
}


int GzScaleMat(GzCoord scale, GzMatrix mat)
{
	if (mat == nullptr) return GZ_FAILURE;

	// Initialize to identity matrix
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			mat[i][j] = (i == j) ? 1.0f : 0.0f;
		}
	}

	mat[0][0] = scale[0];
	mat[1][1] = scale[1];
	mat[2][2] = scale[2];

	return GZ_SUCCESS;
}


int GzPutCamera(GzRender *render, GzCamera *camera) {
	if (render == nullptr || camera == nullptr) {
		// Return failure if either pointer is null
		return GZ_FAILURE;
	}

	// Overwrite the renderer's camera with the new camera definition

	render->camera.position[X] = camera->position[X];
	render->camera.position[Y] = camera->position[Y];
	render->camera.position[Z] = camera->position[Z];
	render->camera.lookat[X] = camera->lookat[X];
	render->camera.lookat[Y] = camera->lookat[Y];
	render->camera.lookat[Z] = camera->lookat[Z];
	render->camera.worldup[X] = camera->worldup[X];
	render->camera.worldup[Y] = camera->worldup[Y];
	render->camera.worldup[Z] = camera->worldup[Z];
	render->camera.FOV = camera->FOV;
	float d;
	d = tan(((render->camera.FOV) / 2) * (3.14159265 / 180));
	render->Xsp[2][2] = INT_MAX * d;

	// Return success
	return GZ_SUCCESS;
}

int GzPushMatrix(GzRender *render, GzMatrix matrix) {
	if (render == nullptr) {
		// Invalid render pointer
		return GZ_FAILURE;
	}

	if (render->matlevel >= (MATLEVELS - 1)) {
		// Stack overflow check
		return GZ_FAILURE;
	}

	// Increment the matrix level to point to a new spot for the incoming matrix
	render->matlevel++;

	// If the stack was empty, directly copy the incoming matrix.
	// Otherwise, multiply the top stack matrix with the incoming matrix and store the result.
	if (render->matlevel == 0) {
		// Directly copying the matrix for the first entry
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				render->Ximage[render->matlevel][i][j] = matrix[i][j];
			}
		}
	}
	else {
		// Multiply the top of the stack with the incoming matrix and store the result
		GzMatrix tempMatrix;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				tempMatrix[i][j] = 0; // Initialize elements of the temp matrix
				for (int k = 0; k < 4; k++) {
					// Perform matrix multiplication
					tempMatrix[i][j] += render->Ximage[render->matlevel - 1][i][k] * matrix[k][j];
				}
			}
		}
		// Copy the multiplication result back to the current top of the stack
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				render->Ximage[render->matlevel][i][j] = tempMatrix[i][j];
			}
		}
	}

	return GZ_SUCCESS;
}

int GzPopMatrix(GzRender *render) {
	if (render == nullptr) {
		// Invalid render pointer
		return GZ_FAILURE;
	}

	if (render->matlevel < 0) {
		// Stack underflow check
		return GZ_FAILURE;
	}

	// Decrement the matrix level to pop the top matrix off the stack
	render->matlevel--;

	return GZ_SUCCESS;
}


void Normalize(GzCoord& vec) {
	float magnitude = sqrt(vec[X] * vec[X] + vec[Y] * vec[Y] + vec[Z] * vec[Z]);
	if (magnitude > 0.0f) {
		vec[X] /= magnitude;
		vec[Y] /= magnitude;
		vec[Z] /= magnitude;
	}
}

float DotProduct(const GzCoord& vec1, const GzCoord& vec2) {
	return vec1[X] * vec2[X] + vec1[Y] * vec2[Y] + vec1[Z] * vec2[Z];
}

void CrossProduct(const GzCoord& vec1, const GzCoord& vec2, GzCoord& outVec) {
	outVec[X] = vec1[Y] * vec2[Z] - vec1[Z] * vec2[Y];
	outVec[Y] = vec1[Z] * vec2[X] - vec1[X] * vec2[Z];
	outVec[Z] = vec1[X] * vec2[Y] - vec1[Y] * vec2[X];
}

void ComputeXiw(GzMatrix& Xiw, const GzCamera& camera) {
	// Calculate the Z axis of camera space
	GzCoord camZ = {
		camera.lookat[X] - camera.position[X],
		camera.lookat[Y] - camera.position[Y],
		camera.lookat[Z] - camera.position[Z]
	};
	Normalize(camZ);

	// Up' = Up - (Up . Z) * Z
	GzCoord upPrime = {
		camera.worldup[X],
		camera.worldup[Y],
		camera.worldup[Z]
	};
	float upZdot = DotProduct(upPrime, camZ);
	for (int i = 0; i < 3; i++) {
		upPrime[i] -= upZdot * camZ[i];
	}
	Normalize(upPrime);

	// The X axis of camera space
	GzCoord camX;
	CrossProduct(upPrime, camZ, camX);
	Normalize(camX);

	// The Y axis is the cross product of Z axis and X axis
	GzCoord camY;
	CrossProduct(camZ, camX, camY);

	// Fill in the rotation part of Xiw
	for (int i = 0; i < 3; i++) {
		Xiw[0][i] = camX[i]; // X axis
		Xiw[1][i] = camY[i]; // Y axis
		Xiw[2][i] = camZ[i]; // Z axis
		Xiw[3][i] = 0.0f;    // Bottom row of the matrix
	}

	// Compute and fill in the translation part of Xiw
	Xiw[0][3] = -DotProduct(camX, camera.position); // Translation along X
	Xiw[1][3] = -DotProduct(camY, camera.position); // Translation along Y
	Xiw[2][3] = -DotProduct(camZ, camera.position); // Translation along Z
	Xiw[3][3] = 1.0f; // Homogeneous coordinate
}

//------------------------------------------------------------------------------------------------

int GzNewRender(GzRender **render, GzRenderClass renderClass, GzDisplay *display) {
	// Check for a valid render class
	if (renderClass != GZ_Z_BUFFER_RENDER) {
		return GZ_FAILURE;
	}

	// Validate input pointers
	if (render == NULL || display == NULL) {
		return GZ_FAILURE;
	}

	// Allocate memory for the new renderer
	GzRender *newRender = (GzRender*)malloc(sizeof(GzRender));
	if (newRender == NULL) {
		// If allocation failed, return failure
		return GZ_FAILURE;
	}

	// Initialize the default camera settings
	newRender->camera.FOV = DEFAULT_FOV;
	newRender->camera.position[X] = DEFAULT_IM_X;
	newRender->camera.position[Y] = DEFAULT_IM_Y;
	newRender->camera.position[Z] = DEFAULT_IM_Z;

	newRender->camera.lookat[X] = 0.0;
	newRender->camera.lookat[Y] = 0.0;
	newRender->camera.lookat[Z] = 0.0;

	newRender->camera.worldup[X] = 0.0;
	newRender->camera.worldup[Y] = 1.0;
	newRender->camera.worldup[Z] = 0.0;

	// The renderer starts closed (not ready for rendering)
	newRender->open = 0;
	// The matrix level is initialized to -1 (no matrices on the stack)
	newRender->matlevel = -1;
	// Assign the display passed to the function
	newRender->display = display;

	// Set up the viewport transformation matrix (Xsp)
	float d = INT_MAX * tan((newRender->camera.FOV / 2) * (Pi / 180));
	newRender->Xsp[0][0] = display->xres / 2.0f;
	newRender->Xsp[0][3] = display->xres / 2.0f;
	newRender->Xsp[1][1] = -display->yres / 2.0f;
	newRender->Xsp[1][3] = display->yres / 2.0f;
	newRender->Xsp[2][2] = d;
	newRender->Xsp[3][3] = 1.0f;
	// All other elements of Xsp should be zero
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			if (i != j && !(i == 2 && j == 2)) { // Skip the diagonal and the special case for Xsp[2][2]
				newRender->Xsp[i][j] = 0.0f;
			}
		}
	}

	// Set the newly created renderer to the output pointer
	*render = newRender;
	return GZ_SUCCESS;
}


int GzFreeRender(GzRender *render)
{
	/*
	-free all renderer resources
	*/
	if (render == 0)
		return GZ_FAILURE;

	GzFreeDisplay(render->display);
	delete[]render->Ximage;
	delete[]render->Xnorm;
	delete[]render->Xsp;
	delete[]render->flatcolor;

	//delete camera
	delete[]render->camera.Xiw;
	delete[]render->camera.Xpi;
	delete[]render->camera.position;
	delete[]render->camera.lookat;
	delete[]render->camera.worldup;

	//delete render
	delete render;
	return GZ_SUCCESS;
}

int GzBeginRender(GzRender *render) {
	if (render == NULL || render->display == NULL) {
		return GZ_FAILURE;
	}

	// Compute Xiw from camera definition
	GzMatrix Xiw;
	ComputeXiw(Xiw, render->camera);

	// Compute Xpi (Projection Transform)
	float d = tan((render->camera.FOV / 2) * (Pi / 180));
	GzMatrix Xpi = {
		{ 1.0, 0.0, 0.0, 0.0 },
		{ 0.0, 1.0, 0.0, 0.0 },
		{ 0.0, 0.0, 1.0, 0.0 },
		{ 0.0, 0.0, d,   1.0 }
	};

	// Initialize Ximage stack
	render->matlevel = -1; // Initialize stack pointer

	GzPushMatrix(render, render->Xsp);
	GzPushMatrix(render, Xpi);
	GzPushMatrix(render, Xiw);
	
	render->open = 1;

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
	render->open = 1;
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

int GzPutTriangle(GzRender *render, int numParts, GzToken *nameList, GzPointer *valueList) {
	if (render == nullptr || nameList == nullptr || valueList == nullptr) {
		return GZ_FAILURE;
	}

	GzCoord* vertices = static_cast<GzCoord*>(valueList[0]);
	if (vertices == nullptr) {
		return GZ_FAILURE;
	}

	// Transform the vertices using the current transformation matrix on the stack
	GzCoord transformedVertices[3];
	for (int i = 0; i < 3; i++) {
		// Homogeneous coordinates
		GzCoord& vertex = vertices[i];
		float w = render->Ximage[render->matlevel][3][0] * vertex[X] +
			render->Ximage[render->matlevel][3][1] * vertex[Y] +
			render->Ximage[render->matlevel][3][2] * vertex[Z] +
			render->Ximage[render->matlevel][3][3];

		for (int j = 0; j < 3; j++) {
			transformedVertices[i][j] = (render->Ximage[render->matlevel][j][0] * vertex[X] +
				render->Ximage[render->matlevel][j][1] * vertex[Y] +
				render->Ximage[render->matlevel][j][2] * vertex[Z] +
				render->Ximage[render->matlevel][j][3]) / w;
		}
	}

	// Clip triangles behind the view plane or off-screen
	// Note: For full clipping, you'd also check against the near and far clipping planes, and clip against the screen edges
	for (int i = 0; i < 3; i++) {
		if (transformedVertices[i][Z] < 0 ||
			transformedVertices[i][X] < 0 || transformedVertices[i][X] > render->display->xres ||
			transformedVertices[i][Y] < 0 || transformedVertices[i][Y] > render->display->yres) {
			return GZ_SUCCESS; // Discard the triangle
		}
	}

	// Sort vertices by Y-coordinate in ascending order
	sortVertices(transformedVertices);

	//Get Bounding box
	int minx = (int)(min(min(transformedVertices[0][0], transformedVertices[1][0]), transformedVertices[2][0]) + 0.5);
	int maxx = (int)(max(max(transformedVertices[0][0], transformedVertices[1][0]), transformedVertices[2][0]) + 0.5);
	int miny = (int)(min(min(transformedVertices[0][1], transformedVertices[1][1]), transformedVertices[2][1]) + 0.5);
	int maxy = (int)(max(max(transformedVertices[0][1], transformedVertices[1][1]), transformedVertices[2][1]) + 0.5);
	
	// Rasterize the triangle
	for (int y = miny; y <= maxy; y++) {
		for (int x = minx; x <= maxx; x++) {
			// Barycentric coordinates for clockwise order
			float f12 = EdgeFunctionCW(transformedVertices[1], transformedVertices[2], x, y);
			float f20 = EdgeFunctionCW(transformedVertices[2], transformedVertices[0], x, y);
			float f01 = EdgeFunctionCW(transformedVertices[0], transformedVertices[1], x, y);

			// Calculate barycentric weights
			float denom = f12 + f20 + f01;
			float alpha = f12 / denom;
			float beta = f20 / denom;
			float gamma = f01 / denom;

			if ((alpha >= 0 && beta >= 0 && gamma >= 0)) {
				// Calculate depth (z) using barycentric coordinates
				int z_inter = alpha * transformedVertices[0][2] + beta * transformedVertices[1][2] 
					+ gamma * transformedVertices[2][2];

				GzIntensity a, b, g, r;
				GzDepth z;
				GzGetDisplay(render->display, x, y, &r, &g, &b, &a, &z);
				if (z_inter < z || z_inter == 0) {
					// Update Z-buffer with the new depth value
					GzPutDisplay(render->display, x, y, ctoi(render->flatcolor[0]),
						ctoi(render->flatcolor[1]), ctoi(render->flatcolor[2]), 1, z_inter);
				}
				
			}
		}
	}

	return GZ_SUCCESS;
}

void Matrixcopy(GzMatrix m1, GzMatrix m2)
{

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			m1[i][j] = m2[i][j];

		}
	}
}

int NormalMatrix(GzRender *render, GzMatrix matrix)
{
	float temp = sqrt(pow(matrix[0][0], 2) + pow(matrix[0][1], 2) + pow(matrix[0][2], 2) + pow(matrix[0][3], 2));
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			matrix[i][j] /= temp;
		}
	}
	Matrixcopy(render->Xnorm[render->matlevel], matrix);
	return GZ_SUCCESS;

}


int Shader(GzRender *render, GzColor color, GzCoord normal) {
	GzColor specular = { 0, 0, 0 };
	GzColor diffuse = { 0, 0, 0 };
	GzColor ambient = { 0, 0, 0 };

	// Normalize the normal vector
	float norm = sqrt(normal[X] * normal[X] + normal[Y] * normal[Y] + normal[Z] * normal[Z]);
	normal[X] /= norm;
	normal[Y] /= norm;
	normal[Z] /= norm;

	// Calculate the ambient component
	ambient[RED] = render->ambientlight.color[RED] * render->Ka[RED];
	ambient[GREEN] = render->ambientlight.color[GREEN] * render->Ka[GREEN];
	ambient[BLUE] = render->ambientlight.color[BLUE] * render->Ka[BLUE];

	// Iterate through all lights to calculate the diffuse and specular components
	for (int i = 0; i < render->numlights; i++) {
		GzCoord lightDirection;
		float dotProduct;

		// Normalize the light direction
		norm = sqrt(render->lights[i].direction[X] * render->lights[i].direction[X] +
			render->lights[i].direction[Y] * render->lights[i].direction[Y] +
			render->lights[i].direction[Z] * render->lights[i].direction[Z]);
		lightDirection[X] = render->lights[i].direction[X] / norm;
		lightDirection[Y] = render->lights[i].direction[Y] / norm;
		lightDirection[Z] = render->lights[i].direction[Z] / norm;

		// Calculate the dot product of the normal and the light direction
		dotProduct = normal[X] * lightDirection[X] +
			normal[Y] * lightDirection[Y] +
			normal[Z] * lightDirection[Z];

		// Calculate the diffuse component
		if (dotProduct > 0) {
			diffuse[RED] += render->lights[i].color[RED] * render->Kd[RED] * dotProduct;
			diffuse[GREEN] += render->lights[i].color[GREEN] * render->Kd[GREEN] * dotProduct;
			diffuse[BLUE] += render->lights[i].color[BLUE] * render->Kd[BLUE] * dotProduct;

			// Calculate the specular component
			GzCoord reflection;
			reflection[X] = 2 * dotProduct * normal[X] - lightDirection[X];
			reflection[Y] = 2 * dotProduct * normal[Y] - lightDirection[Y];
			reflection[Z] = 2 * dotProduct * normal[Z] - lightDirection[Z];

			// Normalize the reflection vector
			norm = sqrt(reflection[X] * reflection[X] + reflection[Y] * reflection[Y] + reflection[Z] * reflection[Z]);
			reflection[X] /= norm;
			reflection[Y] /= norm;
			reflection[Z] /= norm;

			float viewDot = -reflection[Z]; // Assuming the view direction is along the Z-axis

			if (viewDot > 0) {
				float specPower = pow(viewDot, render->spec);
				specular[RED] += render->ambientlight[i].color[RED] * render->Ks[RED] * specPower;
				specular[GREEN] += render->lights[i].color[GREEN] * render->Ks[GREEN] * specPower;
				specular[BLUE] += render->lights[i].color[BLUE] * render->Ks[BLUE] * specPower;
			}
		}
	}

	// Combine the components to get the final color
	color[RED] = specular[RED] + diffuse[RED] + ambient[RED];
	color[GREEN] = specular[GREEN] + diffuse[GREEN] + ambient[GREEN];
	color[BLUE] = specular[BLUE] + diffuse[BLUE] + ambient[BLUE];

	// Clamp the color values to the range [0, 1]
	for (int i = 0; i < 3; i++) {
		if (color[i] < 0) color[i] = 0;
		if (color[i] > 4095) color[i] = 4095;
	}
	GZ_SUCCESS;
}