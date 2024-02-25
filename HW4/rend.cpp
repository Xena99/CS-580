/* CS580 Homework 3 */

#include "stdafx.h"
#include "stdio.h"
#include "math.h"
#include "Gz.h"
#include "rend.h"
#include "disp.h"
#include <algorithm>
#include <iostream>

#define Pi 3.1415926

template <typename T>
T clamp(T value, T min, T max) {
	return (value < min) ? min : (value > max) ? max : value;
}

void ComputeInverseTranspose(const GzMatrix matrix, GzMatrix inverseTranspose);

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

int GzPushMatrix(GzRender *render, GzMatrix matrix, bool isNormalTransform = false) {
	if (render == nullptr) {
		return GZ_FAILURE;
	}

	if (render->matlevel >= (MATLEVELS - 1)) {
		return GZ_FAILURE;
	}

	render->matlevel++;
	if (render->matlevel == 0) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				render->Ximage[render->matlevel][i][j] = matrix[i][j];
			}
		}
	}
	else {
		GzMatrix tempMatrix;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				tempMatrix[i][j] = 0;
				for (int k = 0; k < 4; k++) {
					tempMatrix[i][j] += render->Ximage[render->matlevel - 1][i][k] * matrix[k][j];
				}
			}
		}
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				render->Ximage[render->matlevel][i][j] = tempMatrix[i][j];
			}
		}
	}

	if (isNormalTransform) {
		GzMatrix normalMatrix;
		if (matrix[0][1] == 0 && matrix[0][2] == 0 && matrix[1][0] == 0 &&
			matrix[1][2] == 0 && matrix[2][0] == 0 && matrix[2][1] == 0) {
			// Rotation matrix or uniform scaling matrix
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					normalMatrix[i][j] = matrix[i][j];
				}
				normalMatrix[i][3] = 0.0f; // No translation for normals
			}
		}
		else {
			// Non-uniform scaling matrix
			ComputeInverseTranspose(matrix, normalMatrix);
		}
		normalMatrix[3][0] = normalMatrix[3][1] = normalMatrix[3][2] = 0.0f;
		normalMatrix[3][3] = 1.0f;
		GzPushNormalMatrix(render, normalMatrix);
	}

	return GZ_SUCCESS;
}

int GzPushNormalMatrix(GzRender *render, GzMatrix matrix) {
	if (render == nullptr || render->matlevel >= MATLEVELS - 1) {
		return GZ_FAILURE;
	}

	if (render->matlevel == 0) {
		// If this is the first matrix on the stack, simply copy it
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				render->Xnorm[render->matlevel][i][j] = matrix[i][j];
			}
		}
	}
	else {
		// If there are already matrices on the stack, accumulate the new matrix with the top matrix
		GzMatrix tempMatrix;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				tempMatrix[i][j] = 0;
				for (int k = 0; k < 4; k++) {
					tempMatrix[i][j] += render->Xnorm[render->matlevel - 1][i][k] * matrix[k][j];
				}
			}
		}
		// Copy the accumulated matrix to the top of the stack
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				render->Xnorm[render->matlevel][i][j] = tempMatrix[i][j];
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
	GzMatrix identity =
	{
		1.0,    0.0,    0.0,    0.0,
		0.0,    1.0,    0.0,    0.0,
		0.0,    0.0,    1.0,    0.0,
		0.0,    0.0,    0.0,    1.0
	};
	// Setup Xsp - Screen space transformation
	GzPushMatrix(render, render->Xsp);
	GzPushNormalMatrix(render, identity);
	// No normal matrix needed for Xsp

	// Setup Xpi - Perspective transformation
	GzPushMatrix(render, Xpi);
	GzPushNormalMatrix(render, identity);

	// No normal matrix needed for Xpi

	// Setup Xiw - Camera transformation (modelview matrix)
	GzPushMatrix(render, Xiw);
	GzMatrix inverseTranspose;
	ComputeInverseTranspose(Xiw, inverseTranspose);
	GzPushNormalMatrix(render, inverseTranspose);

	

	render->open = 1;

	return GZ_SUCCESS;
}

int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, GzPointer *valueList) /* void** valuelist */
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
			GzColor* c = (GzColor*)valueList[i];
			render->flatcolor[0] = max(min(c[i][RED], 4095), 0);
			render->flatcolor[1] = max(min(c[i][GREEN], 4095), 0);
			render->flatcolor[2] = max(min(c[i][BLUE], 4095), 0);
		}
		if (nameList[i] == GZ_DIRECTIONAL_LIGHT)
		{
			GzLight* lightDir = (GzLight*)valueList[i];
			if (render->numlights < 0)
				render->numlights = 0;
			render->lights[render->numlights].direction[0] = lightDir->direction[0];
			render->lights[render->numlights].direction[1] = lightDir->direction[1];
			render->lights[render->numlights].direction[2] = lightDir->direction[2];
			render->lights[render->numlights].color[0] = lightDir->color[0];
			render->lights[render->numlights].color[1] = lightDir->color[1];
			render->lights[render->numlights].color[2] = lightDir->color[2];
			render->numlights++;
		}
		if (nameList[i] == GZ_AMBIENT_LIGHT)
		{
			GzLight* lightAmb = (GzLight*)valueList[i];
			render->ambientlight.direction[0] = lightAmb->direction[0];
			render->ambientlight.direction[1] = lightAmb->direction[1];
			render->ambientlight.direction[2] = lightAmb->direction[2];
			render->ambientlight.color[0] = lightAmb->color[0];
			render->ambientlight.color[1] = lightAmb->color[1];
			render->ambientlight.color[2] = lightAmb->color[2];
		}
		if (nameList[i] == GZ_AMBIENT_COEFFICIENT)
		{
			GzColor* colorAmb = (GzColor*)valueList[i];
			render->Ka[0] = (*colorAmb)[0];
			render->Ka[1] = (*colorAmb)[1];
			render->Ka[2] = (*colorAmb)[2];
		}
		if (nameList[i] == GZ_DIFFUSE_COEFFICIENT)
		{
			GzColor* colorDiff = (GzColor*)valueList[0];
			render->Kd[0] = (*colorDiff)[0];
			render->Kd[1] = (*colorDiff)[1];
			render->Kd[2] = (*colorDiff)[2];
		}
		if (nameList[i] == GZ_SPECULAR_COEFFICIENT)
		{
			GzColor* colorSpec = (GzColor*)valueList[i];
			render->Ks[0] = (*colorSpec)[0];
			render->Ks[1] = (*colorSpec)[1];
			render->Ks[2] = (*colorSpec)[2];
		}
		if(nameList[i] == GZ_INTERPOLATE) 
		{
			int* mode = (int*)valueList[i];
			render->interp_mode = *mode;
		}
		if (nameList[i] == GZ_DISTRIBUTION_COEFFICIENT) {
			render->spec = *(float*)valueList[i];
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

int GzShadingEquation(GzRender *render, GzColor color, GzCoord& norm) {
	Normalize(norm);
	
	GzCoord E = {
		0,0,-1
	};
	Normalize(E);

	float NdotE = DotProduct(norm, E);
	GzColor specSum = { 0, 0, 0 }, diffSum = { 0, 0, 0 };

	for (int i = 0; i < render->numlights; ++i) {
		float NdotL = DotProduct(norm, render->lights[i].direction);
		bool validLight = (NdotL >= 0 && NdotE >= 0) || (NdotL < 0 && NdotE < 0);

		if (validLight) {
			if (NdotL < 0) NdotL = (-1) * NdotL;

			GzCoord R = {
				2 * NdotL * norm[X] - render->lights[i].direction[X],
				2 * NdotL * norm[Y] - render->lights[i].direction[Y],
				2 * NdotL * norm[Z] - render->lights[i].direction[Z]
			};
			Normalize(R);

			float RdotE = DotProduct(R, E);
			if (RdotE < 0) RdotE = 0;
			if (RdotE > 1) RdotE = 1;
			for (int j = 0; j < 3; ++j) {
				specSum[j] += render->lights[i].color[j] * pow(RdotE, render->spec);
				diffSum[j] += render->lights[i].color[j] * NdotL;
			}
		}
	}

	for (int i = 0; i < 3; ++i) {
		color[i] = render->Ka[i] * render->ambientlight.color[i] +
			render->Kd[i] * diffSum[i] +
			render->Ks[i] * specSum[i];
	}

	return GZ_SUCCESS;
}

void ComputeInverseTranspose(const GzMatrix matrix, GzMatrix inverseTranspose) {
	// Initialize the last row and column of inverseTranspose matrix to identity
	for (int i = 0; i < 4; ++i) {
		inverseTranspose[i][3] = inverseTranspose[3][i] = 0.0f;
	}
	inverseTranspose[3][3] = 1.0f;

	// Compute the determinant of the 3x3 submatrix
	float det = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
		matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
		matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);

	// Check for non-invertible matrix
	if (fabs(det) < 1e-6) {
		// Matrix is not invertible, return identity matrix
		return;
	}

	// Compute the inverse of the 3x3 submatrix
	GzMatrix inverse;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			// Compute the cofactor of matrix[j][i]
			int i1 = (i + 1) % 3;
			int i2 = (i + 2) % 3;
			int j1 = (j + 1) % 3;
			int j2 = (j + 2) % 3;
			float cofactor = matrix[i1][j1] * matrix[i2][j2] - matrix[i1][j2] * matrix[i2][j1];

			// Apply the sign and divide by the determinant
			inverse[j][i] = cofactor / det;
		}
	}

	// Transpose the inverse of the 3x3 submatrix to get the inverse transpose
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			inverseTranspose[i][j] = inverse[j][i];
		}
	}
}


int GzPutTriangle(GzRender *render, int numParts, GzToken *nameList, GzPointer *valueList) {
	// Error checking
	if (render == nullptr || nameList == nullptr || valueList == nullptr) {
		return GZ_FAILURE;
	}

	GzCoord* vertices = static_cast<GzCoord*>(valueList[0]);
	GzCoord* normals = static_cast<GzCoord*>(valueList[1]);
	if (vertices == nullptr || normals == nullptr) {
		return GZ_FAILURE;
	}

	// Transform the vertices and normals
	GzCoord transformedVertices[3];
	GzCoord transformedNormals[3];
	for (int i = 0; i < 3; i++) {
		// Transform vertices
		GzMatrix& Ximage = render->Ximage[render->matlevel];
		float w = Ximage[3][0] * vertices[i][X] + Ximage[3][1] * vertices[i][Y] + Ximage[3][2] * vertices[i][Z] + Ximage[3][3];
		for (int j = 0; j < 3; j++) {
			transformedVertices[i][j] = (Ximage[j][0] * vertices[i][X] + Ximage[j][1] * vertices[i][Y] + Ximage[j][2] * vertices[i][Z] + Ximage[j][3]) / w;
		}

		// Transform normals
		GzMatrix& Xnorm = render->Xnorm[render->matlevel];
		for (int j = 0; j < 3; j++) {
			transformedNormals[i][j] = Xnorm[j][0] * normals[i][X] + Xnorm[j][1] * normals[i][Y] + Xnorm[j][2] * normals[i][Z];
		}
		Normalize(transformedNormals[i]);
	}
	
	// Sort vertices by Y-coordinate
	sortVertices(transformedVertices);

	// Get Bounding box
	int minx = floor(min(min(transformedVertices[0][0], transformedVertices[1][0]), transformedVertices[2][0]));
	int maxx = ceil(max(max(transformedVertices[0][0], transformedVertices[1][0]), transformedVertices[2][0]));
	int miny = floor(min(min(transformedVertices[0][1], transformedVertices[1][1]), transformedVertices[2][1]));
	int maxy = ceil(max(max(transformedVertices[0][1], transformedVertices[1][1]), transformedVertices[2][1]));

	// Rasterize the triangle
	for (int y = miny; y <= maxy; y++) {
		for (int x = minx; x <= maxx; x++) {
			float f12 = EdgeFunctionCW(transformedVertices[1], transformedVertices[2], x, y);
			float f20 = EdgeFunctionCW(transformedVertices[2], transformedVertices[0], x, y);
			float f01 = EdgeFunctionCW(transformedVertices[0], transformedVertices[1], x, y);

			// Calculate barycentric weights
			float denom = f12 + f20 + f01;
			float alpha = f12 / denom;
			float beta = f20 / denom;
			float gamma = f01 / denom;

			if (alpha >= 0 && beta >= 0 && gamma >= 0) {
				// Interpolate Z-value
				float zInterpolated = alpha * transformedVertices[0][Z] + beta * transformedVertices[1][Z] + gamma * transformedVertices[2][Z];

				// Perform Z-buffer test
				GzIntensity r, g, b, a;
				GzDepth z;
				GzGetDisplay(render->display, x, y, &r, &g, &b, &a, &z);
				if (zInterpolated < z || zInterpolated == INT_MAX) {
					GzColor color;
					if (render->interp_mode == GZ_NORMALS) {
						GzCoord interpolatedNormal = {
							alpha * transformedNormals[0][X] + beta * transformedNormals[1][X] + gamma * transformedNormals[2][X],
							alpha * transformedNormals[0][Y] + beta * transformedNormals[1][Y] + gamma * transformedNormals[2][Y],
							alpha * transformedNormals[0][Z] + beta * transformedNormals[1][Z] + gamma * transformedNormals[2][Z]
						};
						Normalize(interpolatedNormal);

						GzShadingEquation(render, color, interpolatedNormal);
					}
					else if (render->interp_mode == GZ_COLOR) {
						GzColor Color0, Color1, Color2;
						GzShadingEquation(render, Color0, transformedNormals[0]);
						GzShadingEquation(render, Color1, transformedNormals[1]);
						GzShadingEquation(render, Color2, transformedNormals[2]);

						color[0] = alpha * Color0[0] + beta * Color1[0] + gamma * Color2[0]; // Red
						color[1] = alpha * Color0[1] + beta * Color1[1] + gamma * Color2[1]; // Green
						color[2] = alpha * Color0[2] + beta * Color1[2] + gamma * Color2[2]; // Blue
					}
					// Update the pixel value in the display
					GzPutDisplay(render->display, x, y, ctoi(color[RED]), ctoi(color[GREEN]), ctoi(color[BLUE]), 1, zInterpolated);
				}
			}
		}
	}

	return GZ_SUCCESS;
}


