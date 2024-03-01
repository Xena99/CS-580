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


float Determinant3x3(GzMatrix matrix) {
	return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
		matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
		matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}


float Determinant(GzMatrix matrix) {
	float det = 0;
	for (int i = 0; i < 4; i++) {
		GzMatrix submatrix;
		for (int j = 1; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				if (k < i) {
					submatrix[j - 1][k] = matrix[j][k];
				}
				else if (k > i) {
					submatrix[j - 1][k - 1] = matrix[j][k];
				}
			}
		}
		det += (i % 2 == 0 ? 1 : -1) * matrix[0][i] * Determinant3x3(submatrix);
	}
	return det;
}

void Adjoint(GzMatrix matrix, GzMatrix adjoint) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			GzMatrix submatrix;
			for (int k = 0; k < 4; k++) {
				for (int l = 0; l < 4; l++) {
					if (k != i && l != j) {
						submatrix[k < i ? k : k - 1][l < j ? l : l - 1] = matrix[k][l];
					}
				}
			}
			adjoint[j][i] = ((i + j) % 2 == 0 ? 1 : -1) * Determinant3x3(submatrix);
		}
	}
}

int InverseAndTranspose(GzMatrix matrix, GzMatrix &result) {
	// Compute the determinant of the matrix
	float det = Determinant(matrix);
	if (det == 0) {
		return GZ_FAILURE; // Matrix is not invertible
	}

	// Compute the adjoint of the matrix
	GzMatrix adjoint;
	Adjoint(matrix, adjoint);

	// Compute the inverse by dividing the adjoint by the determinant
	GzMatrix inverse;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			inverse[i][j] = adjoint[i][j] / det;
		}
	}

	// Transpose the inverse to get the result
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			result[i][j] = inverse[j][i];
		}
	}

	return GZ_SUCCESS;
}

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

	return GZ_SUCCESS;
}

int GzPushMatrix(GzRender *render, GzMatrix matrix) {
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
	render->normal_level++;
	if (InverseAndTranspose(render->Ximage[render->matlevel], render->Xnorm[render->normal_level]) != GZ_SUCCESS) {
		return GZ_FAILURE;
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
	newRender->normal_level = -1;
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
	render->normal_level = -1;

	GzPushMatrix(render, render->Xsp);

	GzPushMatrix(render, Xpi);

	GzPushMatrix(render, Xiw);
	render->open = 1;

	return GZ_SUCCESS;
}

int GzPutAttribute(GzRender *render, int numAttributes, GzToken	*nameList, GzPointer *valueList)
{
	/* HW 2.1
	-- Set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
	-- In later homeworks set shaders, interpolaters, texture maps, and lights
	*/
	for (int position = 0; position < numAttributes; position++) {
		GzToken currentToken = nameList[position];
		if (currentToken == GZ_RGB_COLOR) {
			float* color = (float*)valueList[position];
			render->flatcolor[0] = color[0];
			render->flatcolor[1] = color[1];
			render->flatcolor[2] = color[2];
		}
		else if (currentToken == GZ_SPECULAR_COEFFICIENT) {
			float* specular = (float*)valueList[position];
			render->Ks[0] = specular[0];
			render->Ks[1] = specular[1];
			render->Ks[2] = specular[2];
		}
		else if (currentToken == GZ_DIFFUSE_COEFFICIENT) {
			float* diffuse = (float*)valueList[position];
			render->Kd[0] = diffuse[0];
			render->Kd[1] = diffuse[1];
			render->Kd[2] = diffuse[2];
		}
		else if (currentToken == GZ_AMBIENT_COEFFICIENT) {
			float* ambient = (float*)valueList[position];
			render->Ka[0] = ambient[0];
			render->Ka[1] = ambient[1];
			render->Ka[2] = ambient[2];
		}
		else if (currentToken == GZ_DISTRIBUTION_COEFFICIENT) {
			render->spec = *(float*)valueList[position];
		}
		else if (currentToken == GZ_DIRECTIONAL_LIGHT) {
			GzLight* dirLight = (GzLight*)valueList[position];
			for (int i = 0; i < 3; i++) {
				render->lights[position].direction[i] = dirLight->direction[i];
				render->lights[position].color[i] = dirLight->color[i];
			}
			render->numlights++;

		}
		else if (currentToken == GZ_AMBIENT_LIGHT) {
			GzLight* ambLight = (GzLight*)valueList[position];
			for (int i = 0; i < 3; i++) {
				render->ambientlight.direction[i] = ambLight->direction[i];
				render->ambientlight.color[i] = ambLight->color[i];
			}

		}
		else if (currentToken == GZ_INTERPOLATE) {
			render->interp_mode = *(int*)valueList[position];
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
int GzShadingEquation(GzRender *render, GzColor color, GzCoord norm) {
	GzCoord E = { 0, 0, -1 };
	Normalize(E);
	GzCoord normal = { norm[X], norm[Y], norm[Z] };
	float NdotE = DotProduct(normal, E);
	if (NdotE < 0) {
		// Flip normal for double-sided surfaces
		normal[X] = -normal[X];
		normal[Y] = -normal[Y];
		normal[Z] = -normal[Z];
		NdotE = -NdotE;
	}

	GzColor specularSum = { 0, 0, 0 };
	GzColor diffuseSum = { 0, 0, 0 };

	for (int i = 0; i < 3; ++i) {
		GzCoord light;
		light[X] = render->lights[i].direction[X];
		light[Y] = render->lights[i].direction[Y];
		light[Z] = render->lights[i].direction[Z];
		Normalize(light);

		float NdotL = DotProduct(normal, light);
		if (NdotL < 0) {
			// Skip light if it's on the wrong side
			continue;
		}

		// Reflect vector R
		GzCoord R = {
			2 * NdotL * norm[X] - light[X],
			2 * NdotL * norm[Y] - light[Y],
			2 * NdotL * norm[Z] - light[Z]
		};
		Normalize(R);

		float RdotE = max(0, DotProduct(R, E));

		// Specular component
		for (int j = 0; j < 3; ++j) {
			specularSum[j] += render->lights[i].color[j] * pow(RdotE, render->spec);
		}

		// Diffuse component
		for (int j = 0; j < 3; ++j) {
			diffuseSum[j] += render->lights[i].color[j] * NdotL;
		}
	}

	// Add in the ambient component
	GzColor ambientContribution = {
		render->ambientlight.color[0] * render->Ka[0],
		render->ambientlight.color[1] * render->Ka[1],
		render->ambientlight.color[2] * render->Ka[2]
	};

	// Combine the components and apply the material's coefficients
	for (int i = 0; i < 3; ++i) {
		color[i] = ambientContribution[i] +
			diffuseSum[i] * render->Kd[i] +
			specularSum[i] * render->Ks[i];
		color[i] = max(0, min(color[i], 1)); // Clamp color to [0, 1]
	}

	return GZ_SUCCESS;
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

	GzCoord transformedNormals[3];
	for (int i = 0; i < 3; i++) {
		// Define the normal vector in homogeneous coordinates
		float normal[4] = { normals[i][0], normals[i][1], normals[i][2], 0 };

		// Transform the normal vector using the combined normal transformation matrix
		float transformedNormal[4] = { 0 };
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				transformedNormal[j] += render->Xnorm[render->normal_level][j][k] * normal[k];
			}
		}

		// Normalize the transformed normal vector
		float magnitude = sqrt(transformedNormal[0] * transformedNormal[0] + transformedNormal[1] * transformedNormal[1] + transformedNormal[2] * transformedNormal[2]);
		transformedNormals[i][0] = transformedNormal[0] / magnitude;
		transformedNormals[i][1] = transformedNormal[1] / magnitude;
		transformedNormals[i][2] = transformedNormal[2] / magnitude;
	}

	// Sort vertices by Y-coordinate
	//sortVertices(transformedVertices);

	//Get Bounding box
	int minx = (int)(min(min(transformedVertices[0][0], transformedVertices[1][0]), transformedVertices[2][0]) + 0.5);
	int maxx = (int)(max(max(transformedVertices[0][0], transformedVertices[1][0]), transformedVertices[2][0]) + 0.5);
	int miny = (int)(min(min(transformedVertices[0][1], transformedVertices[1][1]), transformedVertices[2][1]) + 0.5);
	int maxy = (int)(max(max(transformedVertices[0][1], transformedVertices[1][1]), transformedVertices[2][1]) + 0.5);

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
				if (zInterpolated < z ) {
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
					
					GzPutDisplay(render->display, x, y, (GzIntensity)ctoi(color[RED]), (GzIntensity)ctoi(color[GREEN]), (GzIntensity)ctoi(color[BLUE]), 1, zInterpolated);
				}
			}
		}
	}

	return GZ_SUCCESS;
}

