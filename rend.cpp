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
#include <iostream>
#include <cmath>

#define GZ_SUCCESS 0
#define GZ_FAILURE 1
typedef float GzMatrix[4][4];

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
			// Calculate the submatrix for (i, j)
			GzMatrix submatrix;
			int subi = 0;
			for (int k = 0; k < 4; k++) {
				if (k == i) continue;
				int subj = 0;
				for (int l = 0; l < 4; l++) {
					if (l == j) continue;
					submatrix[subi][subj] = matrix[k][l];
					subj++;
				}
				subi++;
			}
			// Calculate the cofactor for (i, j)
			float cofactor = Determinant3x3(submatrix);
			if ((i + j) % 2 != 0) cofactor = -cofactor;
			adjoint[j][i] = cofactor; // Transpose the cofactor matrix
		}
	}
}

int InverseAndTranspose(GzMatrix matrix, GzMatrix &result) {
	// Compute the determinant of the matrix
	float det = Determinant(matrix);
	if (fabs(det) < 1e-10) {
		return GZ_FAILURE; // Matrix is not invertible
	}

	// Compute the adjoint of the matrix
	GzMatrix adjoint;
	Adjoint(matrix, adjoint);
	GzMatrix invMatrix;
	// Compute the inverse by dividing the adjoint by the determinant
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			invMatrix[i][j] = adjoint[i][j] / det;
		}
	}

	// Transpose only the top-left 3x3 part of the inverse matrix to get the result
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			result[i][j] = invMatrix[j][i];
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

void MatrixMultiply(GzMatrix a, GzMatrix b, GzMatrix& result) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			result[i][j] = 0;
			for (int k = 0; k < 4; k++) {
				result[i][j] += a[i][k] * b[k][j];
			}
		}
	}
}

int GzPushMatrix(GzRender *render, GzMatrix matrix, bool isProjectionMatrix) {
	if (render == nullptr) {
		return GZ_FAILURE;
	}

	if (render->matlevel >= (MATLEVELS - 1)) {
		return GZ_FAILURE;
	}

	// Push the matrix onto the image matrix stack
	render->matlevel++;
	GzMatrix tempMatrix;
	if (render->matlevel == 0) {
		memcpy(render->Ximage[render->matlevel], matrix, sizeof(GzMatrix));
	}
	else {
		MatrixMultiply(render->Ximage[render->matlevel - 1], matrix, render->Ximage[render->matlevel]);
	}

	if (!isProjectionMatrix) {
		render->normal_level++;
		if (render->normal_level == 0) {
			memcpy(render->Xnorm[render->normal_level], matrix, sizeof(GzMatrix));
		}
		else {
			MatrixMultiply(render->Xnorm[render->normal_level - 1], matrix, render->Xnorm[render->normal_level]);
		}

		if (InverseAndTranspose(render->Xnorm[render->normal_level], render->Xnorm[render->normal_level]) != GZ_SUCCESS) {
			return GZ_FAILURE;
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
	newRender->numlights = 0;
	newRender->dx = 0;
	newRender->dy = 0;

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

	// For the screen space transformation matrix
	GzPushMatrix(render, render->Xsp, true);

	// For the perspective projection matrix
	GzPushMatrix(render, Xpi, true);

	// For the world to camera (view) space transformation matrix
	GzPushMatrix(render, Xiw, false);
	render->open = 1;

	return GZ_SUCCESS;
}
int GzPutAttribute(GzRender *render, int numAttributes, GzToken *nameList, GzPointer *valueList)
{

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
		else if (nameList[i] == GZ_AMBIENT_COEFFICIENT) {
			GzColor *ka = (GzColor*)valueList[i];
			render->Ka[RED] = (*ka)[RED];
			render->Ka[GREEN] = (*ka)[GREEN];
			render->Ka[BLUE] = (*ka)[BLUE];

		}

		else if (nameList[i] == GZ_DIFFUSE_COEFFICIENT) {
			GzColor *kd = (GzColor*)valueList[i];
			render->Kd[RED] = (*kd)[RED];
			render->Kd[GREEN] = (*kd)[GREEN];
			render->Kd[BLUE] = (*kd)[BLUE];

		}

		else if (nameList[i] == GZ_SPECULAR_COEFFICIENT) {
			GzColor *ks = (GzColor*)valueList[i];
			render->Ks[RED] = (*ks)[RED];
			render->Ks[GREEN] = (*ks)[GREEN];
			render->Ks[BLUE] = (*ks)[BLUE];

		}

		else if (nameList[i] == GZ_AMBIENT_LIGHT) {
			GzLight *ambientlight = (GzLight*)valueList[i];
			render->ambientlight.color[RED] = ambientlight->color[RED];
			render->ambientlight.color[GREEN] = ambientlight->color[GREEN];
			render->ambientlight.color[BLUE] = ambientlight->color[BLUE];
		}

		else if (nameList[i] == GZ_DIRECTIONAL_LIGHT) {
			if (render->numlights >= MAX_LIGHTS) {
				return GZ_FAILURE;
			}
			GzLight* directionallight = (GzLight*)valueList[i];
			render->lights[render->numlights].color[RED] = directionallight->color[RED];
			render->lights[render->numlights].color[GREEN] = directionallight->color[GREEN];
			render->lights[render->numlights].color[BLUE] = directionallight->color[BLUE];
			render->lights[render->numlights].direction[RED] = directionallight->direction[RED];
			render->lights[render->numlights].direction[GREEN] = directionallight->direction[GREEN];
			render->lights[render->numlights].direction[BLUE] = directionallight->direction[BLUE];
			render->numlights++;
		}

		else if (nameList[i] == GZ_DISTRIBUTION_COEFFICIENT) {
			float *spec = (float*)valueList[i];
			render->spec = *spec;
		}
		else if (nameList[i] == GZ_INTERPOLATE) {
			int *interpolate_mode = (int*)valueList[i];
			render->interp_mode = *interpolate_mode;
		}

		else if (nameList[i] == GZ_TEXTURE_MAP) {
			GzTexture texture = (GzTexture)valueList[i];
			render->tex_fun = texture;
		}

		//hw6

		else if (nameList[i] == GZ_AASHIFTX) {
			float *dx = (float*)valueList[i];
			render->dx = *dx;
		}
		else if (nameList[i] == GZ_AASHIFTY) {
			float *dy = (float*)valueList[i];
			render->dy = *dy;
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
float angleBetween(GzCoord v1, GzCoord v2) {
	float dot = v1[0] * v2[0] + v1[1] * v2[1];
	float det = v1[0] * v2[1] - v1[1] * v2[0];
	return atan2(det, dot);
}

int GzShadingEquation(GzRender *render, GzColor& color, const GzCoord normal) {
	GzCoord E;

	// Compute the view direction vector
	E[X] = render->camera.lookat[X] - render->camera.position[X];
	E[Y] = render->camera.lookat[Y] - render->camera.position[Y];
	E[Z] = render->camera.lookat[Z] - render->camera.position[Z];

	// Normalize the view direction vector
	Normalize(E);

	GzCoord norm = { normal[X], normal[Y], normal[Z] };
	Normalize(norm); // Ensure normal vector is normalized

	GzColor specularSum = { 0, 0, 0 };
	GzColor diffuseSum = { 0, 0, 0 };

	for (int i = 0; i < render->numlights; ++i) {
		GzCoord light = { render->lights[i].direction[X], render->lights[i].direction[Y], render->lights[i].direction[Z] };
		Normalize(light);

		float NdotL = DotProduct(norm, light);
		if (NdotL > 0) {
			GzCoord R;
			for (int j = 0; j < 3; ++j) {
				R[j] = 2 * NdotL * norm[j] - light[j];
			}
			Normalize(R);

			float RdotE = max(0.0f, DotProduct(R, E));
			for (int j = 0; j < 3; ++j) {
				specularSum[j] += render->lights[i].color[j] * std::pow(RdotE, render->spec);
			}

			for (int j = 0; j < 3; ++j) {
				diffuseSum[j] += render->lights[i].color[j] * NdotL;
			}
		}
	}

	// Calculate ambient contribution
	GzColor ambientContribution = {
		render->ambientlight.color[RED] * render->Ka[RED],
		render->ambientlight.color[GREEN] * render->Ka[GREEN],
		render->ambientlight.color[BLUE] * render->Ka[BLUE]
	};

	for (int i = 0; i < 3; ++i) {
		color[i] = ambientContribution[i] +
			diffuseSum[i] * render->Kd[i] +
			specularSum[i] * render->Ks[i];
		color[i] = max(0.0f, min(color[i], 1.0f)); // Clamp the color values to the range [0, 1]
	}

	return GZ_SUCCESS;
}
void sortVerticesAndNormals(GzCoord* vertices, GzCoord* normals, GzTextureIndex* textureIndices) {
	// Find the centroid of the triangle
	GzCoord centroid = {
		(vertices[0][0] + vertices[1][0] + vertices[2][0]) / 3,
		(vertices[0][1] + vertices[1][1] + vertices[2][1]) / 3,
		(vertices[0][2] + vertices[1][2] + vertices[2][2]) / 3
	};

	// Calculate angles of each vertex relative to the centroid
	float angles[3];
	for (int i = 0; i < 3; ++i) {
		GzCoord vec = {
			vertices[i][0] - centroid[0],
			vertices[i][1] - centroid[1],
			0 // We only need X and Y for angle calculation
		};
		angles[i] = atan2(vec[1], vec[0]);
	}

	// Sort vertices based on angles in anti-clockwise order
	for (int i = 0; i < 2; ++i) {
		for (int j = i + 1; j < 3; ++j) {
			if (angles[i] > angles[j]) {
				// Swap vertices
				std::swap(vertices[i], vertices[j]);
				// Swap corresponding normals
				std::swap(normals[i], normals[j]);
				// Swap angles
				std::swap(angles[i], angles[j]);
				// Swap corresponding texture indices
				std::swap(textureIndices[i], textureIndices[j]);
			}
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
	GzTextureIndex* textureIndices = nullptr;
	GzTextureIndex transformedTextureIndices[3];
	float W[3];
	// Find texture coordinates in the valueList
	for (int i = 0; i < numParts; ++i) {
		if (nameList[i] == GZ_TEXTURE_INDEX) {
			textureIndices = static_cast<GzTextureIndex*>(valueList[i]);
			break;
		}
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
		if (textureIndices != nullptr) {
			W[i] = 1.0f / w; 
			transformedTextureIndices[i][U] = textureIndices[i][U] * W[i];
			transformedTextureIndices[i][V] = textureIndices[i][V] * W[i];
		}
	}

	GzCoord transformedNormals[3];
	for (int i = 0; i < 3; i++) {
		// Transform the normal vector using the inverse transpose of the model matrix (upper-left 3x3)
		GzCoord transformedNormal = { 0, 0, 0 };
		for (int j = 0; j < 3; j++) { // Iterate over the first three elements of the normal
			for (int k = 0; k < 3; k++) { // Use only the 3x3 submatrix for transformation
				transformedNormal[j] += render->Xnorm[render->normal_level][j][k] * normals[i][k];
			}
		}
		Normalize(transformedNormal); // Normalize the transformed normal
		for (int j = 0; j < 3; j++) {
			transformedNormals[i][j] = transformedNormal[j];
		}
	}

	// Apply anti-aliasing offsets
	for (int i = 0; i < 3; i++) {
		transformedVertices[i][X] += render->dx;
		transformedVertices[i][Y] += render->dy;
	}

	sortVerticesAndNormals(transformedVertices, transformedNormals, transformedTextureIndices);

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
				float interpolatedW = 1 / (alpha / W[0] + beta / W[1] + gamma / W[2]); // Interpolate 1/w for the pixel

				// Perform Z-buffer test
				GzIntensity r, g, b, a;
				GzDepth z;
				GzGetDisplay(render->display, x, y, &r, &g, &b, &a, &z);
				if (zInterpolated < z) {
					GzColor color;

					float interpolatedUW = alpha * transformedTextureIndices[0][U] + beta * transformedTextureIndices[1][U] + gamma * transformedTextureIndices[2][U];
					float interpolatedVW = alpha * transformedTextureIndices[0][V] + beta * transformedTextureIndices[1][V] + gamma * transformedTextureIndices[2][V];
					float interpolatedW = 1 / (alpha / W[0] + beta / W[1] + gamma / W[2]);

					float correctedU = interpolatedUW / interpolatedW; // Divide by interpolatedW to apply perspective correction
					float correctedV = interpolatedVW / interpolatedW; // Divide by interpolatedW to apply perspective correction

					// Look up texture color at this pixel
					GzColor texColor = { 0.0f, 0.0f, 0.0f }; // Default color if no texture function is provided
					if (render->tex_fun != NULL) {
						render->tex_fun(correctedU, correctedV, texColor);
					}
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

					color[RED] *= texColor[RED];
					color[GREEN] *= texColor[GREEN];
					color[BLUE] *= texColor[BLUE];

					GzPutDisplay(render->display, x, y, (GzIntensity)ctoi(color[RED]), (GzIntensity)ctoi(color[GREEN]), (GzIntensity)ctoi(color[BLUE]), 1, zInterpolated);
				}
			}
		}
	}

	return GZ_SUCCESS;
}
