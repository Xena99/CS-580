// Application5.cpp: implementation of the Application5 class.
//
//////////////////////////////////////////////////////////////////////

/*
 * application test code for homework assignment #5
*/

#include "stdafx.h"
#include "CS580HW.h"
#include "Application5.h"
#include "Gz.h"
#include "disp.h"
#include "rend.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

#define INFILE  "teapot.asc"
#define OUTFILE "output.ppm"


extern int tex_fun(float u, float v, GzColor color); /* image texture function */
extern int ptex_fun(float u, float v, GzColor color); /* procedural texture function */
GzToken nameList1[2];
GzPointer valueList1[2];
float AAFilter[AAKERNEL_SIZE][3] 			/* X, Y, coef */
{
	{ -0.52, 0.38, 0.128 },
	{ 0.41, 0.56, 0.119 },
	{ 0.27, 0.08, 0.294 },
	{ -0.17, -0.29, 0.249 },
	{ 0.58, -0.55, 0.104 },
	{ -0.31, -0.71, 0.106 }
};


void shade(GzCoord norm, GzCoord color);

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Application5::Application5()
{

}

Application5::~Application5()
{
	
}

int Application5::Initialize()
{
 /* to be filled in by the app if it sets camera params */
	GzCamera	camera;
	int		    xRes, yRes, dispClass;	/* display parameters */

	GzToken		nameListShader[9]; 	    /* shader attribute names */
	GzPointer   valueListShader[9];		/* shader attribute pointers */
	GzToken     nameListLights[10];		/* light info */
	GzPointer   valueListLights[10];
	int			shaderType, interpStyle;
	float		specpower;
	int		status;

	status = 0;

	/*
	* Allocate memory for user input
	*/
	m_pUserInput = new GzInput;

	/*
	* initialize the display and the renderer
	*/
	m_nWidth = 800;		// frame buffer and display width
	m_nHeight = 800;    // frame buffer and display height

	for (int i = 0; i < 6; i++) {
		status |= GzNewFrameBuffer(&m_pFrameBuffer, m_nWidth, m_nHeight);

		status |= GzNewDisplay(&m_pDisplay[i], GZ_RGBAZ_DISPLAY, m_nWidth, m_nHeight);

		status |= GzGetDisplayParams(m_pDisplay[i], &xRes, &yRes, &dispClass);

		status |= GzInitDisplay(m_pDisplay[i]);

		status |= GzNewRender(&m_pRender[i], GZ_Z_BUFFER_RENDER, m_pDisplay[i]);

		/* Translation matrix */
		GzMatrix	scale =
		{
			30.0,	0.0,	0.0,	60.5,
			0.0,	30.0,	0.0,	-70.0,
			0.0,	0.0,	30.0,	90,
			0.0,	0.0,	0.0,	1
		};

		GzMatrix	rotateX =
		{
			1.0,	0.0,	0.0,	0.0,
			0.0,	.7071,	.7071,	0.0,
			0.0,	-.7071,	.7071,	0.0,
			0.0,	0.0,	0.0,	1.0
		};

		GzMatrix	rotateY =
		{
			.866,	0.0,	-0.5,	0.0,
			0.0,	1.0,	0.0,	0.0,
			0.5,	0.0,	.866,	0.0,
			0.0,	0.0,	0.0,	1.0
		};

		camera.position[X] = 13.2;
		camera.position[Y] = -8.7;
		camera.position[Z] = -14.8;

		camera.lookat[X] = 0.8;
		camera.lookat[Y] = 0.7;
		camera.lookat[Z] = 4.5;

		camera.worldup[X] = -0.2;
		camera.worldup[Y] = 1.0;
		camera.worldup[Z] = 0.0;

		camera.FOV = 53.7;              /* degrees */

		float pullBackFactor = 210.5; // Increase this number to pull the camera further back

									// Calculate the current view direction vector
		GzCoord viewDirection;
		viewDirection[X] = camera.lookat[X] - camera.position[X];
		viewDirection[Y] = camera.lookat[Y] - camera.position[Y];
		viewDirection[Z] = camera.lookat[Z] - camera.position[Z];

		// Normalize the view direction vector
		float length = sqrt(viewDirection[X] * viewDirection[X] +
			viewDirection[Y] * viewDirection[Y] +
			viewDirection[Z] * viewDirection[Z]);
		viewDirection[X] /= length;
		viewDirection[Y] /= length;
		viewDirection[Z] /= length;

		// Pull the camera back
		camera.position[X] -= viewDirection[X] * pullBackFactor;
		camera.position[Y] -= viewDirection[Y] * pullBackFactor;
		camera.position[Z] -= viewDirection[Z] * pullBackFactor;

		status |= GzPutCamera(m_pRender[i], &camera);

		/* Start Renderer */
		status |= GzBeginRender(m_pRender[i]);

		/* Light */
		//GzLight	light1 = { { -0.7071, 0.7071, 0 },{ 0.6, 0.6, 0.6 } }; // Increase light1 intensity
		//GzLight	light2 = { { 0, -0.7071, -0.7071 },{ 0.8, 0.8, 0.6 } }; // Increase light2 intensity
		//GzLight	light3 = { { 0.7071, 0.0, -0.7071 },{ 0.8, 0.4, 0.6 } }; // Increase light3 intensity
		//GzLight	ambientlight = { { 0, 0, 0 },{ 0.5, 0.5, 0.5 } }; // Increase ambient light intensity

		//														  /* Material property */
		//GzColor specularCoefficient = { 0.8, 0.8, 0.8 }; // Increase specular reflection
		//GzColor ambientCoefficient = { 0.6, 0.6, 0.6 }; // Increase ambient reflection
		//GzColor diffuseCoefficient = { 1.2, 1.2, 1.2 }; // Increase diffuse reflection

		/* Light */
		GzLight	light1 = { { -0.7071, 0.7071, 0 },{ 0.5, 0.4, 0.6 } };
		GzLight	light2 = { { 0, -0.7071, -0.7071 },{ 0.4, 0.4, 0.3 } };
		GzLight	light3 = { { 0.7071, 0.0, -0.7071 },{ 0.4, 0.35, 0.5 } };
		GzLight	ambientlight = { { 0, 0, 0 },{ 0.7, 0.7, 0.7 } };

		/* Material property */
		GzColor specularCoefficient = { 0.4, 0.4, 0.4 };
		GzColor ambientCoefficient = { 0.3, 0.3, 0.3 };
		GzColor diffuseCoefficient = { 0.9, 0.9, 0.9 };
	
		nameListLights[0] = GZ_DIRECTIONAL_LIGHT;
		valueListLights[0] = (GzPointer)&light1;
		nameListLights[1] = GZ_DIRECTIONAL_LIGHT;
		valueListLights[1] = (GzPointer)&light2;
		nameListLights[2] = GZ_DIRECTIONAL_LIGHT;
		valueListLights[2] = (GzPointer)&light3;
		status |= GzPutAttribute(m_pRender[i], 3, nameListLights, valueListLights);

		nameListLights[0] = GZ_AMBIENT_LIGHT;
		valueListLights[0] = (GzPointer)&ambientlight;
		status |= GzPutAttribute(m_pRender[i], 1, nameListLights, valueListLights);
		/*
		* Tokens associated with shading
		*/
		nameListShader[0] = GZ_DIFFUSE_COEFFICIENT;
		valueListShader[0] = (GzPointer)diffuseCoefficient;

		/*
		* Select either GZ_COLOR or GZ_NORMALS as interpolation mode
		*/
		nameListShader[1] = GZ_INTERPOLATE;
		interpStyle = GZ_NORMALS;         /* Phong shading */
		valueListShader[1] = (GzPointer)&interpStyle;

		nameListShader[2] = GZ_AMBIENT_COEFFICIENT;
		valueListShader[2] = (GzPointer)ambientCoefficient;
		nameListShader[3] = GZ_SPECULAR_COEFFICIENT;
		valueListShader[3] = (GzPointer)specularCoefficient;
		nameListShader[4] = GZ_DISTRIBUTION_COEFFICIENT;
		specpower = 32;
		valueListShader[4] = (GzPointer)&specpower;

		nameListShader[5] = GZ_TEXTURE_MAP;
		valueListShader[5] = (GzPointer)(tex_fun);	/* or use ptex_fun */
		//	valueListShader[5] = (GzPointer)(ptex_fun);
		//#endif
		status |= GzPutAttribute(m_pRender[i], 6, nameListShader, valueListShader);

		nameList1[0] = GZ_AASHIFTX;
		valueList1[0] = (GzPointer)&AAFilter[i][0];

		nameList1[1] = GZ_AASHIFTY;
		valueList1[1] = (GzPointer)&AAFilter[i][1];

		status |= GzPutAttribute(m_pRender[i], 2, nameList1, valueList1);

		status |= GzPushMatrix(m_pRender[i], scale, false);
		status |= GzPushMatrix(m_pRender[i], rotateY, false);
		status |= GzPushMatrix(m_pRender[i], rotateX, false);
	}

	if (status) exit(GZ_FAILURE);

	if (status)
		return(GZ_FAILURE);
	else
		return(GZ_SUCCESS);
}

int Application5::Render()
{
	GzToken		nameListTriangle[3]; 	/* vertex attribute names */
	GzPointer	valueListTriangle[3]; 	/* vertex attribute pointers */
	GzCoord		vertexList[3];	/* vertex position coordinates */
	GzCoord		normalList[3];	/* vertex normals */
	GzTextureIndex  	uvList[3];		/* vertex texture map indices */
	char		dummy[256];
	int			status;

	FILE *infile;
	FILE *outfile;

	for (int i = 0; i < 6; i++) {

		/* Initialize Display */
		status |= GzInitDisplay(m_pDisplay[i]);

		/*
		* Tokens associated with triangle vertex values
		*/
		nameListTriangle[0] = GZ_POSITION;
		nameListTriangle[1] = GZ_NORMAL;
		nameListTriangle[2] = GZ_TEXTURE_INDEX;

		// I/O File open

		if ((infile = fopen(INFILE, "r")) == NULL)
		{
			AfxMessageBox("The input file was not opened\n");
			return GZ_FAILURE;
		}


		if ((outfile = fopen(OUTFILE, "wb")) == NULL)
		{
			AfxMessageBox("The output file was not opened\n");
			return GZ_FAILURE;
		}

		/*
		* Walk through the list of triangles, set color
		* and render each triangle
		*/
		while (fscanf(infile, "%s", dummy) == 1) { 	/* read in tri word */
			fscanf(infile, "%f %f %f %f %f %f %f %f",
				&(vertexList[0][0]), &(vertexList[0][1]),
				&(vertexList[0][2]),
				&(normalList[0][0]), &(normalList[0][1]),
				&(normalList[0][2]),
				&(uvList[0][0]), &(uvList[0][1]));
			fscanf(infile, "%f %f %f %f %f %f %f %f",
				&(vertexList[1][0]), &(vertexList[1][1]),
				&(vertexList[1][2]),
				&(normalList[1][0]), &(normalList[1][1]),
				&(normalList[1][2]),
				&(uvList[1][0]), &(uvList[1][1]));
			fscanf(infile, "%f %f %f %f %f %f %f %f",
				&(vertexList[2][0]), &(vertexList[2][1]),
				&(vertexList[2][2]),
				&(normalList[2][0]), &(normalList[2][1]),
				&(normalList[2][2]),
				&(uvList[2][0]), &(uvList[2][1]));

			/*
			* Set the value pointers to the first vertex of the
			* triangle, then feed it to the renderer
			* NOTE: this sequence matches the nameList token sequence
			*/
			valueListTriangle[0] = (GzPointer)vertexList;
			valueListTriangle[1] = (GzPointer)normalList;
			valueListTriangle[2] = (GzPointer)uvList;
			GzPutTriangle(m_pRender[i], 3, nameListTriangle, valueListTriangle);
		}
	}

	for (int i = 0; i < m_pDisplay[0]->xres * m_pDisplay[0]->yres; i++) {

		GzPixel pixel;
		pixel.red = pixel.green = pixel.blue = pixel.z = 0;
		for (int j = 0; j < 6; j++) {
			pixel.red += m_pDisplay[j]->fbuf[i].red * AAFilter[j][2];
			pixel.green += m_pDisplay[j]->fbuf[i].green * AAFilter[j][2];
			pixel.blue += m_pDisplay[j]->fbuf[i].blue * AAFilter[j][2];
			pixel.z += m_pDisplay[j]->fbuf[i].z * AAFilter[j][2];

		}
		m_pDisplay[0]->fbuf[i] = pixel;

	}

	GzFlushDisplay2File(outfile, m_pDisplay[0]);
	GzFlushDisplay2FrameBuffer(m_pFrameBuffer, m_pDisplay[0]);	

	if (fclose(infile))
		AfxMessageBox("The input file was not closed\n");

	if (fclose(outfile))
		AfxMessageBox("The output file was not closed\n");

	if (status)
		return(GZ_FAILURE);
	else
		return(GZ_SUCCESS);
}
int Application5::Clean()
{
	/*
	* Clean up and exit
	*/

	int	status = 0;
	for (int i = 0; i < 6; i++) {
		status |= GzFreeRender(m_pRender[i]);
		status |= GzFreeDisplay(m_pDisplay[i]);
	}
	if (status)
		return(GZ_FAILURE);
	else
		return(GZ_SUCCESS);
}