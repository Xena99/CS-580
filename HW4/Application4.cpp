// Application4.cpp: implementation of the Application4 class.
//
//////////////////////////////////////////////////////////////////////

/*
* application test code for homework assignment
*/

#include "stdafx.h"
#include "CS580HW.h"
#include "Application4.h"
#include "Gz.h"
#include "disp.h"
#include "rend.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#define new DEBUG_NEW
#endif

#define INFILE  "teapot4.asc"
#define OUTFILE "output.ppm"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Application4::Application4()
{

}

Application4::~Application4()
{

}

int Application4::Initialize()
{
	/* to be filled in by the app if it sets camera params */

	GzCamera	camera;
	int		    xRes, yRes, dispClass;	/* display parameters */

	GzToken		nameListShader[9]; 	/* shader attribute names */
	GzPointer   valueListShader[9];		/* shader attribute pointers */
	GzToken     nameListLights[10];		/* light info */
	GzPointer   valueListLights[10];
	int			shaderType, interpStyle;
	float		specpower;
	int		    status;

	status = 0;

	/*
	* Allocate memory for user input
	*/
	m_pUserInput = new GzInput;

	/*
	* initialize the display and the renderer
	*/

	m_nWidth = 400;		// frame buffer and display width
	m_nHeight = 400;    // frame buffer and display height

	status |= GzNewFrameBuffer(&m_pFrameBuffer, m_nWidth, m_nHeight);

	status |= GzNewDisplay(&m_pDisplay, GZ_RGBAZ_DISPLAY, m_nWidth, m_nHeight);

	status |= GzGetDisplayParams(m_pDisplay, &xRes, &yRes, &dispClass);

	status |= GzInitDisplay(m_pDisplay);

	status |= GzNewRender(&m_pRender, GZ_Z_BUFFER_RENDER, m_pDisplay);

	/* Translation matrix */
GzMatrix	scale =
{
	5.5,	0.0,	0.0,	9.5,
	0.0,	5.5,	0.0,	-2,
	0.0,	0.0,	5.5,	45,
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

#if 1 	/* set up app-defined camera if desired, else use camera defaults */
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

	status |= GzPutCamera(m_pRender, &camera);
#endif 

	/* Start Renderer */
	status |= GzBeginRender(m_pRender);

	/* Light */
	GzLight	light1 = { { -0.7071, 0.7071, 0 },{ 0.5, 0.4, 0.6 } };
	GzLight	light2 = { { 0, -0.7071, -0.7071 },{ 0.4, 0.4, 0.3 } };
	GzLight	light3 = { { 0.7071, 0.0, -0.7071 },{ 0.4, 0.35, 0.5 } };
	GzLight	ambientlight = { { 0, 0, 0 },{ 0.7, 0.7, 0.7 } };

	/* Material property */
	GzColor specularCoefficient = { 0.4, 0.4, 0.4 };
	GzColor ambientCoefficient = { 0.3, 0.3, 0.3 };
	GzColor diffuseCoefficient = { 0.9, 0.9, 0.9 };

	/*
	renderer is ready for frame --- define lights and shader at start of frame
	*/

	/*
	* Tokens associated with light parameters
	*/
	nameListLights[0] = GZ_DIRECTIONAL_LIGHT;
	valueListLights[0] = (GzPointer)&light1;
	nameListLights[1] = GZ_DIRECTIONAL_LIGHT;
	valueListLights[1] = (GzPointer)&light2;
	nameListLights[2] = GZ_DIRECTIONAL_LIGHT;
	valueListLights[2] = (GzPointer)&light3;
	status |= GzPutAttribute(m_pRender, 3, nameListLights, valueListLights);

	nameListLights[0] = GZ_AMBIENT_LIGHT;
	valueListLights[0] = (GzPointer)&ambientlight;
	status |= GzPutAttribute(m_pRender, 1, nameListLights, valueListLights);

	/*
	* Tokens associated with shading
	*/
	nameListShader[0] = GZ_DIFFUSE_COEFFICIENT;
	valueListShader[0] = (GzPointer)diffuseCoefficient;

	/*
	* Select either GZ_COLOR or GZ_NORMALS as interpolation mode
	*/
	nameListShader[1] = GZ_INTERPOLATE;
	//interpStyle = GZ_FLAT;			 /*For Flat shading*/
#if 0
	interpStyle = GZ_COLOR;         /* Gouraud shading */
#else 
	interpStyle = GZ_NORMALS;       /* Phong shading */
#endif

	valueListShader[1] = (GzPointer)&interpStyle;
	nameListShader[2] = GZ_AMBIENT_COEFFICIENT;
	valueListShader[2] = (GzPointer)ambientCoefficient;
	nameListShader[3] = GZ_SPECULAR_COEFFICIENT;
	valueListShader[3] = (GzPointer)specularCoefficient;
	nameListShader[4] = GZ_DISTRIBUTION_COEFFICIENT;
	specpower = 32;
	valueListShader[4] = (GzPointer)&specpower;

	status |= GzPutAttribute(m_pRender, 5, nameListShader, valueListShader);

	status |= GzPushMatrix(m_pRender, scale);
	status |= GzPushMatrix(m_pRender, rotateY);
	status |= GzPushMatrix(m_pRender, rotateX);

	if (status) exit(GZ_FAILURE);

	if (status)
		return(GZ_FAILURE);
	else
		return(GZ_SUCCESS);
}

int Application4::Render()
{
	GzToken		nameListTriangle[3]; 	/* vertex attribute names */
	GzPointer	valueListTriangle[3]; 	/* vertex attribute pointers */
	GzCoord		vertexList[3];	/* vertex position coordinates */
	GzCoord		normalList[3];	/* vertex normals */
	GzTextureIndex  	uvList[3];		/* vertex texture map indices */
	char		dummy[256];
	int			status;


	/* Initialize Display */
	status |= GzInitDisplay(m_pDisplay);

	/*
	* Tokens associated with triangle vertex values
	*/
	nameListTriangle[0] = GZ_POSITION;
	nameListTriangle[1] = GZ_NORMAL;

	// I/O File open
	FILE *infile;
	if ((infile = fopen(INFILE, "r")) == NULL)
	{
		AfxMessageBox("The input file was not opened\n");
		return GZ_FAILURE;
	}

	FILE *outfile;
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
		GzPutTriangle(m_pRender, 2, nameListTriangle, valueListTriangle);
	}

	GzFlushDisplay2File(outfile, m_pDisplay); 	/* write out or update display to file*/
	GzFlushDisplay2FrameBuffer(m_pFrameBuffer, m_pDisplay);	// write out or update display to frame buffer

															/*
															* Close file
															*/

	if (fclose(infile))
		AfxMessageBox("The input file was not closed\n");

	if (fclose(outfile))
		AfxMessageBox("The output file was not closed\n");

	if (status)
		return(GZ_FAILURE);
	else
		return(GZ_SUCCESS);
}

int Application4::Clean()
{
	/*
	* Clean up and exit
	*/
	int	status = 0;

	status |= GzFreeRender(m_pRender);
	status |= GzFreeDisplay(m_pDisplay);

	if (status)
		return(GZ_FAILURE);
	else
		return(GZ_SUCCESS);
}
