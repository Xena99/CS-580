#pragma once
#include "disp.h" /* include your own disp.h file (e.g. hw1)*/
#include <map>
/* Camera defaults */
#define	DEFAULT_FOV		35.0
#define	DEFAULT_IM_Z	(-10.0)  /* world coords for image plane origin */
#define	DEFAULT_IM_Y	(5.0)    /* default look-at point = 0,0,0 */
#define	DEFAULT_IM_X	(-10.0)

#define	DEFAULT_AMBIENT	{0.1, 0.1, 0.1}
#define	DEFAULT_DIFFUSE	{0.77, 0.6, 0.5}
#define	DEFAULT_SPECULAR	{0.2, 0.3, 0.4}
#define	DEFAULT_SPEC		32

#define	MATLEVELS	100		/* how many matrix pushes allowed */
#define	MAX_LIGHTS	10		/* how many lights allowed */

#ifndef GZRENDER
#define GZRENDER
typedef struct {			/* define a renderer */
	int normal_level;
	GzRenderClass	renderClass;
	GzDisplay		*display;
	short		    open;
	GzCamera		camera;
	short		    matlevel;	        /* top of stack - current xform */
	GzMatrix		Ximage[MATLEVELS];	/* stack of xforms (Xsm) */
	GzMatrix		Xnorm[MATLEVELS];	/* xforms for norms (Xim) */
	GzMatrix		Xsp;		        /* NDC to screen (pers-to-screen) */
	GzColor		flatcolor;          /* color state for flat shaded triangles */
	int			interp_mode;
	int			numlights;
	GzLight		lights[MAX_LIGHTS];
	GzLight		ambientlight;
	GzColor		Ka, Kd, Ks;
	float		    spec;		/* specular power */
	GzTexture		tex_fun;    /* tex_fun(float u, float v, GzColor color) */
	float dx;
	float dy;
}  GzRender;
#endif

// Function declaration
// HW2
int GzNewRender(GzRender **render, GzRenderClass renderClass, GzDisplay *display);
int GzFreeRender(GzRender *render);
int GzBeginRender(GzRender	*render);
int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList,
	GzPointer *valueList);
int GzPutTriangle(GzRender *render, int	numParts, GzToken *nameList,
	GzPointer *valueList);

// HW3
int GzPutCamera(GzRender *render, GzCamera *camera);
int GzPushMatrix(GzRender *render, GzMatrix	matrix, bool isProjectionMatrix);
int GzPopMatrix(GzRender *render);
void Normalize(GzCoord& vec);

// Object Translation
int GzRotXMat(float degree, GzMatrix mat);
int GzRotYMat(float degree, GzMatrix mat);
int GzRotZMat(float degree, GzMatrix mat);
int GzTrxMat(GzCoord translate, GzMatrix mat);
int GzScaleMat(GzCoord scale, GzMatrix mat);
