// Application.cpp: implementation of the Application class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "CS580HW.h"
#include "Application.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Application::Application()
{
	for (int i = 0; i < AAKERNEL_SIZE; i++) {
		m_pDisplay[i] = NULL;		// the display
		m_pRender[i] = NULL;		// the renderer
		m_pUserInput = NULL;
		m_pFrameBuffer = NULL;
	}
}

Application::~Application()
{
	if(m_pFrameBuffer != NULL)
		delete m_pFrameBuffer;
}

