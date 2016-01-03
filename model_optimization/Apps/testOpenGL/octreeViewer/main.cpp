/*! \file    main.c
    \ingroup demos

    This program is drawing octave space in 3D by using
    OpenGLUT.  It may also be useful to see the design procedure
	in computation fabrication .
 
    Wireframe octave tree space are displayed after mesh model loaded. 
	Some parameters can be adjusted.
 
   Keys:
      -    <tt>Esc &nbsp;</tt> Quit
      -    <tt>q Q &nbsp;</tt> Quit
      -    <tt>i I &nbsp;</tt> Show info
      -    <tt>p P &nbsp;</tt> Toggle perspective or orthographic projection
      -    <tt>r R &nbsp;</tt> Toggle fixed or animated rotation around model X-axis
      -    <tt>s S &nbsp;</tt> Toggle toggle fixed function or shader render path
      -    <tt>n N &nbsp;</tt> Toggle visualization of object's normal vectors
      -    <tt>= + &nbsp;</tt> Increase \a slices
      -    <tt>- _ &nbsp;</tt> Decreate \a slices
      -    <tt>, < &nbsp;</tt> Decreate \a stacks
      -    <tt>. > &nbsp;</tt> Increase \a stacks
      -    <tt>9 ( &nbsp;</tt> Decreate \a depth  (Sierpinski Sponge)
      -    <tt>0 ) &nbsp;</tt> Increase \a depth  (Sierpinski Sponge)
      -    <tt>up&nbsp; &nbsp;</tt> Increase "outer radius"
      -    <tt>down&nbsp;</tt> Decrease "outer radius"
      -    <tt>left&nbsp;</tt> Decrease "inner radius"
      -    <tt>right</tt> Increase "inner radius"
      -    <tt>PgUp&nbsp;</tt> Next shape-drawing function
      -    <tt>PgDn&nbsp;</tt> Prev shape-drawing function

    \author  Written by Yuan Yao Aug. 2015

    \author  Portions Copyright (C) 2004, RMEC.
 
    \image   html boxIndex.png 
    \include octree.cpp
*/

#include <GL/freeglut.h>

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "glmatrix.h"

#include <vector>
#include "../config.h"
#include "suOctree.h"

#ifdef _MSC_VER
/* DUMP MEMORY LEAKS */
#include <crtdbg.h>
#endif
/*
* These global variables control which object is drawn,
* and how it is drawn.  No object uses all of these
* variables.
*/
static int function_index;
static int slices = 16;
static int stacks = 16;
static double yRot = -45;
static double xRot = 45;
static double irad = .25;
static double orad = 1.0;   /* doubles as size for objects other than Torus */
static int depth = 4;
static double offset[3] = { 0, 0, 0 };
static GLboolean show_info = GL_TRUE;
static float ar;
static GLboolean persProject = GL_TRUE;
static GLboolean animateXRot = GL_FALSE;
static GLboolean useShader = GL_FALSE;
static GLboolean visNormals = GL_FALSE;

static int octree_level = 2;//octree level
/*
* Enum to tell drawSizeInfo what to draw for each object
*/
#define GEO_NO_SIZE 0
#define GEO_SIZE 1
#define GEO_SCALE 2
#define GEO_INNER_OUTER_RAD 4
#define GEO_RAD 8
#define GEO_BASE_HEIGHT 16
#define GEO_RAD_HEIGHT 32

/* -- Cube -- */
#define CUBE_NUM_VERT           8
#define CUBE_NUM_FACES          6
#define CUBE_NUM_EDGE_PER_FACE  4
#define CUBE_VERT_PER_OBJ       (CUBE_NUM_FACES*CUBE_NUM_EDGE_PER_FACE)
#define CUBE_VERT_ELEM_PER_OBJ  (CUBE_VERT_PER_OBJ*3)
#define CUBE_VERT_PER_OBJ_TRI   (CUBE_VERT_PER_OBJ+CUBE_NUM_FACES*2)    /* 2 extra edges per face when drawing quads as triangles */

/* Normal Vectors */
static GLfloat cube_n[CUBE_NUM_FACES * 3] =
{
	0.0f, 0.0f, 1.0f,
	1.0f, 0.0f, 0.0f,
	0.0f, 1.0f, 0.0f,
	-1.0f, 0.0f, 0.0f,
	0.0f, -1.0f, 0.0f,
	0.0f, 0.0f, -1.0f
};

/* Vertex indices, as quads, before triangulation */
static GLubyte cube_vi[CUBE_VERT_PER_OBJ] =
{
	7, 4, 5, 6,
	3, 7, 6, 2,
	2, 6, 5, 1,
	5, 4, 0, 1,
	3, 0, 4, 7,	
	2, 1, 0, 3
};

suOctree octree;
Rectangle3f bbox = Rectangle3f(Vec3f(-5, -5, -5), Vec3f(5, 5, 5));

class suCube
{
public:
	suCube(Vec3f minCord, Vec3f maxCord){ box_ = Rectangle3f(minCord, maxCord); }
	suCube(Vec3f minCord, float dx, float dy, float dz){ box_ = Rectangle3f(minCord, minCord + Vec3f(dx, dy, dz)); }
	suCube(){};

	void genVertexs();

public:
	std::vector<Vec3f> vertArr_;    //8 vertex
	Rectangle3f box_;
};

void suCube::genVertexs()
{
	vertArr_.resize(8);
	Vec3f dim = box_.Dimensions();

	vertArr_[0] = box_.Min;
	vertArr_[1] = box_.Min + Vec3f(0.0f, dim.y, 0.0f);
	vertArr_[2] = box_.Min + Vec3f(dim.x, dim.y, 0.0f);
	vertArr_[3] = box_.Min + Vec3f(dim.x, 0.0f, 0.0f);

	vertArr_[4] = box_.Min + Vec3f(0.0f, 0.0f, dim.z);
	vertArr_[5] = box_.Min + Vec3f(0.0f, dim.y, dim.z);
	vertArr_[6] = box_.Max;
	vertArr_[7] = box_.Min + Vec3f(dim.x, 0.0f, dim.z);
}


void drawBox(GLfloat *pVert, GLubyte *pIndice)
{
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);

	glVertexPointer(3, GL_FLOAT, 0, pVert);
	glNormalPointer(GL_FLOAT, 0, cube_n);

	for (int i = 0; i < CUBE_NUM_FACES; i++)
	{
		//std::cout << i*CUBE_NUM_EDGE_PER_FACE << std::endl;
		//glDrawArrays(GL_LINE_LOOP, i*CUBE_NUM_EDGE_PER_FACE, CUBE_NUM_EDGE_PER_FACE);
		glDrawElements(GL_LINE_LOOP, CUBE_NUM_EDGE_PER_FACE, GL_UNSIGNED_BYTE, cube_vi + i*CUBE_NUM_EDGE_PER_FACE);
	}
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
}

void drawCube()
{
	int zSize = 5;
	int xSize = 2;
	int ySize = 2;

	//生成所有六面体空间节点，用vector存储节点数组
	std::vector<Vec3f> points;
	std::vector<suCube> cubes;

	for (int z = -zSize / 2; z <= zSize / 2; z++)
	{
		for (int y = -ySize / 2; y <= ySize / 2; y++)
		{
			for (int x = -xSize / 2; x <= xSize / 2; x++)
			{
				Vec3f p;
				p.x = x; p.y = y; p.z = z;
				points.push_back(p);
			}
		}
	}
	
	for (int i = 0; i < points.size(); i++)
	{
		std::vector<GLfloat> vvert;
		suCube cube(points[i], 1, 1, 1);
		cube.genVertexs();
		cubes.push_back(cube);
		for (int j = 0; j < 8; j++)
		{
			vvert.push_back(cube.vertArr_[j].x);
			vvert.push_back(cube.vertArr_[j].y);
			vvert.push_back(cube.vertArr_[j].z);
		}
		GLfloat *pVert = &vvert[0];

		drawBox(pVert, cube_vi);
	}

}

void initOctree()
{
	octree.DeleteTree();
	octree.SetBox(bbox.Min.x, bbox.Min.y, bbox.Min.z,
		          bbox.Max.x, bbox.Max.y, bbox.Max.z);
	octree.Generate(octree_level);
}

void drawOctree()
{

	//生成所有六面体空间节点，用vector存储节点数组
	std::vector<Vec3f> points;
	std::vector<suCube> cubes;

	//Transverse all leaf node
	std::vector<suCell*>  LeafCellVector;
	suCell *pCell = octree.GetRoot();
	if (pCell == NULL)
		throw suException("No octree exist!");

	queue<suCell*> CellQueue;
	CellQueue.push(pCell);

	while (!CellQueue.empty())
	{
		pCell = CellQueue.front();
		CellQueue.pop();

		if (pCell == NULL)  continue;
		
		if (pCell->pChildren != NULL)
		{
			for (unsigned int i = 0; i < 8; i++)
			{
				suCell *pChildCell = pCell->pChildren + i;
				CellQueue.push(pChildCell);
			}
		}
		else{
			LeafCellVector.push_back(pCell);
		}
	}

	//Bounding box size
	float minx, miny, minz, maxx, maxy, maxz;
	octree.GetBox(minx, miny, minz, maxx, maxy, maxz);
	float xBBSize = maxx - minx;
	float yBBSize = maxy - miny;
	float zBBSize = maxz - minz;
	
	for (int i = 0; i < LeafCellVector.size(); i++)
	{
		std::vector<GLfloat> vvert;
		
		float fx, fy, fz;
		octree.GetLocation(LeafCellVector[i], fx, fy, fz);
		
		float nDevide = pow(2, LeafCellVector[i]->level);
		float xBox_ = xBBSize / nDevide;
		float yBox_ = yBBSize / nDevide;
		float zBox_ = zBBSize / nDevide;

		suCube cube(Vec3f(fx - xBox_ / 2, fy - yBox_ / 2, fz - zBox_ / 2), 
			        Vec3f(fx + xBox_ / 2, fy + yBox_ / 2, fz + zBox_ / 2));
		

		cube.genVertexs();
		cubes.push_back(cube);
		for (int j = 0; j < 8; j++)
		{
			vvert.push_back(cube.vertArr_[j].x);
			vvert.push_back(cube.vertArr_[j].y);
			vvert.push_back(cube.vertArr_[j].z);
		}
		GLfloat *pVert = &vvert[0];

		drawBox(pVert, cube_vi);
	}

}

/* report GL errors, if any, to stderr */
void checkError(const char *functionName)
{
    GLenum error;
    while (( error = glGetError() ) != GL_NO_ERROR) {
        fprintf (stderr, "GL error 0x%X detected in %s\n", error, functionName);
    }
}

/*
 * OpenGL 2+ shader mode needs some function and macro definitions, 
 * avoiding a dependency on additional libraries like GLEW or the
 * GL/glext.h header
 */
#ifndef GL_FRAGMENT_SHADER
#define GL_FRAGMENT_SHADER 0x8B30
#endif

#ifndef GL_VERTEX_SHADER
#define GL_VERTEX_SHADER 0x8B31
#endif

#ifndef GL_COMPILE_STATUS
#define GL_COMPILE_STATUS 0x8B81
#endif

#ifndef GL_LINK_STATUS
#define GL_LINK_STATUS 0x8B82
#endif

#ifndef GL_INFO_LOG_LENGTH
#define GL_INFO_LOG_LENGTH 0x8B84
#endif

typedef ptrdiff_t ourGLsizeiptr;
typedef char ourGLchar;

#ifndef APIENTRY
#define APIENTRY
#endif

#ifndef GL_VERSION_2_0
typedef GLuint (APIENTRY *PFNGLCREATESHADERPROC) (GLenum type);
typedef void (APIENTRY *PFNGLSHADERSOURCEPROC) (GLuint shader, GLsizei count, const ourGLchar **string, const GLint *length);
typedef void (APIENTRY *PFNGLCOMPILESHADERPROC) (GLuint shader);
typedef GLuint (APIENTRY *PFNGLCREATEPROGRAMPROC) (void);
typedef void (APIENTRY *PFNGLATTACHSHADERPROC) (GLuint program, GLuint shader);
typedef void (APIENTRY *PFNGLLINKPROGRAMPROC) (GLuint program);
typedef void (APIENTRY *PFNGLUSEPROGRAMPROC) (GLuint program);
typedef void (APIENTRY *PFNGLGETSHADERIVPROC) (GLuint shader, GLenum pname, GLint *params);
typedef void (APIENTRY *PFNGLGETSHADERINFOLOGPROC) (GLuint shader, GLsizei bufSize, GLsizei *length, ourGLchar *infoLog);
typedef void (APIENTRY *PFNGLGETPROGRAMIVPROC) (GLenum target, GLenum pname, GLint *params);
typedef void (APIENTRY *PFNGLGETPROGRAMINFOLOGPROC) (GLuint program, GLsizei bufSize, GLsizei *length, ourGLchar *infoLog);
typedef GLint (APIENTRY *PFNGLGETATTRIBLOCATIONPROC) (GLuint program, const ourGLchar *name);
typedef GLint (APIENTRY *PFNGLGETUNIFORMLOCATIONPROC) (GLuint program, const ourGLchar *name);
typedef void (APIENTRY *PFNGLUNIFORMMATRIX4FVPROC) (GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
typedef void (APIENTRY *PFNGLUNIFORMMATRIX3FVPROC) (GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
#endif

PFNGLCREATESHADERPROC gl_CreateShader;
PFNGLSHADERSOURCEPROC gl_ShaderSource;
PFNGLCOMPILESHADERPROC gl_CompileShader;
PFNGLCREATEPROGRAMPROC gl_CreateProgram;
PFNGLATTACHSHADERPROC gl_AttachShader;
PFNGLLINKPROGRAMPROC gl_LinkProgram;
PFNGLUSEPROGRAMPROC gl_UseProgram;
PFNGLGETSHADERIVPROC gl_GetShaderiv;
PFNGLGETSHADERINFOLOGPROC gl_GetShaderInfoLog;
PFNGLGETPROGRAMIVPROC gl_GetProgramiv;
PFNGLGETPROGRAMINFOLOGPROC gl_GetProgramInfoLog;
PFNGLGETATTRIBLOCATIONPROC gl_GetAttribLocation;
PFNGLGETUNIFORMLOCATIONPROC gl_GetUniformLocation;
PFNGLUNIFORMMATRIX4FVPROC gl_UniformMatrix4fv;
PFNGLUNIFORMMATRIX3FVPROC gl_UniformMatrix3fv;
 







/*!
    Does printf()-like work using freeglut
    glutBitmapString().  Uses a fixed font.  Prints
    at the indicated row/column position.

    Limitation: Cannot address pixels.
    Limitation: Renders in screen coords, not model coords.
*/
static void shapesPrintf (int row, int col, const char *fmt, ...)
{
    static char buf[256];
    int viewport[4];
    void *font = GLUT_BITMAP_9_BY_15;
    va_list args;

    va_start(args, fmt);
#if defined(WIN32) && !defined(__CYGWIN__)
    (void) _vsnprintf (buf, sizeof(buf), fmt, args);
#else
    (void) vsnprintf (buf, sizeof(buf), fmt, args);
#endif
    va_end(args);

    glGetIntegerv(GL_VIEWPORT,viewport);

    glPushMatrix();
    glLoadIdentity();

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

        glOrtho(0,viewport[2],0,viewport[3],-1,1);

        glRasterPos2i
        (
              glutBitmapWidth(font, ' ') * col,
            - glutBitmapHeight(font) * row + viewport[3]
        );
        glutBitmapString (font, (unsigned char*)buf);

    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}

/* Print info about the about the current shape and render state on the screen */
static void DrawSizeInfo(int *row)
{
	switch (GEO_SIZE)
    {
    case GEO_NO_SIZE:
        break;
    case GEO_SIZE:
        shapesPrintf ((*row)++, 1, "Size  Up  Down : %f", orad);
        break;
    case GEO_SCALE:
        shapesPrintf ((*row)++, 1, "Scale  Up  Down : %f", orad);
        break;
    case GEO_INNER_OUTER_RAD:
        shapesPrintf ((*row)++, 1, "Inner radius Left Right: %f", irad);
        shapesPrintf ((*row)++, 1, "Outer radius  Up  Down : %f", orad);
        break;
    case GEO_RAD:
        shapesPrintf ((*row)++, 1, "Radius  Up  Down : %f", orad);
        break;
    case GEO_BASE_HEIGHT:
        shapesPrintf ((*row)++, 1, "Base   Left Right: %f", irad);
        shapesPrintf ((*row)++, 1, "Height  Up  Down : %f", orad);
        break;
    case GEO_RAD_HEIGHT:
        shapesPrintf ((*row)++, 1, "Radius Left Right: %f", irad);
        shapesPrintf ((*row)++, 1, "Height  Up  Down : %f", orad);
        break;
    }
}

static void drawInfo()
{
    int row = 1;
   // shapesPrintf (row++, 1, "Shape PgUp PgDn: %s", table [function_index].name);
    shapesPrintf (row++, 1, "Octree level +-: %d", octree_level);
    shapesPrintf (row++, 1, "nSides +-: %d   nRings <>: %d", slices, stacks);
    shapesPrintf (row++, 1, "Depth  (): %d", depth);
    DrawSizeInfo(&row);
    if (persProject)
        shapesPrintf (row++, 1, "Perspective projection (p)");
    else
        shapesPrintf (row++, 1, "Orthographic projection (p)");
    if (useShader)
        shapesPrintf (row++, 1, "Using shader (s)");
    else
        shapesPrintf (row++, 1, "Using fixed function pipeline (s)");
    if (animateXRot)
        shapesPrintf (row++, 1, "2D rotation (r)");
    else
        shapesPrintf (row++, 1, "1D rotation (r)");
    shapesPrintf (row++, 1, "visualizing normals: %i (n)",visNormals);
}

/* GLUT callback Handlers */
static void
resize(int width, int height)
{
    ar = (float) width / (float) height;

    glViewport(0, 0, width, height);
}

static void display(void)
{
    const double t = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
    const double a = t*89.0;
    const double b = (animateXRot?t:1)*67.0;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glutSetOption(GLUT_GEOMETRY_VISUALIZE_NORMALS,visNormals);  /* Normals visualized or not? */

	/* fixed function pipeline */
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if (persProject)
		//glFrustum(-ar, ar, -1.0, 1.0, 2.0, 100.0);
		glFrustum(-ar, ar, -1.0, 1.0, 2.0, 100.0);
	else
		glOrtho(-ar * 3, ar * 3, -3.0, 3.0, 2.0, 100.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glEnable(GL_LIGHTING);

	glColor3d(1, 0, 0);

	glPushMatrix();
	glTranslated(0, 1.2, -6);
	glRotated(b, 1, 0, 0);
	glRotated(a, 0, 0, 1);
	// table [function_index].solid ();
	glPopMatrix();

	glPushMatrix();
	glTranslated(0, -1.2, -20);
	glRotated(xRot, 1, 0, 0);
	glRotated(yRot, 0, 1, 0);
	//  table [function_index].wire ();
     	//drawCube();
	drawOctree();
	glPopMatrix();

	glDisable(GL_LIGHTING);
	glColor3d(0.1, 0.1, 0.4);
    
   
    if( show_info )
        /* print info to screen */
        drawInfo();
    else
        /* print to command line instead */
        printf ( "Shape %d slides %d stacks %d\n", function_index, slices, stacks ) ;

    glutSwapBuffers();
}


static void
key(unsigned char key, int x, int y)
{
    switch (key)
    {
    case 27 :
    case 'Q':
    case 'q': glutLeaveMainLoop () ;      break;

    case 'I':
    case 'i': show_info=!show_info;       break;

    case '=':
	case '+': octree_level++;    initOctree();         break;

    case '-':
	case '_': if (octree_level > 0) { octree_level--;  initOctree(); }    break;

    case ',':
    case '<': if( stacks > -1 ) stacks--; break;

    case '.':
    case '>': stacks++;                   break;

    case '9': 
    case '(': if( depth > -1 ) depth--;   break;

    case '0': 
    case ')': ++depth;                    break;

    case 'P':
    case 'p': persProject=!persProject;   break;

    case 'R':
    case 'r': animateXRot=!animateXRot;   break;

    case 'S':
    case 's':
      
        break;

    case 'N':
    case 'n': visNormals=!visNormals;     break;

    default:
        break;
    }

    glutPostRedisplay();
}

static void special (int key, int x, int y)
{
    switch (key)
    {
    case GLUT_KEY_PAGE_UP:    ++function_index; break;
    case GLUT_KEY_PAGE_DOWN:  --function_index; break;
    case GLUT_KEY_UP:         xRot -= 5;        break;
    case GLUT_KEY_DOWN:       xRot += 5;        break;

    case GLUT_KEY_RIGHT:      yRot += 5;        break;
    case GLUT_KEY_LEFT:       yRot -= 5;        break;

    default:
        break;
    }

}


static void
idle(void)
{
    glutPostRedisplay();
}

const GLfloat light_ambient[]  = { 0.0f, 0.0f, 0.0f, 1.0f };
const GLfloat light_diffuse[]  = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_position[] = { 2.0f, 5.0f, 5.0f, 0.0f };

const GLfloat mat_ambient[]    = { 0.7f, 0.7f, 0.7f, 1.0f };
const GLfloat mat_diffuse[]    = { 0.8f, 0.8f, 0.8f, 1.0f };
const GLfloat mat_specular[]   = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat high_shininess[] = { 100.0f };

/* Program entry point */

int main(int argc, char *argv[])
{
	initOctree();

    glutInitWindowSize(800,600);
    glutInitWindowPosition(40,40);
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);

    glutCreateWindow("Octree Viewer");

    glutReshapeFunc(resize);
    glutDisplayFunc(display);
    glutKeyboardFunc(key);
    glutSpecialFunc(special);
    glutIdleFunc(idle);

    glutSetOption ( GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION ) ;

    glClearColor(1,1,1,1);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);

    glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);

    glMaterialfv(GL_FRONT, GL_AMBIENT,   mat_ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,   mat_diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR,  mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, high_shininess);

    glutMainLoop();

#ifdef _MSC_VER
    /* DUMP MEMORY LEAK INFORMATION */
    _CrtDumpMemoryLeaks () ;
#endif

    return EXIT_SUCCESS;
}
