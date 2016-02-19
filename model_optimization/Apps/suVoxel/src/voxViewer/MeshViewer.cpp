//=============================================================================
//                                                                            
//   Example code for the full-day course
//
//   M. Botsch, M. Pauly, L. Kobbelt, P. Alliez, B. Levy,
//   "Geometric Modeling Based on Polygonal Meshes"
//   held at SIGGRAPH 2007, San Diego, USA.
//
//   Copyright (C) 2007 by  Computer Graphics Laboratory, ETH Zurich, 
//                      and Computer Graphics Group,      RWTH Aachen
//
//                                                                            
//-----------------------------------------------------------------------------
//                                                                            
//                                License                                     
//                                                                            
//   This program is free software; you can redistribute it and/or
//   modify it under the terms of the GNU General Public License
//   as published by the Free Software Foundation; either version 2
//   of the License, or (at your option) any later version.
//   
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//   
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 51 Franklin Street, Fifth Floor, 
//   Boston, MA  02110-1301, USA.
//                                                                            
//=============================================================================
//=============================================================================
//
//  CLASS MeshViewer - IMPLEMENTATION
//
//=============================================================================


//== INCLUDES =================================================================

#include "tinyfiledialogs.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "MeshViewer.h"
#include "gl.hh"
#include <iostream>
#include <fstream>
#include "suHyperMesh.h"
#include "tinyfiledialogs.h"
#include <Windows.h>


#define  MENU_FILE_LOAD 1111
#define  MENU_LABEL_BOUNDARY 1112
#define  MENU_LABEL_BNEIGHBOR 1113
#define  MENU_FILE_SAVE 1121
#define  test_ 1
#define  test1_ 2
//bar function;
void TW_CALL getBool(void *value, void *clientData)
{
	*((bool*)value) = *((bool*)clientData);
}

void TW_CALL setBool(const void *value, void *clientData)
{
	*((bool*)clientData) = *((bool*)value);
}


//== IMPLEMENTATION ========================================================== 


MeshViewer::
MeshViewer(const char* _title, int _width, int _height)
: GlutExaminer(_title, _width, _height)
{
	mesh_.request_face_normals();
	mesh_.request_vertex_normals();

	clear_draw_modes();
	add_draw_mode("Wireframe");
	add_draw_mode("Hidden Line");
	add_draw_mode("Solid Flat");
	add_draw_mode("Solid Smooth");
	set_draw_mode(3);

	glutAddMenuEntry("Load...", MENU_FILE_LOAD);
	glutAddMenuEntry("Save...", MENU_FILE_SAVE);
	glutAddMenuEntry("Label IN/OUT", MENU_LABEL_BOUNDARY);
	glutAddMenuEntry("Floodfill", MENU_LABEL_BNEIGHBOR);
	glutAddMenuEntry("test",test_);
	glutAddMenuEntry("test1", test1_);
	octreeForDraw_level_ = 0;
	num_grid = 0;

	init();
}
//-----------------------------------------------------------------------------


bool
MeshViewer::
open_mesh(const char* _filename)
{
	//replace mesh & octave with hyperMesh
	if (open_hyperMesh(_filename))
	{
		mesh_.clear();
		mesh_ = hyperMesh_.mesh;
		// compute face & vertex normals
		mesh_.update_normals();

		// update face indices for faster rendering
		update_face_indices();
		// info
		std::cerr << mesh_.n_vertices() << " vertices, "
			<< mesh_.n_faces() << " faces\n";

		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------

bool 
MeshViewer::open_hyperMesh(const char* _filename)
{
	// load mesh
	hyperMesh_.clear();
	if (!hyperMesh_.LoadMeshFromFile(_filename))
	{
		return false;
	}

	set_scene((OpenMesh::Vec3f)(hyperMesh_.bbMin + hyperMesh_.bbMax)*0.5, 0.5*(hyperMesh_.bbMin - hyperMesh_.bbMax).norm());

	int level = 3;
	hyperMesh_.PartitionSpace(level);
	//hyperMesh_.LabelBoundaryNode();

	hyperMesh_.mesh.update_normals();

	return true;
}



//-----------------------------------------------------------------------------


void
MeshViewer::update_face_indices()
{
	suMesh::ConstFaceIter        f_it(mesh_.faces_sbegin()),
		f_end(mesh_.faces_end());
	suMesh::ConstFaceVertexIter  fv_it;

	indices_.clear();
	indices_.reserve(mesh_.n_faces() * 3);

	for (; f_it != f_end; ++f_it)
	{
		indices_.push_back((fv_it = mesh_.cfv_iter(f_it)).handle().idx());
		indices_.push_back((++fv_it).handle().idx());
		indices_.push_back((++fv_it).handle().idx());
	}
}


//-----------------------------------------------------------------------------


void
MeshViewer::
draw(const std::string& _draw_mode)
{
	//glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

	if (indices_.empty())
	{
		GlutExaminer::draw(_draw_mode);
		return;
	}

	if (_draw_mode == "Wireframe")
	{
		glDisable(GL_LIGHTING);
		glColor3f(0.0, 0.0, 0.0);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

		glEnableClientState(GL_VERTEX_ARRAY);
		GL::glVertexPointer(mesh_.points());

		glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);

		glDisableClientState(GL_VERTEX_ARRAY);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}


	else if (_draw_mode == "Hidden Line")
	{
		glDisable(GL_LIGHTING);
		glColor3f(0.0, 0.0, 0.0);
		glEnableClientState(GL_VERTEX_ARRAY);
		GL::glVertexPointer(mesh_.points());

		glDrawBuffer(GL_NONE);
		glDepthRange(0.01, 1.0);
		glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);

		glDrawBuffer(GL_BACK);
		glDepthRange(0.0, 1.0);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDepthFunc(GL_LEQUAL);
		glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);

		glDisableClientState(GL_VERTEX_ARRAY);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDepthFunc(GL_LESS);
	}


	else if (_draw_mode == "Solid Flat")
	{
		suMesh::ConstFaceIter        f_it(mesh_.faces_begin()),
			f_end(mesh_.faces_end());
		suMesh::ConstFaceVertexIter  fv_it;

		glEnable(GL_LIGHTING);
		glShadeModel(GL_FLAT);

		glBegin(GL_TRIANGLES);
		for (; f_it != f_end; ++f_it)
		{
			GL::glNormal(mesh_.normal(f_it));
			fv_it = mesh_.cfv_iter(f_it.handle());
			GL::glVertex(mesh_.point(fv_it));
			++fv_it;
			GL::glVertex(mesh_.point(fv_it));
			++fv_it;
			GL::glVertex(mesh_.point(fv_it));
		}
		glEnd();
	}


	else if (_draw_mode == "Solid Smooth")
	{
		glEnable(GL_LIGHTING);
		glShadeModel(GL_SMOOTH);

		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);
		GL::glVertexPointer(mesh_.points());
		GL::glNormalPointer(mesh_.vertex_normals());

		glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
	}

	// draw octree
	//octreeForDraw_.draw();
	glColor3f(0.0, 0.0, 1.0);
	hyperMesh_.draw();

	TwDraw();
	glutPostRedisplay();
}

void MeshViewer::processmenu(int i)
{
	if (i == test_)
	{
		hyperMesh_.test();
	}
	if (i == test1_)
	{
		cout << "type the path\n";
		char str[256];
		char *p = str;
		cin >> str;
		if (p)
		{
			if (this->open_mesh(p))
			{
				glutPostRedisplay();
				std::cout << p << " is openned." << std::endl;
				num_grid = hyperMesh_.octree.leafNodesArr.size();
			}
		}
		cout << "type the octree level\n";
		int oct_level;
		cin >> oct_level;
		cout << oct_level;
		hyperMesh_.PartitionSpace(oct_level);

		glutPostRedisplay();
		std::cout << " octree level: " << oct_level << std::endl;
		
		hyperMesh_.LabelBoundaryNode();
		hyperMesh_.LabelBoundaryNeighbors();
		glutPostRedisplay();

		hyperMesh_.floodfill();
		glutPostRedisplay();
		hyperMesh_.test();
		//system("cd /d D:\\oofem\\build2.3\\Debug &oofem.exe&oofem -f test.txt");
		//cout << "oofem done\n";
		//read_point_information_("D:\\oofem\\build2.3\\Debug\\Majnun.out.m0.1.vtu");



	}
	if (i >= MENU_FILE_LOAD)
	{
		//load file				
		if (i == MENU_FILE_LOAD)
		{
			const char* const  extFilename[] = { "*.stl", "*.ply", "*.obj" };
			const char *p = tinyfd_openFileDialog("Open model", "c:\\", 3, extFilename, NULL, 0);
			//char str[256];
			//char *p = str;
			//cin >> str;
			if (p)
			{
				if (this->open_mesh(p))
				{
					glutPostRedisplay();
					std::cout << p << " is openned." << std::endl;
					num_grid = hyperMesh_.octree.leafNodesArr.size();
				}
			}
		}else if (i == MENU_FILE_SAVE && this->hyperMesh_.isLoad())
		{
			if (hyperMesh_.isLoad())
			{
				const char* const  extFilename[] = { "*.vtk" };
				const char *p = tinyfd_saveFileDialog("Save model", "c:\\", 1, extFilename, NULL);

				

				if (p)
				{
					//save to vtk format
					hyperMesh_.saveVTK(p);
				}
			}
			
		}

		//others
		else if (i == MENU_LABEL_BOUNDARY)
		{
			if (hyperMesh_.isLoad())
			{
				hyperMesh_.LabelBoundaryNode();
				hyperMesh_.LabelBoundaryNeighbors();
				glutPostRedisplay();
			}
		}
		else if (i == MENU_LABEL_BNEIGHBOR)
		{
			if (hyperMesh_.isLoad())
			{

				hyperMesh_.floodfill();
				glutPostRedisplay();
			}
		}
	}
	else{
		set_draw_mode(i);
	}

}

void MeshViewer::keyboard(int key, int x, int y)
{
	if (TwEventKeyboardGLUT(key, x, y))  // send event to AntTweakBar
	{
		return;
	}
	// event has not been handled by AntTweakBar
	switch (key)
	{
	case '+':
	{
				if (octreeForDraw_level_ < 8)
				{
					octreeForDraw_level_++;
					//octreeForDraw_.clear();
					//octreeForDraw_.Generate(octreeForDraw_level_);
					hyperMesh_.PartitionSpace(octreeForDraw_level_);
					num_grid = hyperMesh_.octree.leafNodesArr.size();
					TwRefreshBar(pBar);
					glutPostRedisplay();
					std::cout << " octree level: " << octreeForDraw_level_ << std::endl;
				}
				break;

	}
	case  '-':
	{
				 if (octreeForDraw_level_ > 0)
				 {
					 octreeForDraw_level_--;
					 /* octreeForDraw_.clear();
					  octreeForDraw_.Generate(octreeForDraw_level_);*/
					 hyperMesh_.PartitionSpace(octreeForDraw_level_);

					 glutPostRedisplay();
					 std::cout << " octree level: " << octreeForDraw_level_ << std::endl;
				 }
				 break;
	}

	default:
	{
			   GlutExaminer::keyboard(key, x, y);
			   break;
	}
	}

	TwEventKeyboardGLUT(key, x, y);

}

void MeshViewer::init()
{
	pBar = 0;

	//init
	show_grid = false;
	show_boundary = true;
	show_exterior = true;
	show_interior = true;

	TwInit(TW_OPENGL, 0);
	TwGLUTModifiersFunc(glutGetModifiers);
	pBar = TwNewBar("Parameters");
	TwDefine("Parameters size='220 280'");

	TwAddVarRW(pBar, "Lattic", TW_TYPE_BOOL8, &show_grid, "group=Rendering");	
	TwAddVarCB(pBar, "SHOW_BOUNDARY", TW_TYPE_BOOL8, setBool, getBool, &show_boundary, "group=Rendering");
	TwAddVarCB(pBar, "SHOW_EXTERIOR", TW_TYPE_BOOL8, setBool, getBool, &show_exterior, "group=Rendering");
	TwAddVarCB(pBar, "SHOW_INTERIOR", TW_TYPE_BOOL8, setBool, getBool, &show_interior, "group=Rendering");

	TwAddVarRW(pBar, "Level", TW_TYPE_INT32, &octreeForDraw_level_, "group=Grid");
	TwAddVarRO(pBar, "Voxels", TW_TYPE_INT32,&num_grid, "group=Grid");

	TwRefreshBar(pBar);
}

void MeshViewer::reshape(int _w, int _h)
{
	TwWindowSize(_w, _h);
	GlutExaminer::reshape(_w, _h);
}

void MeshViewer::release(void)
{
	hyperMesh_.clear();
	TwTerminate();
}

void MeshViewer::special(int key, int x, int y)
{

	switch (key) {
	case GLUT_KEY_PAGE_UP:
		//cameraHeight = min(800.0f, cameraHeight * 1.1f);
		TwRefreshBar(pBar);
		break;
	case GLUT_KEY_PAGE_DOWN:
		//cameraHeight = max(0.5f, cameraHeight / 1.1f);
		TwRefreshBar(pBar);
		break;
	}
	GlutViewer::special(key, x, y);
}

void MeshViewer::mouse(int button, int state, int x, int y)
{
	if (!TwEventMouseButtonGLUT(button, state, x, y) && state == 0) {

		GlutExaminer::mouse(button, state, x, y);
	}
}
void MeshViewer::motion(int x, int y)
{
	if (!TwMouseMotion(x, y))
	{
		GlutExaminer::motion(x, y);
	}
}

void MeshViewer::passivemotion(int x, int y)
{
	TwMouseMotion(x, y);
}



//=============================================================================
