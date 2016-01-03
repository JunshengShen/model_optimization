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
//                                                                            
//=============================================================================
//=============================================================================
//  CLASS MeshViewerWidget
//
//=============================================================================


#ifndef MESH_VIEWER_WIDGET_HH
#define MESH_VIEWER_WIDGET_HH


//== INCLUDES =================================================================
//#include "../octreeViewer/suMesh.h"
//#include "suObjOctree.h"
#include "suHyperMesh.h"
#include "GlutExaminer.h"
#include <GL/AntTweakBar.h>



//== CLASS DEFINITION =========================================================


class MeshViewer : public GlutExaminer
{
public:
   
  /// default constructor
  MeshViewer(const char* _title, int _width, int _height);

  /// open mesh
  virtual bool open_mesh(const char* _filename);

  bool open_hyperMesh(const char* _filename);

  /// update buffer with face indices
  void update_face_indices();

  /// draw the scene
  virtual void draw(const std::string& _draw_mode);

  /// process menu
  virtual void processmenu(int i);

  /// process event
  virtual void keyboard(int key, int x, int y);
  virtual void special(int key, int x, int y);
  /// mouse click
  virtual void mouse(int button, int state, int x, int y);
  /// mouse motion
  virtual void motion(int x, int y);  
  virtual void passivemotion(int x, int y);

  /// init twbar
  virtual void init();

  /// reshape for twbar
  virtual void reshape(int _w, int _h);

  /// release
  virtual void release(void);

public:
  suHyperMesh                  hyperMesh_;
  suMesh                       mesh_;          //link to hyperMesh.mesh  
  std::vector<unsigned int>    indices_;


private:
	int octreeForDraw_level_;

	//bar
	TwBar *pBar;
	bool show_grid;	
	bool show_boundary ;
	bool show_exterior ;
	bool show_interior ;
	unsigned int num_grid;
};


//=============================================================================
#endif // MESH_VIEWER_WIDGET_HH defined
//=============================================================================

