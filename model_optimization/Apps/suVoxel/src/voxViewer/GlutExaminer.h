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
//  CLASS GlutExaminer
//
//=============================================================================


#ifndef GLUTEXAMINER_HH
#define GLUTEXAMINER_HH


//== INCLUDES =================================================================

#include "GlutViewer.h"
#include <OpenMesh/Core/Geometry/VectorT.hh>

#include <string>
#include <vector>


//== CLASS DEFINITION =========================================================



class GlutExaminer : public GlutViewer
{
public:

  GlutExaminer(const char* _title, int _width, int _height);

  void   set_scene(const OpenMesh::Vec3f& _center, float _radius);
  void   view_all();
  double measure_fps();


protected:

  virtual void init();
  virtual void draw(const std::string& _draw_mode);


  // overloaded glut functions
  virtual void motion(int x, int y);
  virtual void mouse(int button, int state, int x, int y);
  virtual void reshape(int w, int h);
  virtual void keyboard(int key, int x, int y);


  // updates projection matrix
  void update_projection_matrix();
  // translate the scene and update modelview matrix
  void translate(const OpenMesh::Vec3f& _trans);
  // rotate the scene (around its center) and update modelview matrix
  void rotate(const OpenMesh::Vec3f& _axis, float _angle);


  // virtual trackball: map 2D screen point to unit sphere
  bool map_to_sphere(const OpenMesh::Vec2i& _point, OpenMesh::Vec3f& _result);


  // mouse processing functions
  void rotation(int x, int y);
  void translation(int x, int y);
  void zoom(int x, int y);


protected:

  // scene position and dimension
  OpenMesh::Vec3f    center_;
  float    radius_;


  // projection parameters
  float    near_, far_, fovy_;


  // OpenGL matrices
  double   projection_matrix_[16],
           modelview_matrix_[16];


  // trackball helpers
  OpenMesh::Vec2i    last_point_2D_;
  OpenMesh::Vec3f    last_point_3D_;
  bool     last_point_ok_;
  bool     button_down_[10];
};


//=============================================================================
#endif // GLUTEXAMINER_HH defined
//=============================================================================

