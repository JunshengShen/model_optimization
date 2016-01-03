cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

SET(APP_NAME suVoxel)
PROJECT(${APP_NAME})

# ========================================================================
# Common Setting
# ========================================================================
## Add our macro directory to cmake search path
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set (CMAKE_DEBUG_POSTFIX "d")

## Include our cmake files
include (SuCommon)  ##include SuCommon.cmake

## 头文件目录
SET(PROJECT_INCLUDE_DIRS       
     ####这里可以继续添加包含文件目录
	 ####例如 ${PROJECT_SOURCE_DIR}/GLTools/include
        ${PROJECT_SOURCE_DIR}/../../include	  
        ${PROJECT_SOURCE_DIR}/../../include/ThirdParty      
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/../../lib/GitHub/RMECSoft/include/Engines/Core
)

## 库目录
SET(PROJECT_LIBRARY_DIRS 
      ####这里可以继续添加库文件目录
      ####例如：${QT_DIR}/../lib     
       ${PROJECT_SOURCE_DIR}/lib
       ${PROJECT_SOURCE_DIR}/../../lib
       
)



# ========================================================================
# Definitions
# ========================================================================
if (WIN32)
  ##for openMesh
  add_definitions(
      -D_USE_MATH_DEFINES -DNOMINMAX
	  -D_NOMINMAX
      -D_CRT_SECURE_NO_WARNINGS
  )
  ##for glut  
  SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /SAFESEH:NO") 

endif ()

# ========================================================================
# Libs
# ========================================================================
##MeshQuery Lib
FILE(GLOB SRC_MESHQUERY  src/MeshQuery/*.h src/MeshQuery/*.cpp src/MeshQuery/*.cxx)
add_library(MeshQuery ${SRC_MESHQUERY})

##Find Library
FIND_PACKAGE(OpenGL REQUIRED)
LIST(APPEND LIBS ${OPENGL_gl_LIBRARY})

LIST(APPEND LIBS AntTweakBar.lib)


# ========================================================================
# Main Program Settings
# ========================================================================
## Construct libraries
INCLUDE_DIRECTORIES(  
  ${PROJECT_INCLUDE_DIRS}
)
LINK_DIRECTORIES(
  ${PROJECT_LIBRARY_DIRS}
)    


# ========================================================================
# Apps
# ========================================================================
## Test Unit
FILE(GLOB SRC_TESTUNIT  src/*.h src/*.cpp)
SU_ADD_EXECUTABLE(testVoxel src/testUnit/main.cpp ${SRC_TESTUNIT})
target_link_libraries (testVoxel ${LIBS})

## Voxel Viewer
FILE(GLOB SRC_VOXEL_VIEWER  src/voxViewer/*.h src/voxViewer/*.cpp)
SU_ADD_EXECUTABLE(VViewer ${SRC_VOXEL_VIEWER})
target_link_libraries (VViewer ${LIBS} AntTweakBar.lib OpenGL32 glu32 glut32 glew32 MeshQuery)

# ========================================================================
# Comments
# ========================================================================
message(STATUS "##################Todo#####################")
message(STATUS ". Mesh错误检测")
message(STATUS ". 单元选择")
message(STATUS ". Mesh 修复")