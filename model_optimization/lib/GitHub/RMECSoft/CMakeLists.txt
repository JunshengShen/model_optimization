cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
## π”√OpenMesh 3.2

project(RMECSoft)
set(APP_NAME RMECSoft)

###Include directory and Lib directory
set(PROJECT_INCLUDE_DIRS 
         ${${APP_NAME}_SOURCE_DIR}/include
         ${${APP_NAME}_SOURCE_DIR}/include/ThirdParty	
         )
set(PROJECT_LIBRARY_DIRS  ${${APP_NAME}_SOURCE_DIR}/lib)


include_directories(
  ${PROJECT_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}  
  )
  
link_directories( 
  ${PROJECT_LIBRARY_DIRS}     
  )

#### Lib
##MeshQuery Lib
FILE(GLOB SRC_MESHQUERY  Apps/MeshQuery/*.h Apps/MeshQuery/*.cpp Apps/MeshQuery/*.cxx)
add_library(MeshQuery ${SRC_MESHQUERY})


######################################################
## Voxilization program
FILE(GLOB VOX_SOURCES	Apps/voxlize/*.cpp Apps/voxlize/*.c) 
FILE(GLOB VOX_HEADERS	Apps/voxlize/*.h Apps/voxlize/*.hpp)

add_executable (voxel ${VOX_SOURCES} ${VOX_HEADERS})
target_link_libraries (voxel        
         OpenMeshCored
         OpenMeshToolsd		 
	 MeshQuery
)

######################################################
## Voxilization program
FILE(GLOB TEST_SOURCES	Apps/testUnit/*.cpp Apps/testUnit/*.c) 
FILE(GLOB TEST_HEADERS	Apps/testUnit/*.h Apps/testUnit/*.hpp)

add_executable (testUnit ${TEST_SOURCES} ${TEST_HEADERS})
target_link_libraries (testUnit        
         OpenMeshCored
         OpenMeshToolsd		 
	 MeshQuery
)

