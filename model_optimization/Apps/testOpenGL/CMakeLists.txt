cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

SET(APP_NAME OctreeViewer)

PROJECT(${APP_NAME})

###############Find Library#################
FIND_PACKAGE(OpenGL REQUIRED)
LIST(APPEND LIBS ${OPENGL_gl_LIBRARY})

##################Test Unit###################
## 头文件目录
SET(PROJECT_INCLUDE_DIRS       
     ####这里可以继续添加包含文件目录
	 ####例如 ${PROJECT_SOURCE_DIR}/GLTools/include
      ${PROJECT_SOURCE_DIR}/../../include
)

## 库目录
SET(PROJECT_LIBRARY_DIRS 
      ####这里可以继续添加库文件目录
	  ####例如：${QT_DIR}/../lib     
       ${PROJECT_SOURCE_DIR}/lib       
      )

##测试单元作为主程序（也可去掉）
FILE(GLOB MAIN_SRC_SOURCES  *.cpp *.h) 


## Construct libraries
INCLUDE_DIRECTORIES(  
  ${PROJECT_INCLUDE_DIRS}
)
LINK_DIRECTORIES(
  ${PROJECT_LIBRARY_DIRS}
)    
ADD_EXECUTABLE(${APP_NAME} ${MAIN_SRC_SOURCES})

##################MICRO: ADD_DEMO##################
SET(FREEGLUT_BUILD_STATIC_LIBS on)
SET(FREEGLUT_BUILD_SHARED_LIBS off)

MACRO(ADD_DEMO demoDir)
    FILE(GLOB DEMO_SRC_SOURCES   ${demoDir}/*.cpp ${demoDir}/*.tpp ${demoDir}/*.c) 
    FILE(GLOB DEMO_SRC_HEADERS   ${demoDir}/*.hpp ${demoDir}/*.h) 
    
    ##set(_args ${ARGN})  
    
    set(SRC_LIST ${DEMO_SRC_SOURCES} ${DEMO_SRC_HEADERS} )
     
    
    IF(FREEGLUT_BUILD_SHARED_LIBS)
                SET(DEMO_NAME ${demoDir})
                ADD_EXECUTABLE(${DEMO_NAME} ${SRC_LIST} )
                TARGET_LINK_LIBRARIES(${demoDir} ${LIBS} freeglut)
                IF(WIN32 AND MSVC)
                    SET_TARGET_PROPERTIES(${demoDir} PROPERTIES DEBUG_POSTFIX d)
                ENDIF()
                
                
            ENDIF()
            IF(FREEGLUT_BUILD_STATIC_LIBS)
                SET(DEMO_NAME ${demoDir}_static)
                ADD_EXECUTABLE(${DEMO_NAME} ${SRC_LIST} )
                TARGET_LINK_LIBRARIES(${demoDir}_static ${DEMO_LIBS} freeglut_static OpenGL32)
                SET_TARGET_PROPERTIES(${demoDir}_static PROPERTIES COMPILE_FLAGS -DFREEGLUT_STATIC)
                IF(WIN32 AND MSVC)
                    SET_TARGET_PROPERTIES(${demoDir}_static PROPERTIES DEBUG_POSTFIX d)
            ENDIF()
    ENDIF()
    
    ##Copy demo exe to  ${${APP_NAME}_SOURCE_DIR}
    ADD_CUSTOM_COMMAND(TARGET ${DEMO_NAME}
        		  POST_BUILD          
        	          COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${DEMO_NAME}> ${${APP_NAME}_SOURCE_DIR}/bin
    )    
ENDMACRO()


##################Demos###################
message(STATUS "添加测试应用OpenGlut demo：shapes")
ADD_DEMO(shapes)
message(STATUS "添加测试应用OpenGlut demo：octreeViewer")
message(STATUS "OpenGlut库静态链接到exe中，因此不需要dll库")
ADD_DEMO(octreeViewer)