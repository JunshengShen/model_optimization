cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(APP_NAME gicp)
project(${APP_NAME})


#just to avoid the warning
if(COMMAND cmake_policy)
     cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

FILE(GLOB SRC_SOURCES	*.cpp *.c) 
FILE(GLOB SRC_HEADERS	*.h *.hpp)

file(GLOB to_remove scan2ascii.cpp)
list(REMOVE_ITEM SRC_SOURCES ${to_remove})



#Append path to the module path
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

set(PROJECT_INCLUDE_DIRS ${${APP_NAME}_SOURCE_DIR}/../../include)
set(PROJECT_LIBRARY_DIRS ${${APP_NAME}_SOURCE_DIR}/../../lib)


#add boost package
SET(Boost_INCLUDE_DIRS ${CMAKE_INCLUDE_PATH} "E:/DevLib/boost/")
SET(Boost_LIBRARY_DIRS "E:/DevLib/boost/lib32-msvc-12.0")


include_directories(  
  ${Boost_INCLUDE_DIRS}
  ${PROJECT_INCLUDE_DIRS})
  
link_directories(
  ${Boost_LIBRARY_DIRS}  
  ${PROJECT_LIBRARY_DIRS})

add_executable (${APP_NAME}  ${SRC_SOURCES} ${SRC_HEADERS})

target_link_libraries (${APP_NAME}                        
                       ${Boost_LIBRARIES}
                       gsl.lib
                       ann.lib
                       cblas.lib
                       )
                       


FILE(GLOB SRC_SOURCES	*.cpp *.c) 
FILE(GLOB SRC_HEADERS	*.h *.hpp)

file(GLOB to_remove main.cpp)
list(REMOVE_ITEM SRC_SOURCES ${to_remove})


add_executable(scan2ascii ${SRC_SOURCES} ${SRC_HEADERS})
target_link_libraries (scan2ascii 
                       ${Boost_LIBRARIES}
                       gsl.lib
                       ann.lib
                       cblas.lib 
                       )


#####Set the default path for built executables to the "bin" directory
##set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../../bin)

# copy binary to running director
ADD_CUSTOM_COMMAND(TARGET ${APP_NAME} 
          POST_BUILD          
          COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${APP_NAME}> ${PROJECT_SOURCE_DIR}/../../bin
)

message(STATUS "Output exe file in " ${EXECUTABLE_OUTPUT_PATH})


#set optimization level 
set(CMAKE_BUILD_TYPE Release)

