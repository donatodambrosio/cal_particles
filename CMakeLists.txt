cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

Project(cal_particelle_giuseppi)

find_package(OpenCAL 1.0 COMPONENTS  OpenCAL  REQUIRED)
#find_package(Boost COMPONENTS date_time REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

if(NOT MSVC)
  find_library(M_LIB m)
else()
  # not needed on windows
  set(M_LIB "")
endif()

file(GLOB SRC
  ${CMAKE_CURRENT_SOURCE_DIR}/source/*.c
  ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp
  )

file(GLOB HDR
  ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h
  )

link_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/source
  )
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include
#  ${Boost_INCLUDE_DIRS}
  ${OPENGL_INCLUDE_DIRS}
  ${GLUT_INCLUDE_DIR}
  ${OPENCAL_INCLUDE_DIR}
  )


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(${PROJECT_NAME} ${SRC} ${HDR})
#target_link_libraries(${PROJECT_NAME} ${OPENCAL_LIBRARIES} ${Boost_LIBRARIES} ${M_LIB} ${OPENGL_LIBRARIES} ${GLUT_glut_LIBRARY})
target_link_libraries(${PROJECT_NAME} ${OPENCAL_LIBRARIES} ${M_LIB} ${OPENGL_LIBRARIES} ${GLUT_glut_LIBRARY})
