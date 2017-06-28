cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

Project(cal_DEM)

find_package(OpenCAL 1.0 COMPONENTS  OpenCAL  REQUIRED)

OPTION (GUI "GUI" ON)

if(GUI)
    find_package(OpenGL REQUIRED)
    find_package(GLUT REQUIRED)
endif()

if(MSVC)
  set(M_LIB "") # not needed on Microsoft Visual C++
else()
    find_library(M_LIB m)
endif()


file(GLOB HDR
    ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h
    )

file(GLOB SRC
  ${CMAKE_CURRENT_SOURCE_DIR}/source/*.c
  )

file (GLOB CLI_SRC
  ${CMAKE_CURRENT_SOURCE_DIR}/source_cli/*.c
  )

file (GLOB GUI_SRC
  ${CMAKE_CURRENT_SOURCE_DIR}/source_gui/*.c
  )

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${OPENCAL_INCLUDE_DIR}
  ${OPENGL_INCLUDE_DIRS}
  ${GLUT_INCLUDE_DIR}
  )

if(GUI)
include_directories(
  ${OPENGL_INCLUDE_DIRS}
  ${GLUT_INCLUDE_DIR}
  )
endif()


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})


if(GUI)
    message(STATUS "Building with GUI support...")
    add_executable(${PROJECT_NAME} ${GUI_SRC} ${SRC} ${HDR})
    target_link_libraries(${PROJECT_NAME} ${OPENCAL_LIBRARIES} ${M_LIB} ${OPENGL_LIBRARIES} ${GLUT_glut_LIBRARY})
else() # CLI
    message(STATUS "Building with CLI support...")
    add_executable(${PROJECT_NAME} ${CLI_SRC} ${SRC} ${HDR})
    target_link_libraries(${PROJECT_NAME} ${OPENCAL_LIBRARIES} ${M_LIB})
endif()
