cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

Project(cal_DEM)

OPTION (GUI "GUI" OFF)
OPTION (OMP "OpenMP" OFF)
OPTION (VERBOSE "verbose mode" ON)

if(OMP)
    add_definitions(-DOMP=ON)
    find_package(OpenCAL 1.0 COMPONENTS  OpenCAL-OMP REQUIRED)
    FIND_PACKAGE(OpenMP REQUIRED)
    IF(OPENMP_FOUND)
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    ENDIF()
else()
    find_package(OpenCAL 1.0 COMPONENTS  OpenCAL  REQUIRED)
endif()

if(VERBOSE)
    add_definitions(-DVERBOSE=ON)
endif()

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

file(GLOB GUI_HDR
    ${CMAKE_CURRENT_SOURCE_DIR}/include_gui/*.h
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
  ${CMAKE_CURRENT_SOURCE_DIR}/include_gui
  ${OPENGL_INCLUDE_DIRS}
  ${GLUT_INCLUDE_DIR}
  )
endif()


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})


if(GUI)
    message(STATUS "Building with GUI support...")
    add_executable(${PROJECT_NAME} ${GUI_SRC} ${GUI_HDR} ${SRC} ${HDR})
    if(OMP)
        target_link_libraries(${PROJECT_NAME} ${OPENCAL_LIBRARIES} ${M_LIB} ${OPENGL_LIBRARIES} ${GLUT_glut_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
    else()
        target_link_libraries(${PROJECT_NAME} ${OPENCAL_LIBRARIES} ${M_LIB} ${OPENGL_LIBRARIES} ${GLUT_glut_LIBRARY})
    endif()
else() # CLI
    message(STATUS "Building with CLI support...")
    add_executable(${PROJECT_NAME} ${CLI_SRC} ${SRC} ${HDR})
    if(OMP)
        target_link_libraries(${PROJECT_NAME} ${OPENCAL_LIBRARIES} ${M_LIB} ${CMAKE_THREAD_LIBS_INIT})
    else()
        target_link_libraries(${PROJECT_NAME} ${OPENCAL_LIBRARIES} ${M_LIB})
    endif()
endif()
