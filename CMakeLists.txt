#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16)
project(OpNovice2)
find_package(Geant4 REQUIRED)

## Get the ROOT config
include(FindROOT.cmake)
set(INC_ROOT ${ROOT_INCLUDE_DIR})
set(LIB_ROOT ${ROOT_LIBRARY_DIR})

add_definitions(${Geant4_CXX_FLAGS})
set(CMAKE_CXX_FLAGS ${Geant4_CXX_FLAGS})

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(GEANT4_INSTALL_DATA "Data" ON)
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include ${INC_ROOT})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(OpNovice2 OpNovice2.cc ${sources} ${headers})
target_link_libraries(OpNovice2 ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build OpNovice2. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(OpNovice2_SCRIPTS
    OpNovice2.out
    OpNovice2.mac
    boundary.mac
    electron.mac
    fresnel.mac
    scint_by_particle.mac
    vis.mac
    wls.mac
  )

foreach(_script ${OpNovice2_SCRIPTS})
    message(STATUS "Source file: ${CMAKE_CURRENT_SOURCE_DIR}/macros/${_script}")
    message(STATUS "Destination file: ${CMAKE_CURRENT_BINARY_DIR}/${_script}")
    message(STATUS "Configuring file: ${_script}")
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/macros/${_script}
        ${CMAKE_CURRENT_BINARY_DIR}/${_script}
        COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS OpNovice2 DESTINATION bin)
#
#----------------------------------------------------------------------------
ADD_CUSTOM_TARGET(link_target ALL
                  COMMAND ${CMAKE_COMMAND} -E create_symlink ../macros ./macros)

