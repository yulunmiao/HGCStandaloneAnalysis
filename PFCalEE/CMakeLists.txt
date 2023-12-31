
#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(HGCALSS)

#set(CMAKE_MODULE_PATH ${BASEINSTALL}/cmake/modules/FindHepMC.cmake)
set(HepMC_DIR /cvmfs/sft.cern.ch/lcg/releases/HepMC/2.06.10-1a364/x86_64-centos7-gcc8-opt/)
set(HepMC_FOUND 1)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#find_package(HepMC REQUIRED)
#if(HEPMC_FOUND)
#add_subdirectory(HepMC)
#else()
#  message(STATUS "G4 Examples: HepMC package not found.") 
#endif()

#debug symbols
set(CMAKE_BUILD_TYPE Debug)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${HepMC_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/userlib/include)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(PFCalEE PFCalEE.cc ${sources} ${headers})
target_link_libraries(PFCalEE ${Geant4_LIBRARIES} ${HepMC_DIR}/lib/libHepMCfio.so ${HepMC_DIR}/lib/libHepMC.so ${PROJECT_SOURCE_DIR}/userlib/lib/libPFCalEEuserlib.so)
#target_link_libraries(PFCalEE ${Geant4_LIBRARIES} ${PROJECT_SOURCE_DIR}/userlib/lib/libPFCalEEuserlib.so)

#add_executable(TestMaterial TestMaterial.cc ${sources} ${headers})
#target_link_libraries(TestMaterial ${Geant4_LIBRARIES} ${HepMC_DIR}/lib/libHepMCfio.so ${HepMC_DIR}/lib/libHepMC.so ${PROJECT_SOURCE_DIR}/userlib/lib/libPFCalEEuserlib.so)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
#set(EXAMPLEB1_SCRIPTS
#  exampleB1.in
#  exampleB1.out
#  init_vis.mac
#  run1.mac
#  run2.mac
#  vis.mac
#  )

#foreach(_script ${EXAMPLEB1_SCRIPTS})
#  configure_file(
#    ${PROJECT_SOURCE_DIR}/${_script}
#    ${PROJECT_BINARY_DIR}/${_script}
#    COPYONLY
#    )
#endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(HGCALSS DEPENDS PFCalEE)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS PFCalEE DESTINATION bin)


