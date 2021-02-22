cmake_minimum_required(VERSION 3.12.0)

project(PIFO2 VERSION 0.1.0 LANGUAGES CXX)

# Setup CMake behaviour
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY    ${CMAKE_BINARY_DIR}/bin )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY    ${CMAKE_BINARY_DIR}/lib )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY    ${CMAKE_BINARY_DIR}/lib )
set( CMAKE_INSTALL_RPATH               ${CMAKE_INSTALL_PREFIX}/lib )
set( CMAKE_INSTALL_RPATH_USE_LINK_PATH True )
set( CMAKE_CXX_STANDARD 14)
set( CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-O3 -m64 -msse3 -g") # --save-temps -fverbose-asm
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -m64 -msse3")


#############################################
### Optionally add dependencies as subprojects

option( BUNDLE "Bundle dependencies as subprojects" OFF)
  # Can be enabled with `-DBUNDLE=ON`

if( BUNDLE )
    # Note that this case means that the source codes for
    # each of the following dependencies is required to
    # be available inside this directory, which may be
    # achieved via e.g. `copy`, `symbolic links`,
    #   `git submodule`, ...

    # Set features for Atlas and dependencies
    set( ENABLE_MPI   ON  CACHE BOOL "" FORCE )
    set( ENABLE_TESTS OFF CACHE BOOL "" FORCE )

    # Add Atlas and dependencies as subprojects
    add_subdirectory( ecbuild )
    add_subdirectory( eckit )
    add_subdirectory( fckit )
    add_subdirectory( atlas )

endif()

#############################################
### Main project build

find_package( atlas REQUIRED)
find_package( eckit REQUIRED)
find_package( eccodes REQUIRED)

include_directories(AFTER SYSTEM src)

add_subdirectory(src)


#############################################
### Installation

install( TARGETS pifo2
    RUNTIME DESTINATION bin )

#############################################
### Testing

find_package(Boost 1.70 REQUIRED) 

set(TEST_FILES test/test.cpp) 

include(CTest)
enable_testing()
add_executable(pifo_test ${TEST_FILES})
target_link_libraries(pifo_test PUBLIC atlas eckit eccodes app model util boost_unit_test_framework)
add_test(NAME PifoTest COMMAND $<TARGET_FILE:pifo_test>)

# The following is supposed to work but is not...
#find_package(Boost 1.71 REQUIRED COMPONENTS unit_test_framework)
#target_link_libraries(pifo_test PUBLIC atlas eckit eccodes app model util ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
#target_include_directories(pifo_test PRIVATE ${BOOST_INCLUDE_DIRS})



