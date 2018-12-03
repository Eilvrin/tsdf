cmake_minimum_required( VERSION 3.1 FATAL_ERROR)
project(tsdf)

set(CMAKE_CXX_STANDARD 14)

find_package(PCL REQUIRED)
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

include_directories(include ${PCL_INCLUDE_DIRS})

# Download and unpack googletest at configure time
configure_file(cmake/CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
                 ${CMAKE_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()

file(GLOB HDR_FILES include/${PROJECT_NAME}/*.h)
set(SRC_FILES src/main.cpp
              src/marchingCubesTSDF.cpp
              src/TSDF.cpp
              )

add_executable(${PROJECT_NAME} ${SRC_FILES} ${HDR_FILES})
target_link_libraries(${PROJECT_NAME} ${PCL_LIBRARIES})

# Now simply link against gtest or gtest_main as needed. Eg

add_executable(test_TSDFGrid test/test_TSDFGrid.cpp)
target_link_libraries(test_TSDFGrid gtest_main)
add_test(NAME GridMap_test COMMAND test_TSDFGrid)
