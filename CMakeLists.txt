cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(MODES LANGUAGES C CXX VERSION 0.0.1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_VERBOSE_MAKEFILE FALSE)

#if (CMAKE_COMPILER_IS_GNUCXX)
#    set(CMAKE_C_FLAGS_DEBUG "-g -ggdb -pg -DDEBUG")
#    set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb -pg -DDEBUG")
#    set(CMAKE_C_FLAGS_RELEASE "-O4 -DNDEBUG -march=native")
#    set(CMAKE_CXX_FLAGS_RELEASE "-O4 -DNDEBUG -march=native")
#endif ()

# Set CMake build-type. If it not supplied by the user, the default built type is "Release".
if ((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
    set (CMAKE_BUILD_TYPE Release)
    message(STATUS "Building release version!")
elseif (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set (CMAKE_BUILD_TYPE Debug)
    message (STATUS "WARNING: building debug version!")
endif()

#IF (WIN32)
#    SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
#ELSE ()
#    SET(CMAKE_FIND_LIBRARY_SUFFIXES .so .dylib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
#ENDIF ()

# set variables
set (PROJECT_PATH                                "${CMAKE_CURRENT_SOURCE_DIR}")
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}      "${PROJECT_PATH}/cmake/modules/")
set (INCLUDE_PATH                                "${PROJECT_PATH}/include/")
#set (LIB_PATH                                    "${PROJECT_BINARY_DIR}/lib")
set (EXAMPLES_PATH                               "${PROJECT_PATH}/src/")

# options
option(BUILD_STATIC         "Build static library"     OFF)
option(BUILD_EXAMPLES       "Build examples"           ON)
option(ENABLE_UQ            "Use smartuq"              OFF)
option(ENABLE_GSL           "Use GSL"                  OFF)
option(ENABLE_HDF5           "Use HDF5"                  OFF)
option(ENABLE_PAGMO           "Use pagmo"                  OFF)

#find_package(Threads REQUIRED)

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

# Include external libraries
set (MANDATORY_LIBRARIES "")
include(cmake/CMakeLists.txt)

# Include header directories
get_filename_component(HEADERS "include" REALPATH)
include_directories(${HEADERS})

# Include all cpp in src
file(GLOB_RECURSE SOURCES "libraries/*.cpp")

#add library
add_library (MODES SHARED ${SOURCES})
target_link_libraries(MODES ${MANDATORY_LIBRARIES})
set_target_properties(MODES PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION})

if(BUILD_STATIC)
    add_library (MODES-static STATIC ${SOURCES})
    target_link_libraries(MODES-static ${MANDATORY_LIBRARIES})
    set_target_properties(MODES-static PROPERTIES
            VERSION ${PROJECT_VERSION}
            SOVERSION ${PROJECT_VERSION}
            OUTPUT_NAME MODES)
    add_dependencies(MODES-static MODES)
endif()

if ( BUILD_EXAMPLES )
    add_subdirectory ("src")
endif()

# Add main as execturable
add_executable(MODES-main src/main.cpp ${HEADERS})
target_link_libraries(MODES-main MODES ${MANDATORY_LIBRARIES})
set_target_properties(MODES-main PROPERTIES LINKER_LANGUAGE CXX )
add_dependencies(MODES-main MODES)

# user.
install (DIRECTORY "${INCLUDE_PATH}/" DESTINATION include/${PROJECT_NAME} FILES_MATCHING PATTERN "*.h")
if (BUILD_STATIC)
    install (TARGETS MODES-static DESTINATION lib/${CMAKE_PROJECT_NAME})
endif()
install (TARGETS MODES DESTINATION lib/${CMAKE_PROJECT_NAME})
install (TARGETS ${BIN_NAME} DESTINATION bin/${CMAKE_PROJECT_NAME})
