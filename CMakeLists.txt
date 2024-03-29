cmake_minimum_required(VERSION 3.10)

project(minimap2 VERSION 2.24.0 LANGUAGES C)
enable_testing()

option(ONT_MM2_EXE "Build minimap2 exe" ON)

# set default build type to release
if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
  endif()
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# apply the required cross platform patches to minimap2
find_package(Git REQUIRED)
file(GLOB_RECURSE PATCHES *.patch)

foreach(PATCH ${PATCHES})
    # See if we need to apply the patch.
    execute_process(
        COMMAND ${GIT_EXECUTABLE} apply --reverse --check ${PATCH}
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/src/3rdparty/minimap2"
        RESULT_VARIABLE NOT_ALREADY_APPLIED
        ERROR_QUIET
    )
    if (NOT_ALREADY_APPLIED)
        message(STATUS "Applying patch ${PATCH}")
        execute_process(
            COMMAND ${GIT_EXECUTABLE} apply ${PATCH}
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/src/3rdparty/minimap2"
            RESULT_VARIABLE PATCH_FAILED
        )
        if (PATCH_FAILED)
            message(FATAL_ERROR "Failed to apply patch ${PATCH}")
        endif()
    endif()
endforeach()

# handle zlib and pthreads on windows
if(WIN32)
  set(ZLIB_COMPAT ON)
  set(ZLIB_ENABLE_TESTS OFF)
  add_subdirectory(src/3rdparty/zlib-ng)
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/zlib-ng)
  set(ZLIB_LIBRARIES zlib)

  add_compile_definitions(PTW32_STATIC_LIB)
  add_compile_definitions(HAVE_STRUCT_TIMESPEC)
  add_subdirectory(src/3rdparty/pthread-win32)
  include_directories(src/3rdparty/pthread-win32)
else()
  find_package(ZLIB REQUIRED)
  find_package(Threads REQUIRED)
endif()

# compiler settings
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_compile_definitions(HAVE_KALLOC)
add_compile_definitions(__SSE2__)

if(CMAKE_SYSTEM_PROCESSOR MATCHES "^aarch64*|^arm*")
   include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2/sse2neon)
   add_compile_options(-fsigned-char)
   add_compile_definitions(_FILE_OFFSET_BITS=64)
endif()

set(MM2_SOURCES
    src/3rdparty/minimap2/align.c
    src/3rdparty/minimap2/bseq.c
    src/3rdparty/minimap2/esterr.c
    src/3rdparty/minimap2/format.c
    src/3rdparty/minimap2/hit.c
    src/3rdparty/minimap2/index.c
    src/3rdparty/minimap2/kalloc.c
    src/3rdparty/minimap2/ksw2_extd2_sse.c
    src/3rdparty/minimap2/ksw2_exts2_sse.c
    src/3rdparty/minimap2/ksw2_extz2_sse.c
    src/3rdparty/minimap2/ksw2_ll_sse.c
    src/3rdparty/minimap2/kthread.c
    src/3rdparty/minimap2/lchain.c
    src/3rdparty/minimap2/map.c
    src/3rdparty/minimap2/misc.c
    src/3rdparty/minimap2/options.c
    src/3rdparty/minimap2/pe.c
    src/3rdparty/minimap2/sdust.c
    src/3rdparty/minimap2/seed.c
    src/3rdparty/minimap2/sketch.c
    src/3rdparty/minimap2/splitidx.c
)

# Create the minimap2 library
add_library(minimap2 ${MM2_SOURCES})
target_include_directories(minimap2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2)
target_link_libraries(minimap2 ${ZLIB_LIBRARIES})
if(NOT WIN32)
  target_link_libraries(minimap2 m Threads::Threads)
else()
  target_link_libraries(minimap2 pthreadVC3)
endif()
list(APPEND TARGETS minimap2)

# Optionally create the executable
if(ONT_MM2_EXE)
  add_executable(minimap2_exe src/3rdparty/minimap2/main.c)
  target_link_libraries(minimap2_exe minimap2)
  list(APPEND TARGETS minimap2_exe)
  set_target_properties(minimap2_exe PROPERTIES OUTPUT_NAME minimap2)
  add_test(
    NAME minimap2_version_test
    COMMAND minimap2_exe -V
  )
endif()

install(TARGETS ${TARGETS})
