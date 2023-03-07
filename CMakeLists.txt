cmake_minimum_required(VERSION 3.10)

project(minimap2 VERSION 2.22.0 LANGUAGES C)
enable_testing()

option(ONT_MM2_EXE "Build minimap2 exe" ON)

find_package(Git REQUIRED)
file(GLOB_RECURSE PATCHES *.patch)

foreach(PATCH ${PATCHES})
    message(STATUS "Applying patch ${PATCH}")
    exec_program("${GIT_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/minimap2" ARGS "apply ${PATCH}")
endforeach()

if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
  endif()
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

if(WIN32)

  set(ZLIB_COMPAT ON)
  set(ZLIB_ENABLE_TESTS OFF)
  add_subdirectory(zlib-ng) 
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/zlib-ng)
  set(ZLIB_LIBRARIES zlib)

  add_compile_definitions(PTW32_STATIC_LIB)
  add_compile_definitions(HAVE_STRUCT_TIMESPEC)
  add_subdirectory(pthread-win32)
  include_directories(${CMAKE_SOURCE_DIR}/pthread-win32)
  set(CMAKE_THREAD_LIBS_INIT ${CMAKE_CURRENT_BINARY_DIR}/pthread-win32/${CMAKE_BUILD_TYPE}/pthreadVC3.lib)

else()
  find_package(ZLIB REQUIRED)
  find_package(Threads REQUIRED)
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_compile_definitions(HAVE_KALLOC)
add_compile_definitions(__SSE2__)

if(CMAKE_SYSTEM_PROCESSOR MATCHES "^aarch64*|^arm*")
   include_directories(${CMAKE_SOURCE_DIR}/minimap2/sse2neon)
   add_compile_options(-fsigned-char)
   add_compile_definitions(_FILE_OFFSET_BITS=64)
endif()

add_library(mm2
    minimap2/align.c
    minimap2/bseq.c
    minimap2/esterr.c
    minimap2/format.c
    minimap2/hit.c
    minimap2/index.c
    minimap2/kalloc.c
    minimap2/ksw2_extd2_sse.c
    minimap2/ksw2_exts2_sse.c
    minimap2/ksw2_extz2_sse.c
    minimap2/ksw2_ll_sse.c
    minimap2/kthread.c
    minimap2/lchain.c
    minimap2/map.c
    minimap2/misc.c
    minimap2/options.c
    minimap2/pe.c
    minimap2/sdust.c
    minimap2/seed.c
    minimap2/sketch.c
    minimap2/splitidx.c
)

include_directories(${CMAKE_SOURCE_DIR}/minimap2) 
add_library(mm2cpp src/ont_minimap2.c)

target_link_libraries(mm2 ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

if(NOT WIN32)
  target_link_libraries(mm2 m)
endif()

target_link_libraries(mm2cpp mm2)

if(ONT_MM2_EXE)
  add_executable(minimap2 minimap2/main.c)
  target_link_libraries(minimap2 mm2)
  get_target_property(MM2_BIN_DIR minimap2 SOURCE_DIR)
  add_test(
    NAME minimap2_version_test
    WORKING_DIRECTORY ${MM2_BIN_DIR}
    COMMAND minimap2 -V
  )
endif()
