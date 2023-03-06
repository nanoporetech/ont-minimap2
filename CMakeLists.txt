cmake_minimum_required(VERSION 3.10)

project(minimap2 VERSION 2.24.0 LANGUAGES C)

option(ONT_MM2_EXE "Build minimap2 exe" OFF)

find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
        
  if(GIT_SUBMODULE)
    message(STATUS "Submodule update")
    
    file(LOCK ${CMAKE_SOURCE_DIR} DIRECTORY)
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT
    )
    file(LOCK ${CMAKE_SOURCE_DIR} DIRECTORY RELEASE)
    
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()
  endif()

endif()

if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
  endif()
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_compile_definitions(HAVE_KALLOC)

find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

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

target_link_libraries(mm2 PRIVATE ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} m)
target_link_libraries(mm2cpp PRIVATE mm2)

if (ONT_MM2_EXE) 
  add_executable(minimap2 minimap2/main.c)
  target_link_libraries(minimap2 PRIVATE mm2)
endif()
