cmake_minimum_required(VERSION 3.10)

project(minimap2 VERSION 2.24.0 LANGUAGES C)
enable_testing()

if (NOT WIN32)
  add_compile_options(
    # Match the warning options in the Makefile.
    -Wall
    # Always enable assert() macros.
    -UNDEBUG
  )
endif()

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
  find_package(ZLIB QUIET)
  if(NOT TARGET ZLIB::ZLIB)
    set(ZLIB_COMPAT ON)
    set(ZLIB_ENABLE_TESTS OFF)
    add_subdirectory(src/3rdparty/zlib-ng)
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/zlib-ng)
    set(ZLIB_LIBRARIES zlib)
  else()
      set(ZLIB_LIBRARIES ZLIB::ZLIB)
  endif()  
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
add_compile_definitions(
  HAVE_KALLOC
  __SSE2__ # The SSE implementations require this to be set, even for NEON.
  _FILE_OFFSET_BITS=64
)

# KSW is built as both SSE2+4.1 if supported.
if (CMAKE_SYSTEM_PROCESSOR MATCHES "^aarch64*|^arm*")
  # Non-SSE impl.
  add_library(minimap2_ksw2
    src/3rdparty/minimap2/ksw2_extd2_sse.c
    src/3rdparty/minimap2/ksw2_exts2_sse.c
    src/3rdparty/minimap2/ksw2_extz2_sse.c
    src/3rdparty/minimap2/ksw2_ll_sse.c
  )
  target_include_directories(minimap2_ksw2 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2/sse2neon)

else()
  # Dispatch impl.
  function(ont_minimap2_add_ksw2 NAME FLAGS)
    add_library(${NAME}
      src/3rdparty/minimap2/ksw2_extd2_sse.c
      src/3rdparty/minimap2/ksw2_exts2_sse.c
      src/3rdparty/minimap2/ksw2_extz2_sse.c
    )
    target_compile_definitions(${NAME} PRIVATE KSW_CPU_DISPATCH)
    target_compile_options(${NAME} PRIVATE ${${FLAGS}})
  endfunction()

  if (WIN32)
    # Windows gets these by default, but needs the define set to match the code.
    set(FLAGS_SSE2)
    set(FLAGS_SSE41 /D__SSE4_1__)
  else()
    # We build for x64 so the SSE2 flag isn't technically needed.
    set(FLAGS_SSE2 -msse2 -mno-sse4.1)
    set(FLAGS_SSE41 -msse4.1)
  endif()

  # Create both targets and link them into one.
  ont_minimap2_add_ksw2(minimap2_ksw2_sse2 FLAGS_SSE2)
  ont_minimap2_add_ksw2(minimap2_ksw2_sse41 FLAGS_SSE41)
  add_library(minimap2_ksw2
    src/3rdparty/minimap2/ksw2_dispatch.c
    src/3rdparty/minimap2/ksw2_ll_sse.c
  )
  target_compile_definitions(minimap2_ksw2 PRIVATE KSW_CPU_DISPATCH)
  target_link_libraries(minimap2_ksw2 PUBLIC minimap2_ksw2_sse2 minimap2_ksw2_sse41)
endif()

# Create the minimap2 library
add_library(minimap2
    src/3rdparty/minimap2/align.c
    src/3rdparty/minimap2/bseq.c
    src/3rdparty/minimap2/esterr.c
    src/3rdparty/minimap2/format.c
    src/3rdparty/minimap2/hit.c
    src/3rdparty/minimap2/index.c
    src/3rdparty/minimap2/jump.c
    src/3rdparty/minimap2/kalloc.c
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
target_include_directories(minimap2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2)
target_link_libraries(minimap2
  PUBLIC
    ${ZLIB_LIBRARIES}
  PRIVATE
    minimap2_ksw2
)
if(NOT WIN32)
  target_link_libraries(minimap2 PRIVATE m Threads::Threads)
else()
  target_link_libraries(minimap2 PRIVATE pthreadVC3)
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
