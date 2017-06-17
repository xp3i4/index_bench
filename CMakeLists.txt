# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/index_bench/CMakeLists.txt
#
# CMakeLists.txt file for INDEX_BENCH.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_apps_index_bench)
message (STATUS "Configuring apps/index_bench")

set (SEQAN_APP_VERSION "1.1.1")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------
#find_package(CXX11 REQUIRED)
#if (NOT CXX11_FOUND)
#    message (STATUS " C++11 compiler needed to build index_bench")
#    return ()
#endif (NOT CXX11_FOUND)

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

# Update the list of file names below if you add source files to your application.
add_executable (index_bench index_bench.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (index_bench ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${CXX11_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (index_bench)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install index_bench in ${PREFIX}/bin directory
install (TARGETS index_bench
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/index_bench
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
install (FILES example/small.fasta
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (index_bench)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable seqan_tcoffee in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} index_bench CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:index_bench")
  set (CPACK_PACKAGE_NAME "index_bench")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "index_bench - Alignment Free Sequence Comparison")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Chenxu Pan <chenxu.pan@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (index_bench "index_bench")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:index_bench")
