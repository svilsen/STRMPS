# Install script for directory: /home/svilsen/seqan/seqan-src/apps

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Debug")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/fiona/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/splazers/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/razers2/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/seqan_tcoffee/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/tree_recon/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/seqcons/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/dfi/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/insegt/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/rabema/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/seqcons2/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/seqan_flexbar/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/searchjoin/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/pair_align/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/mason2/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/gustaf/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/razers3/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/razers/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/snp_store/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/stellar/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/sam2matrix/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/alf/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/bs_tools/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/samcat/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/ngs_roi/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/yara/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/param_chooser/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/fx_tools/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/sak/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/micro_razers/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/rep_sep/cmake_install.cmake")
  INCLUDE("/home/svilsen/seqan/seqan-build/debug/apps/sgip/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

