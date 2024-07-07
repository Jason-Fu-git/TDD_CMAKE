# Install script for directory: D:/Programming/ClionProjects/TDD_C/code/deps/xtl

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files (x86)/TDD_C")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "C:/Users/Jason Fu/AppData/Local/Programs/CLion/bin/mingw/bin/objdump.exe")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/xtl" TYPE FILE FILES
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xany.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xbasic_fixed_string.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xbase64.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xclosure.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xcompare.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xcomplex.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xcomplex_sequence.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xspan.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xspan_impl.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xdynamic_bitset.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xfunctional.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xhalf_float.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xhalf_float_impl.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xhash.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xhierarchy_generator.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xiterator_base.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xjson.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xmasked_value_meta.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xmasked_value.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xmeta_utils.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xmultimethods.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xoptional_meta.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xoptional.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xoptional_sequence.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xplatform.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xproxy_wrapper.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xsequence.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xsystem.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xtl_config.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xtype_traits.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xvariant.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xvariant_impl.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtl/include/xtl/xvisitor.hpp"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/cmake/xtl" TYPE FILE FILES
    "D:/Programming/ClionProjects/TDD_C/cmake-build-debug-mingw/code/deps/xtl/xtlConfig.cmake"
    "D:/Programming/ClionProjects/TDD_C/cmake-build-debug-mingw/code/deps/xtl/xtlConfigVersion.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/xtl/xtlTargets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/xtl/xtlTargets.cmake"
         "D:/Programming/ClionProjects/TDD_C/cmake-build-debug-mingw/code/deps/xtl/CMakeFiles/Export/2fc63ec57839ed115fc15a5438bb5aec/xtlTargets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/xtl/xtlTargets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/xtl/xtlTargets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/cmake/xtl" TYPE FILE FILES "D:/Programming/ClionProjects/TDD_C/cmake-build-debug-mingw/code/deps/xtl/CMakeFiles/Export/2fc63ec57839ed115fc15a5438bb5aec/xtlTargets.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/pkgconfig" TYPE FILE FILES "D:/Programming/ClionProjects/TDD_C/cmake-build-debug-mingw/code/deps/xtl/xtl.pc")
endif()

