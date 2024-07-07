# Install script for directory: D:/Programming/ClionProjects/TDD_C/code/deps/xtensor

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/xtensor" TYPE FILE FILES
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xaccessible.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xaccumulator.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xadapt.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xarray.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xassign.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xaxis_iterator.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xaxis_slice_iterator.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xblockwise_reducer.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xblockwise_reducer_functors.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xbroadcast.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xbuffer_adaptor.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xbuilder.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xchunked_array.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xchunked_assign.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xchunked_view.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xcomplex.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xcontainer.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xcsv.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xdynamic_view.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xeval.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xexception.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xexpression.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xexpression_holder.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xexpression_traits.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xfixed.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xfunction.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xfunctor_view.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xgenerator.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xhistogram.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xindex_view.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xinfo.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xio.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xiterable.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xiterator.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xjson.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xlayout.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xmanipulation.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xmasked_view.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xmath.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xmime.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xmultiindex_iterator.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xnoalias.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xnorm.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xnpy.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xoffset_view.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xoperation.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xoptional.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xoptional_assembly.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xoptional_assembly_base.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xoptional_assembly_storage.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xpad.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xrandom.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xreducer.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xrepeat.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xscalar.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xsemantic.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xset_operation.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xshape.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xslice.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xsort.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xstorage.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xstrided_view.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xstrided_view_base.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xstrides.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xtensor.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xtensor_config.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xtensor_forward.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xtensor_simd.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xutils.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xvectorize.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xview.hpp"
    "D:/Programming/ClionProjects/TDD_C/code/deps/xtensor/include/xtensor/xview_utils.hpp"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/cmake/xtensor" TYPE FILE FILES
    "D:/Programming/ClionProjects/TDD_C/cmake-build-debug-mingw/code/deps/xtensor/xtensorConfig.cmake"
    "D:/Programming/ClionProjects/TDD_C/cmake-build-debug-mingw/code/deps/xtensor/xtensorConfigVersion.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/xtensor/xtensorTargets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/xtensor/xtensorTargets.cmake"
         "D:/Programming/ClionProjects/TDD_C/cmake-build-debug-mingw/code/deps/xtensor/CMakeFiles/Export/0f2a327e949144b0b747c96acb1eb12f/xtensorTargets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/xtensor/xtensorTargets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/xtensor/xtensorTargets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/cmake/xtensor" TYPE FILE FILES "D:/Programming/ClionProjects/TDD_C/cmake-build-debug-mingw/code/deps/xtensor/CMakeFiles/Export/0f2a327e949144b0b747c96acb1eb12f/xtensorTargets.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/pkgconfig" TYPE FILE FILES "D:/Programming/ClionProjects/TDD_C/cmake-build-debug-mingw/code/deps/xtensor/xtensor.pc")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES "D:/Programming/ClionProjects/TDD_C/cmake-build-debug-mingw/xtensor.hpp")
endif()

