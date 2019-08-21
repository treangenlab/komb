# Generated by Boost 1.71.0

if(TARGET Boost::graph_parallel)
  return()
endif()

if(Boost_VERBOSE OR Boost_DEBUG)
  message(STATUS "Found boost_graph_parallel ${boost_graph_parallel_VERSION} at ${boost_graph_parallel_DIR}")
endif()

# Compute the include and library directories relative to this file.
get_filename_component(_BOOST_CMAKEDIR "${CMAKE_CURRENT_LIST_DIR}/../" ABSOLUTE)
get_filename_component(_BOOST_INCLUDEDIR "${_BOOST_CMAKEDIR}/../../include/" ABSOLUTE)
get_filename_component(_BOOST_LIBDIR "${_BOOST_CMAKEDIR}/../" ABSOLUTE)

# Create imported target Boost::graph_parallel
add_library(Boost::graph_parallel UNKNOWN IMPORTED)

set_target_properties(Boost::graph_parallel PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_BOOST_INCLUDEDIR}"
  INTERFACE_COMPILE_DEFINITIONS "BOOST_ALL_NO_LIB"
)

include(${CMAKE_CURRENT_LIST_DIR}/../BoostDetectToolset-1.71.0.cmake)

if(Boost_DEBUG)
  message(STATUS "Scanning ${CMAKE_CURRENT_LIST_DIR}/libboost_graph_parallel-variant*.cmake")
endif()

file(GLOB __boost_variants "${CMAKE_CURRENT_LIST_DIR}/libboost_graph_parallel-variant*.cmake")

macro(_BOOST_SKIPPED fname reason)
  if(Boost_VERBOSE OR Boost_DEBUG)
    message(STATUS "  [ ] ${fname}")
  endif()
  list(APPEND __boost_skipped "${fname} (${reason})")
endmacro()

foreach(f IN LISTS __boost_variants)
  if(Boost_DEBUG)
    message(STATUS "  Including ${f}")
  endif()
  include(${f})
endforeach()

unset(_BOOST_LIBDIR)
unset(_BOOST_INCLUDEDIR)
unset(_BOOST_CMAKEDIR)

get_target_property(__boost_configs Boost::graph_parallel IMPORTED_CONFIGURATIONS)

if(__boost_variants AND NOT __boost_configs)
  set(__boost_message "No suitable build variant has been found.")
  if(__boost_skipped)
    set(__boost_message "${__boost_message}\nThe following variants have been tried and rejected:")
    foreach(s IN LISTS __boost_skipped)
      set(__boost_message "${__boost_message}\n* ${s}")
    endforeach()
  endif()
  set(boost_graph_parallel_FOUND 0)
  set(boost_graph_parallel_NOT_FOUND_MESSAGE ${__boost_message})
  unset(__boost_message)
  unset(__boost_skipped)
  unset(__boost_configs)
  unset(__boost_variants)
  unset(_BOOST_GRAPH_PARALLEL_DEPS)
  return()
endif()

unset(__boost_skipped)
unset(__boost_configs)
unset(__boost_variants)

if(_BOOST_GRAPH_PARALLEL_DEPS)
  list(REMOVE_DUPLICATES _BOOST_GRAPH_PARALLEL_DEPS)
  if(Boost_VERBOSE OR Boost_DEBUG)
    message(STATUS "Adding boost_graph_parallel dependencies: ${_BOOST_GRAPH_PARALLEL_DEPS}")
  endif()
endif()

foreach(dep_boost_graph_parallel IN LISTS _BOOST_GRAPH_PARALLEL_DEPS)
  set(_BOOST_QUIET)
  if(boost_graph_parallel_FIND_QUIETLY)
    set(_BOOST_QUIET QUIET)
  endif()
  set(_BOOST_REQUIRED)
  if(boost_graph_parallel_FIND_REQUIRED)
    set(_BOOST_REQUIRED REQUIRED)
  endif()
  get_filename_component(_BOOST_CMAKEDIR "${CMAKE_CURRENT_LIST_DIR}/../" ABSOLUTE)
  find_package(boost_${dep_boost_graph_parallel} 1.71.0 EXACT CONFIG ${_BOOST_REQUIRED} ${_BOOST_QUIET} HINTS ${_BOOST_CMAKEDIR})
  set_property(TARGET Boost::graph_parallel APPEND PROPERTY INTERFACE_LINK_LIBRARIES Boost::${dep_boost_graph_parallel})
  unset(_BOOST_QUIET)
  unset(_BOOST_REQUIRED)
  unset(_BOOST_CMAKEDIR)
  if(NOT boost_${dep_boost_graph_parallel}_FOUND)
    set(boost_graph_parallel_FOUND 0)
    set(boost_graph_parallel_NOT_FOUND_MESSAGE "A required dependency, boost_${dep_boost_graph_parallel}, has not been found.")
    unset(_BOOST_GRAPH_PARALLEL_DEPS)
    return()
  endif()
endforeach()

unset(_BOOST_GRAPH_PARALLEL_DEPS)

mark_as_advanced(boost_graph_parallel_DIR)
