cmake_minimum_required(VERSION 3.20)

project(
  standalone_cpp_sasmol
  VERSION 0.1.0
  DESCRIPTION "Standalone C++20 SasMol parity library"
  LANGUAGES C CXX
)

option(SASMOL_ENABLE_SANITIZERS "Build with AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(SASMOL_ENABLE_VMD_ADAPTER "Build optional VMD IMD socket adapter" OFF)

add_library(sasmol_standalone
  src/calculate.cpp
  src/file_io.cpp
  src/linear_algebra.cpp
  src/molecule.cpp
  src/operate.cpp
  src/overlap.cpp
  src/properties.cpp
  src/selection.cpp
  src/subset.cpp
  src/topology.cpp
  src/view.cpp
)

target_include_directories(sasmol_standalone
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_compile_features(sasmol_standalone PUBLIC cxx_std_20)
set_target_properties(sasmol_standalone PROPERTIES
  C_STANDARD 99
  C_STANDARD_REQUIRED YES
  CXX_STANDARD 20
  CXX_STANDARD_REQUIRED YES
  CXX_EXTENSIONS NO
)

if(SASMOL_ENABLE_VMD_ADAPTER)
  target_sources(sasmol_standalone PRIVATE
    src/vmd_legacy/imd.c
    src/vmd_legacy/view_vmd.c
    src/vmd_legacy/vmdsock.c
  )
  target_include_directories(sasmol_standalone PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/vmd_legacy
  )
  target_compile_definitions(sasmol_standalone PUBLIC SASMOL_ENABLE_VMD_ADAPTER=1)
endif()

if(MSVC)
  target_compile_options(sasmol_standalone PRIVATE /W4)
else()
  target_compile_options(sasmol_standalone PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(SASMOL_ENABLE_SANITIZERS)
  if(MSVC)
    message(FATAL_ERROR "SASMOL_ENABLE_SANITIZERS is currently configured for Clang/GCC-style sanitizers.")
  endif()
  target_compile_options(sasmol_standalone PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(sasmol_standalone PUBLIC -fsanitize=address,undefined)
endif()

add_executable(dcd_roundtrip_writer tools/dcd_roundtrip_writer.cpp)
target_link_libraries(dcd_roundtrip_writer PRIVATE sasmol_standalone)

add_executable(pdb_roundtrip_writer tools/pdb_roundtrip_writer.cpp)
target_link_libraries(pdb_roundtrip_writer PRIVATE sasmol_standalone)

include(CTest)

if(BUILD_TESTING)
  add_executable(test_molecule_core tests/test_molecule_core.cpp)
  target_link_libraries(test_molecule_core PRIVATE sasmol_standalone)
  add_test(NAME molecule_core COMMAND test_molecule_core)

  add_executable(test_file_io_contract tests/test_file_io_contract.cpp)
  target_link_libraries(test_file_io_contract PRIVATE sasmol_standalone)
  target_compile_definitions(test_file_io_contract PRIVATE
    SASMOL_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../src/python/test_sasmol/data"
  )
  add_test(NAME file_io_contract COMMAND test_file_io_contract)

  add_executable(test_calculate tests/test_calculate.cpp)
  target_link_libraries(test_calculate PRIVATE sasmol_standalone)
  target_compile_definitions(test_calculate PRIVATE
    SASMOL_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../src/python/test_sasmol/data"
  )
  add_test(NAME calculate COMMAND test_calculate)

  add_executable(test_operate tests/test_operate.cpp)
  target_link_libraries(test_operate PRIVATE sasmol_standalone)
  target_compile_definitions(test_operate PRIVATE
    SASMOL_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../src/python/test_sasmol/data"
  )
  add_test(NAME operate COMMAND test_operate)

  add_executable(test_linear_algebra tests/test_linear_algebra.cpp)
  target_link_libraries(test_linear_algebra PRIVATE sasmol_standalone)
  add_test(NAME linear_algebra COMMAND test_linear_algebra)

  add_executable(test_overlap tests/test_overlap.cpp)
  target_link_libraries(test_overlap PRIVATE sasmol_standalone)
  add_test(NAME overlap COMMAND test_overlap)

  add_executable(test_properties tests/test_properties.cpp)
  target_link_libraries(test_properties PRIVATE sasmol_standalone)
  target_compile_definitions(test_properties PRIVATE
    SASMOL_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../src/python/test_sasmol/data"
  )
  add_test(NAME properties COMMAND test_properties)

  add_executable(test_view tests/test_view.cpp)
  target_link_libraries(test_view PRIVATE sasmol_standalone)
  add_test(NAME view COMMAND test_view)

  add_executable(test_selection tests/test_selection.cpp)
  target_link_libraries(test_selection PRIVATE sasmol_standalone)
  target_compile_definitions(test_selection PRIVATE
    SASMOL_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../src/python/test_sasmol/data"
  )
  add_test(NAME selection COMMAND test_selection)

  add_executable(test_topology tests/test_topology.cpp)
  target_link_libraries(test_topology PRIVATE sasmol_standalone)
  target_compile_definitions(test_topology PRIVATE
    SASMOL_TOPOLOGY_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data/topology"
    SASMOL_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../src/python/test_sasmol/data"
  )
  add_test(NAME topology COMMAND test_topology)

  add_executable(test_subset tests/test_subset.cpp)
  target_link_libraries(test_subset PRIVATE sasmol_standalone)
  target_compile_definitions(test_subset PRIVATE
    SASMOL_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../src/python/test_sasmol/data"
  )
  add_test(NAME subset COMMAND test_subset)

  add_executable(test_biomt_fixture_semantics tests/test_biomt_fixture_semantics.cpp)
  target_link_libraries(test_biomt_fixture_semantics PRIVATE sasmol_standalone)
  add_test(NAME biomt_fixture_semantics COMMAND test_biomt_fixture_semantics)

  add_executable(test_biomt_metadata_fixture_semantics
    tests/test_biomt_metadata_fixture_semantics.cpp)
  target_link_libraries(test_biomt_metadata_fixture_semantics PRIVATE sasmol_standalone)
  add_test(NAME biomt_metadata_fixture_semantics
    COMMAND test_biomt_metadata_fixture_semantics)

  add_executable(test_dcd_reader tests/test_dcd_reader.cpp)
  target_link_libraries(test_dcd_reader PRIVATE sasmol_standalone)
  target_compile_definitions(test_dcd_reader PRIVATE
    SASMOL_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../src/python/test_sasmol/data"
    SASMOL_TEST_BUILD_DIR="${CMAKE_CURRENT_BINARY_DIR}"
  )
  add_test(NAME dcd_reader COMMAND test_dcd_reader)

  add_executable(test_pdb_helpers tests/test_pdb_helpers.cpp)
  target_link_libraries(test_pdb_helpers PRIVATE sasmol_standalone)
  target_compile_definitions(test_pdb_helpers PRIVATE
    SASMOL_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../src/python/test_sasmol/data"
  )
  add_test(NAME pdb_helpers COMMAND test_pdb_helpers)
endif()
