cmake_minimum_required(VERSION 3.20)

project(
  standalone_cpp_sasmol
  VERSION 0.1.0
  DESCRIPTION "Standalone C++20 SasMol parity library"
  LANGUAGES CXX
)

option(SASMOL_ENABLE_SANITIZERS "Build with AddressSanitizer and UndefinedBehaviorSanitizer" OFF)

add_library(sasmol_standalone
  src/calculate.cpp
  src/file_io.cpp
  src/molecule.cpp
  src/operate.cpp
  src/properties.cpp
)

target_include_directories(sasmol_standalone
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_compile_features(sasmol_standalone PUBLIC cxx_std_20)
set_target_properties(sasmol_standalone PROPERTIES
  CXX_STANDARD 20
  CXX_STANDARD_REQUIRED YES
  CXX_EXTENSIONS NO
)

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
