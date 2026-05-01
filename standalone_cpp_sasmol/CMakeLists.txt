cmake_minimum_required(VERSION 3.20)

project(
  standalone_cpp_sasmol
  VERSION 0.1.0
  DESCRIPTION "Standalone C++20 SasMol parity library"
  LANGUAGES CXX
)

add_library(sasmol_standalone
  src/molecule.cpp
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

include(CTest)

if(BUILD_TESTING)
  add_executable(test_molecule_core tests/test_molecule_core.cpp)
  target_link_libraries(test_molecule_core PRIVATE sasmol_standalone)
  add_test(NAME molecule_core COMMAND test_molecule_core)
endif()
