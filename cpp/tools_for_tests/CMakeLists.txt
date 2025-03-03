add_library(tools_for_tests STATIC tiffwrappers.cc tiffwrappers.h)
target_link_libraries(tools_for_tests ${TIFF_LIBRARY})
target_include_directories(tools_for_tests PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/.."
  "${PROJECT_BINARY_DIR}/include/"
)
target_include_directories(tools_for_tests SYSTEM PUBLIC ${TIFF_INCLUDE_DIR})
if(spdlog_INCLUDE_DIR)
  target_include_directories(tools_for_tests SYSTEM PUBLIC ${spdlog_INCLUDE_DIR})
endif()
if(eigen_INCLUDE_DIR)
  target_include_directories(tools_for_tests SYSTEM PUBLIC ${eigen_INCLUDE_DIR})
endif()
# On different platforms the CMakeDeps generator in conan seems to install eigen
# as either "eigen" or "Eigen3" because the recipe does not explicitly define the
# name (yet). To work around this we have to check for both.
if(Eigen3_INCLUDE_DIR)
  target_include_directories(tools_for_tests SYSTEM PUBLIC ${Eigen3_INCLUDE_DIR})
endif()


# Simple manual tester of read/write tiff
add_executable(copy_tiff copy_tiff.cc)
target_link_libraries(copy_tiff tools_for_tests sopt)

configure_file("directories.in.h" "${PROJECT_BINARY_DIR}/include/tools_for_tests/directories.h")
#file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/outputs")
#file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/outputs/sdmm")

