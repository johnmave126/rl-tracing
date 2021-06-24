cmake_minimum_required (VERSION 2.8.3)
project(rl-tracing)

add_subdirectory(ext ext_build)

include_directories(
  # include files
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  # tinyformat string formatting library
  ${TFM_INCLUDE_DIR}
  # Eigen linear algebra library
  ${EIGEN_INCLUDE_DIR}
  # OpenEXR high dynamic range bitmap library
  ${OPENEXR_INCLUDE_DIRS}
  # Intel Thread Building Blocks
  ${TBB_INCLUDE_DIR}
  # Pseudorandom number generator
  ${PCG32_INCLUDE_DIR}
  # PugiXML parser
  ${PUGIXML_INCLUDE_DIR}
  # Helper functions for statistical hypothesis tests
  ${HYPOTHESIS_INCLUDE_DIR}
  # GLFW library for OpenGL context creation
  ${GLFW_INCLUDE_DIR}
  # GLEW library for accessing OpenGL functions
  ${GLEW_INCLUDE_DIR}
  # NanoVG drawing library
  ${NANOVG_INCLUDE_DIR}
  # NanoGUI user interface library
  ${NANOGUI_INCLUDE_DIR}
  ${NANOGUI_EXTRA_INCS}
  # Portable filesystem API
  ${FILESYSTEM_INCLUDE_DIR}
)

# The following lines build the main executable. If you add a source
# code file to tracer, be sure to include it in this list.
add_executable(rl-tracer

  # Header files
  include/rl-tracer/bbox.h
  include/rl-tracer/bitmap.h
  include/rl-tracer/block.h
  include/rl-tracer/bsdf.h
  include/rl-tracer/accel.h
  include/rl-tracer/camera.h
  include/rl-tracer/color.h
  include/rl-tracer/common.h
  include/rl-tracer/dpdf.h
  include/rl-tracer/frame.h
  include/rl-tracer/integrator.h
  include/rl-tracer/emitter.h
  include/rl-tracer/mesh.h
  include/rl-tracer/object.h
  include/rl-tracer/parser.h
  include/rl-tracer/proplist.h
  include/rl-tracer/ray.h
  include/rl-tracer/rfilter.h
  include/rl-tracer/sampler.h
  include/rl-tracer/scene.h
  include/rl-tracer/timer.h
  include/rl-tracer/transform.h
  include/rl-tracer/vector.h
  include/rl-tracer/warp.h
  include/rl-tracer/lightprobe.h
  include/rl-tracer/guider.h

  # Source code files
  src/bitmap.cpp
  src/block.cpp
  src/accel.cpp
  src/chi2test.cpp
  src/common.cpp
  src/diffuse.cpp
  src/gui.cpp
  src/independent.cpp
  src/main.cpp
  src/mesh.cpp
  src/obj.cpp
  src/object.cpp
  src/parser.cpp
  src/perspective.cpp
  src/proplist.cpp
  src/rfilter.cpp
  src/scene.cpp
  src/ttest.cpp
  src/warp.cpp
  src/microfacet.cpp
  src/mirror.cpp
  src/dielectric.cpp
  src/normals.cpp
  src/simple.cpp
  src/ao.cpp
  src/lightprobe.cpp
  src/area.cpp
  src/whitted.cpp
  src/path_simple.cpp
  src/path.cpp
  src/roughdielectric.cpp
  src/qtable.cpp
  src/path_guided.cpp
  src/path_naive.cpp
  src/path_guided_simple.cpp
  src/path_guided_mis.cpp
  src/qtable_sphere.cpp
  src/probe.cpp
)

add_definitions(${NANOGUI_EXTRA_DEFS})

# The following lines build the warping test application
add_executable(warptest
  include/rl-tracer/color.h
  include/rl-tracer/warp.h
  include/rl-tracer/vector.h
  include/rl-tracer/bitmap.h
  include/rl-tracer/lightprobe.h
  src/bitmap.cpp
  src/warp.cpp
  src/warptest.cpp
  src/microfacet.cpp
  src/object.cpp
  src/proplist.cpp
  src/common.cpp
  src/lightprobe.cpp
)

target_link_libraries(rl-tracer tbb_static pugixml IlmImf nanogui ${NANOGUI_EXTRA_LIBS})
target_link_libraries(warptest tbb_static IlmImf nanogui ${NANOGUI_EXTRA_LIBS})

# vim: set et ts=2 sw=2 ft=cmake nospell:
