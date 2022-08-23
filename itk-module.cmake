# the top-level README is used for describing this module, just
# re-used it for documentation here
get_filename_component(MY_CURRENT_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
file(READ "${MY_CURRENT_DIR}/README.md" DOCUMENTATION)

itk_module(Ransac
  ENABLE_SHARED
  COMPILE_DEPENDS
    ITKCommon
    ITKMesh
    ITKIOMeshBase
    ITKRegistrationCommon
    ITKMetricsv4
    ITKRegistrationMethodsv4
    ITKOptimizersv4
    ITKOptimizers
  TEST_DEPENDS
    ITKTestKernel
    ITKMetaIO
  DESCRIPTION
    "${DOCUMENTATION}"
  EXCLUDE_FROM_DEFAULT
  ENABLE_SHARED
)
