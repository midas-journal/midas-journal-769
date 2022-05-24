itk_module_test()

set(RansacTests
	itkRansacTest_PlaneEstimation.cxx
  )

CreateTestDriver(Ransac "${Ransac-Test_LIBRARIES}" "${RansacTests}")

itk_add_test(NAME itkRansacTest_PlaneEstimation
	COMMAND RansacTestDriver
  itkRansacTest_PlaneEstimation
  )

