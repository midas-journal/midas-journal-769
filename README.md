ITKRANSAC
=================================

[<img src="https://github.com/InsightSoftwareConsortium/ITKRANSAC/actions/workflows/build-test-package.yml/badge.svg">](https://github.com/InsightSoftwareConsortium/ITKRANSAC/actions)
[<img src="https://img.shields.io/pypi/v/itk-ransac.svg">](https://pypi.python.org/pypi/itk-ransac)

Overview
--------
This is the source code for a C++ templated implementation of the RANSAC
algorithm. The implementation is multi-threaded. This repository is only for 
pointset registratation and differs slightly from the original generic implementation due to optimization. 

For implementation related to plane and sphere estimation
please refer https://github.com/midas-journal/midas-journal-769.

The code is "in the style of ITK". That is, it is very similar to the official
ITK style but does not follow all of the required conventions.

**Manifest:**

1. [RANSAC {h,txx}](./include/itkRANSAC.hxx) - Multi-threaded implementation of the generic RANSAC algorithm.
2. [ParametersEstimator.{h,hxx}](./include/itkParametersEstimator.hxx) - Super class of all parameter estimation objects
that can be used with the RANSAC algorithm. This is an abstract class that
defines an interface.
3. [itkLandmarkRegistrationEstimator.{h,hxx}](./include/itkLandmarkRegistrationEstimator.hxx) - Estimation code for landmark based pointset registration.
4. [Testing/*.cxx](./test/itkRansacTest_LandmarkRegistration) - Test for the PointSet registration using landmark points.

<br/><br/>
**Sample Usage in Python for 3D PointSet is shown here:**
```python
data = itk.vector[itk.Point[itk.D, 6]]()
agreeData = itk.vector[itk.Point[itk.D, 6]]()
GenerateData(data, agreeData)

transformParameters = itk.vector.D()

maximumDistance = inlier_value
RegistrationEstimatorType = itk.Ransac.LandmarkRegistrationEstimator[6]
registrationEstimator = RegistrationEstimatorType.New()
registrationEstimator.SetMinimalForEstimate(number_of_ransac_points)
registrationEstimator.SetAgreeData(agreeData)
registrationEstimator.SetDelta(maximumDistance)
registrationEstimator.LeastSquaresEstimate(data, transformParameters)

desiredProbabilityForNoOutliers = 0.99
RANSACType = itk.RANSAC[itk.Point[itk.D, 6], itk.D]
ransacEstimator = RANSACType.New()
ransacEstimator.SetData(data)
ransacEstimator.SetAgreeData(agreeData)
ransacEstimator.SetMaxIteration(number_of_iterations)
ransacEstimator.SetNumberOfThreads(8)
ransacEstimator.SetParametersEstimator(registrationEstimator)
    
percentageOfDataUsed = ransacEstimator.Compute( transformParameters, desiredProbabilityForNoOutliers )
for i in transformParameters:
  print(i)
```

<br/><br/>
**Landmarks can be obtained by performing feature matching.**

For this one can use the [ITKFPFH](https://github.com/InsightSoftwareConsortium/ITKFPFH) library.
