import itk
import numpy as np

def GenerateData(data):
    dimension = 6
    mesh1 = itk.meshread("/home/pranjal.sahu/ethicon_aws/deep_learning/deep_learning/train/fixed.vtk", itk.F)
    mesh2 = itk.meshread("/home/pranjal.sahu/ethicon_aws/deep_learning/deep_learning/train/moving.vtk", itk.F)
    data.reserve(mesh1.GetNumberOfPoints())
    for i in range(mesh1.GetNumberOfPoints()):
        point1 = mesh1.GetPoint(i)
        point2 = mesh2.GetPoint(i)
        input_data = [point1[0], point1[1], point1[2], point2[0], point2[1], point2[2]]
        input_data = [float(x) for x in input_data]
        data.push_back(input_data)
    return

data = itk.vector.itkPointD6()
GenerateData(data)


DIMENSION = 3
INLIERS = 90
OUTLIERS = 10
outlierDistance = 20.0

RANSACType = itk.RANSAC[itk.Point[itk.D, 6], itk.D]

transformParameters = itk.vector.D()
maximalDistanceFromPlane = 6
RegistrationEstimatorType = itk.Ransac.LandmarkRegistrationEstimator[6]
registrationEstimator = RegistrationEstimatorType.New()
registrationEstimator.SetMinimalForEstimate(3)
registrationEstimator.SetDelta(maximalDistanceFromPlane)
registrationEstimator.LeastSquaresEstimate(data, transformParameters)

print("Least squares hyper(plane) parameters: [n,a]")
print("[", ", ".join([str(np.round(x, 3)) for x in transformParameters]), "]")

# Performing RANSAC based estimation
desiredProbabilityForNoOutliers = 0.999
ransacEstimator = RANSACType.New()
ransacEstimator.SetData(data)
ransacEstimator.SetNumberOfThreads(8)
ransacEstimator.SetParametersEstimator(registrationEstimator)
percentageOfDataUsed = ransacEstimator.Compute( transformParameters, desiredProbabilityForNoOutliers )

print("Percentage of points used ", percentageOfDataUsed)
print("RANSAC parameters: [n,a]")
print("[", ", ".join([str(np.round(x, 3)) for x in transformParameters]), "]")


