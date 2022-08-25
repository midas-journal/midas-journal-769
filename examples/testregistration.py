import itk
import numpy as np
import random
import copy

def GenerateData(data, agreeData):
    dimension = 6
    mesh1 = itk.meshread("/data/Apedata/Slicer-cli-outputs/EJEJHH_moving_corr.vtk", itk.F)
    mesh2 = itk.meshread("/data/Apedata/Slicer-cli-outputs/EJEJHH_fixed_corr.vtk", itk.F)
    data.reserve(mesh1.GetNumberOfPoints())
    for i in range(mesh1.GetNumberOfPoints()):
        point1 = mesh1.GetPoint(i)
        point2 = mesh2.GetPoint(i)
        input_data = [point1[0], point1[1], point1[2], point2[0], point2[1], point2[2]]
        input_data = [float(x) for x in input_data]
        data.push_back(input_data)
    
    mesh1 = itk.meshread("/data/Apedata/Slicer-cli-outputs/EJEJHH_movingMeshPointsBefore.vtk", itk.F)
    mesh2 = itk.meshread("/data/Apedata/Slicer-cli-outputs/EJEJHH_fixedMeshPoints.vtk", itk.F)
    count_min = int(np.min([mesh1.GetNumberOfPoints(), mesh2.GetNumberOfPoints()]))

    mesh1_points  = copy.deepcopy(itk.array_from_vector_container(mesh1.GetPoints()))
    mesh2_points  = copy.deepcopy(itk.array_from_vector_container(mesh2.GetPoints()))

    agreeData.reserve(count_min)
    for i in range(count_min):
        point1 = mesh1_points[i]
        point2 = mesh2_points[i]
        input_data = [point1[0], point1[1], point1[2], point2[0], point2[1], point2[2]]
        input_data = [float(x) for x in input_data]
        agreeData.push_back(input_data)
    return

data = itk.vector.itkPointD6()
agreeData = itk.vector.itkPointD6()
GenerateData(data, agreeData)


DIMENSION = 3
INLIERS = 90
OUTLIERS = 10
outlierDistance = 20.0

RANSACType = itk.RANSAC[itk.Point[itk.D, 6], itk.D]

transformParameters = itk.vector.D()
maximalDistanceFromPlane = 1.7
RegistrationEstimatorType = itk.Ransac.LandmarkRegistrationEstimator[6]
registrationEstimator = RegistrationEstimatorType.New()
registrationEstimator.SetMinimalForEstimate(100)
registrationEstimator.SetDelta(maximalDistanceFromPlane)
registrationEstimator.LeastSquaresEstimate(data, transformParameters)

print("Least squares hyper(plane) parameters: [n,a]")
print("[", ", ".join([str(np.round(x, 3)) for x in transformParameters]), "]")

# Performing RANSAC based estimation
desiredProbabilityForNoOutliers = 0.999
ransacEstimator = RANSACType.New()
ransacEstimator.SetData(data)
ransacEstimator.SetAgreeData(agreeData)
ransacEstimator.SetMaxIteration(10000)
ransacEstimator.SetNumberOfThreads(8)
ransacEstimator.SetParametersEstimator(registrationEstimator)
percentageOfDataUsed = ransacEstimator.Compute( transformParameters, desiredProbabilityForNoOutliers )

print("Percentage of points used ", percentageOfDataUsed)
print("RANSAC parameters: [n,a]")
print("[", ", ".join([str(np.round(x, 3)) for x in transformParameters]), "]")

transform = itk.Similarity3DTransform.D.New()
p = transform.GetParameters()
f = transform.GetFixedParameters()
for i in range(7):
    p.SetElement(i, transformParameters[i])
counter = 0
for i in range(7, 10):
    f.SetElement(counter, transformParameters[i])
    counter = counter + 1
transform.SetParameters(p)
transform.SetFixedParameters(f)

mesh1 = itk.meshread("/data/Apedata/Slicer-cli-outputs/EJEJHH_movingMesh.vtk", itk.F)
output = itk.transform_mesh_filter(mesh1, transform=transform)
itk.meshwrite(output, '/data/Apedata/Slicer-cli-outputs/EJEJHH_movingMesh_transformed.vtk')


