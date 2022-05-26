import itk
import numpy as np

def GenerateData(numInliers, numOutliers, outlierDistance, data, planeParameters ):
    dimension = 3
    normal = np.zeros(3)
    noise = np.zeros(3)
    tmp = np.zeros(3)
    
    pointOnPlane = np.zeros(3)
    randomPoint = np.zeros(3)
    
    noiseStandardDeviation = 0.4 # noise standard deviation
    coordinateMax = 1000.0
    
    planeParameters.clear()

    # generate points on random (hyper) plane
    for i in range(dimension):
        normal[i] = np.random.uniform()
        pointOnPlane[i] = np.random.uniform(-coordinateMax, coordinateMax)
  
    normal = normal / np.sqrt(np.sum(normal**2))

    #normal.Normalize()
    print(normal)
    print(pointOnPlane)

    for i in range(dimension):
        planeParameters.push_back( normal[i])
    
    for i in range(dimension):
        planeParameters.push_back( pointOnPlane[i])

    # generate inliers
    for i in range(numInliers):
        for j in range(dimension):
            randomPoint[j] = np.random.uniform( -coordinateMax, coordinateMax )
            noise[j] = np.random.normal( noiseStandardDeviation )
    
        # project random point onto the plane and add noise
        tmp = randomPoint - pointOnPlane
        randomPoint = pointOnPlane + noise + (tmp - (tmp*normal)*normal)
        data.push_back( randomPoint)
    
    # generate outliers (via rejection)
    count_Outliers = 0
    while(count_Outliers < numOutliers):
        for j in range(dimension):
            randomPoint[j] = np.random.uniform( -coordinateMax, coordinateMax )
        tmp = randomPoint - pointOnPlane
        if (np.abs(np.sum(np.multiply(tmp, normal))) >= outlierDistance ):
            data.push_back(randomPoint)
            count_Outliers = numOutliers+1
    
DIMENSION = 3
INLIERS = 90
OUTLIERS = 10
outlierDistance = 20.0

maximalDistanceFromPlane = 0.5

PlaneEstimatorType = itk.PlaneParametersEstimator[3]
#RANSACType = itk.RANSAC[itk.Point[itk.D, 3], itk.D]

true_array = [0.460999, 0.721842, 0.516164, -874.896, 897.281, 101.569]

data = itk.vector.itkPointD3()
planeParameters = itk.vector.D()
truePlaneParameters = itk.vector.D()

for i, k in enumerate(true_array):
    truePlaneParameters.push_back(true_array[i])

print(truePlaneParameters)


GenerateData(INLIERS, OUTLIERS, outlierDistance, data, truePlaneParameters)
print('Data Generation Done')
#data_np = itk.array_from_vector_container(data)
#print(data_np.shape)

planeEstimator = PlaneEstimatorType.New()
planeEstimator.SetDelta(maximalDistanceFromPlane)
planeEstimator.LeastSquaresEstimate(data, planeParameters)


print("Least squares hyper(plane) parameters: [n,a]")
for i in range(2*DIMENSION):
    print(planeParameters[i])

dotProduct = 0
for i in range(DIMENSION):
    dotProduct += planeParameters[i]*truePlaneParameters[i]
	  
print("dotProduct  ", dotProduct)