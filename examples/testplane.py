import itk
import numpy as np

def GenerateData(numInliers, numOutliers, outlierDistance, data, planeParameters ):
    dimension = 3
    normal = np.zeros(dimension)
    noise = np.zeros(dimension)
    tmp = np.zeros(dimension)
    
    pointOnPlane = np.zeros(dimension)
    randomPoint = np.zeros(dimension)
    
    # noise standard deviation
    noiseStandardDeviation = 0.4
    coordinateMax = 1000.0
    
    planeParameters.clear()

    # generate points on random (hyper) plane
    for i in range(dimension):
        normal[i] = np.random.uniform()
        pointOnPlane[i] = np.random.uniform(-coordinateMax, coordinateMax)
  
    normal = normal / np.sqrt(np.sum(normal**2))

    for i in range(dimension):
        planeParameters.push_back( normal[i])
    
    for i in range(dimension):
        planeParameters.push_back( pointOnPlane[i])

    # generate inliers
    for i in range(numInliers):
        for j in range(dimension):
            randomPoint[j] = np.random.uniform( -coordinateMax, coordinateMax )
            noise[j] = np.random.normal(0, noiseStandardDeviation)
    
        # project random point onto the plane and add noise
        tmp = randomPoint - pointOnPlane
        
        randomPoint = pointOnPlane + noise + (tmp - (np.sum(tmp*normal))*normal)
        data.push_back([randomPoint[0], randomPoint[1], randomPoint[2]])
    
    # generate outliers (via rejection)
    count_Outliers = 0
    while(count_Outliers < numOutliers):
        for j in range(dimension):
            randomPoint[j] = np.random.uniform( -coordinateMax, coordinateMax )
        
        tmp = randomPoint - pointOnPlane
        if (np.abs(np.sum(tmp * normal)) >= outlierDistance ):
            data.push_back([randomPoint[0], randomPoint[1], randomPoint[2]])
            count_Outliers = count_Outliers+1
    
DIMENSION = 3
INLIERS = 90
OUTLIERS = 10
outlierDistance = 20.0

maximalDistanceFromPlane = 0.5

PlaneEstimatorType = itk.PlaneParametersEstimator[3]
RANSACType = itk.RANSAC[itk.Point[itk.D, 3], itk.D]

data = itk.vector.itkPointD3()
planeParameters = itk.vector.D()
truePlaneParameters = itk.vector.D()

GenerateData(INLIERS, OUTLIERS, outlierDistance, data, truePlaneParameters)

print("True plane parameters")
print("[", ", ".join([str(np.round(x, 3)) for x in truePlaneParameters]), "]")


print('Data Generation Done')
print('Number of samples generated ', data.size())


planeEstimator = PlaneEstimatorType.New()
planeEstimator.SetDelta(maximalDistanceFromPlane)
planeEstimator.LeastSquaresEstimate(data, planeParameters)

print("Least squares hyper(plane) parameters: [n,a]")
print("[", ", ".join([str(np.round(x, 3)) for x in planeParameters]), "]")


dotProduct = 0
for i in range(DIMENSION):
    dotProduct += truePlaneParameters[i]*planeParameters[i]

print("dotProduct  ", dotProduct)

dotProduct = 0
print("Check if computed point is on known plane [0=correct]")
for i in range(DIMENSION):
    dotProduct += (planeParameters[DIMENSION+i] - truePlaneParameters[DIMENSION+i])*truePlaneParameters[i]
print(dotProduct)



print('Size of data after LEast Squares ', data.size())
# Performing RANSAC based estimation
desiredProbabilityForNoOutliers = 0.999
ransacEstimator = RANSACType.New()
ransacEstimator.SetData( data )
ransacEstimator.SetNumberOfThreads(8)
ransacEstimator.SetParametersEstimator(planeEstimator)
percentageOfDataUsed = ransacEstimator.Compute( planeParameters, desiredProbabilityForNoOutliers )



print('percentageOfDataUsed  is ', percentageOfDataUsed)

print("RANSAC hyper(plane) parameters: [n,a]")
print("[", ", ".join([str(np.round(x, 3)) for x in planeParameters]), "]")


dotProduct = 0
for i in range(DIMENSION):
    dotProduct += truePlaneParameters[i]*planeParameters[i]
print("dotProduct  ", dotProduct)


dotProduct = 0
print("Check if computed point is on known plane [0=correct]")
for i in range(DIMENSION):
    dotProduct += (planeParameters[DIMENSION+i] - truePlaneParameters[DIMENSION+i])*truePlaneParameters[i]
print(dotProduct)
