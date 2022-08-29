ITKRANSAC
=================================

Overview
--------
This is the source code for a C++ templated implementation of the RANSAC
algorithm. The implementation is multi-threaded. This repository is only for 
pointset registratation and differs slightly from the original generic implementation due to optimization. 

For implementation related to plane and sphere estimation
please refer https://github.com/midas-journal/midas-journal-769.

The code is "in the style of ITK". That is, it is very similar to the official
ITK style but does not follow all of the required conventions.

Manifest:

RANSAC.{h,txx} - Multi-threaded implementation of the generic RANSAC algorithm.

ParametersEstimator.{h,txx} - Super class of all parameter estimation objects
that can be used with the RANSAC algorithm. This is an abstract class that
defines an interface.

itkLandmarkRegistrationEstimator.{h,hxx} - Estimation code for landmark based pointset registration.

Testing/*.cxx - Test for the PointSet registration using landmark points.
