//// //////////////////////////////////////////////////////////////////////////////////////////////////////
//// Cpp functions 2015-3-11 Thong Pham
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <fstream>
#include <cmath>
#include <random>
#include <cstdlib>
#include <stdio.h>
#include <float.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]


double my_zeroin(double, double, std::function <double (double)>, double, long);
