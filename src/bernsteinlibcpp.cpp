#include <iostream>
#include <Eigen/Dense>
#include <bernsteinlibcpp.h>
#include <algorithm>
#include <string>
#include <boost/math/special_functions/binomial.hpp>

// TODO: try putting everything on the .h maybe

using namespace Eigen;

// namespace BernsteinLib;
// using binom = typename boost::math::binomial_coefficient<double>;

int sum() {
Matrix2d a;
  a << 1, 2,
	   3, 4;
	Vector3d v(1,2,3);
	std::cout << "ola" << std::endl;
	std::cout << "a * 2.5 =\n" << a * 2.5 << std::endl;
	std::cout << "0.1 * v =\n" << 0.1 * v << std::endl;
	std::cout << "Doing v *= 2;" << std::endl;
	v *= 2;
	std::cout << "Now v =\n" << v << std::endl;
	std::cout <<   boost::math::binomial_coefficient<double>(4, 2)  << std::endl; 
	std::cout << "test\n";
	return 0;
}

float basis(int k, int x, float tau) {
	return 0;
}

// vim:set tabstop=4

