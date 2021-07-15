#include <iostream>
#include <Eigen/Dense>
#include <bernsteinlibcpp.h>
#include <algorithm>
#include <string>
#include <boost/math/special_functions/binomial.hpp>
#include <iostream>

// TODO: try putting everything on the .h maybe

namespace BernsteinLib {

// https://stackoverflow.com/a/33786772


template <typename T>
const auto comb = boost::math::binomial_coefficient<T>;

double basis(int k, int n, double tau) {
	// return boost::math::binomial_coefficient<double>(n, k)*std::pow(1-tau, n-k)*std::pow(tau,k);
	return comb<double>(n, k)*std::pow(1-tau, n-k)*std::pow(tau,k);
}


/* Returns a matrix for degree elevation */
Eigen::MatrixXd degrelevmat(int n, int m) {
	Eigen::MatrixXd mat(m+1, n+1);
	for (int i = 0; i < m+1; i++) {
		for (int j = 0; j < n+1; j++) {
			mat(i, j) = (std::max(0, i - m + n) <= j and j <= std::min(n, i))? comb<double>(i, j) * comb<double>(m-i, n-j)/comb<double>(m, n): 0;
		}
	}
	return mat;
}


/*Performs degree elevation*/
Eigen::MatrixXd degrelev(Eigen::MatrixXd p, int m) {
	return (p.rows() -1 > m)? p : degrelevmat(p.rows() -1, m)*p;
}


/*Returns an integral/antiderivative matrix
 *Returns an integral/antiderivative matrix that must be multiplied by a column vector
 *to obtain the control points of the integral function"""
 */
Eigen::MatrixXd antiderivmat(int n, double t) {
	Eigen::MatrixXd mat(n+2, n+1);
	double val = t/(n+1);
	for(int i=0; i<n+2; i++) {
		for(int j=0; j<n+1; j++) {
			mat(i, j) = (i > j)? val: 0;
		}
	}
	return mat;
}


/*Returns the control points of the integral function of for the control points of p*/
Eigen::MatrixXd antideriv(Eigen::MatrixXd p, double t, Eigen::MatrixXd p0) {
	std::cout <<"p0.replicate(p.rows(),1): " << p0.replicate(p.rows()+1,1);
	std::cout <<"\nantiderivmat(p.rows()-1, t): " << antiderivmat(p.rows()-1, t);
	std::cout <<"\np: " << p << std::endl;
	return p0.replicate(p.rows()+1,1) + antiderivmat(p.rows()-1, t) * p;
}


/*Returns a matrix to perform derivation*/
Eigen::MatrixXd derivmat(int n, double t) {
	Eigen::MatrixXd mat1(n, n+1), mat2(n, n+1);
	mat1 << Eigen::MatrixXd::Identity(n, n), Eigen::MatrixXd::Zero(n, 1);
	mat2 << Eigen::MatrixXd::Zero(n, 1), Eigen::MatrixXd::Identity(n, n);
	return (n/t)*(mat2 - mat1);
}


/* Performs derivation of control points p with final time t*/
Eigen::MatrixXd deriv(Eigen::MatrixXd p, double t) {
	return derivmat(p.rows()-1, t)*p;
}


/*Returns a matrix to perform derivation but preseves order*/
Eigen::MatrixXd deivelevmat(int n, double t){
	return derivmat(n+1, t)*degrelevmat(n, n+1);
}


/*Returns a matrix to perform evaluation on times times*/
Eigen::MatrixXd evalmat(int n, double t, Eigen::MatrixXd times){
	Eigen::MatrixXd mat(times.rows(), n+1);
	for (int i=0; i < mat.rows(); i++) {
		for (int j=0; j < n+1; j++) {
			mat(i, j) = basis(j, n, times(i, 0)/t);
		}
	}
	return mat;
}

Eigen::MatrixXd evalmat(int n, double tf, double t){
	Eigen::MatrixXd mat(1, n+1);
	for (int j=0; j < n+1; j++) {
		mat(0, j) = basis(j, n, t/tf);
	}
	return mat;
}

/*Performs evaluation on times times*/
Eigen::MatrixXd eval(Eigen::MatrixXd p, double tf, Eigen::MatrixXd times) {
	return evalmat(p.rows()-1, tf, times)*p;
}

Eigen::MatrixXd eval(Eigen::MatrixXd p, double tf, double t) {
	return evalmat(p.rows()-1, tf, t)*p;
}

/*Calculates the integral*/
Eigen::MatrixXd integr(Eigen::MatrixXd p, double t) {
	return (t / p.rows()) * p.colwise().sum();
}


/*Control points for multiplication*/
Eigen::MatrixXd mul(Eigen::MatrixXd p1, Eigen::MatrixXd p2) {
	int m = p1.rows()-1, n = p2.rows()-1;
	Eigen::MatrixXd new_p = Eigen::MatrixXd::Zero(m+n+1, p1.cols());
	for (int i = 0; i < m+n+1; i++) {
		for (int j=std::max(0, i-n); j < std::min(m, i)+1; j++) {
			new_p.row(i).array() += comb<double>(i, j) * comb<double>(m+n-i,m-j)* p1.row(j).array()*p2.row(i-j).array() / comb<double>(m+n,n);
		}
	}
	return new_p;
}


/*Control points for power*/
Eigen::MatrixXd pow(Eigen::MatrixXd p, int y) {
	if ( y == 0 ) {
		return Eigen::MatrixXd::Ones(1, p.cols());
	}
	Eigen::MatrixXd temp_p = pow(p, y/2);
	if (y % 2 == 0) {
		return mul(temp_p, temp_p);
	} else {
		return mul(p, mul(temp_p, temp_p));
	}
}


/*Calulates the addition of two polynomials*/
Eigen::MatrixXd add(Eigen::MatrixXd p1, Eigen::MatrixXd p2) {
	if (p1.rows() > p2.rows()) {
		p2 = degrelev(p2, p1.rows()-1);
	} else if (p2.rows() > p1.rows()) {
		p1 = degrelev(p1, p2.rows()-1);
	}
	return p1 + p2;
}



/*Calculates a matrix to convert the control points to monomial polynomails*/
// TODO: Didn't do the tomonmat functions
Eigen::MatrixXd tomonmat(Eigen::MatrixXd n, double t) {
	return Eigen::MatrixXd(1, 1);
	
}

};

// vim:set tabstop=4

