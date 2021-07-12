#include <iostream>
#include <Eigen/Dense>
#include <bernsteinlibcpp.h>

// https://stackoverflow.com/a/33786772


using namespace BernsteinLib;

int main()
{
	Eigen::MatrixXf p(5, 1);
	p << 2, 4, -1, 2, 5;
	Eigen::MatrixXf p2(5,2);
	// p2 << p,  4.0,1.0,7.0,2.0,6;
	p2 << p, p;
	double T = 3;
	Eigen::MatrixXf times(3,1);
	times << 1.2, 1.3, 1.8;
	std::cout << "test the bernstein basis.\n";
	std::cout << basis(2,5,0.7) << std::endl;
	std::cout << "test the degree elevation matrix\n";
	std::cout << degrelevmat(3,7) << std::endl;
	std::cout << "test the anti derivative matrix\n";
	std::cout << antiderivmat(3, 6) << std::endl;
	std::cout << "test the antiderivative\n";
	Eigen::MatrixXf p0(1,1);
	p0 << 0;
	std::cout << antideriv(p, T, p0);
	std::cout << "test degree elevation\n";
	std::cout << degrelev(p, 8) << std::endl;
	std::cout << derivmat(4, 1) << std::endl;
	std::cout << evalmat(p.rows()-1,T,  times) << std::endl;
	std::cout << integr(p2, T) << std::endl;
	std::cout << "test" << std::endl;
	std::cout << mul(p2, p2) << std::endl;
	std::cout << pow(p2, 4) << std::endl;

	Eigen::Vector3d vec = Eigen::MatrixXd::Ones(3,1);
	std::cout << vec << std::endl;
}

