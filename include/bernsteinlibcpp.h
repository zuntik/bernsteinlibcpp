#pragma once

namespace BernsteinLib {

	double basis(int k, int n, double tau);
	Eigen::MatrixXd degrelevmat(int n, int m);
	Eigen::MatrixXd degrelev(Eigen::MatrixXd p, int m);
	Eigen::MatrixXd antiderivmat(int n, double t);
	Eigen::MatrixXd antideriv(Eigen::MatrixXd p, double t, Eigen::MatrixXd p0);
	Eigen::MatrixXd derivmat(int n, double t);
	Eigen::MatrixXd deriv(Eigen::MatrixXd p, double t);
	Eigen::MatrixXd deivelevmat(int n, double t);
	Eigen::MatrixXd evalmat(int n, double tf, double t);
	Eigen::MatrixXd evalmat(int n, double tf, Eigen::MatrixXd times);
	Eigen::MatrixXd eval(Eigen::MatrixXd p, double t, Eigen::MatrixXd times);
	Eigen::MatrixXd eval(Eigen::MatrixXd p, double tf, double t);
	Eigen::MatrixXd integr(Eigen::MatrixXd p, double t);
	Eigen::MatrixXd mul(Eigen::MatrixXd p1, Eigen::MatrixXd p2);
	Eigen::MatrixXd pow(Eigen::MatrixXd p, int y);
	Eigen::MatrixXd add(Eigen::MatrixXd p1, Eigen::MatrixXd p2);
	Eigen::MatrixXd tomonmat(Eigen::MatrixXd n, double t);

};
