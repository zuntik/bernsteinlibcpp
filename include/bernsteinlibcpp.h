#ifndef BERNSTEINLIB_H
#define BERNSTEINLIB_H


namespace BernsteinLib {

	double basis(int k, int n, double tau);
	Eigen::MatrixXf degrelevmat(int n, int m);
	Eigen::MatrixXf degrelev(Eigen::MatrixXf p, int m);
	Eigen::MatrixXf antiderivmat(int n, double t);
	Eigen::MatrixXf antideriv(Eigen::MatrixXf p, double t, Eigen::MatrixXf p0);
	Eigen::MatrixXf derivmat(int n, double t);
	Eigen::MatrixXf deriv(Eigen::MatrixXf p, double t);
	Eigen::MatrixXf deivelevmat(int n, double t);
	Eigen::MatrixXf evalmat(int n, double t, Eigen::MatrixXf times);
	Eigen::MatrixXf eval(Eigen::MatrixXf p, double t, Eigen::MatrixXf times);
	Eigen::MatrixXf integr(Eigen::MatrixXf p, double t);
	Eigen::MatrixXf mul(Eigen::MatrixXf p1, Eigen::MatrixXf p2);
	Eigen::MatrixXf pow(Eigen::MatrixXf p, int y);
	Eigen::MatrixXf add(Eigen::MatrixXf p1, Eigen::MatrixXf p2);
	Eigen::MatrixXf tomonmat(Eigen::MatrixXf n, double t);

};

#endif
