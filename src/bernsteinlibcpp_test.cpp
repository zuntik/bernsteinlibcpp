#include <iostream>
#include <Eigen/Dense>
#include <bernsteinlibcpp.h>
#include <boost/math/special_functions/binomial.hpp>

// https://stackoverflow.com/a/33786772


template <typename T>
const auto comb = boost::math::binomial_coefficient<T>;

double basis(int k, int n, double tau) {
	// return boost::math::binomial_coefficient<double>(n, k)*std::pow(1-tau, n-k)*std::pow(tau,k);
	return comb<double>(n, k)*std::pow(1-tau, n-k)*std::pow(tau,k);
}


/* Returns a matrix for degree elevation */
Eigen::MatrixXf degrelevmat(int n, int m) {
	Eigen::MatrixXf mat(m+1, n+1);
	for (int i = 0; i < n+1; i++) {
		for (int j = 0; i < m+1; j++) {
			mat(i, j) = (std::max(0, i - m + n) <= j <= std::min(n, i) + 1)? comb<double>(i, j) * comb<double>(m-i, n-j): 0;
		}
	}
	return mat;
}


#if 0

degrelevmat(n, m):
    """Returns a matrix for degree elevation"""
    return np.array([[comb(i, j) * comb(m - i, n - j) if max(0, i - m + n) <= j <= min(n, i) + 1 else 0
                      for j in range(n + 1)] for i in range(m+1)]) / comb(m, n)
#endif


#if 0

/*Returns an integral/antiderivative matrix
 *Returns an integral/antiderivative matrix that must be multiplied by a column vector
 *to obtain the control points of the integral function"""
 */
double antiderivmat(n, t) {
    return np.tril(np.ones((n + 2, n + 1)) * t / (n + 1), -1)
}


double antideriv(p, t, p0) {
    """Returns the control points of the integral function of for the control points of p"""
    return p0 + antiderivmat(p.shape[0] - 1, t) @ p
}


def degrelev(p, m):
    """Performs degree elevation"""
    return p if p.shape[0] - 1 > m else degrelevmat(p.shape[0] - 1, m) @ p


def derivmat(n, t):
    """Returns a matrix to perform derivation"""
    return n / t * (np.hstack((np.zeros((n, 1)), np.eye(n))) - np.hstack((np.eye(n), np.zeros((n, 1)))))


def deriv(p, t):
    """Performs derivation of control points p with final time t"""
    return derivmat(p.shape[0] - 1, t) @ p


def derivelevmat(n, t):
    """Performs derivation while preserving order"""
    return derivmat(n + 1, t) @ degrelevmat(n, n + 1)


def derivelev(p, t):
    """Performs degree elevation"""
    return derivelevmat(p.shape[0] - 1, t) @ p


def evalmat(n, t, times):
    """Returns a matrix to perform evaluation on times times"""
    return np.array([[basis(j, n, ti / t) for j in range(n + 1)] for ti in np.array(times).flatten()])


def eval(p, t, times):
    """Performs evaluation on times times"""
    return evalmat(p.shape[0] - 1, t, times) @ p


def integr(p, t):
    """Calculates the integral"""
    return t / p.shape[0] * np.sum(p, 0)


def mul(p1, p2):
    """Control points for multiplication"""
    if p1.shape[0] < p2.shape[0]:
        p1, p2 = p2, p1
    m, n = p1.shape[0]-1, p2.shape[0]-1
    return np.array([np.sum(
        [comb(i, j) * comb(m + n - i, m - j) * p1[j, :] * p2[i - j, :] for j in range(max(0, i - n), min(m, i) + 1)], 0)
                     for i in range(m + n + 1)]) / comb(m + n, n)


def pow(p, y):
    """Control points for power"""
    if y == 0:
        return np.ones((1, p.shape[1]))
    temp_p = pow(p, y // 2)
    if y % 2 == 0:
        return mul(temp_p, temp_p)
    else:
        return mul(p, mul(temp_p, temp_p))


def add(p1, p2):
    """Calulates the addition of two polynomials"""
    if p1.shape[0] > p2.shape[0]:
        p2 = degrelev(p2, p1.shape[0] - 1)
    elif p2.shape[0] > p1.shape[0]:
        p1 = degrelev(p1, p2.shape[0] - 1)
    return p1 + p2


def tomonmat(n, t):
    """Calculates a matrix to convert the control points to monomial polynomails"""
    return np.flipud([[0 if i > k else comb(n, k) * comb(k, i) * (-1) ** (k - i) for i in range(n + 1)] for k in
                      range(n + 1)]) / np.array([t ** i for i in range(n, -1, -1)]).reshape((-1, 1))


def tomon(p, t):
    """Convers the control points to monomial coefficients"""
    return tomonmat(p.shape[0] - 1, t) @ p
#endif

int main()
{
	sum();
	std::cout << std::pow(3,2) << std::endl;
	std::cout << basis(2,5,0.7) << std::endl;
}
