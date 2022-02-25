
#include <iostream>
#include <fstream>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <unsupported/Eigen/Splines>
#include <vector>
using namespace Eigen;

//typedef Eigen::Spline<double, 1, 2> Spline1D;
//typedef Eigen::SplineFitting<Spline1D> SplineFitting1D;

//sparse
typedef Triplet<int> Trip;

class mat 
{
	public:
		double LBCY, UBCY, LBCX, UBCX, dx, dy;
		int x, y, edgenum;
		int init_var(double lx, double ux, double ly, double uy, int xx, int yy)
		{this->LBCY=ly;this->LBCX=lx;this->UBCY=uy;this->UBCX=ux;this->x=xx;this->y=yy;this->dx=(UBCX-LBCX)/(double)xx;this->dy=(UBCY-LBCY)/(double)yy;return 1;}
		double interp(double r1, double f1, double r2, double f2, double r3);
		int init_mat(void);
		int init_fvm_mat(void);
		int print_fvm_result(void);
		double pi=3.14159265358979;
		// center
		Eigen::Matrix<double, Dynamic, Dynamic> f1;
		Eigen::Matrix<double, Dynamic, Dynamic> x1;
		Eigen::Matrix<double, Dynamic, Dynamic> y1;

		// 
		Eigen::Matrix<double, Dynamic, Dynamic> matT2;
		//sparse 

};

