
#include <iostream>
#include <fstream>
#include <Eigen/Eigen>
#include "mat.hpp"
#include <vector>

using namespace Eigen;
typedef Triplet<int> Trip;

int mat::init_mat(void)
{
	int i, j, num=0;
	double tmp1, tmp2;
	std::ofstream centerfile;
	centerfile.open("center.txt");
	f1.setZero(this->x,this->y);
	x1.setZero(this->x,this->y);
	y1.setZero(this->x,this->y);
	//std::cout<<"SIZEX:"<<this->x<<',';
	//std::cout<<"SIZEY:"<<this->y<<',';
	for(i=0;i<=this->x-1;i++)
	{
		for(j=0;j<=this->y-1;j++)
			{
			num++;
			tmp1 = this->LBCX + (double)i * this->dx; 
			tmp2 = this->LBCY + (double)j * this->dy; 
                        x1(i,j) = tmp1;
                        y1(i,j) = tmp2;
			f1(i,j) = sin(pi*tmp1)*sin(pi*tmp2);
			centerfile << tmp1 << " " << tmp2 << " " <<f1(i,j) << "\n";
			}
	}
	centerfile.close();
	return 1;
}
int mat::print_fvm_result(void)
{
	int i, j;
	//const char* filename="fvm_result";
	//std::ofstream file((filename+num+".txt").c_str());
	for(i=0;i<=this->x-1;i++)
	{
		for(j=0;j<=this->y-1;j++)
		{
			std::cout << this->x1(i,j) << " " << this->y1(i,j) << " " << this->matT2(i,j) << "\n";
		}
	}
	//file.close();
	return 1;
}




int main()
{
	int test_op=0;
	double tmp;
	mat test3;
	test3.init_var(0.0,1.0,0.0,1.0, 200, 200);
	test_op = test3.init_mat();
	test_op = test3.init_fvm_mat();
	//std::cout<<test_op;
	return 0;

}
