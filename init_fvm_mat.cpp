#include <iostream>
#include <fstream>
#include <Eigen/Eigen>
#include "mat.hpp"

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double

int mat::init_fvm_mat(void)
{
	double dt, dx, dy;
	double dxw, dxe, dyn, dys;
	double tmp;
	int i, j, k, kk;

	std::vector<Trip> trp, val;
	std::vector<Trip> trp2, val2;
	SparseMatrix<double> A(this->x*this->y,this->x*this->y);
	SparseMatrix<double> B(this->x*this->y,this->x*this->y);
	Eigen::VectorXd matT(this->x*this->y);
	std::string str_i;
	Eigen::VectorXd matC(this->x*this->y);
        //should modified
	dx = (this->UBCX - this->LBCX) / (double)this->x;
	dy = (this->UBCY - this->LBCY) / (double)this->y;
	dt = 0.1;
	dxw = dx;
	dxe = dx;
	dyn = dy;
	dys = dy;
	matT.setZero();
	this->matT2.setZero(this->x*this->y, 1);

	for (i=0;i<=this->x-1;i++)
	{
		for(j = 0;j<=this->y-1;j++)
		{
			k = i*this->y + j;
			if(i == 0)
			{
            			// East
            			kk = (i+1)*this->y + j;
				tmp = dy / dxe;
				trp.push_back(Trip(k,kk,tmp));
            			// North
            			if(j != 0)
				{
					kk = i*this->y + j - 1;
					tmp = dx / dyn;
				trp.push_back(Trip(k,kk,tmp));
				}
            			// South
            			if(j != this->y - 1)
				{
					kk = i*this->y + j + 1;
					tmp = dx / dys;
				trp.push_back(Trip(k,kk,tmp));
				}
            			// Present
            			tmp = -(dy/dxe + dy/dxw + dx/dyn + dx/dys) + 1/dt;
				trp.push_back(Trip(k,k,tmp));
				trp2.push_back(Trip(k,k,1/dt));
			}else if(i == this->x - 1)
			{
            			// West
            			kk = (i-1)*this->y + j;
				tmp = dy / dxw;
				trp.push_back(Trip(k,kk,tmp));
            			// North
            			if(j != 0)
				{
					kk = i*this->y + j - 1;
					tmp = dx / dyn;
				trp.push_back(Trip(k,kk,tmp));
				}
            			// South
            			if(j != this->y - 1)
				{
					kk = i*this->y + j + 1;
					tmp = dx / dys;
				trp.push_back(Trip(k,kk,tmp));
				}
            			// Present
            			tmp = -(dy/dxe + dy/dxw + dx/dyn + dx/dys) + 1/dt;
				trp.push_back(Trip(k,k,tmp));
				trp2.push_back(Trip(k,k,1/dt));
			}else if(j == 0)
			{
            			// West
            			if(i != 0)
				{
            				kk = (i-1)*this->y + j;
					tmp = dy / dxw;
				trp.push_back(Trip(k,kk,tmp));
				}
            			// East
            			if(i != this->x - 1)
				{
            				kk = (i+1)*this->y + j;
					tmp = dy / dxe;
				trp.push_back(Trip(k,kk,tmp));
				}
            			// North
					matC(k) = dx/dyn*sin((double)pi*((double)i/(double)this->x));
					matT(k) = sin(pi*(double)i/(double)this->x);
            			// South
					kk = i*this->y + j + 1;
					tmp = dx / dys;
				trp.push_back(Trip(k,kk,tmp));
            			// Present
            			tmp = -(dy/dxe + dy/dxw + dx/dyn + dx/dys) + 1/dt;
				trp.push_back(Trip(k,k,tmp));
				trp2.push_back(Trip(k,k,1/dt));
			}else if(j == this->y - 1)
			{
            			// West
            			if(i != 0)
				{
	            			kk = (i-1)*this->y + j;
					tmp = dy / dxw;
				trp.push_back(Trip(k,kk,tmp));
				}
            			// East
            			if(i != this->x - 1)
				{
	            			kk = (i+1)*this->y + j;
					tmp = dy / dxe;
				trp.push_back(Trip(k,kk,tmp));
				}
            			// North
					kk = i*this->y + j - 1;
					tmp = dx / dyn;
					trp.push_back(Trip(k,kk,tmp));
            			// Present
            			tmp = -(dy/dxe + dy/dxw + dx/dyn + dx/dys) + 1/dt;
				trp.push_back(Trip(k,k,tmp));
				trp2.push_back(Trip(k,k,1/dt));
			}else
			{
            			// West
            			kk = (i-1)*this->y + j;
				tmp = dy / dxw;
				trp.push_back(Trip(k,kk,tmp));
            			// East
            			kk = (i+1)*this->y + j;
				tmp = dy / dxe;
				trp.push_back(Trip(k,kk,tmp));
            			// North
				kk = i*this->y + j - 1;
				tmp = dx / dyn;
				trp.push_back(Trip(k,kk,tmp));
            			// South
				kk = i*this->y + j + 1;
				tmp = dx / dys;
				trp.push_back(Trip(k,kk,tmp));
            			// Present
            			tmp = -(dy/dxe + dy/dxw + dx/dyn + dx/dys) + 1/dt;
				trp.push_back(Trip(k,k,tmp));
				trp2.push_back(Trip(k,k,1/dt));
			}
		
		}
	
	}

	
	A.setFromTriplets(trp.begin(), trp.end());
	B.setFromTriplets(trp2.begin(), trp2.end());

	Eigen::SimplicialCholesky<SpMat> chol(B);  // performs a Cholesky factorization of A
	Eigen::SimplicialCholesky<SpMat> chol2(B);  // performs a Cholesky factorization of A
	for( i=0; i<=2000;i++)
	{ 
		matT = chol.solve(A*matT+matC);         // use the factorization to solve for the given right hand side
		matT = chol.solve(A*matT+matC);         // use the factorization to solve for the given right hand side
		if(ceil((double)i/100) == (double)i/100)
		{
			std::cout<<matT<<"\n";
			//str_i << i;
			//this->print_fvm_result();
		}
	}

//	this->matT2 = this->matA.fullPivHouseholderQr().solve(matT); 
	//this->matT2.resize(this->y, this->x);
	return 1;
}


