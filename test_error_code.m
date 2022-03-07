#include <iostream>
#include <fstream>
#include <Eigen/Eigen>
#include "mat.hpp"
#include<Eigen/SparseCholesky>
//#include<Eigen/IterativeLinearSolvers>
#include<Eigen/SparseLU>

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double

int mat::init_fvm_mat(void)
{
	double dt, dx, dy, dV, PE;
	double dxw, dxe, dyn, dys;
	double tmp, tmp2;
	int i, j, k, kk;

	std::vector<Trip> trp;
	std::vector<Trip> trp2;
	SparseMatrix<double> A(this->x*this->y,this->x*this->y);
	SparseMatrix<double> B(this->x*this->y,this->x*this->y);
	Eigen::VectorXd matT(this->x*this->y);
	Eigen::VectorXd matT2(this->x*this->y);
        //should modified
	dx = (this->UBCX - this->LBCX) / (double)this->x;
	dy = (this->UBCY - this->LBCY) / (double)this->y;
	dxw = dx;
	dxe = dx;
	dyn = dy;
	dys = dy;
	dV = dx * dy; 
	PE = 0.01;
	dt = 1/(double)this->x*PE/10;
	matT.setZero();
	matT2.setZero();

	for (i=1;i<=this->x;i++)
	{
		for(j = 1;j<=this->y;j++)
		{
			k = (i-1)*this->y + j;
            		// West
			if(i==1)
			{
				kk = (this->x-1) * this->y + j;
			}
			else
			{
	            		kk = (i-2)*this->y + j;
			}
			tmp = dy / dxw;
			tmp2=-dy / dxw;
			trp2.push_back(Trip(k-1,kk-1,tmp));
			trp.push_back(Trip(k-1,kk-1, tmp2));

            		// East
			if(i == this->x)
			{
				kk = j;
			}
			else
			{
	            		kk = i*this->y + j;
			}
			tmp = dy / dxe - dV * PE /dxe;
			tmp2=-dy / dxe + dV * PE /dxe;
			trp2.push_back(Trip(k-1,kk-1,tmp));
			trp.push_back(Trip(k-1,kk-1, tmp2));
            		// North
			if(j != 1)
			{
				kk = (i-1)*this->y + j - 1;
				tmp = dx / dyn - dV * PE / dyn;
				tmp2=-dx / dyn + dV * PE / dyn;
				trp2.push_back(Trip(k-1,kk-1,tmp));
				trp.push_back(Trip(k-1,kk-1,tmp2));
			}
            		//South
			if(j != this->y)
			{
				kk = (i-1)*this->y + j + 1;
				tmp = dx / dys;
				tmp2=-dx / dys;
				trp2.push_back(Trip(k-1,kk-1,tmp));
				trp.push_back(Trip(k-1,kk-1,tmp2));
			}
            		// Present
			if(j == 1)
			{
            			tmp = -(dy/dxe + dy/dxw + dx/dys) + 2*dV/dt + dV*PE/dxe;
	            		tmp2 = (dy/dxe + dy/dxw + dx/dys) + 2*dV/dt - dV*PE/dxe;
			}
			else if(j == this->y )
			{
        	    		tmp = -(dy/dxe + dy/dxw + dx/dyn) + 2*dV/dt + dV*PE/dyn + dV*PE/dxe;
            			tmp2 = (dy/dxe + dy/dxw + dx/dyn) + 2*dV/dt - dV*PE/dyn - dV*PE/dxe;	
			}
			else
			{
	            		tmp = -(dy/dxe + dy/dxw + dx/dyn + dx/dys) + 2*dV/dt + dV*PE/dyn + dV*PE/dxe;
	            		tmp2 = (dy/dxe + dy/dxw + dx/dyn + dx/dys) + 2*dV/dt - dV*PE/dyn - dV*PE/dxe;
			}
			trp2.push_back(Trip(k-1,k-1,tmp));
			trp.push_back(Trip(k-1,k-1,tmp2));
		
		}
	
	}

	A.setFromTriplets(trp2.begin(), trp2.end());
	B.setFromTriplets(trp.begin(), trp.end());

	//Eigen::SimplicialCholesky<SpMat> chol(A);  // performs a Cholesky factorization of A
	SparseLU<SparseMatrix<double> >solver(A);
	
	solver.analyzePattern(A);   // for this step the numerical values of A are not used
	solver.factorize(A);	
//	solver.setTolerance(pow(10, -10));

	for( i=0; i<=1;i++)
	{ 
		for(j=0;j<=this->y-1;j++)
		{
			matT(j) = 1;
			matT(this->x*this->y - this->x+j) = 0;
		}
		matT2 = solver.solve(B*matT);         // use the factorization to solve for the given right hand side
		
		if(ceil((double)i/1) == (double)i/1)
		{
			////std::cout<<B *matT - A*matT2<<"\n";
			std::cout<<B<<"\n";
		}
		matT = matT2;
	}

//	this->matT2 = this->matA.fullPivHouseholderQr().solve(matT); 
	//this->matT2.resize(this->y, this->x);
	return 1;
}


