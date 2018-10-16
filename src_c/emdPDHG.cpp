// An implementation of the EMD example in the following paper:
// Jacobs, M.; Leger, F.; Li, W.; Osher, S. 
// "Solving Large-Scale Optimization Problems with a Convergence Rate Independent of Grid Size."
// CAM report, May 2018.
// This file is mainly authored by Matthew Jacobs @ UCLA (majaco@math.ucla.edu).
// The cases of p=1 and p=inf are edited by Jialin Liu @ UCLA (liujl11@math.ucla.edu).

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctime>
#include <float.h>
#include <stdio.h>
#include <fftw3.h>


typedef struct{
    fftw_plan dctIn;
    fftw_plan dctOut;
    double *kernel;
    double *workspace;
}poisson_solver;


double *create_inverse_laplace_kernel(int n1, int n2){
    double *kernel = (double *)calloc(n1*n2,sizeof(double));
    for(int i=0;i<n2;i++){
        for(int j=0;j<n1;j++){
            double x=M_PI*j/(n1*1.0);
            double y=M_PI*i/(n2*1.0);
            
            double negativeLaplacian=2*n1*n1*(1-cos(x))+2*n2*n2*(1-cos(y));
            if(i>0||j>0){
                kernel[i*n1+j]=-1/(4*n1*n2*negativeLaplacian);
                
            }
            
        }
    }
    return kernel;
}


poisson_solver create_poisson_solver_workspace(int n1, int n2, double * ffttime){
    clock_t b,e;
    b=clock();

    poisson_solver fftps;
    fftps.workspace=(double *)calloc(n1*n2,sizeof(double));
    fftps.kernel=create_inverse_laplace_kernel(n1,n2);
    
    fftps.dctIn=fftw_plan_r2r_2d(n2, n1, fftps.workspace, fftps.workspace,
                                 FFTW_REDFT10, FFTW_REDFT10,
                                 FFTW_MEASURE);
    fftps.dctOut=fftw_plan_r2r_2d(n2, n1, fftps.workspace, fftps.workspace,
                                  FFTW_REDFT01, FFTW_REDFT01,
                                  FFTW_MEASURE);
    
    e=clock();
    *ffttime = (e-b)/((double)CLOCKS_PER_SEC); // record fft set up time
    
    return fftps;
}



void destroy_poisson_solver(poisson_solver fftps){
    free(fftps.kernel);
    free(fftps.workspace);
    fftw_destroy_plan(fftps.dctIn);
    fftw_destroy_plan(fftps.dctOut);
}



void invert_laplacian(poisson_solver fftps, double *phi, int n1, int n2){
    
    memcpy(fftps.workspace,phi,n1*n2*sizeof(double));
    
    fftw_execute(fftps.dctIn);
    
    for(int i=0;i<n1*n2;i++){
        fftps.workspace[i]*=fftps.kernel[i];
    }
    
    fftw_execute(fftps.dctOut);
    
    memcpy(phi,fftps.workspace,n1*n2*sizeof(double));
    
}




double calc_divergence(double *ux, double *uy, double *divergence, int n1, int n2){
    double divTot=0;
    for(int i=0;i<n2-1;i++){
        for(int j=0;j<n1-1;j++){
            
            int xp=(j+1);
            int yp=(i+1);
            
            divergence[i*n1+j]=n1*(ux[i*n1+xp]-ux[i*n1+j])+n2*(uy[yp*n1+j]-uy[i*n1+j]);
            divTot+=pow(divergence[i*n1+j],2);
        }
        divergence[i*n1+n1-1]=-n1*ux[i*n1+n1-1]+n2*(uy[(i+1)*n1+n1-1]-uy[i*n1+n1-1]);
        divTot+=pow(divergence[i*n1+n1-1],2);
    }
    for(int j=0;j<n1-1;j++){
        divergence[(n2-1)*n1+j]=n1*(ux[(n2-1)*n1+j+1]-ux[(n2-1)*n1+j])-n2*uy[(n2-1)*n1+j];
        divTot+=pow(divergence[(n2-1)*n1+j],2);
        
    }
    divergence[n2*n1-1]=-(n1*ux[n1*n2-1]+n2*uy[n1*n2-1]);
    divTot+=pow(divergence[n2*n1-1],2);
    
    divTot/=n1*n2;
    divTot=sqrt(divTot);
    return divTot;
}

void calc_gradient(double *ux, double *uy, double *potential, int n1, int n2){
    
    ux[0]=0;
    memset(uy,0,n1*sizeof(double));
    for(int j=1;j<n1;j++){
        ux[j]=n1*(potential[j]-potential[j-1]);
    }
    
    
    for(int i=1;i<n2;i++){
        ux[n1*i]=0;
        uy[n1*i]=n2*(potential[i*n1]-potential[(i-1)*n1]);
        for(int j=1;j<n1;j++){
            int xm=j-1;
            int ym=i-1;
            ux[i*n1+j]=n1*(potential[i*n1+j]-potential[i*n1+xm]);
            uy[i*n1+j]=n2*(potential[i*n1+j]-potential[ym*n1+j]);
        }
    }
    
}


void proj_l1(double * px, double *py, double dx, double dy){

int case1 = 0;
if (dx + dy < -1.0) case1 = 0;
if (dx + dy <=1.0 && dx + dy >= -1.0) case1 = 1;
if (dx + dy > 1.0) case1 = 2;

int case2 = 0;
if (-dx + dy < -1.0) case2 = 0;
if (-dx + dy <=1.0 && -dx + dy >= -1.0) case2 = 1;
if (-dx + dy > 1.0) case2 = 2;

//area 0
if (case1 == 1 && case2 == 1) {px[0] = dx; py[0] = dy;}
//area 1
if (case1 == 2 && case2 == 1) {px[0] = (dx-dy+1.0)/2.0; py[0] = (dy-dx+1.0)/2.0;}
//area 2
if (case1 == 1 && case2 == 2) {px[0] = -(-dx-dy+1.0)/2.0; py[0] = (dy+dx+1.0)/2.0;}
//area 3
if (case1 == 0 && case2 == 1) {px[0] = -(-dx+dy+1.0)/2.0; py[0] = -(-dy+dx+1.0)/2.0;}
//area 4
if (case1 == 1 && case2 == 0) {px[0] = (dx+dy+1.0)/2.0; py[0] = -(-dy-dx+1.0)/2.0;}
//area 5
if (case1 == 2 && case2 == 2) {px[0] = 0.0; py[0] = 1.0;}
//area 6
if (case1 == 0 && case2 == 2) {px[0] = -1.0; py[0] = 0.0;}
//area 7
if (case1 == 0 && case2 == 0) {px[0] = 0.0; py[0] = -1.0;}
//area 8
if (case1 == 2 && case2 == 0) {px[0] = 1.0; py[0] = 0.0;}
}


double update_p(double *px, double *py, double *ux, double *uy, double *uxold, double *uyold, double *psix, double *psiy,  double sigma, double theta, int n1, int n2, int p){
    double change=0;
    int pcount=n1*n2;

	if(p==1) {


	for(int i=0;i<pcount;i++){
        	double pxold=px[i];
        	double pyold=py[i];      
		double factor=0.0;
        	px[i]+=sigma*(psix[i]+ux[i]+theta*(ux[i]-uxold[i]));
		factor = fmax(fabs(px[i]),1.0);
		px[i]/=factor;
        	py[i]+=sigma*(psiy[i]+uy[i]+theta*(uy[i]-uyold[i]));
        	factor = fmax(fabs(py[i]),1.0);
        	py[i]/=factor;
        	change+=fabs(px[i]-pxold)+fabs(py[i]-pyold);
    	}
    	change/=pcount;

	}
	else if (p == 2) {

	for(int i=0;i<pcount;i++){
        	double pxold=px[i];
        	double pyold=py[i];      
        	px[i]+=sigma*(psix[i]+ux[i]+theta*(ux[i]-uxold[i]));
        	py[i]+=sigma*(psiy[i]+uy[i]+theta*(uy[i]-uyold[i]));        
        	double norm=sqrt(px[i]*px[i]+py[i]*py[i]);        
        	double factor=fmax(norm,1);       
        	px[i]/=factor;
        	py[i]/=factor;
		change+=pow(px[i]-pxold,2)+pow(py[i]-pyold,2);        
    	}
    	change/=pcount;
    	change=sqrt(change);

	}
	else if (p==3) {

	for(int i=0;i<pcount;i++){
       		double pxold=px[i];
        	double pyold=py[i];      
		double factor=0.0;
        	px[i]+=sigma*(psix[i]+ux[i]+theta*(ux[i]-uxold[i]));
		py[i]+=sigma*(psiy[i]+uy[i]+theta*(uy[i]-uyold[i]));        
		proj_l1(px+i,py+i,px[i],py[i]);
        	change+=fmax( fabs(px[i]-pxold), fabs(py[i]-pyold) );
    	}
    	change/=pcount;

	}
    
    return change;
}

void leray_projection(poisson_solver fftps, double *ux, double *uy, double *divergence, int n1, int n2){
    
    int pcount=n1*n2;
    
    double *tempx=(double *)calloc(pcount,sizeof(double));
    double *tempy=(double *)calloc(pcount,sizeof(double));
    
   
    memset(fftps.workspace,0,pcount*sizeof(double));
    
    calc_divergence(ux, uy, fftps.workspace, n1, n2);
    
    fftw_execute(fftps.dctIn);
    
    for(int i=0;i<pcount;i++){
        fftps.workspace[i]*=fftps.kernel[i];
    }
    
    fftw_execute(fftps.dctOut);
    
 
    calc_gradient(tempx, tempy, fftps.workspace , n1, n2);
    
    for(int i=0;i<pcount; i++){
        ux[i]-=tempx[i];
        uy[i]-=tempy[i];
        
    }
    
    free(tempx);
    free(tempy);
    
}




void update_u(poisson_solver fftps, double *ux, double *uy, double *px, double *py, double *divergence, double tau, int n1, int n2){
    int pcount=n1*n2;
    for(int i=0;i<pcount;i++){
        ux[i]-=tau*px[i];
        uy[i]-=tau*py[i];
    }
    leray_projection(fftps, ux, uy, divergence, n1, n2); 
}




void check_divergence(double *density, double *divergence, double *psix, double *psiy, int n1, int n2){
    
    calc_divergence(psix, psiy, divergence, n1, n2);
    double sum=0;
    int pcount=n1*n2;
    for(int i=0;i<pcount;i++){
        sum+=pow(density[i]-divergence[i],2);
    }
    sum/=pcount;
    sum=sqrt(sum);
    
}


void get_grad_psi(poisson_solver fftps, double *density, double *psix, double *psiy, int n1, int n2){
    
    int pcount=n1*n2;
    
    double *potential=(double *)calloc(pcount,sizeof(double));
    
    memcpy(potential, density,pcount*sizeof(double));
    
    invert_laplacian(fftps, potential, n1, n2);
    
    calc_gradient(psix, psiy, potential, n1, n2);
    
    check_divergence(density, potential, psix, psiy, n1, n2);
    
    free(potential);
}


double calc_primal(double *ux, double *uy, double *psix, double *psiy, int n1, int n2, int p){
    
if (p==1) {

    int pcount=n1*n2;
    double value=0;
    for(int i=0;i<pcount;i++){
        double norm=fabs(ux[i]+psix[i])+fabs(uy[i]+psiy[i]);
        value+=norm;
    }
    value/=pcount;
    return value;

}
else if (p==2) {

    int pcount=n1*n2;
    double value=0;
    for(int i=0;i<pcount;i++){
        double norm=sqrt(pow(ux[i]+psix[i],2)+pow(uy[i]+psiy[i],2));
        value+=norm;
    }
    value/=pcount;
    return value;

}
else if (p==3) {

    int pcount=n1*n2;
    double value=0;
    for(int i=0;i<pcount;i++){
        double norm=fmax( fabs(ux[i]+psix[i]) , fabs(uy[i]+psiy[i]) );
        value+=norm;
    }
    value/=pcount;
    return value;
}
    
}

double calc_dual(double *px, double *py, double *psix, double *psiy,  int n1, int n2){
    
    int pcount=n1*n2;
    
    
    double value=0;
    for(int i=0;i<pcount;i++){
        value+=px[i]*psix[i]+py[i]*psiy[i];
        
    }
    value/=pcount;
    return value;
    
}


double calc_residual(double *ux, double *uy, double *uxold, double *uyold, double *px, double *py, double *pxold, double *pyold, double tau, double sigma, int n1, int n2){
   int pcount=n1*n2;
    double res=0;
    for(int i=0;i<pcount;i++){
        res += ((pow(ux[i]-uxold[i],2) + pow(uy[i]-uyold[i],2))/tau);
	res += ((pow(px[i]-pxold[i],2) + pow(py[i]-pyold[i],2))/sigma);
	res += (( (ux[i]-uxold[i])*(px[i]-pxold[i]) )*2.0);
	res += (( (uy[i]-uyold[i])*(py[i]-pyold[i]) )*2.0);
    }
    res /= pcount;
    return res;
}


void emdPDHG(int p, double *density, double tau, double sigma, double error, int maxIters, int n1, int n2, double *ux, double *uy, double *px, double *py, double *psix, double *psiy, double *kpotential, double * iters, double * ffttime, double * caltime){
// int p = 1 means we use 1-norm as the ground metric in EMD calculation.
// p could be 1,2,3. p=3 means sup-norm.
// [ux uy] is the primal flux; [px py] is the gradient of the Kantorovich potential. 
    int pcount=n1*n2;
    poisson_solver fftps=create_poisson_solver_workspace(n1, n2, ffttime);
    
    double *uxold=(double *)calloc(pcount,sizeof(double));
    double *uyold=(double *)calloc(pcount,sizeof(double));

    memcpy(uxold,ux,pcount*sizeof(double));
    memcpy(uyold,uy,pcount*sizeof(double));
    
    double *pxold=(double *)calloc(pcount,sizeof(double));
    double *pyold=(double *)calloc(pcount,sizeof(double));

    memcpy(pxold,px,pcount*sizeof(double));
    memcpy(pyold,py,pcount*sizeof(double));
    
    double *divergence=(double *)calloc(pcount,sizeof(double));

    clock_t b,e;
    b=clock();// start timer
    
    get_grad_psi(fftps, density, psix, psiy, n1, n2);
    
    double theta=1.0;
    
    for(int i=0;i<maxIters;i++){
	update_u(fftps, ux, uy, px, py, divergence, tau, n1, n2);
        update_p(px, py, ux, uy, uxold, uyold, psix, psiy,  sigma, theta, n1, n2, p);

	//double primal=calc_primal(ux, uy, psix, psiy, n1, n2, p);
	//double dual=calc_dual(px, py, psix, psiy, n1, n2);
	double res = calc_residual(ux, uy, uxold, uyold, px, py, pxold, pyold, tau, sigma, n1, n2);

        memcpy(uxold,ux,pcount*sizeof(double));
        memcpy(uyold,uy,pcount*sizeof(double));
	memcpy(pxold,px,pcount*sizeof(double));
        memcpy(pyold,py,pcount*sizeof(double));

	if(fabs(res)<error){//stop
            *iters = (double)i;
            if (i==0) {*iters = 1.0;}
            break;
        }       
    }

    calc_divergence(px, py, kpotential, n1, n2);
    invert_laplacian(fftps, kpotential, n1, n2);

    e=clock();
    *caltime = (e-b)/((double)CLOCKS_PER_SEC); // record calculation time
    
    destroy_poisson_solver(fftps);
    free(uxold);
    free(uyold);
    free(pxold);
    free(pyold);
    free(divergence); 
}
















