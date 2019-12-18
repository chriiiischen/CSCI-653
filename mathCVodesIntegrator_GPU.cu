//! @file CVodesIntegrator.cpp

#include "mathCVodesIntegrator.h"

#include <iostream>
using namespace std;



// XG: Type of tolorance
#define CV_SS 1  // all the equations have the same tolerance
#define CV_SV 2  // every equation has its own tolerance

typedef long int sd_size_t;

namespace math
{
/**
 * Function called by cvodes to evaluate ydot given y.  The CVODE integrator
 * allows passing in a void* pointer to access external data.
 */
double *dev_transport, *dev_molWeights, *dev_density;
N_Vector m_cache_netProdRates=NULL;

__global__ void cvodes_rhs_cudaKernel(double* molWeights,double* ydotPointer,
                                      double* transport,double* density, int nEq)
{
  int iId = blockDim.x * blockIdx.x + threadIdx.x;

  if(iId < nEq)
  {
    ydotPointer[iId] = (molWeights[iId]*ydotPointer[iId]+transport[iId])/(*density) ;
  } 

  return ;
        
}

__global__ void cvodes_rhs_cudaKernel2(double* molWeights, double* ydotPointerDst,
				       double* ydotPointerSrc,
                                      double* transport,double* density, int nEq)
{
  int iId = blockDim.x * blockIdx.x + threadIdx.x;

  if(iId < nEq)
  {
    ydotPointerDst[iId] = (molWeights[iId]*ydotPointerSrc[iId]+transport[iId])/(*density) ;
  }

  return ;

}


__global__ void vectorReduce(double* pVectorArray, int iNEq)
{
  pVectorArray[iNEq-1]=1.0;
  for(int i=0 ; i<iNEq-1 ; i++){
          pVectorArray[iNEq-1] -= pVectorArray[i] ;
  }
  return;
}

void setDevMem(double **pTarget, int iSize)
{
    cudaError_t cudaStatus;
    if(pTarget!=NULL)
    {
        cudaFree(pTarget);
    }
    cudaStatus = cudaMalloc((void**)pTarget, sizeof(double)*iSize);
    //cout<<pTarget<<' '<<endl;
 		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
                        exit(0); 
		}
    return;
}

void host2dev(double *pDst, double *pSrc, int iSize)
{
   cudaError_t cudaStatus;
    cudaStatus = cudaMemcpy(pDst,pSrc, sizeof(double) * iSize, cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
                        cout<<"cudaMemcpy error: "<<cudaGetErrorString(cudaStatus)<<" "<<pDst<<endl;
                        exit(0);
		}
    return;

}


static int cvodes_rhs(realtype t, N_Vector y, N_Vector ydot, void* data)//JD: cuda version
{
    //cout<<"cvode_cuda "<<endl;
    std::pair<Cantera::IdealGasMix*, double*>* f = (std::pair<Cantera::IdealGasMix*, double*>*) data;
    Cantera::IdealGasMix* gas = f->first ;
    int nEq = gas->nSpecies() ;
    // XG: first nEq elements are molecular weights
    // XG: second nEq elements are transport terms
    double* transport  = &(f->second[nEq]) ;
    double* molWeights = &(f->second[0]) ;
    //gas->setState_TRY(temperature, density, N_VGetDeviceArrayPointer_Cuda(y));
   // N_VCopyFromDevice_Cuda(y);   
   
    //double* ydotPointer_h =  N_VGetHostArrayPointer_Cuda(m_cache_netProdRates) ;
    //gas->getNetProductionRates(ydotPointer_h) ;
    //N_VCopyToDevice_Cuda(m_cache_netProdRates);
    
    //ydotPointer_h =  N_VGetHostArrayPointer_Cuda(ydot) ;
    //gas->getNetProductionRates(ydotPointer_h) ;
    //N_VCopyToDevice_Cuda(ydot);
    
    double * pm_cache_netProdRates = N_VGetDeviceArrayPointer_Cuda(m_cache_netProdRates);
    double* ydotPointer_dev =  N_VGetDeviceArrayPointer_Cuda(ydot) ;
    //cudaMemcpy(ydotPointer_dev, pm_cache_netProdRates, sizeof(double)*nEq, cudaMemcpyDeviceToDevice);
    //double dDensity = gas->density();

    unsigned block = 32;
    unsigned grid = (nEq + block - 1)/block;

    // prepare data for device
    // init device memory
    //double *dev_transport, *dev_molWeights, *dev_density;
    /*cudaError_t cudaStatus;
    cudaStatus = cudaMalloc((void**)&dev_transport, sizeof(double)*nEq);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
		}
    cudaStatus = cudaMalloc((void**)&dev_molWeights,sizeof(double)*nEq);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
		}
    cudaStatus = cudaMalloc((void**)&dev_density,sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
		}

    //copy data to device
    cudaStatus = cudaMemcpy(dev_transport,transport, sizeof(double) * nEq, cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
		}

    cudaStatus = cudaMemcpy(dev_molWeights, molWeights, sizeof(double) * nEq, cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
		}

    cudaStatus = cudaMemcpy(dev_density, &dDensity, sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
		}

    */

    //cvodes_rhs_cudaKernel<<<grid,block>>>(dev_molWeights,ydotPointer_dev,
    //                                   dev_transport, dev_density, nEq);
     
    cvodes_rhs_cudaKernel2<<<grid,block>>>(dev_molWeights,ydotPointer_dev,
 					pm_cache_netProdRates,
                                       dev_transport, dev_density, nEq);
   
    //N_VCopyFromDevice_Cuda(ydot);
    //cudaDeviceSynchronize();	
    //cudaFree(dev_transport);
    //cudaFree(dev_molWeights);
    //cudaFree(dev_density);
    return 0;

    
}

int test_cnt=0;

static int cvodes_rhs_ori(realtype t, N_Vector y, N_Vector ydot, void* data)

{
    cout<<"ori";
    N_VCopyFromDevice_Cuda(y);
    N_VCopyFromDevice_Cuda(ydot);
    test_cnt++;
    std::pair<Cantera::IdealGasMix*, double*>* f = (std::pair<Cantera::IdealGasMix*, double*>*) data;
    Cantera::IdealGasMix* gas = f->first ;
    int nEq = gas->nSpecies() ;
    //double temperature = gas->temperature() ;
    //double density = gas->density() ;
    // XG: first nEq elements are molecular weights
    // XG: second nEq elements are transport terms
    double* transport  = &(f->second[nEq]) ;
    double* molWeights = &(f->second[0]) ;
    //gas->setState_TRY(temperature, density, N_VGetDeviceArrayPointer_Cuda(y));
    double* ydotPointer =  N_VGetHostArrayPointer_Cuda(ydot) ;
    /*cout<<"0: ";
    for(int m=0; m<nEq;m++)
    {
	cout<<ydotPointer[m]<<' ';
    }
    cout<<endl;*/
    gas->getNetProductionRates(ydotPointer) ;
    /*cout<<"1: ";
    for(int m=0; m<nEq;m++)
    {
        cout<<ydotPointer[m]<<' ';
    }
    cout<<endl;*/ 
    //cout<<"2: ";
    for(int m=0 ; m<nEq ; m++){
        ydotPointer[m] = (molWeights[m]*ydotPointer[m]+transport[m])/gas->density() ;
   	//cout<<ydotPointer[m]<<' ';
    }
    //cout<<endl;
    
    N_VCopyToDevice_Cuda(ydot);
    //cout<<"3: ";test_printVector(y,nEq);
   /*double *temp =  N_VGetHostArrayPointer_Cuda(ydot);
    cout<<"3: ";
    for(int m=0; m <nEq ; m++)
    {
       cout<<temp[m]<<' ';
    }
    cout<<endl;*/
    N_VCopyToDevice_Cuda(ydot);
    if(test_cnt==2){
    //exit(0);
    }
    return 0 ; 
}

static int cvodes_jac(realtype t, N_Vector u, N_Vector fu, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    int flag = SUNMatZero(J) ;
    return 0 ;
}

CVodesIntegrator::CVodesIntegrator() :
    m_neq(0),
    m_cvode_mem(0),
    m_t0(0.0),
    //m_y(0),i
    //m_cache_netProdRates(NULL),
    m_abstol(0),
    m_type(DENSE+NOJAC),
    m_itol(CV_SS),
    m_method(CV_BDF),
    m_maxord(0),
    m_reltol(1.e-9),
    m_abstols(1.e-15),
    m_nabs(0),
    m_hmax(0.0),
    m_hmin(0.0),
    m_maxsteps(20000),
    m_maxErrTestFails(0),
    m_mupper(1), m_mlower(1),
    m_A(NULL),
    m_LS(NULL),
    gas(NULL),
    //gas_user(NULL),
    m_temperature(0.0),
    m_density(0.0),
    m_intEnergy(0.0),
    m_userData(NULL),
    m_time_old(0.0),
    max_rhseval(0)
{
    m_y = NULL;
    m_test_Log = new CSciLog("RunTimeRecord");
}

CVodesIntegrator::~CVodesIntegrator()
{

    if (this->m_cvode_mem) {
        CVodeFree(&this->m_cvode_mem);
    }

    if (m_y) {
        N_VDestroy_Cuda(m_y);
    }
    if (m_abstol) {
        N_VDestroy_Serial(m_abstol);
    }
    
    if(m_LS)
    {
      /* Free the linear solver memory */
      SUNLinSolFree(m_LS);
    }

    /* Free the matrix memory */
    
    if(this->m_A)
    {
      SUNMatDestroy(this->m_A);
    }

    if (this->m_userData != NULL)
    {
      delete [] this->m_userData ;
      this->m_userData = NULL ;
    }

    //if (this->gas_user != NULL)
    //{
    //  delete this->gas_user;
    //  this->gas_user = NULL;
    //}

    //ivanComment this->gas should probably not be set to NULL here, as there
    //ivanComment should be another Delete() method from the object who
    //ivanComment instantiated it responsible for deleting it (if not NULL),
    //ivanComment because this->gas is passed at initialization, not created
    //ivanComment within this class
    //gas = NULL ;
}

double& CVodesIntegrator::solution(size_t k)
{
    realtype* pm_y_host = N_VGetHostArrayPointer_Cuda(m_y);
    return pm_y_host[k];
    //return NV_Ith_S(m_y, k);
}

double* CVodesIntegrator::solution()
{
    return N_VGetHostArrayPointer_Cuda(m_y);
}

void CVodesIntegrator::setTolerances(double reltol, size_t n, double* abstol)
{
    m_itol = CV_SV;
    m_nabs = n;
    if (n != m_neq) {
        if (m_abstol) {
            N_VDestroy_Serial(m_abstol);
        }
        m_abstol = N_VNew_Serial(static_cast<sd_size_t>(n));
    }

    //realtype* pm_abstol_host = N_VGetHostArrayPointer_Cuda(m_abstol);

    for (size_t i=0; i<n; i++) {
        NV_Ith_S(m_abstol, i) = abstol[i];
        //pm_abstol_host[i] = abstol[i];
    }
    m_reltol = reltol;
}

void CVodesIntegrator::setTolerances(double reltol, double abstol)
{
    m_itol = CV_SS;
    m_reltol = reltol;
    m_abstols = abstol;
}

void CVodesIntegrator::setProblemType(int probtype)
{
    m_type = probtype;
}

void CVodesIntegrator::setMethod(MethodType t)
{
    if (t == BDF_Method)
    {
        m_method = CV_BDF;
    }
    else if (t == Adams_Method)
    {
        m_method = CV_ADAMS;
    }
    else
    {
        common::errorHandler->SendError("Unknow integration method for CVodesIntegrator");
        common::errorHandler->AbortExecution();
    }
}

void CVodesIntegrator::setMaxStepSize(double hmax)
{
    m_hmax = hmax;
    if (this->m_cvode_mem) {
        CVodeSetMaxStep(this->m_cvode_mem, hmax);
    }
}

void CVodesIntegrator::setMinStepSize(double hmin)
{
    m_hmin = hmin;
    if (this->m_cvode_mem) {
        CVodeSetMinStep(this->m_cvode_mem, hmin);
    }
}

void CVodesIntegrator::setMaxSteps(int nmax)
{
    m_maxsteps = nmax;
    if (this->m_cvode_mem) {
        CVodeSetMaxNumSteps(this->m_cvode_mem, m_maxsteps);
    }
}

void CVodesIntegrator::setMaxErrTestFails(int n)
{
    m_maxErrTestFails = n;
    if (this->m_cvode_mem) {
        CVodeSetMaxErrTestFails(this->m_cvode_mem, n);
    }
}

void CVodesIntegrator::initialize(double t0, Cantera::IdealGasMix *gas,
                                  const std::vector<double> molweights, 
                                  const std::vector<double> transport)
{   
    this->gas = gas;
    //this->gas_user = new Cantera::IdealGasMix(gas[0]);
    m_neq = gas->nSpecies();
    m_temperature = gas->temperature() ;
    m_density = gas->density() ;
    m_intEnergy = gas->intEnergy_mass() ;
    m_t0 = t0;
    m_time = t0;

    if (m_userData != NULL)
    {
      delete [] m_userData;
      m_userData = NULL;
    }

    int size = static_cast<int>(m_neq) ;
    m_userData = new double [2*size] ; 

    if(molweights.size()!=size or transport.size()!=size)
    {
      common::errorHandler->SendError("CVodesIntegrator::initialize failed, molweights or transport do not have the right size");
      common::errorHandler->AbortExecution();
    }  
    else
    {
      for(int i=0 ; i<size ; i++)
      {
        m_userData[i] = molweights[i] ;
        m_userData[i+size] = transport[i] ;
      }
    }

    if (m_y)
    {
        N_VDestroy(m_y); // free solution std::vector if already allocated
    }
    if(m_cache_netProdRates)
    {
 	N_VDestroy(m_cache_netProdRates);
    }

    double * testDev;

    int codeTemp;
    
    cudaError_t cudaStatus = cudaGetDevice(&codeTemp);
    if (cudaStatus != cudaSuccess) {
        cout<<cudaGetErrorString(cudaStatus)<<endl;
        fprintf(stderr, "cudaGetdevice failed!");
    } 	
    else
    {
	cout<<"using device: "<<codeTemp<<endl;
    }

    std::cout<<"prepared to N_VNewManaged: "<<m_neq<<std::endl;
    const realtype neq = m_neq;
    
    m_cache_netProdRates = N_VNew_Cuda(neq);
   
    m_y = N_VNew_Cuda(neq); // allocate solution std::vector
    if(check_retval((void*)m_y, "N_VNew_Cuda", 0))
    {
         std::cout<<"managed memory failed"<<std::endl;
    }
    else
    {
	cout<<"managed memory allocated success"<<endl;
    }
    gas->getMassFractions(N_VGetHostArrayPointer_Cuda(m_y)) ;
    N_VCopyToDevice_Cuda(m_y);    
    // check abs tolerance array size
    if (m_itol == CV_SV && m_nabs < m_neq) {
        common::errorHandler->SendError("not enough tolerances given for option CV_SV");
        common::errorHandler->AbortExecution();
    }

    if (this->m_cvode_mem) {
        CVodeFree(&this->m_cvode_mem);
    }

    //! Specify the method and the iteration type. Cantera Defaults:
    //!        CV_BDF  - Use BDF methods
    //!        CV_NEWTON - use Newton's method
    this->m_cvode_mem = CVodeCreate(m_method);
    if (!this->m_cvode_mem)
    {
        common::errorHandler->SendError("CVodeCreate failed");
        common::errorHandler->AbortExecution();
    }

    int flag = CVodeInit(this->m_cvode_mem, cvodes_rhs, m_t0, m_y);
    if (flag != CV_SUCCESS) {
        if (flag == CV_MEM_FAIL)
        {
          common::errorHandler->SendError("CVodeInit memory allocation failed");
          common::errorHandler->AbortExecution();
        }
        else if (flag == CV_ILL_INPUT)
        {
          common::errorHandler->SendError("Wrong input argument for CVodeInit");
          common::errorHandler->AbortExecution();
        }
        else 
        {
          common::errorHandler->SendError("CVodeInit failed");
          common::errorHandler->AbortExecution();
        }
    }

    if (m_itol == CV_SV)
    {
        flag = CVodeSVtolerances(this->m_cvode_mem, m_reltol, m_abstol);
    }
    else
    {
        flag = CVodeSStolerances(this->m_cvode_mem, m_reltol, m_abstols);
    }
    if (flag != CV_SUCCESS)
    {
        if (flag == CV_MEM_FAIL)
        {
          common::errorHandler->SendError("CVodeSStolerances memory allocation failed");
          common::errorHandler->AbortExecution();
        }
        else if (flag == CV_ILL_INPUT)
        {
          common::errorHandler->SendError("Wrong input argument for CVodeSStolerances");
          common::errorHandler->AbortExecution();
        }
        else 
        {
          common::errorHandler->SendError("CVodeSStolerances failed");
          common::errorHandler->AbortExecution();
        }
    }

    //this->gas_user->setState_TRY(m_temperature, m_density, this->solution());
    
    m_userData_pair = std::make_pair(this->gas,this->m_userData) ;
   
    //std::pair<Cantera::IdealGasMix*, double*>* jinf
    // = (std::pair<Cantera::IdealGasMix*, double*>*) &m_userData_pair;
    //Cantera::IdealGasMix* jingas = jinf->first ;
    //int jinnEq = jingas->nSpecies() ;
    // XG: first nEq elements are molecular weights
    // XG: second nEq elements are transport terms
    //double* jintransport  = &(jinf->second[jinnEq]) ;
    //double* jinmolWeights = &(jinf->second[0]) ;

    setDevMem(&dev_transport,m_neq);
    setDevMem(&dev_molWeights,m_neq);
    setDevMem(&dev_density,1);
    //setDevMem(&dev_netProdRates,m_neq);
    host2dev(dev_transport,&(this->m_userData[m_neq]),m_neq);
    host2dev(dev_molWeights,&(this->m_userData[0]),m_neq);
    double dDensity = (this->gas)->density();
    host2dev(dev_density,&dDensity,1);    
    flag = CVodeSetUserData(this->m_cvode_mem, &m_userData_pair);
    if (flag != CV_SUCCESS)
    {
      common::errorHandler->SendError("CVodeSetUserData failed");
      common::errorHandler->AbortExecution();
    }
    applyOptions();
    cout<<"initialize finished"<<endl;

}



void CVodesIntegrator::reinitialize(double t0, Cantera::IdealGasMix *gas,
                                    const std::vector<double> molweights,
                                    const std::vector<double> transport)
{
    m_t0 = t0;
    m_time = t0;
    this->gas = gas ;
    m_temperature = gas->temperature() ;
    m_density = gas->density() ;
    m_intEnergy = gas->intEnergy_mass() ;

    if(m_userData!=NULL)
    {
      delete [] m_userData ;
      m_userData = NULL ;
    }

    int size = static_cast<int>(m_neq) ;
    m_userData = new double [2*size] ; 

    if(molweights.size()!=size or transport.size()!=size)
    {
      common::errorHandler->SendError("CVodesIntegrator::reinitialize failed, molweights or transport do not have the right size");
      common::errorHandler->AbortExecution();
    }  
    else
    {
      for(int i=0 ; i<size ; i++)
      {
        m_userData[i] = molweights[i] ;
        m_userData[i+size] = transport[i] ;
      }
    }

    gas->getMassFractions(N_VGetHostArrayPointer_Cuda(m_y)) ;
    N_VCopyToDevice_Cuda(m_y);
    //cout<<"4: ";test_printVector(m_y,m_neq);
    int result = CVodeReInit(this->m_cvode_mem, m_t0, m_y);
    if (result != CV_SUCCESS)
    {
        common::errorHandler->SendError("CVodesIntegrator reinitialize failed.");
        common::errorHandler->AbortExecution();
    }

    
    //this->gas_user->setState_TRY(m_temperature, m_density, this->solution());
    m_userData_pair = std::make_pair(this->gas,this->m_userData) ;
    host2dev(dev_transport,&(this->m_userData[m_neq]),m_neq);
    host2dev(dev_molWeights,&(this->m_userData[0]),m_neq);
    double dDensity = (this->gas)->density();
    host2dev(dev_density,&dDensity,1);
    result = CVodeSetUserData(this->m_cvode_mem, &m_userData_pair);
    if (result != CV_SUCCESS)
    {
      common::errorHandler->SendError("CVodeSetUserData failed");
      common::errorHandler->AbortExecution();
    }
}
void test_printVector(N_Vector a, int size)
{
    double * temp = N_VGetHostArrayPointer_Cuda(a);
    for (int i=0;i<size;i++)
    {
	std::cout<<temp[i]<<' ';
    }
    std::cout<<std::endl;
}
void CVodesIntegrator::applyOptions()
{
    int flag = 0 ;
    if (m_type == DENSE + NOJAC or m_type == DENSE + JAC)
    {
        sd_size_t N = static_cast<sd_size_t>(m_neq);

        /* Create dense SUNMatrix for use in linear solves */

        this->m_A = SUNDenseMatrix(N, N);
  
        /* Create dense SUNLinearSolver object for use by CVode */
        m_LS = SUNLinSol_Dense(m_y, this->m_A);
        flag = CVodeSetLinearSolver(this->m_cvode_mem, m_LS, this->m_A);
    } 
    else if (m_type == DIAG + JAC)
    {
        flag = CVDiag(this->m_cvode_mem);
    } 
    else if (m_type == GMRES)
    {
        m_LS = SUNLinSol_SPGMR(m_y, PREC_NONE, 0);
        flag = CVodeSetLinearSolver(this->m_cvode_mem, m_LS, NULL);
    }
    else if (m_type == BAND + NOJAC)
    {
        sd_size_t N = static_cast<sd_size_t>(m_neq);
        long int nu = m_mupper;
        long int nl = m_mlower;

        /* Create banded SUNMatrix for the forward problem */
        this->m_A = SUNBandMatrix(N, nu, nl);

        /* Create banded SUNLinearSolver for the forward problem */
        m_LS = SUNLinSol_Band(m_y, this->m_A);
        flag = CVodeSetLinearSolver(this->m_cvode_mem, m_LS, this->m_A);
    }
    else
    {
      common::errorHandler->SendError("CVodesIntegrator::applyOptions, unsupported option");
      common::errorHandler->AbortExecution();
    }
    
   // m_LS = SUNLinSol_SPGMR(m_y, PREC_NONE, 0);
     //   flag = CVodeSetLinearSolver(this->m_cvode_mem, m_LS, NULL);

    if(flag == CVLS_MEM_NULL)
    {
      common::errorHandler->SendError("CVodesIntegrator::applyOptions, this->m_cvode_mem points to NULL");
      common::errorHandler->AbortExecution();
    }
    else if(flag == CVLS_ILL_INPUT)
    {
      common::errorHandler->SendError("CVodesIntegrator::applyOptions, CVodeSetLinearSolver has wrong input argument types");
      common::errorHandler->AbortExecution();
    }
    else if(flag == CVLS_SUNLS_FAIL)
    {
      common::errorHandler->SendError("CVodesIntegrator::applyOptions, CVodeSetLinearSolver fails to call m_LS");
      common::errorHandler->AbortExecution();
    }
    else if(flag == CVLS_MEM_FAIL)
    {
      common::errorHandler->SendError("CVodesIntegrator::applyOptions, CVodeSetLinearSolver fails to allocate memory");
      common::errorHandler->AbortExecution();
    }
    

    if(m_type == DENSE+NOJAC or m_type == BAND+NOJAC)
      flag = CVodeSetJacFn(this->m_cvode_mem, cvodes_jac);

    if(flag == CVLS_LMEM_NULL)
    {
      common::errorHandler->SendError("CVodesIntegrator::applyOptions, CVLS is not initialized");
      common::errorHandler->AbortExecution();
    }
    else if(flag == CVLS_MEM_NULL)
    {
      common::errorHandler->SendError("CVodesIntegrator::applyOptions, m_cvode_mem points to null");
      common::errorHandler->AbortExecution();
    }



    if (m_maxord > 0) {
        CVodeSetMaxOrd(this->m_cvode_mem, m_maxord);
    }
    if (m_maxsteps > 0) {
        CVodeSetMaxNumSteps(this->m_cvode_mem, m_maxsteps);
    }
    if (m_hmax > 0) {
        CVodeSetMaxStep(this->m_cvode_mem, m_hmax);
    }
    if (m_hmin > 0) {
        CVodeSetMinStep(this->m_cvode_mem, m_hmin);
    }
    if (m_maxErrTestFails > 0) {
        CVodeSetMaxErrTestFails(this->m_cvode_mem, m_maxErrTestFails);
    }
}

void CVodesIntegrator::integrate(double tout)
{
    bool readbaleFlag = false;
    m_test_Log->recordTimeFast("integrate start",true,true);
   //N_VCopyToDevice_Cuda(m_y);
    if (tout == m_time) {
        return;
    }
    while(m_time<tout)
    {
      m_time_old = m_time ;
      m_test_Log->recordTimeFast("start step",false,readbaleFlag);
      double dummy = this->step(tout) ;
      m_test_Log->recordTimeFast("step end",false,readbaleFlag);
      if(m_time<=tout)
      {
        //gas->getMassFractions(N_VGetDeviceArrayPointer_Cuda(m_y)) ;

        // XG: we solve n-1 species in the main function, the nth element of the
        // transport array passed from the main function cannot be trusted
        
        /*N_VCopyFromDevice_Cuda(m_y);// JD add
        double *y = N_VGetHostArrayPointer_Cuda(m_y) ;
        y[m_neq-1] = 1.0 ;
        for(int i=0 ; i<m_neq-1 ; i++)
          y[m_neq-1] -= y[i] ;
        y = NULL;*/
        double *y_dev = N_VGetDeviceArrayPointer_Cuda(m_y);
        vectorReduce<<<1,1>>>(y_dev,m_neq);
        y_dev = NULL ;
        m_test_Log->recordTimeFast("start gas related",false,readbaleFlag);
        N_VCopyFromDevice_Cuda(m_y);// JD add   
        m_test_Log->recordTimeFast("m_y host to device end",false,readbaleFlag);
        this->gas->setState_TRY(m_temperature, m_density, this->solution()) ; // JD bottleneck
        this->gas->setState_UV(m_intEnergy, 1.0/m_density);
        m_temperature = this->gas->temperature() ;
        m_userData_pair = std::make_pair(this->gas,this->m_userData) ;
        m_test_Log->recordTimeFast("gas related end",true,readbaleFlag);
        int flag = CVodeSetUserData(this->m_cvode_mem, &m_userData_pair);
        if (flag != CV_SUCCESS)
        {
          common::errorHandler->SendError("CVodeSetUserData failed");
          common::errorHandler->AbortExecution();
        }
      }
    }
    //N_VCopyToDevice_Cuda(m_y);
    if(m_time>tout)
    {
     // N_VCopyFromDevice_Cuda(m_y);// JD add 
      m_time = m_time_old ;
      gas->getMassFractions(N_VGetHostArrayPointer_Cuda(m_y)) ;
      N_VCopyToDevice_Cuda(m_y); // JD add
      int flag = CVode(this->m_cvode_mem, tout, m_y, &m_time, CV_NORMAL);
      if (flag != CV_SUCCESS)
      {
        common::errorHandler->SendError("CVodeIntegrator integration failed");
        common::errorHandler->AbortExecution();
      }
      // XG: we solve n-1 species in the main function, the nth element of the
      // transport array passed from the main function cannot be trusted
      /*double *y = N_VGetHostArrayPointer_Cuda(m_y) ;
      y[m_neq-1] = 1.0 ;
      for(int i=0 ; i<m_neq-1 ; i++)
        y[m_neq-1] -= y[i] ;
      y = NULL ;*/
      double *y_dev = N_VGetDeviceArrayPointer_Cuda(m_y);
      vectorReduce<<<1,1>>>(y_dev,m_neq);
      y_dev = NULL ;
      N_VCopyFromDevice_Cuda(m_y);// JD add
      this->gas->setState_TRY(m_temperature, m_density, this->solution());
      this->gas->setState_UV(m_intEnergy, 1.0/m_density);
      m_temperature = this->gas->temperature() ;
      m_userData_pair = std::make_pair(this->gas,this->m_userData) ;
      flag = CVodeSetUserData(this->m_cvode_mem, &m_userData_pair);
      if (flag != CV_SUCCESS)
      {
        common::errorHandler->SendError("CVodeSetUserData failed");
        common::errorHandler->AbortExecution();
      }
    }
    if(m_time>tout)
    {
      common::errorHandler->SendError("CVodeIntegrator integration exceeded ending time");
      common::errorHandler->AbortExecution();
    }
    //long int RhsEval = 0 ;
    //CVodeGetNumRhsEvals(this->m_cvode_mem, &RhsEval);
    //if(RhsEval>max_rhseval)
    //{
    //  double dummy = 0.0 ;
    //  CVodeGetCurrentStep(this->m_cvode_mem,&dummy) ;
    //  long int JacEval = 0 ;
    //  CVodeGetNumJacEvals(this->m_cvode_mem,&JacEval) ;
    //  max_rhseval = RhsEval ;
    //  cout<<"number of JacEval is "<<JacEval<<endl ;
    //  cout<<"number of RhsEval is "<<RhsEval<<endl ;
    //  cout<<"current time step is "<<dummy<<endl ;
    //}
    //N_VCopyToDevice_Cuda(m_y);
    //N_VPrint_Cuda(m_y);
    //N_VCopyFromDevice_Cuda(m_y);
    m_test_Log->recordTimeFast("integrate end",true,true);
    m_test_Log->save();
}

double CVodesIntegrator::step(double tout)
{
    //cout<<"5: ";
    //test_printVector(m_y,m_neq);
    
    double* ydotPointer_h =  N_VGetHostArrayPointer_Cuda(m_cache_netProdRates) ;
    this->gas->getNetProductionRates(ydotPointer_h);
  
    N_VCopyToDevice_Cuda(m_cache_netProdRates);
    int flag = CVode(this->m_cvode_mem, tout+1.0, m_y, &m_time, CV_ONE_STEP);
    if (flag != CV_SUCCESS)
    {
      common::errorHandler->SendError("CVodeIntegrator stepping failed");
      common::errorHandler->AbortExecution();
    }
    //cout<<"6: ";
    //test_printVector(m_y,m_neq);
    //exit(0);  N_VCopyFromDevice_Cuda(y);
    // N_VCopyFromDevice_Cuda(m_y);
    return m_time;
}

int CVodesIntegrator::nEvals() const
{
    long int ne;
    CVodeGetNumRhsEvals(this->m_cvode_mem, &ne);
    return ne;
}

static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */

  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  /* Check if retval < 0 */

  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */

  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  return(0);
}

}
