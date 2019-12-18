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
static int cvodes_rhs(realtype t, N_Vector y, N_Vector ydot, void* data)
{
    std::pair<Cantera::IdealGasMix*, double*>* f = (std::pair<Cantera::IdealGasMix*, double*>*) data;
    Cantera::IdealGasMix* gas = f->first ;
    int nEq = gas->nSpecies() ;
    //double temperature = gas->temperature() ;
    //double density = gas->density() ;
    // XG: first nEq elements are molecular weights
    // XG: second nEq elements are transport terms
    double* transport  = &(f->second[nEq]) ;
    double* molWeights = &(f->second[0]) ;
    //gas->setState_TRY(temperature, density, N_VGetArrayPointer(y));
    double* ydotPointer =  N_VGetArrayPointer(ydot) ;
    gas->getNetProductionRates(ydotPointer) ;
    for(int m=0 ; m<nEq ; m++)
        ydotPointer[m] = (molWeights[m]*ydotPointer[m]+transport[m])/gas->density() ;
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
    m_y(0),
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
{}

CVodesIntegrator::~CVodesIntegrator()
{

    if (this->m_cvode_mem) {
        CVodeFree(&this->m_cvode_mem);
    }

    if (m_y) {
        N_VDestroy_Serial(m_y);
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
    return NV_Ith_S(m_y, k);
}

double* CVodesIntegrator::solution()
{
    return N_VGetArrayPointer(m_y);
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
    for (size_t i=0; i<n; i++) {
        NV_Ith_S(m_abstol, i) = abstol[i];
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
        N_VDestroy_Serial(m_y); // free solution std::vector if already allocated
    }
    m_y = N_VNew_Serial(static_cast<sd_size_t>(m_neq)); // allocate solution std::vector

    gas->getMassFractions(N_VGetArrayPointer(m_y)) ;
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
    flag = CVodeSetUserData(this->m_cvode_mem, &m_userData_pair);
    if (flag != CV_SUCCESS)
    {
      common::errorHandler->SendError("CVodeSetUserData failed");
      common::errorHandler->AbortExecution();
    }
    applyOptions();
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

    gas->getMassFractions(N_VGetArrayPointer(m_y)) ;

    int result = CVodeReInit(this->m_cvode_mem, m_t0, m_y);
    if (result != CV_SUCCESS)
    {
        common::errorHandler->SendError("CVodesIntegrator reinitialize failed.");
        common::errorHandler->AbortExecution();
    }

    
    //this->gas_user->setState_TRY(m_temperature, m_density, this->solution());
    m_userData_pair = std::make_pair(this->gas,this->m_userData) ;
    result = CVodeSetUserData(this->m_cvode_mem, &m_userData_pair);
    if (result != CV_SUCCESS)
    {
      common::errorHandler->SendError("CVodeSetUserData failed");
      common::errorHandler->AbortExecution();
    }
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
    if (tout == m_time) {
        return;
    }
    while(m_time<tout)
    {
      m_time_old = m_time ;
      double dummy = this->step(tout) ;
      if(m_time<=tout)
      {
        //gas->getMassFractions(N_VGetArrayPointer(m_y)) ;

        // XG: we solve n-1 species in the main function, the nth element of the
        // transport array passed from the main function cannot be trusted
        double *y = N_VGetArrayPointer(m_y) ;
        y[m_neq-1] = 1.0 ;
        for(int i=0 ; i<m_neq-1 ; i++)
          y[m_neq-1] -= y[i] ;
        y = NULL ;
        this->gas->setState_TRY(m_temperature, m_density, this->solution());
        this->gas->setState_UV(m_intEnergy, 1.0/m_density);
        m_temperature = this->gas->temperature() ;
        m_userData_pair = std::make_pair(this->gas,this->m_userData) ;
        int flag = CVodeSetUserData(this->m_cvode_mem, &m_userData_pair);
        if (flag != CV_SUCCESS)
        {
          common::errorHandler->SendError("CVodeSetUserData failed");
          common::errorHandler->AbortExecution();
        }
      }
    }
    if(m_time>tout)
    {
      m_time = m_time_old ;
      gas->getMassFractions(N_VGetArrayPointer(m_y)) ;
      int flag = CVode(this->m_cvode_mem, tout, m_y, &m_time, CV_NORMAL);
      if (flag != CV_SUCCESS)
      {
        common::errorHandler->SendError("CVodeIntegrator integration failed");
        common::errorHandler->AbortExecution();
      }
      // XG: we solve n-1 species in the main function, the nth element of the
      // transport array passed from the main function cannot be trusted
      double *y = N_VGetArrayPointer(m_y) ;
      y[m_neq-1] = 1.0 ;
      for(int i=0 ; i<m_neq-1 ; i++)
        y[m_neq-1] -= y[i] ;
      y = NULL ;
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
}

double CVodesIntegrator::step(double tout)
{
    int flag = CVode(this->m_cvode_mem, tout+1.0, m_y, &m_time, CV_ONE_STEP);
    if (flag != CV_SUCCESS)
    {
      common::errorHandler->SendError("CVodeIntegrator stepping failed");
      common::errorHandler->AbortExecution();
    }
    return m_time;
}

int CVodesIntegrator::nEvals() const
{
    long int ne;
    CVodeGetNumRhsEvals(this->m_cvode_mem, &ne);
    return ne;
}

}
