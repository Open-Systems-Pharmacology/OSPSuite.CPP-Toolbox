/*
 * Integrator.cpp
 *
 */

#include "Integrator.h"

#include<map>

#ifndef PARALLEL
#include "nvector/nvector_serial.h"
inline N_Vector myN_VMake(sunindextype N, realtype *ws) {
    return N_VMake_Serial(N, ws);
}
inline N_Vector* myN_VMakeArray(sunindextype N, N_Vector x, realtype *ws) {
    N_Vector *arr = N_VCloneVectorArrayEmpty_Serial(N, x);
    for(sunindextype i=0; i<N; i++)
        N_VSetArrayPointer_Serial(&ws[i*NV_LENGTH_S(x)], arr[i]);
    return arr;
}
inline void myN_VDestroy(N_Vector x) {
    N_VDestroy_Serial(x);
}
inline void myN_VDestroyArray(int N, N_Vector* arr) {
    N_VDestroyVectorArray_Serial(arr, N);
}
#else
#include "nvector/nvector_pthreads.h"
constexpr int num_threads = 4;
inline N_Vector myN_VMake(int N, double *ws) {
    return N_VMake_Pthreads(N, num_threads, ws);
}
inline void myN_VDestroy(N_Vector x) {
    N_VDestroy_Pthreads(x);
}
#endif

#include "cvodes/cvodes.h"
#include "cvodes/cvodes_direct.h"

#ifdef DER_DENSE
#include "sunmatrix/sunmatrix_dense.h"
#include "sunlinsol/sunlinsol_dense.h"
//int ODEJac_dense (long int N, realtype Time, N_Vector vec_y, N_Vector fy, DlsMat J, void *userdata, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
#ifdef DER_KLU
//#include "sunlinsol/sunlinsol_sparse.h"
#include "sunmatrix/sunmatrix_sparse.h"
#include "sunlinsol/sunlinsol_klu.h"
//int ODEJac_sparse(realtype Time, N_Vector vec_y, N_Vector fy, SlsMat J, void *userdata, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
#ifdef DER_SUPERLU
//#include "sunlinsol/sunlinsol_sparse.h"
#include "sunmatrix/sunmatrix_sparse.h"
#include "sunlinsol/sunlinsol_superlumt.h"
//int ODEJac_sparse(realtype Time, N_Vector vec_y, N_Vector fy, SlsMat J, void *userdata, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif

#ifdef DER_CPPAD
#include "cppad/cppad.h"
int JacDenseCppAD(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
//int JacSparseCppAD(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif

int(*SimInterface::printfPtr)(const char *c,...) = printf;
//unsigned int (*SimInterface::getDimModels)() = nullptr;
//Model* (*SimInterface::getModelPtr)(unsigned int) = nullptr;

struct UserData {
    Model *model;
    double *P;
    unsigned int *S;

    sunindextype *jacCols;
    sunindextype *jacRows;

    // dont call it _P here: some problems with MACRO expansion in Android
    UserData(Model *ptrModel, double *ptrP, unsigned int *ptrS, sunindextype *ptrJacCols, sunindextype *ptrJacRows) : model(ptrModel), P(ptrP), S(ptrS), jacCols(ptrJacCols), jacRows(ptrJacRows) {};
};

int f(realtype t, N_Vector y, N_Vector dy, void *userData) {
    UserData *d = (UserData*)userData;
    d->model->ODERHSFunction( t, N_VGetArrayPointer(y), d->P, d->S, N_VGetArrayPointer(dy) );
    return 0;

    //	bool error = false;
    //	double *ptr = N_VGetArrayPointer(y);
    //	for(int i=0;i<d->model->dimY;i++)
    //		error |= ptr[i] < 0.0;
    //	return error;
}

int g(realtype t, N_Vector y, realtype *gout, void *userData) {
    UserData *d = (UserData*)userData;
    d->model->ODEImplicitSwitch(t, N_VGetArrayPointer(y), d->P, d->S, gout);
    return 0;
}

//int fS(int Ns, realtype t, N_Vector y, N_Vector ydot, N_Vector *yS, N_Vector *ySdot, void *userData, N_Vector tmp1, N_Vector tmp2) {
//	UserData *d = (UserData*)userData;
//
//	d->model->ODEJacSparse   (t, N_VGetArrayPointer(y), d->P, d->S, d->jac_fx->data);
//	d->model->ODEJacSparse_fp(t, N_VGetArrayPointer(y), d->P, d->S, d->jac_fp->data);
//
//	for(int i=0; i<Ns; i++) {
//		double *ptr_ySdot = N_VGetArrayPointer(ySdot[i]);
//		SlsMatvec(d->jac_fx, N_VGetArrayPointer(y), ptr_ySdot);
//		std::fill(ptr_ySdot, ptr_ySdot+d->model->dimY, 0.0);
//		for (int j=d->jac_fp->colptrs[i]; j<d->jac_fp->colptrs[i+1]; j++)
//			ptr_ySdot[d->jac_fp->rowvals[j]] += d->jac_fp->data[j];
//	}
//
//	static int counter = 0;
//	std::cout << "Eval " << ++counter << std::endl;
//	return 0;
//}

//int g_max(realtype t, N_Vector y, realtype *gout, void *userData) {
//	UserData *d = (UserData*)userData;
//}

int JacDense(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *userData, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    UserData *d = (UserData*)userData;
    d->model->ODEJacDense(t, N_VGetArrayPointer(y), d->P, d->S, SM_DATA_D(J));

    return 0;
}

#if defined(DER_KLU) || defined(DER_SUPERLU)
int JacSparse(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *userData, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    UserData *d = (UserData*)userData;
    // this seems to be called with several SlsMat instances, because copying the
    // sparsity pattern only the first time for each model does not work for superLU
    // SimInterface::printfPtr("[DEBUG] JacO1: %e, JacO2: %e\n", SM_INDEXVALS_S(J)[3], SM_INDEXPTRS_S(J)[3]);

    memcpy(SUNSparseMatrix_IndexPointers(J), d->jacCols, (d->model->dimY+1)*sizeof(sunindextype));
    memcpy(SUNSparseMatrix_IndexValues(J), d->jacRows, (d->model->dimNNZ)*sizeof(sunindextype));

    d->model->ODEJacSparse(t, N_VGetArrayPointer(y), d->P, d->S, SUNSparseMatrix_Data(J));
    return 0;
}
#endif

void SimInterface::InitializeMemory() {
    dimT = model->dimOutputTimePoints;

    wsY = new double[model->dimY];
    wsP = new double[model->dimP];
    wsT = new double[model->dimOutputTimePoints];
    //	wsSens = new double[model->dimY*model->dimP_free];

    wsSwitch = new unsigned int[model->dimS];

    jacCols = new sunindextype[model->dimY+1];
    jacRows = new sunindextype[model->dimNNZ];

    userData = (void*) new UserData(model, wsP, wsSwitch, jacCols, jacRows);

    if(model->dimImplicitSwitches)
        isImplicitSwitched = new unsigned int[model->dimImplicitSwitches];
    else
        isImplicitSwitched = nullptr;

    nvecY = myN_VMake( model->dimY, wsY );
    //	nvecArrSens = myN_VMakeArray( model->dimP_free, nvecY, wsSens );

    memcpy(wsP, model->P_init, model->dimP*sizeof(double));
    memcpy(wsT, model->outputTimePoints, model->dimOutputTimePoints*sizeof(double));
}

SimInterface::SimInterface(const SimInterface &s) {

    model = s.model;
    InitializeMemory();

    memcpy(wsY, s.wsY, model->dimY*sizeof(double));
    memcpy(wsP, s.wsP, model->dimP*sizeof(double));
    memcpy(wsT, s.wsT, dimT*sizeof(double));
    //	memcpy(wsSens, s.wsSens, model->dimY*model->dimP_free*sizeof(double));
    memcpy(wsSwitch, s.wsSwitch, model->dimS*sizeof(double));
    memcpy(jacCols, s.jacCols, (model->dimY+1)*sizeof(sunindextype));
    memcpy(jacRows, s.jacRows, (model->dimNNZ)*sizeof(sunindextype));

    absTol   = s.absTol;
    relTol   = s.relTol;
    hInit    = s.hInit;
    hMin     = s.hMin;
    hMax     = s.hMax;
    maxSteps = s.maxSteps;
    maxOrder = s.maxOrder;

    useJac = s.useJac;
    reduceOrder = s.reduceOrder;
    isScaled = s.isScaled;

    jacType = s.jacType;

    progress = 0.0;
    errorFlag = InitializeCVODE();
}

void SimInterface::swap(SimInterface &s) {
    using std::swap;
    swap(wsY, s.wsY);
    swap(wsP, s.wsP);
    swap(wsT, s.wsT);
    swap(wsSens, s.wsSens);
    swap(wsSwitch, s.wsSwitch);
    swap(jacCols, s.jacCols);
    swap(jacRows, s.jacRows);
    swap(dimT, s.dimT);
    swap(isImplicitSwitched, s.isImplicitSwitched);
    swap(nvecY, s.nvecY);
    swap(nvecArrSens, s.nvecArrSens);
    swap(model, s.model);
    swap(cvodeMem, s.cvodeMem);
    swap(userData, s.userData);
    swap(absTol, s.absTol);
    swap(relTol, s.relTol);
    swap(hInit, s.hInit);
    swap(hMin, s.hMin);
    swap(hMax, s.hMax);
    swap(maxSteps, s.maxSteps);
    swap(maxOrder, s.maxOrder);
    swap(useJac, s.useJac);
    swap(reduceOrder, s.reduceOrder);
    swap(isScaled, s.isScaled);
    swap(progress, s.progress);
    swap(errorFlag, s.errorFlag);
    swap(jacType, s.jacType);
}

// copy and swap assignment
SimInterface& SimInterface::operator=(SimInterface s) {
    s.swap(*this);
    return *this;
}

int SimInterface::InitializeCVODE() {

    int cvodeFlag;

    // default parameters that are not set in model file
    maxOrder = 5;
    reduceOrder = false; // false
    progress = 0.0;

    cvodeMem = CVodeCreate(CV_BDF); // , CV_NEWTONinitialize integrator (CV_ADAMS and CV_BDF, CV_FUNCTIONAL and CV_NEWTON)
    assert(cvodeMem!=NULL);

    cvodeFlag = CVodeInit(cvodeMem, f, 0.0, nvecY); // reserve memory
    if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;

    if(model->dimImplicitSwitches) {
        cvodeFlag = CVodeRootInit(cvodeMem, model->dimImplicitSwitches, g); // set implicit root finding function for implicit switches
        if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;
    }

    cvodeFlag = CVodeSetUserData(cvodeMem, userData); // set user data
    if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;

    switch(jacType) {
        case DENSE_FD:
            // only set dense linear solver, obtain derivatives by built-in finite differences
        JMatrix = SUNDenseMatrix(model->dimY, model->dimY);
        LSolver = SUNLinSol_Dense(nvecY, JMatrix); //SUNDenseLinearSolver(nvecY, JMatrix);
        cvodeFlag = CVodeSetLinearSolver(cvodeMem, LSolver, JMatrix);// CVDlsSetLinearSolver
            if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;
            break;
#ifdef DER_DENSE
        case DENSE_CPP:
        JMatrix = SUNDenseMatrix(model->dimY, model->dimY);
        LSolver = SUNDenseLinearSolver(nvecY, JMatrix);
        cvodeFlag = CVodeSetLinearSolver(cvodeMem, LSolver, JMatrix);// CVDlsSetLinearSolver
            if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;
        cvodeFlag = CVodeSetJacFn(cvodeMem, JacDense);
            if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;
            break;
#endif
#ifdef DER_CPPAD
        case DENSE_CPPAD:
            CppadRecord(); // record DAG
        JMatrix = SUNDenseMatrix(model->dimY, model->dimY);
        LSolver = SUNDenseLinearSolver(nvecY, JMatrix);
        cvodeFlag = CVodeSetLinearSolver(cvodeMem, LSolver, JMatrix);// CVDlsSetLinearSolver
            if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;
        cvodeFlag = CVodeSetJacFn(cvodeMem, JacDenseCppAD);
            if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;
            break;
#endif
#ifdef DER_KLU
        case SPARSE_KLU:
        JMatrix = SUNSparseMatrix(model->dimY, model->dimY, model->dimNNZ, CSC_MAT);

        for(unsigned int i=0; i<model->dimY+1; i++)
            jacCols[i] = model->jacCols[i];
        for(unsigned int i=0; i<model->dimNNZ; i++)
            jacRows[i] = model->jacRows[i];

        //memcpy(SM_INDEXVALS_S(JMatrix), model->jacRows, (model->dimNNZ)*sizeof(int));
        //memcpy(SM_INDEXPTRS_S(JMatrix), model->jacCols, (model->dimY+1)*sizeof(int));
        //this->printfPtr("[DEBUG] JacO1: %e, JacO2: %e\n", SM_INDEXVALS_S(JMatrix)[3], SM_INDEXPTRS_S(JMatrix)[3]);
        LSolver = SUNLinSol_KLU(nvecY, JMatrix); // initialize KLU
        cvodeFlag = SUNLinSol_KLUSetOrdering(LSolver, 0); // AMD instead of standard COLAMD
        if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;
        cvodeFlag = CVodeSetLinearSolver(cvodeMem, LSolver, JMatrix);// CVDlsSetLinearSolver
            if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;
        cvodeFlag = CVodeSetJacFn(cvodeMem, JacSparse);
            if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;
            break;
#endif
#ifdef DER_SUPERLU
        case SPARSE_SUPERLU:
            cvodeFlag = CVSuperLUMT(cvodeMem, 1, model->dimY, model->dimNNZ); // initialize SuperLUMT
            if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;
        cvodeFlag = CVodeSetJacFn(cvodeMem, JacSparse);
            if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;
            break;
#endif
        default:
            return -1; // unknown type
    }

    // set up remaining optional options

    // default parameters from model
    model->ODEOptions(relTol, absTol, hInit, hMin, hMax, maxSteps, useJac, model->P_init);

    // check if model is scaled by custom vector
    memcpy(wsY, model->Y_sca, model->dimY*sizeof(double));
    isScaled = false;
    for(unsigned int i=0; i<model->dimY; i++)
        isScaled = isScaled || wsY[i]!=1;

    if(isScaled) {
        for(unsigned int i=0; i<model->dimY; i++) {
            //printfPtr("[DEBUG] %e\n", wsY[i]*absTol);
            wsY[i]*=absTol;
        }
        cvodeFlag = CVodeSVtolerances(cvodeMem, relTol, nvecY);
        printfPtr("[INFO] Using custom abs tol scaling from model file.\n");
    }
    else {
        cvodeFlag = CVodeSStolerances(cvodeMem, relTol, absTol); if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;
    }
//    cvodeFlag = CVodeSStolerances(cvodeMem, relTol, absTol); if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;


    cvodeFlag = CVodeSetMaxNumSteps(cvodeMem, maxSteps);   if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;
    cvodeFlag = CVodeSetMinStep(cvodeMem, hMin);           if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;
    cvodeFlag = CVodeSetInitStep(cvodeMem, hInit);         if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;
    cvodeFlag = CVodeSetMaxStep(cvodeMem, hMax);           if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;
    cvodeFlag = CVodeSetMaxOrd(cvodeMem, maxOrder);        if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;
    cvodeFlag = CVodeSetStabLimDet(cvodeMem, reduceOrder); if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;
    cvodeFlag = CVodeSetMaxErrTestFails(cvodeMem, 15);     if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;


    //      cvodesFlag = CVodeWFtolerances(cvode_mem, ewt); // custom error function
    //      if( cvodesFlag!=CV_SUCCESS ) return cvodesFlag;

    //	// set sensitivity options
    //	cvodeFlag = CVodeSensInit(cvodeMem, model->dimP_free, CV_STAGGERED, fS, nvecArrSens); if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;
    ////	cvodeFlag = CVodeSetSensParams(cvodeMem, ((UserData*)userData)->P, NULL, NULL);         if( cvodeFlag!=CV_SUCCESS ) return cvodeFlag;
    //	std::vector<double> vecAbsTol(model->dimP_free, 1e2*absTol);
    //	cvodeFlag = CVodeSensSStolerances(cvodeMem, 1e2*relTol, vecAbsTol.data());
    //	//cvodeFlag = CVodeSetSensDQMethod(cvode mem, DQtype, DQrhomax);

    return 0;
}

SimInterface::SimInterface(const std::string &name, enumJacTypes type) :
			                        wsY(nullptr), wsP(nullptr), wsT(nullptr), wsSwitch(nullptr), isImplicitSwitched(nullptr),
			                        nvecY(nullptr), nvecArrSens(nullptr), model(nullptr), cvodeMem(nullptr), userData(nullptr),
			                        JMatrix(nullptr), LSolver(nullptr), jacCols(nullptr), jacRows(nullptr) {
    model = getModelByName(name);
    if(model) {
        jacType = type;

        InitializeMemory();
        errorFlag = InitializeCVODE();
        if(errorFlag != 0) {
            printfPtr("[ERROR] CVODE: %d", errorFlag);
        }
    }
    else
        errorFlag = -1;
}

SimInterface::~SimInterface() {
    delete[] wsY;
    delete[] wsP;
    delete[] wsT;
    delete[] wsSwitch;

    delete[] isImplicitSwitched;

    delete (UserData*) userData;
    if(nvecY!=nullptr) myN_VDestroy(nvecY);
    if(nvecArrSens!=nullptr) myN_VDestroyArray( model->dimP_free, nvecArrSens );
    CVodeFree(&cvodeMem);
    if(LSolver!=nullptr) SUNLinSolFree(LSolver);
    if(JMatrix!=nullptr) SUNMatDestroy(JMatrix);
    delete[] jacCols;
    delete[] jacRows;
}

bool SimInterface::Simulate(double *obs, double *traj, double *sens) {

    if(errorFlag != 0) {
        printfPtr("[ERROR] Model is not initialized correctly, aborting simulation.\n");
        return false;
    }

    unsigned int maxRestarts = 5;
    unsigned int noRestarts = 0;

    int cvodeFlag;

    progress = 0.0;

    memcpy(&wsP[model->dimP_free],&model->P_init[model->dimP_free],(model->dimP - model->dimP_free)*sizeof(double));
    memcpy( wsSwitch             , model->S_init                  , model->dimS                    *sizeof(unsigned int));

    model->ODEInitialParameters(wsP);
    model->ODEInitialValues(wsY, wsP, wsSwitch);

//    model->ODEOptions(relTol, absTol, hInit, hMin, hMax, maxSteps, useJac, wsP);
//    CVodeMem test = (CVodeMem) cvodeMem;
//    double rel = test->cv_reltol;
//    double abs = test->cv_Sabstol;
//    this->printfPtr("Rel: %e, Abs: %e\n", rel, abs);

    unsigned int posOutput = 0;
    std::set<double> explicitSwitches = model->explicitSwitches(wsP);
    std::set<double>::iterator iterSwitch = explicitSwitches.begin();
    double t_s = wsT[posOutput];
    double t_f;

    model->ODEExplicitSwitch(t_s, wsY, wsP, wsSwitch);

    if( iterSwitch!=explicitSwitches.end() && t_s==*iterSwitch )
        ++iterSwitch;

    if(traj!=nullptr) memcpy(traj,wsY,model->dimY*sizeof(double));
    if(obs!=nullptr) model->ODEObservers(t_s, wsY, wsP, wsSwitch, obs);

    cvodeFlag = CVodeReInit(cvodeMem, t_s, nvecY);
    if( cvodeFlag!=CV_SUCCESS ) return false;

    //	if(sens==nullptr)
    //		CVodeSensToggleOff(cvodeMem);
    //	else
    //		CVodeSensReInit(cvodeMem, CV_STAGGERED, nvecArrSens);

    while(posOutput<dimT-1) {
        if(iterSwitch==explicitSwitches.end()||*iterSwitch>wsT[posOutput+1]) {
            // next stop is output
            ++posOutput;
            t_f = wsT[posOutput];
        }
        else if(*iterSwitch==wsT[posOutput+1]) {
            // next stop is explicit switch and output
            ++posOutput;
            ++iterSwitch;
            t_f = wsT[posOutput];
        }
        else { //*iterSwitch<wsT[posOutput+1]
            // next stop is switch only
            t_f = *iterSwitch;
            ++iterSwitch;
        }

        cvodeFlag = CVode(cvodeMem, t_f, nvecY, &t_s, CV_NORMAL);

        // error codes drom the CVode manual
        // >= 0:
        // CV_SUCCESS: CVode succeeded and no roots were found.
        // CV_TSTOP_RETURN:  CVode succeeded by reaching the stopping point specified through the optional input function CVodeSetStopTime (see 4.5.6.1).
        // CV_ROOT_RETURN: CVode succeeded and found one or more roots. In this case, tret is the location of the root. If nrtfn > 1,
        //                 call CVodeGetRootInfo to see which gi were found to have a root.
        if(cvodeFlag<0) {
            switch(cvodeFlag) {
                case CV_MEM_NULL:
                    printfPtr("[ERROR] CVODE returned CV_MEM_NULL: The cvode mem argument was NULL. t_s = %e\n", t_s);
                    return false;
                case CV_NO_MALLOC:
                    printfPtr("[ERROR] CVODE returned CV_NO_MALLOC: The cvodes memory was not allocated by a call to CVodeInit. t_s = %e\n", t_s);
                    return false;
                case CV_ILL_INPUT:
                    printfPtr("[ERROR] CVODE returned CV_ILL_INPUT: One of the inputs to CVode was illegal. t_s = %e\n", t_s);
                    return false;
                case CV_TOO_CLOSE:
                    if(hMin == 0.0) {
                        printfPtr("[ERROR] CVODE returned CV_TOO_CLOSE: The initial time t0 and the final time tout are too close to each other. t_s = %e\n", t_s);
                        return false;
                    }
                    else {
                        printfPtr("[WARNING] CVODE returned CV_TOO_CLOSE: The initial time t0 and the final time tout are too close to each other. t_s = %e\n", t_s);
                        CVodeSetMinStep(cvodeMem, 0.0); // try with 0.0
                    }
                    break;
                case CV_TOO_MUCH_WORK:
                    printfPtr("[ERROR] CVODE returned CV_TOO_MUCH_WORK: The solver took mxstep internal steps but still could not reach tout. t_s = %e\n", t_s);
                    return false;
                case CV_TOO_MUCH_ACC:
                    printfPtr("[WARNING] CVODE returned CV_TOO_MUCH_ACC: The solver could not satisfy the accuracy demanded by the user for some internal step. t_s = %e\n", t_s);
                    //CVodeSStolerances(cvodeMem, 10*relTol, 10*absTol); // TODO: really change or just restart?
                    break;
                case CV_ERR_FAILURE:
                    printfPtr("[WARNING] CVODE returned CV_ERR_FAILURE: Either error test failures occurred too many times (MXNEF = 7) during one internal time step, or with |h| = hmin. t_s = %e\n", t_s);
                    break;
                case CV_CONV_FAILURE:
                    printfPtr("[WARNING] CVODE returned CV_CONV_FAILURE: Either convergence test failures occurred too many times (MXNCF = 10) during one internal time step, or with |h| = hmin. t_s = %e\n", t_s);
                    break;
                case CV_LINIT_FAIL:
                    printfPtr("[WARNING] CVODE returned CV_LINIT_FAIL: The linear solver's initialization function failed. t_s = %e\n", t_s);
                    break;
                case CV_LSETUP_FAIL:
                    printfPtr("[WARNING] CVODE returned CV_LSETUP_FAIL: The linear solver's setup function failed in an unrecoverable manner. t_s = %e\n", t_s);
                    break;
                case CV_LSOLVE_FAIL:
                    printfPtr("[WARNING] CVODE returned CV_LSOLVE_FAIL: The linear solver's solve function failed in an unrecoverable manner. t_s = %e\n", t_s);
                    break;
                case CV_RHSFUNC_FAIL:
                    printfPtr("[WARNING] CVODE returned CV_RHSFUNC_FAIL: The right-hand side function failed in an unrecoverable manner. t_s = %e\n", t_s);
                    break;
                case CV_FIRST_RHSFUNC_ERR:
                    printfPtr("[WARNING] CVODE returned CV_FIRST_RHSFUNC_ERR: The right-hand side function had a recoverable error at the first call. t_s = %e\n", t_s);
                    break;
                case CV_REPTD_RHSFUNC_ERR:
                    printfPtr("[WARNING] CVODE returned CV_REPTD_RHSFUNC_ERR: Convergence test failures occurred too many times due to repeated recoverable errors in the RHS function. t_s = %e\n", t_s);
                    break;
                case CV_UNREC_RHSFUNC_ERR:
                    printfPtr("[WARNING] CVODE returned CV_UNREC_RHSFUNC_ERR: The right-hand function had a recoverable error, but no recovery was possible. t_s = %e\n", t_s);
                    break;
                case CV_RTFUNC_FAIL:
                    printfPtr("[WARNING] CVODE returned CV_RTFUNC_FAIL: The rootfinding function failed. t_s = %e\n", t_s);
                    break;
                default:
                    printfPtr("[ERROR] CVODE: Unknown Error in CVODE. t_s = %e\n", t_s);
                    return false;
            }
            if(noRestarts<maxRestarts) {
                ++noRestarts;
                CVodeReInit(cvodeMem, t_s, nvecY);
                cvodeFlag = CVode(cvodeMem, t_f, nvecY, &t_s, CV_NORMAL);
                if(cvodeFlag<0) {
                    printfPtr("[ERROR] CVODE: Restart failed.\n");
                    return false; // failed once more
                }
            }
            else {
                printfPtr("[ERROR] CVODE: Maximum number of restarts exceeded.\n");
                return false;
            }
        }

        if( t_f == wsT[posOutput] ) {
            if(traj!=nullptr) memcpy(&traj[posOutput*model->dimY],wsY,model->dimY*sizeof(double));
            if(obs!=nullptr) model->ODEObservers(t_s, wsY, wsP, wsSwitch, &obs[posOutput*model->dimObservers]);
        }

        progress = (wsT[posOutput] - wsT[0]) / (wsT[dimT-1] - wsT[0]);

        if(posOutput<dimT-1) {
            bool reinit = model->ODEExplicitSwitch(t_s, wsY, wsP, wsSwitch);

            if(reinit) {
                //				std::stringstream s;
                //				s << "Switch at " << t_s << std::endl;
                //				outputString(s.str());

                CVodeReInit(cvodeMem, t_s, nvecY);
            }
        }

        //		if(model->dimImplicitSwitches) {
        //			CVodeGetRootInfo(cvodeMem,isImplicitSwitched);
        //			for(unsigned int i=0; i<model->dimImplicitSwitches; i++)
        //				if(isImplicitSwitched[i]) {
        //					model->ODEExplicitSwitch(t_s, wsY, wsP, wsS);
        //					//memcpy(res,ws_y,dimY*sizeof(double));
        //					CVodeRootInit(cvodeMem, 0, g);
        //					CVodeReInit(cvodeMem, t_s, nvecY);
        //					break;
        //				}
        //		}
    }

    //	if(reduceOrder||noRestarts) {
    //		long int no_reductions;
    //		CVodeGetNumStabLimOrderReds(cvodeMem, &no_reductions);
    //		std::stringstream str;
    //		str << "No of reductions: " << no_reductions << ", no of restarts: " << noRestarts;
    //		outputString(str.str());
    //	}

    //	static long int maxSteps = 0;
    //	long int noSteps;
    //	CVodeGetNumSteps(cvodeMem, &noSteps);
    //	if(noSteps>maxSteps) maxSteps = noSteps;
    //	std::stringstream sstr;
    //	sstr << "No steps: " << noSteps << "  max: " << maxSteps;
    //	outputString(sstr.str());

    return true;
}

void SimInterface::setT(unsigned int newDimT, const double *t) {
    if( newDimT != dimT ) {
        delete[] wsT;
        wsT=new double[newDimT];
        dimT=newDimT;
    }
    memcpy(wsT,t,dimT*sizeof(double));
    //	UpdateTimePoints();
}

void SimInterface::setRelTol(double _relTol) { relTol = _relTol; CVodeSStolerances(cvodeMem, relTol, absTol); };
void SimInterface::setAbsTol(double _absTol) { absTol = _absTol; CVodeSStolerances(cvodeMem, relTol, absTol); };
void SimInterface::setHInit (double _hInit)  { hInit  = _hInit; CVodeSetInitStep(cvodeMem, hInit); };
void SimInterface::setHMin  (double _hMin)   { hMin   = _hMin;  CVodeSetMinStep(cvodeMem, hMin); };
void SimInterface::setHMax  (double _hMax)   { hMax   = _hMax;  CVodeSetMaxStep(cvodeMem, hMax); };
void SimInterface::setMaxSteps(long int _maxSteps) { maxSteps = _maxSteps; CVodeSetMaxNumSteps(cvodeMem, maxSteps); };

// load/free shared library
#ifdef _WIN32
#include "windows.h"
//const char * LIBNAME = "ModelLibrary.dll";
//HINSTANCE libHandle = nullptr;
const char* LIBEXTENSION = ".dll";
typedef HINSTANCE libHandleType;
#define LOADLIBRARY(x) LoadLibrary(x)
#define GETFUNCTION(x) GetProcAddress(libHandle, x)
#define FREELIBRARY(x) FreeLibrary(x)
#define GETLIBRARYERROR GetLastError()
#else
#include <dlfcn.h>
//const char* LIBNAME = "ModelLibrary.so";
//void* libHandle;
const char* LIBEXTENSION = ".so";
typedef void* libHandleType;
#define LOADLIBRARY(x) dlopen(x,RTLD_NOW)
#define GETFUNCTION(x) dlsym(libHandle, x)
#define FREELIBRARY(x) dlclose(x)
#define GETLIBRARYERROR dlerror()
#endif

std::map<std::string, libHandleType> vecLibHandles;
std::map<std::string, Model*> vecModels;

Model* SimInterface::getModelByName(const std::string &nameModel) {
    std::map<std::string, Model*>::iterator mIter = vecModels.find(nameModel);
    if(mIter==vecModels.end()) {
        std::string nameLib(nameModel + LIBEXTENSION);
        SimInterface::printfPtr("[INFO] Model name not found, trying to load shared library: %s\n", nameLib.c_str());
        if(loadModelLibrary(nameLib)) {
            mIter = vecModels.find(nameModel);
            if(mIter==vecModels.end()) {
                SimInterface::printfPtr("[ERROR] Loaded shared library, but model name was not found.\n");
				std::stringstream sstr;
				for(std::map<std::string, Model*>::iterator mIter = vecModels.begin(); mIter !=vecModels.end(); ++mIter)
					sstr << " " << mIter->first;
				SimInterface::printfPtr("[DEBUG] Available model names:%s\n", sstr.str().c_str());
                return nullptr;
            }
        }
        else {
            SimInterface::printfPtr("[ERROR] Failed to load shared library: %s\n", nameLib.c_str());
            return nullptr;
        }
    }
    return mIter->second;
}

bool SimInterface::loadModelLibrary(const std::string &nameLib) {
    std::map<std::string, libHandleType>::iterator libHandleIter = vecLibHandles.find(nameLib);

    if(libHandleIter!=vecLibHandles.end())
        return true; // already loaded

    libHandleType libHandle = LOADLIBRARY(nameLib.c_str());

    if (libHandle!=nullptr) {
        unsigned int (*getDimModels)() = (unsigned int(*)()) GETFUNCTION("getDimModels");
        Model* (*getModelPtr)(unsigned int) = (Model*(*)(unsigned int)) GETFUNCTION("getModelPtr");
        if(getDimModels==nullptr || getModelPtr==nullptr) {
            FREELIBRARY(libHandle);
            SimInterface::printfPtr("[ERROR] Failed to get function addresses.\n");
            return false;
        }

        unsigned int dimModels = getDimModels();
        SimInterface::printfPtr("[INFO] Loaded shared library %s with %u models.\n", nameLib.c_str(), dimModels);
        for(unsigned int i=0; i<dimModels;i++) {
            Model* m = getModelPtr(i);
            std::pair<std::map<std::string, Model*>::iterator, bool> infoPair = vecModels.insert(std::pair<std::string, Model*>(m->name, m));
            if(!infoPair.second)
                SimInterface::printfPtr("[WARNING] Duplicate model name: %s. Not loading second version.\n", m->name.c_str());
            else
                SimInterface::printfPtr("[INFO] Loaded model %s.\n", m->name.c_str());
        }
        vecLibHandles.insert(std::pair<std::string, libHandleType>(nameLib, libHandle));
        //SimInterface::getDimModels = (unsigned int(*)()) GETFUNCTION("getDimModels");
        //SimInterface::getModelPtr = (Model*(*)(unsigned int)) GETFUNCTION("getModelPtr");
//        if(SimInterface::getDimModels==nullptr || SimInterface::getModelPtr==nullptr) {
//            FREELIBRARY;
//            SimInterface::printfPtr("[ERROR] Failed to get function addresses.\n");
//            return false;
//        }
    }
    else {
		// only a single call to GETLIBRARYERROR as it returns the last error only once
		char* e = NULL;
#ifdef _WIN32
      ::FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM,
         NULL, GetLastError(),
         MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
         (LPTSTR)&e, 0, NULL);
#else
      e = GETLIBRARYERROR;
#endif
		if (e!=nullptr)
			SimInterface::printfPtr("[ERROR] Failed to load library (Error: %s).\n", e);
		else
			SimInterface::printfPtr("[ERROR] Failed to load library (unknown error).\n");
        return false;
    }
    return true;
}

void SimInterface::freeModelLibrary() {
    for(std::map<std::string, libHandleType>::iterator mIter=vecLibHandles.begin(); mIter!=vecLibHandles.end(); ++mIter)
        FREELIBRARY(mIter->second);
    vecLibHandles.clear();
    vecModels.clear();
    //libHandle = nullptr;
}
//////

#ifdef DER_CPPAD

// TODO: not thread safe yet

CppAD::ADFun<double> ffcn;

std::vector<double> yad_pos;
std::vector<double> yad_dir;
std::vector<double> yad_res;

int JacDenseCppAD(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {

    UserData *d = (UserData*) user_data;
    double *y_ptr = N_VGetArrayPointer(y);

    memcpy(yad_pos.data(), y_ptr, d->model->dimY*sizeof(double));
    yad_pos[d->model->dimY] = t;

    yad_res = ffcn.Forward(0 , yad_pos);

    for( int i = 0; i < d->model->dimY; i++)  {
        yad_dir[i] = 1.0;
        yad_res = ffcn.Forward( 1, yad_dir );

        memcpy(SUNDenseMatrix_Column(J, i), yad_res.data(), d->model->dimY*sizeof(double));
        yad_dir[i] = 0.0;
    }

    return 0;
}

void SimInterface::CppadRecord() {
    //	model->ODEInitialValues(wsY, wsP);

    // domain
    std::vector< CppAD::AD<double> > y(model->dimY+1);

    for(int i=0; i<model->dimY; i++)
        y[i] = wsY[i];
    y[model->dimY] = 0.0;

    // start recording
    CppAD::Independent(y);

    // range
    std::vector< CppAD::AD<double> > dy(model->dimY);
    model->ODERHSFunction(y[model->dimY], y.data(), wsP, dy.data());

    // stop recording
    ffcn.Dependent(y, dy);

    ffcn.optimize();

    yad_pos.resize(model->dimY+1,0.0);
    yad_dir.resize(model->dimY+1,0.0);
    yad_res.resize(model->dimY);
}

void SimInterface::CppADJac(const double t, const double *y_ptr, double *jac) {

    memcpy( yad_pos.data(), y_ptr, model->dimY*sizeof(double) );
    yad_pos[ model->dimY ] = t;

    yad_res = ffcn.Forward( 0 , yad_pos );

    for( int i = 0; i < model->dimY; i++)  {
        yad_dir[i] = 1.0;
        yad_res = ffcn.Forward( 1, yad_dir );
        memcpy( &jac[model->dimY*i], yad_res.data(), model->dimY*sizeof(double));
        yad_dir[i] = 0.0;
    }
}

#endif
