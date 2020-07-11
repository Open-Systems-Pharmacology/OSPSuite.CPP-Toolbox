/*
 * Integrator.h
 *
 *  Created on: 12.03.2016
 *      Author: GGKMT
 */

#ifndef SRC_INTEGRATOR_H_
#define SRC_INTEGRATOR_H_

// glibc workaround for old glibc version on cluster
#ifndef _WIN32
__asm__(".symver memcpy,memcpy@GLIBC_2.2.5");
#endif

//#include <cstdint>
#include <cassert>
#include <cstring>
#include <set>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <ctime>
#include <climits>

#include "Model.h"
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>

// flags to disable certain parts if libraries/sources are not available
//#define PARALLEL
//#define DER_CPPAD
#define DER_DENSE
#define DER_KLU
//#define DER_SUPERLU

#if __cplusplus < 201103L
#define nullptr NULL
#endif

//struct N_Vector;
typedef struct _generic_N_Vector *N_Vector;

class SimInterface {
public:
    double* wsY;
    double* wsP;
    double* wsT;
    double* wsSens;
    unsigned int *wsSwitch;

    enum enumJacTypes { DENSE_FD, DENSE_CPP, DENSE_CPPAD, SPARSE_KLU, SPARSE_SUPERLU };

private:
    unsigned int dimT;
    unsigned int *isImplicitSwitched;

    N_Vector nvecY;
    N_Vector *nvecArrSens;
    Model *model;

    void* cvodeMem;
    void* userData;

    SUNMatrix JMatrix;
    SUNLinearSolver LSolver;

    sunindextype *jacCols;
    sunindextype *jacRows;

    // marked options are accessed and modified by linked model file
    double absTol; // *
    double relTol; // *
    double hInit; // *
    double hMin; // *
    double hMax; // *
    long int maxSteps; // *
    int maxOrder;

    bool useJac; // * // -> removed
    bool reduceOrder;
    bool isScaled;

    double progress;
    int errorFlag;

    enumJacTypes jacType;

    void InitializeMemory();
    int InitializeCVODE();
    void swap(SimInterface &s);

public:
//    SimInterface(unsigned int index = 0, enumJacTypes type = SPARSE_KLU );
    SimInterface(const std::string &name, enumJacTypes type = SPARSE_KLU);
    SimInterface(const SimInterface& );
    SimInterface& operator=(SimInterface s);
    //	SimInterface(SimInterface&&) = default;
    ~SimInterface();

    bool Simulate(double *obs, double *traj, double *sens);
    void setT(unsigned int newDimT, const double *t);
    void setRelTol(double _relTol);
    void setAbsTol(double _absTol);
    void setHInit (double _hInit);
    void setHMin  (double _hMin);
    void setHMax  (double _hMax);
    void setMaxSteps(long int _maxSteps);

    inline double getProgress() const { return progress; }
    inline unsigned int getDimT() const { return dimT; }
    inline unsigned int getDimY() const { return model->dimY; }
    inline unsigned int getDimP() const { return model->dimP; }
    inline unsigned int getDimP_free() const { return model->dimP_free; }
    inline unsigned int getDimO() const { return model->dimObservers; }
    inline std::string getName() const { return model->name; }
    inline Model* getModel() const { return model; }
    inline int getErrorFlag() const { return errorFlag; }

    // static
    static int(*printfPtr)(const char *c,...);
    static inline void setPrintfPtr( int(*ptr)(const char *,...) ) { printfPtr = ptr; }

//    static unsigned int (*getDimModels)();
//    static Model* (*getModelPtr)(unsigned int);
    static Model* getModelByName(const std::string &name);
    static bool loadModelLibrary(const std::string &name);
    static void freeModelLibrary();

    inline unsigned int getParameterIndex(const unsigned int id) const {
        for(unsigned int pos=0; pos<model->dimP_free; pos++)
            if(model->P_map[pos]==id)
                return pos;
        return -1;
    }

    inline unsigned int getObserverIndex(const unsigned int id) const {
        for(unsigned int pos=0; pos<model->dimObservers; pos++)
            if(model->O_map[pos]==id)
                return pos;
        return -1;
    }

    inline unsigned int getStateIndex(const unsigned int id) const {
        for(unsigned int pos=0; pos<model->dimY; pos++)
            if(model->Y_map[pos]==id)
                return pos;
        return -1;
    }

#ifdef DER_CPPAD
void CppadRecord();
void CppADJac(const double t, const double *y_ptr, double *jac); // for testing purposes
#endif
};

// helper function to round double indices from matlab to nearest int
//inline int roundMex(const double dblIndex) {
//	return static_cast<int>( dblIndex>=0 ? dblIndex+0.5 : dblIndex-0.5 );
//}

// fast version for indices (using > 0)
inline unsigned int roundMex(const double dblIndex) {
    return static_cast<unsigned int>(dblIndex+0.5);
}

#endif /* SRC_INTEGRATOR_H_ */
