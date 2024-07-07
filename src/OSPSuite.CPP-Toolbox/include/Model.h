/*
 * Model.h
 *
 */

#ifndef SRC_MODEL_H_
#define SRC_MODEL_H_

#include<string>
#include<set>

#ifdef _MSC_VER
#define EXPORTSHARED __declspec(dllexport)
#else
#define EXPORTSHARED
#endif

struct Model;

extern "C" {
    EXPORTSHARED unsigned int getDimModels();
    EXPORTSHARED Model* getModelPtr(unsigned int index);
}

struct Model {
    const std::string name;

    const unsigned int dimY;
    const unsigned int dimP;
    const unsigned int dimP_free;
    const unsigned int dimS;
    const unsigned int dimOutputTimePoints;
    const unsigned int dimImplicitSwitches;
    const unsigned int dimObservers;
    const unsigned int dimNNZ;
    const unsigned int dimNNZ_fp;
    const unsigned int dimNNZ_ox;
    const unsigned int dimNNZ_op;

    const unsigned int* P_map;
    const unsigned int* O_map;
    const unsigned int* Y_map;
    const double* Y_sca;
    const double* P_init;
    const unsigned int* S_init;
    const double* outputTimePoints;

    // sparse format for cvodes requires ints
    const int* jacRows;
    const int* jacCols;
    const int* jacRows_fp;
    const int* jacCols_fp;
    const int* jacRows_ox;
    const int* jacCols_ox;
    const int* jacRows_op;
    const int* jacCols_op;

    const unsigned int index;
    const size_t hash;

    std::set<double> (*explicitSwitches)(const double*);

    virtual void ODEOptions(double &absTol, double &relTol, double &hInit, double &hMin, double &hMax, long int &maxSteps, bool &useJac, const double * P) = 0;
    virtual void ODEInitialValues(double *y, double *P, unsigned int *S) = 0;
    virtual void ODEInitialParameters(double *P) = 0;
    virtual void ODERHSFunction(const double Time, const double *y, const double *P, const unsigned int *S, double *dy) = 0;
    virtual void ODEObservers(const double Time, const double *y, const double *P, const unsigned int *S, double *obs) = 0;
    virtual bool ODEExplicitSwitch(const double Time, double *y, double *P, unsigned int *S) = 0;
    virtual void ODEImplicitSwitch(const double Time, double *y, double *P, unsigned int *S, double *gout) = 0;
    virtual void ODEJacDense    (const double Time, const double *y, const double *P, const unsigned int *S, double *J) = 0;
    virtual void ODEJacDense_fp (const double Time, const double *y, const double *P, const unsigned int *S, double *J) = 0;
    virtual void ODEJacDense_ox (const double Time, const double *y, const double *P, const unsigned int *S, double *J) = 0;
    virtual void ODEJacDense_op (const double Time, const double *y, const double *P, const unsigned int *S, double *J) = 0;
    virtual void ODEJacSparse   (const double Time, const double *y, const double *P, const unsigned int *S, double *J) = 0;
    virtual void ODEJacSparse_fp(const double Time, const double *y, const double *P, const unsigned int *S, double *J) = 0;
    virtual void ODEJacSparse_ox(const double Time, const double *y, const double *P, const unsigned int *S, double *J) = 0;
    virtual void ODEJacSparse_op(const double Time, const double *y, const double *P, const unsigned int *S, double *J) = 0;

#ifdef DER_CPPAD
virtual void ODERHSFunction(CppAD::AD<double> &Time, CppAD::AD<double> *y, const double *P, CppAD::AD<double> *dy) = 0;
#endif

Model(const std::string &_name, unsigned int _dimY, unsigned int _dimP, unsigned int _dimP_free, unsigned int _dimS, unsigned int _dimOutputTimePoints,
        unsigned int _dimImplicitSwitches, unsigned int _dimObservers, unsigned int _dimNNZ, unsigned int _dimNNZ_fp, unsigned int _dimNNZ_ox, unsigned int _dimNNZ_op,
        const unsigned int* _P_map, const unsigned int* _O_map, const unsigned int* _Y_map,
        const double* _Y_sca, const double* _P_init, const unsigned int* _S_init,
        const double *_outputTimePoints, const int *_jacRows, const int *_jacCols,
        const int *_jacRows_fp, const int *_jacCols_fp, const int *_jacRows_ox, const int *_jacCols_ox,
        const int *_jacRows_op, const int *_jacCols_op, uint64_t _hash,
        std::set<double> (*_explicitSwitches)(const double*));
virtual ~Model() { }
};

#endif /* SRC_MODEL_H_ */
