/*
 * Model.cpp
 *
 */

#include "Model.h"

#include <vector>
#include <iostream>

namespace {
std::vector<Model*>& getVecSimModels() { static std::vector<Model*> vecSimModels; return vecSimModels; };
}

// should only be called by ctor to register singleton static object
unsigned int registerModel(Model *m) {
    getVecSimModels().push_back(m);
    return getVecSimModels().size()-1;
}

unsigned int getDimModels() {
    return getVecSimModels().size();
}

Model* getModelPtr(unsigned int index) {
    return getVecSimModels()[index];
}

Model::Model(const std::string &_name, unsigned int _dimY, unsigned int _dimP,
        unsigned int _dimP_free, unsigned int _dimS, unsigned int _dimOutputTimePoints,
        unsigned int _dimImplicitSwitches, unsigned int _dimObservers, unsigned int _dimNNZ,
        unsigned int _dimNNZ_fp, unsigned int _dimNNZ_ox, unsigned int _dimNNZ_op,
        const unsigned int* _P_map, const unsigned int* _O_map, const unsigned int* _Y_map,
        const double *_Y_sca, const double* _P_init, const unsigned int* _S_init,
        const double *_outputTimePoints, const int *_jacRows, const int *_jacCols,
        const int *_jacRows_fp, const int *_jacCols_fp, const int *_jacRows_ox, const int *_jacCols_ox,
        const int *_jacRows_op, const int *_jacCols_op, uint64_t _hash,
        std::set<double> (*_explicitSwitches)(const double*)) :
                    name(_name), dimY(_dimY), dimP(_dimP), dimP_free(_dimP_free), dimS(_dimS), dimOutputTimePoints(_dimOutputTimePoints),
                    dimImplicitSwitches(_dimImplicitSwitches), dimObservers(_dimObservers),
                    dimNNZ(_dimNNZ), dimNNZ_fp(_dimNNZ_fp), dimNNZ_ox(_dimNNZ_ox), dimNNZ_op(_dimNNZ_op),
                    P_map(_P_map), O_map(_O_map), Y_map(_Y_map), Y_sca(_Y_sca), P_init(_P_init), S_init(_S_init), outputTimePoints(_outputTimePoints),
                    jacRows(_jacRows), jacCols(_jacCols), jacRows_fp(_jacRows_fp), jacCols_fp(_jacCols_fp),
                    jacRows_ox(_jacRows_ox), jacCols_ox(_jacCols_ox), jacRows_op(_jacRows_op), jacCols_op(_jacCols_op),
                    index(registerModel(this)), hash(_hash) { explicitSwitches = _explicitSwitches; }
