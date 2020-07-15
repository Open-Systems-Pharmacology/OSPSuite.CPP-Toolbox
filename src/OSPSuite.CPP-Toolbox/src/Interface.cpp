#include "Integrator.h"
#include "ThreadPool.h"

static std::vector<SimInterface*> vecSimulations;
// defined below: static ThreadPool* pool = nullptr;

// attention: access is safeguarded only through getDim by default
inline bool checkIndex(unsigned int index) {
    return index<vecSimulations.size() && vecSimulations[index]!=nullptr;
}

void setOutputPtr(int(*outputPtr)(const char*, ...)) {
    SimInterface::setPrintfPtr(outputPtr);
}

unsigned int initSimulation(const unsigned int index, const char * name ) {
#ifdef DER_KLU
	auto jacobianType = SimInterface::SPARSE_KLU;
#else
	auto jacobianType = SimInterface::DENSE_CPP;
#endif

    SimInterface* sim = new SimInterface(std::string(name), jacobianType);//DENSE_FD, DENSE_CPP, SPARSE_KLU

    if(sim->getErrorFlag()<0) {
//        if(SimInterface::getDimModels!=nullptr) {
//            std::stringstream sstr;
//            for (unsigned int i=0; i<SimInterface::getDimModels(); i++)
//                sstr << SimInterface::getModelPtr(i)->name << "|";
//            sstr << std::endl;
//
//            SimInterface::printfPtr("[Error] Model name \"%s\" not found. Models available: %s\n", name, sstr.str().c_str());
//        }
//        else {
//            SimInterface::printfPtr("[Error] Model library not loaded properly.\n");
//        }
        delete sim;
        return -1;
    }

    if(index==UINT_MAX) {
        // try to find free index
        unsigned int freeIndex;
        for(freeIndex=0; freeIndex<vecSimulations.size(); freeIndex++)
            if(vecSimulations[freeIndex]==nullptr) {
                vecSimulations[freeIndex] = sim;
                return freeIndex;
            }
        // no free index -> append
        vecSimulations.push_back(sim);
        freeIndex = vecSimulations.size()-1;
        return freeIndex;
    }

    // increase size if necessary
    if( index>=vecSimulations.size() ) {
        vecSimulations.resize(index+1, nullptr);
    }

    if( vecSimulations[index]!=nullptr ) {
        SimInterface::printfPtr("[Warning] Replacing simulation.\n");
        delete vecSimulations[index];
    }
    vecSimulations[index] = sim;
    return index;
}

void freeSimulation(const unsigned int index) {
    if(index==UINT_MAX) { // clear all
        for(size_t i=0; i<vecSimulations.size();i++)
            delete vecSimulations[i];
        vecSimulations.clear();
    }
    else if( index<vecSimulations.size() ){
        delete vecSimulations[index];  vecSimulations[index] = nullptr;
    }
}

#ifdef USE_EXPERIMENTAL_FEATURES
#include "ImportXML.h"
void exportSimulation(const char *xmlFile, const char *outputDir, const char *modelName, unsigned int N_obs, const double *observer, unsigned int N_par, const double *parameter, bool reduceVariables) {
    std::set<unsigned int> setObserver;
    for(unsigned int i=0; i<N_obs; i++)
        setObserver.emplace(roundMex(observer[i]));
    std::set<unsigned int> setParameter;
    for(unsigned int i=0; i<N_par; i++)
        setParameter.emplace(roundMex(parameter[i]));

    convertXMLtoCpp(std::string(xmlFile), std::string(outputDir), std::string(modelName), setObserver, setParameter, reduceVariables);
}

bool checkSimulation(const unsigned int index, const char *xmlFile, unsigned int N_obs, const double *observer, unsigned int N_par, const double *parameter) {
    SimInterface * s = vecSimulations[index];

    std::set<unsigned int> setObserver;
    for(unsigned int i=0; i<N_obs; i++)
        setObserver.emplace(roundMex(observer[i]));
    std::set<unsigned int> setParameter;
    for(unsigned int i=0; i<N_par; i++)
        setParameter.emplace(roundMex(parameter[i]));

    uint64_t hash = getModelHash(std::string(xmlFile), setObserver, setParameter);

    if(hash==s->getModel()->hash) {
        SimInterface::printfPtr("[INFO] Model %u hash matches.\n", index);
        return true;
    }
    else {
        SimInterface::printfPtr("[ERROR] Model %u hash mismatch: %llu vs %llu. Please recompile.\n", index, hash, s->getModel()->hash);
        return false;
    }
}
#else
void exportSimulation(const char *xmlFile, const char *outputDir, const char *modelName, unsigned int N_obs, const double *observer, unsigned int N_par, const double *parameter, bool reduceVariables) {
    SimInterface::printfPtr("[ERROR] Experimental feature export simulation currently disabled.\n");
}
bool checkSimulation(const unsigned int index, const char *xmlFile, unsigned int N_obs, const double *observer, unsigned int N_par, const double *parameter) {
    SimInterface::printfPtr("[ERROR] Experimental feature check simulation currently disabled.\n");
    return false;
}
#endif

bool getDim(const unsigned int index, unsigned int* dimT, unsigned int* dimY, unsigned int* dimO, unsigned int* dimP_free) {
    if(!checkIndex(index))
        return false;

    *dimT = vecSimulations[index]->getDimT();
    *dimY = vecSimulations[index]->getDimY();
    *dimO = vecSimulations[index]->getDimO();
    *dimP_free = vecSimulations[index]->getDimP_free();
    return true;
}

unsigned int getDimNNZ(const unsigned int index) {
    return vecSimulations[index]->getModel()->dimNNZ;
}

void getTimePoints(const unsigned int index, double* t) {
    memcpy( t, vecSimulations[index]->wsT, vecSimulations[index]->getDimT()*sizeof(double) );
}

void setTimePoints(const unsigned int index, const unsigned int no, const double* t) {
    vecSimulations[index]->setT(no, t);
}

void setAllParameters(const unsigned int index, const double* p) {
    memcpy(vecSimulations[index]->wsP, p, vecSimulations[index]->getDimP_free()*sizeof(double));
}

void getAllParameters(const unsigned int index, double* p) {
    memcpy(p, vecSimulations[index]->wsP, vecSimulations[index]->getDimP_free()*sizeof(double));
}

void setParameters(const unsigned int index, const unsigned int ids_size, const double* ids, const double* p, bool isRelativeToInitialValue ) {
    SimInterface * s = vecSimulations[index];
    //	if(byPathID) {
    //		if(vecPermutations[index].empty()) {
    //			vecPermutations[index].resize(ids_size);
    //			for(unsigned int i=0; i<ids_size; i++) {
    //				vecPermutations[index][i] = s->getParameter(ids[i]);
    //				if(vecPermutations[index][i]==-1) {
    //					s->outputString("(ERROR) setParameters: unknown id!");
    //					vecPermutations[index].clear();
    //					return;
    //				}
    //			}
    //		}
    //
    //		const unsigned int * m = vecPermutations[index].data();
    //		const unsigned int size = vecPermutations[index].size();
    //
    //		for(unsigned int i=0; i<size; i++)
    //			s->wsP[m[i]]=p[i];
    //	}
    //	else {
    if(isRelativeToInitialValue) {
        for(unsigned int i=0; i<ids_size; i++) {
            unsigned int indexP = roundMex(ids[i]);
            s->wsP[indexP]=p[i]*s->getModel()->P_init[indexP];
        }
    }
    else {
        for(unsigned int i=0; i<ids_size; i++)
            s->wsP[roundMex(ids[i])]=p[i];
    }
    //	}
}

void getParameterIndex(const unsigned int index, const unsigned int ids_size, const double *ids, double *indices) {
    SimInterface * s = vecSimulations[index];

    for(unsigned int i=0; i<ids_size; i++) {
        unsigned int returnedIndex = s->getParameterIndex(roundMex(ids[i]));
        if(returnedIndex==UINT_MAX) {
            SimInterface::printfPtr("[ERROR] getParameterIndex - unknown id: %u\n", roundMex(ids[i]));
            indices[i] = -1;
        }
        else
            indices[i] = returnedIndex;
    }
}

void getObserverIndex(const unsigned int index, const unsigned int ids_size, const double *ids, double *indices) {
    SimInterface * s = vecSimulations[index];

    for(unsigned int i=0; i<ids_size; i++) {
        unsigned int returnedIndex = s->getObserverIndex(ids[i]);
        if(returnedIndex==UINT_MAX) {
            SimInterface::printfPtr("[ERROR] getObserverIndex - unknown id: %u\n", roundMex(ids[i]));
            indices[i] = -1;
        }
        else
            indices[i] = returnedIndex;
    }
}

void getStateIndex(const unsigned int index, const unsigned int ids_size, const double *ids, double *indices) {
    SimInterface * s = vecSimulations[index];

    for(unsigned int i=0; i<ids_size; i++) {
        unsigned int returnedIndex = s->getStateIndex(ids[i]);
        if(returnedIndex==UINT_MAX) {
            SimInterface::printfPtr("[ERROR] getStateIndex - unknown id: %u\n", roundMex(ids[i]));
            indices[i] = -1;
        }
        else
            indices[i] = returnedIndex;
    }
}

bool simulate(const unsigned int index, double* obs, double* traj, double* sens) {
    return vecSimulations[index]->Simulate(obs,traj,sens);
}

void getRHS(const unsigned int index, const double *y, double *rhs) {
    SimInterface * s = vecSimulations[index];
    s->getModel()->ODEInitialParameters(s->wsP);

    if(y==nullptr)
        s->getModel()->ODEInitialValues(s->wsY,s->wsP,s->wsSwitch);
    else
        memcpy(s->wsY, y, s->getDimY()*sizeof(double));

    s->getModel()->ODERHSFunction(1.0, s->wsY, s->wsP, s->wsSwitch, rhs);
}

void getJacInfo(const unsigned int index, int* nnz,  const int** rows, const int** cols) {
    *nnz = vecSimulations[index]->getModel()->dimNNZ;
    *rows = vecSimulations[index]->getModel()->jacRows;
    *cols = vecSimulations[index]->getModel()->jacCols;
}

void getJacSparse(const unsigned int index, const double *y, double *jac) {
    SimInterface * s = vecSimulations[index];
    s->getModel()->ODEInitialParameters(s->wsP);

    if(y==nullptr)
        s->getModel()->ODEInitialValues(s->wsY,s->wsP,s->wsSwitch);
    else
        memcpy(s->wsY, y, s->getDimY()*sizeof(double));

    s->getModel()->ODEJacSparse(0.0, s->wsY, s->wsP, s->wsSwitch, jac);
}

void getJacDense(const unsigned int index, const double *y, double *jac) {
    SimInterface * s = vecSimulations[index];
    s->getModel()->ODEInitialParameters(s->wsP);

    if(y==nullptr)
        s->getModel()->ODEInitialValues(s->wsY,s->wsP,s->wsSwitch);
    else
        memcpy(s->wsY, y, s->getDimY()*sizeof(double));

    s->getModel()->ODEJacDense(0.0, s->wsY, s->wsP, s->wsSwitch, jac);
}

void getJacDenseFD(const unsigned int index, const double *y, double *jac) {
    SimInterface * s = vecSimulations[index];
    s->getModel()->ODEInitialParameters(s->wsP);

    if(y==nullptr)
        s->getModel()->ODEInitialValues(s->wsY,s->wsP,s->wsSwitch);
    else
        memcpy(s->wsY, y, s->getDimY()*sizeof(double));

    double h=1e-8;
    double* rhs = new double[s->getDimY()];
    double* rhs_h = new double[s->getDimY()];

    s->getModel()->ODERHSFunction(0.0, s->wsY, s->wsP, s->wsSwitch, rhs);

    for(unsigned int i=0; i<s->getDimY(); i++) {
        double temp = s->wsY[i];
        s->wsY[i] += h;
        s->getModel()->ODERHSFunction(0.0, s->wsY, s->wsP, s->wsSwitch, rhs_h);
        for(unsigned int j=0; j<s->getDimY(); j++) {
            jac[i*s->getDimY()+j] = (rhs_h[j]-rhs[j])/h;
        }
        s->wsY[i] = temp;
    }

    delete[] rhs;
    delete[] rhs_h;
}

// negative values indicate: do not change
void setOptions(const unsigned int index, const double relTol, const double absTol, const double hInit,
                const double hMin, const double hMax, const long int maxSteps) {//, const double maxTime
    if(relTol  >=0) vecSimulations[index]->setRelTol(relTol);
    if(absTol  >=0) vecSimulations[index]->setAbsTol(absTol);
    if(hInit   >=0) vecSimulations[index]->setHInit(hInit);
    if(hMin    >=0) vecSimulations[index]->setHMin(hMin);
    if(hMax    >=0) vecSimulations[index]->setHMax(hMax);
    if(maxSteps>=0) vecSimulations[index]->setMaxSteps(maxSteps);
}

#if __cplusplus > 201103L
static ThreadPool* pool = nullptr;
void parallelSim(const unsigned int simIndex, const unsigned int resIndex, bool *success, double* obs, double* traj) {
    pool->pushTask(simIndex, resIndex, success, obs, traj);
}

void initThreads(unsigned int noThreads) {
    pool = new ThreadPool(noThreads, vecSimulations);
}

void freeThreads() {
    delete pool; pool = nullptr;
}
#else
void parallelSim(const unsigned int simIndex, const unsigned int resIndex, bool *success, double* obs, double* traj) {
    SimInterface::printfPtr("[ERROR] cppToolbox was compiled without thread support. Do not use this function or recompile.\n");
}

void initThreads(unsigned int noThreads) {
    SimInterface::printfPtr("[ERROR] cppToolbox was compiled without thread support. Do not use this function or recompile.\n");
}

void freeThreads() {
    SimInterface::printfPtr("[ERROR] cppToolbox was compiled without thread support. Do not use this function or recompile.\n");
}
#endif
