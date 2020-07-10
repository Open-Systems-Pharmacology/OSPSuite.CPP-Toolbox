#include "Integrator.h"
#include "Interface.h"

#include "mex.h"

/*
 *
 *
 *
 */


const unsigned int NAMELENGTH = 1024;
bool isInitialized = false;
unsigned int simIndex = -1;
unsigned int dimT, dimY, dimO, dimP_free;

struct resultFuture {
    mxArray* success;
    mxArray* obs;
    mxArray* traj;
    double *p;

    resultFuture() : success(nullptr), obs(nullptr), traj(nullptr), p(nullptr) {}
    ~resultFuture() { mxDestroyArray(success); mxDestroyArray(obs); mxDestroyArray(traj); delete[] p; }
    void get(mxArray* _success, mxArray* _obs, mxArray* _traj) { // clear without releasing the memory that is returned to matlab
        _success = success; success = nullptr;
        _obs     = obs;     obs     = nullptr;
        _traj    = traj;    traj    = nullptr;
        delete[] p;
    }
};
std::vector<resultFuture*> vecResults;

int nnz;
const int *rows, *cols;

extern "C" {

void mexCleanup() {
    // serial simulations
    freeSimulation(-1); // UINT_MAX to clear all simulations

    // parallel memory
    for(unsigned int i=0; i<vecResults.size(); i++)
        delete vecResults[i];
    vecResults.clear();
    freeThreads();

    // free library itself
    SimInterface::freeModelLibrary();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    // determine action by control char
    if(nrhs<2 || !mxIsChar(prhs[0]) )
        mexErrMsgIdAndTxt("MATLAB:SimMatlab", "SimMatlab requires at least two input arguments: one control char and one simulation index.");

    const mxChar *c = mxGetChars(prhs[0]);
    simIndex = roundMex(mxGetScalar(prhs[1]));

    // check initialization
    if(!isInitialized) {
        setOutputPtr(mexPrintf);
        mexAtExit(mexCleanup);
//        if(c[0]!='e') {// for export only set outputPtr
//            isInitialized = SimInterface::loadModelLibrary();
//            if(isInitialized) {
//                mexAtExit(mexCleanup);
//            }
//            else {
//                //         mexErrMsgIdAndTxt("MATLAB:SimMatlab", "Failed to load model shared library.\n");
//                mexErrMsgTxt("Failed to load model shared library.\n");
//                if(nlhs>0) {
//                    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
//                    mxGetPr(plhs[0])[0] = -1;
//                }
//                return;
//            }
//        }
    }

    //	if(simIndex!=oldSimIndex) {
    if( !getDim(simIndex, &dimT, &dimY, &dimO, &dimP_free) )
        if(c[0]!='i' && c[0]!='f' && c[0]!='e') { // init and free may be called with invalid indices
            //mexErrMsgIdAndTxt("MATLAB:SimMatlab", "Invalid simulation index.");
            mexWarnMsgTxt("Invalid simulation index.");
            return;
        }
    //	}

    //	oldSimIndex = simIndex;

    switch(c[0]) {
        case 'i': // init
            if(nrhs<3) {
                simIndex = initSimulation(simIndex, "Standard");
            }
            else {
                char simName[NAMELENGTH];
                mxGetString(prhs[2], simName, NAMELENGTH);
                simIndex = initSimulation(simIndex, simName);
            }

            if(simIndex==UINT_MAX)
                mexWarnMsgTxt("Invalid simulation name.");

            //if(nlhs>0) { // always return a status, because nlhs is 0 if only "ans" is used
                plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
                mxGetPr(plhs[0])[0] = (simIndex==UINT_MAX ? (double)-1 : (double)simIndex); // negative index for Matlab to indicate error, explicit cast to prevent conversion to unsigned int
            //}
            break;
        case 'f': // free
            freeSimulation(simIndex);
            break;
        case 'e': // export
        {
            if(nrhs<4||nrhs>7) {
                mexPrintf("%u\n", nrhs);
                mexErrMsgIdAndTxt("MATLAB:SimMatlab","Export requires four to seven input arguments.");
            }

            unsigned int N_obs = 0;
            double *observers = nullptr;
            if(nrhs>=5 && mxIsDouble(prhs[4]) && mxGetM(prhs[4])*mxGetN(prhs[4]) > 0) {
                N_obs = mxGetM(prhs[4])*mxGetN(prhs[4]);
                observers = mxGetPr(prhs[4]);
            }
            unsigned int N_par = 0;
            double *parameters = nullptr;
            if(nrhs>=6 && mxIsDouble(prhs[5]) && mxGetM(prhs[5])*mxGetN(prhs[5]) > 0) {
                N_par = mxGetM(prhs[5])*mxGetN(prhs[5]);
                parameters = mxGetPr(prhs[5]);
            }
            bool reduceVariables = true;
            if(nrhs>=7 && mxIsLogical(prhs[6])) {
                reduceVariables = mxGetLogicals(prhs[6])[0];
            }

            char xmlFile[NAMELENGTH];
            char outputDir[NAMELENGTH];
            char modelName[NAMELENGTH];
            mxGetString(prhs[1], xmlFile, NAMELENGTH);
            mxGetString(prhs[2], outputDir, NAMELENGTH);
            mxGetString(prhs[3], modelName, NAMELENGTH);

            exportSimulation(xmlFile, outputDir, modelName, N_obs, observers, N_par, parameters, reduceVariables);
            break;
        }
        case 'c': // check
        {
            if(nrhs<3||nrhs>5) {
                mexPrintf("%u\n", nrhs);
                mexErrMsgIdAndTxt("MATLAB:SimMatlab","Export requires three to five input arguments.");
            }

            unsigned int N_obs = 0;
            double *observers = nullptr;
            if(nrhs>=4 && mxIsDouble(prhs[3]) && mxGetM(prhs[3])*mxGetN(prhs[3]) > 0) {
                N_obs = mxGetM(prhs[3])*mxGetN(prhs[3]);
                observers = mxGetPr(prhs[3]);
            }
            unsigned int N_par = 0;
            double *parameters = nullptr;
            if(nrhs>=5 && mxIsDouble(prhs[4]) && mxGetM(prhs[4])*mxGetN(prhs[4]) > 0) {
                N_par = mxGetM(prhs[4])*mxGetN(prhs[4]);
                parameters = mxGetPr(prhs[4]);
            }

            char xmlFile[NAMELENGTH];
            mxGetString(prhs[2], xmlFile, NAMELENGTH);

            plhs[0] = mxCreateLogicalMatrix(1, 1);
            mxGetLogicals(plhs[0])[0] = checkSimulation(simIndex, xmlFile, N_obs, observers, N_par, parameters);
            break;
        }
        case 's': // simulate
        {
            if (nlhs < 1||nlhs > 3) {
                // output validation
                mexErrMsgIdAndTxt("MATLAB:SimMatlab", "Simulation requires one to three output arguments.");
            }

            plhs[0] = mxCreateLogicalMatrix(1, 1);

            if(nlhs<=1) {
                // nothing more to do
                mxGetLogicals(plhs[0])[0] = true;
                break;
            }

            double *observers = nullptr;
            double *trajectories = nullptr;

            if(nlhs>=2) {
                // include observers
                plhs[1] = mxCreateDoubleMatrix(dimO,dimT,mxREAL);
                observers = mxGetPr(plhs[1]);
            }
            if(nlhs>=3) {
                // include full trajectories
                plhs[2] = mxCreateDoubleMatrix(dimY,dimT,mxREAL);
                trajectories = mxGetPr(plhs[2]);
            }

            mxGetLogicals(plhs[0])[0] = simulate(simIndex,observers,trajectories,nullptr);
            break;
        }
        case 't': // set time
            if( nrhs < 3 ) {
                if( nlhs < 1 )
                    mexErrMsgIdAndTxt("MATLAB:SimMatlab","Get time requires output argument.");

                plhs[0] = mxCreateDoubleMatrix(dimT,1,mxREAL);
                getTimePoints(simIndex, mxGetPr(plhs[0]));
                break;
            }

            if( nrhs < 3 || !mxIsDouble(prhs[2]) )
                mexErrMsgIdAndTxt("MATLAB:SimMatlab:notDouble","Set time requires double vector with new time grid.");

            setTimePoints( simIndex, mxGetM(prhs[2])*mxGetN(prhs[2]), mxGetPr(prhs[2]) );
            break;
        case 'p': // set parameter
            if (nrhs == 2) { // get parameter
                if(nlhs<1)
                    mexErrMsgIdAndTxt("MATLAB:SimMatlab","Get all parameters requires one output argument.");
                plhs[0] = mxCreateDoubleMatrix(dimP_free,1,mxREAL);
                getAllParameters(simIndex, mxGetPr(plhs[0]));
            }
            else if (nrhs == 3) { // set parameter (fast)
                if( !mxIsDouble(prhs[2]) )
                    mexErrMsgIdAndTxt("MATLAB:SimMatlab:notDouble","Set parameter requires double vector with new parameter.");

                if( mxGetM(prhs[2])*mxGetN(prhs[2]) != dimP_free )
                    mexErrMsgIdAndTxt("MATLAB:SimMatlab:notVector","Set parameter input vector has wrong dimension.");

                setAllParameters(simIndex, mxGetPr(prhs[2]));
            }
            else if(nrhs == 4) { // set parameter (slow)
                if( !(mxIsDouble(prhs[2]) && mxIsDouble(prhs[3])) )
                    mexErrMsgIdAndTxt("MATLAB:SimMatlab:notDouble","Set parameter requires double vectors with indices and new parameters.");

                unsigned int dimArg2 = mxGetM(prhs[2])*mxGetN(prhs[2]);
                unsigned int dimArg3 = mxGetM(prhs[3])*mxGetN(prhs[3]);

                if( dimArg2 != dimArg3 )
                    mexErrMsgIdAndTxt("MATLAB:SimMatlab:notVector","Set parameter input vector has wrong dimension.");

                if(mxGetN(prhs[0])<2)
                    mexErrMsgIdAndTxt("MATLAB:SimMatlab","Set parameter: specification of value type missing.");

                if(c[1]=='a') // absolute
                    setParameters(simIndex, dimArg2, mxGetPr(prhs[2]), mxGetPr(prhs[3]), false);
                else if(c[1]=='r') // relative
                    setParameters(simIndex, dimArg2, mxGetPr(prhs[2]), mxGetPr(prhs[3]), true);
                else
                    mexErrMsgIdAndTxt("MATLAB:SimMatlab","Set parameter: unknown specification of index type.");

            }
            else {
                mexErrMsgIdAndTxt("MATLAB:SimMatlab","Setting p requires one or two additional input vectors.");
            }
            break;
        case 'g': // get parameter or observer index by id
        {
            if( nlhs != 1 || nrhs != 3 ) mexErrMsgIdAndTxt("MATLAB:SimMatlab", "Get parameter and observer index by id require three inputs and one output.");
            unsigned int dim = mxGetM(prhs[2])*mxGetN(prhs[2]);
            plhs[0] = mxCreateDoubleMatrix(dim,1,mxREAL);

            if(mxGetN(prhs[0])<2)
                mexErrMsgIdAndTxt("MATLAB:SimMatlab","Get parameter or observer index: specification of type missing.");

            if(c[1]=='p') // parameter
                getParameterIndex(simIndex, dim, mxGetPr(prhs[2]), mxGetPr(plhs[0]));
            else if(c[1]=='o') // observer
                getObserverIndex(simIndex, dim, mxGetPr(prhs[2]), mxGetPr(plhs[0]));
            else if(c[1]=='s') // states
                getStateIndex(simIndex, dim, mxGetPr(prhs[2]), mxGetPr(plhs[0]));
            else
                mexErrMsgIdAndTxt("MATLAB:SimMatlab","Get parameter or observer index: unknown specification of type.");
        }
        break;
        case 'o': // set options
            if(nrhs<4)
                mexErrMsgIdAndTxt("MATLAB:SimMatlab","Too few input arguments.");

            {
                double hInit    = nrhs < 5 ? -1 : mxGetScalar(prhs[4]);
                double hMin     = nrhs < 6 ? -1 : mxGetScalar(prhs[5]);
                double hMax     = nrhs < 7 ? -1 : mxGetScalar(prhs[6]);
                double maxSteps = nrhs < 8 ? -1 : mxGetScalar(prhs[7]);

                setOptions(simIndex, mxGetScalar(prhs[2]), mxGetScalar(prhs[3]), hInit, hMin, hMax, maxSteps);
            }

            break;
            //	case 'f': // get free params
            //		plhs[0] = mxCreateDoubleMatrix(dimP_free,1,mxREAL);
            //		memcpy( mxGetPr(plhs[0]), sim.wsP, dimP_free*sizeof(double) );
            //		break;
            //	case 'm': // get p_map for free params
            //		plhs[0] = mxCreateDoubleMatrix(dimP_free,1,mxREAL);
            //		ptr = mxGetPr(plhs[0]);
            //		for(unsigned int i=0; i<dimP_free; i++)
            //			ptr[i] = sim.getModel()->P_map[i];
            //		break;
        case 'r': // get rhs
            if(nlhs<1||nrhs<2)
                mexErrMsgIdAndTxt("MATLAB:SimMatlab","Too few input or output arguments.");

            plhs[0] = mxCreateDoubleMatrix(dimY,1,mxREAL);

            // supply nullptr to evaluate at initial value
            getRHS(simIndex, ( nrhs<3 ? nullptr : mxGetPr(prhs[2]) ), mxGetPr(plhs[0]));

            break;
        case 'd': // get derivative
            if(mxGetN(prhs[0])<2)
                mexErrMsgIdAndTxt("MATLAB:SimMatlab","Specification of Jacobian type missing.");
            if(nlhs<1||nrhs<2)
                mexErrMsgIdAndTxt("MATLAB:SimMatlab","Too few input or output arguments.");
            if(c[1]=='d') { // dense
                plhs[0] = mxCreateDoubleMatrix(dimY,dimY,mxREAL);
                getJacDense(simIndex, ( nrhs<3 ? nullptr : mxGetPr(prhs[2]) ), mxGetPr(plhs[0]));
            }
            else if(c[1]=='s') { // get sparse derivative
                getJacInfo(simIndex, &nnz, &rows, &cols);
                plhs[0] = mxCreateSparse(dimY,dimY,nnz,mxREAL);
                getJacSparse(simIndex, ( nrhs<3 ? nullptr : mxGetPr(prhs[2]) ), mxGetPr(plhs[0]));
                for(int i=0; i<nnz; i++)
                    mxGetIr(plhs[0])[i] = rows[i];
                for(unsigned int i=0; i<dimY+1; i++)
                    mxGetJc(plhs[0])[i] = cols[i];
                //		memcpy(mxGetIr(plhs[0]), sim.getModel()->jacRows, (sim.getModel()->dimNNZ)*sizeof(int));
                //		memcpy(mxGetJc(plhs[0]), sim.getModel()->jacCols, (sim.getModel()->dimY+1)*sizeof(int));
            }
            else if(c[1]=='f') { // finite difference dense
                plhs[0] = mxCreateDoubleMatrix(dimY,dimY,mxREAL);
                getJacDenseFD(simIndex, ( nrhs<3 ? nullptr : mxGetPr(prhs[2]) ), mxGetPr(plhs[0]));
            }
            else
                mexErrMsgIdAndTxt("MATLAB:SimMatlab","Specification of Jacobian type unknown.");
            break;
        case 'a': { // threaded simulations

            if(mxGetN(prhs[0])<2)
                mexErrMsgIdAndTxt("MATLAB:SimMatlab","Specification of async operation missing.");

            unsigned int resIndex = roundMex(mxGetScalar(prhs[2]));

            switch(c[1]) {
                case 'o':
                case 't':
                case 'b':
                {
                    resultFuture* res = new resultFuture();
                    res->success = mxCreateLogicalMatrix(1, 1);
                    if(c[1]=='o' || c[1]=='b')  {
                        // include observers
                        res->obs = mxCreateDoubleMatrix(dimO,dimT,mxREAL);
                    }
                    if(c[1]=='t' || c[1]=='b') {
                        // include full trajectories
                        res->traj = mxCreateDoubleMatrix(dimY,dimT,mxREAL);
                    }
                    if(vecResults[resIndex]!=nullptr) {
                        mexWarnMsgTxt("Overwriting result future, please make sure that the former task was not in the queue anymore.\n");
                        delete vecResults[resIndex];
                    }

                    vecResults[resIndex] = res;
                    parallelSim(simIndex, resIndex, mxGetLogicals(res->success), mxGetPr(res->obs), mxGetPr(res->traj));
                }
                break;
                case 'i':
                {
                    unsigned int resIndex  = roundMex(mxGetScalar(prhs[2]));
                    unsigned int noThreads = roundMex(mxGetScalar(prhs[3]));
                    for(unsigned int i=0; i<vecResults.size(); i++)
                        delete vecResults[i];
                    vecResults = std::vector<resultFuture*>(resIndex,nullptr);
                    initThreads(noThreads);
                }
                break;
                case 'f':
                    for(unsigned int i=0; i<vecResults.size(); i++)
                        delete vecResults[i];
                    vecResults.clear();
                    freeThreads();
                    break;
            }
        }
        break;
        case 'j':
        {
            unsigned int resIndex = roundMex(mxGetScalar(prhs[2]));

            plhs[0] = vecResults[resIndex]->success;
            plhs[1] = vecResults[resIndex]->obs;
            plhs[2] = vecResults[resIndex]->traj;
            vecResults.clear();
        }
        break;
        default:
            mexErrMsgIdAndTxt("MATLAB:SimMatlab","Unknown control char.");
    }
}

}
