/*
 * ExportCpp.cpp
 *
 */

#include "ImportXML.h"
#include "Expression.h"
#include "Integrator.h" // outputString

#include <iostream>
#include <vector>
#include <limits>
#include <sstream>
#include <fstream>
#include <numeric> // partial_sum

inline void formatStream(std::ostream &os) {
    os.precision(std::numeric_limits<double>::max_digits10);
    os << std::scientific;
}

// reuseable output for nan and inf
inline std::string toString(double value) {
    if(std::isnan(value))
        return "NAN";
    else if(std::isinf(value))
        if(value>0)
            return "INFINITY";
        else
            return "-INFINITY";
    else {
        std::stringstream sstr;
        formatStream(sstr);
        sstr << value;
        return sstr.str();
    }

}

void writeParamLocal(const Simulation &sim, std::ostream &os) {
    for(auto iter = sim.parameters.begin(); iter != sim.parameters.end(); ++iter) {
        if( sim.isLocal[iter->id] && sim.isDependent[iter->id] ) {//!iter->canBeVaried &&
            os << "    const double " << iter->name << " = " << iter->e->toString() << ";" << std::endl;
        }
    }
}

void writeODEInitialValues(Simulation &sim, std::ostream &os) {
    // ODE initial values
    os << "void CLASSNAME(MODELNAME)::ODEInitialValues(double *__restrict__ y, double *__restrict__ P, unsigned int *__restrict__ S) {" << std::endl;

    sim.resetDependencies();
    for (auto iter = sim.variables.begin(); iter!=sim.variables.end(); ++iter) {
//        if( iter->e ) // formula is optional here
            iter->e->args[0]->getDependencies(sim.isDependent);
    }
    writeParamLocal(sim, os);

    unsigned int indexY = 0;
    for(auto iter = sim.variables.begin(); iter!=sim.variables.end(); ++iter) {
        os << "    y[" << indexY++ << "] = ";
        if( iter->formulaId == UINT_MAX )
            os << toString(iter->value); // toString for correct nan/inf output
        else {
            os << iter->e->args[0]->toString();//iter->scalefactor << "*("<< ")"
        }
        os << "; // " << iter->path << std::endl;
    }
    os << std::endl << "}" << std::endl << std::endl;
}

void writeODEParameterAndObserver(Simulation &sim, std::ostream &os) {

    os << "void CLASSNAME(MODELNAME)::ODEObservers(const double Time, const double *__restrict__ y, const double *__restrict__ P, const unsigned int *__restrict__ S, double *__restrict__ obs) {" << std::endl;

    // collect dependencies on parameters
    sim.resetDependencies();
    for (auto iter = sim.observers.begin(); iter!=sim.observers.end(); ++iter)
        iter->e->getDependencies(sim.isDependent);
    writeParamLocal(sim, os);

    unsigned int i = 0;
    for (auto iter = sim.observers.begin(); iter!=sim.observers.end(); ++iter)
        os << "    obs[" << i++ << "] = " << iter->e->toString() << "; // " << iter->path << std::endl;
    os << "}" << std::endl << std::endl;

    os << "void CLASSNAME(MODELNAME)::ODEInitialParameters(double *__restrict__ P) {" << std::endl;
    sim.resetDependencies();
    for (auto iter = sim.parameters.begin(); iter!=sim.parameters.end(); ++iter) {
        if(iter->e && !sim.isLocal[iter->id] && !iter->canBeVaried) // might not have a formula attached
            iter->e->getDependencies(sim.isDependent);
    }
    writeParamLocal(sim, os);
    for (auto iter = sim.parameters.begin(); iter!=sim.parameters.end(); ++iter) {
        if(iter->canBeSwitched) // TODO: merge?
            os << "    " << iter->name << " = " << iter->e->toString() << "; // " << sim.vecNodes[iter->id]->path << std::endl;
        else if(iter->e && !sim.isLocal[iter->id] && !iter->canBeVaried) // might not have a formula attached
            os << "    " << iter->name << " = " << iter->e->toString() << "; // " << sim.vecNodes[iter->id]->path << std::endl;
    }
    os << "}" << std::endl << std::endl;

    // initialize on first use and anonymous namespace
    os << "namespace {" << std::endl;
    os << "unsigned int _dimY() { static const unsigned int dimY = " << sim.variables.size() << "; return dimY; }" << std::endl << std::endl;

    // write output times vector
    os << "const double* _outputTimePoints() { static const double outputTimePoints[] = { ";
    unsigned int dimOutput = 0;
    for(auto iter = sim.outputTimepoints.begin(); iter!=sim.outputTimepoints.end(); ++iter) {
        if(dimOutput!=0)
            os << ", ";
        if (++dimOutput % 100 == 0)
            os << std::endl << "                                                                               ";
        os << *iter;
    }
    os << " }; return outputTimePoints; }" << std::endl;
    os << "unsigned int _dimOutputTimePoints() { static const unsigned int dimOutputTimePoints = " << dimOutput << "; return dimOutputTimePoints; }" << std::endl << std::endl;

    os << "const unsigned int* _positiveStateIndices() { static const unsigned int positiveStateIndices[] = { ";
    unsigned int dimPositiveStateIndices = 0;
    for(auto iter = sim.variables.begin(); iter != sim.variables.end(); ++iter)	{
        if (iter->negativeValuesAllowed) {
            if (dimPositiveStateIndices!=0)
                os << ", ";
            if (dimPositiveStateIndices % 100 == 0)
                os << std::endl << "                                                                                 ";
            os << dimPositiveStateIndices; //os << iter->second->GetODEIndex();
            ++dimPositiveStateIndices;
        }
    }
    os << "}; return positiveStateIndices; }" << std::endl;
    os << "const unsigned int _dimPositiveStateIndices() { static const unsigned int dimPositiveStateIndices = " << dimPositiveStateIndices << "; return dimPositiveStateIndices; }" << std::endl << std::endl;

    sim.resetDependencies();
    for(auto iter = sim.events.begin(); iter != sim.events.end(); ++iter)
        if (iter->e->isExplicit())
            iter->e->getDependencies(sim.isDependent);

    os << "std::set<double> _explicitSwitches(const double *P) {" << std::endl;
    os << "    std::set<double> s;" << std::endl;
    //os << writeParamLocal(usedIDs, os);// TODO: include

    // add evaluated formulas
    for(auto iter = sim.events.begin(); iter != sim.events.end(); ++iter) {
        std::cout << "Switch: " << iter->e->isExplicit() << std::endl;
        std::cout << iter->e->toString() << std::endl;
        if (iter->e->isExplicit()) {
            os << "    s.insert(";
            os << ( iter->e->args[0]->isTime() ? iter->e->args[1]->toString() : iter->e->args[0]->toString() );
            os << ");" << std::endl;
        }
    }
    os << "    return s;" << std::endl;
    os << "}" << std::endl << std::endl;

    // TODO: fix
    os << "unsigned int _dimImplicitSwitches() { static const unsigned int dimImplicitSwitches = " << "0; return dimImplicitSwitches; }" << std::endl;

    // ODE initial values and initial time to evaluate global parameter initial values
    //double * InitialValues = sim->GetDEInitialValues();

    os << "const double *_P_init() { static const double P_init[] = { ";
    unsigned int dimP = 0;
    for(auto iter = sim.parameters.begin(); iter!=sim.parameters.end(); ++iter) {
        if(dimP!=0)
            os << ", ";
        if (++dimP % 100 == 0)
            os << std::endl << "                                                           ";
        os << toString(sim.vecNodes[iter->id]->value); // use node here, because it is evaluated;  toString for correct nan/inf output
    }
    os << "}; return P_init; }" << std::endl;

    os << "const double *_Y_sca() { static const double Y_sca[] = { ";
    unsigned int dimY = 0;
    for(auto iter = sim.variables.begin(); iter!=sim.variables.end(); ++iter) {
        if(dimY!=0)
            os << ", ";
        if (++dimY % 100 == 0)
            os << std::endl << "                                                           ";
        os << toString(iter->scalefactor); // toString for correct nan/inf output
    }
    os << "}; return Y_sca; }" << std::endl;

    // TODO: use real map instead to improve ODEParamValues?
    os << "const unsigned int *_P_map() { static const unsigned int P_map[] = {";

    dimP = 0;
    for(auto iter = sim.parameters.begin(); iter!=sim.parameters.end(); ++iter) {
        if(dimP!=0)
            os << ", ";
        if (++dimP % 100 == 0)
            os << std::endl << "                                                  ";
        os << sim.vecNodes[iter->id]->id;
    }
    os << "}; return P_map; }" << std::endl << std::endl;

    os << "const unsigned int *_O_map() { static const unsigned int O_map[] = {";
    unsigned int dimO = 0;
    for(auto iter = sim.observers.begin(); iter!=sim.observers.end(); ++iter) {
        if(dimO!=0)
            os << ", ";
        if (++dimO % 100 == 0)
            os << std::endl << "                                                  ";
        os << sim.vecNodes[iter->id]->id;
    }
    os << "}; return O_map; }" << std::endl << std::endl;

    os << "const unsigned int *_Y_map() { static const unsigned int Y_map[] = {";
    dimY = 0;
    for(auto iter = sim.variables.begin(); iter!=sim.variables.end(); ++iter) {
        if(dimY!=0)
            os << ", ";
        if (++dimY % 100 == 0)
            os << std::endl << "                                                  ";
        os << sim.vecNodes[iter->id]->id;
    }
    os << "}; return Y_map; }" << std::endl << std::endl;

    os << "const unsigned int *_S_init() { static const unsigned int S_init[] = {";
    unsigned int dimS = 0;
    for(auto iter = sim.events.begin(); iter!=sim.events.end(); ++iter) {//sim->Switches().size()
        if(dimS!=0)
            os << ", ";
        if (++dimS % 100 == 0)
            os << std::endl << "                                                    ";
        os << 1;
    }
    //TODO: fix
    //	os << std::endl << "                                                    ";
    //	for (std::map<unsigned int, formulaParameterInfo >::const_iterator iter = formulaParameterIDs.begin(); iter != formulaParameterIDs.end(); ++iter)
    //	{
    //		os << ", " << iter->second.initialIndex;
    //	}
    os << "}; return S_init; }" << std::endl << std::endl;

    os << "unsigned int _dimP() { static const unsigned int dimP = " << dimP << "; return dimP; }" << std::endl;
    os << "unsigned int _dimP_free() { static const unsigned int dimP_free = " << sim.dimP_free << "; return dimP_free; }" << std::endl;
    os << "unsigned int _dimS() { static const unsigned int dimS = " << sim.events.size() << "; return dimS; }" << std::endl;//sim->Switches().size() + formulaParameterIDs.size()
    os << "unsigned int _dimObservers() { static const unsigned int dimObservers = " << sim.observers.size() << "; return dimObservers; }" << std::endl;
    os << "uint64_t _hash() { static const uint64_t hash = " << sim.hash << "U; return hash; }"  << std::endl  << std::endl;

    //	os << jacConstants;

    os << "}" << std::endl << std::endl;

    // function ODEOptions
    //const DESolverProperties & solverProperties = sim->GetSolver().GetSolverProperties();
    os << "void CLASSNAME(MODELNAME)::ODEOptions(double &relTol, double &absTol, double &hInit, double &hMin, double &hMax, long int &maxSteps, bool &useJac, const double *__restrict__ P) {" << std::endl;
    if(sim.solver.relTolID>=0)
        os << "    relTol = "   << sim.vecNodes[sim.solver.relTolID]->e->toString() << ";" << std::endl;
    if(sim.solver.absTolID>=0)
        os << "    absTol = "   << sim.vecNodes[sim.solver.absTolID]->e->toString() << ";" << std::endl;
    if(sim.solver.hInitID>=0)
        os << "    hInit  = "   << sim.vecNodes[sim.solver.hInitID]->e->toString()  << ";" << std::endl;
    if(sim.solver.hMinID>=0)
        os << "    hMin   = "   << sim.vecNodes[sim.solver.hMinID]->e->toString()   << ";" << std::endl;
    if(sim.solver.hMaxID>=0)
        os << "    hMax   = "   << sim.vecNodes[sim.solver.hMaxID]->e->toString()   << ";" << std::endl;
    if(sim.solver.maxStepsID>=0)
        os << "    maxSteps = " << sim.vecNodes[sim.solver.maxStepsID]->e->toString() << ";" << std::endl; // (int)(sim.vecNodes[sim.solver.maxStepsID]->e->getValue()+0.5)
    if(sim.solver.useJacID>=0)
        os << "    useJac = "   << sim.vecNodes[sim.solver.useJacID]->e->toString()   << ";" << std::endl; // (int)(sim.vecNodes[sim.solver.useJacID]->e->getValue()+0.5)
    //TODO
    //os << "                 'MaxOrder'," << 5 << ", ..." << endl;//solverProperties.GetMaxOrd()<<", ..."<<endl;
    os << "}" << std::endl << std::endl;
}

void writeODERHSFunction(Simulation &sim, std::ostream &os) {
    os << "void CLASSNAME(MODELNAME)::ODERHSFunction(const double Time, const double *__restrict__ y, const double *__restrict__ P, const unsigned int *__restrict__ S, double *__restrict__ dy) {" << std::endl;

    // collect dependencies on parameters
    sim.resetDependencies();
    for (auto iter = sim.variables.begin(); iter!=sim.variables.end(); ++iter)
        iter->e->args[1]->getDependencies(sim.isDependent);
    writeParamLocal(sim, os);

    unsigned int indexY = 0;
    for(auto iter = sim.variables.begin(); iter!=sim.variables.end(); ++iter) {
        os << "    dy[" << indexY++ << "] = ";
        os << iter->e->args[1]->toString();
        // or to string
        os << "; // " << std::endl;
    }
    os << "}" << std::endl << std::endl;
}

void writeODESwitches(Simulation &sim, std::ostream &os) {
    // explicit switches
    os << "bool CLASSNAME(MODELNAME)::ODEExplicitSwitch(const double Time, double *__restrict__ y, double *__restrict__ P, unsigned int *__restrict__ S) {" << std::endl;
    os << "    bool switchUpdate = false;" << std::endl;

    // collect dependencies on parameters
    sim.resetDependencies();
    for(auto iter = sim.events.begin(); iter != sim.events.end(); ++iter) {
        if (iter->e->isExplicit()) {
            iter->e->getDependencies(sim.isDependent);
        }
    }
    writeParamLocal(sim, os);

    unsigned int switchIndex = 0;
    for(auto eIter = sim.events.begin(); eIter!=sim.events.end(); ++eIter) {
        os << "    if(";
        if(eIter->oneTime)
            os << "S[" << switchIndex << "] && ";
        os << eIter->e->toString() << ") {" << std::endl;
        if(eIter->oneTime)
            os << "        S[" << switchIndex << "] = 0;" << std::endl;

        for(auto aIter=eIter->listAssignments.begin(); aIter!=eIter->listAssignments.end(); ++aIter) {
//            if(sim.vecNodes[aIter->oldid]->rhs!=nullptr) {
//                // state
//            }
            os << "        " << sim.vecNodes[aIter->oldid]->name << " = " << sim.vecNodes[aIter->newid]->e->toString() << ";" << std::endl;
        }

        os << "        switchUpdate = true;" << std::endl;
        os << "    }" << std::endl << std::endl;
        ++switchIndex;
    }

    os << "    return switchUpdate;" << std::endl;
    os << "}" << std::endl << std::endl;

    //	// explicit timepoints
    //	os << "std::set<double> _explicitSwitches(const double *P) {" << std::endl;
    //	os << "    std::set<double> s;" << std::endl;
    //	//os << writeParamLocal(usedIDs, os);// TODO: include
    //
    //	for(auto iter = sim.events.begin(); iter != sim.events.end(); ++iter) {
    //		//write
    //	}
    //
    //	os << "    return s;" << std::endl;
    //	os << "}" << std::endl << std::endl;

    os << "void CLASSNAME(MODELNAME)::ODEImplicitSwitch(const double Time, double *__restrict__ y, double *__restrict__ P, unsigned int *__restrict__ S, double *__restrict__ gout) {" << std::endl;
    os << "}" << std::endl << std::endl;
}

void writeODETableParameters(Simulation &sim, std::ostream &os) {
    //	unsigned int i;
    //	vector <Parameter *> tableParameters;
    //
    //	//oss << "function SetupTableParameters" << endl << endl;
    //
    //	//---- assign table formula to every table parameter
    //	for (i = 0; i < sim->Parameters().size(); i++)
    //	{
    //		Parameter * param = sim->Parameters()[i];
    //		if (!param->IsTable())
    //			continue;
    //
    //		// TODO: check removal
    //		//string paramName = param->GetShortUniqueName();
    //
    //		//declare as global
    //		//oss << "    global " << paramName << ";" << endl;
    //
    //		//assign table function
    //		//outfile << "    " << paramName << " = @(Time,y) " << param->TableFunctionNameForCpp() << "(Time,y);" << endl << endl;
    //
    //		//cache table parameter
    //		tableParameters.push_back(param);
    //	}
    //
    //	//---- now write table functions itself
    //	for (size_t i = 0; i < tableParameters.size(); i++)
    //		tableParameters[i]->writeTableFunctionForCpp(os);
}

struct JacFormula {
    int row;
    int col;
    std::string str;
    JacFormula(int _row, int _col, std::string _str) : row(_row), col(_col), str(_str) {}
};

void writeSparseStructure(Simulation &sim, const std::string &nameExt, const std::vector<JacFormula> &vecJacFormulas,
        const int noRows, const int noCols, std::ostream &dense, std::ostream &sparse, std::ostream &osConsts) {
    // write dense matrix
    dense << "void CLASSNAME(MODELNAME)::ODEJacDense" << nameExt << "(const double Time, const double *__restrict__ y, const double *__restrict__ P, const unsigned int *__restrict__ S, double *__restrict__ J) {" << std::endl;
    writeParamLocal(sim, dense);
    for (std::vector<JacFormula>::const_iterator iter = vecJacFormulas.begin(); iter != vecJacFormulas.end(); ++iter)
        dense << "    J[" << iter->col*noRows + iter->row << "] = " << iter->str << ";" << std::endl;
    dense << "}" << std::endl << std::endl;

    // write sparse matrix
    sparse << "void CLASSNAME(MODELNAME)::ODEJacSparse" << nameExt << "(const double Time, const double *__restrict__ y, const double *__restrict__ P, const unsigned int *__restrict__ S, double *__restrict__ J) {" << std::endl;
    writeParamLocal(sim, sparse);
    for (size_t i = 0; i < vecJacFormulas.size(); ++i)
        sparse << "    J[" << i << "] = " << vecJacFormulas[i].str << ";" << std::endl;
    sparse << "}" << std::endl << std::endl;

    // write sparse structure
    osConsts << "int _dimNNZ" << nameExt << "() { static int dimNNZ" << nameExt << " = " << vecJacFormulas.size() << "; return dimNNZ" << nameExt << "; }" << std::endl;

    std::vector<int> sparseI(noCols+1);
    sparseI[0] = 0;
    for (size_t i = 0; i < vecJacFormulas.size(); i++)
        ++sparseI[vecJacFormulas[i].col+1];
    std::partial_sum(sparseI.begin(), sparseI.end(), sparseI.begin());

    osConsts << "const int* _jacRows" << nameExt << "() { static const int jacRows" << nameExt << "[] = {";
    for (size_t i=0; i < vecJacFormulas.size();) // increment below
    {
        if (i==0)
            osConsts << vecJacFormulas[i].row; // first without comma
        else
            osConsts << ", " << vecJacFormulas[i].row;
        if ((++i) % 100 == 0)
            osConsts << std::endl << "                                                      ";
    }
    osConsts << "}; return jacRows" << nameExt << "; }" << std::endl;

    osConsts << "const int* _jacCols" << nameExt << "() { static const int jacCols" << nameExt << "[] = {";
    for (size_t i = 0; i < sparseI.size();)  // increment below
    {
        if (i == 0)
            osConsts << sparseI[i]; // first without comma
        else
            osConsts << ", " << sparseI[i];
        if ((++i) % 100 == 0)
            osConsts << std::endl << "                                                      ";
    }
    osConsts << "}; return jacCols" << nameExt << "; }" << std::endl;
}

void writeODEJacobian(Simulation &sim, std::ostream &dense, std::ostream &sparse, std::ostream &osConsts) {
    std::ostringstream ossFormula;
    formatStream(ossFormula);
    std::vector<JacFormula> vecJacFormulas;
    const int noRows = sim.variables.size();
    const int noCols = sim.variables.size();
    vecJacFormulas.reserve(noRows*noCols);

    // collect dependencies on parameters
    sim.resetDependencies();
    unsigned int col = 0;
    for(auto colIter = sim.variables.begin(); colIter != sim.variables.end(); ++colIter) { // ! col wise !
        unsigned int row = 0;
        for(auto rowIter = sim.variables.begin(); rowIter != sim.variables.end(); ++rowIter) {
//            if(rowIter->rhs==nullptr) {
//                SimInterface::printfPtr("[ERROR] RHS not initialized. Id: %u\n", rowIter->id);
//                exit(0);
//            }

            Expression *e = rowIter->e->args[1]->sd(colIter->id);//rowIter->e->args[1]->newPtr();
            e = e->simplify();

//            if(col==210 && row==211) {
//                Expression *e = rowIter->e->args[1]->sd(colIter->id);//iter->diff(col);
//                std::cout << e->toString() << std::endl;;
//                e = e->simplify();
//                std::cout << e->toString() << std::endl;;
//
//                std::cout << rowIter->id << "|" << colIter->id << std::endl;
//                std::cout << rowIter->e->args[1]->toString() << std::endl;
//                exit(0);
//            }

            // even the dense Jacobian need non-zero elements only
            if (!e->isZero()) {
                e->getDependencies(sim.isDependent);
                //if (CheckNewParameter(&f, mapNewP))
                //	f->InsertNewParameters(mapNewP);

                ossFormula.clear();
                ossFormula.str("");
                ossFormula << e->toString();
                vecJacFormulas.push_back(JacFormula(row, col, ossFormula.str()));
            }
            ++row;

            // clean up
            e->delPtr();
        }
        ++col;
    }

    // nzz structure, write rows directly, use ossFormula to cache column positions
    writeSparseStructure(sim, "", vecJacFormulas, noRows, noCols, dense, sparse, osConsts);
}

void writeSensJacobian(Simulation &sim, std::ostream &dense, std::ostream &sparse, std::ostream &osConsts) {
    // fp
    std::ostringstream ossFormula;
    formatStream(ossFormula);
    std::vector<JacFormula> vecJacFormulas;
    int noCols = sim.dimP_free;
    int noRows = sim.variables.size();
    vecJacFormulas.reserve(noCols*noRows);

    // collect dependencies on parameters
    sim.resetDependencies();
    for (int col = 0; col < noCols; col++) { // ! col wise !
        int row = 0;
        for(auto iter = sim.variables.begin(); iter != sim.variables.end(); ++iter) {
            if(iter->e==nullptr) {
                SimInterface::printfPtr("[ERROR] Id %u not initialized.\n", iter->id);
                exit(0);
            }
            //			Expression *e = sim.vecExpressions[iter->id]->sd(col);//iter->second->DE_Jacobian(-vecParameters[col]);
            //
            ////			e = e->simplify();
            //
            //			// even the dense Jacobian need non-zero elements only
            //			if (!e->isZero()) {
            //				e->getDependencies(deps);
            //				//if (CheckNewParameter(&f, mapNewP))
            //				//	f->InsertNewParameters(mapNewP);
            //
            //				ossFormula.clear();
            //				ossFormula.str("");
            //				ossFormula << e->toString();
            //				vecJacFormulas.push_back(JacFormula(row, col, ossFormula.str()));
            //			}
            //			++row;

            // clean up
            //			delete e;
        }
    }

    // nzz structure, write rows directly, use ossFormula to cache column positions
    writeSparseStructure(sim, "_fp", vecJacFormulas, noRows, noCols, dense, sparse, osConsts);

    // ox
    ossFormula.clear();
    ossFormula.str("");
    vecJacFormulas.clear();
    noCols = sim.variables.size();
    noRows = sim.observers.size();
    vecJacFormulas.reserve(noCols*noRows);

    sim.resetDependencies();
    for (int col = 0; col < noCols; col++) // ! col wise !
    {
        int row = 0;
        for (auto iter = sim.observers.begin(); iter != sim.observers.end(); ++iter) {
            if(iter->e==nullptr) {
                SimInterface::printfPtr("[ERROR] Id %u not initialized.\n", iter->id);
                exit(0);
            }
            //			Expression *e = sim.vecExpressions[iter->id]->sd(col); // obs[row]->DE_Jacobian(col);
            //
            ////			e = e->simplify();
            //
            //			// even the dense Jacobian need non-zero elements only
            //			if (!e->isZero()) {
            //				e->getDependencies(deps);
            //				//if (CheckNewParameter(&f, mapNewP))
            //				//	f->InsertNewParameters(mapNewP);
            //
            //				ossFormula.clear();
            //				ossFormula.str("");
            //				ossFormula << e->toString();
            //				vecJacFormulas.push_back(JacFormula(row, col, ossFormula.str()));
            //			}
            //			++row;

            // clean up
            //			delete e;
        }
    }

    // nzz structure, write rows directly, use ossFormula to cache column positions
    writeSparseStructure(sim, "_ox", vecJacFormulas, noRows, noCols, dense, sparse, osConsts);

    // op
    ossFormula.clear();
    ossFormula.str("");
    vecJacFormulas.clear();
    noCols = sim.dimP_free;
    noRows = sim.observers.size();
    vecJacFormulas.reserve(noCols*noRows);

    sim.resetDependencies();
    for (int col = 0; col < noCols; col++) // ! col wise !
    {
        int row = 0;
        for (auto iter = sim.observers.begin(); iter != sim.observers.end(); ++iter) {
            if(iter->e==nullptr) {
                SimInterface::printfPtr("[ERROR] Id %u not initialized.\n", iter->id);
                exit(0);
            }
            //			Expression *e = sim.vecExpressions[iter->id]->sd(col); //obs[j]->DE_Jacobian(-vecParameters[col]);
            //
            ////			e = e->simplify();
            //
            //			// even the dense Jacobian need non-zero elements only
            //			if (!e->isZero()) {
            //				e->getDependencies(deps);
            //				//if (CheckNewParameter(&f, mapNewP))
            //				//	f->InsertNewParameters(mapNewP);
            //
            //				ossFormula.clear();
            //				ossFormula.str("");
            //				ossFormula << e->toString();
            //				vecJacFormulas.push_back(JacFormula(row, col, ossFormula.str()));
            //			}
            //			++row;

            // clean up
            //			delete e;
        }
    }

    // nzz structure, write rows directly, use ossFormula to cache column positions
    writeSparseStructure(sim, "_op", vecJacFormulas, noRows, noCols, dense, sparse, osConsts);
}

void cppExport(Simulation &sim, const std::string &exportDir) {
    std::ostringstream ossBufferExport;
    std::ostringstream ossBufferJacDense;
    std::ostringstream ossBufferJacSparse;
    formatStream(ossBufferExport);
    formatStream(ossBufferJacDense);
    formatStream(ossBufferJacSparse);

    ossBufferExport << "#define MODELNAME " << sim.name << std::endl;
    ossBufferExport << "#include \"ModelDerived.h\"" << std::endl << std::endl;

//    ossBufferJacDense << "#define MODELNAME " << sim.name << std::endl;
//    ossBufferJacDense << "#include \"ModelDerived.h\"" << std::endl << std::endl;
//
//    ossBufferJacSparse << "#define MODELNAME " << sim.name << std::endl;
//    ossBufferJacSparse << "#include \"ModelDerived.h\"" << std::endl << std::endl;

    //writeODEDefinitions(sim, ossBufferExport);
    writeODEInitialValues(sim, ossBufferExport);
    writeODERHSFunction(sim, ossBufferExport);
    writeODESwitches(sim, ossBufferExport);
    writeODETableParameters(sim, ossBufferExport);

    writeODEParameterAndObserver(sim, ossBufferExport);

    ossBufferExport << "namespace {" << std::endl;
    writeODEJacobian(sim, ossBufferJacDense, ossBufferJacSparse, ossBufferExport);
    writeSensJacobian(sim, ossBufferJacDense, ossBufferJacSparse, ossBufferExport);
    ossBufferExport << "} // namespace" << std::endl << std::endl;

    ossBufferExport << "#define CONVERT(x) #x" << std::endl;
    ossBufferExport << "#define STRING(x) CONVERT(x)" << std::endl;
    ossBufferExport << "CLASSNAME(MODELNAME)::CLASSNAME(MODELNAME)() : Model(STRING(MODELNAME), _dimY(), _dimP(), _dimP_free(), _dimS(), _dimOutputTimePoints()," << std::endl;
    ossBufferExport << "                                              _dimImplicitSwitches(), _dimObservers(), _dimNNZ(), _dimNNZ_fp(), _dimNNZ_ox(), _dimNNZ_op()," << std::endl;
    ossBufferExport << "                                              _P_map(), _O_map(), _Y_map(), _Y_sca(), _P_init(), _S_init(), _outputTimePoints()," << std::endl;
    ossBufferExport << "                                              _jacRows(), _jacCols(), _jacRows_fp(), _jacCols_fp(), _jacRows_op(), _jacCols_op()," << std::endl;
    ossBufferExport << "                                              _jacRows_op(), _jacCols_op(), _hash(), _explicitSwitches) {" << std::endl;
    ossBufferExport << "}" << std::endl;
    ossBufferExport << "static CLASSNAME(MODELNAME) MODELNAME;" << std::endl;

    // buffered output finished, write buffer
    std::string filenameBase = exportDir + sim.name;

    std::ofstream outfile;
    std::string filename;

    filename = filenameBase + ".cpp";
    outfile.open(filename.c_str());
    outfile << ossBufferExport.str();
    outfile << ossBufferJacDense.str();
    outfile << ossBufferJacSparse.str();
    outfile.close();

//    filename = filenameBase + "Base.cpp";
//    outfile.open(filename.c_str());
//    outfile << ossBufferExport.str();
//    outfile.close();
//    filename = filenameBase + "JacDense.cpp";
//    outfile.open(filename.c_str());
//    outfile << ossBufferJacDense.str();
//    outfile.close();
//    filename = filenameBase + "JacSparse.cpp";
//    outfile.open(filename.c_str());
//    outfile << ossBufferJacSparse.str();
}


