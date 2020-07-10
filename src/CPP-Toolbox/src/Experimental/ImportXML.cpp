/*
 * ImportXML.cpp
 *
 */

#include <iostream>
#include <list>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <iterator> // std::distance

#include "ImportXML.h"
#include "Expression.h"
#include "Integrator.h" // outputString

//#include "FuncParser/FuncParser.h"
//using namespace FuncParserNative;

using namespace rapidxml;


struct parameter {

};

struct event {

};

Node Simulation::Time;

void Node::print() {
    SimInterface::printfPtr("Id:       %u\n", id);
    SimInterface::printfPtr("Expr:     %p\n", e);
    SimInterface::printfPtr("Name:     %s\n", name.c_str());
    SimInterface::printfPtr("Path:     %s\n", path.c_str());
    SimInterface::printfPtr("Unit:     %s\n", unit.c_str());
    SimInterface::printfPtr("Equation: %s\n", equation.c_str());
    SimInterface::printfPtr("fId:      %u\n", formulaId);
    SimInterface::printfPtr("Value:    %f\n", value);
    SimInterface::printfPtr("Scaling:  %f\n", scalefactor);
    SimInterface::printfPtr("Persist:  %u\n", persistent);
    SimInterface::printfPtr("Variable: %u\n", canBeVaried);
    SimInterface::printfPtr("noNegVal: %u\n", negativeValuesAllowed);
    SimInterface::printfPtr("OneTime:  %u\n", oneTime);
    SimInterface::printfPtr("RHS size: %u\n", listRHS.size());
    SimInterface::printfPtr("EqRef s.: %u\n", mapEqRef.size());
    SimInterface::printfPtr("Assign s: %u\n", listAssignments.size());;
}

void Simulation::getTimepoints(xml_node<> *n, std::set<double> &setT) {
    for (xml_node<> *child = n->first_node(); child; child = child->next_sibling()) {
        if( strcmp(child->name(), "OutputIntervalList" )==0 ) {
            for (xml_node<> *grandchild = child->first_node(); grandchild; grandchild = grandchild->next_sibling()) {
                if( strcmp( grandchild->name(), "OutputInterval" )==0 ) {
                    if ( strcmp(grandchild->first_attribute("distribution")->value(), "Uniform")==0 ) {
                        double t_start=0;
                        double t_end=0;
                        unsigned int N=1;
                        for (xml_node<> *ggc = grandchild->first_node(); ggc; ggc = ggc->next_sibling()) {
                            if      ( strcmp( ggc->name(), "StartTime" )==0 ) {
                                t_start = atof(ggc->value());
                            }
                            else if ( strcmp( ggc->name(), "EndTime" )==0 ) {
                                t_end = atof(ggc->value());
                            }
                            else if ( strcmp( ggc->name(), "NumberOfTimePoints" )==0 ) {
                                N = atoi(ggc->value());
                            }
                        }
                        double d = (t_end - t_start)/(N-1);
                        for(unsigned int i=0;i<N;i++)
                            setT.emplace(t_start+i*d);
                    }
                    else
                        SimInterface::printfPtr("[ERROR] Unknown output interval distribution: %s\n", grandchild->first_attribute("distribution")->value());
                }
                else
                    SimInterface::printfPtr("[ERROR] Expected output interval node, got %s\n", grandchild->name());
            }
        }
        else
            SimInterface::printfPtr("[ERROR] Expected output interval list node, got %s\n", child->name());
    }
}

void Simulation::getSolver(xml_node<> *n, SolverOptions &solver) {
    xml_attribute<> * att = n->first_attribute();
    if(att) {
        if(strcmp(att->name(), "name")==0)
            solver.name = att->value();
    }
    for (xml_node<> *child = n->first_node(); child; child = child->next_sibling()) {
        if( strcmp(child->name(), "RelTol" )==0 )
            solver.relTolID = atol(child->first_attribute()->value());
        else if( strcmp(child->name(), "AbsTol" )==0 )
            solver.absTolID = atol(child->first_attribute()->value());
        else if( strcmp(child->name(), "H0" )==0 )
            solver.hInitID = atol(child->first_attribute()->value());
        else if( strcmp(child->name(), "HMin" )==0 )
            solver.hMinID = atol(child->first_attribute()->value());
        else if( strcmp(child->name(), "HMax" )==0 )
            solver.hMaxID = atol(child->first_attribute()->value());
        else if( strcmp(child->name(), "MxStep" )==0 )
            solver.maxStepsID = atol(child->first_attribute()->value());
        else if( strcmp(child->name(), "UseJacobian" )==0 )
            solver.useJacID = atol(child->first_attribute()->value());
        else
            SimInterface::printfPtr("[ERROR] Expected solver attribute, got %s\n", child->name());
    }
}

void Simulation::getList(xml_node<> *n, std::list<Node> &l) {
    for (xml_node<> *child = n->first_node(); child; child = child->next_sibling()) {
        Node newNode;
        for (xml_attribute<> *att = child->first_attribute(); att; att = att->next_attribute()) {
            if      (strcmp(att->name(), "id")==0) {
                newNode.id = atoi(att->value());
                if(newNode.id>maxID)
                    maxID = newNode.id;
            }
            else if (strcmp(att->name(), "entityId")==0) {
                // ignore for now
            }
            else if (strcmp(att->name(), "name")==0) {
                newNode.name = att->value();
            }
            else if (strcmp(att->name(), "path")==0) {
                newNode.path = att->value();
            }
            else if (strcmp(att->name(), "unit")==0) {
                newNode.unit = att->value();
            }
            else if (strcmp(att->name(), "formulaId")==0) {
                newNode.formulaId = atoi(att->value());
            }
            else if (strcmp(att->name(), "conditionFormulaId")==0) {
                newNode.formulaId = atoi(att->value()); // treat equally to formulaId
            }
            else if (strcmp(att->name(), "initialValueFormulaId")==0) {
                newNode.formulaId = atoi(att->value()); // treat equally to formulaId
            }
            else if (strcmp(att->name(), "value")==0) {
                newNode.value = atof(att->value());
            }
            else if (strcmp(att->name(), "persistable")==0) {
                newNode.persistent = atoi(att->value());
            }
            else if (strcmp(att->name(), "canBeVaried")==0) { // parameter
                newNode.canBeVaried = atoi(att->value());
            }
            else if (strcmp(att->name(), "negativeValuesAllowed")==0) { // variable
                newNode.negativeValuesAllowed = atoi(att->value());
            }
            else if (strcmp(att->name(), "oneTime")==0) { // event
                newNode.oneTime = atoi(att->value());
            }
            else {
                SimInterface::printfPtr("[ERROR] Unknown attribute %s\n", att->name());
            }
        }
        for (xml_node<> *grandchild = child->first_node(); grandchild; grandchild = grandchild->next_sibling()) {
            if( strcmp( grandchild->name(), "Equation" )==0 ) {
                newNode.equation = grandchild->value();
            }
            else if( strcmp( grandchild->name(), "ReferenceList" )==0 ) {
                for (xml_node<> *ggc = grandchild->first_node(); ggc; ggc = ggc->next_sibling()) {
                    newNode.mapEqRef.emplace(ggc->first_attribute("alias")->value(), atoi(ggc->first_attribute("id")->value()));
                }
            }
            else if( strcmp( grandchild->name(), "ScaleFactor" )==0 ) {
                newNode.scalefactor = atof(grandchild->value());
            }
            else if( strcmp( grandchild->name(), "RHSFormulaList" )==0 ) {
                for (xml_node<> *ggc = grandchild->first_node(); ggc; ggc = ggc->next_sibling()) {
                    newNode.listRHS.push_back(atoi(ggc->first_attribute("id")->value()));
                }
            }
            else if( strcmp( grandchild->name(), "AssignmentList" )==0 ) {
                for (xml_node<> *ggc = grandchild->first_node(); ggc; ggc = ggc->next_sibling()) {
                    unsigned int oldid = atoi(ggc->first_attribute("objectId")->value());
                    unsigned int newid = atoi(ggc->first_attribute("newFormulaId")->value());
                    bool useAsValue = atoi(ggc->first_attribute("useAsValue")->value());
                    newNode.listAssignments.emplace_back(oldid, newid, useAsValue);
                }
            }
            else{
                SimInterface::printfPtr("[ERROR] Unknown node %s\n", grandchild->name());
            }
        }
        if(newNode.id==UINT_MAX) {
            newNode.print();
            exit(0);
        }
        l.push_back(newNode);
    }
}

Simulation::Simulation(xml_node<> *n) {
    std::list<Node> listNodes;

    for (xml_node<> *child = n->first_node(); child; child = child->next_sibling()) {
        if      (strcmp(child->name(), "ObserverList")==0) {
            getList(child, observers);
        }
        else if (strcmp(child->name(), "EventList")==0) {
            getList(child, events);
        }
        else if (strcmp(child->name(), "FormulaList")==0) {
            getList(child, formulas);
        }
        else if (strcmp(child->name(), "VariableList")==0) {
            getList(child, variables);
        }
        else if (strcmp(child->name(), "ParameterList")==0) {
            getList(child, parameters);
        }
        else if (strcmp(child->name(), "Solver")==0) {
            getSolver(child, solver);
        }
        else if (strcmp(child->name(), "OutputSchema")==0) {
            getTimepoints(child, outputTimepoints);
        }
        else {
            SimInterface::printfPtr("[ERROR] Unknown node %s\n", child->name());
        }
    }

    // consolidate into single vector for access by id
    vecNodes.assign(maxID+1, nullptr);
    vecNodes[0] = &Time;// explicitly insert time (id==0)
    vecNodes[0]->id = 0;
    vecNodes[0]->name = "Time";
    vecNodes[0]->value = 0.0;
    vecNodes[0]->e = Expression::fromString(&(vecNodes[0]->name));
    vecNodes[0]->e->id = 0;
    vecNodes[0]->e->args.clear();
    vecNodes[0]->e->isSimplified = true;

    for(auto iter = observers.begin(); iter != observers.end(); ++iter)
        if(vecNodes[iter->id]==nullptr)
            vecNodes[iter->id] = &(*iter);
        else
            SimInterface::printfPtr("[ERROR] Duplicate id in xml: %u.\n", iter->id);
    for(auto iter = events.begin(); iter != events.end(); ++iter)
        if(vecNodes[iter->id]==nullptr)
            vecNodes[iter->id] = &(*iter);
        else
            SimInterface::printfPtr("[ERROR] Duplicate id in xml: %u.\n", iter->id);
    for(auto iter = formulas.begin(); iter != formulas.end(); ++iter)
        if(vecNodes[iter->id]==nullptr)
            vecNodes[iter->id] = &(*iter);
        else
            SimInterface::printfPtr("[ERROR] Duplicate id in xml: %u.\n", iter->id);
    for(auto iter = variables.begin(); iter != variables.end(); ++iter)
        if(vecNodes[iter->id]==nullptr)
            vecNodes[iter->id] = &(*iter);
        else
            SimInterface::printfPtr("[ERROR] Duplicate id in xml: %u.\n", iter->id);
    for(auto iter = parameters.begin(); iter != parameters.end(); ++iter)
        if(vecNodes[iter->id]==nullptr)
            vecNodes[iter->id] = &(*iter);
        else
            SimInterface::printfPtr("[ERROR] Duplicate id in xml: %u.\n", iter->id);
}

void Simulation::printStats() {
    SimInterface::printfPtr("[INFO] Observers:         %u\n", observers.size());
    SimInterface::printfPtr("[INFO] Variables:         %u\n", variables.size());
    SimInterface::printfPtr("[INFO] Parameters:        %u\n", parameters.size());
    SimInterface::printfPtr("[INFO] Events:            %u\n", events.size());
    SimInterface::printfPtr("[INFO] Formulas:          %u\n", formulas.size());
    //	SimInterface::printfPtr("[INFO] Solvers:           %u\n", solvers.size());
    SimInterface::printfPtr("[INFO] Output timepoints: %u\n", outputTimepoints.size());
}

void mark(Simulation &sim, std::vector<unsigned int> &marked, Node &n) {
    if(marked[n.id]==2) {
        SimInterface::printfPtr("[ERROR] Not a DAG. Id: ", n.id);;
        exit(0);
    }
    if(!marked[n.id]) {
        marked[n.id] = 2; // current recursion

        if(sim.vecNodes[n.id]->e->id!=n.id) { // variable only
            std::deque<unsigned int> deps;
            sim.vecNodes[n.id]->e->getDependencies(deps);
            if(sim.vecNodes[n.id]->canBeSwitched)
                for(auto aIter = sim.vecNodes[n.id]->listAssignments.begin();  aIter!=sim.vecNodes[n.id]->listAssignments.end(); ++aIter)
                    sim.vecNodes[aIter->newid]->e->getDependencies(deps);

            for(auto iter = deps.begin(); iter!=deps.end(); ++iter)
                mark(sim, marked, *sim.vecNodes[*iter]);
        }

        marked[n.id] = 1; // permanent mark
        n.order = sim.orderIndex;
        ++sim.orderIndex;
    }
}

void sortByDependencies(Simulation &sim) {
    // sim.orderIndex = 0; // do not reset here, instead sort already spefified free params to the beginning
    std::vector<unsigned int> marked(sim.maxID+1, 0);
    sim.orderIndex = 2*sim.maxID;

    sim.parameters.sort(); // move user-defined parameters to the beginning
    auto iter = sim.parameters.begin();
    for(; iter->order != UINT_MAX && iter != sim.parameters.end(); ++iter)
        marked[iter->id] = 1; // flag those as already processed to prevent overwriting order

    for(; iter != sim.parameters.end(); ++iter) {
        if(!marked[iter->id])
            mark(sim, marked, *iter);
    }
    //	for(auto iter = sim.parameters.begin(); iter != sim.parameters.end(); ++iter) {
    //		auto element = std::find(sim.order.begin(), sim.order.end(), iter->id);
    //		size_t index = std::distance(element,sim.order.begin());
    //		iter->order = index;
    //	}
    sim.parameters.sort();
}

// 64bit fnFNV1a, documentation @ http://www.isthe.com/chongo/tech/comp/fnv/#FNV-1a
static uint64_t fnFNV1a(const char* begin, const char* const end) {
    const uint64_t FNV_prime = 1099511628211U; // prime value for 64-bit hash
    uint64_t hash = 14695981039346656037U; // offset value for 64-bit hash

    for (; begin < end; ++begin)
        hash = (hash ^ *begin) * FNV_prime; // bit-wise XOR and multiplication

    return hash;
}

void convertXMLtoCpp(const std::string &xmlFile, const std::string &_outputDir, const std::string &_modelName,
        const std::set<unsigned int> &setObservers, const std::set<unsigned int> &_setParameters, bool reduceVariables) {
    std::list<Simulation> simulations;
    std::string outputDir(_outputDir);
    std::string modelName(_modelName);
    std::set<unsigned int> setParameters(_setParameters);

    // merge for hash
    std::vector<unsigned int> vecParametersAndObservers(_setParameters.size()+setObservers.size());
    std::copy(_setParameters.begin(), _setParameters.end(), std::back_inserter(vecParametersAndObservers));
    std::copy(setObservers.begin(), setObservers.end(), std::back_inserter(vecParametersAndObservers));
    int idMemorySize = vecParametersAndObservers.size()*sizeof(unsigned int);

    // normalize export dir
    if(!outputDir.empty() && !(outputDir.back() == '\\' || outputDir.back() == '/'))
        outputDir.push_back('/');

    // derive modelname if not supplied
    if(modelName.empty()) {
        size_t firstIndex = xmlFile.find_last_of("\\/");
        firstIndex = firstIndex!=std::string::npos ? firstIndex+1 : 0;
        size_t lastIndex = xmlFile.find_last_of(".");
        lastIndex = lastIndex!=std::string::npos ? lastIndex : xmlFile.size();
        if(firstIndex<lastIndex) {
            modelName = xmlFile.substr(firstIndex, lastIndex-firstIndex);
            std::replace(modelName.begin(), modelName.end(), ' ', '_');
        }
        else {
            modelName = "Standard";
        }
    }
    SimInterface::printfPtr("[INFO] Model name: %s\n", modelName.c_str());
    SimInterface::printfPtr("[INFO] Output dir: %s\n", outputDir.c_str());

    // c style non-constant char array for rapidxml
    FILE *f = fopen(xmlFile.c_str(), "rb");
    if(f==nullptr) {
        SimInterface::printfPtr("File not found: %s\n", xmlFile.c_str());
        return;
    }
    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    fseek(f, 0, SEEK_SET);

    char *xmlContent = new char[fsize+idMemorySize+1];// (char*)malloc(fsize + 1);
    fread(xmlContent, fsize, 1, f);
    fclose(f);

    memcpy(&xmlContent[fsize], vecParametersAndObservers.data(), idMemorySize);
    uint64_t hash = fnFNV1a(xmlContent, xmlContent+idMemorySize+1);

    xmlContent[fsize] = 0; // null terminate string
    xml_document<> doc;
    doc.parse<parse_trim_whitespace | parse_no_data_nodes>(xmlContent);
    for (xml_node<> *child = doc.first_node(); child; child = child->next_sibling()) {
        if( strcmp( child->name(), "Simulation" )==0 ) {
            // simulations
            simulations.push_back(Simulation(child));
            simulations.back().name = modelName;
        }
        else {
            SimInterface::printfPtr("[ERROR] Unknown node %s\n", child->name());
        }
    }

    for(auto simIter = simulations.begin(); simIter!=simulations.end(); ++simIter) {
        simIter->printStats();
        simIter->hash = hash;

        // modify canBeVaried
        if(!setParameters.empty()) {
            for(auto pIter = simIter->parameters.begin(); pIter!=simIter->parameters.end(); ++pIter) {
                auto p = setParameters.find(pIter->id);
                if( p!=setParameters.end() ) {
                    pIter->canBeVaried = true;
                    SimInterface::printfPtr("Initializing id %u : %s : %s\n", pIter->id, pIter->name.c_str(), pIter->path.c_str());
                    setParameters.erase(p);
                }
                else
                    pIter->canBeVaried = false;
            }
        }
        SimInterface::printfPtr("Size: %u\n", setParameters.size());

        //		simIter->vecExpressions = std::vector<Expression*>(simIter->maxID+1, nullptr);
        // parse equations
        //FuncParser p;

        for(auto fIter = simIter->formulas.begin(); fIter!=simIter->formulas.end(); ++fIter) {
            simIter->vecNodes[fIter->id]->e = parseExpression(simIter->vecNodes[fIter->id]->equation);

//            std::vector<std::string> VariableNames;
//            std::vector<std::string> ParameterNames;
//            bool CaseSensitive = true;
//            bool LogicOperatorsAllowed = true;
//            double ComparisonTolerance = 1e-20;
//            bool LogicalNumericMixAllowed = true;
//            FuncParserErrorData ED;
//            std::cout << simIter->vecNodes[fIter->id]->equation << std::endl;
//            FuncNode *f = p.Parse(simIter->vecNodes[fIter->id]->equation, VariableNames, ParameterNames, CaseSensitive, LogicOperatorsAllowed, ComparisonTolerance, LogicalNumericMixAllowed, ED);
//            std::cout << ED.GetNumber() << ", " << ED.GetDescription() << std::endl;
//            for(auto vIter = VariableNames.begin(); vIter!=VariableNames.end(); ++vIter)
//                std::cout << *vIter << std::endl;
//            std::cout << "***" << std::endl;
//            for(auto pIter = ParameterNames.begin(); pIter!=ParameterNames.end(); ++pIter)
//                std::cout << *pIter << std::endl;
//            std::cout << "---" << std::endl;
//            return;
        }

        // parse events
        int offsetSwitch = simIter->events.size();
        for(auto eIter = simIter->events.begin(); eIter!=simIter->events.end(); ++eIter) {
            for(auto aIter=eIter->listAssignments.begin(); aIter!=eIter->listAssignments.end(); ++aIter) {
                auto pIterOld = simIter->vecNodes[aIter->oldid];
                auto pIterNew = simIter->vecNodes[aIter->newid];

                //if(pIterOld->e->isConstant() && (pIterNew->e->isConstant() || aIter->useAsValue))

                if(!pIterOld->canBeSwitched) {
                    pIterOld->canBeSwitched = true;
                    //Expression *e = new Expression(OP_SWITCH);
                    //e->name = "S[" + std::to_string(offsetSwitch++) + "]";// set invalid name for debugging
                    //e->delPtr();

                    pIterOld->e = Expression::fromString(&(pIterOld->name));
                    pIterOld->e->id = pIterOld->id;
                    pIterOld->e->type = OP_SWITCH;
                    pIterOld->e->isSimplified = true; // can't simplify those
                    if (pIterOld->formulaId==UINT_MAX) // set initial value
                        pIterOld->e->args[0] = Expression::fromConstant(pIterOld->value);
                    else
                        pIterOld->e->args[0] = simIter->vecNodes[pIterOld->formulaId]->e->newPtr();
                }
                else {
                    //if(!(pIterNew->e->isConstant() || aIter->useAsValue))
                    if (pIterNew->formulaId==UINT_MAX) // set initial value
                        pIterOld->e->args.push_back(Expression::fromConstant(pIterNew->value));
                    else
                        pIterOld->e->args.push_back(simIter->vecNodes[pIterNew->formulaId]->e->newPtr());
                }
            }
            eIter->e = simIter->vecNodes[eIter->formulaId]->e->newPtr();
        }

        std::list<unsigned int> global;
        // parse parameter
        for(auto pIter = simIter->parameters.begin(); pIter!=simIter->parameters.end(); ++pIter) {
            if(pIter->persistent) {
                simIter->observers.emplace_back();//*pIter
                simIter->observers.back().id = pIter->id;
                simIter->observers.back().formulaId = pIter->formulaId;
            }
            if(pIter->canBeVaried) { //|| (pIter->canBeSwitched && pIter->e->args.size()==1)
                // global free parameter if it can be varied or if it can be switched to constants only
                pIter->name = "P[" + std::to_string(simIter->dimP_free++) + "]";
                pIter->order = simIter->orderIndex++; // sort user specified parameters to the beginning later
                pIter->e = Expression::fromString(&(pIter->name));
                pIter->e->id = pIter->id;
                pIter->e->args.clear();
                pIter->e->isSimplified = true; // can't simplify those
            }
            else if(pIter->canBeSwitched) {
                // full switch with new expressions
                //pIter->name
//                pIter->e->setName(&(pIter->name));
//                for(unsigned int i=0; i<pIter->e->args.size(); i++) {
//                    pIter->e->args[i]->args[0] = simIter->vecNodes[i]->e->newPtr();
//                }
//                pIter->order = simIter->orderIndex++ + simIter->maxID; // sort user specified parameters to the beginning later
//                pIter->name = "P[]";// set invalid name for debugging
//                pIter->e = Expression::fromString(&(pIter->name));
//                pIter->e->id = pIter->id;
//                pIter->e->isSimplified = true; // can't simplify those
//                if (pIter->formulaId==UINT_MAX) {// set initial value
//                    std::cout << "Constant" << std::endl;
//                    pIter->e->args[0] = Expression::fromConstant(pIter->value);
//                }
//                else
//                    pIter->e->args[0] = simIter->vecNodes[pIter->formulaId]->e->newPtr();
            }
            else if (pIter->formulaId==UINT_MAX) {
                pIter->name = "P_" + std::to_string(pIter->id) + "_";
                pIter->e = Expression::fromConstant(pIter->value);
            }
            else {
                pIter->name = "P_" + std::to_string(pIter->id) + "_";
                pIter->e = simIter->vecNodes[pIter->formulaId]->e->newPtr();
            }
        }
        SimInterface::printfPtr("[INFO] User defined parameters: %u\n", simIter->dimP_free);

        // observers first
        for(auto oIter = simIter->observers.begin(); oIter!=simIter->observers.end(); ++oIter) {
            oIter->e = simIter->vecNodes[oIter->formulaId]->e->newPtr();
        }

        // parse variables
        for(auto vIter = simIter->variables.begin(); vIter!=simIter->variables.end(); ++vIter) {
            vIter->name = "y[]";
            Expression *rhsExpression = Expression::fromConstant(0);
            for(auto rhs = vIter->listRHS.begin(); rhs!=vIter->listRHS.end(); ++rhs) {
                rhsExpression = new Expression(OP_SUM, rhsExpression, simIter->vecNodes[*rhs]->e->newPtr());
            }
            Expression *initExpression;
            if(vIter->formulaId!=UINT_MAX)
                initExpression = simIter->vecNodes[vIter->formulaId]->e->newPtr();
            else if(!std::isnan(vIter->value))
                initExpression = Expression::fromConstant(vIter->value);
            else {
                SimInterface::printfPtr("[ERROR] Neither initial value nor formula for state. Id: %u.\n", vIter->id);
                return;
            }
            simIter->vecNodes[vIter->id]->e = new Expression(OP_VARIABLE, initExpression, rhsExpression);
            simIter->vecNodes[vIter->id]->e->setName(&(vIter->name));
            simIter->vecNodes[vIter->id]->e->id = vIter->id;
            simIter->vecNodes[vIter->id]->e->isSimplified = true;
        }

        // update names and references
        for(auto fIter = simIter->formulas.begin(); fIter!=simIter->formulas.end(); ++fIter) {
            simIter->vecNodes[fIter->id]->e->replaceReferences(simIter->vecNodes[fIter->id]->mapEqRef, simIter->vecNodes);
//            simIter->vecNodes[fIter->id]->value = simIter->vecNodes[fIter->id]->e->evaluate(&simIter->vecNodes);
        }

        // sort parameters by dependencies and assign missing values to initialize parameters correctly
        sortByDependencies(*simIter);

        for(auto pIter = simIter->parameters.begin(); pIter!=simIter->parameters.end(); ++pIter) {
            if(std::isnan(pIter->value)) {
                if(pIter->formulaId<simIter->vecNodes.size()) {
                    pIter->value = simIter->vecNodes[pIter->formulaId]->e->evaluate(&simIter->vecNodes);
                }
                else
                    pIter->value = pIter->e->evaluate(&simIter->vecNodes);
                if(std::isnan(pIter->value))
                    SimInterface::printfPtr("[WARNING] Parameter is nan. Id: %u\n", pIter->id);
            }
        }

        // simplify
        for(auto iter = simIter->vecNodes.begin(); iter!=simIter->vecNodes.end(); ++iter) {// all expressions
            if((*iter)->e) {
                (*iter)->e = (*iter)->e->simplify();
                if((*iter)->e->type==OP_VARIABLE) {
                    (*iter)->e->args[0] = (*iter)->e->args[0]->simplify();
                    (*iter)->e->args[1] = (*iter)->e->args[1]->simplify();
                }
            }
        }
//        for(auto vIter = simIter->variables.begin(); vIter!=simIter->variables.end(); ++vIter) // all rhs
//            vIter->rhs = vIter->rhs->simplify();

        // reduce and collect dependencies
        unsigned int originalSize = simIter->observers.size();
        unsigned int reducedSize = 0;
        simIter->resetDependencies();
        if(setObservers.empty()) {
            for(auto oIter = simIter->observers.begin(); oIter!=simIter->observers.end(); ++oIter) {
                oIter->e->getDependencies(simIter->isDependent);
            }
            SimInterface::printfPtr("[INFO] Not reducing observers.\n");
        }
        else {
            for(auto oIter = simIter->observers.begin(); oIter!=simIter->observers.end();) {
                if( setObservers.count(oIter->id) ) {
                    oIter->e->getDependencies(simIter->isDependent);
                    ++oIter;
                }
                else {
                    ++reducedSize;
                    simIter->filteredObservers.splice(simIter->filteredObservers.end(), simIter->observers, oIter++);
                }
            }
            SimInterface::printfPtr("[INFO] Reduced observers: %u / %u\n", reducedSize, originalSize);
        }

        originalSize = simIter->variables.size();
        reducedSize = 0;
        simIter->dimY = 0; // reset names
        if(reduceVariables) {
            for(auto vIter = simIter->variables.begin(); vIter!=simIter->variables.end();++vIter) { // mark deps
                if(vIter->canBeSwitched)
                    simIter->isDependent[vIter->id] = true; // don't remove potentially switched states
                if(simIter->isDependent[vIter->id]) {
                    vIter->e->getDependenciesWithRHS(simIter->isDependent);
                }
            }
        }
        else {
            for(auto vIter = simIter->variables.begin(); vIter!=simIter->variables.end();++vIter) { // mark deps
               simIter->isDependent[vIter->id] = true;
               vIter->e->getDependenciesWithRHS(simIter->isDependent);
            }
        }

        for(auto vIter = simIter->variables.begin(); vIter!=simIter->variables.end();) { // remove without deps
//            if(vIter->e->args[1]->isZero()) { // MoBi removes those always
//                vIter->e->type = OP_PARAMETER;
//                vIter->e->args[1]->delPtr();
//                vIter->e->args.pop_back();
//                vIter->canBeVaried = false;
//                vIter->value = vIter->e->args[0]->evaluate(&simIter->vecNodes);
////                if(vIter->canBeVaried) {
////                    vIter->name = "P[" + std::to_string(simIter->dimP_free++) + "]";
////                    vIter->order = simIter->orderIndex++; // sort user specified parameters to the beginning later
////                    vIter->e->setName(&(vIter->name));// = Expression::fromString(&(pIter->name));
////                    vIter->e->args.clear();
////                    vIter->e->isSimplified = true; // can't simplify those
////                }
//                ++reducedSize;
//                simIter->parameters.splice(simIter->parameters.end(), simIter->variables, vIter++); // move to parameters
//            }
//            else
            if(simIter->isDependent[vIter->id]) {
                vIter->name = "y[" + std::to_string(simIter->dimY++) + "]"; // set new names
                ++vIter;
            }
            else {
                ++reducedSize;
                simIter->filteredVariables.splice(simIter->filteredVariables.end(), simIter->variables, vIter++);// simIter->variables.erase(vIter++);
            }
        }
        if(reduceVariables)
            SimInterface::printfPtr("[INFO] Reduced variables: %u / %u\n", reducedSize, originalSize);
        else
            SimInterface::printfPtr("[INFO] Not reducing variables.\n");


        //TODO: filter switches?
        //originalSize = simIter->events.size();
        //reducedSize = 0;

        // mark used parameters in switches
//        for(auto eIter = simIter->events.begin(); eIter != simIter->events.end(); ++eIter) {
//            eIter->e->getDependencies(simIter->isDependent);
//            // potential new dependencies after switches
//            for(auto aIter=eIter->listAssignments.begin(); aIter!=eIter->listAssignments.end(); ++aIter)
//                simIter->vecNodes[aIter->newid]->e->getDependencies(simIter->isDependent);
//        }

        originalSize = simIter->parameters.size();
        reducedSize = 0;
        for(auto pIter = simIter->parameters.begin(); pIter!=simIter->parameters.end();) {
            if(!simIter->isDependent[pIter->id] && !pIter->canBeVaried) {
                ++reducedSize;
                simIter->filteredParameters.splice(simIter->filteredParameters.end(), simIter->parameters, pIter++);// simIter->parameters.erase(pIter++)
            }
            else {
                ++pIter;
            }
        }
        SimInterface::printfPtr("[INFO] Reduced parameters: %u / %u\n", reducedSize, originalSize);

        // populate isLocal
        unsigned int globalP = simIter->dimP_free;
        simIter->isLocal.assign(simIter->maxID+1, false);
//        for(auto pIter = simIter->parameters.begin(); pIter!= simIter->parameters.end(); ++pIter) { // put switched parameters first
//            if(pIter->canBeSwitched) {
//                pIter->name = "P[" + std::to_string(globalP++) + "]";
//                std::cout << pIter->name << std::endl;
//            }
//        }
        for(auto pIter = simIter->parameters.begin(); pIter!= simIter->parameters.end(); ++pIter) { // global ones afterwards
            if(!pIter->canBeVaried) {
                bool isLocal = pIter->e->isLocal();
                simIter->isLocal[pIter->id] = isLocal;
                if(!isLocal) {
                    pIter->name = "P[" + std::to_string(globalP++) + "]";
                }
            }
        }
        //SimInterface::printfPtr("[INFO] Free parameters: %u\n", simIter->dimP_free);
        //return;

        cppExport(*simIter, outputDir);

        // free memory
        for(unsigned int i=0; i<simIter->vecNodes.size(); i++) { // all expressions
            if(simIter->vecNodes[i] && simIter->vecNodes[i]->e!=nullptr) {
                simIter->vecNodes[i]->e->delPtr();
                simIter->vecNodes[i]->e = nullptr;
            }
//            if(simIter->vecNodes[i]->rhs!=nullptr) {
//                simIter->vecNodes[i]->rhs->delPtr();
//                simIter->vecNodes[i]->rhs = nullptr;
//            }
        }
    }

    delete[] xmlContent;
}

uint64_t getModelHash(const std::string &xmlFile, const std::set<unsigned int> &setObservers, const std::set<unsigned int> &_setParameters) {
    // merge sets into vec for hash
    std::vector<unsigned int> vecParametersAndObservers(_setParameters.size()+setObservers.size());
    std::copy(_setParameters.begin(), _setParameters.end(), std::back_inserter(vecParametersAndObservers));
    std::copy(setObservers.begin(), setObservers.end(), std::back_inserter(vecParametersAndObservers));
    int idMemorySize = vecParametersAndObservers.size()*sizeof(unsigned int);

    // c style non-constant char array for rapidxml
    FILE *f = fopen(xmlFile.c_str(), "rb");
    if(f==nullptr) {
        SimInterface::printfPtr("File not found: %s\n", xmlFile.c_str());
        return 0;
    }
    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    fseek(f, 0, SEEK_SET);

    char *xmlContent = new char[fsize+idMemorySize];
    fread(xmlContent, fsize, 1, f);
    fclose(f);

    memcpy(&xmlContent[fsize], vecParametersAndObservers.data(), idMemorySize);
    uint64_t hash = fnFNV1a(xmlContent, xmlContent+idMemorySize+1);

    delete[] xmlContent;

    return hash;
}
