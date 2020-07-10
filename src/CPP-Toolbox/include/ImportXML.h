/*
 * ImportXML.h
 *
 *
 */

#ifndef SRC_IMPORTXML_H_
#define SRC_IMPORTXML_H_

#include "rapidxml-1.13/rapidxml.hpp"
#include <cmath>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <functional>
#include <climits>
#include <string>

struct assignment {
    unsigned int oldid;
    unsigned int newid;
    bool useAsValue;
    assignment(unsigned int _oldid, unsigned int _newid, bool _useAsValue) : oldid(_oldid), newid(_newid), useAsValue(_useAsValue) {}
};

enum { XML_OBSERVER, XML_EVENT, XML_FORMULA, XML_VARIABLE, XML_PARAMETER, XML_SOLVER };

class Expression;

struct Node {
    unsigned int id = -1;
    Expression *e = nullptr;
//    Expression *rhs = nullptr;
    std::string name = "";
    std::string path = "";
    std::string unit = "";
    std::string equation = "";
    unsigned int formulaId = -1;
    double value = NAN;
    double scalefactor = NAN;
    bool persistent = false;
    bool canBeVaried = true;
    bool canBeSwitched = false;
    bool negativeValuesAllowed = true;
    bool oneTime = false;
    std::list<unsigned int> listRHS;
    std::map< std::string, unsigned int > mapEqRef;
    std::list<assignment> listAssignments;
    bool allAssignmentsConstant = true;
//unsigned int switchIndex = 0;

    // internal
    unsigned int order = -1;
    unsigned int switchIndex = 0;
    bool switchConstant = true;
    void print();

    bool operator<(const Node &other) const {
        return order < other.order;
    }
};

//struct variable {
//	unsigned int id = -1;
//	double value = NAN;
//	double scaleFactor = NAN;
//	Expression *rhs = nullptr;
//	Expression *init = nullptr;
//
//	bool persistable = false;
//	bool negativeValuesAllowed = true;
//
//	std::string name = "";
//	std::string path = "";
//	std::string unit = "";
//
//	void parse();
//	void appendDependencies();
//
//private:
//	unsigned int formulaId = -1;
//	std::list<unsigned int> listRHS;
//};

struct SolverOptions {
    std::string name = "";
    unsigned int relTolID = -1;
    unsigned int absTolID = -1;
    unsigned int hInitID = -1;
    unsigned int hMinID = -1;
    unsigned int hMaxID = -1;
    unsigned int maxStepsID = -1;
    unsigned int useJacID = -1;
    //	double relTol;
    //	double absTol;
    //	double hInit;
    //	double hMin;
    //	double hMax;
    //	long int maxSteps;
    //	bool useJac;
};

class Simulation {
public:
    unsigned int maxID = 0;
    unsigned int dimP_free = 0;
    unsigned int dimP = 0;
    unsigned int dimY = 0;
    std::string name;
    uint64_t hash;

    static Node Time; // explicit node for time (id=0)
    std::vector<Node*> vecNodes;

    std::list<Node> observers;
    std::list<Node> events;
    std::list<Node> formulas;
    std::list<Node> variables;
    std::list<Node> parameters;
    SolverOptions solver;
    std::set<double> outputTimepoints;

    // lists of filtered nodes, dont erase directly, because
    // some parameters such as solver settings may not have
    // an obvious dependency
    std::list<Node> filteredObservers;
    std::list<Node> filteredVariables;
    std::list<Node> filteredParameters;

    //    std::vector<Expression*> vecExpressions;
    unsigned int orderIndex = 0;
    std::vector<bool> isDependent;
    std::vector<bool> isLocal;

    Simulation(rapidxml::xml_node<> *n);
    void printStats();
    void resetDependencies() {
        isDependent.assign(maxID+1,false);
    }

private:
    void getTimepoints(rapidxml::xml_node<> *n, std::set<double> &setT);
    void getSolver(rapidxml::xml_node<> *n, SolverOptions &solver);
    void getList(rapidxml::xml_node<> *n, std::list<Node> &l);
};

void convertXMLtoCpp(const std::string &xmlFile, const std::string &_exportDir, const std::string &_modelName,
                     const std::set<unsigned int> &setObservers, const std::set<unsigned int> &setParameters, bool reduceVariables);
uint64_t getModelHash(const std::string &xmlFile, const std::set<unsigned int> &setObservers, const std::set<unsigned int> &_setParameters);
void cppExport(Simulation &sim, const std::string &exportDir);

#endif /* SRC_SIMIMPORTXML_H_ */
