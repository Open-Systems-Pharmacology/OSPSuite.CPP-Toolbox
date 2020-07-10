/*
 * Expression.cpp
 *
 */

#include "Expression.h"
#include "ImportXML.h" // Node

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>

#ifndef M_PI
constexpr double M_PI = 3.14159265358979323846;
#endif
#ifndef M_E
constexpr double M_E = 2.7182818284590452354;
#endif

const std::vector<std::string> defaultNames = {
        "",         // OP_CONSTANT -- skip
        "",         // OP_PARAMETER -- skip
        "",         // OP_VARIABLE -- skip
        "",         // OP_SWITCH -- skip
        " + ",      // OP_SUM
        " - ",      // OP_DIFF
        " * ",      // OP_PROD
        " / ",      // OP_DIV
        " ^ ",      // OP_POW
        " && ",     // OP_AND
        " || ",     // OP_OR
        " ! ",      // OP_UNARYNOT
        " u+ ",     // OP_UNARYPLUS
        " u- ",     // OP_UNARYMINUS
        " ( ",      // OP_LEFTBRACKET
        " ) ",      // OP_RIGHTBRACKET
        " t? ",     // OP_IFOPEN
        " ? ",      // OP_IF
        " == ",     // OP_EQUAL
        " != ",     // OP_NOTEQUAL
        " <= ",     // OP_LESSEQUAL
        " >= ",     // OP_GREATEREQUAL
        " < ",      // OP_LESS
        " > ",      // OP_GREATER
        "min",    // OP_MIN
        "max",    // OP_MAX
        "acos",   // OP_ACOS
        "asin",   // OP_ASIN
        "atan",   // OP_ATAN
        "cos",    // OP_COS
        "cosh",   // OP_COSH
        "exp",    // OP_EXP
        "log",    // OP_LOG
        "log10",  // OP_LOG10
        "sin",    // OP_SIN
        "sinh",   // OP_SINH
        "sqrt",   // OP_SQRT
        "tan",    // OP_TAN
        "tanh" }; //OP_TANH

//Expression::Expression(const Expression &other) : type(other.type), id(other.id), value(other.value), name(other.name) {
//	for(unsigned int i=0; i<3; i++)
//		if(other.args[i]!=nullptr)
//			args[i] = other.args[i]->newPtr();
//}

bool Expression::operator<(const Expression &other) const {
    bool result = type < other.type && value < other.value && args.size() < other.args.size();
    for(unsigned int i=0; i<args.size(); i++)
        result = result && args[i] < other.args[i];
    return result;
}

//const Expression& Expression::operator=(const Expression &other) {
//    return other;
//}

bool Expression::isLocal() const {

    if(type==OP_VARIABLE || (id==0 && type==OP_PARAMETER)) //args[1]!=nullptr,(id!=UINT_MAX && getName()[0]=='y'),  || vecNodes[id]->canBeSwitched
        return true;

    bool isLocalVar = false;
    for(unsigned int i=0; i<args.size();i++)
        isLocalVar = isLocalVar || args[i]->isLocal();
    return isLocalVar;
}

void Expression::getDependencies(std::deque<unsigned int> &dependencies) const {
    if(type>=OP_PARAMETER && type<=OP_SWITCH && id!=UINT_MAX) {
        // TODO: add circular dependence protection?
//        if(dependencies[id])
//            return;
//        else
        dependencies.emplace_back(id);
    }
    else {
        for(unsigned int i=0; i<args.size();i++)
            args[i]->getDependencies(dependencies);
    }
}

void Expression::getDependencies(std::vector<bool> &dependencies) const {
    if(type>=OP_PARAMETER && type<=OP_SWITCH && id!=UINT_MAX) {//check
        if(dependencies[id])
            return;
        else {
            dependencies[id]=true;
            if(type==OP_VARIABLE)
                return;
        }
    }
    for(unsigned int i=0; i<args.size();i++)
        args[i]->getDependencies(dependencies);
}

void Expression::getDependenciesWithRHS(std::vector<bool> &dependencies) const {
    if(type!=OP_VARIABLE)
        std::cout << "Not a variable!" << std::endl;
    if(args.size()!=2)
        std::cout << "Args wrong size!" << std::endl;

    args[0]->getDependencies(dependencies);
    args[1]->getDependencies(dependencies);
}

unsigned int Expression::getPrecedence() const {
    switch(type) {
        case OP_IFOPEN:
        case OP_IF:
            return 1;
        case OP_OR:
            return 2;
        case OP_AND:
            return 3;
        case OP_EQUAL:
        case OP_NOTEQUAL:
            return 4;
        case OP_LESSEQUAL:
        case OP_GREATEREQUAL:
        case OP_LESS:
        case OP_GREATER:
            return 5;
        case OP_SUM:
        case OP_DIFF:
            return 6;
        case OP_PROD:
        case OP_DIV:
            return 7;
        case OP_UNARYNOT:
        case OP_UNARYPLUS:
        case OP_UNARYMINUS:
            return 8;
        case OP_POW:
            return 9;
        default:
            //		if(type!=OP_LEFTBRACKET)
            //			std::cout << "Error: getPrecedence called for " << type << std::endl;
            return 0;
    }
}

Expression* Expression::ad(Expression *x, Expression *dx) {
    switch(type) {
        case OP_CONSTANT:
        case OP_PARAMETER:
        case OP_VARIABLE:
        case OP_SWITCH:
        case OP_SUM:
            return new Expression(OP_SUM, &dx[0], &dx[1]);
        case OP_DIFF:
            return new Expression(OP_DIFF, &dx[0], &dx[1]);
        case OP_PROD:
            return new Expression(OP_SUM, new Expression(OP_PROD, &dx[0], &x[1]), new Expression(OP_PROD, &x[0], &dx[1]));
        case OP_DIV:
            return new Expression(OP_DIFF, new Expression(OP_PROD, &dx[0], &x[1]), new Expression(OP_PROD, &x[0], &dx[1]));
        case OP_POW:
            return new Expression(OP_DIV, new Expression(OP_SUM, new Expression(OP_PROD, &dx[0], &x[1]), new Expression(OP_PROD, &x[0], &dx[1])), new Expression(OP_PROD, &x[1],&x[1]));
        case OP_AND:
        case OP_OR:
        case OP_UNARYNOT:
            // should not happen
            return nullptr;
        case OP_UNARYPLUS:
            return new Expression(OP_UNARYPLUS,dx);
        case OP_UNARYMINUS:
            return new Expression(OP_UNARYMINUS,dx);
        case OP_LEFTBRACKET:
        case OP_RIGHTBRACKET:
        case OP_IFOPEN:
            // should not happen
            return nullptr;
        case OP_IF:
            return new Expression(OP_IF, &x[0], &dx[1], &dx[2]);
        case OP_EQUAL:
        case OP_NOTEQUAL:
        case OP_LESSEQUAL:
        case OP_GREATEREQUAL:
        case OP_LESS:
        case OP_GREATER:
        case OP_MIN:
            return new Expression(OP_MIN, &dx[0], &dx[1]);
        case OP_MAX:
            return new Expression(OP_MAX, &dx[0], &dx[1]);
        case OP_ACOS:
            return new Expression(OP_DIV, Expression::fromConstant(1), new Expression(OP_SQRT, new Expression(OP_DIFF, Expression::fromConstant(1),new Expression(OP_PROD, &x[0], &x[0]))));
        case OP_ASIN:
            return new Expression(OP_DIV, Expression::fromConstant(-1), new Expression(OP_SQRT, new Expression(OP_DIFF, Expression::fromConstant(1),new Expression(OP_PROD, &x[0], &x[0]))));
        case OP_ATAN:
            return new Expression(OP_DIV, Expression::fromConstant(-1), new Expression(OP_SUM, Expression::fromConstant(1),new Expression(OP_PROD, &x[0], &x[0])));
        case OP_COS:
            return new Expression(OP_UNARYMINUS, new Expression(OP_SIN, x));
        case OP_COSH:
            return new Expression(OP_COS, x);
        case OP_EXP:
            return new Expression(OP_EXP, x);
        case OP_LOG:
            return new Expression(OP_DIV, Expression::fromConstant(1), x);
        case OP_LOG10:
            return new Expression(OP_DIV, Expression::fromConstant(1), new Expression(OP_PROD, x, Expression::fromConstant(log(10))));
        case OP_SIN:
            return new Expression(OP_DIV, x);
        case OP_SINH:
            return new Expression(OP_COSH, x);
        case OP_SQRT:
            return new Expression(OP_DIV, Expression::fromConstant(1), new Expression(OP_PROD, Expression::fromConstant(2), new Expression(OP_SQRT, x)));
        case OP_TAN:
            return new Expression(OP_DIV, Expression::fromConstant(1), new Expression(OP_PROD, new Expression(OP_COS, x), new Expression(OP_COS, x)));
        case OP_TANH:
            return new Expression(OP_DIV, Expression::fromConstant(1), new Expression(OP_PROD, new Expression(OP_COS, x), new Expression(OP_COS, x)));
    }
    // should not happen
    return nullptr;
}

Expression* Expression::sd(unsigned int c) {
//    std::cout << type << std::endl;
    switch(type) {
        case OP_CONSTANT:
            return Expression::fromConstant(0);
        case OP_PARAMETER:
            if(id==c)
                return Expression::fromConstant(1);
            else if(args.size())
                return args[0]->sd(c);
            else
                return Expression::fromConstant(0);
        case OP_VARIABLE:
            if(id==c)
                return Expression::fromConstant(1);
            else
                return Expression::fromConstant(0);
        case OP_SWITCH:
        {
            Expression *e = new Expression(OP_SWITCH);
            e->args.assign(args.size(), nullptr);
            for(unsigned int i=0; i<args.size(); i++)
                e->args[i] = args[i]->sd(c);
            return e;
//            return this->newPtr();
        }
        case OP_SUM:
            return new Expression(OP_SUM, args[0]->sd(c), args[1]->sd(c));
        case OP_DIFF:
            return new Expression(OP_DIFF, args[0]->sd(c), args[1]->sd(c));
        case OP_PROD:
            return new Expression(OP_SUM, new Expression(OP_PROD, args[0]->sd(c), args[1]->newPtr()), new Expression(OP_PROD, args[0]->newPtr(), args[1]->sd(c)));
        case OP_DIV:
            return new Expression(OP_DIV,
                                  new Expression(OP_DIFF, new Expression(OP_PROD, args[0]->sd(c), args[1]->newPtr()), new Expression(OP_PROD, args[0]->newPtr(), args[1]->sd(c))),
                                  new Expression(OP_PROD, args[1]->newPtr(),args[1]->newPtr()));
        case OP_POW:
            //d/dx(f(x)^(g(x))) = f(x)^(g(x) - 1) (g(x) f'(x) + f(x) log(f(x)) g'(x))
            return new Expression(OP_PROD,
                                  new Expression(OP_POW, args[0]->newPtr(), new Expression(OP_SUM, args[1]->newPtr(), fromConstant(-1))),
                                  new Expression(OP_SUM,
                                                 new Expression(OP_PROD, args[0]->sd(c), args[1]->newPtr()),
                                                 new Expression(OP_PROD,
                                                                new Expression(OP_PROD, args[0]->newPtr(), args[1]->sd(c)),
                                                                new Expression(OP_LOG, args[0]->newPtr()))));
        case OP_AND:
        case OP_OR:
        case OP_UNARYNOT:
            // should not happen
            return nullptr;
        case OP_UNARYPLUS:
            return new Expression(OP_UNARYPLUS,args[0]->sd(c));
        case OP_UNARYMINUS:
            return new Expression(OP_UNARYMINUS,args[0]->sd(c));
        case OP_LEFTBRACKET:
        case OP_RIGHTBRACKET:
        case OP_IFOPEN:
            // should not happen
            return nullptr;
        case OP_IF:
            return new Expression(OP_IF, args[0]->newPtr(), args[1]->sd(c), args[2]->sd(c));
        case OP_EQUAL:
        case OP_NOTEQUAL:
        case OP_LESSEQUAL:
        case OP_GREATEREQUAL:
        case OP_LESS:
        case OP_GREATER:
            // should not happen
            std::cout << "sd called for boolean operation! " << id << "|" << toString() << std::endl;
            return nullptr;
        case OP_MIN:
            return new Expression(OP_IF, new Expression(OP_LESSEQUAL, args[0]->newPtr(), args[1]->newPtr()), args[0]->sd(c), args[1]->sd(c));
        case OP_MAX:
            return new Expression(OP_IF, new Expression(OP_GREATEREQUAL, args[0]->newPtr(), args[1]->newPtr()), args[0]->sd(c), args[1]->sd(c));
        case OP_ACOS:
            return new Expression(OP_PROD, new Expression(OP_DIV, Expression::fromConstant(1), new Expression(OP_SQRT, new Expression(OP_DIFF, Expression::fromConstant(1),new Expression(OP_PROD, args[0]->newPtr(), args[0]->newPtr())))), args[0]->sd(c));
        case OP_ASIN:
            return new Expression(OP_PROD, new Expression(OP_DIV, Expression::fromConstant(-1), new Expression(OP_SQRT, new Expression(OP_DIFF, Expression::fromConstant(1),new Expression(OP_PROD, args[0]->newPtr(), args[0]->newPtr())))), args[0]->sd(c));
        case OP_ATAN:
            return new Expression(OP_PROD, new Expression(OP_DIV, Expression::fromConstant(-1), new Expression(OP_SUM, Expression::fromConstant(1),new Expression(OP_PROD, args[0]->newPtr(), args[0]->newPtr()))), args[0]->sd(c));
        case OP_COS:
            return new Expression(OP_PROD, new Expression(OP_UNARYMINUS, new Expression(OP_SIN, args[0]->newPtr())), args[0]->sd(c));
        case OP_COSH:
            return new Expression(OP_PROD, new Expression(OP_COS, args[0]->newPtr()), args[0]->sd(c));
        case OP_EXP:
            return new Expression(OP_PROD, new Expression(OP_EXP, args[0]->newPtr()), args[0]->sd(c));
        case OP_LOG:
            return new Expression(OP_PROD, new Expression(OP_DIV, Expression::fromConstant(1), args[0]->newPtr()), args[0]->sd(c));
        case OP_LOG10:
            return new Expression(OP_PROD, new Expression(OP_DIV, Expression::fromConstant(1), new Expression(OP_PROD, args[0]->newPtr(), Expression::fromConstant(log(10)))), args[0]->sd(c));
        case OP_SIN:
            return new Expression(OP_PROD, new Expression(OP_DIV, args[0]->newPtr()), args[0]->sd(c));
        case OP_SINH:
            return new Expression(OP_PROD, new Expression(OP_COSH, args[0]->newPtr()), args[0]->sd(c));
        case OP_SQRT:
            return new Expression(OP_PROD, new Expression(OP_DIV, Expression::fromConstant(1), new Expression(OP_PROD, Expression::fromConstant(2), new Expression(OP_SQRT, args[0]))), args[0]->sd(c));
        case OP_TAN:
            return new Expression(OP_PROD, new Expression(OP_DIV, Expression::fromConstant(1), new Expression(OP_PROD, new Expression(OP_COS, args[0]->newPtr()), new Expression(OP_COS, args[0]->newPtr()))), args[0]->sd(c));
        case OP_TANH:
            return new Expression(OP_PROD, new Expression(OP_DIV, Expression::fromConstant(1), new Expression(OP_PROD, new Expression(OP_COS, args[0]->newPtr()), new Expression(OP_COS, args[0]->newPtr()))), args[0]->sd(c));
    }
    // should not happen
    return nullptr;
}

Expression* Expression::simplify() {
    Expression *e = this;
    if(isSimplified)
        return e;
//    isSimplified = true;

    bool isArgsConstant = args.size();//true; // variable can store ref in arg[0], cant be evaluated below without vecNodes
    for(unsigned int i=0; i<args.size(); i++) {
        args[i] = args[i]->simplify();
        isArgsConstant = isArgsConstant && args[i]->isConstant();
    }

    if(isArgsConstant) {
        e = Expression::fromConstant(evaluate());
        delPtr();
    }
    else
    {
        switch(type) {
            case OP_CONSTANT:
                break;
            case OP_PARAMETER:
                if(args[0]->isConstant()) {
                    e = args[0]->newPtr();
                    delPtr();
                }
                else if(args[0]->type==OP_PARAMETER || args[0]->type==OP_VARIABLE) {
                    e = args[0]->newPtr();
                    delPtr();
                }
                break;
            case OP_VARIABLE:
                break;
            case OP_SWITCH:
                break;
            case OP_SUM:
                if(args[0]->isZero()) {
                    e = args[1]->newPtr();
                    delPtr();
                }
                else if(args[1]->isZero()) {
                    e = args[0]->newPtr();
                    delPtr();
                }
                else if(args[1]->type==OP_UNARYMINUS) {
                    e = new Expression(OP_DIFF, args[0]->newPtr(), args[1]->args[0]->newPtr());
                    e = e->simplify();
                    delPtr();
                }
                else if(args[1]->isConstant()) { // commute to ease further simplifications
                    Expression *temp = args[1];
                    args[1] = args[0];
                    args[0] = temp;
                }
                break;
            case OP_DIFF:
                if(args[0]->isZero()) {
                    e=new Expression(OP_UNARYMINUS, args[1]->newPtr());
                    e = e->simplify();
                    delPtr();
                }
                else if(args[1]->isZero()) {
                    e=args[0]->newPtr();
                    delPtr();
                }
                break;
            case OP_PROD:
                if(args[0]->isZero() || args[1]->isZero()) {
                    e = Expression::fromConstant(0);
                    delPtr();
                }
                else if(args[0]->isOne()) {
                    e = args[1]->newPtr();
                    delPtr();
                }
                else if(args[1]->isOne()) {
                    e = args[0]->newPtr();
                    delPtr();
                }
                else if(args[1]->isConstant()) { // commute to ease further simplifications
                    Expression *temp = args[1];
                    args[1] = args[0];
                    args[0] = temp;
                }
                break;
            case OP_DIV:
                if(args[1]->isOne()) {
                    e = args[0]->newPtr();
                    delPtr();
                }
                else if(args[0]->isZero()) {
                    e = Expression::fromConstant(0);
                    delPtr();
                }
                else if(args[1]->isConstant() && !args[1]->isZero()) { // replace by product
                    e = new Expression(OP_PROD, fromConstant(1.0/args[1]->value), args[0]->newPtr());
                    delPtr();
                }
                break;
            case OP_POW:
                if(args[0]->isZero()) {
                    e = Expression::fromConstant(0);
                    delPtr();
                }
                else if(args[1]->isZero()) {
                    e = Expression::fromConstant(1);
                    delPtr();
                }
                else if(args[1]->isOne()) {
                    e = args[0]->newPtr();
                    delPtr();
                }
                break;
            case OP_AND:
            case OP_OR:
            case OP_UNARYNOT:
                break;
            case OP_UNARYPLUS:
                e = args[0]->newPtr();
                delPtr();
                break;
            case OP_UNARYMINUS:
                if (args[0]->type==OP_UNARYMINUS) {
                    e = args[0]->args[0]->newPtr();
                    delPtr();
                }
                else if(args[0]->isConstant()) {
                    e = Expression::fromConstant(-args[0]->value);
                    delPtr();
                }
                break;
            case OP_LEFTBRACKET:
            case OP_RIGHTBRACKET:
            case OP_IFOPEN:
                break;
            case OP_IF:
                if(args[0]->isConstant()) {
                    e = args[0]->value ? args[1]->newPtr() : args[2]->newPtr();
                    delPtr();
                }
                else if(args[1]->isConstant() && args[2]->isConstant() && args[1]->value == args[2]->value) {
                    e = args[1]->newPtr();
                    delPtr();
                }
                break;
            case OP_EQUAL:
            case OP_NOTEQUAL:
            case OP_LESSEQUAL:
            case OP_GREATEREQUAL:
            case OP_LESS:
            case OP_GREATER:
            case OP_MIN:
            case OP_MAX:
            case OP_ACOS:
            case OP_ASIN:
            case OP_ATAN:
            case OP_COS:
            case OP_COSH:
            case OP_EXP:
            case OP_LOG:
            case OP_LOG10:
            case OP_SIN:
            case OP_SINH:
            case OP_SQRT:
            case OP_TAN:
            case OP_TANH:
                break;
        }
    }

    return e;
}

unsigned int Expression::getNoExpressions() const {
    switch(type) {
        case OP_CONSTANT:
            return 0;
        case OP_PARAMETER:  // may contain a reference, but is zero during buildExpressionTree
            return 1;
        case OP_VARIABLE:
            return 2;
        case OP_SWITCH: // reserve for old_id only first, then push new ids
            return 1;
        case OP_SUM:
        case OP_DIFF:
        case OP_PROD:
        case OP_DIV:
        case OP_POW:
        case OP_AND:
        case OP_OR:
            return 2;
        case OP_UNARYNOT:
        case OP_UNARYPLUS:
        case OP_UNARYMINUS:
            return 1;
        case OP_LEFTBRACKET:
        case OP_RIGHTBRACKET:
        case OP_IFOPEN:
            // should not happen
            return 0;
        case OP_IF:
            return 3;
        case OP_EQUAL:
        case OP_NOTEQUAL:
        case OP_LESSEQUAL:
        case OP_GREATEREQUAL:
        case OP_LESS:
        case OP_GREATER:
            return 2;
        case OP_MIN:
        case OP_MAX:
            return 2;
        case OP_ACOS:
        case OP_ASIN:
        case OP_ATAN:
        case OP_COS:
        case OP_COSH:
        case OP_EXP:
        case OP_LOG:
        case OP_LOG10:
        case OP_SIN:
        case OP_SINH:
        case OP_SQRT:
        case OP_TAN:
        case OP_TANH:
            return 1;
    }
    // should not happen
    return 0;
}

double Expression::evaluate(std::vector<Node*>* vecNodes) const {
    switch(type) {
        case OP_CONSTANT:
            return value;
        case OP_PARAMETER:
        case OP_VARIABLE:
        case OP_SWITCH:
            if(args.size())
                return args[0]->evaluate(vecNodes);
            else if(vecNodes==nullptr)
                return NAN;
            else {
                //          std::cout << "id: " << id << "   " << name << "   " << value << std::endl;
                if(id>=vecNodes->size()) {
                    std::cout << "wrong size: " << id << "   " << getName() << std::endl;
                    return NAN;
                }
                if(std::isnan((*vecNodes)[id]->value))//(*vecNodes)[id]->canBeVaried
                    if((*vecNodes)[id]->formulaId!=UINT_MAX && (*vecNodes)[(*vecNodes)[id]->formulaId]->e!=nullptr)
                        return (*vecNodes)[(*vecNodes)[id]->formulaId]->e->evaluate(vecNodes);
                //          std::cout << "Returns: " << (*vecNodes)[id]->value << std::endl;
                return (*vecNodes)[id]->value;
            }
        case OP_SUM:
            return args[0]->evaluate(vecNodes)+args[1]->evaluate(vecNodes);
        case OP_DIFF:
            return args[0]->evaluate(vecNodes)-args[1]->evaluate(vecNodes);
        case OP_PROD:
            return args[0]->evaluate(vecNodes)*args[1]->evaluate(vecNodes);
        case OP_DIV:
            return args[0]->evaluate(vecNodes)/args[1]->evaluate(vecNodes);
        case OP_POW:
            return pow(args[0]->evaluate(vecNodes),args[1]->evaluate(vecNodes));
        case OP_AND:
            return args[0]->evaluate(vecNodes)&&args[1]->evaluate(vecNodes);
        case OP_OR:
            return args[0]->evaluate(vecNodes)||args[1]->evaluate(vecNodes);
        case OP_UNARYNOT:
            return !args[0]->evaluate(vecNodes);
        case OP_UNARYPLUS:
            return args[0]->evaluate(vecNodes);
        case OP_UNARYMINUS:
            return -args[0]->evaluate(vecNodes);
        case OP_LEFTBRACKET:
        case OP_RIGHTBRACKET:
        case OP_IFOPEN:
            // should not happen
            return NAN;
        case OP_IF:
            return args[0]->evaluate(vecNodes) ? args[1]->evaluate(vecNodes) : args[2]->evaluate(vecNodes);
        case OP_EQUAL:
            return args[0]->evaluate(vecNodes)==args[1]->evaluate(vecNodes);
        case OP_NOTEQUAL:
            return args[0]->evaluate(vecNodes)!=args[1]->evaluate(vecNodes);
        case OP_LESSEQUAL:
            return args[0]->evaluate(vecNodes)<=args[1]->evaluate(vecNodes);
        case OP_GREATEREQUAL:
            return args[0]->evaluate(vecNodes)>=args[1]->evaluate(vecNodes);
        case OP_LESS:
            return args[0]->evaluate(vecNodes)<args[1]->evaluate(vecNodes);
        case OP_GREATER:
            return args[0]->evaluate(vecNodes)>args[1]->evaluate(vecNodes);
        case OP_MIN:
            return std::min(args[0]->evaluate(vecNodes),args[1]->evaluate(vecNodes));
        case OP_MAX:
            return std::max(args[0]->evaluate(vecNodes),args[1]->evaluate(vecNodes));
        case OP_ACOS:
            return acos(args[0]->evaluate(vecNodes));
        case OP_ASIN:
            return asin(args[0]->evaluate(vecNodes));
        case OP_ATAN:
            return atan(args[0]->evaluate(vecNodes));
        case OP_COS:
            return cos(args[0]->evaluate(vecNodes));
        case OP_COSH:
            return cosh(args[0]->evaluate(vecNodes));
        case OP_EXP:
            return exp(args[0]->evaluate(vecNodes));
        case OP_LOG:
            return log(args[0]->evaluate(vecNodes));
        case OP_LOG10:
            return log10(args[0]->evaluate(vecNodes));
        case OP_SIN:
            return sin(args[0]->evaluate(vecNodes));
        case OP_SINH:
            return sinh(args[0]->evaluate(vecNodes));
        case OP_SQRT:
            return sqrt(args[0]->evaluate(vecNodes));
        case OP_TAN:
            return tan(args[0]->evaluate(vecNodes));
        case OP_TANH:
            return tanh(args[0]->evaluate(vecNodes));
    }
    // should not happen
    return NAN;
}

Expression* Expression::fromString(const std::string &s) {
    // upper and lower case is not consistent for function names
    std::string s_lower(s);
    std::transform(s_lower.begin(), s_lower.end(), s_lower.begin(), ::tolower);

    auto iter = std::find(defaultNames.begin(), defaultNames.end(), s_lower);
    if(iter!=defaultNames.end())
        return new Expression(static_cast<Types>(iter-defaultNames.begin()));
    else {
        if(s_lower=="pi")
            return fromConstant(M_PI);
        else if(s_lower=="e")
            return fromConstant(M_E);
        else if(s_lower=="ln")
            return new Expression(OP_LOG);
        else if(s_lower=="and") //"AND"
            return new Expression(OP_AND);
        else if(s_lower=="not") //"NOT"
            return new Expression(OP_UNARYNOT);
        else if(s_lower=="or") //"OR"
            return new Expression(OP_OR);
        else if(s_lower=="eq") //"EQ"
            return new Expression(OP_EQUAL);
        else if(s_lower=="neq") //"NEQ"
            return new Expression(OP_NOTEQUAL);
        else if(s_lower=="leq") //"LEQ"
            return new Expression(OP_LESSEQUAL);
        else if(s_lower=="geq") //"GEQ"
            return new Expression(OP_GREATEREQUAL);
        else if(s_lower=="lt") // "LT"
            return new Expression(OP_LESS);
        else if(s_lower=="gt") //"GT"
            return new Expression(OP_GREATER);
    }
    Expression* e = new Expression(OP_PARAMETER);
    e->name = s; // use unmodified name
    return e;
}

Expression* Expression::fromString(const std::string *s) {
    Expression* e = new Expression(OP_PARAMETER);
    e->ptrName = s;
    return e;
}

//Expression* Expression::fromID(const unsigned int id) {
//	return nullptr;
//}

std::string Expression::toString(unsigned int parentPrecedence) const {
    unsigned int p = getPrecedence();
    bool useBraces = true;//type==OP_DIFF || type==OP_DIV || ( parentPrecedence!=UINT_MAX && parentPrecedence > p ? true : false );
    switch(type) {
        case OP_CONSTANT:
        case OP_PARAMETER:
        case OP_VARIABLE:
        case OP_SWITCH:
            return getName();
        case OP_SUM:
        case OP_DIFF:
        case OP_PROD:
        case OP_DIV:
            if(useBraces)
                return "(" + args[0]->toString(p) + defaultNames[type] + args[1]->toString(p) + ")";
            else
                return args[0]->toString(p) + defaultNames[type] + args[1]->toString(p);
        case OP_POW:
            return "pow(" + args[0]->toString(p) + ", " + args[1]->toString(p) + ")";
        case OP_AND:
        case OP_OR:
            return "(" + args[0]->toString(p) + defaultNames[type] + args[1]->toString(p) + ")";
        case OP_UNARYNOT:
            return " ! " + args[0]->toString(parentPrecedence);
        case OP_UNARYPLUS:
            return " + " + args[0]->toString(parentPrecedence);
        case OP_UNARYMINUS:
            return " - " + args[0]->toString(parentPrecedence);
        case OP_LEFTBRACKET:
        case OP_RIGHTBRACKET:
        case OP_IFOPEN:
            // should not happen
            return "err";
        case OP_IF:
            return "(" + args[0]->toString(p) + " ? " + args[1]->toString(p) + " : " + args[2]->toString(p) + ")";
            //	case OP_EQUAL:
            //		return "(" + args[0]->toString(p) + " == " + args[1]->toString(p) + ")";
            //	case OP_NOTEQUAL:
            //		return "(" + args[0]->toString(p) + " != " + args[1]->toString(p) + ")";
            //	case OP_LESSEQUAL:
            //		return "(" + args[0]->toString(p) + " <= " + args[1]->toString(p) + ")";
            //	case OP_GREATEREQUAL:
            //		return "(" + args[0]->toString(p) + " >= " + args[1]->toString(p) + ")";
            //	case OP_LESS:
            //		return "(" + args[0]->toString(p) + " < " + args[1]->toString(p) + ")";
            //	case OP_GREATER:
            //		return "(" + args[0]->toString(p) + " > " + args[1]->toString(p) + ")";
        case OP_EQUAL:
        case OP_NOTEQUAL:
        case OP_LESSEQUAL:
        case OP_GREATEREQUAL:
        case OP_LESS:
        case OP_GREATER:
            return "(" + args[0]->toString(p) + defaultNames[type] + args[1]->toString(p) + ")";
        case OP_MIN:
        case OP_MAX:
            return defaultNames[type] + "(" + args[0]->toString(p) + ", " + args[1]->toString(p) + ")";
        case OP_ACOS:
        case OP_ASIN:
        case OP_ATAN:
        case OP_COS:
        case OP_COSH:
        case OP_EXP:
        case OP_LOG:
        case OP_LOG10:
        case OP_SIN:
        case OP_SINH:
        case OP_SQRT:
        case OP_TAN:
        case OP_TANH:
            return defaultNames[type] + "(" + args[0]->toString(p) + ")";
    }
    // should not happen
    return std::string();
}

void Expression::setValue(double val) {
    value = val;
    // stringstream outputs nan instead of NAN and inf instead of Inf
    if(std::isnan(val))
        name = "NAN";
    else if(std::isinf(val))
        if(val>0)
            name = "INFINITY";
        else
            name = "-INFINITY";
    else {
        // update name
        std::stringstream sstr;
        sstr << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
        name = sstr.str();
    }
}

//std::string Expression::getName() const {
//    if(name.empty())
//        return defaultNames[type];
//    else
//        return name;
//}

const std::string& Expression::getName() const {
    return *ptrName;
//    if(ptrName)
//        return *ptrName;
//    else if(name.empty())
//        return defaultNames[type];
//    else
//        return name;
}

void Expression::replaceReferences(const std::map<std::string, unsigned int> &refIds, const std::vector<Node*> &vecNodes) {
    if(type>=OP_PARAMETER && type<=OP_SWITCH) {
        auto iter = refIds.find(name);
        if(iter==refIds.end())
            std::cout << "Error equation ref not found." << std::endl;
        else {
            id = iter->second;
            if(id>=vecNodes.size())// && vecNodes[iter->second]->e!=nullptr
                std::cout << "Error with id " << id << std::endl;
            else {
//                Expression *e = vecNodes[iter->second]->e->newPtr();
//                delPtr();
//                return e;
                if(vecNodes[id]->name.empty())
                    std::cout << "Empty name! id: " << id << std::endl;
                ptrName = &vecNodes[id]->name;
                if(vecNodes[id]->e->type==OP_VARIABLE) {
                    type = OP_VARIABLE;
                    isSimplified = true;
                    args[0] = vecNodes[id]->e->args[0]->newPtr();
                    args.push_back(vecNodes[id]->e->args[1]->newPtr());
//                    if(args[0]) {
//                    args[0]->delPtr();
//                    args[0] = nullptr;
//                    }
                }
                else
                    args[0] = vecNodes[id]->e->newPtr();
            }
        }
    }
    else {
        for(unsigned int i=0; i<args.size(); i++)
            args[i]->replaceReferences(refIds, vecNodes);
    }
//    return this;
}
