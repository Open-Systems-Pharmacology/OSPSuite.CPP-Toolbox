/*
 * Expression.h
 *
 *
 */

#ifndef SRC_EXPRESSION_H_
#define SRC_EXPRESSION_H_

#include <string>
#include <deque>
#include <vector>
#include <map>
#include <cmath>
#include <climits>

#include "ImportXML.h"

enum Types { OP_CONSTANT, OP_PARAMETER, OP_VARIABLE, OP_SWITCH,
    OP_SUM, OP_DIFF, OP_PROD, OP_DIV, OP_POW, OP_AND, OP_OR, OP_UNARYNOT, OP_UNARYPLUS, OP_UNARYMINUS,
    OP_LEFTBRACKET, OP_RIGHTBRACKET, OP_IFOPEN, // helper ops that should never be on output stack
    OP_IF,
    OP_EQUAL, OP_NOTEQUAL, OP_LESSEQUAL, OP_GREATEREQUAL, OP_LESS, OP_GREATER, // bool operators
    OP_MIN, OP_MAX,
    OP_ACOS, OP_ASIN, OP_ATAN, OP_COS, OP_COSH, OP_EXP, OP_LOG, OP_LOG10, OP_SIN, OP_SINH, OP_SQRT, OP_TAN, OP_TANH };

// + pi, e, RND, SRND
// + LT, GT, NEQ, GEQ, LEQ, EQ as functions
// + AND, OR, NOT & |

struct Node;

class Expression {
public:
    Types type;
    unsigned int id = -1;
    unsigned int refs = 1;
    bool isSimplified = false;
    std::vector<Expression*> args;//[3] = { nullptr, nullptr, nullptr };

private:
    double value = std::nan(""); // protect value, because changes have to be written to name as well
    std::string name = "";
    const std::string *ptrName = &name;

    ~Expression() {
        for(unsigned int i=0; i<args.size(); i++)
            args[i]->delPtr();
    }

public:
    Expression(Types opType) : type(opType), args(getNoExpressions(),nullptr) {};
    Expression(Types opType, Expression* e) : args(1,nullptr) { type = opType; args[0] = e; }
    Expression(Types opType, Expression* el, Expression* er) : args(2,nullptr) { type = opType; args[0] = el; args[1] = er; }
    Expression(Types opType, Expression* el, Expression* ec, Expression* er) : args(3,nullptr) { type = opType; args[0] = el; args[1] = ec; args[2] = er; }
    Expression(const Expression &other) = delete;

    static Expression* fromString(const std::string &s);
    static Expression* fromString(const std::string *s);
    static Expression* fromConstant(double val) { Expression* e = new Expression(OP_CONSTANT); e->setValue(val); return e; }
    bool operator<(const Expression &other) const;
//    Expression& operator=(const Expression &other);

    Expression* ad(Expression *x, Expression *dx);
    Expression* sd(unsigned int componentId);
    Expression* simplify();
//    Expression* sd_verbose(unsigned int c);
    Expression* newPtr() { ++refs; return this; };
    void delPtr() { --refs; if(refs==0) delete this; };

    double evaluate(std::vector<Node*>* vecNodes = nullptr) const;
    inline bool isConstant() const { return type == OP_CONSTANT; }
    inline bool isParameter() const { return type == OP_PARAMETER; }
    inline bool isFunction() const { return type >= OP_MIN; }
    inline bool isTime() const { return id==0; }
    inline bool isZero() const { return isConstant() && value == 0.0; }//std::fpclassify(value) == FP_ZERO
    inline bool isOne() const { return isConstant() && value == 1; }
    inline bool isExplicit() const { return type == OP_EQUAL && args[0]!=nullptr && args[1]!=nullptr
                                                             && ((args[0]->isTime() && args[1]->isConstant())
                                                             || (args[0]->isConstant() && args[1]->isTime())); }
    bool isLocal() const;
    void getDependencies(std::deque<unsigned int> &dependencies) const;
    void getDependencies(std::vector<bool> &dependencies) const;
    void getDependenciesWithRHS(std::vector<bool> &dependencies) const;
    unsigned int getPrecedence() const;
    unsigned int getNoExpressions() const;
//    std::string getName() const;
    const std::string& getName() const;
    void setName(const std::string *n) { ptrName = n; };
    std::string toString(unsigned int parentPrecedence = -1) const;
    void setValue(double val);
    double getValue() { return value; }
    void replaceReferences(const std::map<std::string, unsigned int> &refs, const std::vector<Node*> &vecNodes);
};


Expression* parseExpression(std::string &stringExpression);
Expression* buildAST(const std::deque<Expression*> &out);

#endif /* SRC_SIMEXPRESSION_H_ */
