#include "Expression.h"

#include <vector>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>

// debug functions
void printStack(std::deque<Expression*> &out, std::deque<Expression*> &ops) {
    for(std::deque<Expression*>::iterator iter = out.begin(); iter != out.end(); ++iter)
        std::cout << (*iter)->getName() << " ";
    std::cout << " | ";
    for(std::deque<Expression*>::iterator iter = ops.begin(); iter != ops.end(); ++iter)
        std::cout << (*iter)->getName() << " ";
    std::cout << std::endl;
}

// parse error -> clean-up memory
Expression* parseError(std::deque<Expression*> &out, std::deque<Expression*> &ops, const std::string &str) {
    // parse error
    if(!str.empty())
        std::cout << str << std::endl;
    for(std::deque<Expression*>::iterator iter = out.begin(); iter != out.end(); ++iter)
        (*iter)->delPtr();
    for(std::deque<Expression*>::iterator iter = ops.begin(); iter != ops.end(); ++iter)
        (*iter)->delPtr();
    return nullptr;
}

// parse error in build expression tree -> clean-up memory, but remove already established cross references first
Expression* parseError(std::deque<Expression*> &out, const std::string &str) {
    if(!str.empty())
        std::cout << str << std::endl;
    for(std::deque<Expression*>::iterator iter = out.begin(); iter != out.end(); ++iter)
        (*iter)->delPtr();
    return nullptr;
}

Expression* buildExpressionTree(std::deque<Expression*> &out) { // &out not const to be able to clean up memory in case of errors
    std::deque<Expression*> val;
    for(std::deque<Expression*>::const_iterator iter=out.begin(); iter!=out.end(); ++iter) {
        if( (*iter)->type<=OP_PARAMETER ) {
            val.push_back(*iter);
        }
        else {
            size_t n = (*iter)->getNoExpressions();
            if(n>val.size())
                return parseError(out, "buildExpressionTree error: expected more expressions left.");
            if(n>(*iter)->args.size()) {
                std::cout << (*iter)->type << "|" << (*iter)->args.size() << std::endl;
                return parseError(out, "ARG.");
            }

            for(int i=n-1; i>=0; i--) { // keep original order
                (*iter)->args[i] = val.back();
                val.pop_back();
            }
            val.push_back(*iter);
        }
        //		for(std::deque<expr*>::iterator iter = val.begin(); iter != val.end(); ++iter)
        //			std::cout << (*iter)->getName() << " ";
        //		std::cout << std::endl;
    }
    // val should be a single element
    if(val.size()!=1)
        return parseError(out, "buildExpressionTree error: expected only a single expression left.");

    return val.back();
}

// shunting yard algorithm, string not const due to in-place-conversion for comma separated floats
Expression* parseExpression(std::string &stringExpression ) {
    std::deque<Expression*> out;
    std::deque<Expression*> ops;

    const char *pos = stringExpression.c_str();
    Expression* lastExpr = nullptr; // last expression for unary minus

    while(*pos != '\0') {
        // ignore white spaces
        while ( isspace(*pos) )
            ++pos;

        // conversion, because MoBi allows comma separated floats in equations
        if(*pos==',') {
            std::replace( stringExpression.begin(), stringExpression.end(), ',', '.');
            parseError(out, ops, ""); // clean memory and restart parsing
            out.clear();
            ops.clear();
            pos = stringExpression.c_str();
            continue;
        }

        // number
        if (isdigit(*pos) || *pos == '.')  {
            char *end;
            double value = strtod(pos, &end);
            if(end!=NULL) {
                pos = end-1; // next char at end of loop
                out.push_back(lastExpr = Expression::fromConstant(value));
            }
            else {
                return parseError(out, ops, "parseExpression error: reading constant failed.");
            }
        }
        // function or variable name
        else if ( isalpha(*pos) ) {
            const char *start = pos;
            while(isalnum(*pos) || *pos=='_')
                ++pos;
            lastExpr = Expression::fromString( std::string(start, pos-start) );
            if(lastExpr->type<=OP_SWITCH) { // can be parameter, variable or named constant
                out.push_back(lastExpr); // variable name
            }
            else
                ops.push_back(lastExpr); // function name
            --pos; // next char at end of loop
        }
        // operator
        else if (*pos==';') { // || *pos==','
            while(!ops.empty() && ops.back()->type!=OP_LEFTBRACKET) {
                out.push_back(ops.back());
                ops.pop_back();
            }
            if(ops.empty())
                return parseError(out, ops, "parseExpression error: unexpected \";\" or \",\" .)");
            lastExpr = ops.back(); // point to the left bracket to be prepared for unary operators
        }
        else if (*pos=='+') {
            // distinguish between unary and binary -
            if(lastExpr ==nullptr  || lastExpr->type>=OP_SUM || lastExpr->type==OP_LEFTBRACKET) { // beginning or
                ops.push_back(lastExpr = new Expression(OP_UNARYPLUS));
            }
            else {
                while(!ops.empty() && 6 <= ops.back()->getPrecedence() ) {
                    out.push_back(ops.back());
                    ops.pop_back();
                }
                ops.push_back(lastExpr = new Expression(OP_SUM));
            }
        }
        else if (*pos=='-') {
            // distinguish between unary and binary -
            if(lastExpr == nullptr || lastExpr->type>=OP_SUM || lastExpr->type==OP_LEFTBRACKET) { // beginning or
                ops.push_back(lastExpr = new Expression(OP_UNARYMINUS));
            }
            else {
                while(!ops.empty() && 6 <= ops.back()->getPrecedence() ) {
                    out.push_back(ops.back());
                    ops.pop_back();
                }
                ops.push_back(lastExpr = new Expression(OP_DIFF));
            }
        }
        else if (*pos=='*') {
            while(!ops.empty() && 7 <= ops.back()->getPrecedence() ) {
                out.push_back(ops.back());
                ops.pop_back();
            }
            ops.push_back(lastExpr = new Expression(OP_PROD));
        }
        else if (*pos=='/') {
            while(!ops.empty() && 7 <= ops.back()->getPrecedence() ) {
                out.push_back(ops.back());
                ops.pop_back();
            }
            ops.push_back(lastExpr = new Expression(OP_DIV));
        }
        else if (*pos=='^') {
            while(!ops.empty() && 9 <= ops.back()->getPrecedence() ) {
                out.push_back(ops.back());
                ops.pop_back();
            }
            ops.push_back(lastExpr = new Expression(OP_POW));
        }
        else if (*pos=='(') {
            ops.push_back(lastExpr = new Expression(OP_LEFTBRACKET));
        }
        else if (*pos==')') {
            while(!ops.empty() && ops.back()->type!=OP_LEFTBRACKET ) {
                out.push_back(ops.back());
                ops.pop_back();
            }
            if(ops.empty()) {
                // mismatch error - no left bracket
                return parseError(out, ops, "parseExpression error: bracket mismatch.");
            }
            ops.back()->delPtr(); // free memory -> bracket not necessary on output stack
            ops.pop_back(); // remove bracket
            if(!ops.empty() && ops.back()->isFunction()) {
                out.push_back(ops.back());
                ops.pop_back();
            }
        }
        else if (*pos=='?') {
            while(!ops.empty() && 1 < ops.back()->getPrecedence() ) {
                out.push_back(ops.back());
                ops.pop_back();
            }
            ops.push_back(lastExpr = new Expression(OP_IFOPEN));
        }
        else if (*pos==':') {
            // pop till ?
            while(!ops.empty() && ops.back()->type!=OP_IFOPEN ) {
                out.push_back(ops.back());
                ops.pop_back();
            }
            if(ops.empty()) {
                // mismatch error - no ? operator
                return parseError(out, ops, "parseExpression error: \":\" without \"?\"");
            }
            ops.back()->type = OP_IF;
            ops.back()->args.assign(3,nullptr);
        }
        else if (*pos=='&') {
            while(!ops.empty() && 3 <= ops.back()->getPrecedence() ) {
                out.push_back(ops.back());
                ops.pop_back();
            }
            ops.push_back(lastExpr = new Expression(OP_AND));
        }
        else if (*pos=='|') {
            while(!ops.empty() && 2 <= ops.back()->getPrecedence() ) {
                out.push_back(ops.back());
                ops.pop_back();
            }
            ops.push_back(lastExpr = new Expression(OP_OR));
        }
        else if (*pos=='!') {
            while(!ops.empty() && 8 <= ops.back()->getPrecedence() ) {
                out.push_back(ops.back());
                ops.pop_back();
            }
            ops.push_back(lastExpr = new Expression(OP_UNARYNOT));
        }
        else if (*pos=='=') {
            while(!ops.empty() && 4 <= ops.back()->getPrecedence() ) {
                out.push_back(ops.back());
                ops.pop_back();
            }
            ops.push_back(lastExpr = new Expression(OP_EQUAL));
        }
        else if (*pos=='<') {
            if(*(pos+1)=='=') {
                lastExpr = new Expression(OP_LESSEQUAL);
                ++pos;
            }
            else if(*(pos+1)=='>') {
                lastExpr = new Expression(OP_NOTEQUAL);
                ++pos;
            }
            else {
                lastExpr = new Expression(OP_LESS);
            }
            while(!ops.empty() && 5 <= ops.back()->getPrecedence() ) {
                out.push_back(ops.back());
                ops.pop_back();
            }
            ops.push_back(lastExpr);
        }
        else if (*pos=='>') {
            if(*(pos+1)=='=') {
                lastExpr = new Expression(OP_GREATEREQUAL);
                ++pos;
            }
            else
                lastExpr = new Expression(OP_GREATER);
            while(!ops.empty() && 5 <= ops.back()->getPrecedence() ) {
                out.push_back(ops.back());
                ops.pop_back();
            }
            ops.push_back(lastExpr);
        }
        else {
            std::cout << *pos << std::endl;
            return parseError(out, ops, "parseExpression error: unknown token.");
        }
        ++pos; // next character

//        printStack(out,ops);
    }

    while(!ops.empty()) {
        if(ops.back()->type==OP_LEFTBRACKET || ops.back()->type==OP_RIGHTBRACKET) {
            // mismatch
            return parseError(out, ops, "parseExpression error: bracket mismatch.");
        }
        out.push_back(ops.back());
        ops.pop_back();
    }

    //	printStack(out,ops);

    return buildExpressionTree(out);
}

//int main() {
//	// tests
//	const std::vector<std::pair<std::string, double> > tests = {
//			{"1+2", 1+2}, // 0
//			{"   1    -    2", 1-2},
//			{"1 + 2 * 3", 1+2*3},
//			{"sin(1)", sin(1)},
//			{"cos(1 + tan(2))", cos(1+tan(2))},
//			{"--1", - -1},
//			{"1---2^2", -3},
//			{"2/-1*-4", 2/-1*-4 },
//			{"sqrt(-5.3+1e20/3)", sqrt(-5.3+1e20/3)},
//			{"5^2^5+2^2^(2*3)", 9769721},
//			{"max(1;3)*min(5;3)", 9}, // 10
//			{"1<2", 1<2},
//			{"1>=2", 1>=2},
//			{"!1", !1},
//			{"NEQ(5+5^2 ; 5+5^2)", 0},
//			{"1 ? 1*5 : 0+1", 5}, // 15
//			{"0 ? 0 ? 2 : 3 : 4", 0 ? 0 ? 2 : 3 : 4},
//			{"1 ? 0 ? 2 : 3 : 4", 1 ? 0 ? 2 : 3 : 4},
//			{"0 ? 1 ? 2 : 3 : 4", 0 ? 1 ? 2 : 3 : 4},
//			{"1 ? 1 ? 2 : 3 : 4", 1 ? 1 ? 2 : 3 : 4},
//			{"- (NOT 1 ? 0 : 5)", -5 }, // 20
//			{"- NOT ( 1 ? 0 : 5)", -1 }
//	};
//	//Expression* e = parseExpression("3 + 4 * 2 / ( 1 - 5 ) ^ 2 ^ 3");//sin(10)^2+cos(5)
//	//Expression* e = parseExpression("sin ( max ( 2, 3 ) / 3 * 3.1415 )");
//	//Expression* e = parseExpression("OralApplicationsEnabled ? P_int_para*Aeff*(min(Solubility/ MW ;DrugLiquid/(Liquid))-(NOT SinkCondition_para ? fu*DrugMucosa : 0)) : 0");
//	//Expression* e = parseExpression("pi*(r1*r1+r1*r2+r2*r2)/3*L");
////	std::string str("X+Y+Z<=10"); // double free
////	Expression* e = parseExpression(str);
////	if(e!=nullptr)
////		std::cout << e->toString() << std::endl;
////	e->delPtr();
//    std::string str("x+y");
//    Expression* e = parseExpression(str);
//    if(e!=nullptr) {
//        std::cout << e->toString() << std::endl;
//        std::cout << e->sd(0) << std::endl;
//    }
////    e->delPtr();
////	for(size_t i=0; i<tests.size(); i++) {
////		std::string str(tests[i].first);
////		Expression *e = parseExpression(str);
////		if(e!=nullptr) {
////				std::cout << e->toString() << std::endl;
////				double res = e->evaluate();
////				std::cout << res << " == " << tests[i].second << " : " << ((res == tests[i].second) ? "[OK]" : "[FAILED]" )<< std::endl;
////		}
////		else {
////			std::cout << "Parse error: " << i << std::endl;
////		}
////		e->delPtr();
////		std::cout << "*******************" << std::endl;
////	}
//	return 0;
//}
