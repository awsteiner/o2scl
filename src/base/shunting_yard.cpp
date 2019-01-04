/*
  -------------------------------------------------------------------
  
  Original author: Jesse Brown
  Modifications: Brandon Amos, redpois0n
  Modifications for O2scl copyright (C) 2017-2019, Andrew W. Steiner
  
  This file is part of O2scl.
  
  O2scl is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  O2scl is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with O2scl. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/
/*
  Original cpp-expresion-parser license from Brandon Amos:

  The MIT License (MIT)
  
  Copyright (c) 2013 Brandon Amos <http://bamos.github.io>
  
  Permission is hereby granted, free of charge, to any person
  obtaining a copy of this software and associated documentation files
  (the "Software"), to deal in the Software without restriction,
  including without limitation the rights to use, copy, modify, merge,
  publish, distribute, sublicense, and/or sell copies of the Software,
  and to permit persons to whom the Software is furnished to do so,
  subject to the following conditions:

  The above copyright notice and this permission notice shall be
  included in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
  BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
  ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
  CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cmath>

#include <o2scl/shunting_yard.h>

using namespace o2scl;

std::map<std::string, int> calculator::buildOpPrecedence() {
  std::map<std::string, int> opp;

  // Create the operator precedence map based on C++ default
  // precedence order as described on cppreference website:
  // http://en.cppreference.com/w/c/language/operator_precedence
  opp["sin"]=2;
  opp["cos"]=2;
  opp["tan"]=2;
  opp["sqrt"]=2;
  opp["log"]=2;
  opp["exp"]=2;
  opp["abs"]=2;
  opp["log10"]=2;
  opp["asin"]=2;
  opp["acos"]=2;
  opp["atan"]=2;
  opp["sinh"]=2;
  opp["cosh"]=2;
  opp["tanh"]=2;
  opp["asinh"]=2;
  opp["acosh"]=2;
  opp["atanh"]=2;
  opp["^"]=2;
  opp["*"]=3;
  opp["/"]=3;
  opp["%"]=3;
  opp["+"]=4;
  opp["-"]=4;
  opp["<<"]=5;
  opp[">>"]=5;
  opp["<"]=6;
  opp["<="]=6;
  opp[">="]=6;
  opp[">"]=6;
  opp["=="]=7;
  opp["!="]=7;
  opp["&&"]=11;
  opp["||"]=12;
  opp["("]=16;

  return opp;
}

// Builds the opPrecedence map only once:
std::map<std::string, int> calculator::opPrecedence=
  calculator::buildOpPrecedence();

bool calculator::isvariablechar(const char c) {
  return (isalpha(c) || c == '_');
}

TokenQueue_t calculator::toRPN(const char* expr,
			       std::map<std::string, double>* vars,
			       bool debug,
			       std::map<std::string, int> opPrec) {

  TokenQueue_t rpnQueue;
  std::stack<std::string> operatorStack;
  bool lastTokenWasOp = true;

  // In one pass, ignore whitespace and parse the expression into RPN
  // using Dijkstra's Shunting-yard algorithm.
  while (*expr && isspace(*expr)) ++expr;

  while (*expr) {
    
    if (isdigit(*expr)) {

      // If the token is a number, add it to the output queue.
      char* nextChar = 0;
      double digit = strtod(expr,&nextChar);
      if (debug) {
	std::cout << "digit: " << digit << std::endl;
      }
      rpnQueue.push(new Token<double>(digit, NUM));
      expr = nextChar;
      lastTokenWasOp = false;

    } else if (isvariablechar(*expr)) {

      // If the function is a variable, resolve it and
      // add the parsed number to the output queue.
      std::stringstream ss;
      ss << *expr;
      ++expr;
      while (isvariablechar(*expr) || isdigit(*expr)) {
        ss << *expr;
        ++expr;
      }

      bool found = false;
      double val;

      std::string key = ss.str();

      if (key=="sin") {
	operatorStack.push("sin");
	lastTokenWasOp=true;
      } else if (key=="cos") {
	operatorStack.push("cos");
	lastTokenWasOp=true;
      } else if (key=="tan") {
	operatorStack.push("tan");
	lastTokenWasOp=true;
      } else if (key=="sqrt") {
	operatorStack.push("sqrt");
	lastTokenWasOp=true;
      } else if (key=="log") {
	operatorStack.push("log");
	lastTokenWasOp=true;
      } else if (key=="exp") {
	operatorStack.push("exp");
	lastTokenWasOp=true;
      } else if (key=="abs") {
	operatorStack.push("abs");
	lastTokenWasOp=true;
      } else if (key=="log10") {
	operatorStack.push("log10");
	lastTokenWasOp=true;
      } else if (key=="asin") {
	operatorStack.push("asin");
	lastTokenWasOp=true;
      } else if (key=="acos") {
	operatorStack.push("acos");
	lastTokenWasOp=true;
      } else if (key=="atan") {
	operatorStack.push("atan");
	lastTokenWasOp=true;
      } else if (key=="sinh") {
	operatorStack.push("sinh");
	lastTokenWasOp=true;
      } else if (key=="cosh") {
	operatorStack.push("cosh");
	lastTokenWasOp=true;
      } else if (key=="tanh") {
	operatorStack.push("tanh");
	lastTokenWasOp=true;
      } else if (key=="asinh") {
	operatorStack.push("asinh");
	lastTokenWasOp=true;
      } else if (key=="acosh") {
	operatorStack.push("acosh");
	lastTokenWasOp=true;
      } else if (key=="atanh") {
	operatorStack.push("atanh");
	lastTokenWasOp=true;
      } else {
	
	if (key == "true") {
	  found = true;
	  val = 1;
	} else if (key == "false") {
	  found = true;
	  val = 0;
	} else if (vars) {
	  std::map<std::string, double>::iterator it = vars->find(key);
	  if(it != vars->end()) {
	    found = true;
	    val = it->second;
	  }
	}
	
	if (found) {
	  // Save the number
	  if (debug) {
	    std::cout << "val: " << val << std::endl;
	  }
	  rpnQueue.push(new Token<double>(val, NUM));;
	} else {
	  // Save the variable name:
	  if (debug) {
	    std::cout << "key: " << key << std::endl;
	  }
	  rpnQueue.push(new Token<std::string>(key, VAR));
	}
	
	lastTokenWasOp = false;
	
      }

    } else {

      // Otherwise, the variable is an operator or parenthesis.
      switch (*expr) {
      case '(':
	operatorStack.push("(");
	++expr;
	break;
      case ')':
	while (operatorStack.top().compare("(")) {
	  rpnQueue.push(new Token<std::string>(operatorStack.top(),OP));
	  operatorStack.pop();
	}
	operatorStack.pop();
	++expr;
	break;

      default:
	{
	  // The token is an operator.
	  //
	  // Let p(o) denote the precedence of an operator o.
	  //
	  // If the token is an operator, o1, then
	  //   While there is an operator token, o2, at the top
	  //       and p(o1) <= p(o2), then
	  //     pop o2 off the stack onto the output queue.
	  //   Push o1 on the stack.
	  std::stringstream ss;
	  ss << *expr;
	  ++expr;
	  while (*expr && !isspace(*expr) && !isdigit(*expr)
		 && !isvariablechar(*expr) && *expr != '(' &&
		 *expr != ')') {
	    ss << *expr;
	    ++expr;
	  }
	  ss.clear();
	  std::string str;
	  ss >> str;
	  if (debug) {
	    std::cout << "str: " << str << std::endl;
	  }

	  if (lastTokenWasOp) {
	    // Convert unary operators to binary in the RPN.
	    if (!str.compare("-") || !str.compare("+")) {
	      rpnQueue.push(new Token<double>(0, NUM));
	    } else {
	      throw std::domain_error("Unrecognized unary operator: '" +
				      str + "'.");
	    }
	  }

	  while (!operatorStack.empty() &&
		 opPrec[str] >= opPrec[operatorStack.top()]) {
	    rpnQueue.push(new Token<std::string>(operatorStack.top(), OP));
	    operatorStack.pop();
	  }
	  operatorStack.push(str);
	  lastTokenWasOp = true;
	}
      }
    } while (*expr && isspace(*expr)) ++expr;

    // End of while (*expr)
  }
  
  while (!operatorStack.empty()) {
    rpnQueue.push(new Token<std::string>(operatorStack.top(), OP));
    operatorStack.pop();
  }

  return rpnQueue;
}

double calculator::calculate(const char* expr,
			     std::map<std::string, double>* vars,
			     bool debug) {

  // Convert to RPN with Dijkstra's Shunting-yard algorithm.
  TokenQueue_t rpn = toRPN(expr,vars,debug,opPrecedence);

  double ret = calculate(rpn);

  cleanRPN(rpn);

  return ret;
}

double calculator::calculate(TokenQueue_t rpn,
			     std::map<std::string, double>* vars) {

  // Evaluate the expression in RPN form.
  std::stack<double> evaluation;
  while (!rpn.empty()) {
    TokenBase* base = rpn.front();
    rpn.pop();

    // Operator:
    if (base->type == OP) {
      Token<std::string>* strTok = static_cast<Token<std::string>*>(base);
      std::string str = strTok->val;
      /*
	if (evaluation.size() < 2) {
        throw std::domain_error("Invalid equation.");
	}
      */
      double right = evaluation.top();
      evaluation.pop();
      if (!str.compare("sin")) {
	evaluation.push(sin(right));
      } else if (!str.compare("cos")) {
	evaluation.push(cos(right));
      } else if (!str.compare("tan")) {
	evaluation.push(tan(right));
      } else if (!str.compare("sqrt")) {
	evaluation.push(sqrt(right));
      } else if (!str.compare("log")) {
	evaluation.push(log(right));
      } else if (!str.compare("exp")) {
	evaluation.push(exp(right));
      } else if (!str.compare("abs")) {
	evaluation.push(std::abs(right));
      } else if (!str.compare("log10")) {
	evaluation.push(log10(right));
      } else if (!str.compare("asin")) {
	evaluation.push(asin(right));
      } else if (!str.compare("acos")) {
	evaluation.push(acos(right));
      } else if (!str.compare("atan")) {
	evaluation.push(atan(right));
      } else if (!str.compare("sinh")) {
	evaluation.push(sinh(right));
      } else if (!str.compare("cosh")) {
	evaluation.push(cosh(right));
      } else if (!str.compare("tanh")) {
	evaluation.push(tanh(right));
      } else if (!str.compare("asinh")) {
	evaluation.push(asinh(right));
      } else if (!str.compare("acosh")) {
	evaluation.push(acosh(right));
      } else if (!str.compare("atanh")) {
	evaluation.push(atanh(right));
      } else {
	double left  = evaluation.top();
	evaluation.pop();
	if (!str.compare("+")) {
	  evaluation.push(left + right);
	} else if (!str.compare("*")) {
	  evaluation.push(left * right);
	} else if (!str.compare("-")) {
	  evaluation.push(left - right);
	} else if (!str.compare("/")) {
	  evaluation.push(left / right);
	} else if (!str.compare("<<")) {
	  evaluation.push((int) left << (int) right);
	} else if (!str.compare("^")) {
	  evaluation.push(pow(left, right));
	} else if (!str.compare(">>")) {
	  evaluation.push((int) left >> (int) right);
	} else if (!str.compare("%")) {
	  evaluation.push((int) left % (int) right);
	} else if (!str.compare("<")) {
	  evaluation.push(left < right);
	} else if (!str.compare(">")) {
	  evaluation.push(left > right);
	} else if (!str.compare("<=")) {
	  evaluation.push(left <= right);
	} else if (!str.compare(">=")) {
	  evaluation.push(left >= right);
	} else if (!str.compare("==")) {
	  evaluation.push(left == right);
	} else if (!str.compare("!=")) {
	  evaluation.push(left != right);
	} else if (!str.compare("&&")) {
	  evaluation.push((int) left && (int) right);
	} else if (!str.compare("||")) {
	  evaluation.push((int) left || (int) right);
	} else {
	  throw std::domain_error("Unknown operator: '" + str + "'.");
	}
      }
    } else if (base->type == NUM) { // Number
      Token<double>* doubleTok = static_cast<Token<double>*>(base);
      evaluation.push(doubleTok->val);
    } else if (base->type == VAR) { // Variable
      if (!vars) {
        throw std::domain_error
	  ("Detected variable, but the variable map is null.");
      }
      
      Token<std::string>* strTok = static_cast<Token<std::string>*>(base);
      
      std::string key = strTok->val;
      std::map<std::string, double>::iterator it = vars->find(key);
      
      if (it == vars->end()) {
        throw std::domain_error("Unable to find the variable '" + key + "'.");
      }
      evaluation.push(it->second);
    } else {
      throw std::domain_error("Invalid token.");
    }
  }
  return evaluation.top();
}

void calculator::cleanRPN(TokenQueue_t& rpn) {
  while (rpn.size()) {
    delete rpn.front();
    rpn.pop();
  }
  return;
}

calculator::~calculator() {
  cleanRPN(this->RPN);
}

calculator::calculator(const char* expr,
		       std::map<std::string, double>* vars,
		       bool debug,
		       std::map<std::string, int> opPrec) {
  compile(expr,vars,debug,opPrec);
}

void calculator::compile(const char* expr,
			 std::map<std::string, double>* vars,
			 bool debug,
			 std::map<std::string, int> opPrec) {

  // Make sure it is empty:
  cleanRPN(this->RPN);

  this->RPN = calculator::toRPN(expr,vars,debug,opPrec);
}

double calculator::eval(std::map<std::string, double>* vars) {
  return calculate(this->RPN, vars);
}

std::string calculator::RPN_to_string() {
  std::stringstream ss;
  TokenQueue_t rpn = this->RPN;

  ss << "calculator { RPN: [ ";
  while (rpn.size()) {
    TokenBase* base = rpn.front();

    Token<double>* doubleTok = dynamic_cast<Token<double>*>(base);
    if(doubleTok) {
      ss << doubleTok->val;
    }

    Token<std::string>* strTok = dynamic_cast<Token<std::string>*>(base);
    if(strTok) {
      ss << "'" << strTok->val << "'";
    }

    rpn.pop();

    ss << (rpn.size() ? ", ":"");
  }
  ss << " ] }";
  return ss.str();
}
