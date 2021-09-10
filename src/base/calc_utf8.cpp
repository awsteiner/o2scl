/*
  -------------------------------------------------------------------
  
  Original author: Jesse Brown
  Modifications: Brandon Amos, redpois0n
  Modifications for O2scl copyright (C) 2017-2021, Andrew W. Steiner
  
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

#include <o2scl/calc_utf8.h>
#include <o2scl/err_hnd.h>
#include <o2scl/string_conv.h>

using namespace o2scl;

typedef std::queue<token_base*> token_queue_t;

calc_utf8::calc_utf8(const std::u32string &expr,
                     const std::map<std::u32string, double> *vars,
                     bool debug,
                     std::map<std::string, int> op_prec) {
  compile(expr,vars,debug,op_prec);
}

calc_utf8::~calc_utf8() {
  cleanRPN(this->RPN);
}

token_queue_t calc_utf8::get_RPN() {
  return this->RPN;
}

std::vector<std::u32string> calc_utf8::get_var_list() {
  std::vector<std::u32string> list;
  token_queue_t rpn=get_RPN();
  while (rpn.size()) {
    token_base *tb=rpn.front();
    if (tb->type==2) {
      token32<std::u32string>* str=dynamic_cast<token32<std::u32string>*>(tb);
      if (str!=0) {
	list.push_back(str->val);
      } else {
	O2SCL_ERR2("Token was of type var but could not convert to ",
		   "string",o2scl::exc_einval);
      }
    }
    rpn.pop();
  }
  return list;
}  

std::map<std::string, int> calc_utf8::build_op_precedence() {
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
  opp["floor"]=2;
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

// Builds the op_precedence map only once:
std::map<std::string, int> calc_utf8::op_precedence=
  calc_utf8::build_op_precedence();

bool calc_utf8::is_variable_char(const char32_t c) {
  // Presume it's a variable if it's not an operator or digit
  if (c=='^' || c=='*' || c=='/' || c=='%' || c=='+' || c=='-' ||
      c=='<' || c=='=' || c=='>' || c=='!' || c=='&' || c=='|' ||
      c=='(' || c==')' || isdigit(c) || c=='.' || c==',') {
    return false;
  }
  return true;
}

int calc_utf8::toRPN_nothrow(const std::u32string &expr,
                             const std::map<std::u32string, double> *vars,
                             bool debug,
                             std::map<std::string, int> op_prec,
                             token_queue_t &rpn_queue2) {

  token_queue_t rpn_queue;
  std::stack<std::string> operator_stack;
  bool lastTokenWasOp = true;

  // In one pass, ignore whitespace and parse the expression into RPN
  // using Dijkstra's Shunting-yard algorithm.
  
  size_t i=0;
  while (i<expr.length()) {
    
    if (isdigit(expr[i])) {
      
      std::u32string part;
      part=part+expr[i];

      bool exponent=false;
      // The idea here is that this will handle both commas or
      // periods, depending on locale, but this hasn't been tested
      while (i+1<expr.length() &&
             (expr[i+1]=='.' || expr[i+1]==',' ||
              isdigit(expr[i+1]) || expr[i+1]=='e' ||
              expr[i+1]=='E' || (exponent && (expr[i+1]=='+' ||
                                              expr[i+1]=='-')))) {
        if (expr[i+1]=='e' || expr[i+1]=='E') exponent=true;
        part=part+expr[i+1];
        i++;
      }

      double value;
      int iret=s32tod_nothrow(part,value);
      if (iret!=0) {
        O2SCL_ERR("Sanity 1 in calc_utf8.",o2scl::exc_einval);
      }
      if (debug) {
	std::cout << "value: " << value << std::endl;
      }
      rpn_queue.push(new token32<double>(value,token_num));

    } else if (is_variable_char(expr[i])) {

      std::u32string key;
      key+=expr[i];
      
      while (i+1<expr.length() &&
             (is_variable_char(expr[i+1]) || isdigit(expr[i+1]))) {
        key=key+expr[i+1];
        i++;
      }
        
      bool found=false;
      double val;

      if (key.length()==3 &&
          key[0]=='s' && key[1]=='i' && key[2]=='n') {
	operator_stack.push("sin");
	lastTokenWasOp=true;
      } else if (key.length()==3 &&
                 key[0]=='c' && key[1]=='o' && key[2]=='s') {
	operator_stack.push("cos");
	lastTokenWasOp=true;
      } else if (key.length()==3 &&
                 key[0]=='t' && key[1]=='a' && key[2]=='n') {
	operator_stack.push("tan");
	lastTokenWasOp=true;
      } else if (key.length()==4 && key[0]=='s' &&
                 key[1]=='q' && key[2]=='r' && key[3]=='t') {
	operator_stack.push("sqrt");
	lastTokenWasOp=true;
      } else if (key.length()==3 &&
                 key[0]=='l' && key[1]=='o' && key[2]=='g') {
	operator_stack.push("log");
	lastTokenWasOp=true;
      } else if (key.length()==3 &&
                 key[0]=='e' && key[1]=='x' && key[2]=='p') {
	operator_stack.push("exp");
	lastTokenWasOp=true;
      } else if (key.length()==3 &&
                 key[0]=='a' && key[1]=='b' && key[2]=='s') {
	operator_stack.push("abs");
	lastTokenWasOp=true;
      } else if (key.length()==5 && key[0]=='l' && key[1]=='o' &&
                 key[2]=='g' && key[3]=='1' && key[4]=='0') {
	operator_stack.push("log10");
	lastTokenWasOp=true;
      } else if (key.length()==4 && key[0]=='a' &&
                 key[1]=='s' && key[2]=='i' && key[3]=='n') {
	operator_stack.push("asin");
	lastTokenWasOp=true;
      } else if (key.length()==4 && key[0]=='a' &&
                 key[1]=='c' && key[2]=='o' && key[3]=='s') {
	operator_stack.push("acos");
	lastTokenWasOp=true;
      } else if (key.length()==4 && key[0]=='a' &&
                 key[1]=='t' && key[2]=='a' && key[3]=='n') {
	operator_stack.push("atan");
	lastTokenWasOp=true;
      } else if (key.length()==4 && key[0]=='s' &&
                 key[1]=='i' && key[2]=='n' && key[3]=='h') {
	operator_stack.push("sinh");
	lastTokenWasOp=true;
      } else if (key.length()==4 && key[0]=='c' &&
                 key[1]=='o' && key[2]=='s' && key[3]=='h') {
	operator_stack.push("cosh");
	lastTokenWasOp=true;
      } else if (key.length()==4 && key[0]=='t' &&
                 key[1]=='a' && key[2]=='n' && key[3]=='h') {
	operator_stack.push("tanh");
	lastTokenWasOp=true;
      } else if (key.length()==5 && key[0]=='a' && key[1]=='s' &&
                 key[2]=='i' && key[3]=='n' && key[4]=='h') {
	operator_stack.push("asinh");
	lastTokenWasOp=true;
      } else if (key.length()==5 && key[0]=='a' && key[1]=='c' &&
                 key[2]=='o' && key[3]=='s' && key[4]=='h') {
	operator_stack.push("acosh");
	lastTokenWasOp=true;
      } else if (key.length()==5 && key[0]=='a' && key[1]=='t' &&
                 key[2]=='a' && key[3]=='n' && key[4]=='h') {
	operator_stack.push("atanh");
	lastTokenWasOp=true;
      } else if (key.length()==5 && key[0]=='f' && key[1]=='l' &&
                 key[2]=='o' && key[3]=='o' && key[4]=='r') {
	operator_stack.push("floor");
	lastTokenWasOp=true;
      } else {
        
        if (key.length()==4 && key[0]=='t' &&
                 key[1]=='r' && key[2]=='u' && key[3]=='e') {
	  found = true;
	  val = 1;
        } else if (key.length()==5 && key[0]=='f' && key[1]=='a' &&
                   key[2]=='l' && key[3]=='s' && key[4]=='e') {
	  found = true;
	  val = 0;
	} else if (vars) {
	  std::map<std::u32string,double>::const_iterator it=vars->find(key);
	  if (it != vars->end()) {
	    found = true;
	    val = it->second;
	  }
	}
	
	if (found) {
          
	  // Save the number
	  if (debug) {
	    std::cout << "val: " << val << std::endl;
	  }
	  rpn_queue.push(new token32<double>(val,token_num));
          
	} else {
          
	  // Save the variable name:
	  if (debug) {
	    std::cout << "key: ";
            for(size_t ik=0;ik<key.length();ik++) {
              std::cout << key[ik];
            }
            std::cout << std::endl;
	  }
	  rpn_queue.push(new token32<std::u32string>(key,token_var));
	}
	
	lastTokenWasOp=false;
	
      }

    } else {

      // Otherwise, the variable is an operator or parenthesis.
      switch (expr[i]) {
        
      case '(':
        
	operator_stack.push("(");
        i++;
	break;
        
      case ')':
        
	while (operator_stack.top().compare("(")) {
	  rpn_queue.push(new token32<std::string>(operator_stack.top(),
                                                  token_op));
	  operator_stack.pop();
	}
	operator_stack.pop();
	i++;
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
	  ss << expr[i];
          i++;
	  while (i<expr.length() && !isspace(expr[i]) &&
                 !isdigit(expr[i]) && !is_variable_char(expr[i]) &&
                 expr[i] != '(' && expr[i] != ')') {
	    ss << expr[i];
            i++;
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
	      rpn_queue.push(new token32<double>(0,token_num));
	    } else {
	      //"Unrecognized unary operator: '" +
              //str + "'.");
              cleanRPN(rpn_queue);
              return 1;
	    }
	  }

	  while (!operator_stack.empty() &&
		 op_prec[str] >= op_prec[operator_stack.top()]) {
	    rpn_queue.push(new token32<std::string>(operator_stack.top(),
                                                    token_op));
	    operator_stack.pop();
	  }
          
	  operator_stack.push(str);
	  lastTokenWasOp = true;
	}
      }
    }
    
    while (i+1<expr.length() && isspace(expr[i+1])) i++;
    
    // End of while (i<expr.length()) {
  }
  
  while (!operator_stack.empty()) {
    rpn_queue.push(new token32<std::string>(operator_stack.top(), token_op));
    operator_stack.pop();
  }

  rpn_queue2=rpn_queue;
  
  return 0;
}

double calc_utf8::calculate(const std::u32string &expr,
                            const std::map<std::u32string, double> *vars,
                            bool debug) {

  // Convert to RPN with Dijkstra's Shunting-yard algorithm.
  token_queue_t rpn;
  int retx=toRPN_nothrow(expr,vars,debug,op_precedence,rpn);
  if (retx!=0) {
    O2SCL_ERR("Failed.",o2scl::exc_efailed);
  }

  double ret = calculate(rpn);

  cleanRPN(rpn);

  return ret;
}

int calc_utf8::calculate_nothrow(const std::u32string &expr,
                                 const std::map<std::u32string, double> *vars,
                                 bool debug, double &result) {

  // Convert to RPN with Dijkstra's Shunting-yard algorithm.
  token_queue_t rpn;
  int ret1=toRPN_nothrow(expr,vars,debug,op_precedence,rpn);
  if (ret1!=0) return ret1;

  double ret;
  int ret2=calculate_nothrow(rpn,0,ret);
  if (ret2!=0) return 10+ret2;

  cleanRPN(rpn);
  result=ret;

  return 0;
}

double calc_utf8::calculate(token_queue_t rpn,
                            const std::map<std::u32string, double> *vars) {

  // Evaluate the expression in RPN form.
  std::stack<double> evaluation;
  while (!rpn.empty()) {
    token_base *base=rpn.front();
    rpn.pop();

    // Operator:
    if (base->type == token_op) {
      token32<std::string>* strTok = static_cast<token32<std::string>*>(base);
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
      } else if (!str.compare("floor")) {
	evaluation.push(floor(right));
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
    } else if (base->type == token_num) { // Number
      token32<double>* doubleTok = static_cast<token32<double>*>(base);
      evaluation.push(doubleTok->val);
    } else if (base->type == token_var) { // Variable
      if (!vars) {
        throw std::domain_error
	  ("Detected variable, but the variable map is null.");
      }
      
      token32<std::u32string>* strTok =
        static_cast<token32<std::u32string>*>(base);
      
      std::u32string key = strTok->val;
      std::map<std::u32string, double>::const_iterator it = vars->find(key);
      
      if (it == vars->end()) {
        O2SCL_ERR("Unable to find variable.",o2scl::exc_efailed);
        //throw std::domain_error("Unable to find the variable '"
        //+ key + "'.");
      }
      evaluation.push(it->second);
    } else {
      throw std::domain_error("Invalid token.");
    }
  }
  return evaluation.top();
}

int calc_utf8::calculate_nothrow(token_queue_t rpn,
                                 const std::map<std::u32string, double> *vars,
                                 double &result) {

  // Evaluate the expression in RPN form.
  std::stack<double> evaluation;
  while (!rpn.empty()) {
    token_base* base = rpn.front();
    rpn.pop();

    // Operator:
    if (base->type == token_op) {
      token32<std::string>* strTok =
        static_cast<token32<std::string>*>(base);
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
      } else if (!str.compare("floor")) {
	evaluation.push(floor(right));
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
	  cleanRPN(rpn);
	  return 1;
	  //throw std::domain_error("Unknown operator: '" + str + "'.");
	}
      }
    } else if (base->type == token_num) { // Number
      token32<double>* doubleTok = static_cast<token32<double>*>(base);
      evaluation.push(doubleTok->val);
    } else if (base->type == token_var) { // Variable
      if (!vars) {
	cleanRPN(rpn);
	return 2;
        //throw std::domain_error
	//("Detected variable, but the variable map is null.");
      }
      
      token32<std::u32string>* strTok =
        static_cast<token32<std::u32string>*>(base);
      
      std::u32string key = strTok->val;
      std::map<std::u32string, double>::const_iterator it = vars->find(key);
      
      if (it == vars->end()) {
	cleanRPN(rpn);
	return 3;
        //throw std::domain_error("Unable to find the variable '" + key + "'.");
      }
      evaluation.push(it->second);
    } else {
      cleanRPN(rpn);
      return 4;
      //throw std::domain_error("Invalid token.");
    }
  }
  result=evaluation.top();
  return 0;
}

void calc_utf8::cleanRPN(token_queue_t& rpn) {
  while (rpn.size()) {
    delete rpn.front();
    rpn.pop();
  }
  return;
}

void calc_utf8::compile(const std::u32string &expr,
                        const std::map<std::u32string, double> *vars,
                        bool debug,
                        std::map<std::string, int> op_prec) {

  // Make sure it is empty:
  cleanRPN(this->RPN);

  int retx=calc_utf8::toRPN_nothrow(expr,vars,debug,op_prec,this->RPN);
  if (retx!=0) {
    O2SCL_ERR("Failed.",o2scl::exc_einval);
  }

  return;
}

int calc_utf8::compile_nothrow(const std::u32string &expr,
                               const std::map<std::u32string, double> *vars,
                               bool debug,
                               std::map<std::string, int> op_prec) {

  // Make sure it is empty:
  cleanRPN(this->RPN);

  int ret=calc_utf8::toRPN_nothrow(expr,vars,debug,op_prec,this->RPN);
  return ret;
}

double calc_utf8::eval(const std::map<std::u32string, double> *vars) {
  return calculate(this->RPN, vars);
}

int calc_utf8::eval_nothrow(const std::map<std::u32string, double> *vars,
                            double &result) {
  int ret=calculate_nothrow(this->RPN,vars,result);
  return ret;
}

std::string calc_utf8::RPN_to_string() {
  std::stringstream ss;
  token_queue_t rpn = this->RPN;

  ss << "calc_utf8 { RPN: [ ";
  while (rpn.size()) {
    token_base* base = rpn.front();

    token32<double>* doubleTok = dynamic_cast<token32<double>*>(base);
    if(doubleTok) {
      ss << doubleTok->val;
    }

    token32<std::string>* strTok = dynamic_cast<token32<std::string>*>(base);
    if(strTok) {
      ss << "'" << strTok->val << "'";
    }

    rpn.pop();

    ss << (rpn.size() ? ", ":"");
  }
  ss << " ] }";
  return ss.str();
}
