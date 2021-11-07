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
#ifndef O2SCL_CALC_UTF8_H
#define O2SCL_CALC_UTF8_H

/** \file calc_utf8.h
    \brief Definitions for \ref o2scl::calc_utf8
*/

#include <map>
#include <stack>
#include <string>
#include <queue>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cmath>

#include <o2scl/rng.h>
#include <o2scl/err_hnd.h>
#include <o2scl/string_conv.h>
#include <o2scl/find_constants.h>

namespace o2scl {

  /** \brief Evaluate a mathematical expression stored in a UTF8 string
      
      This is based on Brandon Amos' code at 
      https://github.com/bamos/cpp-expression-parser
      in turn based on Jesse Brown's code at
      http://www.daniweb.com/software-development/cpp/code/427500/calculator-using-shunting-yard-algorithm .

      The original code has been modified for use in \o2 .
      
      \verbatim embed:rst

      .. todo::

      In class calc_utf8:

      - Future: Add functions atan2, cot, csc, ceil, floor, max, min,
      and maybe if?

      - Future: There is some code duplication across the
      functions, especially with regard to conversion between UTF8
      and char32, which could be removed.

      \endverbatim
  */
  template<class fp_t=double> class calc_utf8 {

  public:
    
    /** \brief Token base data type for \ref o2scl::calc_utf8
     */
    struct token_base {
      
      /** \brief The token type
       */
      int type;
      
      virtual ~token_base() {}
      
    };
    
    /** \brief Token class for \ref o2scl::calc_utf8
     */
    template<class T> class token32 : public token_base {
      
    public:
      
      /** \brief Create a token of type \c type with value \c t
          
          \comment
          AWS: renamed type to typex to prevent shadow warnings
          \endcomment
      */
      token32(T t, int typex) : val(t) { this->type=typex; }
      
      /** \brief The actual value stored
       */
      T val;
      
    };
    
  protected:

    rng<> r;
    
    /** \brief A typedef for a queue of tokens for \ref o2scl::calculator
     */
    typedef std::queue<token_base *> token_queue_t;
    
    /// \name Token types
    //@{
    static const int token_none=0;
    static const int token_op=1;
    static const int token_var=2;
    static const int token_num=3;
    //@}

    /** \brief A map denoting operator precedence
     */
    std::map<std::string, int> op_precedence;
    
    /** \brief Build the operator precedence map
     */
    std::map<std::string, int> build_op_precedence() {
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

    /** \brief Return true if \c is a variable
     */
    bool is_variable_char(const char32_t c) {
      // Presume it's a variable if it's not an operator or digit
      if (c=='^' || c=='*' || c=='/' || c=='%' || c=='+' || c=='-' ||
          c=='<' || c=='=' || c=='>' || c=='!' || c=='&' || c=='|' ||
          c=='(' || c==')' || isdigit(c) || c=='.' || c==',' || isspace(c)) {
        return false;
      }
      return true;
    }      
    
    /** \brief Compile and evaluate the expression in \c RPN using
	definitions in \c vars and return an integer to indicate
	success or failure
    */
    int calc_RPN_nothrow(token_queue_t rpn,
                         const std::map<std::u32string, fp_t> *vars,
                         fp_t &result) {

      // Evaluate the expression in RPN form.
      std::stack<fp_t> evaluation;
      while (!rpn.empty()) {
        token_base *base=rpn.front();
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
          fp_t right = evaluation.top();
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
            evaluation.push(abs(right));
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
            fp_t left  = evaluation.top();
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
              //throw std::domain_error("Unknown operator: '" + str + "'.");
              return 2;
            }
          }
        } else if (base->type == token_num) { // Number
          token32<fp_t> *doubleTok = static_cast<token32<fp_t>*>(base);
          evaluation.push(doubleTok->val);
        } else if (base->type == token_var) { // Variable
          if (!vars) {
            //throw std::domain_error
            //("Detected variable, but the variable map is null.");
            return 3;
          }
      
          token32<std::u32string> *strTok=
            static_cast<token32<std::u32string>*>(base);
      
          std::u32string key=strTok->val;
          if (key.length()==4 && ((char)key[0])=='r' &&
              ((char)key[1])=='a' && ((char)key[2])=='n' &&
              ((char)key[3])=='d') {
            evaluation.push(r.random());
          } else {
        
            typename std::map<std::u32string, fp_t>::const_iterator it=
              vars->find(key);
        
            if (it == vars->end()) {
              //O2SCL_ERR("Unable to find variable.",o2scl::exc_efailed);
              //throw std::domain_error("Unable to find the variable '"
              //+ key + "'.");
              return 4;
            }
            evaluation.push(it->second);
          }
        } else {
          //throw std::domain_error("Invalid token.");
          return 5;
        }
      }
      result=evaluation.top();
      return 0;
    }      
    
    /** \brief Empty and free memory associated with \c rpn
        
	\note This is called by the destructor to free the memory in
	\ref RPN .
    */
    void cleanRPN(token_queue_t &rpn) {
      while (rpn.size()) {
        delete rpn.front();
        rpn.pop();
      }
      return;
    }
    
    /** \brief Convert the expression in \c expr to RPN and return an
	integer to indicate success or failure
    */
    int toRPN_nothrow(const std::u32string &expr,
                      const std::map<std::u32string, fp_t> *vars,
                      std::map<std::string, int> op_prec,
                      token_queue_t &rpn_queue2) {
      
      token_queue_t rpn_queue;
      std::stack<std::string> operator_stack;
      bool last_token_was_op = true;

      if (verbose>=2) {
        std::cout << "In toRPN_nothrow(), processing expression: ";
        std::string stmp;
        char32_to_utf8(expr,stmp);
        std::cout << stmp << std::endl;
      }
  
      // In one pass, ignore whitespace and parse the expression into RPN
      // using Dijkstra's Shunting-yard algorithm.
  
      size_t i=0;
      while (i<expr.length()) {
    
        if (isdigit(expr[i])) {
      
          if (verbose>=2) {
            std::cout << "In toRPN_nothrow(), found digit at character "
                 << i << std::endl;
          }
      
          std::u32string str_num;
          str_num=str_num+expr[i];

          bool exponent=false;
          // The idea here is that this will handle both commas or
          // periods, depending on locale, but this hasn't been tested
          while (i+1<expr.length() &&
                 (expr[i+1]=='.' || expr[i+1]==',' ||
                  isdigit(expr[i+1]) || expr[i+1]=='e' ||
                  expr[i+1]=='E' || (exponent && (expr[i+1]=='+' ||
                                                  expr[i+1]=='-')))) {
            if (expr[i+1]=='e' || expr[i+1]=='E') exponent=true;
            str_num=str_num+expr[i+1];
            i++;
          }
      
          if (verbose>=2) {
            std::cout << "In toRPN_nothrow(), str_num: ";
            std::string stmp;
            char32_to_utf8(str_num,stmp);
            std::cout << stmp << std::endl;
          }
      
          fp_t value;
          int iret=s32tod_nothrow(str_num,value);
          if (verbose>=1) {
            std::cout << "In toRPN_nothrow(), value: " << value << std::endl;
          }
      
          if (iret!=0) {
            O2SCL_ERR("Sanity 1 in calc_utf8.",o2scl::exc_einval);
          }
          rpn_queue.push(new token32<fp_t>(value,token_num));

          last_token_was_op=false;
      
        } else if (is_variable_char(expr[i])) {

          if (verbose>=2) {
            std::cout << "In toRPN_nothrow(), found variable at character "
                 << i << std::endl;
          }
      
          std::u32string key;
          key+=expr[i];
      
          while (i+1<expr.length() &&
                 (is_variable_char(expr[i+1]) || isdigit(expr[i+1]))) {
            key=key+expr[i+1];
            i++;
          }
        
          if (verbose>=2) {
            std::cout << "In toRPN_nothrow(), key: ";
            std::string stmp;
            char32_to_utf8(key,stmp);
            std::cout << stmp << std::endl;
          }
      
          bool found=false;
          fp_t val;

          if (key.length()==3 &&
              key[0]=='s' && key[1]=='i' && key[2]=='n') {
            operator_stack.push("sin");
            last_token_was_op=true;
          } else if (key.length()==3 &&
                     key[0]=='c' && key[1]=='o' && key[2]=='s') {
            operator_stack.push("cos");
            last_token_was_op=true;
          } else if (key.length()==3 &&
                     key[0]=='t' && key[1]=='a' && key[2]=='n') {
            operator_stack.push("tan");
            last_token_was_op=true;
          } else if (key.length()==4 && key[0]=='s' &&
                     key[1]=='q' && key[2]=='r' && key[3]=='t') {
            operator_stack.push("sqrt");
            last_token_was_op=true;
          } else if (key.length()==3 &&
                     key[0]=='l' && key[1]=='o' && key[2]=='g') {
            operator_stack.push("log");
            last_token_was_op=true;
          } else if (key.length()==3 &&
                     key[0]=='e' && key[1]=='x' && key[2]=='p') {
            operator_stack.push("exp");
            last_token_was_op=true;
          } else if (key.length()==3 &&
                     key[0]=='a' && key[1]=='b' && key[2]=='s') {
            operator_stack.push("abs");
            last_token_was_op=true;
          } else if (key.length()==5 && key[0]=='l' && key[1]=='o' &&
                     key[2]=='g' && key[3]=='1' && key[4]=='0') {
            operator_stack.push("log10");
            last_token_was_op=true;
          } else if (key.length()==4 && key[0]=='a' &&
                     key[1]=='s' && key[2]=='i' && key[3]=='n') {
            operator_stack.push("asin");
            last_token_was_op=true;
          } else if (key.length()==4 && key[0]=='a' &&
                     key[1]=='c' && key[2]=='o' && key[3]=='s') {
            operator_stack.push("acos");
            last_token_was_op=true;
          } else if (key.length()==4 && key[0]=='a' &&
                     key[1]=='t' && key[2]=='a' && key[3]=='n') {
            operator_stack.push("atan");
            last_token_was_op=true;
          } else if (key.length()==4 && key[0]=='s' &&
                     key[1]=='i' && key[2]=='n' && key[3]=='h') {
            operator_stack.push("sinh");
            last_token_was_op=true;
          } else if (key.length()==4 && key[0]=='c' &&
                     key[1]=='o' && key[2]=='s' && key[3]=='h') {
            operator_stack.push("cosh");
            last_token_was_op=true;
          } else if (key.length()==4 && key[0]=='t' &&
                     key[1]=='a' && key[2]=='n' && key[3]=='h') {
            operator_stack.push("tanh");
            last_token_was_op=true;
          } else if (key.length()==5 && key[0]=='a' && key[1]=='s' &&
                     key[2]=='i' && key[3]=='n' && key[4]=='h') {
            operator_stack.push("asinh");
            last_token_was_op=true;
          } else if (key.length()==5 && key[0]=='a' && key[1]=='c' &&
                     key[2]=='o' && key[3]=='s' && key[4]=='h') {
            operator_stack.push("acosh");
            last_token_was_op=true;
          } else if (key.length()==5 && key[0]=='a' && key[1]=='t' &&
                     key[2]=='a' && key[3]=='n' && key[4]=='h') {
            operator_stack.push("atanh");
            last_token_was_op=true;
          } else if (key.length()==5 && key[0]=='f' && key[1]=='l' &&
                     key[2]=='o' && key[3]=='o' && key[4]=='r') {
            operator_stack.push("floor");
            last_token_was_op=true;
          } else {

            if (key.length()==4 && ((char)key[0])=='t' &&
                ((char)key[1])=='r' && ((char)key[2])=='u' &&
                ((char)key[3])=='e') {
              found=true;
              val=1;
            } else if (key.length()==5 && ((char)key[0])=='f' &&
                       ((char)key[1])=='a' && ((char)key[2])=='l' &&
                       ((char)key[3])=='s' && ((char)key[4])=='e') {
              found=true;
              val=0;
            } else if (key.length()==4 && ((char)key[0])=='r' &&
                       ((char)key[1])=='a' && ((char)key[2])=='n' &&
                       ((char)key[3])=='d') {
              found=false;
            } else if (vars) {
              typename std::map<std::u32string,fp_t>::const_iterator it=
                vars->find(key);
              if (it != vars->end()) {
                found = true;
                val = it->second;
              }
            }
	
            if (found) {
          
              // Save the number
              if (verbose>=1) {
                std::cout << "In toRPN_nothrow(), value: " << val << std::endl;
              }
              rpn_queue.push(new token32<fp_t>(val,token_num));
          
            } else {
          
              // Save the variable name:
              if (verbose>=1) {
                std::cout << "In toRPN_nothrow(), variable: ";

                std::string stmp;
                char32_to_utf8(key,stmp);
                std::cout << stmp << std::endl;
            
              }
              rpn_queue.push(new token32<std::u32string>(key,token_var));
            }
	
            last_token_was_op=false;
	
          }

        } else {

          if (verbose>=1) {
            std::cout << "In toRPN_nothrow(), operator or parenthesis: "
                      << ((char)expr[i]) << std::endl;
          }
      
          // Otherwise, the variable is an operator or parenthesis.
          switch (expr[i]) {
        
          case '(':
        
            operator_stack.push("(");
            break;
        
          case ')':
        
            while (operator_stack.top().compare("(")) {
              rpn_queue.push(new token32<std::string>(operator_stack.top(),
                                                      token_op));
              operator_stack.pop();
            }
            operator_stack.pop();
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
              ss << ((char)expr[i]);
              while (i+1<expr.length() && !isspace(expr[i+1]) &&
                     !isdigit(expr[i+1]) && !is_variable_char(expr[i+1]) &&
                     expr[i+1] != '(' && expr[i+1] != ')') {
                ss << ((char)expr[i+1]);
                i++;
              }
              ss.clear();
              std::string str;
              ss >> str;
          
              if (verbose>=1) {
                std::cout << "In toRPN_nothrow(), processing operator: "
                          << str << std::endl;
              }

              if (last_token_was_op) {
                if (verbose>=1) {
                  std::cout << "In toRPN_nothrow(), last token was an operator: "
                       << str.compare("-") << " " << str.compare("+") << std::endl;
                }
                // Convert unary operators to binary in the RPN.
                if (!str.compare("-") || !str.compare("+")) {
                  rpn_queue.push(new token32<fp_t>(0,token_num));
                } else {
                  //"Unrecognized unary operator: '" +
                  //str + "'.");
                  if (verbose>=2) {
                    std::cout << "In toRPN_nothrow(), "
                         << "unrecognized unary operator." << std::endl;
                  }
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
              last_token_was_op=true;
            }
          }
        }

        //char ch;
        //cin >> ch;
    
        // Skip spaces if necessary
        while (i+1<expr.length() && isspace(expr[i+1])) i++;
    
        if (verbose>=1) {
          std::cout << "In toRPN_nothrow(), "
                    << "moving to next character." << std::endl;
        }
    
        // Move to the next character
        i++;

        // End of while (i<expr.length()) {
      }
  
      while (!operator_stack.empty()) {
        rpn_queue.push(new token32<std::string>(operator_stack.top(),
                                                token_op));
        operator_stack.pop();
      }
      
      rpn_queue2=rpn_queue;
      
      return 0;
    }      
    
    /** \brief The current expression in RPN
     */
    token_queue_t RPN;

  public:

    /** \brief Create an empty calc_utf8 object
     */
    calc_utf8() {
      // Create the operator precedence object
      op_precedence=calc_utf8::build_op_precedence();
      
      verbose=0;
      
      r.clock_seed();
    }      

    /** \brief Compile expression \c expr using variables 
	specified in \c vars
    */
    calc_utf8(const std::u32string &expr,
              const std::map<std::u32string, fp_t> *vars=0) {
      verbose=0;
  
      // Create the operator precedence object
      op_precedence=calc_utf8::build_op_precedence();
  
      compile(expr,vars);

      r.clock_seed();
    }      
    
    ~calc_utf8() {
      // Clear memory associated with the RPN object
      cleanRPN(this->RPN);
    }      

    /// \name Compile and evaluate 
    //@{
    /** \brief Compile and evaluate \c expr using definitions in 
	\c vars
    */
    fp_t calculate(const std::u32string &expr,
                     const std::map<std::u32string, fp_t> *vars=0) {

      // Convert to RPN with Dijkstra's Shunting-yard algorithm.
      token_queue_t rpn;
      int retx=toRPN_nothrow(expr,vars,op_precedence,rpn);
      if (retx!=0) {
        O2SCL_ERR("Failed.",o2scl::exc_efailed);
      }

      fp_t ret;
      int iret=calc_RPN_nothrow(rpn,0,ret);

      cleanRPN(rpn);

      return ret;
    }
    
    /** \brief Compile and evaluate \c expr using definitions in 
	\c vars
    */
    fp_t calculate(const std::string &expr,
                     const std::map<std::string, fp_t> *vars=0) {
      fp_t ret;
      std::u32string expr2;
      utf8_to_char32(expr,expr2);
      if (vars==0) {
        ret=calculate(expr2,0);
        return ret;
      }
      std::map<std::u32string, fp_t> vars2;
      for(typename std::map<std::string, fp_t>::const_iterator it=
            vars->begin();
          it!=vars->end();it++) {
        std::u32string tmp;
        utf8_to_char32(it->first,tmp);
        vars2.insert(std::make_pair(tmp,it->second));
      }
      ret=calculate(expr2,&vars2);
      return ret;
    }
    
    /** \brief Compile and evaluate \c expr using definitions in 
	\c vars and return an integer to indicate success or failure
    */
    int calculate_nothrow(const std::u32string &expr,
                          const std::map<std::u32string, fp_t> *vars,
                          fp_t &result) {

      // Convert to RPN with Dijkstra's Shunting-yard algorithm.
      token_queue_t rpn;
      int ret1=toRPN_nothrow(expr,vars,op_precedence,rpn);
      if (ret1!=0) return ret1;

      fp_t ret;
      int ret2=calc_RPN_nothrow(rpn,0,ret);
      if (ret2!=0) return 10+ret2;

      cleanRPN(rpn);
      result=ret;

      return 0;
    }      
    
    /** \brief Compile and evaluate \c expr using definitions in 
	\c vars and return an integer to indicate success or failure
    */
    int calculate_nothrow(const std::string &expr,
                          const std::map<std::string, fp_t> *vars,
                          fp_t &result) {

      std::u32string expr2;
      utf8_to_char32(expr,expr2);
  
      std::map<std::u32string, fp_t> vars2;
      if (vars!=0) {
        for(typename std::map<std::string, fp_t>::const_iterator it=
              vars->begin();
            it!=vars->end();it++) {
          std::u32string tmp;
          utf8_to_char32(it->first,tmp);
          vars2.insert(std::make_pair(tmp,it->second));
        }
      }
  
      // Convert to RPN with Dijkstra's Shunting-yard algorithm.
      token_queue_t rpn;
      int ret1;
      if (vars==0) {
        ret1=toRPN_nothrow(expr2,0,op_precedence,rpn);
      } else {
        ret1=toRPN_nothrow(expr2,&vars2,op_precedence,rpn);
      }
      if (ret1!=0) return ret1;

      fp_t ret;
      int ret2=calc_RPN_nothrow(rpn,0,ret);
      if (ret2!=0) return 10+ret2;

      cleanRPN(rpn);
      result=ret;

      return 0;
    }
    //@}

    /// \name Compile functions
    //@{
    /** \brief Compile expression \c expr using variables 
	specified in \c vars and return an
	integer to indicate success or failure
    */
    void compile(const std::u32string &expr,
		 const std::map<std::u32string, fp_t> *vars=0) {

      // Make sure it is empty:
      cleanRPN(this->RPN);

      int retx=calc_utf8::toRPN_nothrow(expr,vars,op_precedence,this->RPN);
      if (retx!=0) {
        O2SCL_ERR("Failed.",o2scl::exc_einval);
      }

      return;
    }      
    
    /** \brief Compile expression \c expr using variables 
	specified in \c vars and return an
	integer to indicate success or failure
    */
    int compile_nothrow(const std::u32string &expr,
			const std::map<std::u32string, fp_t> *vars=0) {

      // Make sure it is empty:
      cleanRPN(this->RPN);

      int ret=calc_utf8::toRPN_nothrow(expr,vars,op_precedence,this->RPN);
      return ret;
    }      
    
    /** \brief Compile expression \c expr using variables 
	specified in \c vars and return an
	integer to indicate success or failure
    */
    void compile(const std::string &expr,
		 const std::map<std::string, fp_t> *vars=0) {

      std::u32string expr2;
      utf8_to_char32(expr,expr2);
  
      std::map<std::u32string, fp_t> vars2;
      if (vars!=0) {
        for(typename std::map<std::string, fp_t>::const_iterator it=
              vars->begin();
            it!=vars->end();it++) {
          std::u32string tmp;
          utf8_to_char32(it->first,tmp);
          vars2.insert(std::make_pair(tmp,it->second));
        }
      }
  
      // Make sure it is empty:
      cleanRPN(this->RPN);

      int retx;
      if (vars==0) {
        retx=calc_utf8::toRPN_nothrow(expr2,0,op_precedence,this->RPN);
      } else {
        retx=calc_utf8::toRPN_nothrow(expr2,&vars2,op_precedence,this->RPN);
      }
      if (retx!=0) {
        O2SCL_ERR("Failed.",o2scl::exc_einval);
      }

      return;
    }      
    
    /** \brief Compile expression \c expr using variables 
	specified in \c vars and return an
	integer to indicate success or failure
    */
    int compile_nothrow(const std::string &expr,
			const std::map<std::string, fp_t> *vars=0) {

      std::u32string expr2;
      utf8_to_char32(expr,expr2);
  
      std::map<std::u32string, fp_t> vars2;
      if (vars!=0) {
        for(typename std::map<std::string, fp_t>::const_iterator it=
              vars->begin();
            it!=vars->end();it++) {
          std::u32string tmp;
          utf8_to_char32(it->first,tmp);
          vars2.insert(std::make_pair(tmp,it->second));
        }
      }
  
      // Make sure it is empty:
      cleanRPN(this->RPN);

      int ret;
      if (vars==0) {
        ret=calc_utf8::toRPN_nothrow(expr2,0,
                                     op_precedence,this->RPN);
      } else {
        ret=calc_utf8::toRPN_nothrow(expr2,&vars2,
                                     op_precedence,this->RPN);
      }
      return ret;
    }
    //@}

    /// \name Evaluate functions
    //@{
    /** \brief Evalate the previously compiled expression using
	variables specified in \c vars
    */
    fp_t eval_char32(const std::map<std::u32string, fp_t> *vars=0) {
      fp_t ret;
      int iret=calc_RPN_nothrow(this->RPN, vars,ret);
      return ret;
    }

    /** \brief Evalate the previously compiled expression using
	variables specified in \c vars
    */
    int eval_char32_nothrow(const std::map<std::u32string, fp_t> *vars,
                            fp_t &result) {
      int ret=calc_RPN_nothrow(this->RPN,vars,result);
      return ret;
    }
    
    /** \brief Evalate the previously compiled expression using
	variables specified in \c vars
    */
    fp_t eval(const std::map<std::string, fp_t> *vars=0) {

      std::map<std::u32string, fp_t> vars2;
      if (vars!=0) {
        for(typename std::map<std::string, fp_t>::const_iterator it=
              vars->begin();
            it!=vars->end();it++) {
          std::u32string tmp;
          utf8_to_char32(it->first,tmp);
          vars2.insert(std::make_pair(tmp,it->second));
        }
      }
  
      fp_t ret;
      int iret;
      if (vars==0) {
        iret=calc_RPN_nothrow(this->RPN,0,ret);
      } else {
        iret=calc_RPN_nothrow(this->RPN,&vars2,ret);
      }
      return ret;
    }

    /** \brief Evalate the previously compiled expression using
	variables specified in \c vars
    */
    int eval_nothrow(const std::map<std::string, fp_t> *vars,
		     fp_t &result) {

      std::map<std::u32string, fp_t> vars2;
      if (vars!=0) {
        for(typename std::map<std::string, fp_t>::const_iterator it=
              vars->begin();
            it!=vars->end();it++) {
          std::u32string tmp;
          utf8_to_char32(it->first,tmp);
          vars2.insert(std::make_pair(tmp,it->second));
        }
      }
  
      int ret;
      if (vars==0) {
        ret=calc_RPN_nothrow(this->RPN,0,result);
      } else {
        ret=calc_RPN_nothrow(this->RPN,&vars2,result);
      }
  
      return ret;
    }
    //@}
    
    /// Verbosity parameter
    int verbose;

    /// \name Other functions
    //@{
    /** \brief Convert the RPN expression to a string

	\note This is mostly useful for debugging
    */
    std::string RPN_to_string() {
  
      std::stringstream ss;
      token_queue_t rpn = this->RPN;

      ss << "calc_utf8 { RPN: [ ";
  
      while (rpn.size()) {
    
        token_base *base=rpn.front();

        token32<fp_t> *doubleTok=dynamic_cast<token32<fp_t>*>(base);
        if (doubleTok) {
          ss << doubleTok->val;
        }
    
        token32<std::string> *strTok=dynamic_cast<token32<std::string>*>(base);
        if (strTok) {
          ss << "'" << strTok->val << "'";
        }

        token32<std::u32string> *str32Tok=
          dynamic_cast<token32<std::u32string>*>(base);
        if (str32Tok) {
          std::string stmp;
          char32_to_utf8(str32Tok->val,stmp);
          ss << "'" << stmp << "'";
        }

        rpn.pop();

        // Output ", " only if we're not at the end
        ss << (rpn.size() ? ", ":"");
      }
  
      ss << " ] }";
  
      return ss.str();
    }
  
    /** \brief Get a copy of the RPN version

	\warning This copy contains raw pointers which are invalid
	if changed.
    */
    token_queue_t get_RPN() {
      return this->RPN;
    }      

    /** \brief Get the variable list
     */
    std::vector<std::u32string> get_var_list() {
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
    //@}
    
  };
  
}

// End of "#ifndef O2SCL_CALC_UTF8_H"
#endif
