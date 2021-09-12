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
#ifndef O2SCL_SHUNTING_YARD_H
#define O2SCL_SHUNTING_YARD_H

/** \file shunting_yard.h
    \brief Definitions for \ref o2scl::calculator
*/

#include <map>
#include <stack>
#include <string>
#include <queue>

namespace o2scl {

  /** \brief Token list for \ref o2scl::calculator
   */
  enum tokType {NONE,OP,VAR,NUM};

  /** \brief Token base data type for \ref o2scl::calculator
   */
  struct TokenBase {
    /** \brief The token type
     */
    tokType type;
    virtual ~TokenBase() {}
  };

  /** \brief Token class for \ref o2scl::calculator
   */
  template<class T> class Token : public TokenBase {

  public:
    
    /** \brief Create a token of type \c type with value \c t

	\comment
	AWS: renamed type to typex to prevent shadow warnings
	\endcomment
     */
  Token(T t, tokType typex) : val(t) { this->type=typex; }
    
    /** \brief The actual value stored
     */
    T val;
  };

  /** \brief A typedef for a queue of tokens for \ref o2scl::calculator
   */
  typedef std::queue<TokenBase*> TokenQueue_t;

  /** \brief Evaluate a mathematical expression in a string

      \note This class will eventually be deprecated in favor
      of \ref calc_utf8 .

      This is based on Brandon Amos' code at 
      https://github.com/bamos/cpp-expression-parser
      in turn based on Jesse Brown's code at
      http://www.daniweb.com/software-development/cpp/code/427500/calculator-using-shunting-yard-algorithm .

      The original code has been modified for use in \o2 .

      \future Add functions atan2, cot, csc, ceil, floor, int, max, min,
      and maybe if?

      \warning The nothrow() functions are a naive attempt at 
      more detailed error handling than Amos' original code. I have
      tried to make sure they don't create any memory leaks, but 
      I have not fully tested this. 
   */
  class calculator {

  private:

    /** \brief A map denoting operator precedence
     */
    static std::map<std::string, int> opPrecedence;

    /** \brief Build the operator precedence map
     */
    static std::map<std::string, int> buildOpPrecedence();

    /** \brief Return true if \c is a variable
     */
    static bool isvariablechar(char c);
    
  public:

    static int nt_debug;
    
    /** \brief Compile and evaluate \c expr using definitions in 
	\c vars
     */
    static double calculate(const char* expr,
			    const std::map<std::string, double>* vars = 0,
			    bool debug=false);
    
    /** \brief Compile and evaluate \c expr using definitions in 
	\c vars and return an integer to indicate success or failure
    */
    static int calculate_nothrow(const char* expr,
				 const std::map<std::string, double>* vars,
				 bool debug, double &result);
    
  private:

    /** \brief Compile and evaluate the expression in \c RPN
	using definitions in \c vars
     */
    static double calculate(TokenQueue_t RPN,
			    const std::map<std::string, double>* vars = 0);
    
    /** \brief Compile and evaluate the expression in \c RPN using
	definitions in \c vars and return an integer to indicate
	success or failure
     */
    static int calculate_nothrow(TokenQueue_t RPN,
				 const std::map<std::string, double>* vars,
				 double &result);

    /** \brief Empty and free memory associated with \c rpn

	\note This is called by the destructor to free the memory in
	\ref RPN .
     */
    static void cleanRPN(TokenQueue_t &rpn);
    
    /** \brief Convert the expression in \c expr to RPN 
     */
    static TokenQueue_t toRPN(const char* expr,
			      const std::map<std::string, double>* vars,
			      bool debug=false,
			      std::map<std::string, int> opPrec=opPrecedence);

    /** \brief Convert the expression in \c expr to RPN and return an
	integer to indicate success or failure
     */
    static int toRPN_nothrow(const char* expr,
			     const std::map<std::string, double>* vars,
			     bool debug,
			     std::map<std::string, int> opPrec,
			     TokenQueue_t &queue2);

  private:

    /** \brief The current expression in RPN
     */
    TokenQueue_t RPN;

  public:

    ~calculator();
    
    /** \brief Create an empty calculator object
     */
    calculator(){
      nt_debug=0;
    }
    
    /** \brief Compile expression \c expr using variables 
	specified in \c vars
     */
    calculator(const char* expr,
	       const std::map<std::string, double> *vars=0,
	       bool debug=false,
	       std::map<std::string, int> opPrec=opPrecedence);
    
    /** \brief Compile expression \c expr using variables 
	specified in \c vars and return an
	integer to indicate success or failure
    */
    void compile(const char* expr,
		 const std::map<std::string, double> *vars=0,
		 bool debug=false,
		 std::map<std::string, int> opPrec=opPrecedence);
    
    /** \brief Compile expression \c expr using variables 
	specified in \c vars and return an
	integer to indicate success or failure
    */
    int compile_nothrow(const char* expr,
			const std::map<std::string, double> *vars=0,
			bool debug=false,
			std::map<std::string, int> opPrec=opPrecedence);
    
    /** \brief Evalate the previously compiled expression using
	variables specified in \c vars
     */
    double eval(const std::map<std::string, double> *vars=0);

    /** \brief Evalate the previously compiled expression using
	variables specified in \c vars
     */
    int eval_nothrow(const std::map<std::string, double> *vars,
		     double &result);
    
    /** \brief Convert the RPN expression to a string

	\note This is mostly useful for debugging
    */
    std::string RPN_to_string();

    /** \brief Get a copy of the RPN version

	\warning This copy contains raw pointers which are invalid
	if the changed.
     */
    TokenQueue_t get_RPN();

    /** \brief Get the variable list
     */
    std::vector<std::string> get_var_list();
    
  };

}

// End of "#ifndef O2SCL_SHUNTING_YARD_H"
#endif
