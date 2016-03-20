// Original author: Jesse Brown
// Modifications: Brandon Amos, redpois0n
// Additional modifications for use in O2scl by Andrew W. Steiner

/*
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

#include <map>
#include <stack>
#include <string>
#include <queue>

namespace o2scl {

  /** \brief Desc
   */
  enum tokType {NONE,OP,VAR,NUM};

  /** \brief Desc
   */
  struct TokenBase {
    /** \brief Desc
     */
    tokType type;
    virtual ~TokenBase() {}
  };

  /** \brief Desc
   */
  template<class T> class Token : public TokenBase {

  public:
    
    /** \brief Desc
     */
  Token(T t, tokType type) : val(t) { this->type=type; }
    
    /** \brief Desc
     */
    T val;
  };

  typedef std::queue<TokenBase*> TokenQueue_t;

  /** \brief Desc
   */
  class calculator {

  private:

    /** \brief Desc
     */
    static std::map<std::string, int> opPrecedence;

    /** \brief Desc
     */
    static std::map<std::string, int> buildOpPrecedence();

  public:

    /** \brief Desc
     */
    static double calculate(const char* expr,
			    std::map<std::string, double>* vars = 0);

  private:

    /** \brief Desc
     */
    static double calculate(TokenQueue_t RPN,
			    std::map<std::string, double>* vars = 0);

    /** \brief Desc
     */
    static void cleanRPN(TokenQueue_t& rpn);

    /** \brief Desc
     */
    static TokenQueue_t toRPN(const char* expr,
			      std::map<std::string, double>* vars);

  private:

    /** \brief Desc
     */
    TokenQueue_t RPN;

  public:

    ~calculator();
    
    /** \brief Desc
     */
    calculator(){}
    
    /** \brief Desc
     */
    calculator(const char* expr,
	       std::map<std::string, double> *vars=0);
    
    /** \brief Desc
     */
    void compile(const char* expr,
		 std::map<std::string, double> *vars=0);
    
    /** \brief Desc
     */
    double eval(std::map<std::string, double> *vars=0);
    
    /** \brief Desc
     */
    std::string str();
  };

}

// End of "#ifndef O2SCL_SHUNTING_YARD_H"
#endif
