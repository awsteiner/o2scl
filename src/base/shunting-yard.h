// Source: http://www.daniweb.com/software-development/cpp/code/427500/calculator-using-shunting-yard-algorithm#
// Author: Jesse Brown
// Modifications: Brandon Amos

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

#ifndef _SHUNTING_YARD_H
#define _SHUNTING_YARD_H

#include <map>
#include <stack>
#include <string>
#include <queue>

enum tokType { NONE, OP, VAR, NUM };

struct TokenBase {
  tokType type;
  virtual ~TokenBase() {}
};

template<class T> class Token : public TokenBase {
public:
  Token (T t, tokType type) : val(t) { this->type=type; }
  T val;
};

typedef std::queue<TokenBase*> TokenQueue_t;

class calculator {
private:
  static std::map<std::string, int> opPrecedence;
  static std::map<std::string, int> buildOpPrecedence();

public:
  static double calculate(const char* expr,
      std::map<std::string, double>* vars = 0);

private:
  static double calculate(TokenQueue_t RPN,
      std::map<std::string, double>* vars = 0);
  static void cleanRPN(TokenQueue_t& rpn);
  static TokenQueue_t toRPN(const char* expr,
			    std::map<std::string, double>* vars);

private:
  TokenQueue_t RPN;
public:
  ~calculator();
  calculator(){}
  calculator(const char* expr,
	     std::map<std::string, double>* vars = 0);
  void compile(const char* expr,
	       std::map<std::string, double>* vars = 0);
  double eval(std::map<std::string, double>* vars = 0);
  std::string str();
};

#endif // _SHUNTING_YARD_H
