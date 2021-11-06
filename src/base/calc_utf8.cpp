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

/*
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cmath>

#include <o2scl/calc_utf8.h>
#include <o2scl/err_hnd.h>
#include <o2scl/string_conv.h>

using namespace std;
using namespace o2scl;

calc_utf8::calc_utf8() {
}
  
calc_utf8::calc_utf8(const std::u32string &expr,
                     const std::map<std::u32string, double> *vars) {
}

calc_utf8::~calc_utf8() {
}

std::queue<calc_utf8::token_base *> calc_utf8::get_RPN() {
}

std::vector<std::u32string> calc_utf8::get_var_list() {
}  

std::map<std::string, int> calc_utf8::build_op_precedence() {
}

bool calc_utf8::is_variable_char(const char32_t c) {
}

int calc_utf8::toRPN_nothrow(const std::u32string &expr,
                             const std::map<std::u32string, double> *vars,
                             std::map<std::string, int> op_prec,
                             token_queue_t &rpn_queue2) {
}

double calc_utf8::calculate(const std::string &expr,
                            const std::map<std::string, double> *vars) {
}

double calc_utf8::calculate(const std::u32string &expr,
                            const std::map<std::u32string, double> *vars) {
}

int calc_utf8::calculate_nothrow(const std::u32string &expr,
                                 const std::map<std::u32string, double> *vars,
                                 double &result) {
}

int calc_utf8::calculate_nothrow(const std::string &expr,
                                 const std::map<std::string, double> *vars,
                                 double &result) {
}

int calc_utf8::calc_RPN_nothrow(token_queue_t rpn,
                                const std::map<std::u32string,
                                double> *vars,
                                double &result) {
}

void calc_utf8::cleanRPN(token_queue_t& rpn) {

void calc_utf8::compile(const std::u32string &expr,
                        const std::map<std::u32string, double> *vars) {
}

int calc_utf8::compile_nothrow(const std::u32string &expr,
                               const std::map<std::u32string,
                               double> *vars) {
}

void calc_utf8::compile(const std::string &expr,
                        const std::map<std::string, double> *vars) {
}

int calc_utf8::compile_nothrow(const std::string &expr,
                               const std::map<std::string, double> *vars) {
}

double calc_utf8::eval_char32(const std::map<std::u32string, double> *vars) {

int calc_utf8::eval_char32_nothrow(const std::map<std::u32string, double> *vars,
                            double &result) {

double calc_utf8::eval(const std::map<std::string, double> *vars) {

int calc_utf8::eval_nothrow(const std::map<std::string, double> *vars,
                            double &result) {

std::string calc_utf8::RPN_to_string() {
*/
