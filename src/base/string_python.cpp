/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2020-2025, Andrew W. Steiner

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

  ───────────────────────────────────────────────────────────────────
*/

#include <o2scl/string_python.h>

using namespace std;
using namespace o2scl;

void *o2scl_create_std_string() {
  std::string *ptr=new std::string;
  return ptr;
}

void o2scl_free_std_string(void *vptr) {
  std::string *ptr=(std::string *)vptr;
  delete ptr;
  return;
}

void o2scl_copy_std_string(void *vsrc, void *vdest) {
  std::string *src=(std::string *)vsrc;
  std::string *dest=(std::string *)vdest;
  *dest=*src;
}

size_t o2scl_std_string_length(void *vptr) {
  std::string *ptr=(std::string *)vptr;
  size_t ret=ptr->length();
  return ret;
}

char o2scl_std_string_getitem(void *vptr, size_t n) {
  std::string *ptr=(std::string *)vptr;
  /* tag 4 */ char ret=ptr->operator[](n);
  return ret;
}

void o2scl_std_string_setitem(void *vptr, size_t i, char val) {
  std::string *ptr=(std::string *)vptr;
  (*ptr)[i]=val;
  return;
}

void o2scl_std_string_resize(void *vptr, size_t n) {
  std::string *ptr=(std::string *)vptr;
  ptr->resize(n);
  return;
}

