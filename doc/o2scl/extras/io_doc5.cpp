virtual int object_in_mem(cinput *cin, in_file_format *ins,
			  object *&op, std::string &name);
virtual int object_in_mem(cinput *cin, in_file_format *ins, object *&op,
			  int &sz, std::string &name);
virtual int object_in_mem(cinput *cin, in_file_format *ins,
			  object **&op, int &sz, int &sz2,
			  std::string &name);
template<size_t N>
int object_in_mem(cinput *co, in_file_format *ins,
		  object op[][N], int &sz, std::string &name);
