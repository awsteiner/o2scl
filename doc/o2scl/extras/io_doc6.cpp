virtual int object_out(coutput *cout, out_file_format *outs, object *op,
		       int sz=0, std::string name="");
virtual int object_out(coutput *cout, out_file_format *outs, object **op,
		       int sz, int sz2, std::string name="");
template<size_t N>
int object_out(coutput *cout, out_file_format *outs,
	       object op[][N], int sz, std::string name="");
