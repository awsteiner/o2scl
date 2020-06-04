class my_type_t {
  virtual member_func();
};
my_type_t my_instance;
class my_derived_type_t : public my_type_t {
  virtual member_func();
};
my_derived_type_t my_inst2;
mm_funct_mfptr<my_type_t> func(&my_inst2,&my_instance::member_func);
