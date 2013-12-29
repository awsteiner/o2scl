class my_class {
protected:
  ovector a;
public:
  int afunction() {
    a.allocate(1);
    // do something with a
    return 0;
  }
  int get_result(const ovector_view &av) { av=a; return 0; }
};
