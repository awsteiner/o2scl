class my_class {
  int afunction(ovector &a) {
    if (a.get_size()>0 && a.is_owner()==true) {
      O2SCL_ERR("Unallocated vector not sent to afunction().",1);
      return 1;
    } else {
      a.allocate(1);
      // do something with a
      return 0;
    }
  }
};
