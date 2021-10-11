// Simple code to write the function to a file
int write_file(double x1) {
  my_class c;
  c.set_parameter();
  ofstream fout;
  double p=1.1;
  fout.open("ex_fptr.out");
  fout << "x y" << endl;
  for(double x=-1.0;x<=2.00001;x+=0.001) {
    fout << x << " " << c.member_function_parameter(x,p) << endl;
  }
  fout.close();
  return 0;
}
