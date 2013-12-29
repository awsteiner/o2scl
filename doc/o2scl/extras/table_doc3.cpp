for(size_t i=0;i<t.get_ncolumns();i++) {
  cout << t.get_column_name(i) << " ";
}
cout << endl;
for(size_t i=0;i<t.get_ncolumns();i++) {
  cout << t.get(i,1) << " ";
}
cout << endl;
