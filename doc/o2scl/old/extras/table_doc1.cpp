// Set the 1st row of column "colx" to 1.0
t.set("colx",0,1.0);
// Set the 2nd row of column "colz" to 2.0
t.set("colz",1,2.0);
// Set the 3rd row of column "coly" to 4.0
t.set("coly",2,4.0);
// This will print out 2.0
cout << t.get("colz",1) << endl;
