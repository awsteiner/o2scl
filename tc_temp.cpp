#include <hdf5.h>

#include <o2scl/test_mgr.h>
#include <o2scl/hdf_file.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

int main(void) {

  system("h5dump /usr/local/share/o2scl/fermion_cal.o2");
  system("h5dump /usr/local/share/o2scl//fermion_cal.o2");
  
  hid_t file;
  file=H5Fopen("/usr/local/share/o2scl//fermion_cal.o2",
	       H5F_ACC_RDWR,H5P_DEFAULT);
  cout << "file: " << file << endl;
  hid_t group=H5Gopen1(file,"acol");
  cout << "group: " << group << endl;
  hid_t dset=H5Dopen1(group,"o2scl_type");
  cout << "dset: " << dset << endl;

  return 0;
}

