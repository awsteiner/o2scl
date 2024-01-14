/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2011-2024, Andrew W. Steiner
  
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
#include <o2scl/test_mgr.h>
#include <o2scl/rng.h>
#include <o2scl/expval.h>
#include <o2scl/constants.h>

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

using namespace std;
using namespace o2scl;

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);
  cout.precision(4);

  test_mgr t;
  t.set_output_level(1);

  rng<> gr;
  int seed=10;
  if (argc>=2) {
    seed=o2scl::stoi(argv[1]);
  }
  gr.set_seed(seed);

  double avg, sd, avge;

  size_t i_block, i_per_block;
  size_t m_block, m_per_block;

  // ------------------------------------------------------------------
  // Test expval_scalar with an odd number of blocks
  // ------------------------------------------------------------------

  if (true) {

    expval_scalar se(9,5);

    cout << "n_per_block=5: " << endl;
    cout << "x          i+1 ib ipb mb mpb prog       f "
	 << "avg        sd         avge         data" << endl;

    for(size_t i=0;i<100;i++) {
      double x=((double)(i/10))+1.0;

      // Add measurement and get info
      se.add(x);
      se.get_block_indices(i_block,i_per_block);

      // Output
      se.current_avg_stats(avg,sd,avge,m_block,m_per_block);
    
      {
	cout << x << " ";
	cout.width(3);
	cout << i+1 << " ";
	cout.width(2);
	cout << i_block << " ";
	cout.width(2);
	cout << i_per_block << " ";
	cout.width(2);
	cout << m_block << "   ";
	cout.width(2);
	cout << m_per_block << " ";
	cout << se.progress() << " " << se.finished() << " "
	     << avg << " " << sd << " " << avge << " : ";
	cout.unsetf(ios::scientific);
	cout << endl;
	//cout << se.get_data() << endl;
	cout.setf(ios::scientific);
      }

      if (i==44 || i==89) {
	t.test_gen(se.finished()==true,"fin");
      } else {
	t.test_gen(se.finished()==false,"fin");
      }

      if (i<=44) {
	for(size_t k=0;k<9;k++) {
	  if (se.get_data()[k]>0.0) {
	    t.test_rel(se.get_data()[k],((double)(k/2+1)),1.0e-12,"data 1");
	  }
	}
	t.test_rel(se.progress(),((double)(i+1))/45.0,1.0e-12,"progress 1");
      } else if (i<=89) {
	for(size_t k=0;k<9;k++) {
	  if (se.get_data()[k]>0.0) {
	    t.test_rel(se.get_data()[k],((double)(k+1)),1.0e-12,"data 2");
	  }
	}
	t.test_rel(se.progress(),((double)(i+1))/90.0,1.0e-12,"progress 2");
      } else {
	for(size_t k=0;k<9;k++) {
	  if (se.get_data()[k]>0.0) {
	    t.test_rel(se.get_data()[k],((double)(k*2))+1.5,1.0e-12,"data 3");
	  }
	}
	t.test_rel(se.progress(),((double)(i+1))/180.0,1.0e-12,"progress 3");
      }

      if (i==89) {
	se.reblock_avg_stats(3,avg,sd,avge,m_per_block);
	cout << m_per_block << " ";
	cout << avg << " " << sd << " " << avge << endl;
	t.test_gen(m_per_block==30,"mpb 1");
	t.test_rel(avg,5.0,1.0e-12,"avg 1");
      
	se.reblock_avg_stats(2,avg,sd,avge,m_per_block);
	cout << m_per_block << " ";
	cout << avg << " " << sd << " " << avge << endl;
	t.test_gen(m_per_block==40,"mpb 2");
	t.test_rel(avg,4.5,1.0e-12,"avg 2");
      }

    }

    se.reblock_avg_stats(2,avg,sd,avge,m_per_block);
    cout << m_per_block << " ";
    cout << avg << " " << sd << " " << avge << endl;
    t.test_gen(m_per_block==40,"mpb 3");
    t.test_rel(avg,4.5,1.0e-12,"avg 3");

    cout << endl;

  }

  // ------------------------------------------------------------------
  // Test expval_scalar with an even number of blocks
  // ------------------------------------------------------------------

  if (true) {

    expval_scalar se2(10,5);

    cout << "n_per_block=5: " << endl;
    cout << "x          i+1 ib ipb mb mpb prog       f "
	 << "avg        sd         avge         data" << endl;

    for(size_t i=0;i<100;i++) {
      double x=((double)(i/10))+1.0;

      // Add measurement and get info
      se2.add(x);
      se2.get_block_indices(i_block,i_per_block);

      // Output
      se2.current_avg_stats(avg,sd,avge,m_block,m_per_block);

      cout << x << " ";
      cout.width(3);
      cout << i+1 << " ";
      cout.width(2);
      cout << i_block << " ";
      cout.width(2);
      cout << i_per_block << " ";
      cout.width(2);
      cout << m_block << "   ";
      cout.width(2);
      cout << m_per_block << " ";
      cout << se2.progress() << " " << se2.finished() << " "
	   << avg << " " << sd << " " << avge << " : ";
      cout.unsetf(ios::scientific);
      //cout << se2.get_data() << endl;
      cout << endl;
      cout.setf(ios::scientific);

      if (i<=49) {
	for(size_t k=0;k<9;k++) {
	  if (se2.get_data()[k]>0.0) {
	    t.test_rel(se2.get_data()[k],((double)(k/2+1)),1.0e-12,"data 4");
	  }
	}
	t.test_rel(se2.progress(),((double)(i+1))/50.0,1.0e-12,"progress 4");
      } else {
	for(size_t k=0;k<9;k++) {
	  if (se2.get_data()[k]>0.0) {
	    t.test_rel(se2.get_data()[k],((double)(k))+1.0,1.0e-12,"data 5");
	  }
	}
	t.test_rel(se2.progress(),((double)(i+1))/100.0,1.0e-12,"progress 5");
      }
    
      if (i<50) {
	if (i>3) {
	  t.test_gen(m_per_block==5,"mpb 1");
	}
	t.test_gen((i+1)%5==i_per_block,"ipb 1");
      } else {
	t.test_gen(m_per_block==10,"mpb 2");
	t.test_gen((i+1)%10==i_per_block,"ipb 2");
      }

      if (i==49 || i==99) {
	t.test_gen(se2.finished()==true,"fin");
      } else {
	t.test_gen(se2.finished()==false,"fin");
      }

    }
    se2.reblock_avg_stats(3,avg,sd,avge,m_per_block);
    cout << m_per_block << " ";
    cout << avg << " " << sd << " " << avge << endl;
    t.test_gen(m_per_block==30,"mpb 1");
    t.test_rel(avg,5.0,1.0e-12,"avg 4");
    cout << endl;

  }

  // ------------------------------------------------------------------
  // Test expval_vector with an odd number of blocks
  // ------------------------------------------------------------------

  if (true) {

    expval_vector ve(10,9,5);
    ubvector vtmp(10);
    ubvector vavg(10), vsd(10), vavge(10);

    cout << "n_per_block=5: " << endl;
    cout << "x          i+1 ib ipb mb mpb prog       f "
	 << "avg        sd         avge         data" << endl;
  
    for(size_t i=0;i<100;i++) {
      for(size_t il=0;il<vtmp.size();il++) vtmp[il]=((double)(i/10))+1.0;
    
      // Add measurement and get info
      ve.add(vtmp);
      ve.get_block_indices(i_block,i_per_block);

      // Output
      ve.current_avg_stats(vavg,vsd,vavge,m_block,m_per_block);
    
      size_t ix=((size_t)(fabs(sin(i*o2scl_const::pi*100.0)*10.0)));

      {
	cout << vtmp[ix] << " ";
	cout.width(3);
	cout << i+1 << " ";
	cout.width(2);
	cout << i_block << " ";
	cout.width(2);
	cout << i_per_block << " ";
	cout.width(2);
	cout << m_block << "   ";
	cout.width(2);
	cout << m_per_block << " ";
	cout << ve.progress() << " " << ve.finished() << " "
	     << vavg[ix] << " " << vsd[ix] << " " << vavge[ix] << " ";
	cout.unsetf(ios::scientific);
	for(size_t k=0;k<9;k++) {
	  //cout << ve.get_data()(ix,k) << " ";
	}
	cout.setf(ios::scientific);
	cout << endl;
      }

      if (i==44 || i==89) {
	t.test_gen(ve.finished()==true,"fin");
      } else {
	t.test_gen(ve.finished()==false,"fin");
      }

      if (i<=44) {
	for(size_t k=0;k<9;k++) {
	  if (ve.get_data()(ix,k)>0.0) {
	    t.test_rel(ve.get_data()(ix,k),((double)(k/2+1)),
		       1.0e-12,"data 6");
	  }
	}
	t.test_rel(ve.progress(),((double)(i+1))/45.0,1.0e-12,"progress 6");
      } else if (i<=89) {
	for(size_t k=0;k<9;k++) {
	  if (ve.get_data()(ix,k)>0.0) {
	    t.test_rel(ve.get_data()(ix,k),((double)(k+1)),
		       1.0e-12,"data 7");
	  }
	}
	t.test_rel(ve.progress(),((double)(i+1))/90.0,1.0e-12,"progress 7");
      } else {
	for(size_t k=0;k<9;k++) {
	  if (ve.get_data()(ix,k)>0.0) {
	    t.test_rel(ve.get_data()(ix,k),((double)(k*2))+1.5,
		       1.0e-12,"data 8");
	  }
	}
	t.test_rel(ve.progress(),((double)(i+1))/180.0,1.0e-12,"progress 8");
      }

      if (i==89) {
	ve.reblock_avg_stats(3,vavg,vsd,vavge,m_per_block);
	cout << m_per_block << " ";
	cout << vavg[ix] << " " << vsd[ix] << " " << vavge[ix] << endl;
	t.test_gen(m_per_block==30,"mpb 1");
	t.test_rel(vavg[ix],5.0,1.0e-12,"avg 5");
      
	ve.reblock_avg_stats(2,vavg,vsd,vavge,m_per_block);
	cout << m_per_block << " ";
	cout << vavg[ix] << " " << vsd[ix] << " " << vavge[ix] << endl;
	t.test_gen(m_per_block==40,"mpb 2");
	t.test_rel(vavg[ix],4.5,1.0e-12,"avg 6");
      }

    }

    ve.reblock_avg_stats(2,vavg,vsd,vavge,m_per_block);
    cout << m_per_block << " ";
    cout << vavg[0] << " " << vsd[0] << " " << vavge[0] << endl;
    t.test_gen(m_per_block==40,"mpb 3");
    t.test_rel(vavg[0],4.5,1.0e-12,"avg 7");

    cout << endl;

  }

  // ------------------------------------------------------------------
  // Test expval_vector with an even number of blocks
  // ------------------------------------------------------------------

  if (true) {

    expval_vector ve2(10,10,5);
    ubvector vtmp(10);
    ubvector vavg(10), vsd(10), vavge(10);

    cout << "n_per_block=5: " << endl;
    cout << "x          i+1 ib ipb mb mpb prog       f "
	 << "avg        sd         avge         data" << endl;

    for(size_t i=0;i<100;i++) {
      for(size_t il=0;il<vtmp.size();il++) vtmp[il]=((double)(i/10))+1.0;

      // Add measurement and get info
      ve2.add(vtmp);
      ve2.get_block_indices(i_block,i_per_block);

      // Output
      ve2.current_avg_stats(vavg,vsd,vavge,m_block,m_per_block);

      size_t ix=((size_t)(fabs(sin(i*o2scl_const::pi*100.0)*10.0)));

      cout << vtmp[ix] << " ";
      cout.width(3);
      cout << i+1 << " ";
      cout.width(2);
      cout << i_block << " ";
      cout.width(2);
      cout << i_per_block << " ";
      cout.width(2);
      cout << m_block << "   ";
      cout.width(2);
      cout << m_per_block << " ";
      cout << ve2.progress() << " " << ve2.finished() << " "
	   << vavg[ix] << " " << vsd[ix] << " " << vavge[ix] << " : ";
      cout.unsetf(ios::scientific);
      for(size_t k=0;k<10;k++) {
	//cout << ve2.get_data()(ix,k) << " ";
      }
      cout.setf(ios::scientific);
      cout << endl;

      if (i<=49) {
	for(size_t k=0;k<9;k++) {
	  if (ve2.get_data()(ix,k)>0.0) {
	    t.test_rel(ve2.get_data()(ix,k),((double)(k/2+1)),
		       1.0e-12,"data 9");
	  }
	}
	t.test_rel(ve2.progress(),((double)(i+1))/50.0,
		   1.0e-12,"progress 9");
      } else {
	for(size_t k=0;k<9;k++) {
	  if (ve2.get_data()(ix,k)>0.0) {
	    t.test_rel(ve2.get_data()(ix,k),((double)(k))+1.0,
		       1.0e-12,"data 10");
	  }
	}
	t.test_rel(ve2.progress(),((double)(i+1))/100.0,
		   1.0e-12,"progress 10");
      }
    
      if (i<50) {
	if (i>3) {
	  t.test_gen(m_per_block==5,"mpb 1");
	}
	t.test_gen((i+1)%5==i_per_block,"ipb 1");
      } else {
	t.test_gen(m_per_block==10,"mpb 2");
	t.test_gen((i+1)%10==i_per_block,"ipb 2");
      }

      if (i==49 || i==99) {
	t.test_gen(ve2.finished()==true,"fin");
      } else {
	t.test_gen(ve2.finished()==false,"fin");
      }

    }
    ve2.reblock_avg_stats(3,vavg,vsd,vavge,m_per_block);
    cout << m_per_block << " ";
    cout << vavg[0] << " " << vsd[0] << " " << vavge[0] << endl;
    t.test_gen(m_per_block==30,"mpb 1");
    t.test_rel(vavg[0],5.0,1.0e-12,"avg 8");

    cout << endl;

  }

  // ------------------------------------------------------------------
  // Test expval_matrix with an odd number of blocks
  // ------------------------------------------------------------------

  if (true) {

    expval_matrix me(10,10,9,5);
    ubmatrix mtmp(10,10);
    ubmatrix mavg(10,10), msd(10,10), mavge(10,10);

    cout << "n_per_block=5: " << endl;
    cout << "x          i+1 ib ipb mb mpb prog       f "
	 << "avg        sd         avge         data" << endl;
  
    for(size_t i=0;i<100;i++) {

      for(size_t il=0;il<mtmp.size1();il++) {
	for(size_t jl=0;jl<mtmp.size2();jl++) {
	  mtmp(il,jl)=((double)(i/10))+1.0;
	}
      }

      // Add measurement and get info
      me.add(mtmp);
      me.get_block_indices(i_block,i_per_block);

      // Output
      me.current_avg_stats(mavg,msd,mavge,m_block,m_per_block);
    
      size_t ix=((size_t)(fabs(sin(i*o2scl_const::pi*100.0)*10.0)));
      size_t jx=((size_t)(fabs(sin((i+1)*o2scl_const::pi*100.0)*10.0)));

      {
	cout << mtmp(ix,jx) << " ";
	cout.width(3);
	cout << i+1 << " ";
	cout.width(2);
	cout << i_block << " ";
	cout.width(2);
	cout << i_per_block << " ";
	cout.width(2);
	cout << m_block << "   ";
	cout.width(2);
	cout << m_per_block << " ";
	cout << me.progress() << " " << me.finished() << " "
	     << mavg(ix,jx) << " " << msd(ix,jx) << " " 
	     << mavge(ix,jx) << " ";
	cout.unsetf(ios::scientific);
	for(size_t k=0;k<9;k++) {
	  //cout << me.get_data().get(ix,jx,k) << " ";
	}
	cout.setf(ios::scientific);
	cout << endl;
      }

      if (i==44 || i==89) {
	t.test_gen(me.finished()==true,"fin");
      } else {
	t.test_gen(me.finished()==false,"fin");
      }

      if (i<=44) {
	for(size_t k=0;k<9;k++) {
	  if (me.get_data().get(ix,jx,k)>0.0) {
	    t.test_rel(me.get_data().get(ix,jx,k),
		       ((double)(k/2+1)),1.0e-12,"data 11");
	  }
	}
	t.test_rel(me.progress(),((double)(i+1))/45.0,
		   1.0e-12,"progress 11");
      } else if (i<=89) {
	for(size_t k=0;k<9;k++) {
	  if (me.get_data().get(ix,jx,k)>0.0) {
	    t.test_rel(me.get_data().get(ix,jx,k),
		       ((double)(k+1)),1.0e-12,"data 12");
	  }
	}
	t.test_rel(me.progress(),((double)(i+1))/90.0,
		   1.0e-12,"progress 12");
      } else {
	for(size_t k=0;k<9;k++) {
	  if (me.get_data().get(ix,jx,k)>0.0) {
	    t.test_rel(me.get_data().get(ix,jx,k),
		       ((double)(k*2))+1.5,1.0e-12,"data 13");
	  }
	}
	t.test_rel(me.progress(),((double)(i+1))/180.0,
		   1.0e-12,"progress 13");
      }

      if (i==89) {
	me.reblock_avg_stats(3,mavg,msd,mavge,m_per_block);
	cout << m_per_block << " ";
	cout << mavg(ix,jx) << " " << msd(ix,jx) << " " 
	     << mavge(ix,jx) << endl;
	t.test_gen(m_per_block==30,"mpb 9");
	t.test_rel(mavg(ix,jx),5.0,1.0e-12,"avg 9");
      
	me.reblock_avg_stats(2,mavg,msd,mavge,m_per_block);
	cout << m_per_block << " ";
	cout << mavg(ix,jx) << " " << msd(ix,jx) << " " 
	     << mavge(ix,jx) << endl;
	t.test_gen(m_per_block==40,"mpb 10");
	t.test_rel(mavg(ix,jx),4.5,1.0e-12,"avg 10");
      }

    }

    me.reblock_avg_stats(2,mavg,msd,mavge,m_per_block);
    cout << m_per_block << " ";
    cout << mavg(0,0) << " " << msd(0,0) << " " << mavge(0,0) << endl;
    t.test_gen(m_per_block==40,"mpb 3");
    t.test_rel(mavg(0,0),4.5,1.0e-12,"avg 11");

    cout << endl;

  }

  // ------------------------------------------------------------------
  // Test expval_matrix with an even number of blocks
  // ------------------------------------------------------------------

  if (true) {

    expval_matrix me2(10,10,10,5);
    ubmatrix mtmp(10,10);
    ubmatrix mavg(10,10), msd(10,10), mavge(10,10);

    cout << "n_per_block=5: " << endl;
    cout << "x          i+1 ib ipb mb mpb prog       f "
	 << "avg        sd         avge         data" << endl;

    for(size_t i=0;i<100;i++) {

      for(size_t il=0;il<mtmp.size1();il++) {
	for(size_t jl=0;jl<mtmp.size2();jl++) {
	  mtmp(il,jl)=((double)(i/10))+1.0;
	}
      }

      // Add measurement and get info
      me2.add(mtmp);
      me2.get_block_indices(i_block,i_per_block);

      // Output
      me2.current_avg_stats(mavg,msd,mavge,m_block,m_per_block);

      size_t ix=((size_t)(fabs(sin(i*o2scl_const::pi*100.0)*10.0)));
      size_t jx=((size_t)(fabs(sin((i+1)*o2scl_const::pi*100.0)*10.0)));

      cout << mtmp(ix,jx) << " ";
      cout.width(3);
      cout << i+1 << " ";
      cout.width(2);
      cout << i_block << " ";
      cout.width(2);
      cout << i_per_block << " ";
      cout.width(2);
      cout << m_block << "   ";
      cout.width(2);
      cout << m_per_block << " ";
      cout << me2.progress() << " " << me2.finished() << " "
	   << mavg(ix,jx) << " " << msd(ix,jx) << " " 
	   << mavge(ix,jx) << " : ";
      cout.unsetf(ios::scientific);
      for(size_t k=0;k<10;k++) {
	//cout << me2.get_data().get(ix,jx,k) << " ";
      }
      cout.setf(ios::scientific);
      cout << endl;

      if (i<=49) {
	for(size_t k=0;k<9;k++) {
	  if (me2.get_data().get(ix,jx,k)>0.0) {
	    t.test_rel(me2.get_data().get(ix,jx,k),
		       ((double)(k/2+1)),1.0e-12,"data 14");
	  }
	}
	t.test_rel(me2.progress(),((double)(i+1))/50.0,
		   1.0e-12,"progress 14");
      } else {
	for(size_t k=0;k<9;k++) {
	  if (me2.get_data().get(ix,jx,k)>0.0) {
	    t.test_rel(me2.get_data().get(ix,jx,k),
		       ((double)(k))+1.0,1.0e-12,"data 15");
	  }
	}
	t.test_rel(me2.progress(),((double)(i+1))/100.0,
		   1.0e-12,"progress 15");
      }
    
      if (i<50) {
	if (i>3) {
	  t.test_gen(m_per_block==5,"mpb 11");
	}
	t.test_gen((i+1)%5==i_per_block,"ipb 11");
      } else {
	t.test_gen(m_per_block==10,"mpb 12");
	t.test_gen((i+1)%10==i_per_block,"ipb 12");
      }

      if (i==49 || i==99) {
	t.test_gen(me2.finished()==true,"fin");
      } else {
	t.test_gen(me2.finished()==false,"fin");
      }

    }
    me2.reblock_avg_stats(3,mavg,msd,mavge,m_per_block);
    cout << m_per_block << " ";
    cout << mavg(0,0) << " " << msd(0,0) << " " << mavge(0,0) << endl;
    t.test_gen(m_per_block==30,"mpb 1");
    t.test_rel(mavg(0,0),5.0,1.0e-12,"avg 13");

    cout << endl;

  }

  // ------------------------------------------------------------------
  // Done
  // ------------------------------------------------------------------

  t.report();

  return 0;
}
