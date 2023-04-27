/**
 * ExchCXX Copyright (c) 2020-2022, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * (1) Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 
 * (2) Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * 
 * (3) Neither the name of the University of California, Lawrence Berkeley
 * National Laboratory, U.S. Dept. of Energy nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 * 
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * 
 * You are under no obligation whatsoever to provide any bug fixes, patches,
 * or upgrades to the features, functionality or performance of the source
 * code ("Enhancements") to anyone; however, if you choose to make your
 * Enhancements available either publicly, or directly to Lawrence Berkeley
 * National Laboratory, without imposing a separate written license agreement
 * for such Enhancements, then you hereby grant the following license: a
 * non-exclusive, royalty-free perpetual license to install, use, modify,
 * prepare derivative works, incorporate into other computer software,
 * distribute, and sublicense such enhancements or derivative works thereof,
 * in binary and source code form.
 */

#include "reference_values.hpp"
#include <exchcxx/exchcxx.hpp>

std::vector<double> rho = {0.1, 0.2, 0.3, 0.4, 0.5};
std::vector<double> sigma = {0.2, 0.3, 0.4, 0.5, 0.6};

std::vector<double> rho_polarized =
  {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
std::vector<double> sigma_polarized =
  {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
   1.1, 1.2, 1.3, 1.5, 1.5 };

std::vector<double> exc_xc_lda_x_ref_unp = {
  -0.342808612301, -0.431911786723,
  -0.494415573788, -0.544174751790,
  -0.586194481348,
};

std::vector<double> vxc_xc_lda_x_ref_unp = {
  -0.457078149734, -0.575882382297,
  -0.659220765051, -0.725566335720,
  -0.781592641797
};

std::vector<double> exc_xc_lda_x_ref_pol = {
  -0.506753763434, -0.658748952120,
  -0.763800785778, -0.846274084184,
  -0.915314307811
};

std::vector<double> vxc_xc_lda_x_ref_pol = {
-0.575882382297, -0.725566335720,
-0.830566118415, -0.914156299468,
-0.984745021843, -1.046447735921,
-1.101623366705, -1.151764764594,
-1.197883627397, -1.240700981799
};



std::vector<double> exc_xc_gga_c_lyp_ref_unp = {
  -0.007040306272, -0.031424640440,
  -0.037479119388, -0.040429224120,
  -0.042290563929
};

std::vector<double> vrho_xc_gga_c_lyp_ref_unp = {
  -0.081854247031, -0.055198496086,
  -0.051617025994, -0.050995654065,
  -0.051084686930
};

std::vector<double> vsigma_xc_gga_c_lyp_ref_unp = {
  0.013598460611,
  0.004629650473,
  0.002429957976,
  0.001529632674,
  0.001065244937
};

std::vector<double> exc_xc_gga_c_lyp_ref_pol = {
  -0.031543975366, -0.043113613690,
  -0.046604883008, -0.048519647105,
  -0.049807583631
};

std::vector<double> vrho_xc_gga_c_lyp_ref_pol = {
  -0.089983823779, -0.035745759262,
  -0.062361000975, -0.045947114249,
  -0.059003605615, -0.049400798274,
  -0.058191535482, -0.051370405717,
  -0.058012961950, -0.052706534799
};

std::vector<double> vsigma_xc_gga_c_lyp_ref_pol = {
0.008447669161 , 0.006467154082,
-0.000638497084, 0.001421914705,
0.001031651601 , 0.000257537600,
0.000581699649 , 0.000435910598,
0.000202132738 , 0.000321427269,
0.000246773907 , 0.000146744820,
0.000206495563 , 0.000160996240,
0.000110118944
};


template <typename T, typename OutIt>
void copy_iterable( const T& src, OutIt&& dest ) {
  std::copy( src.begin(), src.end(), dest );
}


lda_reference load_lda_reference_values(ExchCXX::Kernel k, ExchCXX::Spin p) {

  using namespace ExchCXX; 

  lda_reference ref_vals;

  if( p == Spin::Unpolarized ) {

    ref_vals.npts = rho.size();
    copy_iterable( rho, std::back_inserter(ref_vals.rho) );

    switch(k) {
      case Kernel::SlaterExchange:
        copy_iterable( exc_xc_lda_x_ref_unp, std::back_inserter(ref_vals.exc) );
        copy_iterable( vxc_xc_lda_x_ref_unp, std::back_inserter(ref_vals.vrho) );
        break;
      default: 
        throw std::runtime_error("No Reference Values for Specified Kernel");
    }

  } else {

    copy_iterable( rho_polarized, std::back_inserter(ref_vals.rho) );
    ref_vals.npts = rho_polarized.size() / 2;

    switch(k) {
      case Kernel::SlaterExchange:
        copy_iterable( exc_xc_lda_x_ref_pol, std::back_inserter(ref_vals.exc) );
        copy_iterable( vxc_xc_lda_x_ref_pol, std::back_inserter(ref_vals.vrho) );
        break;
      default: 
        throw std::runtime_error("No Reference Values for Specified Kernel");
    }

  }

  return ref_vals;
}




gga_reference load_gga_reference_values(ExchCXX::Kernel k, ExchCXX::Spin p) {

  using namespace ExchCXX; 

  gga_reference ref_vals;

  if( p == Spin::Unpolarized ) {

    ref_vals.npts = rho.size();
    copy_iterable( rho, std::back_inserter(ref_vals.rho) );
    copy_iterable( sigma, std::back_inserter(ref_vals.sigma) );

    switch(k) {
      case Kernel::LYP:
        copy_iterable(exc_xc_gga_c_lyp_ref_unp, std::back_inserter(ref_vals.exc) );
        copy_iterable(vrho_xc_gga_c_lyp_ref_unp, std::back_inserter(ref_vals.vrho) );
        copy_iterable(vsigma_xc_gga_c_lyp_ref_unp, std::back_inserter(ref_vals.vsigma) );
        break;
      default: 
        throw std::runtime_error("No Reference Values for Specified Kernel");
    }

  } else {

    ref_vals.npts = rho_polarized.size() / 2;
    copy_iterable(  rho_polarized, std::back_inserter(ref_vals.rho) );
    copy_iterable(  sigma_polarized, std::back_inserter(ref_vals.sigma) );

    switch(k) {
      case Kernel::LYP:
        copy_iterable(exc_xc_gga_c_lyp_ref_pol, std::back_inserter(ref_vals.exc) );
        copy_iterable(vrho_xc_gga_c_lyp_ref_pol, std::back_inserter(ref_vals.vrho) );
        copy_iterable(vsigma_xc_gga_c_lyp_ref_pol, std::back_inserter(ref_vals.vsigma) );
        break;
      default: 
        throw std::runtime_error("No Reference Values for Specified Kernel");
    }

  }

  return ref_vals;
}



lda_reference gen_lda_reference_values(ExchCXX::Backend backend, 
  ExchCXX::Kernel kern, ExchCXX::Spin polar) {
  
  using namespace ExchCXX;
  lda_reference ref_vals;
  std::tie( ref_vals.npts, ref_vals.rho ) = load_reference_density( polar );

  XCKernel func( backend, kern, polar );

  ref_vals.exc.resize( func.exc_buffer_len( ref_vals.npts ) );
  ref_vals.vrho.resize( func.vrho_buffer_len( ref_vals.npts ) );

  func.eval_exc_vxc( ref_vals.npts, ref_vals.rho.data(), ref_vals.exc.data(),
    ref_vals.vrho.data() );

  return ref_vals;

}

gga_reference gen_gga_reference_values(ExchCXX::Backend backend, 
  ExchCXX::Kernel kern, ExchCXX::Spin polar) {
  
  using namespace ExchCXX;
  gga_reference ref_vals;
  std::tie( ref_vals.npts, ref_vals.rho )   = load_reference_density( polar );
  std::tie( ref_vals.npts, ref_vals.sigma ) = load_reference_sigma( polar );

  XCKernel func( backend, kern, polar );

  ref_vals.exc.resize( func.exc_buffer_len( ref_vals.npts ) );
  ref_vals.vrho.resize( func.vrho_buffer_len( ref_vals.npts ) );
  ref_vals.vsigma.resize( func.vsigma_buffer_len( ref_vals.npts ) );

  func.eval_exc_vxc( ref_vals.npts, ref_vals.rho.data(), ref_vals.sigma.data(), 
    ref_vals.exc.data(), ref_vals.vrho.data(), ref_vals.vsigma.data() );

  return ref_vals;

}






std::pair<int,std::vector<double>> load_reference_density(ExchCXX::Spin s) {

  std::vector<double> ref_rho;
  int npts;

  if( s == ExchCXX::Spin::Polarized ) {
    npts = rho_polarized.size() / 2;
    copy_iterable( rho_polarized, std::back_inserter(ref_rho) );
  } else {
    npts = rho.size();
    copy_iterable( rho, std::back_inserter(ref_rho) );
  }
  return std::make_pair(npts, ref_rho);

}


std::pair<int,std::vector<double>> load_reference_sigma(ExchCXX::Spin s) { 

  std::vector<double> ref_sigma;
  int npts;

  if( s == ExchCXX::Spin::Polarized ) {
    npts = sigma_polarized.size() / 3;
    copy_iterable( sigma_polarized, std::back_inserter(ref_sigma) );
  } else {
    npts = sigma.size();
    copy_iterable( sigma, std::back_inserter(ref_sigma) );
  }
  return std::make_pair(npts, ref_sigma);

}


double load_reference_exx( ExchCXX::Kernel kern ) {

  using namespace ExchCXX;

  double exx;
  switch(kern) {
    case Kernel::B3LYP: exx = 0.2;  break;
    case Kernel::PBE0:  exx = 0.25; break;
    default:            exx = 0;
  }

  return exx;
}
