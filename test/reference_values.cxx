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
std::vector<double> lapl  = {0.3, 0.4, 0.5, 0.6, 0.7};
std::vector<double> tau   = {0.2, 0.3, 0.4, 0.5, 0.6};

std::vector<double> rho_polarized =
  {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
std::vector<double> sigma_polarized =
  {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
   1.1, 1.2, 1.3, 1.4, 1.5 };
std::vector<double> lapl_polarized =
  {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
std::vector<double> tau_polarized =
  {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

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





std::vector<double> exc_xc_mgga_c_scan_ref_unp = {
-0.396376974817, -0.471540080123,
-0.533396348869, -0.591369866023,
-0.642807787274
};

std::vector<double> vrho_xc_mgga_c_scan_ref_unp = {
-0.460152650566, -0.623118535281,
-0.764470500527, -0.862470512705,
-0.936094993961
};

std::vector<double> vsigma_xc_mgga_c_scan_ref_unp = {
-0.049357276976, -0.056332014377,
-0.040207846026, -0.028050533162,
-0.020140871184
};

std::vector<double> vtau_xc_mgga_c_scan_ref_unp = {
0.042672319773, 0.087890594383,
0.088306469557, 0.080389984903,
0.071734366025
};

std::vector<double> exc_xc_mgga_c_scan_ref_pol = {
-0.592191935274, -0.756007107013,
-0.875625350878, -0.972027141138,
-1.053486673118
};

std::vector<double> vrho_xc_mgga_c_scan_ref_pol = {
-0.642368591061, -0.822909425557,
-0.958799773865, -1.063088704863,
-1.153559630345, -1.228156943455,
-1.296083819945, -1.355697783203,
-1.411287187583, -1.461796825163
};

std::vector<double> vsigma_xc_mgga_c_scan_ref_pol = {
-0.044495462460, 0.000000000000,
-0.020457936958, -0.018239002704,
 0.000000000000, -0.012892377595,
-0.010154834731, 0.000000000000,
-0.007876478998, -0.006455830014,
 0.000000000000, -0.005319202936,
-0.004524387110, 0.000000000000,
-0.003875648399
};

std::vector<double> vtau_xc_mgga_c_scan_ref_pol = {
0.038229318445, 0.034287395774,
0.042250078467, 0.039339723043,
0.038151864869, 0.035489486994,
0.033824550829, 0.031903228786,
0.030524916152, 0.029113290061
};




std::vector<double> exc_xc_mgga_c_r2scanl_ref_unp = {
-0.394823861318, -0.491171995952,
-0.528380290237, -0.565443056965,
-0.600480870171
};

std::vector<double> vrho_xc_mgga_c_r2scanl_ref_unp = {
-0.546609816979, -0.546087438942,
-0.621955980174, -0.694599669881,
-0.762595835512
};

std::vector<double> vsigma_xc_mgga_c_r2scanl_ref_unp = {
 0.003783375354, -0.036415775717,
-0.033904222674, -0.027285116327,
-0.020495017477
};

std::vector<double> vlapl_xc_mgga_c_r2scanl_ref_unp = {
0.000000000000, 0.011056364163,
0.013679018572, 0.012650392673,
0.011802313376
};

std::vector<double> vtau_xc_mgga_c_r2scanl_ref_unp = {
0.000000000000, 0.000000000000,
0.000000000000, 0.000000000000,
0.000000000000
};

std::vector<double> exc_xc_mgga_c_r2scanl_ref_pol = {
-0.575725446547, -0.692588892216,
-0.784198364004, -0.862584626416,
-0.930239925125
};

std::vector<double> vrho_xc_mgga_c_r2scanl_ref_pol = {
-0.681495044795, -0.688746989901,
-0.790625972596, -0.882861971126,
-0.971394010990, -1.040800712793,
-1.103937287811, -1.157509381409,
-1.207429868663, -1.252191434665
};

std::vector<double> vsigma_xc_mgga_c_r2scanl_ref_pol = {
 0.003296172548,  0.000000000000,
-0.035982064895, -0.031040505262,
 0.000000000000, -0.022323946859,
-0.015645099288,  0.000000000000,
-0.011675746617, -0.008965945413,
 0.000000000000, -0.007315682423,
-0.006084573039,  0.000000000000,
-0.005226777712
};

std::vector<double> vlapl_xc_mgga_c_r2scanl_ref_pol = {
0.000011563079, 0.010540407693,
0.011426280664, 0.010546991623,
0.009746880821, 0.009123409167,
0.008581982177, 0.008176059558,
0.007814849266, 0.007529117778
};

std::vector<double> vtau_xc_mgga_c_r2scanl_ref_pol = {
0.000000000000, 0.000000000000,
0.000000000000, 0.000000000000,
0.000000000000, 0.000000000000,
0.000000000000, 0.000000000000,
0.000000000000, 0.000000000000
};


template <typename T, typename OutIt>
void copy_iterable( const T& src, OutIt&& dest ) {
  std::copy( src.begin(), src.end(), dest );
}


lda_reference load_lda_reference_values(ExchCXX::Kernel k, ExchCXX::Spin p) {

  using namespace ExchCXX; 

  lda_reference ref_vals;

  if( p == Spin::Unpolarized ) {

    copy_iterable( rho, std::back_inserter(ref_vals.rho) );
    ref_vals.npts = rho.size();

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

mgga_reference load_mgga_reference_values(ExchCXX::Kernel k, ExchCXX::Spin p, bool need_lap) {

  using namespace ExchCXX; 

  mgga_reference ref_vals;

  if( p == Spin::Unpolarized ) {

    copy_iterable( rho, std::back_inserter(ref_vals.rho) );
    copy_iterable( sigma, std::back_inserter(ref_vals.sigma) );
    copy_iterable( tau, std::back_inserter(ref_vals.tau) );
    if(need_lap) copy_iterable( lapl, std::back_inserter(ref_vals.lapl) );
    ref_vals.npts = rho.size();

    switch(k) {
      case Kernel::SCAN_X:
        copy_iterable(exc_xc_mgga_c_scan_ref_unp, std::back_inserter(ref_vals.exc) );
        copy_iterable(vrho_xc_mgga_c_scan_ref_unp, std::back_inserter(ref_vals.vrho) );
        copy_iterable(vsigma_xc_mgga_c_scan_ref_unp, std::back_inserter(ref_vals.vsigma) );
        copy_iterable(vtau_xc_mgga_c_scan_ref_unp, std::back_inserter(ref_vals.vtau) );
        break;
      case Kernel::R2SCANL_X:
        copy_iterable(exc_xc_mgga_c_r2scanl_ref_unp, std::back_inserter(ref_vals.exc) );
        copy_iterable(vrho_xc_mgga_c_r2scanl_ref_unp, std::back_inserter(ref_vals.vrho) );
        copy_iterable(vsigma_xc_mgga_c_r2scanl_ref_unp, std::back_inserter(ref_vals.vsigma) );
        copy_iterable(vlapl_xc_mgga_c_r2scanl_ref_unp, std::back_inserter(ref_vals.vlapl) );
        copy_iterable(vtau_xc_mgga_c_r2scanl_ref_unp, std::back_inserter(ref_vals.vtau) );
        break;
      default: 
        throw std::runtime_error("No Reference Values for Specified Kernel");
    }

  } else {

    copy_iterable( rho_polarized, std::back_inserter(ref_vals.rho) );
    copy_iterable( sigma_polarized, std::back_inserter(ref_vals.sigma) );
    copy_iterable( tau_polarized, std::back_inserter(ref_vals.tau) );
    if(need_lap) copy_iterable( lapl_polarized, std::back_inserter(ref_vals.lapl) );
    ref_vals.npts = rho_polarized.size() / 2;

    switch(k) {
      case Kernel::SCAN_X:
        copy_iterable(exc_xc_mgga_c_scan_ref_pol, std::back_inserter(ref_vals.exc) );
        copy_iterable(vrho_xc_mgga_c_scan_ref_pol, std::back_inserter(ref_vals.vrho) );
        copy_iterable(vsigma_xc_mgga_c_scan_ref_pol, std::back_inserter(ref_vals.vsigma) );
        copy_iterable(vtau_xc_mgga_c_scan_ref_pol, std::back_inserter(ref_vals.vtau) );
        break;
      case Kernel::R2SCANL_X:
        copy_iterable(exc_xc_mgga_c_r2scanl_ref_pol, std::back_inserter(ref_vals.exc) );
        copy_iterable(vrho_xc_mgga_c_r2scanl_ref_pol, std::back_inserter(ref_vals.vrho) );
        copy_iterable(vsigma_xc_mgga_c_r2scanl_ref_pol, std::back_inserter(ref_vals.vsigma) );
        copy_iterable(vlapl_xc_mgga_c_r2scanl_ref_pol, std::back_inserter(ref_vals.vlapl) );
        copy_iterable(vtau_xc_mgga_c_r2scanl_ref_pol, std::back_inserter(ref_vals.vtau) );
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

mgga_reference gen_mgga_reference_values(ExchCXX::Backend backend, 
  ExchCXX::Kernel kern, ExchCXX::Spin polar) {
  
  using namespace ExchCXX;
  XCKernel func( backend, kern, polar );

  mgga_reference ref_vals;
  std::tie( ref_vals.npts, ref_vals.rho )   = load_reference_density( polar );
  std::tie( ref_vals.npts, ref_vals.sigma ) = load_reference_sigma( polar );
  std::tie( ref_vals.npts, ref_vals.tau )   = load_reference_tau( polar );
  if(func.needs_laplacian()) {
    std::tie( ref_vals.npts, ref_vals.lapl )  = load_reference_lapl( polar );
  }


  ref_vals.exc.resize( func.exc_buffer_len( ref_vals.npts ) );
  ref_vals.vrho.resize( func.vrho_buffer_len( ref_vals.npts ) );
  ref_vals.vsigma.resize( func.vsigma_buffer_len( ref_vals.npts ) );
  ref_vals.vlapl.resize( func.vlapl_buffer_len( ref_vals.npts ) );
  ref_vals.vtau.resize( func.vtau_buffer_len( ref_vals.npts ) );

  func.eval_exc_vxc( ref_vals.npts, ref_vals.rho.data(), ref_vals.sigma.data(), 
    ref_vals.lapl.data(), ref_vals.tau.data(), ref_vals.exc.data(), 
    ref_vals.vrho.data(), ref_vals.vsigma.data(), ref_vals.vlapl.data(),
    ref_vals.vtau.data() );

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

std::pair<int,std::vector<double>> load_reference_lapl(ExchCXX::Spin s) {

  std::vector<double> ref_lapl;
  int npts;

  if( s == ExchCXX::Spin::Polarized ) {
    npts = lapl_polarized.size() / 2;
    copy_iterable( lapl_polarized, std::back_inserter(ref_lapl) );
  } else {
    npts = lapl.size();
    copy_iterable( lapl, std::back_inserter(ref_lapl) );
  }
  return std::make_pair(npts, ref_lapl);

}

std::pair<int,std::vector<double>> load_reference_tau(ExchCXX::Spin s) {

  std::vector<double> ref_tau;
  int npts;

  if( s == ExchCXX::Spin::Polarized ) {
    npts = tau_polarized.size() / 2;
    copy_iterable( tau_polarized, std::back_inserter(ref_tau) );
  } else {
    npts = tau.size();
    copy_iterable( tau, std::back_inserter(ref_tau) );
  }
  return std::make_pair(npts, ref_tau);

}


double load_reference_exx( ExchCXX::Kernel kern ) {

  using namespace ExchCXX;

  double exx;
  switch(kern) {
    case Kernel::B3LYP:    exx = 0.2;  break;
    case Kernel::PBE0:     exx = 0.25; break;
    case Kernel::M062X_X:  exx = 0.54; break;
    default:               exx = 0;
  }

  return exx;
}
