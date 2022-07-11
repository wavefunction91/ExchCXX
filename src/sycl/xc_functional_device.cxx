#include <exchcxx/xc_functional.hpp>
#include <string>

template <typename T> class scal_device_tag;
template <typename T> class add_scal_device_tag;


void scal_device( const int N, const double fact, const double* X_device, double* Y_device, sycl::queue* queue ) {
  queue->parallel_for<scal_device_tag<double>>( sycl::range<1>(N),
    [=]( sycl::id<1> idx ) { Y_device[idx] = fact * X_device[idx]; });
}

void add_scal_device( const int N, const double fact, const double* X_device, double* Y_device, sycl::queue* queue ) {
  queue->parallel_for<add_scal_device_tag<double>>( sycl::range<1>(N),
    [=]( sycl::id<1> idx ) { Y_device[idx] += fact * X_device[idx]; });
}


template <typename T = double>
T* safe_sycl_malloc( size_t N, sycl::queue* queue ) {
  return sycl::malloc_device<T>( N, *queue );
}

template <typename T>
void safe_zero( size_t len, T* ptr, sycl::queue* queue ) {
  queue->memset( ptr, 0, len*sizeof(T) );
}

namespace ExchCXX {






LDA_EXC_GENERATOR_DEVICE( XCFunctional::eval_exc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT LDA",  is_lda() );

  size_t len_exc_buffer = exc_buffer_len( N );

  double* eps_scr = nullptr;
  if( kernels_.size() > 1 and not supports_inc_interface() )
    eps_scr = safe_sycl_malloc( len_exc_buffer, queue );

  safe_zero( len_exc_buffer, eps, queue );


  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      kernels_[i].second.eval_exc_inc_device(
        kernels_[i].first, N, rho, eps, queue
      );

    } else {

      double* eps_eval = i ? eps_scr : eps;
      kernels_[i].second.eval_exc_device(N, rho, eps_eval, queue);

      if( i )
        add_scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, queue );
      else
        scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, queue );

    }

  }

  if( eps_scr ) sycl::free( eps_scr, *queue );

}


LDA_EXC_VXC_GENERATOR_DEVICE( XCFunctional::eval_exc_vxc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT LDA",  is_lda() );

  size_t len_exc_buffer = exc_buffer_len( N );
  size_t len_vxc_buffer = vrho_buffer_len( N );

  double* eps_scr(nullptr), *vxc_scr(nullptr);
  if( kernels_.size() > 1 and not supports_inc_interface() ) {
    eps_scr = safe_sycl_malloc( len_exc_buffer, queue );
    vxc_scr = safe_sycl_malloc( len_vxc_buffer, queue );
  }

  safe_zero( len_exc_buffer, eps, queue );
  safe_zero( len_vxc_buffer, vxc, queue );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      kernels_[i].second.eval_exc_vxc_inc_device(
        kernels_[i].first, N, rho, eps, vxc, queue
      );

    } else {

      double* eps_eval = i ? eps_scr : eps;
      double* vxc_eval = i ? vxc_scr : vxc;
      kernels_[i].second.eval_exc_vxc_device(N, rho, eps_eval, vxc_eval, queue);

      if( i ) {

        add_scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, queue );
        add_scal_device( len_vxc_buffer, kernels_[i].first, vxc_eval, vxc, queue );

      } else {

        scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, queue );
        scal_device( len_vxc_buffer, kernels_[i].first, vxc_eval, vxc, queue );

      }

    }

  }

  if( eps_scr ) sycl::free( eps_scr, *queue );
  if( vxc_scr ) sycl::free( vxc_scr, *queue );

}



// GGA Interfaces

GGA_EXC_GENERATOR_DEVICE( XCFunctional::eval_exc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT GGA",  is_gga() );

  size_t len_exc_buffer = exc_buffer_len( N );

  double* eps_scr = nullptr;
  if( kernels_.size() > 1 and not supports_inc_interface() )
    eps_scr = safe_sycl_malloc( len_exc_buffer, queue );

  safe_zero( len_exc_buffer, eps, queue );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_inc_device(
          kernels_[i].first, N, rho, sigma, eps, queue
        );
      else
        kernels_[i].second.eval_exc_inc_device(
          kernels_[i].first, N, rho, eps, queue
        );

    } else {

      double* eps_eval = i ? eps_scr : eps;

      if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_device(N, rho, sigma, eps_eval, queue);
      else
        kernels_[i].second.eval_exc_device(N, rho, eps_eval, queue);

      if( i )
        add_scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, queue );
      else
        scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, queue );

    }
  }

  if( eps_scr ) sycl::free( eps_scr, *queue );

}


GGA_EXC_VXC_GENERATOR_DEVICE( XCFunctional::eval_exc_vxc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT GGA",  is_gga() );

  size_t len_exc_buffer    = exc_buffer_len(N);
  size_t len_vrho_buffer   = vrho_buffer_len(N);
  size_t len_vsigma_buffer = vsigma_buffer_len(N);

  double* eps_scr(nullptr), *vrho_scr(nullptr), *vsigma_scr(nullptr);
  if( kernels_.size() > 1 and not supports_inc_interface() ) {
    eps_scr    = safe_sycl_malloc( len_exc_buffer, queue );
    vrho_scr   = safe_sycl_malloc( len_vrho_buffer, queue );
    vsigma_scr = safe_sycl_malloc( len_vsigma_buffer, queue );
  }

  safe_zero( len_exc_buffer,    eps,    queue );
  safe_zero( len_vrho_buffer,   vrho,   queue );
  safe_zero( len_vsigma_buffer, vsigma, queue );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_vxc_inc_device(
          kernels_[i].first, N, rho, sigma, eps, vrho,
          vsigma, queue
        );
      else
        kernels_[i].second.eval_exc_vxc_inc_device(
          kernels_[i].first, N, rho, eps, vrho, queue
        );

    } else {

      double* eps_eval    = i ? eps_scr    : eps;
      double* vrho_eval   = i ? vrho_scr   : vrho;
      double* vsigma_eval = i ? vsigma_scr : vsigma;

      if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_vxc_device(N, rho, sigma, eps_eval, vrho_eval,
          vsigma_eval, queue );
      else
        kernels_[i].second.eval_exc_vxc_device(N, rho, eps_eval, vrho_eval, queue);

      if( i ) {

        add_scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, queue );
        add_scal_device( len_vrho_buffer, kernels_[i].first, vrho_eval, vrho, queue);
        if( kernels_[i].second.is_gga() )
          add_scal_device( len_vsigma_buffer, kernels_[i].first, vsigma_eval, vsigma, queue );

      } else {

        scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, queue );
        scal_device( len_vrho_buffer, kernels_[i].first, vrho_eval, vrho, queue );
        if( kernels_[i].second.is_gga() )
          scal_device( len_vsigma_buffer, kernels_[i].first, vsigma_eval, vsigma, queue );

      }

    }
  }

  if( eps_scr ) sycl::free( eps_scr, *queue );
  if( vrho_scr ) sycl::free( vrho_scr, *queue );
  if( vsigma_scr ) sycl::free( vsigma_scr, *queue );
}




// mGGA Interfaces

MGGA_EXC_GENERATOR_DEVICE( XCFunctional::eval_exc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT MGGA",  is_mgga() );

  size_t len_exc_buffer = exc_buffer_len( N );

  double* eps_scr = nullptr;
  if( kernels_.size() > 1 and not supports_inc_interface() )
    eps_scr = safe_sycl_malloc( len_exc_buffer, queue );

  safe_zero( len_exc_buffer,    eps,    queue );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      if( kernels_[i].second.is_mgga() )
        kernels_[i].second.eval_exc_inc_device(
          kernels_[i].first, N, rho, sigma, lapl, tau, eps, queue
        );
      else if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_inc_device(
          kernels_[i].first, N, rho, sigma, eps, queue
        );
      else
        kernels_[i].second.eval_exc_inc_device(
          kernels_[i].first, N, rho, eps, queue
        );

    } else {

      double* eps_eval = i ? eps_scr : eps;

      if( kernels_[i].second.is_mgga() )
        kernels_[i].second.eval_exc_device(N, rho, sigma, lapl, tau, eps_eval, queue);
      else if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_device(N, rho, sigma, eps_eval, queue);
      else
        kernels_[i].second.eval_exc_device(N, rho, eps_eval, queue);

      if( i )
        add_scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, queue );
      else
        scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, queue );

    }
  }

  if( eps_scr ) sycl::free( eps_scr, *queue );

}


MGGA_EXC_VXC_GENERATOR_DEVICE( XCFunctional::eval_exc_vxc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT GGA",  is_gga() );

  size_t len_exc_buffer    = exc_buffer_len(N);
  size_t len_vrho_buffer   = vrho_buffer_len(N);
  size_t len_vsigma_buffer = vsigma_buffer_len(N);
  size_t len_vlapl_buffer   = vlapl_buffer_len(N);
  size_t len_vtau_buffer   = vtau_buffer_len(N);

  double* eps_scr(nullptr), *vrho_scr(nullptr), *vsigma_scr(nullptr),
    *vlapl_scr(nullptr), *vtau_scr(nullptr);
  if( kernels_.size() > 1 and not supports_inc_interface() ) {
    eps_scr    = safe_sycl_malloc( len_exc_buffer, queue );
    vrho_scr   = safe_sycl_malloc( len_vrho_buffer, queue );
    vsigma_scr = safe_sycl_malloc( len_vsigma_buffer, queue );
    vlapl_scr  = safe_sycl_malloc( len_vlapl_buffer, queue );
    vtau_scr   = safe_sycl_malloc( len_vtau_buffer, queue );
  }

  safe_zero( len_exc_buffer, eps, queue );
  safe_zero( len_vrho_buffer, vrho, queue );
  safe_zero( len_vsigma_buffer, vsigma, queue );
  safe_zero( len_vlapl_buffer, vlapl, queue );
  safe_zero( len_vtau_buffer, vtau, queue );


  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      if( kernels_[i].second.is_mgga() )
        kernels_[i].second.eval_exc_vxc_inc_device(
          kernels_[i].first, N, rho, sigma, lapl, tau, eps,
          vrho, vsigma, vlapl, vtau, queue
        );
      else if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_vxc_inc_device(
          kernels_[i].first, N, rho, sigma, eps, vrho,
          vsigma, queue
        );
      else
        kernels_[i].second.eval_exc_vxc_inc_device(
          kernels_[i].first, N, rho, eps, vrho, queue
        );

    } else {

      double* eps_eval    = i ? eps_scr    : eps;
      double* vrho_eval   = i ? vrho_scr   : vrho;
      double* vsigma_eval = i ? vsigma_scr : vsigma;
      double* vlapl_eval  = i ? vlapl_scr  : vlapl;
      double* vtau_eval   = i ? vtau_scr   : vtau;

      if( kernels_[i].second.is_mgga() )
        kernels_[i].second.eval_exc_vxc_device(N, rho, sigma, lapl, tau, eps_eval,
          vrho_eval, vsigma_eval, vlapl_eval, vtau_eval, queue );
      else if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_vxc_device(N, rho, sigma, eps_eval, vrho_eval,
          vsigma_eval, queue );
      else
        kernels_[i].second.eval_exc_vxc_device(N, rho, eps_eval, vrho_eval, queue);

      if( i ) {

        add_scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, queue );
        add_scal_device( len_vrho_buffer, kernels_[i].first, vrho_eval, vrho, queue );

        if( kernels_[i].second.is_gga() )
          add_scal_device( len_vsigma_buffer, kernels_[i].first, vsigma_eval, vsigma, queue );

        if( kernels_[i].second.is_mgga() ) {
          add_scal_device( len_vlapl_buffer, kernels_[i].first, vlapl_eval, vlapl, queue );
          add_scal_device( len_vtau_buffer,  kernels_[i].first, vtau_eval,  vtau, queue  );
        }

      } else {

        scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, queue );
        scal_device( len_vrho_buffer, kernels_[i].first, vrho_eval, vrho, queue );

        if( kernels_[i].second.is_gga() or kernels_[i].second.is_mgga() )
          scal_device( len_vsigma_buffer, kernels_[i].first, vsigma_eval, vsigma, queue );

        if( kernels_[i].second.is_mgga() ) {
          scal_device( len_vlapl_buffer, kernels_[i].first, vlapl_eval, vlapl, queue );
          scal_device( len_vtau_buffer,  kernels_[i].first, vtau_eval,  vtau, queue  );
        }

      }
    }
  }

  if( eps_scr ) sycl::free( eps_scr, *queue );
  if( vrho_scr ) sycl::free( vrho_scr, *queue );
  if( vsigma_scr ) sycl::free( vsigma_scr, *queue );
  if( vlapl_scr ) sycl::free( vlapl_scr, *queue );
  if( vtau_scr ) sycl::free( vtau_scr, *queue );
}

}
