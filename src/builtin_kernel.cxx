#include <exchcxx/impl/builtin/kernel.hpp>
#include <exchcxx/impl/builtin/util.hpp>

namespace ExchCXX {


// LDA interface
LDA_EXC_GENERATOR( BuiltinKernel::eval_exc ) const {
  disabled_lda_interface();
};
LDA_EXC_VXC_GENERATOR( BuiltinKernel::eval_exc_vxc ) const {
  disabled_lda_interface();
};

// GGA interface
GGA_EXC_GENERATOR( BuiltinKernel::eval_exc ) const {
  disabled_gga_interface();
};
GGA_EXC_VXC_GENERATOR( BuiltinKernel::eval_exc_vxc ) const {
  disabled_gga_interface();
};

// MGGA interface
MGGA_EXC_GENERATOR( BuiltinKernel::eval_exc ) const {
  disabled_mgga_interface();
};
MGGA_EXC_VXC_GENERATOR( BuiltinKernel::eval_exc_vxc ) const {
  disabled_mgga_interface();
};

#ifdef EXCHCXX_ENABLE_DEVICE

// LDA interface
LDA_EXC_GENERATOR( BuiltinKernel::eval_exc_device ) const {
  disabled_lda_interface();
};
LDA_EXC_VXC_GENERATOR( BuiltinKernel::eval_exc_vxc_device ) const {
  disabled_lda_interface();
};
LDA_EXC_GENERATOR_DEVICE( BuiltinKernel::eval_exc_device_async ) const {
  disabled_lda_interface();
};
LDA_EXC_VXC_GENERATOR_DEVICE( BuiltinKernel::eval_exc_vxc_device_async ) const {
  disabled_lda_interface();
};

// GGA interface
GGA_EXC_GENERATOR( BuiltinKernel::eval_exc_device ) const {
  disabled_gga_interface();
};
GGA_EXC_VXC_GENERATOR( BuiltinKernel::eval_exc_vxc_device ) const {
  disabled_gga_interface();
};
GGA_EXC_GENERATOR_DEVICE( BuiltinKernel::eval_exc_device_async ) const {
  disabled_gga_interface();
};
GGA_EXC_VXC_GENERATOR_DEVICE( BuiltinKernel::eval_exc_vxc_device_async ) const {
  disabled_gga_interface();
};

// MGGA interface
MGGA_EXC_GENERATOR( BuiltinKernel::eval_exc_device ) const {
  disabled_mgga_interface();
};
MGGA_EXC_VXC_GENERATOR( BuiltinKernel::eval_exc_vxc_device ) const {
  disabled_mgga_interface();
};
MGGA_EXC_GENERATOR_DEVICE( BuiltinKernel::eval_exc_device_async ) const {
  disabled_mgga_interface();
};
MGGA_EXC_VXC_GENERATOR_DEVICE( BuiltinKernel::eval_exc_vxc_device_async ) const {
  disabled_mgga_interface();
};

#endif

}

