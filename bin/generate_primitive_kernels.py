import os,sys
import textwrap

def indent(text, amount, ch=' '):
  return textwrap.indent(text, amount * ch)

def _is_conv_to_float( x ):
  try:
    x = float(x)
    return True
  except ValueError:
    return False

class XCFunc:

  command_map = { 
    'EXC' : "g++ -P -DXC_DONT_COMPILE_{{V,F,K,L}}XC -E {}",
    'VXC' : "g++ -P -DXC_DONT_COMPILE_{{E,F,K,L}}XC -E {}", 
    'EXC_VXC' : "g++ -P -DXC_DONT_COMPILE_{{F,K,L}}XC -E {}" 
  }

  unpol_vars = {
    'rho[0]' : 'rho',
    'sigma[0]' : 'sigma',
    'zk[0]' : 'eps'
  }

  unpol_outputs = {
    'EXC' : ['eps'],
    'VXC' : ['vrho', 'vsigma'],
    'EXC_VXC' : ['eps','vrho', 'vsigma']
  }

  pol_vars = {
    'rho[0]'   : 'rho_a',
    'rho[1]'   : 'rho_b',
    'sigma[0]' : 'sigma_aa',
    'sigma[1]' : 'sigma_ab',
    'sigma[2]' : 'sigma_bb',
    'zk[0]' : 'eps'
  }

  pol_outputs = {
    'EXC' : ['eps'],
    'VXC' : ['vrho_a', 'vrho_b', 'vsigma_aa', 'vsigma_ab', 'vsigma_bb'],
    'EXC_VXC' : ['eps','vrho_a', 'vrho_b', 'vsigma_aa', 'vsigma_ab', 'vsigma_bb']
  }

  arith_ops = [ '+', '-', '*', '/' ]

  def __init__(self, fname, xc_type = 'EXC' ):
    self.fname = fname
    self.xc_type = xc_type

    # Generate Preprocessed output
    self.cmd = self.command_map[xc_type].format(fname)
    self.xc_out = os.popen(self.cmd).read()

    # Sanitize output
    self.xc_lines = self.sanitize_xc_out()

    # Separate out different functions
    self.separate_xc_out()

    # Handle special cases
    self.unpol_func_lines = self.special_cases( self.unpol_func_lines )
    self.ferr_func_lines  = self.special_cases( self.ferr_func_lines )
    self.pol_func_lines   = self.special_cases( self.pol_func_lines )


    self.const_params = []
    tmp_params, self.unpol_xc_body = self.finalize_lines( self.unpol_func_lines, self.unpol_vars, self.unpol_outputs[xc_type] )
    self.const_params = self.const_params + tmp_params

    tmp_params, self.ferr_xc_body = self.finalize_lines( self.ferr_func_lines, self.unpol_vars, self.unpol_outputs[xc_type] )
    self.const_params = self.const_params + tmp_params

    tmp_params, self.pol_xc_body = self.finalize_lines( self.pol_func_lines, self.pol_vars, self.pol_outputs[xc_type] )
    self.const_params = self.const_params + tmp_params

    self.const_params = list(set(self.const_params))

  def sanitize_xc_out( self ):
    
    self.xc_out = self.xc_out.replace('(','( ')
    self.xc_out = self.xc_out.replace(')',' )')
    #self.xc_out = self.xc_out.replace('params->','')
    self.xc_out = self.xc_out.replace('POW_1_3','cbrt')
    self.xc_out = self.xc_out.replace('POW_3_2','pow_3_2')
    self.xc_out = self.xc_out.replace('POW_2','square')
    self.xc_out = self.xc_out.replace('M_CBRT','constants::m_cbrt_')
    self.xc_out = self.xc_out.replace('  t','t')
    self.xc_out = self.xc_out.replace('    z','z')
    self.xc_out = self.xc_out.replace('    v','v')
    self.xc_out = self.xc_out.replace('0.1e1 / M_PI','constants::m_one_ov_pi')
    self.xc_out = self.xc_out.replace('M_PI * M_PI','constants::m_pi_sq')
    self.xc_out = self.xc_out.replace('cbrt( constants::m_one_ov_pi )','constants::m_cbrt_one_ov_pi')
    self.xc_out = self.xc_out.replace('cbrt( constants::m_pi_sq )','constants::m_cbrt_pi_sq')
    self.xc_out = self.xc_out.replace('my_piecewise3','piecewise_functor_3')
    self.xc_out = self.xc_out.replace('my_piecewise5','piecewise_functor_5')

    xc_lines = self.xc_out.splitlines()
    xc_lines = list(filter( lambda x: not x.startswith('  double'), xc_lines ))
    xc_lines = list(filter( lambda x: not x.startswith('  if'), xc_lines ))
    xc_lines = list(filter( lambda x: not x.startswith('  assert'), xc_lines ))
    xc_lines = list(filter( lambda x: '*params' not in x, xc_lines ))

    return xc_lines

  def separate_xc_out( self ):
    unpol_i = next( i for i,v in enumerate( self.xc_lines ) if 'func_unpol' in v )
    ferr_i  = next( i for i,v in enumerate( self.xc_lines ) if 'func_ferr'  in v )
    pol_i   = next( i for i,v in enumerate( self.xc_lines ) if 'func_pol'   in v )

    self.unpol_func_lines = self.xc_lines[(unpol_i+2):(ferr_i-2)]
    self.ferr_func_lines  = self.xc_lines[(ferr_i+2):(pol_i-2)]
    self.pol_func_lines   = self.xc_lines[(pol_i+2):-1]
  

  def special_cases( self, xc_lines ):
    one_ov_pi_var = [ line.split(' = ')[0].strip() for line in xc_lines if 'constants::m_one_ov_pi' in line ]
    if len(one_ov_pi_var) > 0:
      xc_lines = [ x.replace( 'cbrt( ' + one_ov_pi_var[0] + ' )', 'constants::m_cbrt_one_ov_pi' ) for x in xc_lines]
    pi_sq_var = [ line.split(' = ')[0].strip() for line in xc_lines if 'constants::m_pi_sq' in line ]
    if len(pi_sq_var) > 0:
      xc_lines = [ x.replace( 'cbrt( ' + pi_sq_var[0] + ' )', 'constants::m_cbrt_pi_sq' ) for x in xc_lines]
    return xc_lines


  def finalize_lines( self, _xc_lines, xc_vars, xc_output ):

    xc_lines = _xc_lines.copy()
    
    # this replaces the variables to sane text 
    for k,v in xc_vars.items():
      xc_lines = [ x.replace(k,v) for x in xc_lines ]

    # This changes all parameter arrays to fixed values
    xc_lines = [x.replace('[','_') for x in xc_lines ]
    xc_lines = [x.replace(']','')  for x in xc_lines ]

    res_lines = [ x for x in xc_lines if x.split(' = ')[0].strip() in xc_output ]
    #print(res_lines)

    for line in res_lines: xc_lines.remove( line )

    # Remove all output lines that are not in the specified xc_output
    xc_lines = [x for x in xc_lines if 't' in x.split(' = ')[0] or x.split(' = ')[0].strip() in xc_output ]

    # TODO: transverse the dependency tree to remove unneeded statements
    #req_vars = [ x.split(' = ')[1].strip().replace(';','') for x in res_lines ]
    #req_vars = " ".join(req_vars)

    #for op in self.arith_ops:
    #  req_vars = req_vars.replace(op,'')
    #req_vars = req_vars.split(' ')
    #req_vars = [x for x in req_vars if len(x) > 0 and 't' in x ]
    #print(req_vars)

    const_lines, xc_lines = self.get_const_lines( xc_lines, xc_vars )
    #print(const_lines)

    const_lines = [ 'constexpr double ' + x for x in const_lines ]
    xc_lines    = [ 'const double '    + x for x in xc_lines    ]


    # Get parameter lines
    const_params = []
    for line in const_lines:
      var, expr = line.split(' = ')
      expr.strip()
      for op in self.arith_ops:
        expr = expr.replace(op,'')
      expr = expr.split(' ')

      param_terms = [x.replace('params>','').replace(';','') for x in expr if 'param' in x ]
      const_params = const_params + param_terms
        
    const_params = list(set(const_params))
    #print(const_params)

    const_lines.append('\n')
    xc_lines.append('\n')

    xc_body_lines = const_lines + xc_lines + res_lines

    # Combine and perform final sanitization of parameters
    xc_body = "\n".join(xc_body_lines).replace('params->','')

    return const_params, xc_body


  def get_const_lines( self, _xc_lines, xc_vars ):
    xc_lines = _xc_lines.copy()
    const_lines = [x for x in xc_lines if 'constant' in x ]
    for line in const_lines: xc_lines.remove(line)

    const_vars = [ line.split(' = ')[0].strip() for line in const_lines ]
    for line in xc_lines:
      var, expr = line.split(' = ')
      var = var.strip()
      expr = expr.strip()

      expr = expr.replace(';','')
      for op in self.arith_ops:
        expr = expr.replace(op,'')
      expr = expr.split(' ')
      expr = [x for x in expr if len(x) > 0]

      special_functions = [x for x in expr if '(' in x or ')' in x]
      if len(special_functions) > 0: continue

      dep_vars = [x for x in expr if x in xc_vars.values()]
      if len(dep_vars) > 0: continue


      #numbers = [x for x in expr if 't' not in x]
      numbers = [ x for x in expr if _is_conv_to_float(x) ]
      for x in numbers: expr.remove(x)

      const_in_expr = [ x for x in expr if x in const_vars or 'params>' in x]
      for x in const_in_expr: expr.remove(x)

      if len(expr) > 0: continue

      const_lines.append(line)
      const_vars.append(var)
     
    for line in const_lines: 
      if line in xc_lines: xc_lines.remove(line)

    return const_lines, xc_lines



class GenMetaData:
  def __init__( self, local_name, libxc_file, ofname, xc_type, dens_tol, exx_coeff,
    params = {} ):
    self.local_name = local_name
    self.libxc_file = libxc_file
    self.ofname     = ofname
    self.xc_type    = xc_type
    self.dens_tol   = dens_tol
    self.exx_coeff  = exx_coeff
    self.params     = params

    self.is_hyb = abs(float(exx_coeff)) > 1e-10



libxc_prefix = '/ccs/home/dbwy/cscratch_chm136/xc_proxy_app_bin/external/libxc/src/maple2c/' 
kernel_prefix = './include/exchcxx/impl/builtin/kernels/'
gen_table = {

  'SlaterExchange' : GenMetaData( 'BuiltinSlaterExchange', 
    libxc_prefix + 'lda_exc/lda_x.c', 
    kernel_prefix + 'slater_exchange.hpp',
    'LDA', 1e-24, 0., { 'alpha' : '1.0' } 
    ),

  'VWN3' : GenMetaData( 'BuiltinVWN3', 
    libxc_prefix + 'lda_exc/lda_c_vwn_3.c', 
    kernel_prefix + 'vwn3.hpp',
    'LDA', 1e-24, 0. 
    ),

  'VWN_RPA' : GenMetaData( 'BuiltinVWN_RPA', 
    libxc_prefix + 'lda_exc/lda_c_vwn_rpa.c', 
    kernel_prefix + 'vwn_rpa.hpp',
    'LDA', 1e-24, 0. 
    ),

  'PW91' : GenMetaData( 'BuiltinPW91_LDA',
     libxc_prefix + 'lda_exc/lda_c_pw.c',
    kernel_prefix + 'pw91_lda.hpp',
     'LDA', 1e-24, 0., 
     { 'pp'     : '{1., 1., 1.}',
       'a'      : '{0.031091,  0.015545,   0.016887}',
       'alpha1' : '{0.21370,  0.20548,  0.11125}',
       'beta1'  : '{7.5957, 14.1189, 10.357}',
       'beta2'  : '{3.5876, 6.1977, 3.6231}',
       'beta3'  : '{1.6382, 3.3662, 0.88026}',
       'beta4'  : '{0.49294, 0.62517, 0.49671}',
       'fz20'   : '1.709921' }),

  'PW91_MOD' : GenMetaData( 'BuiltinPW91_LDA_MOD',
     libxc_prefix + 'lda_exc/lda_c_pw.c',
    kernel_prefix + 'pw91_lda_mod.hpp',
     'LDA', 1e-24, 0., 
     { 'pp'     : '{1., 1., 1.}',
       'a'      : '{0.0310907, 0.01554535, 0.0168869}',
       'alpha1' : '{0.21370,  0.20548,  0.11125}',
       'beta1'  : '{7.5957, 14.1189, 10.357}',
       'beta2'  : '{3.5876, 6.1977, 3.6231}',
       'beta3'  : '{1.6382, 3.3662, 0.88026}',
       'beta4'  : '{0.49294, 0.62517, 0.49671}',
       'fz20'   : '1.709920934161365617563962776245' }),

  'PW91_RPA' : GenMetaData( 'BuiltinPW91_LDA_RPA',
     libxc_prefix + 'lda_exc/lda_c_pw.c',
    kernel_prefix + 'pw91_lda_rpa.hpp',
     'LDA', 1e-24, 0., 
     { 
       'pp'     : '{0.75, 0.75, 1.0}',
       'a'      : '{0.031091,  0.015545,   0.016887}',
       'alpha1' : '{0.082477, 0.035374, 0.028829}',
       'beta1'  : '{ 5.1486, 6.4869, 10.357}',
       'beta2'  : '{1.6483, 1.3083, 3.6231}',
       'beta3'  : '{0.23647, 0.15180, 0.47990}',
       'beta4'  : '{0.20614, 0.082349, 0.12279}',
       'fz20'   : '1.709921' }),


  'PZ81' : GenMetaData( 'BuiltinPZ81',
     libxc_prefix + 'lda_exc/lda_c_pz.c',
    kernel_prefix + 'pz81.hpp',
     'LDA', 1e-24, 0., 
    {
      'gamma' : '{-0.1423, -0.0843}' ,  
      'beta1' : '{ 1.0529,  1.3981}' ,  
      'beta2' : '{ 0.3334,  0.2611}' ,  
      'a'     : '{ 0.0311,  0.01555}', 
      'b'     : '{-0.048,  -0.0269}' ,  
      'c'     : '{ 0.0020,  0.0007}' ,  
      'd'     : '{-0.0116, -0.0048}'   
    }),
       

  'PZ81_MOD' : GenMetaData( 'BuiltinPZ81_MOD',
     libxc_prefix + 'lda_exc/lda_c_pz.c',
    kernel_prefix + 'pz81_mod.hpp',
     'LDA', 1e-24, 0., 
    {
      'gamma' : '{-0.1423, -0.0843}'                          , 
      'beta1' : '{ 1.0529,  1.3981}'                          , 
      'beta2' : '{ 0.3334,  0.2611}'                          , 
      'a'     : '{ 0.0311,  0.01555}'                         ,
      'b'     : '{-0.048,  -0.0269}'                          , 
      'c'     : '{ 0.0020191519406228,  0.00069255121311694}' , 
      'd'     : '{-0.0116320663789130, -0.00480126353790614}'   
    }),


  'B88' : GenMetaData( 'BuiltinB88', 
    libxc_prefix + 'gga_exc/gga_x_b88.c', 
    kernel_prefix + 'b88.hpp',
    'GGA', 1e-25, 0., {'beta':'0.0042', 'gamma':'6.0'} 
    ),

  'LYP' : GenMetaData( 'BuiltinLYP', 
    libxc_prefix + 'gga_exc/gga_c_lyp.c', 
    kernel_prefix + 'lyp.hpp',
    'GGA', 1e-32, 0., {'A':'0.04918', 'B':'0.132', 'c':'0.2533', 'd':'0.349'} 
    ),

  'PBE_X' : GenMetaData( 'BuiltinPBE_X', 
    libxc_prefix + 'gga_exc/gga_x_pbe.c', 
    kernel_prefix + 'pbe_x.hpp',
    'GGA', 1e-32, 0., {'kappa':'0.8040', 'mu':' 0.2195149727645171'} 
    ),

  'PBE_C' : GenMetaData( 'BuiltinPBE_C', 
    libxc_prefix + 'gga_exc/gga_c_pbe.c', 
    kernel_prefix + 'pbe_c.hpp',
    'GGA', 1e-12, 0., {'beta':'0.06672455060314922', 'gamma':'0.031090690869654895034', 'BB':'1.'} 
    )

}








file_prefix = """#pragma once
#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

#include <exchcxx/impl/builtin/kernels/screening_interface.hpp>

"""

struct_prefix = """
namespace ExchCXX {{

template <>
struct kernel_traits< {0} > :
  public {1}_screening_interface< {0} > {{

  static constexpr bool is_lda  = {2};
  static constexpr bool is_gga  = {3};
  static constexpr bool is_mgga = {4};

  static constexpr double dens_tol  = {5};

  static constexpr bool is_hyb  = {6};
  static constexpr double exx_coeff = {7};
"""

lda_exc_args_unpolar     = "double rho, double& eps"
lda_exc_vxc_args_unpolar = "double rho, double& eps, double& vrho" 
gga_exc_args_unpolar     = "double rho, double sigma, double& eps"
gga_exc_vxc_args_unpolar = "double rho, double sigma, double& eps, double& vrho, double& vsigma" 

lda_exc_args_polar     = "double rho_a, double rho_b, double& eps"
lda_exc_vxc_args_polar = "double rho_a, double rho_b, double& eps, double& vrho_a, double& vrho_b" 
gga_exc_args_polar     = "double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps"
gga_exc_vxc_args_polar = "double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb" 


exc_args_unpolar = {
  'LDA' : lda_exc_args_unpolar,
  'GGA' : gga_exc_args_unpolar
}
exc_args_polar = {
  'LDA' : lda_exc_args_polar,
  'GGA' : gga_exc_args_polar
}
exc_vxc_args_unpolar = {
  'LDA' : lda_exc_vxc_args_unpolar,
  'GGA' : gga_exc_vxc_args_unpolar
}
exc_vxc_args_polar = {
  'LDA' : lda_exc_vxc_args_polar,
  'GGA' : gga_exc_vxc_args_polar
}



func_prefix = """
BUILTIN_KERNEL_EVAL_RETURN
  eval_{0}_{1}_impl( {2} ) {{
"""


builtin_xc_type_format = """
struct {0} : detail::BuiltinKernelImpl< {0} > {{

  {0}( Spin p ) :
    detail::BuiltinKernelImpl< {0} >(p) {{ }}
  
  virtual ~{0}() = default;

}};
"""


for name, meta_data in gen_table.items():
  xc_exc = XCFunc( meta_data.libxc_file, 'EXC' )
  xc_exc_vxc = XCFunc( meta_data.libxc_file, 'EXC_VXC' )

  xc_type = meta_data.xc_type

  is_lda  = xc_type == 'LDA'
  is_gga  = xc_type == 'GGA'
  is_mgga = xc_type == 'MGGA'

  xc_struct_prefix = struct_prefix.format(
    meta_data.local_name, xc_type.lower(), str(is_lda).lower(), 
    str(is_gga).lower(), str(is_mgga).lower(), meta_data.dens_tol,
    str(meta_data.is_hyb).lower(), meta_data.exx_coeff )

  xc_param_lines = []
  for pname, valstr in meta_data.params.items():
    if '{' in valstr:
      origvalstr = valstr
      valstr = valstr.replace('{','')
      valstr = valstr.replace('}','')
      valstr = valstr.split(',')
      valstr = [x.strip() for x in valstr]
      #nitem = len(valstr)
      #xc_param_lines.append('  static constexpr std::array<double,' + str(nitem) + '> ' + pname + ' = ' + origvalstr + ";")
      for i,v in enumerate(valstr):
        xc_param_lines.append('  static constexpr double ' + pname + '_' + str(i) + ' = ' + v + ';' )
    else:
      xc_param_lines.append('  static constexpr double ' + pname + ' = ' + valstr + ";")
  xc_param_lines = '\n'.join(xc_param_lines)

  exc_prefix_unpolar = func_prefix.format( 'exc', 'unpolar', 
    exc_args_unpolar[xc_type] )
  exc_vxc_prefix_unpolar = func_prefix.format( 'exc_vxc', 'unpolar', 
    exc_vxc_args_unpolar[xc_type] )

  exc_prefix_ferr = func_prefix.format( 'exc', 'ferr', 
    exc_args_unpolar[xc_type] )
  exc_vxc_prefix_ferr = func_prefix.format( 'exc_vxc', 'ferr', 
    exc_vxc_args_unpolar[xc_type] )

  exc_prefix_polar = func_prefix.format( 'exc', 'polar', 
    exc_args_polar[xc_type] )
  exc_vxc_prefix_polar = func_prefix.format( 'exc_vxc', 'polar', 
    exc_vxc_args_polar[xc_type] )


  exc_func_unpolar_body = "\n".join( [
    exc_prefix_unpolar, 
    indent(xc_exc.unpol_xc_body,2), 
    "\n}"
  ] )
  exc_func_ferr_body = "\n".join( [
    exc_prefix_ferr, 
    indent(xc_exc.ferr_xc_body,2), 
    "\n}"
  ] )
  exc_func_polar_body = "\n".join( [
    exc_prefix_polar, 
    indent(xc_exc.pol_xc_body,2), 
    "\n}"
  ] )

  exc_func_unpolar_body = indent(exc_func_unpolar_body, 2)
  exc_func_ferr_body = indent(exc_func_ferr_body, 2)
  exc_func_polar_body = indent(exc_func_polar_body, 2)

  exc_vxc_func_unpolar_body = "\n".join( [
    exc_vxc_prefix_unpolar, 
    indent(xc_exc_vxc.unpol_xc_body,2), 
    "\n}"
  ] )
  exc_vxc_func_ferr_body = "\n".join( [
    exc_vxc_prefix_ferr, 
    indent(xc_exc_vxc.ferr_xc_body,2), 
    "\n}"
  ] )
  exc_vxc_func_polar_body = "\n".join( [
    exc_vxc_prefix_polar, 
    indent(xc_exc_vxc.pol_xc_body,2), 
    "\n}"
  ] )

  exc_vxc_func_unpolar_body = indent(exc_vxc_func_unpolar_body, 2)
  exc_vxc_func_ferr_body = indent(exc_vxc_func_ferr_body, 2)
  exc_vxc_func_polar_body = indent(exc_vxc_func_polar_body, 2)


  xc_struct_str = "\n".join([
    xc_struct_prefix, 
    xc_param_lines, 
    exc_func_unpolar_body,
    exc_vxc_func_unpolar_body,
    exc_func_ferr_body,
    exc_vxc_func_ferr_body,
    exc_func_polar_body,
    exc_vxc_func_polar_body,
    "\n\n};"
  ])

  
  builtin_impl_str = builtin_xc_type_format.format(meta_data.local_name)

  xc_file_body = "\n".join([
    file_prefix,
    xc_struct_str,
    builtin_impl_str,
    "\n\n} // namespace ExchCXX"
  ])

  with open( meta_data.ofname, 'w' ) as f:
    f.write( xc_file_body )
  
