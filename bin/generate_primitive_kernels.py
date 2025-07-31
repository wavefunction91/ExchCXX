import os,sys,re
import textwrap
import subprocess

def indent(text, amount, ch=' '):
  return textwrap.indent(text, amount * ch)

def _is_conv_to_float( x ):
  try:
    x = float(x)
    return True
  except ValueError:
    return False

def add_double_to_piecewise_functor_3_args(xc_lines):
  pattern = re.compile(r'piecewise_functor_3\(([^,]+),\s*([^,]+),\s*([^)]+)\)')
  def replacer(match):
    arg1 = match.group(1).strip()
    arg2 = match.group(2).strip()
    arg3 = match.group(3).strip()
    
    if _is_conv_to_float(arg2):
      arg2 = str(float(arg2))
    if _is_conv_to_float(arg3):
      arg3 = str(float(arg3))

    return f'piecewise_functor_3( {arg1}, {arg2}, {arg3} )'
  return [pattern.sub(replacer, line) for line in xc_lines]

class XCFunc:

  command_map = { 
    'EXC' : "g++ -P -DXC_DONT_COMPILE_{{V,F,K,L}}XC -E {}",
    'VXC' : "g++ -P -DXC_DONT_COMPILE_{{E,F,K,L}}XC -E {}", 
    'FXC' : "g++ -P -DXC_DONT_COMPILE_{{E,V,K,L}}XC -E {}",
    'EXC_VXC' : "g++ -P -DXC_DONT_COMPILE_{{E,F,K,L}}XC -E {}",
    'VXC_FXC' : "g++ -P -DXC_DONT_COMPILE_{{E,V,K,L}}XC -E {}"
  }

  unpol_lda_vars = {
    'rho[0]' : 'rho',
  }
  unpol_gga_vars = {
    'rho[0]' : 'rho',
    'sigma[0]' : 'sigma',
  }

  unpol_mgga_vars = {
    'rho[0]' : 'rho',
    'sigma[0]' : 'sigma',
    'lapl[0]' : 'lapl',
    'tau[0]' : 'tau',
  }

  unpol_vars = {
    'LDA' : unpol_lda_vars, 
    'GGA' : unpol_gga_vars,
    'MGGA' : unpol_mgga_vars
  }








  unpol_lda_outputs = {
    'EXC' : ['eps'],
    'VXC' : ['vrho'],
    'FXC' : ['v2rho2'],
    'EXC_VXC' : ['eps','vrho'],
    'VXC_FXC' : ['vrho', 'v2rho2']
  }

  unpol_gga_outputs = {
    'EXC' : ['eps'],
    'VXC' : ['vrho', 'vsigma'],
    'FXC' : ['v2rho2', 'v2rhosigma', 'v2sigma2'],
    'EXC_VXC' : ['eps','vrho', 'vsigma'],
    'VXC_FXC' : ['vrho', 'vsigma', 'v2rho2', 'v2rhosigma', 'v2sigma2']
  }

  unpol_mgga_outputs = {
    'EXC' : ['eps'],
    'VXC' : ['vrho', 'vsigma', 'vlapl', 'vtau'],
    'FXC' : ['v2rho2', 'v2rhosigma', 'v2rholapl', 'v2rhotau', 'v2sigma2', 'v2sigmalapl', 'v2sigmatau', 'v2lapl2', 'v2lapltau', 'v2tau2'],
    'EXC_VXC' : ['eps', 'vrho', 'vsigma', 'vlapl', 'vtau'],
    'VXC_FXC' : ['vrho', 'vsigma', 'vlapl', 'vtau', 'v2rho2', 'v2rhosigma', 'v2rholapl', 'v2rhotau', 'v2sigma2', 'v2sigmalapl', 'v2sigmatau', 'v2lapl2', 'v2lapltau', 'v2tau2']
  }

  unpol_outputs = {
    'LDA' : unpol_lda_outputs, 
    'GGA' : unpol_gga_outputs,
    'MGGA' : unpol_mgga_outputs
  }
  unpol_substitutions = {
    'tzk0' : 'eps',
    'tvrho0' : 'vrho',
    'tvsigma0': 'vsigma',
    'tvlapl0' : 'vlapl',
    'tvtau0': 'vtau',
    'tv2rho20': 'v2rho2',
    'tv2rhosigma0': 'v2rhosigma',
    'tv2rholapl0': 'v2rholapl',
    'tv2rhotau0': 'v2rhotau',
    'tv2sigma20': 'v2sigma2',
    'tv2sigmalapl0': 'v2sigmalapl',
    'tv2sigmatau0': 'v2sigmatau',
    'tv2lapl20': 'v2lapl2',
    'tv2lapltau0': 'v2lapltau',
    'tv2tau20': 'v2tau2'
  }




  pol_lda_vars = {
    'rho[0]'   : 'rho_a',
    'rho[1]'   : 'rho_b',
  }

  pol_gga_vars = {
    'rho[0]'   : 'rho_a',
    'rho[1]'   : 'rho_b',
    'sigma[0]' : 'sigma_aa',
    'sigma[1]' : 'sigma_ab',
    'sigma[2]' : 'sigma_bb',
  }

  pol_mgga_vars = {
    'rho[0]'   : 'rho_a',
    'rho[1]'   : 'rho_b',
    'sigma[0]' : 'sigma_aa',
    'sigma[1]' : 'sigma_ab',
    'sigma[2]' : 'sigma_bb',
    'lapl[0]' : 'lapl_a',
    'lapl[1]' : 'lapl_b',
    'tau[0]'  : 'tau_a',
    'tau[1]'  : 'tau_b',
  }

  pol_vars = {
    'LDA' : pol_lda_vars, 
    'GGA' : pol_gga_vars,
    'MGGA': pol_mgga_vars
  }







  pol_lda_outputs = {
    'EXC' : ['eps'],
    'VXC' : ['vrho_a', 'vrho_b'],
    'FXC' : ['v2rho2_aa', 'v2rho2_ab', 'v2rho2_bb'],
    'EXC_VXC' : ['eps','vrho_a', 'vrho_b'],
    'VXC_FXC' : ['vrho_a', 'vrho_b', 'v2rho2_aa', 'v2rho2_ab', 'v2rho2_bb']
  }

  pol_gga_outputs = {
    'EXC' : ['eps'],
    'VXC' : ['vrho_a', 'vrho_b', 'vsigma_aa', 'vsigma_ab', 'vsigma_bb'],
    'FXC' : ['v2rho2_aa', 'v2rho2_ab', 'v2rho2_bb', 
             'v2rhosigma_a_aa', 'v2rhosigma_a_ab', 'v2rhosigma_a_bb',
             'v2rhosigma_b_aa', 'v2rhosigma_b_ab', 'v2rhosigma_b_bb', 
             'v2sigma2_aa_aa', 'v2sigma2_aa_ab', 'v2sigma2_aa_bb', 
             'v2sigma2_ab_ab', 'v2sigma2_ab_bb', 'v2sigma2_bb_bb'],
    'EXC_VXC' : ['eps','vrho_a', 'vrho_b', 'vsigma_aa', 'vsigma_ab', 'vsigma_bb'],
    'VXC_FXC' : ['vrho_a', 'vrho_b', 'vsigma_aa', 'vsigma_ab', 'vsigma_bb', 
                'v2rho2_aa', 'v2rho2_ab', 'v2rho2_bb', 
                'v2rhosigma_a_aa', 'v2rhosigma_a_ab', 'v2rhosigma_a_bb',
                'v2rhosigma_b_aa', 'v2rhosigma_b_ab', 'v2rhosigma_b_bb', 
                'v2sigma2_aa_aa', 'v2sigma2_aa_ab', 'v2sigma2_aa_bb', 
                'v2sigma2_ab_ab', 'v2sigma2_ab_bb', 'v2sigma2_bb_bb']
  }

  pol_mgga_outputs = {
    'EXC' : ['eps'],
    'VXC' : ['vrho_a', 'vrho_b', 'vsigma_aa', 'vsigma_ab', 'vsigma_bb', 'vlapl_a', 'vlapl_b', 'vtau_a', 'vtau_b'],
    'FXC' : ['v2rho2_aa', 'v2rho2_ab', 'v2rho2_bb', 
             'v2rhosigma_a_aa', 'v2rhosigma_a_ab', 'v2rhosigma_a_bb',
             'v2rhosigma_b_aa', 'v2rhosigma_b_ab', 'v2rhosigma_b_bb', 
             'v2rholapl_a_a', 'v2rholapl_a_b', 'v2rholapl_b_a', 'v2rholapl_b_b',
             'v2rhotau_a_a', 'v2rhotau_a_b', 'v2rhotau_b_a', 'v2rhotau_b_b', 
             'v2sigma2_aa_aa', 'v2sigma2_aa_ab', 'v2sigma2_aa_bb', 
             'v2sigma2_ab_ab', 'v2sigma2_ab_bb', 'v2sigma2_bb_bb', 
             'v2sigmalapl_aa_a', 'v2sigmalapl_aa_b', 'v2sigmalapl_ab_a', 'v2sigmalapl_ab_b', 'v2sigmalapl_bb_a', 'v2sigmalapl_bb_b',
             'v2sigmatau_aa_a', 'v2sigmatau_aa_b', 'v2sigmatau_ab_a', 'v2sigmatau_ab_b', 'v2sigmatau_bb_a', 'v2sigmatau_bb_b',
             'v2lapl2_aa', 'v2lapl2_ab', 'v2lapl2_bb', 
             'v2lapltau_a_a', 'v2lapltau_a_b', 'v2lapltau_b_a', 'v2lapltau_b_b',
             'v2tau2_aa', 'v2tau2_ab', 'v2tau2_bb'],
    'EXC_VXC' : ['eps','vrho_a', 'vrho_b', 'vsigma_aa', 'vsigma_ab', 'vsigma_bb', 'vlapl_a', 'vlapl_b', 'vtau_a', 'vtau_b'],
    'VXC_FXC' : ['vrho_a', 'vrho_b', 'vsigma_aa', 'vsigma_ab', 'vsigma_bb', 'vlapl_a', 'vlapl_b', 'vtau_a', 'vtau_b',
                'v2rho2_aa', 'v2rho2_ab', 'v2rho2_bb', 
                'v2rhosigma_a_aa', 'v2rhosigma_a_ab', 'v2rhosigma_a_bb',
                'v2rhosigma_b_aa', 'v2rhosigma_b_ab', 'v2rhosigma_b_bb', 
                'v2rholapl_a_a', 'v2rholapl_a_b', 'v2rholapl_b_a', 'v2rholapl_b_b',
                'v2rhotau_a_a', 'v2rhotau_a_b', 'v2rhotau_b_a', 'v2rhotau_b_b', 
                'v2sigma2_aa_aa', 'v2sigma2_aa_ab', 'v2sigma2_aa_bb', 
                'v2sigma2_ab_ab', 'v2sigma2_ab_bb', 'v2sigma2_bb_bb', 
                'v2sigmalapl_aa_a', 'v2sigmalapl_aa_b', 'v2sigmalapl_ab_a', 'v2sigmalapl_ab_b', 'v2sigmalapl_bb_a', 'v2sigmalapl_bb_b',
                'v2sigmatau_aa_a', 'v2sigmatau_aa_b', 'v2sigmatau_ab_a', 'v2sigmatau_ab_b', 'v2sigmatau_bb_a', 'v2sigmatau_bb_b',
                'v2lapl2_aa', 'v2lapl2_ab', 'v2lapl2_bb', 
                'v2lapltau_a_a', 'v2lapltau_a_b', 'v2lapltau_b_a', 'v2lapltau_b_b',
                'v2tau2_aa', 'v2tau2_ab', 'v2tau2_bb']
  }

  zero_out_spin_cross_terms = ['v2rho2_ab', 'v2rhosigma_a_ab', 'v2rhosigma_a_bb', 'v2rhosigma_b_aa', 'v2rhosigma_b_ab', \
             'v2rholapl_a_b', 'v2rholapl_b_a', 'v2rhotau_a_b', 'v2rhotau_b_a',  \
             'v2sigma2_aa_ab', 'v2sigma2_aa_bb', 'v2sigma2_ab_ab', 'v2sigma2_ab_bb', \
             'v2sigmalapl_aa_b', 'v2sigmalapl_ab_a', 'v2sigmalapl_ab_b', 'v2sigmalapl_bb_a', \
             'v2sigmatau_aa_b', 'v2sigmatau_ab_a', 'v2sigmatau_ab_b', 'v2sigmatau_bb_a', \
              'v2lapl2_ab', 'v2lapltau_a_b', 'v2lapltau_b_a', 'v2tau2_ab']

  pol_outputs = {
    'LDA' : pol_lda_outputs, 
    'GGA' : pol_gga_outputs,
    'MGGA': pol_mgga_outputs
  }

  pol_substitutions = {
    'tzk0' : 'eps',
    'tvrho0' : 'vrho_a',
    'tvrho1' : 'vrho_b',
    'tvsigma0' : 'vsigma_aa',
    'tvsigma1' : 'vsigma_ab',
    'tvsigma2' : 'vsigma_bb',
    'tvlapl0'  : 'vlapl_a',
    'tvlapl1'  : 'vlapl_b',
    'tvtau0'   : 'vtau_a',
    'tvtau1'   : 'vtau_b',
    'tv2rho20' : 'v2rho2_aa',
    'tv2rho21' : 'v2rho2_ab',
    'tv2rho22' : 'v2rho2_bb',
    'tv2rhosigma0' : 'v2rhosigma_a_aa',
    'tv2rhosigma1' : 'v2rhosigma_a_ab',
    'tv2rhosigma2' : 'v2rhosigma_a_bb',
    'tv2rhosigma3' : 'v2rhosigma_b_aa',
    'tv2rhosigma4' : 'v2rhosigma_b_ab',
    'tv2rhosigma5' : 'v2rhosigma_b_bb',
    'tv2rholapl0' : 'v2rholapl_a_a',
    'tv2rholapl1' : 'v2rholapl_a_b',
    'tv2rholapl2' : 'v2rholapl_b_a',
    'tv2rholapl3' : 'v2rholapl_b_b',
    'tv2rhotau0' : 'v2rhotau_a_a',
    'tv2rhotau1' : 'v2rhotau_a_b',
    'tv2rhotau2' : 'v2rhotau_b_a',
    'tv2rhotau3' : 'v2rhotau_b_b',
    'tv2sigma20' : 'v2sigma2_aa_aa',
    'tv2sigma21' : 'v2sigma2_aa_ab',
    'tv2sigma22' : 'v2sigma2_aa_bb',
    'tv2sigma23' : 'v2sigma2_ab_ab',
    'tv2sigma24' : 'v2sigma2_ab_bb',
    'tv2sigma25' : 'v2sigma2_bb_bb',
    'tv2sigmalapl0' : 'v2sigmalapl_aa_a',
    'tv2sigmalapl1' : 'v2sigmalapl_aa_b',
    'tv2sigmalapl2' : 'v2sigmalapl_ab_a',
    'tv2sigmalapl3' : 'v2sigmalapl_ab_b',
    'tv2sigmalapl4' : 'v2sigmalapl_bb_a',
    'tv2sigmalapl5' : 'v2sigmalapl_bb_b',
    'tv2sigmatau0' : 'v2sigmatau_aa_a',
    'tv2sigmatau1' : 'v2sigmatau_aa_b',
    'tv2sigmatau2' : 'v2sigmatau_ab_a',
    'tv2sigmatau3' : 'v2sigmatau_ab_b',
    'tv2sigmatau4' : 'v2sigmatau_bb_a',
    'tv2sigmatau5' : 'v2sigmatau_bb_b',
    'tv2lapl20' : 'v2lapl2_aa',
    'tv2lapl21' : 'v2lapl2_ab',
    'tv2lapl22' : 'v2lapl2_bb',
    'tv2lapltau0' : 'v2lapltau_a_a',
    'tv2lapltau1' : 'v2lapltau_a_b',
    'tv2lapltau2' : 'v2lapltau_b_a',
    'tv2lapltau3' : 'v2lapltau_b_b',
    'tv2tau20' : 'v2tau2_aa',
    'tv2tau21' : 'v2tau2_ab',
    'tv2tau22' : 'v2tau2_bb'
  }

  arith_ops = [ '+', '-', '*', '/' ]

  def __init__(self, fname, xc_approx, xc_type):
    self.fname = fname
    self.xc_approx = xc_approx
    self.xc_type = xc_type

    # Generate Preprocessed output
    self.cmd = self.command_map[xc_type].format(fname)
    self.xc_out = subprocess.check_output(self.cmd, shell=True, executable="/bin/bash").decode()

    # Sanitize output
    self.xc_lines = self.sanitize_xc_out()

    # Separate out different functions
    self.separate_xc_out()

    # Handle special cases
    self.unpol_func_lines = self.special_cases( self.unpol_func_lines )
    self.pol_func_lines   = self.special_cases( self.pol_func_lines )



    fname_end = self.fname.split('/')[-1]
    zero_out_spin_cross = ("_x.c" in fname_end or "_x_" in fname_end) and self.xc_type in ["FXC", "VXC_FXC"]

    self.const_params = []
    tmp_params, self.unpol_xc_body = self.finalize_lines( self.unpol_func_lines, self.unpol_vars[xc_approx], self.unpol_outputs[xc_approx][xc_type], self.unpol_substitutions)
    self.const_params = self.const_params + tmp_params

    tmp_params, self.pol_xc_body = self.finalize_lines( self.pol_func_lines, self.pol_vars[xc_approx], self.pol_outputs[xc_approx][xc_type], self.pol_substitutions, zero_out_spin_cross)
    self.const_params = self.const_params + tmp_params

    self.const_params = list(set(self.const_params))

  def sanitize_xc_out( self ):
    
    self.xc_out = self.xc_out.replace('(','( ')
    self.xc_out = self.xc_out.replace(')',' )')
    #self.xc_out = self.xc_out.replace('params->','')
    self.xc_out = self.xc_out.replace('POW_1_3','cbrt')
    self.xc_out = self.xc_out.replace('POW_3_2','pow_3_2')
    self.xc_out = self.xc_out.replace('POW_1_4','pow_1_4')
    self.xc_out = self.xc_out.replace('POW_2','square')
    self.xc_out = self.xc_out.replace('M_CBRT','constants::m_cbrt_')
    self.xc_out = self.xc_out.replace('  t','t')
    self.xc_out = self.xc_out.replace('    z','z')
    self.xc_out = self.xc_out.replace('    v','v')
    self.xc_out = self.xc_out.replace('0.1e1 / M_PI','constants::m_one_ov_pi')
    self.xc_out = self.xc_out.replace('M_PI * M_PI','constants::m_pi_sq')
    #self.xc_out = self.xc_out.replace('DBL_EPSILON', 'std::numeric_limits<double>::epsilon()')
    self.xc_out = self.xc_out.replace(' 0,',' 0.0,')
    self.xc_out = self.xc_out.replace(' 0 ',' 0.0 ')
    self.xc_out = self.xc_out.replace(' 1 ',' 1.0 ')
    self.xc_out = self.xc_out.replace(' 1,',' 1.0,')
    self.xc_out = self.xc_out.replace('cbrt( constants::m_one_ov_pi )','constants::m_cbrt_one_ov_pi')
    self.xc_out = self.xc_out.replace('cbrt( constants::m_pi_sq )','constants::m_cbrt_pi_sq')
    self.xc_out = self.xc_out.replace('constants::m_cbrt_PI','constants::m_cbrt_pi')
    self.xc_out = self.xc_out.replace('my_piecewise3','piecewise_functor_3')
    self.xc_out = self.xc_out.replace('my_piecewise5','piecewise_functor_5')
    self.xc_out = self.xc_out.replace('p->zeta_threshold','zeta_tol')
    self.xc_out = self.xc_out.replace('p->dens_threshold','dens_tol')
    self.xc_out = self.xc_out.replace('p->sigma_threshold', 'sigma_tol')
    self.xc_out = self.xc_out.replace('p->tau_threshold', 'tau_tol')
    self.xc_out = self.xc_out.replace('p->cam_omega', 'omega')

    xc_lines = self.xc_out.splitlines()
    xc_lines = list(filter( lambda x: not x.startswith('  double'), xc_lines ))
    xc_lines = list(filter( lambda x: not x.startswith('  if'), xc_lines ))
    xc_lines = list(filter( lambda x: not x.startswith('  assert'), xc_lines ))
    xc_lines = list(filter( lambda x: not x.startswith('    out->'), xc_lines ))
    xc_lines = list(filter( lambda x: '*params' not in x, xc_lines ))

    # # If piecewise_functor_3() is used, convert any integer to a double float
    xc_lines = add_double_to_piecewise_functor_3_args(xc_lines)

    return xc_lines

  def separate_xc_out( self ):
    # Map XC type to appropriate function name in output
    if self.xc_type == "EXC":
        string = "exc"
    elif self.xc_type == "VXC" or self.xc_type == "EXC_VXC":
        string = "vxc"
    elif self.xc_type == "FXC" or self.xc_type == "VXC_FXC":
        string = "fxc"
    else:
        raise ValueError(f"Unknown XC type: {self.xc_type}")
      
    try:
        unpol_i = next( i for i,v in enumerate( self.xc_lines ) if f'func_{string}_unpol' in v )
        pol_i   = next( i for i,v in enumerate( self.xc_lines ) if f'func_{string}_pol'   in v )
        
        self.unpol_func_lines = self.xc_lines[(unpol_i+2):(pol_i-2)]
        self.pol_func_lines   = self.xc_lines[(pol_i+2):-1]
    except StopIteration:
        # If patterns aren't found, initialize with empty lists and print a warning
        print(f"Warning: Could not find expected function patterns for {self.xc_type}")
        print(self.xc_lines)
        sys.exit()
        
  

  def special_cases( self, xc_lines ):
    one_ov_pi_var = [ line.split(' = ')[0].strip() for line in xc_lines if 'constants::m_one_ov_pi' in line ]
    if len(one_ov_pi_var) > 0:
      xc_lines = [ x.replace( 'cbrt( ' + one_ov_pi_var[0] + ' )', 'constants::m_cbrt_one_ov_pi' ) for x in xc_lines]
    pi_sq_var = [ line.split(' = ')[0].strip() for line in xc_lines if 'constants::m_pi_sq' in line ]
    if len(pi_sq_var) > 0:
      xc_lines = [ x.replace( 'cbrt( ' + pi_sq_var[0] + ' )', 'constants::m_cbrt_pi_sq' ) for x in xc_lines]

    # replace cbrt
    xc_lines = [ x.replace('cbrt(','safe_math::cbrt(') for x in xc_lines ]
    xc_lines = [ x.replace('sqrt(','safe_math::sqrt(') for x in xc_lines ]
    xc_lines = [ x.replace('log(','safe_math::log(') for x in xc_lines ]
    xc_lines = [ x.replace('exp(','safe_math::exp(') for x in xc_lines ]
    xc_lines = [ x.replace('pow(','safe_math::pow(') for x in xc_lines ]
    xc_lines = [ x.replace('atan(','safe_math::atan(') for x in xc_lines ]
    xc_lines = [ x.replace('erf(','safe_math::erf(') for x in xc_lines ]
    xc_lines = [ x.replace('xc_E1_scaled(','safe_math::xc_E1_scaled(') for x in xc_lines ]
    xc_lines = [ x.replace('xc_erfcx(','safe_math::xc_erfcx(') for x in xc_lines ]
    return xc_lines


  def finalize_lines( self, _xc_lines, xc_vars, xc_output, substitutions, zero_out_spin_cross = False  ):

    xc_lines = _xc_lines.copy()
    
    # this replaces the variables to sane text 
    for k,v in xc_vars.items():
      xc_lines = [ x.replace(k,v) for x in xc_lines ]

    # Create a new list rather than modifying while iterating
    if zero_out_spin_cross:
      zero_terms = [output for output in xc_output if output in self.zero_out_spin_cross_terms]
      xc_output = [output for output in xc_output if output not in self.zero_out_spin_cross_terms]
    else:
      zero_terms = []

    # print("zero terms ", zero_terms)
    # print("xc_output ", xc_output)

    for k,v in substitutions.items():
      if v in xc_output:
        xc_lines = [ x.replace(k,v) for x in xc_lines ]

    # This changes all parameter arrays to fixed values
    xc_lines = [x.replace('[','_') for x in xc_lines ]
    xc_lines = [x.replace(']','')  for x in xc_lines ]

    res_lines = [ x for x in xc_lines if x.split(' = ')[0].strip() in xc_output ]

    for line in res_lines: xc_lines.remove( line )

    # Remove all output lines that are not in the specified xc_output
    xc_lines = [x for x in xc_lines if 't' in x.split(' = ')[0] or x.split(' = ')[0].strip() in xc_output ]

    # TODO: transverse the dependency tree to remove unneeded statements

    # Get RHS of result lines
    req_vars = [ x.split(' = ')[1].strip().replace(';','').replace(',','') for x in res_lines ]
    req_vars = " ".join(req_vars)

    # Strip RHS into individual variables
    for op in self.arith_ops:
      req_vars = req_vars.replace(op,'')
    req_vars = req_vars.split(' ')
    # additional_req_vars = [x for x in req_vars if x in substitutions.values()]
    req_vars = [x for x in req_vars if len(x) > 0 and 't' in x and 'tau' not in x]
    req_vars = [x for x in req_vars if 'functor' not in x ]
    req_vars = list(req_vars)

    # Get a list of all variables and their dependencies
    all_vars = dict()
    for line in xc_lines:
      tvar, rhs = line.split(' = ')
      tvar = tvar.strip()
      rhs  = rhs.strip().replace(';','')
      rhs  = rhs.strip().replace(',','')
      for op in self.arith_ops:
        rhs = rhs.replace(op,'')
      rhs = rhs.split(' ')
      rhs = [ x for x in rhs if len(x) > 0 and 't' in x and 'tau' not in x ]
      rhs = [ x for x in rhs if 'numeric' not in x ]
      rhs = [ x for x in rhs if 'constant' not in x ]
      rhs = [ x for x in rhs if 'cbrt' not in x ]
      rhs = [ x for x in rhs if 'atan' not in x ]
      rhs = [ x for x in rhs if 'sqrt' not in x ]
      rhs = [ x for x in rhs if 'param' not in x ]
      rhs = [ x for x in rhs if 'functor' not in x ]
      rhs = [ x for x in rhs if 'log' not in x ]
      rhs = [ x for x in rhs if 'pow' not in x ]
      rhs = [ x for x in rhs if 'exp' not in x ]
      rhs = [ x for x in rhs if 'dens_tol' not in x ]
      rhs = [ x for x in rhs if 'zeta_tol' not in x ]
      rhs = [ x for x in rhs if 'tau_tol' not in x ]
      rhs = [ x for x in rhs if 'sigma_tol' not in x ]
      rhs = [ x for x in rhs if 'xc_E1_scaled' not in x ]
      rhs = [ x for x in rhs if 'xc_erfcx' not in x ]
      rhs = [ x for x in rhs if 'erf' not in x ]
      
      #if len(rhs) > 0:
      all_vars[tvar] = set(rhs)

    # Transverse the dependency tree until list of required variables is unchanging
    have_all_req_vars = False
    refine_iter = 0
    while not have_all_req_vars:
      _req_vars_old = req_vars.copy()
      new_req_vars = set(req_vars) 
      
      for v in req_vars:
        #print( v, all_vars[v] )
        new_req_vars = new_req_vars.union( all_vars[v] )
      #print('new vars', new_req_vars)
      req_vars = list(new_req_vars)
      have_all_req_vars = new_req_vars == set(_req_vars_old)
      refine_iter = refine_iter + 1
      if refine_iter > 20: raise RuntimeError('Too many Refinements')
    
    #print( 'all req vars', req_vars )
    #sys.exit()

    # Remove unused lines
    unused_tvars = set(all_vars.keys()).difference( set(req_vars) )
    #if len(unused_vars): print('unused vars',unused_vars)
    unused_xc_lines = [ x for x in xc_lines if x.split(' = ')[0].strip() in unused_tvars ]
    for line in unused_xc_lines: xc_lines.remove( line )
    


    # Determine if there are unused xc_vars in the RHS 
    unused_xc_vars = []
    for v in xc_vars.values():
      dep_lines = [x for x in list(set(xc_lines).union(set(res_lines))) if v in x.split(' = ')[1].strip() ]
      if len(dep_lines) == 0: unused_xc_vars.append(v)
    #if len(unused_xc_vars): print( 'unused vars',unused_xc_vars )
   
    unused_xc_var_lines = [ '(void)('+x+');' for x in unused_xc_vars ]
    #if len(unused_xc_var_lines): print( 'unused vars',unused_xc_var_lines )
    




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
    if zero_out_spin_cross:
      zero_output_lines = [ f'{x} = 0.0;' for x in zero_terms ]
      zero_output_lines.append('\n')
    else:
      zero_output_lines = []
    

    xc_body_lines = unused_xc_var_lines + const_lines + xc_lines + res_lines + zero_output_lines
    #xc_body_lines = const_lines + xc_lines + res_lines

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

      const_in_expr = [ x for x in expr if x in const_vars or 'params>' in x ]
      for x in const_in_expr: expr.remove(x)

      if len(expr) > 0: continue

      const_lines.append(line)
      const_vars.append(var)
     
    for line in const_lines: 
      if line in xc_lines: xc_lines.remove(line)

    return const_lines, xc_lines



class GenMetaData:
  def __init__( self, local_name, libxc_file, ofname, xc_type, dens_tol,
    params = {}, needs_laplacian = False, is_kedf = False, is_epc = False ):
    self.local_name = local_name
    self.libxc_file = libxc_file
    self.ofname     = ofname
    self.xc_type    = xc_type
    self.dens_tol   = dens_tol
    self.params     = params
    self.needs_laplacian = needs_laplacian
    self.is_kedf = is_kedf
    self.is_epc  = is_epc




# libxc_prefix = '/home/davidwillia/Development/ExchCXX/libxc-7.0.0/src/maple2c/'
libxc_prefix = '/home/jsliang/madft-exchcxx/build/_deps/libxc-src/src/maple2c/'
kernel_prefix = '/home/jsliang/madft-exchcxx/include/exchcxx/impl/builtin/kernels/'
gen_table = {

  'B88' : GenMetaData( 'BuiltinB88', 
    libxc_prefix + 'gga_exc/gga_x_b88.c', 
    kernel_prefix + 'b88.hpp',
    'GGA', 1e-15, {'beta':'0.0042', 'gamma':'6.0'} 
    ),

  'M062X_X' : GenMetaData( 'BuiltinM062X_X',
    libxc_prefix + 'mgga_exc/hyb_mgga_x_m05.c',
    kernel_prefix + 'm06_2x_x.hpp',
    'MGGA', 1e-15,
    {'a_0' : '0.46',
     'a_1' : '-0.2206052',
     'a_2' : '-9.431788e-02',  
     'a_3' : '2.164494e+00', 
     'a_4' : '-2.556466e+00', 
     'a_5' : '-1.422133e+01',
     'a_6' : '1.555044e+01',  
     'a_7' : '3.598078e+01', 
     'a_8' : '-2.722754e+01', 
     'a_9' : '-3.924093e+01',  
     'a_10': '1.522808e+01',  
     'a_11': '1.522227e+01',
     'csi_HF': '1.0',
     'cx':     '0.54'}
    ), 

  'M05_X': GenMetaData( 'BuiltinM05_X',
    libxc_prefix + 'mgga_exc/hyb_mgga_x_m05.c',
    kernel_prefix + 'm05_x.hpp',
    'MGGA', 1e-15,
    { 'a_0': '1.0', 'a_1': '0.08151', 'a_2': '-0.43956', 'a_3': '-3.22422', 'a_4': '2.01819', 'a_5': '8.79431', 'a_6': '-0.00295',
      'a_7': '9.82029', 'a_8': '-4.82351', 'a_9': '-48.17574', 'a_10': '3.64802', 'a_11': '34.02248',
      'csi_HF': '0.72', 'cx': '0.28'}
    ),

  'M05_2X_X': GenMetaData( 'BuiltinM05_2X_X',
    libxc_prefix + 'mgga_exc/hyb_mgga_x_m05.c',
    kernel_prefix + 'm05_2x_x.hpp',
    'MGGA', 1e-15,
    { 'a_0': '1.0', 'a_1': '-0.56833', 'a_2': '-1.30057', 'a_3': '5.50070', 'a_4': '9.06402', 'a_5': '-32.21075', 'a_6': '-23.73298',
      'a_7': '70.22996', 'a_8': '29.88614', 'a_9': '-60.25778', 'a_10': '-13.22205', 'a_11': '15.23694',
      'csi_HF': '0.44', 'cx': '0.56'}
    ),


  'EPC17_1' : GenMetaData( 'BuiltinEPC17_1', 
    libxc_prefix + 'lda_exc/lda_c_epc17.c', 
    kernel_prefix + 'epc17_1.hpp',
    'LDA', 1e-24,
    { 'a': '2.35',
      'b': '2.40',
      'c': '3.20' },
    False, False, True 
    ),

  'EPC17_2' : GenMetaData( 'BuiltinEPC17_2', 
    libxc_prefix + 'lda_exc/lda_c_epc17.c', 
    kernel_prefix + 'epc17_2.hpp',
    'LDA', 1e-24, 
    { 'a': '2.35',
      'b': '2.40',
      'c': '6.60' },
    False, False, True 
    ),

  'EPC18_1' : GenMetaData( 'BuiltinEPC18_1', 
    libxc_prefix + 'lda_exc/lda_c_epc18.c', 
    kernel_prefix + 'epc18_1.hpp',
    'LDA', 1e-24, 
    { 'a': '1.80',
      'b': '0.10',
      'c': '0.03' },
    False, False, True 
    ),

  'EPC18_2' : GenMetaData( 'BuiltinEPC18_2', 
    libxc_prefix + 'lda_exc/lda_c_epc18.c', 
    kernel_prefix + 'epc18_2.hpp',
    'LDA', 1e-24,
    { 'a': '3.90',
      'b': '0.50',
      'c': '0.06' },
    False, False, True 
    ),

  'SlaterExchange' : GenMetaData( 'BuiltinSlaterExchange', 
    libxc_prefix + 'lda_exc/lda_x.c', 
    kernel_prefix + 'slater_exchange.hpp',
    'LDA', 1e-24, { 'alpha' : '1.0' } 
    ),

  'VWN3' : GenMetaData( 'BuiltinVWN3', 
    libxc_prefix + 'lda_exc/lda_c_vwn_3.c', 
    kernel_prefix + 'vwn3.hpp',
    'LDA', 1e-24, 
    ),

  'VWN_RPA' : GenMetaData( 'BuiltinVWN_RPA', 
    libxc_prefix + 'lda_exc/lda_c_vwn_rpa.c', 
    kernel_prefix + 'vwn_rpa.hpp',
    'LDA', 1e-24,
    ),

  'PW91' : GenMetaData( 'BuiltinPW91_LDA',
     libxc_prefix + 'lda_exc/lda_c_pw.c',
    kernel_prefix + 'pw91_lda.hpp',
     'LDA', 1e-24,
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
     'LDA', 1e-24,
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
     'LDA', 1e-24,
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
     'LDA', 1e-24,
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
     'LDA', 1e-24,
    {
      'gamma' : '{-0.1423, -0.0843}'                          , 
      'beta1' : '{ 1.0529,  1.3981}'                          , 
      'beta2' : '{ 0.3334,  0.2611}'                          , 
      'a'     : '{ 0.0311,  0.01555}'                         ,
      'b'     : '{-0.048,  -0.0269}'                          , 
      'c'     : '{ 0.0020191519406228,  0.00069255121311694}' , 
      'd'     : '{-0.0116320663789130, -0.00480126353790614}'   
    }),



  'LYP' : GenMetaData( 'BuiltinLYP', 
    libxc_prefix + 'gga_exc/gga_c_lyp.c', 
    kernel_prefix + 'lyp.hpp',
    'GGA', 1e-32, {'a':'0.04918', 'b':'0.132', 'c':'0.2533', 'd':'0.349'} 
    ),

  'PBE_X' : GenMetaData( 'BuiltinPBE_X', 
    libxc_prefix + 'gga_exc/gga_x_pbe.c', 
    kernel_prefix + 'pbe_x.hpp',
    'GGA', 1e-32, {'kappa':'0.8040', 'mu':' 0.2195149727645171'} 
    ),

  'revPBE_X' : GenMetaData( 'BuiltinRevPBE_X', 
    libxc_prefix + 'gga_exc/gga_x_pbe.c', 
    kernel_prefix + 'rev_pbe_x.hpp',
    'GGA', 1e-32, {'kappa':'1.245', 'mu':' 0.2195149727645171'} 
    ),

  'PBE_C' : GenMetaData( 'BuiltinPBE_C', 
    libxc_prefix + 'gga_exc/gga_c_pbe.c', 
    kernel_prefix + 'pbe_c.hpp',
    'GGA', 1e-12, {'beta':'0.06672455060314922', 'gamma':'0.031090690869654895034', 'BB':'1.'} 
    ),

  'SCAN_X' : GenMetaData( 'BuiltinSCAN_X',
    libxc_prefix + 'mgga_exc/mgga_x_scan.c',
    kernel_prefix + 'scan_x.hpp',
    'MGGA', 1e-15, 
    {'c1': '0.667',
     'c2': '0.8',
     'd' : '1.24',
     'k1': '0.065'}
    ),

  'SCAN_C' : GenMetaData( 'BuiltinSCAN_C',
    libxc_prefix + 'mgga_exc/mgga_c_scan.c',
    kernel_prefix + 'scan_c.hpp',
    'MGGA', 1e-15,
    ),

  'R2SCAN_X' : GenMetaData( 'BuiltinR2SCAN_X',
    libxc_prefix + 'mgga_exc/mgga_x_r2scan.c',
    kernel_prefix + 'r2scan_x.hpp',
    'MGGA', 1e-11,
    {'c1': '0.667',
     'c2': '0.8',
     'd': '1.24',
     'k1': '0.065',
     'eta': '0.001',
     'dp2': '0.361'}
    ),

  'R2SCAN_C' : GenMetaData( 'BuiltinR2SCAN_C',
    libxc_prefix + 'mgga_exc/mgga_c_r2scan.c',
    kernel_prefix + 'r2scan_c.hpp',
    'MGGA', 1e-15,
    {'eta': '0.001'}
    ),


  'M062X_C' : GenMetaData( 'BuiltinM062X_C',
    libxc_prefix + 'mgga_exc/mgga_c_m06l.c',
    kernel_prefix + 'm06_2x_c.hpp',
    'MGGA', 1e-12,
    {
      'gamma_ss':  '0.06', 
      'gamma_ab':  '0.0031', 
      'alpha_ss':  '0.00515088', 
      'alpha_ab':  '0.00304966',
      'css_0':     '3.097855e-01', 
      'css_1':    '-5.528642e+00',  
      'css_2':     '1.347420e+01', 
      'css_3':    '-3.213623e+01',  
      'css_4':     '2.846742e+01',
      'cab_0':     '8.833596e-01',  
      'cab_1':     '3.357972e+01', 
      'cab_2':    '-7.043548e+01',  
      'cab_3':     '4.978271e+01', 
      'cab_4':    '-1.852891e+01',
      'dss_0':     '6.902145e-01',  
      'dss_1':     '9.847204e-02',  
      'dss_2':     '2.214797e-01', 
      'dss_3':    '-1.968264e-03', 
      'dss_4':    '-6.775479e-03',  
      'dss_5':     '0.000000e+00',
      'dab_0':     '1.166404e-01', 
      'dab_1':    '-9.120847e-02', 
      'dab_2':    '-6.726189e-02',  
      'dab_3':     '6.720580e-05',  
      'dab_4':     '8.448011e-04',  
      'dab_5':     '0.000000e+00',
      'Fermi_D_cnst': '1e-10'
    }
    ), 

  'PKZB_X' : GenMetaData( 'BuiltinPKZB_X',
    libxc_prefix + 'mgga_exc/mgga_x_pkzb.c',
    kernel_prefix + 'pkzb_x.hpp',
    'MGGA', 1e-15, {}),

  'PKZB_C' : GenMetaData( 'BuiltinPKZB_C',
    libxc_prefix + 'mgga_exc/mgga_c_pkzb.c',
    kernel_prefix + 'pkzb_c.hpp',
    'MGGA', 1e-13, {}),

  'FT98' : GenMetaData( 'BuiltinFT98_X',
    libxc_prefix + 'mgga_exc/mgga_x_ft98.c',
    kernel_prefix + 'ft98_x.hpp',
    'MGGA', 1e-12,
    { 'a': '0.00528014',
      'b': '0.00003904539', 
      'a1': '2.816049', 
      'a2': '0.879058', 
      'b1': '0.398773', 
      'b2': '66.364138'},
    True
    ),

  'PC07' : GenMetaData( 'BuiltinPC07_K',
    libxc_prefix + 'mgga_exc/mgga_k_pc07.c',
    kernel_prefix + 'pc07_k.hpp',
    'MGGA', 1e-15,
    { 'a': '0.5389',
      'b': '3.0' },
    True, True
    ),

  'PC07OPT' : GenMetaData( 'BuiltinPC07OPT_K',
    libxc_prefix + 'mgga_exc/mgga_k_pc07.c',
    kernel_prefix + 'pc07opt_k.hpp',
    'MGGA', 1e-15,
    { 'a': '1.784720',
      'b': '0.258304' },
    True, True
    ),

   'TPSS_X' : GenMetaData( 'BuiltinTPSS_X',
     libxc_prefix + 'mgga_exc/mgga_x_tpss.c',
     kernel_prefix + 'tpss_x.hpp',
     'MGGA', 1e-15,
     { 'b': '0.40', 
       'c': '1.59096', 
       'e': '1.537', 
       'kappa' : '0.804', 
       'mu' : '0.21951',
       'BLOC_a' : '2.0',
       'BLOC_b' : '0.0' }
    ),

    'revTPSS_X' : GenMetaData( 'BuiltinRevTPSS_X',
      libxc_prefix + 'mgga_exc/mgga_x_tpss.c',
      kernel_prefix + 'revtpss_x.hpp',
      'MGGA', 1e-15,
      { 'b': '0.40', 
        'c': '2.35203946',
        'e': '2.16769874',
        'kappa' : '0.804',
        'mu' : '0.14',
        'BLOC_a' : '3.0',
        'BLOC_b' : '0.0' }
    ),

  'B97_D' : GenMetaData( 'BuiltinB97_D',
    libxc_prefix + 'gga_exc/gga_xc_b97.c',
    kernel_prefix + 'b97_d.hpp',
    'GGA', 1e-14, 
    { 'c_x_0' : '1.08662', 'c_x_1' : '-0.52127', 'c_x_2' : '3.25429', 'c_x_3' : '0.0', 'c_x_4' : '0.0',
      'c_ss_0' : '0.2234', 'c_ss_1' : '-1.56208', 'c_ss_2' : '1.94293', 'c_ss_3' : '0.0', 'c_ss_4' : '0.0',
      'c_ab_0' : '0.69041', 'c_ab_1' : '6.3027', 'c_ab_2' : '-14.9712', 'c_ab_3' : '0.0', 'c_ab_4' : '0.0',
      'c_xx'   : '0.0' }
    ),

  'ITYH_X' : GenMetaData( 'BuiltinITYH_X',
    libxc_prefix + 'gga_exc/gga_x_ityh.c',
    kernel_prefix + 'ityh_x.hpp',
    'GGA', 1e-14, 
    { 'omega' : '0.2' }
    ),

  'ITYH_X_033' : GenMetaData( 'BuiltinITYH_X_033',
    libxc_prefix + 'gga_exc/gga_x_ityh.c',
    kernel_prefix + 'ityh_x_033.hpp',
    'GGA', 1e-14, 
    { 'omega' : '0.33' }
    ),
    
  'ITYH_X_015' : GenMetaData( 'BuiltinITYH_X_015',
    libxc_prefix + 'gga_exc/gga_x_ityh.c',
    kernel_prefix + 'ityh_x_015.hpp',
    'GGA', 1e-14, 
    { 'omega' : '0.15' }
    ),

  'VWN' : GenMetaData( 'BuiltinVWN',
    libxc_prefix + 'lda_exc/lda_c_vwn.c',
    kernel_prefix + 'vwn.hpp',
    'LDA', 1e-15,
    ),

  'M06_L_X' : GenMetaData( 'BuiltinM06_L_X',
    libxc_prefix + 'mgga_exc/mgga_x_m06l.c',
    kernel_prefix + 'm06_l_x.hpp',
    'MGGA', 1e-15,
    { 'a_0' : '0.3987756', 'a_1' : '0.2548219', 'a_2' : '0.3923994', 'a_3' : '-2.103655', 'a_4' : '-6.302147',
      'a_5' : '10.97615', 'a_6' : '30.97273', 'a_7' : '-23.18489', 'a_8' : '-56.73480', 'a_9' : '21.60364',
      'a_10': '34.21814', 'a_11': '-9.049762', 'd_0' : '0.6012244', 'd_1' : '0.004748822', 'd_2' : '-0.008635108',
      'd_3' : '-0.000009308062', 'd_4' : '0.00004482811', 'd_5' : '0.0'}
    ),
    
  'M06_X' : GenMetaData( 'BuiltinM06_X',
    libxc_prefix + 'mgga_exc/mgga_x_m06l.c',
    kernel_prefix + 'm06_x.hpp',
    'MGGA', 1e-15,
    { 'a_0' : '0.5877943', 'a_1' : '-0.1371776', 'a_2' : '0.2682367', 'a_3' : '-2.515898', 'a_4' : '-2.978892',
      'a_5' : '8.710679', 'a_6' : '16.88195', 'a_7' : '-4.489724', 'a_8' : '-32.99983', 'a_9' : '-14.49050',
      'a_10': '20.43747', 'a_11': '12.56504', 'd_0' : '0.1422057', 'd_1' : '0.0007370319', 'd_2' : '-0.01601373',
      'd_3' : '0.0', 'd_4' : '0.0', 'd_5' : '0.0'}
    ),

    'revM06_L_X' : GenMetaData( 'BuiltinRevM06_L_X',
    libxc_prefix + 'mgga_exc/mgga_x_m06l.c',
    kernel_prefix + 'rev_m06_l_x.hpp',
    'MGGA', 1e-15,
    { 'a_0' : '1.423227252', 'a_1' : '0.471820438', 'a_2' : '-0.167555701', 'a_3' : '-0.250154262', 'a_4' : '0.062487588',
      'a_5' : '0.733501240', 'a_6' : '-2.359736776', 'a_7' : '-1.436594372', 'a_8' : '0.444643793', 'a_9' : '1.529925054',
      'a_10': '2.053941717', 'a_11' : '-0.036536031',
      'd_0' : '-0.423227252', 'd_1' : '0.0', 'd_2' : '0.003724234', 'd_3' : '0.0', 'd_4' : '0.0', 'd_5' : '0.0'}
    ),

    'M06_SX_X' : GenMetaData( 'BuiltinM06_SX_X',
    libxc_prefix + 'mgga_exc/mgga_x_m06l.c',
    kernel_prefix + 'm06_sx_x.hpp',
    'MGGA', 1e-15, 
    { 'a_0' : '0.996501680264007', 'a_1' : '0.0301264933631367', 'a_2' : '-0.103366758333673', 'a_3' : '-0.155653062500239', 'a_4' : '0.00795768051149902',
      'a_5' : '0.0871986277454856', 'a_6' : '-0.816152625764469', 'a_7' : '0.672773006612420', 'a_8' : '0.521127186174968', 'a_9' : '0.399466945122217',
      'a_10': '0.519400018999204', 'a_11': '-0.965261552636835', 'd_0' : '-0.347792307472902', 'd_1' : '0.0', 'd_2' : '-0.00270366787478266', 'd_3' : '0.0', 
      'd_4' : '0.0', 'd_5' : '0.0'}
    ),

  'M06_HF_X' : GenMetaData( 'BuiltinM06_HF_X',
    libxc_prefix + 'mgga_exc/mgga_x_m06l.c',
    kernel_prefix + 'm06_hf_x.hpp',
    'MGGA', 1e-15,
    { 'a_0' : '0.1179732', 'a_1' : '-1.066708', 'a_2' : '-0.1462405', 'a_3' : '7.481848', 'a_4' : '3.776679',
      'a_5' : '-44.36118', 'a_6' : '-18.30962', 'a_7' : '100.3903', 'a_8' : '38.64360', 'a_9' : '-98.06018',
      'a_10': '-25.57716', 'a_11': '35.90404', 'd_0' : '-0.1179732', 'd_1' : '-0.0025', 'd_2' : '-0.01180065',
      'd_3' : '0.0', 'd_4' : '0.0', 'd_5' : '0.0'}
    ),
                     
  'M06_L_C' : GenMetaData( 'BuiltinM06_L_C',
    libxc_prefix + 'mgga_exc/mgga_c_m06l.c',
    kernel_prefix + 'm06_l_c.hpp',
    'MGGA', 1e-12,
    {
      'gamma_ss':  '0.06', 'gamma_ab':  '0.0031', 'alpha_ss':  '0.00515088', 'alpha_ab':  '0.00304966',
      'css_0': '5.349466e-01', 'css_1': '5.396620e-01', 'css_2': '-3.161217e+01', 'css_3': '5.149592e+01', 'css_4': '-2.919613e+01',
      'cab_0': '6.042374e-01', 'cab_1': '1.776783e+02', 'cab_2': '-2.513252e+02', 'cab_3': '7.635173e+01', 'cab_4': '-1.255699e+01',
      'dss_0': '4.650534e-01', 'dss_1': '1.617589e-01', 'dss_2': '1.833657e-01', 'dss_3': '4.692100e-04', 'dss_4': '-4.990573e-03', 'dss_5': '0.000000e+00',
      'dab_0': '3.957626e-01', 'dab_1': '-5.614546e-01', 'dab_2': '1.403963e-02', 'dab_3': '9.831442e-04', 'dab_4': '-3.577176e-03', 'dab_5': '0.000000e+00',
      'Fermi_D_cnst': '1e-10'}
    ),
  'M06_HF_C' : GenMetaData( 'BuiltinM06_HF_C',
    libxc_prefix + 'mgga_exc/mgga_c_m06l.c',
    kernel_prefix + 'm06_hf_c.hpp',
    'MGGA', 1e-12,
    {
      'gamma_ss':  '0.06', 'gamma_ab':  '0.0031', 'alpha_ss':  '0.00515088', 'alpha_ab':  '0.00304966',
      'css_0': '1.023254e-01', 'css_1': '-2.453783e+00', 'css_2': '2.913180e+01', 'css_3': '-3.494358e+01', 'css_4': '2.315955e+01',
      'cab_0': '1.674634e+00', 'cab_1': '5.732017e+01', 'cab_2': '5.955416e+01', 'cab_3': '-2.311007e+02', 'cab_4': '1.255199e+02',
      'dss_0': '8.976746e-01', 'dss_1': '-2.345830e-01', 'dss_2': '2.368173e-01', 'dss_3': '-9.913890e-04', 'dss_4': '-1.146165e-02', 'dss_5': '0.000000e+00',
      'dab_0': '-6.746338e-01', 'dab_1': '-1.534002e-01', 'dab_2': '-9.021521e-02', 'dab_3': '-1.292037e-03', 'dab_4': '-2.352983e-04', 'dab_5': '0.000000e+00',
      'Fermi_D_cnst': '1e-10'}
    ),
  'M06_C' : GenMetaData( 'BuiltinM06_C',
    libxc_prefix + 'mgga_exc/mgga_c_m06l.c',
    kernel_prefix + 'm06_c.hpp',
    'MGGA', 1e-12,
    {
      'gamma_ss':  '0.06', 'gamma_ab':  '0.0031', 'alpha_ss':  '0.00515088', 'alpha_ab':  '0.00304966',
      'css_0': '5.094055e-01', 'css_1': '-1.491085e+00', 'css_2': '1.723922e+01', 'css_3': '-3.859018e+01', 'css_4': '2.845044e+01',
      'cab_0': '3.741539e+00', 'cab_1': '2.187098e+02', 'cab_2': '-4.531252e+02', 'cab_3': '2.936479e+02', 'cab_4': '-6.287470e+01',
      'dss_0': '4.905945e-01', 'dss_1': '-1.437348e-01', 'dss_2': '2.357824e-01', 'dss_3': '1.871015e-03', 'dss_4': '-3.788963e-03', 'dss_5': '0.000000e+00',
      'dab_0': '-2.741539e+00', 'dab_1': '-6.720113e-01', 'dab_2': '-7.932688e-02', 'dab_3': '1.918681e-03', 'dab_4': '-2.032902e-03', 'dab_5': '0.000000e+00',
      'Fermi_D_cnst': '1e-10'}
    ),
  'M06_SX_C' : GenMetaData( 'BuiltinM06_SX_C',
    libxc_prefix + 'mgga_exc/mgga_c_m06l.c',
    kernel_prefix + 'm06_sx_c.hpp',
    'MGGA', 1e-14,
    {
      'gamma_ss':  '0.06', 'gamma_ab':  '0.0031', 'alpha_ss':  '0.00515088', 'alpha_ab':  '0.00304966',
      'css_0': '1.17575011057022E+00', 'css_1': '6.58083496678423E-01', 'css_2': '-2.78913774852905E+00', 'css_3': '-1.18597601856255E+00', 'css_4': '1.16439928209688E+00',
      'cab_0': '1.63738167314691E-01', 'cab_1': '-4.36481171027951E-01', 'cab_2': '-1.90232628449712E+00', 'cab_3': '-1.42432902881841E+00', 'cab_4': '-9.05909137360893E-01',
      'dss_0': '8.17322574473352E-02', 'dss_1': '-2.88531085759385E-02', 'dss_2': '9.05917734868130E-02', 'dss_3': '0.0', 'dss_4': '0.0', 'dss_5': '-4.86297499082106E-04',
      'dab_0': '7.40594619832397E-01', 'dab_1': '1.23306511345974E-02', 'dab_2': '-1.88253421850249E-02', 'dab_3': '0.0', 'dab_4': '0.0', 'dab_5': '4.87276242162303E-04',
      'Fermi_D_cnst': '1e-10'}
    ),
  'revM06_L_C' : GenMetaData( 'BuiltinRevM06_L_C',
    libxc_prefix + 'mgga_exc/mgga_c_m06l.c',
    kernel_prefix + 'rev_m06_l_c.hpp',
    'MGGA', 1e-12,
    {
      'gamma_ss':  '0.06', 'gamma_ab':  '0.0031', 'alpha_ss':  '0.00515088', 'alpha_ab':  '0.00304966',
      'css_0': '1.227659748', 'css_1': '0.855201283', 'css_2': '-3.113346677', 'css_3': '-2.239678026', 'css_4': '0.354638962',
      'cab_0': '0.344360696', 'cab_1': '-0.557080242', 'cab_2': '-2.009821162', 'cab_3': '-1.857641887', 'cab_4': '-1.076639864',
      'dss_0': '-0.538821292', 'dss_1': '-0.028296030', 'dss_2': '0.023889696', 'dss_3': '0.0', 'dss_4': '0.0', 'dss_5': '-0.002437902',
      'dab_0': '0.400714600', 'dab_1': '0.015796569', 'dab_2': '-0.032680984', 'dab_3': '0.0', 'dab_4': '0.0', 'dab_5': '0.001260132',
      'Fermi_D_cnst': '1e-10'}
    ),

  'M05_C' : GenMetaData( 'BuiltinM05_C',
    libxc_prefix + 'mgga_exc/mgga_c_m05.c',
    kernel_prefix + 'm05_c.hpp',
    'MGGA', 1e-15,
    { 'gamma_ss':  '0.06', 'gamma_ab':  '0.0031',
      'css_0': '1.00000e0', 'css_1': '3.77344e0', 'css_2': '-26.04463e0', 'css_3': '30.69913e0', 'css_4': '-9.22695e0',
      'cab_0': '1.00000e0', 'cab_1': '3.78569e0', 'cab_2': '-14.15261e0', 'cab_3': '-7.46589e0', 'cab_4': '17.94491e0',
      'Fermi_D_cnst': '1e-10'}
    ),
  'M05_2X_C' : GenMetaData( 'BuiltinM05_2X_C',
    libxc_prefix + 'mgga_exc/mgga_c_m05.c',
    kernel_prefix + 'm05_2x_c.hpp',
    'MGGA', 1e-15,
    { 'gamma_ss':  '0.06', 'gamma_ab':  '0.0031',
      'css_0': '1.00000e0', 'css_1': '-3.05430e0', 'css_2': '7.61854e0', 'css_3': '1.47665e0', 'css_4': '-11.92365e0',
      'cab_0': '1.00000e0', 'cab_1': '1.09297e0', 'cab_2': '-3.79171e0', 'cab_3': '2.82810e0', 'cab_4': '-10.58909e0',
      'Fermi_D_cnst': '1e-10'}
    ),
                  
  'M08_HX_C' : GenMetaData( 'BuiltinM08_HX_C',
    libxc_prefix + 'mgga_exc/mgga_c_m08.c',
    kernel_prefix + 'm08_hx_c.hpp',
    'MGGA', 1e-15,
    { 'm08_a_0': '1.0000000e+00', 'm08_a_1': '-4.0661387e-01', 'm08_a_2': '-3.3232530e+00', 'm08_a_3': '1.5540980e+00', 'm08_a_4': '4.4248033e+01', 'm08_a_5': '-8.4351930e+01',
      'm08_a_6': '-1.1955581e+02', 'm08_a_7': '3.9147081e+02', 'm08_a_8': '1.8363851e+02', 'm08_a_9': '-6.3268223e+02', 'm08_a_10': '-1.1297403e+02', 'm08_a_11': '3.3629312e+02',
      'm08_b_0': '1.3812334e+00', 'm08_b_1': '-2.4683806e+00', 'm08_b_2': '-1.1901501e+01', 'm08_b_3': '-5.4112667e+01', 'm08_b_4': '1.0055846e+01', 'm08_b_5': '1.4800687e+02',
      'm08_b_6': '1.1561420e+02', 'm08_b_7': '2.5591815e+02', 'm08_b_8': '2.1320772e+02', 'm08_b_9': '-4.8412067e+02', 'm08_b_10': '-4.3430813e+02', 'm08_b_11': '5.6627964e+01'}
    ),
  'M08_SO_C' : GenMetaData( 'BuiltinM08_SO_C',
    libxc_prefix + 'mgga_exc/mgga_c_m08.c',
    kernel_prefix + 'm08_so_c.hpp',
    'MGGA', 1e-15,
    { 'm08_a_0': '1.0000000e+00', 'm08_a_1': '0.0000000e+00', 'm08_a_2': '-3.9980886e+00', 'm08_a_3': '1.2982340e+01', 'm08_a_4': '1.0117507e+02', 'm08_a_5': '-8.9541984e+01',
      'm08_a_6': '-3.5640242e+02', 'm08_a_7': '2.0698803e+02', 'm08_a_8': '4.6037780e+02', 'm08_a_9': '-2.4510559e+02', 'm08_a_10': '-1.9638425e+02', 'm08_a_11': '1.1881459e+02',
      'm08_b_0': '1.0000000e+00', 'm08_b_1': '-4.4117403e+00', 'm08_b_2': '-6.4128622e+00', 'm08_b_3': '4.7583635e+01', 'm08_b_4': '1.8630053e+02', 'm08_b_5': '-1.2800784e+02',
      'm08_b_6': '-5.5385258e+02', 'm08_b_7': '1.3873727e+02', 'm08_b_8': '4.1646537e+02', 'm08_b_9': '-2.6626577e+02', 'm08_b_10': '5.6676300e+01', 'm08_b_11': '3.1673746e+02'}
    ),
  'M11_C' : GenMetaData( 'BuiltinM11_C',
    libxc_prefix + 'mgga_exc/mgga_c_m08.c',
    kernel_prefix + 'm11_c.hpp',
    'MGGA', 1e-15,
    { 'm08_a_0': '1.0000000e+00', 'm08_a_1': '0.0000000e+00', 'm08_a_2': '-3.8933250e+00', 'm08_a_3': '-2.1688455e+00', 'm08_a_4': '9.3497200e+00', 'm08_a_5': '-1.9845140e+01',
      'm08_a_6': '2.3455253e+00', 'm08_a_7': '7.9246513e+01', 'm08_a_8': '9.6042757e+00', 'm08_a_9': '-6.7856719e+01', 'm08_a_10': '-9.1841067e+00', 'm08_a_11': '0.0000000e+00',
      'm08_b_0': '7.2239798e-01', 'm08_b_1': '4.3730564e-01', 'm08_b_2': '-1.6088809e+01', 'm08_b_3': '-6.5542437e+01', 'm08_b_4': '3.2057230e+01', 'm08_b_5': '1.8617888e+02',
      'm08_b_6': '2.0483468e+01', 'm08_b_7': '-7.0853739e+01', 'm08_b_8': '4.4483915e+01', 'm08_b_9': '-9.4484747e+01', 'm08_b_10': '-1.1459868e+02', 'm08_b_11': '0.0000000e+00'}
    ),
  'MN12_L_C' : GenMetaData( 'BuiltinMN12_L_C',
    libxc_prefix + 'mgga_exc/mgga_c_m08.c',
    kernel_prefix + 'mn12_l_c.hpp',
    'MGGA', 1e-15,
    { 'm08_a_0': '8.844610e-01', 'm08_a_1': '-2.202279e-01', 'm08_a_2': '5.701372e+00', 'm08_a_3': '-2.562378e+00', 'm08_a_4': '-9.646827e-01', 'm08_a_5': '1.982183e-01',
      'm08_a_6': '1.019976e+01', 'm08_a_7': '9.789352e-01', 'm08_a_8': '-1.512722e+00', 'm08_a_9': '0.000000e+00', 'm08_a_10': '0.000000e+00', 'm08_a_11': '0.000000e+00',
      'm08_b_0': '5.323948e-01', 'm08_b_1': '-5.831909e+00', 'm08_b_2': '3.882386e+00', 'm08_b_3': '5.878488e+00', 'm08_b_4': '1.493228e+01', 'm08_b_5': '-1.374636e+01',
      'm08_b_6': '-8.492327e+00', 'm08_b_7': '-2.486548e+00', 'm08_b_8': '-1.822346e+01', 'm08_b_9': '0.000000e+00', 'm08_b_10': '0.000000e+00', 'm08_b_11': '0.000000e+00'}
    ),
  'MN12_SX_C' : GenMetaData( 'BuiltinMN12_SX_C',
    libxc_prefix + 'mgga_exc/mgga_c_m08.c',
    kernel_prefix + 'mn12_sx_c.hpp',
    'MGGA', 1e-15,
    { 'm08_a_0': '7.171161e-01', 'm08_a_1': '-2.380914e+00', 'm08_a_2': '5.793565e+00', 'm08_a_3': '-1.243624e+00', 'm08_a_4': '1.364920e+01', 'm08_a_5': '-2.110812e+01',
      'm08_a_6': '-1.598767e+01', 'm08_a_7': '1.429208e+01', 'm08_a_8': '6.149191e+00', 'm08_a_9': '0.000000e+00', 'm08_a_10': '0.000000e+00', 'm08_a_11': '0.000000e+00',
      'm08_b_0': '4.663699e-01', 'm08_b_1': '-9.110685e+00', 'm08_b_2': '8.705051e+00', 'm08_b_3': '-1.813949e+00', 'm08_b_4': '-4.147211e-01', 'm08_b_5': '-1.021527e+01',
      'm08_b_6': '8.240270e-01', 'm08_b_7': '4.993815e+00', 'm08_b_8': '-2.563930e+01', 'm08_b_9': '0.000000e+00', 'm08_b_10': '0.000000e+00', 'm08_b_11': '0.000000e+00'}
    ),
  'MN15_L_C' : GenMetaData( 'BuiltinMN15_L_C',
    libxc_prefix + 'mgga_exc/mgga_c_m08.c',
    kernel_prefix + 'mn15_l_c.hpp',
    'MGGA', 1e-15,
    { 'm08_a_0': '0.952058087', 'm08_a_1': '-0.756954364', 'm08_a_2': '5.677396094', 'm08_a_3': '-5.017104782', 'm08_a_4': '-5.10654071', 'm08_a_5': '-4.812053335',
      'm08_a_6': '3.397640087', 'm08_a_7': '1.980041517', 'm08_a_8': '10.1231046', 'm08_a_9': '0.0', 'm08_a_10': '0.0', 'm08_a_11': '0.0',
      'm08_b_0': '0.819504932', 'm08_b_1': '-7.689358913', 'm08_b_2': '-0.70532663', 'm08_b_3': '-0.600096421', 'm08_b_4': '11.03332527', 'm08_b_5': '5.861969337',
      'm08_b_6': '8.913865465', 'm08_b_7': '5.74529876', 'm08_b_8': '4.254880837', 'm08_b_9': '0.0', 'm08_b_10': '0.0', 'm08_b_11': '0.0'}
    ),
  'MN15_C' : GenMetaData( 'BuiltinMN15_C',
    libxc_prefix + 'mgga_exc/mgga_c_m08.c',
    kernel_prefix + 'mn15_c.hpp',
    'MGGA', 1e-15,
    { 'm08_a_0': '1.093250748', 'm08_a_1': '-0.269735037', 'm08_a_2': '6.368997613', 'm08_a_3': '-0.245337101', 'm08_a_4': '-1.587103441', 'm08_a_5': '0.124698862',
      'm08_a_6': '1.605819855', 'm08_a_7': '0.466206031', 'm08_a_8': '3.484978654', 'm08_a_9': '0.0', 'm08_a_10': '0.0', 'm08_a_11': '0.0',
      'm08_b_0': '1.427424993', 'm08_b_1': '-3.57883682', 'm08_b_2': '7.398727547', 'm08_b_3': '3.927810559', 'm08_b_4': '2.789804639', 'm08_b_5': '4.988320462',
      'm08_b_6': '3.079464318', 'm08_b_7': '3.521636859', 'm08_b_8': '4.769671992', 'm08_b_9': '0.0', 'm08_b_10': '0.0', 'm08_b_11': '0.0'}
    ),
  'CF22D_C' : GenMetaData( 'BuiltinCF22D_C',
    libxc_prefix + 'mgga_exc/mgga_c_m08.c',
    kernel_prefix + 'cf22d_c.hpp',
    'MGGA', 1e-15,
    { 'm08_a_0': '0.873863376', 'm08_a_1': '0.078066142', 'm08_a_2': '6.576550257', 'm08_a_3': '-1.126030147', 'm08_a_4': '-3.244797887', 'm08_a_5': '-2.186090839',
      'm08_a_6': '-3.489135041', 'm08_a_7': '3.090689716', 'm08_a_8': '3.866592474', 'm08_a_9': '0.0', 'm08_a_10': '0.0', 'm08_a_11': '0.0',
      'm08_b_0': '0.828203832', 'm08_b_1': '-2.518707202', 'm08_b_2': '10.436806314', 'm08_b_3': '3.588267084', 'm08_b_4': '-5.789404145', 'm08_b_5': '3.353560215',
      'm08_b_6': '-2.432384384', 'm08_b_7': '-1.147183331', 'm08_b_8': '2.991316045', 'm08_b_9': '0.0', 'm08_b_10': '0.0', 'm08_b_11': '0.0'}
    ),

  'P86_C' : GenMetaData( 'BuiltinP86_C',
    libxc_prefix + 'gga_exc/gga_c_p86.c',
    kernel_prefix + 'p86_c.hpp',
    'GGA', 1e-15,
    { 'malpha': '0.023266', 'mbeta': '7.389e-6', 'mgamma': '8.723', 'mdelta': '0.472', 'aa': '0.001667', 
     'bb': '0.002568', 'ftilde': '0.19195'}
    ),

  'PW91_C' : GenMetaData( 'BuiltinPW91_C',
    libxc_prefix + 'gga_exc/gga_c_pw91.c',
    kernel_prefix + 'pw91_c.hpp',
    'GGA', 1e-12),
    
  'PBE_SOL_C' : GenMetaData( 'BuiltinPBE_SOL_C',
    libxc_prefix + 'gga_exc/gga_c_pbe.c',
    kernel_prefix + 'pbe_sol_c.hpp',
    'GGA', 1e-12, {'beta':'0.046', 'gamma':'0.031090690869654895034', 'BB':'1.'} 
    ),

  'BMK_C' : GenMetaData( 'BuiltinBMK_C',
    libxc_prefix + 'gga_exc/gga_c_bmk.c',
    kernel_prefix + 'bmk_c.hpp',
    'GGA', 1e-14,
    {'c_ss_0':'-2.19098', 'c_ss_1':'23.8939', 'c_ss_2':'-44.3303', 'c_ss_3':'22.5982', 'c_ss_4':'0.0',
      'c_ab_0':'1.22334', 'c_ab_1':'-3.4631', 'c_ab_2':'10.0731', 'c_ab_3':'-11.1974', 'c_ab_4':'0.0'}
      ),  
  'N12_C' : GenMetaData( 'BuiltinN12_C',
    libxc_prefix + 'gga_exc/gga_c_bmk.c',
    kernel_prefix + 'n12_c.hpp',
    'GGA', 1e-14,
    {'c_ss_0':'1.00000e+00', 'c_ss_1':'-5.53170e+00', 'c_ss_2':'3.07958e+01', 'c_ss_3':'-5.64196e+01', 'c_ss_4':'3.21250e+01',
      'c_ab_0':'1.00000e+00', 'c_ab_1':'3.24511e+00', 'c_ab_2':'-2.52893e+01', 'c_ab_3':'1.44407e+01', 'c_ab_4':'1.96870e+01'}
    ),
  'N12_SX_C' : GenMetaData( 'BuiltinN12_SX_C',
    libxc_prefix + 'gga_exc/gga_c_bmk.c',
    kernel_prefix + 'n12_sx_c.hpp',
    'GGA', 1e-14,
    {'c_ss_0':'2.63373e+00', 'c_ss_1':'-1.05450e+00', 'c_ss_2':'-7.29853e-01', 'c_ss_3':'4.94024e+00', 'c_ss_4':'-7.31760e+00',
      'c_ab_0':'8.33615e-01', 'c_ab_1':'3.24128e+00', 'c_ab_2':'-1.06407e+01', 'c_ab_3':'-1.60471e+01', 'c_ab_4':'2.51047e+01'}
    ),

  'P86VWN_FT_C' : GenMetaData( 'BuiltinP86VWN_FT_C',
    libxc_prefix + 'gga_exc/gga_c_p86vwn.c',
    kernel_prefix + 'p86vwn_ft_c.hpp',
    'GGA', 1e-15,
    { 'malpha': '0.023266', 'mbeta': '7.389e-6', 'mgamma': '8.723', 'mdelta': '0.472', 'aa': '0.001667', 
     'bb': '0.002568', 'ftilde': '1.7454151061251239789*0.11'}
    ),

  'SOGGA11_X_C' : GenMetaData( 'BuiltinSOGGA11_X_C',
    libxc_prefix + 'gga_exc/gga_c_sogga11.c',
    kernel_prefix + 'sogga11_x_c.hpp',
    'GGA', 1e-14,
    { 'sogga11_a_0': '0.50000', 'sogga11_a_1': '78.2439', 'sogga11_a_2': '25.7211', 'sogga11_a_3': '-13.8830', 'sogga11_a_4': '-9.87375', 'sogga11_a_5': '-14.1357',
      'sogga11_b_0': '0.50000', 'sogga11_b_1': '-79.2439', 'sogga11_b_2': '16.3725', 'sogga11_b_3': '2.08129', 'sogga11_b_4': '7.50769', 'sogga11_b_5': '-10.1861'}
    ),

  'TPSS_C' : GenMetaData( 'BuiltinTPSS_C',
    libxc_prefix + 'mgga_exc/mgga_c_tpss.c',
    kernel_prefix + 'tpss_c.hpp',
    'MGGA', 1e-15,
    { 'beta': '0.06672455060314922', 'd': '2.8', 'C0_c_0': '0.53', 'C0_c_1': '0.87', 'C0_c_2': '0.50', 'C0_c_3': '2.26'}
    ),
  'BC95_C' : GenMetaData( 'BuiltinBC95_C',
    libxc_prefix + 'mgga_exc/mgga_c_bc95.c',
    kernel_prefix + 'bc95_c.hpp',
    'MGGA', 1e-14,
    { 'css': '0.038', 'copp': '0.0031'}
    ),

  'revTPSS_C' : GenMetaData( 'BuiltinRevTPSS_C',
    libxc_prefix + 'mgga_exc/mgga_c_revtpss.c',
    kernel_prefix + 'revtpss_c.hpp',
    'MGGA', 1e-13,
    { 'd': '2.8', 'C0_c_0': '0.59', 'C0_c_1': '0.9269', 'C0_c_2': '0.6225', 'C0_c_3': '2.1540'}
    ),

  'RSCAN_C' : GenMetaData( 'BuiltinRSCAN_C',
    libxc_prefix + 'mgga_exc/mgga_c_rscan.c',
    kernel_prefix + 'rscan_c.hpp',
    'MGGA', 1e-15),

  'PW91_X': GenMetaData( 'BuiltinPW91_X',
    libxc_prefix + 'gga_exc/gga_x_pw91.c',
    kernel_prefix + 'pw91_x.hpp',
    'GGA', 1e-15,
    { 'a': '0.19645', 'b': '7.7956', 'c': '0.2743', 'd': '-0.1508', 'f': '0.004', 'alpha': '100.0', 'expo': '4.0'}
    ),

  'MPW91_X': GenMetaData( 'BuiltinMPW91_X',
    libxc_prefix + 'gga_exc/gga_x_pw91.c',
    kernel_prefix + 'mpw91_x.hpp',
    'GGA', 1e-15,
    { 'a': '6.0*0.00426/constants::X2S', 'b': '1.0/constants::X2S', 
      'c': '0.00426/(constants::X_FACTOR_C*constants::X2S*constants::X2S)', 
      'd': '-(0.00426 - 5.0*0.0003780762333399851)/(constants::X_FACTOR_C*constants::X2S*constants::X2S)', 
      'f': '1.0e-6/(constants::X_FACTOR_C*0.00048120394750740677)', 
      'alpha': '100.0', 'expo': '3.72'}
    ),
    
  'OPTX_X': GenMetaData( 'BuiltinOPTX_X',
    libxc_prefix + 'gga_exc/gga_x_optx.c',
    kernel_prefix + 'optx_x.hpp',
    'GGA', 1e-15,
    { 'a': '1.05151', 'b': '1.43169/constants::X_FACTOR_C', 'gamma': '0.006'}
    ),
    
  'RPBE_X': GenMetaData( 'BuiltinRPBE_X',
    libxc_prefix + 'gga_exc/gga_x_rpbe.c',
    kernel_prefix + 'rpbe_x.hpp',
    'GGA', 1e-15,
    { 'rpbe_kappa': '0.8040', 'rpbe_mu': '0.2195149727645171'}
    ),
    
  'SOGGA11_X_X': GenMetaData( 'BuiltinSOGGA11_X_X',
    libxc_prefix + 'gga_exc/gga_x_sogga11.c',
    kernel_prefix + 'sogga11_x_x.hpp',
    'GGA', 1e-15,
    { 'kappa': '0.552', 'mu': '0.1234567901234567901234567901234567901235',
      'a_0': '0.29925', 'a_1': '3.21638', 'a_2': '-3.55605', 'a_3': '7.65852', 'a_4': '-11.2830', 'a_5': '5.25813',
      'b_0': '0.29925', 'b_1': '-2.88595', 'b_2': '3.23617', 'b_3': '-2.45393', 'b_4': '-3.75495', 'b_5': '3.96613',
      'cx': '0.4015'}
    ),
    
  'PW86_X': GenMetaData( 'BuiltinPW86_X',
    libxc_prefix + 'gga_exc/gga_x_pw86.c',
    kernel_prefix + 'pw86_x.hpp',
    'GGA', 1e-15,
    { 'aa': '1.296', 'bb': '14.0', 'cc': '0.2'}
    ),

  'mBEEF_X': GenMetaData( 'BuiltinMBEEF_X',
    libxc_prefix + 'mgga_exc/mgga_x_mbeef.c',
    kernel_prefix + 'mbeef_x.hpp',
    'MGGA', 1e-15,
  ),

  'RSCAN_X': GenMetaData( 'BuiltinRSCAN_X',
    libxc_prefix + 'mgga_exc/mgga_x_rscan.c',
    kernel_prefix + 'rscan_x.hpp',
    'MGGA', 1e-11,
    { 'c2': '0.8', 'd': '1.24', 'k1': '0.065', 'taur': '1.0e-4', 'alphar': '1.0e-3'}
    ),

  'BMK_X': GenMetaData( 'BuiltinBMK_X',
    libxc_prefix + 'mgga_exc/mgga_x_tau_hcth.c',
    kernel_prefix + 'bmk_x.hpp',
    'MGGA', 1e-15,
    { 'cx_local_0': '0.474302', 'cx_local_1': '2.77701', 'cx_local_2': '-11.4230', 'cx_local_3': '11.7167',
      'cx_nlocal_0': '-0.192212', 'cx_nlocal_1': '4.73936', 'cx_nlocal_2': '-26.6188', 'cx_nlocal_3': '22.4891', 'ax': '0.42'}
    ),

  'M08_HX_X': GenMetaData( 'BuiltinM08_HX_X',
    libxc_prefix + 'mgga_exc/mgga_x_m08.c',
    kernel_prefix + 'm08_hx_x.hpp',
    'MGGA', 1e-15,
    { 'a_0': '1.3340172e+00', 'a_1': '-9.4751087e+00', 'a_2': '-1.2541893e+01', 'a_3': '9.1369974e+00', 'a_4': '3.4717204e+01', 'a_5': '5.8831807e+01',
      'a_6': '7.1369574e+01', 'a_7': '2.3312961e+01', 'a_8': '4.8314679e+00', 'a_9': '-6.5044167e+00', 'a_10': '-1.4058265e+01', 'a_11': '1.2880570e+01',
      'b_0': '-8.5631823e-01', 'b_1': '9.2810354e+00', 'b_2': '1.2260749e+01', 'b_3': '-5.5189665e+00', 'b_4': '-3.5534989e+01', 'b_5': '-8.2049996e+01',
      'b_6': '-6.8586558e+01', 'b_7': '3.6085694e+01', 'b_8': '-9.3740983e+00', 'b_9': '-5.9731688e+01', 'b_10': '1.6587868e+01', 'b_11': '1.3993203e+01',
      'ax': '0.5223'}
    ),

  'M08_SO_X': GenMetaData( 'BuiltinM08_SO_X',
    libxc_prefix + 'mgga_exc/mgga_x_m08.c',
    kernel_prefix + 'm08_so_x.hpp',
    'MGGA', 1e-15,
    { 'a_0': '-3.4888428e-01', 'a_1': '-5.8157416e+00', 'a_2': '3.7550810e+01', 'a_3': '6.3727406e+01', 'a_4': '-5.3742313e+01', 'a_5': '-9.8595529e+01',
      'a_6': '1.6282216e+01', 'a_7': '1.7513468e+01', 'a_8': '-6.7627553e+00', 'a_9': '1.1106658e+01', 'a_10': '1.5663545e+00', 'a_11': '8.7603470e+00',
      'b_0': '7.8098428e-01', 'b_1': '5.4538178e+00', 'b_2': '-3.7853348e+01', 'b_3': '-6.2295080e+01', 'b_4': '4.6713254e+01', 'b_5': '8.7321376e+01',
      'b_6': '1.6053446e+01', 'b_7': '2.0126920e+01', 'b_8': '-4.0343695e+01', 'b_9': '-5.8577565e+01', 'b_10': '2.0890272e+01', 'b_11': '1.0946903e+01',
      'ax': '0.5679'}
    ),

  'MN12_L_X': GenMetaData( 'BuiltinMN12_L_X',
    libxc_prefix + 'mgga_exc/mgga_x_mn12.c',
    kernel_prefix + 'mn12_l_x.hpp',
    'MGGA', 1e-15,
    { 'c_0': '6.735981e-01', 'c_1': '-2.270598e+00', 'c_2': '-2.613712e+00', 'c_3': '3.993609e+00', 'c_4': '4.635575e+00', 'c_5': '1.250676e+00',
      'c_6': '8.444920e-01', 'c_7': '-1.301173e+01', 'c_8': '-1.777730e+01', 'c_9': '-4.627211e+00', 'c_10': '5.976605e+00', 'c_11': '1.142897e+00',
      'c_12': '-2.040226e+01', 'c_13': '-2.382843e+01', 'c_14': '7.119109e+00', 'c_15': '-2.335726e+01', 'c_16': '-1.622633e+01', 'c_17': '1.482732e+01',
      'c_18': '1.449285e+00', 'c_19': '1.020598e+01', 'c_20': '4.407450e+00', 'c_21': '-2.008193e+01', 'c_22': '-1.253561e+01', 'c_23': '-5.435031e+00',
      'c_24': '1.656736e+01', 'c_25': '2.000229e+01', 'c_26': '-2.513105e+00', 'c_27': '9.658436e+00', 'c_28': '-3.825281e+00', 'c_29': '-2.500000e+01',
      'c_30': '-2.070080e+00', 'c_31': '-9.951913e+00', 'c_32': '8.731211e-01', 'c_33': '2.210891e+01', 'c_34': '8.822633e+00', 'c_35': '2.499949e+01',
      'c_36': '2.500000e+01', 'c_37': '6.851693e-01', 'c_38': '-7.406948e-02', 'c_39': '-6.788000e-01'}
    ),

  'MN15_L_X': GenMetaData( 'BuiltinMN15_L_X',
    libxc_prefix + 'mgga_exc/mgga_x_mn12.c',
    kernel_prefix + 'mn15_l_x.hpp',
    'MGGA', 1e-15,
    { 'c_0': '0.670864162', 'c_1': '-0.822003903', 'c_2': '-1.022407046', 'c_3': '1.689460986', 'c_4': '-0.00562032', 'c_5': '-0.110293849',
      'c_6': '0.972245178', 'c_7': '-6.697641991', 'c_8': '-4.322814495', 'c_9': '-6.786641376', 'c_10': '-5.687461462', 'c_11': '9.419643818',
      'c_12': '11.83939406', 'c_13': '5.086951311', 'c_14': '4.302369948', 'c_15': '-8.07344065', 'c_16': '2.429988978', 'c_17': '11.09485698',
      'c_18': '1.247333909', 'c_19': '3.700485291', 'c_20': '0.867791614', 'c_21': '-0.591190518', 'c_22': '-0.295305435', 'c_23': '-5.825759145',
      'c_24': '2.537532196', 'c_25': '3.143390933', 'c_26': '2.939126332', 'c_27': '0.599342114', 'c_28': '2.241702738', 'c_29': '2.035713838',
      'c_30': '-1.525344043', 'c_31': '-2.325875691', 'c_32': '1.141940663', 'c_33': '-1.563165026', 'c_34': '7.882032871', 'c_35': '11.93400684',
      'c_36': '9.852928303', 'c_37': '0.584030245', 'c_38': '-0.720941131', 'c_39': '-2.836037078'}
    ),

  'MN15_X': GenMetaData( 'BuiltinMN15_X',
    libxc_prefix + 'mgga_exc/mgga_x_mn12.c',
    kernel_prefix + 'mn15_x.hpp',
    'MGGA', 1e-15,
    { 'c_0': '0.073852235', 'c_1': '-0.839976156', 'c_2': '-3.082660125', 'c_3': '-1.02881285', 'c_4': '-0.811697255', 'c_5': '-0.063404387',
      'c_6': '2.54805518', 'c_7': '-5.031578906', 'c_8': '0.31702159', 'c_9': '2.981868205', 'c_10': '-0.749503735',
      'c_11': '0.231825661', 'c_12': '1.261961411', 'c_13': '1.665920815', 'c_14': '7.483304941', 'c_15': '-2.544245723', 'c_16': '1.384720031',
      'c_17': '6.902569885', 'c_18': '1.657399451', 'c_19': '2.98526709', 'c_20': '6.89391326', 'c_21': '2.489813993', 'c_22': '1.454724691',
      'c_23': '-5.054324071', 'c_24': '2.35273334', 'c_25': '1.299104132', 'c_26': '1.203168217', 'c_27': '0.121595877', 'c_28': '8.048348238',
      'c_29': '21.91203659', 'c_30': '-1.852335832', 'c_31': '-3.4722735', 'c_32': '-1.564591493', 'c_33': '-2.29578769', 'c_34': '3.666482991',
      'c_35': '10.87074639', 'c_36': '9.696691388', 'c_37': '0.630701064', 'c_38': '-0.505825216', 'c_39': '-3.562354535', 'ax': '0.44'}
    ),

  'CF22D_X': GenMetaData( 'BuiltinCF22D_X',
    libxc_prefix + 'mgga_exc/mgga_x_mn12.c',
    kernel_prefix + 'cf22d_x.hpp',
    'MGGA', 1e-15,
    { 'c_0': '0.24416116', 'c_1': '-0.389728151', 'c_2': '-1.829675858', 'c_3': '1.396044771', 'c_4': '2.315047133', 'c_5': '0.397552547',
      'c_6': '1.082144406', 'c_7': '-7.894560034', 'c_8': '-3.656253030', 'c_9': '2.574496508', 'c_10': '4.031038406',
      'c_11': '-3.931389433', 'c_12': '0.333519075', 'c_13': '-3.032270318', 'c_14': '3.673752289', 'c_15': '3.005997956', 'c_16': '-6.463733874',
      'c_17': '-4.596755225', 'c_18': '0.964839180', 'c_19': '0.363791944', 'c_20': '1.646506623', 'c_21': '-3.504641550', 'c_22': '-3.922228074',
      'c_23': '0.843718076', 'c_24': '10.779373313', 'c_25': '2.293612669', 'c_26': '7.088363286', 'c_27': '2.598770741', 'c_28': '-0.088522116',
      'c_29': '7.180809030', 'c_30': '-1.017514009', 'c_31': '1.735020310', 'c_32': '3.499241561', 'c_33': '0.922224945', 'c_34': '-2.212903920',
      'c_35': '0.243080429', 'c_36': '17.306321840', 'c_37': '0.311402396', 'c_38': '-3.257126009', 'c_39': '-3.372399742', 'ax': '0.462806'}
    ),
    
  'MN12_SX_X': GenMetaData( 'BuiltinMN12_SX_X',
    libxc_prefix + 'mgga_exc/mgga_x_mn12.c',
    kernel_prefix + 'mn12_sx_x.hpp',
    'MGGA', 1e-15,
    { 'c_0': '5.226556e-01', 'c_1': '-2.681208e-01', 'c_2': '-4.670705e+00', 'c_3': '3.067320e+00', 'c_4': '4.095370e+00', 'c_5': '2.653023e+00',
      'c_6': '5.165969e-01', 'c_7': '-2.035442e+01', 'c_8': '-9.946472e+00', 'c_9': '2.938637e+00', 'c_10': '1.131100e+01', 'c_11': '4.752452e+00',
      'c_12': '-3.061331e+00', 'c_13': '-2.523173e+01', 'c_14': '1.710903e+01', 'c_15': '-2.357480e+01', 'c_16': '-2.727754e+01', 'c_17': '1.603291e+01',
      'c_18': '1.842503e+00', 'c_19': '1.927120e+00', 'c_20': '1.107987e+01', 'c_21': '-1.182087e+01', 'c_22': '-1.117768e+01', 'c_23': '-5.821000e+00',
      'c_24': '2.266545e+01', 'c_25': '8.246708e+00', 'c_26': '-4.778364e+00', 'c_27': '5.329122e-01', 'c_28': '-6.666755e+00', 'c_29': '1.671429e+00',
      'c_30': '-3.311409e+00', 'c_31': '3.415913e-01', 'c_32': '-6.413076e+00', 'c_33': '1.038584e+01', 'c_34': '9.026277e+00', 'c_35': '1.929689e+01',
      'c_36': '2.669232e+01', 'c_37': '1.517278e+00', 'c_38': '-3.442503e+00', 'c_39': '1.100161e+00', 'ax': '0.0', 'sx': '0.25', 'omega': '0.11'}
    ),

  'M11_X': GenMetaData( 'BuiltinM11_X',
    libxc_prefix + 'mgga_exc/mgga_x_m11.c',
    kernel_prefix + 'm11_x.hpp',
    'MGGA', 1e-11,
    { 'a_0': '-0.18399900e+00', 'a_1': '-1.39046703e+01', 'a_2': '1.18206837e+01', 'a_3': '3.10098465e+01', 'a_4': '-5.19625696e+01', 'a_5': '1.55750312e+01',
      'a_6': '-6.94775730e+00', 'a_7': '-1.58465014e+02', 'a_8': '-1.48447565e+00', 'a_9': '5.51042124e+01', 'a_10': '-1.34714184e+01', 'a_11': '0.00000000e+00',
      'b_0': '0.75599900e+00', 'b_1': '1.37137944e+01', 'b_2': '-1.27998304e+01', 'b_3': '-2.93428814e+01', 'b_4': '5.91075674e+01', 'b_5': '-2.27604866e+01',
      'b_6': '-1.02769340e+01', 'b_7': '1.64752731e+02', 'b_8': '1.85349258e+01', 'b_9': '-5.56825639e+01', 'b_10': '7.47980859e+00', 'b_11': '0.00000000e+00',
      'alpha': '1.0', 'beta': '-0.572', 'omega': '0.25'}
    ),

  'wB97_XC': GenMetaData( 'BuiltinWB97_XC',
    libxc_prefix + 'gga_exc/hyb_gga_xc_wb97.c',
    kernel_prefix + 'wb97_xc.hpp',
    'GGA', 1e-14,
    { 'c_x_0': '1.00000e+00', 'c_x_1': '1.13116e+00', 'c_x_2': '-2.74915e+00', 'c_x_3': '1.20900e+01', 'c_x_4': '-5.71642e+00',
      'c_ss_0': '1.00000e+00', 'c_ss_1': '-2.55352e+00', 'c_ss_2': '1.18926e+01', 'c_ss_3': '-2.69452e+01', 'c_ss_4': '1.70927e+01',
      'c_ab_0': '1.00000e+00', 'c_ab_1': '3.99051e+00', 'c_ab_2': '-1.70066e+01', 'c_ab_3': '1.07292e+00', 'c_ab_4': '8.88211e+00',
      'alpha': '1.0', 'beta': '-1.0', 'omega': '0.4'}
    ),

  'wB97X_XC': GenMetaData( 'BuiltinWB97X_XC',
    libxc_prefix + 'gga_exc/hyb_gga_xc_wb97.c',
    kernel_prefix + 'wb97x_xc.hpp',
    'GGA', 1e-14,
    { 'c_x_0': '8.42294e-01', 'c_x_1': '7.26479e-01', 'c_x_2': '1.04760e+00', 'c_x_3': '-5.70635e+00', 'c_x_4': '1.32794e+01',
      'c_ss_0': '1.00000e+00', 'c_ss_1': '-4.33879e+00', 'c_ss_2': '1.82308e+01', 'c_ss_3': '-3.17430e+01', 'c_ss_4': '1.72901e+01',
      'c_ab_0': '1.00000e+00', 'c_ab_1': '2.37031e+00', 'c_ab_2': '-1.13995e+01', 'c_ab_3': '6.58405e+00', 'c_ab_4': '-3.78132e+00',
      'alpha': '1.0', 'beta': '-0.842294', 'omega': '0.3'}
    ),

  'wB97X_V_XC': GenMetaData( 'BuiltinWB97X_V_XC',
    libxc_prefix + 'gga_exc/hyb_gga_xc_wb97.c',
    kernel_prefix + 'wb97x_v_xc.hpp',
    'GGA', 1e-14,
    { 'c_x_0': '0.833', 'c_x_1': '0.603', 'c_x_2': '1.194', 'c_x_3': '0.0', 'c_x_4': '0.0',
      'c_ss_0': '0.556', 'c_ss_1': '-0.257', 'c_ss_2': '0.0', 'c_ss_3': '0.0', 'c_ss_4': '0.0',
      'c_ab_0': '1.219', 'c_ab_1': '-1.850', 'c_ab_2': '0.0', 'c_ab_3': '0.0', 'c_ab_4': '0.0',
      'alpha': '1.0', 'beta': '-0.833', 'omega': '0.3'}
    ),

  'wB97X_D_XC': GenMetaData( 'BuiltinWB97X_D_XC',
    libxc_prefix + 'gga_exc/hyb_gga_xc_wb97.c',
    kernel_prefix + 'wb97x_d_xc.hpp',
    'GGA', 1e-14,
    { 'c_x_0': '7.77964e-01', 'c_x_1': '6.61160e-01', 'c_x_2': '5.74541e-01', 'c_x_3': '-5.25671e+00', 'c_x_4': '1.16386e+01',
      'c_ss_0': '1.00000e+00', 'c_ss_1': '-6.90539e+00', 'c_ss_2': '3.13343e+01', 'c_ss_3': '-5.10533e+01', 'c_ss_4': '2.64423e+01',
      'c_ab_0': '1.00000e+00', 'c_ab_1': '1.79413e+00', 'c_ab_2': '-1.20477e+01', 'c_ab_3': '1.40847e+01', 'c_ab_4': '-8.50809e+00',
      'alpha': '1.0', 'beta': '-0.777964', 'omega': '0.2'}
    ),

  'wB97X_D3_XC': GenMetaData( 'BuiltinWB97X_D3_XC',
    libxc_prefix + 'gga_exc/hyb_gga_xc_wb97.c',
    kernel_prefix + 'wb97x_d3_xc.hpp',
    'GGA', 1e-14,
    { 'c_x_0': '0.804272', 'c_x_1': '0.698900', 'c_x_2': '0.508940', 'c_x_3': '-3.744903', 'c_x_4': '10.060790',
      'c_ss_0': '1.00000', 'c_ss_1': '-4.868902', 'c_ss_2': '21.295726', 'c_ss_3': '-36.020866', 'c_ss_4': '19.177018',
      'c_ab_0': '1.00000', 'c_ab_1': '2.433266', 'c_ab_2': '-15.446008', 'c_ab_3': '17.644390', 'c_ab_4': '-8.879494',
      'alpha': '1.0', 'beta': '-0.804272', 'omega': '0.25'}
    ),

  'HJS_PBE_X': GenMetaData( 'BuiltinHJS_PBE_X',
    libxc_prefix + 'gga_exc/gga_x_hjs.c',
    kernel_prefix + 'hjs_pbe_x.hpp',
    'GGA', 5e-12,
    { 'a_0': '0.0159941', 'a_1': '0.0852995', 'a_2': '-0.160368', 'a_3': '0.152645', 'a_4': '-0.0971263', 'a_5': '0.0422061',
      'b_0': '5.33319', 'b_1': '-12.4780', 'b_2': '11.0988', 'b_3': '-5.11013', 'b_4': '1.71468', 
      'b_5': '-0.610380', 'b_6': '0.307555', 'b_7': '-0.0770547', 'b_8': '0.0334840', 'omega': '0.11'}
    ),

  'LRCwPBE_HJS_PBE_X': GenMetaData( 'BuiltinLRCwPBE_HJS_PBE_X',
    libxc_prefix + 'gga_exc/gga_x_hjs.c',
    kernel_prefix + 'lrcwpbe_hjs_pbe_x.hpp',
    'GGA', 5e-12,
    { 'a_0': '0.0159941', 'a_1': '0.0852995', 'a_2': '-0.160368', 'a_3': '0.152645', 'a_4': '-0.0971263', 'a_5': '0.0422061',
      'b_0': '5.33319', 'b_1': '-12.4780', 'b_2': '11.0988', 'b_3': '-5.11013', 'b_4': '1.71468', 
      'b_5': '-0.610380', 'b_6': '0.307555', 'b_7': '-0.0770547', 'b_8': '0.0334840', 'omega': '0.3'}
    ),

  'LRCwPBEh_HJS_PBE_X': GenMetaData( 'BuiltinLRCwPBEh_HJS_PBE_X',
    libxc_prefix + 'gga_exc/gga_x_hjs.c',
    kernel_prefix + 'lrcwpbeh_hjs_pbe_x.hpp',
    'GGA', 5e-12,
    { 'a_0': '0.0159941', 'a_1': '0.0852995', 'a_2': '-0.160368', 'a_3': '0.152645', 'a_4': '-0.0971263', 'a_5': '0.0422061',
      'b_0': '5.33319', 'b_1': '-12.4780', 'b_2': '11.0988', 'b_3': '-5.11013', 'b_4': '1.71468', 
      'b_5': '-0.610380', 'b_6': '0.307555', 'b_7': '-0.0770547', 'b_8': '0.0334840', 'omega': '0.2'}
    ),                              
                                  
  'wPBEh_X_default0': GenMetaData( 'BuiltinWPBEh_X_default0',
    libxc_prefix + 'gga_exc/gga_x_wpbeh.c',
    kernel_prefix + 'wpbeh_x_default0.hpp',
    'GGA', 1e-14,
    { 'omega': '0.0'}
    ),

  'LCwPBE_wPBEh_X': GenMetaData( 'BuiltinLCwPBE_wPBEh_X',
    libxc_prefix + 'gga_exc/gga_x_wpbeh.c',
    kernel_prefix + 'lcwpbe_wpbeh_x.hpp',
    'GGA', 5e-12,
    { 'a_0': '0.0159941', 'a_1': '0.0852995', 'a_2': '-0.160368', 'a_3': '0.152645', 'a_4': '-0.0971263', 'a_5': '0.0422061',
      'b_0': '5.33319', 'b_1': '-12.4780', 'b_2': '11.0988', 'b_3': '-5.11013', 'b_4': '1.71468', 
      'b_5': '-0.610380', 'b_6': '0.307555', 'b_7': '-0.0770547', 'b_8': '0.0334840', 'omega': '0.4'}
    ),

  'HSE03_wPBEh_X': GenMetaData( 'BuiltinHSE03_wPBEh_X',
    libxc_prefix + 'gga_exc/gga_x_wpbeh.c',
    kernel_prefix + 'hse03_wpbeh_x.hpp',
    'GGA', 1e-14,
    { 'omega': '0.188988157484230974715081591092'}
    ),

  'HSE06_wPBEh_X': GenMetaData( 'BuiltinHSE06_wPBEh_X',
    libxc_prefix + 'gga_exc/gga_x_wpbeh.c',
    kernel_prefix + 'hse06_wpbeh_x.hpp',
    'GGA', 1e-14,
    { 'omega': '0.11'}
    ),



  #'R2SCANL_X' : GenMetaData( 'BuiltinR2SCANL_X',
  #  libxc_prefix + 'mgga_exc/mgga_x_r2scanl.c',
  #  kernel_prefix + 'r2scanl_x.hpp',
  #  'MGGA', 1e-11,
  #  {'c1': '0.667',
  #   'c2': '0.8',
  #   'd': '1.24',
  #   'k1': '0.065',
  #   'eta': '0.001',
  #   'dp2': '0.361',
  #   'a': '1.784720',
  #   'b': '0.258304'}
  #  ),

  #'R2SCANL_C' : GenMetaData( 'BuiltinR2SCANL_C',
  #  libxc_prefix + 'mgga_exc/mgga_c_r2scanl.c',
  #  kernel_prefix + 'r2scanl_c.hpp',
  #  'MGGA', 1e-11,
  #  {'eta': '0.001',
  #   'a': '1.784720',
  #   'b': '0.258304'}
  #  )

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
  static constexpr bool needs_laplacian = {5};
  static constexpr bool is_kedf = {6};
  static constexpr bool is_epc  = {7};

  static constexpr double dens_tol  = {8};
  static constexpr double zeta_tol  = 1e-15;
  static constexpr double sigma_tol  = {9};
  static constexpr double tau_tol = is_kedf ? 0.0 : 1e-20;

"""

lda_exc_args_unpolar     = "double rho, double& eps"
lda_exc_vxc_args_unpolar = "double rho, double& eps, double& vrho" 
gga_exc_args_unpolar     = "double rho, double sigma, double& eps"
gga_exc_vxc_args_unpolar = "double rho, double sigma, double& eps, double& vrho, double& vsigma"
mgga_exc_args_unpolar    = "double rho, double sigma, double lapl, double tau, double& eps"
mgga_exc_vxc_args_unpolar = "double rho, double sigma, double lapl, double tau, double& eps, double& vrho, double& vsigma, double& vlapl, double& vtau" 

lda_exc_args_polar     = "double rho_a, double rho_b, double& eps"
lda_exc_vxc_args_polar = "double rho_a, double rho_b, double& eps, double& vrho_a, double& vrho_b" 
gga_exc_args_polar     = "double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps"
gga_exc_vxc_args_polar = "double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb" 
mgga_exc_args_polar    = "double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double lapl_a, double lapl_b, double tau_a, double tau_b, double& eps"
mgga_exc_vxc_args_polar = "double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double lapl_a, double lapl_b, double tau_a, double tau_b, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb, double& vlapl_a, double& vlapl_b, double& vtau_a, double& vtau_b"

# Define VXC arguments (no energy output)
lda_vxc_args_unpolar = "double rho, double& vrho"
gga_vxc_args_unpolar = "double rho, double sigma, double& vrho, double& vsigma"
mgga_vxc_args_unpolar = "double rho, double sigma, double lapl, double tau, double& vrho, double& vsigma, double& vlapl, double& vtau"

lda_vxc_args_polar = "double rho_a, double rho_b, double& vrho_a, double& vrho_b"
gga_vxc_args_polar = "double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb"
mgga_vxc_args_polar = "double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double lapl_a, double lapl_b, double tau_a, double tau_b, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb, double& vlapl_a, double& vlapl_b, double& vtau_a, double& vtau_b"

# Create dictionaries for VXC argument maps
vxc_args_unpolar = {
  'LDA' : lda_vxc_args_unpolar,
  'GGA' : gga_vxc_args_unpolar,
  'MGGA': mgga_vxc_args_unpolar
}
vxc_args_polar = {
  'LDA' : lda_vxc_args_polar,
  'GGA' : gga_vxc_args_polar,
  'MGGA': mgga_vxc_args_polar
}

# LDA FXC argument definitions
lda_fxc_args_unpolar = "double rho, double& v2rho2"
lda_vxc_fxc_args_unpolar = "double rho, double& vrho, double& v2rho2"
lda_fxc_args_polar = "double rho_a, double rho_b, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb"
lda_vxc_fxc_args_polar = "double rho_a, double rho_b, double& vrho_a, double& vrho_b, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb"

# GGA FXC argument definitions
gga_fxc_args_unpolar = "double rho, double sigma, double& v2rho2, double& v2rhosigma, double& v2sigma2"
gga_vxc_fxc_args_unpolar = "double rho, double sigma, double& vrho, double& vsigma, double& v2rho2, double& v2rhosigma, double& v2sigma2"
gga_fxc_args_polar = "double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb, double& v2rhosigma_a_aa, double& v2rhosigma_a_ab, double& v2rhosigma_a_bb, double& v2rhosigma_b_aa, double& v2rhosigma_b_ab, double& v2rhosigma_b_bb, double& v2sigma2_aa_aa, double& v2sigma2_aa_ab, double& v2sigma2_aa_bb, double& v2sigma2_ab_ab, double& v2sigma2_ab_bb, double& v2sigma2_bb_bb"
gga_vxc_fxc_args_polar = "double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb, double& v2rhosigma_a_aa, double& v2rhosigma_a_ab, double& v2rhosigma_a_bb, double& v2rhosigma_b_aa, double& v2rhosigma_b_ab, double& v2rhosigma_b_bb, double& v2sigma2_aa_aa, double& v2sigma2_aa_ab, double& v2sigma2_aa_bb, double& v2sigma2_ab_ab, double& v2sigma2_ab_bb, double& v2sigma2_bb_bb"

# MGGA FXC argument definitions
mgga_fxc_args_unpolar = "double rho, double sigma, double lapl, double tau, double& v2rho2, double& v2rhosigma, double& v2rholapl, double& v2rhotau, double& v2sigma2, double& v2sigmalapl, double& v2sigmatau, double& v2lapl2, double& v2lapltau, double& v2tau2"
mgga_vxc_fxc_args_unpolar = "double rho, double sigma, double lapl, double tau, double& vrho, double& vsigma, double& vlapl, double& vtau, double& v2rho2, double& v2rhosigma, double& v2rholapl, double& v2rhotau, double& v2sigma2, double& v2sigmalapl, double& v2sigmatau, double& v2lapl2, double& v2lapltau, double& v2tau2"
mgga_fxc_args_polar = "double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double lapl_a, double lapl_b, double tau_a, double tau_b, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb, double& v2rhosigma_a_aa, double& v2rhosigma_a_ab, double& v2rhosigma_a_bb, double& v2rhosigma_b_aa, double& v2rhosigma_b_ab, double& v2rhosigma_b_bb, double& v2rholapl_a_a, double& v2rholapl_a_b, double& v2rholapl_b_a, double& v2rholapl_b_b, double& v2rhotau_a_a, double& v2rhotau_a_b, double& v2rhotau_b_a, double& v2rhotau_b_b, double& v2sigma2_aa_aa, double& v2sigma2_aa_ab, double& v2sigma2_aa_bb, double& v2sigma2_ab_ab, double& v2sigma2_ab_bb, double& v2sigma2_bb_bb, double& v2sigmalapl_aa_a, double& v2sigmalapl_aa_b, double& v2sigmalapl_ab_a, double& v2sigmalapl_ab_b, double& v2sigmalapl_bb_a, double& v2sigmalapl_bb_b, double& v2sigmatau_aa_a, double& v2sigmatau_aa_b, double& v2sigmatau_ab_a, double& v2sigmatau_ab_b, double& v2sigmatau_bb_a, double& v2sigmatau_bb_b, double& v2lapl2_aa, double& v2lapl2_ab, double& v2lapl2_bb, double& v2lapltau_a_a, double& v2lapltau_a_b, double& v2lapltau_b_a, double& v2lapltau_b_b, double& v2tau2_aa, double& v2tau2_ab, double& v2tau2_bb"
mgga_vxc_fxc_args_polar = "double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double lapl_a, double lapl_b, double tau_a, double tau_b, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb, double& vlapl_a, double& vlapl_b, double& vtau_a, double& vtau_b, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb, double& v2rhosigma_a_aa, double& v2rhosigma_a_ab, double& v2rhosigma_a_bb, double& v2rhosigma_b_aa, double& v2rhosigma_b_ab, double& v2rhosigma_b_bb, double& v2rholapl_a_a, double& v2rholapl_a_b, double& v2rholapl_b_a, double& v2rholapl_b_b, double& v2rhotau_a_a, double& v2rhotau_a_b, double& v2rhotau_b_a, double& v2rhotau_b_b, double& v2sigma2_aa_aa, double& v2sigma2_aa_ab, double& v2sigma2_aa_bb, double& v2sigma2_ab_ab, double& v2sigma2_ab_bb, double& v2sigma2_bb_bb, double& v2sigmalapl_aa_a, double& v2sigmalapl_aa_b, double& v2sigmalapl_ab_a, double& v2sigmalapl_ab_b, double& v2sigmalapl_bb_a, double& v2sigmalapl_bb_b, double& v2sigmatau_aa_a, double& v2sigmatau_aa_b, double& v2sigmatau_ab_a, double& v2sigmatau_ab_b, double& v2sigmatau_bb_a, double& v2sigmatau_bb_b, double& v2lapl2_aa, double& v2lapl2_ab, double& v2lapl2_bb, double& v2lapltau_a_a, double& v2lapltau_a_b, double& v2lapltau_b_a, double& v2lapltau_b_b, double& v2tau2_aa, double& v2tau2_ab, double& v2tau2_bb"

# Create dictionaries for argument maps
fxc_args_unpolar = {
  'LDA' : lda_fxc_args_unpolar,
  'GGA' : gga_fxc_args_unpolar,
  'MGGA': mgga_fxc_args_unpolar
}
fxc_args_polar = {
  'LDA' : lda_fxc_args_polar,
  'GGA' : gga_fxc_args_polar,
  'MGGA': mgga_fxc_args_polar
}
vxc_fxc_args_unpolar = {
  'LDA' : lda_vxc_fxc_args_unpolar,
  'GGA' : gga_vxc_fxc_args_unpolar,
  'MGGA': mgga_vxc_fxc_args_unpolar
}
vxc_fxc_args_polar = {
  'LDA' : lda_vxc_fxc_args_polar,
  'GGA' : gga_vxc_fxc_args_polar,
  'MGGA': mgga_vxc_fxc_args_polar
}

exc_args_unpolar = {
  'LDA' : lda_exc_args_unpolar,
  'GGA' : gga_exc_args_unpolar,
  'MGGA': mgga_exc_args_unpolar
}
exc_args_polar = {
  'LDA' : lda_exc_args_polar,
  'GGA' : gga_exc_args_polar,
  'MGGA': mgga_exc_args_polar
}
exc_vxc_args_unpolar = {
  'LDA' : lda_exc_vxc_args_unpolar,
  'GGA' : gga_exc_vxc_args_unpolar,
  'MGGA': mgga_exc_vxc_args_unpolar
}
exc_vxc_args_polar = {
  'LDA' : lda_exc_vxc_args_polar,
  'GGA' : gga_exc_vxc_args_polar,
  'MGGA': mgga_exc_vxc_args_polar
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
  print(name)
  xc_exc = XCFunc( meta_data.libxc_file, meta_data.xc_type, 'EXC' )
  # xc_vxc = XCFunc( meta_data.libxc_file, meta_data.xc_type, 'VXC' )
  xc_exc_vxc = XCFunc( meta_data.libxc_file, meta_data.xc_type, 'EXC_VXC' )
  xc_fxc = XCFunc( meta_data.libxc_file, meta_data.xc_type, 'FXC' )
  xc_vxc_fxc = XCFunc( meta_data.libxc_file, meta_data.xc_type, 'VXC_FXC' )

  xc_type = meta_data.xc_type

  is_lda  = xc_type == 'LDA'
  is_gga  = xc_type == 'GGA'
  is_mgga = xc_type == 'MGGA'
  needs_laplacian = (xc_type == 'MGGA') and meta_data.needs_laplacian
  is_kedf = meta_data.is_kedf
  is_epc  = meta_data.is_epc

  xc_struct_prefix = struct_prefix.format(
    meta_data.local_name, xc_type.lower(), str(is_lda).lower(), 
    str(is_gga).lower(), str(is_mgga).lower(), str(needs_laplacian).lower(), 
    str(is_kedf).lower(), str(is_epc).lower(), meta_data.dens_tol, 
    meta_data.dens_tol**(4.0/3.0) )

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

  exc_prefix_polar = func_prefix.format( 'exc', 'polar', 
    exc_args_polar[xc_type] )
  exc_vxc_prefix_polar = func_prefix.format( 'exc_vxc', 'polar', 
    exc_vxc_args_polar[xc_type] )
  
  # # Generate VXC-only function bodies
  # vxc_func_unpolar_body = "\n".join( [
  #   vxc_prefix_unpolar, 
  #   indent(xc_vxc.unpol_xc_body,2), 
  #   "\n}"
  # ] )
  # vxc_func_polar_body = "\n".join( [
  #   vxc_prefix_polar, 
  #   indent(xc_vxc.pol_xc_body,2), 
  #   "\n}"
  # ] )

  # vxc_func_unpolar_body = indent(vxc_func_unpolar_body, 2)
  # vxc_func_polar_body = indent(vxc_func_polar_body, 2)

  # Add functions for FXC and VXC_FXC
  fxc_prefix_unpolar = func_prefix.format('fxc', 'unpolar', 
    fxc_args_unpolar[xc_type])
  vxc_fxc_prefix_unpolar = func_prefix.format('vxc_fxc', 'unpolar', 
    vxc_fxc_args_unpolar[xc_type])

  fxc_prefix_polar = func_prefix.format('fxc', 'polar', 
    fxc_args_polar[xc_type])
  vxc_fxc_prefix_polar = func_prefix.format('vxc_fxc', 'polar', 
    vxc_fxc_args_polar[xc_type])

  # Add VXC-only function prefixes
  vxc_prefix_unpolar = func_prefix.format('vxc', 'unpolar', 
    vxc_args_unpolar[xc_type])
  vxc_prefix_polar = func_prefix.format('vxc', 'polar', 
    vxc_args_polar[xc_type])

  exc_func_unpolar_body = "\n".join( [
    exc_prefix_unpolar, 
    indent(xc_exc.unpol_xc_body,2), 
    "\n}"
  ] )
  exc_func_polar_body = "\n".join( [
    exc_prefix_polar, 
    indent(xc_exc.pol_xc_body,2), 
    "\n}"
  ] )

  exc_func_unpolar_body = indent(exc_func_unpolar_body, 2)
  exc_func_polar_body = indent(exc_func_polar_body, 2)

  exc_vxc_func_unpolar_body = "\n".join( [
    exc_vxc_prefix_unpolar, 
    indent(xc_exc_vxc.unpol_xc_body,2), 
    "\n}"
  ] )
  exc_vxc_func_polar_body = "\n".join( [
    exc_vxc_prefix_polar, 
    indent(xc_exc_vxc.pol_xc_body,2), 
    "\n}"
  ] )

  exc_vxc_func_unpolar_body = indent(exc_vxc_func_unpolar_body, 2)
  exc_vxc_func_polar_body = indent(exc_vxc_func_polar_body, 2)

  # Generate function bodies for FXC and VXC_FXC
  fxc_func_unpolar_body = "\n".join( [
    fxc_prefix_unpolar, 
    indent(xc_fxc.unpol_xc_body,2), 
    "\n}"
  ] )
  fxc_func_polar_body = "\n".join( [
    fxc_prefix_polar, 
    indent(xc_fxc.pol_xc_body,2), 
    "\n}"
  ] )

  fxc_func_unpolar_body = indent(fxc_func_unpolar_body, 2)
  fxc_func_polar_body = indent(fxc_func_polar_body, 2)

  vxc_fxc_func_unpolar_body = "\n".join( [
    vxc_fxc_prefix_unpolar, 
    indent(xc_vxc_fxc.unpol_xc_body,2), 
    "\n}"
  ] )
  vxc_fxc_func_polar_body = "\n".join( [
    vxc_fxc_prefix_polar, 
    indent(xc_vxc_fxc.pol_xc_body,2), 
    "\n}"
  ] )

  vxc_fxc_func_unpolar_body = indent(vxc_fxc_func_unpolar_body, 2)
  vxc_fxc_func_polar_body = indent(vxc_fxc_func_polar_body, 2)


  xc_struct_str = "\n".join([
    xc_struct_prefix, 
    xc_param_lines, 
    exc_func_unpolar_body,
    # vxc_func_unpolar_body,
    exc_vxc_func_unpolar_body,
    fxc_func_unpolar_body,
    vxc_fxc_func_unpolar_body,
    exc_func_polar_body,
    # vxc_func_polar_body,
    exc_vxc_func_polar_body,
    fxc_func_polar_body,
    vxc_fxc_func_polar_body,
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

