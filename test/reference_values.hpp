#pragma once
#include <exchcxx/enums/enums.hpp>
#include <array>
#include <vector>


struct lda_reference {
  int npts;
  std::vector<double> rho;
  std::vector<double> exc;
  std::vector<double> vrho;
};

struct gga_reference {
  int npts;
  std::vector<double> rho;
  std::vector<double> sigma;
  std::vector<double> exc;
  std::vector<double> vrho;
  std::vector<double> vsigma;
};

lda_reference load_lda_reference_values(ExchCXX::Kernel, ExchCXX::Spin);
gga_reference load_gga_reference_values(ExchCXX::Kernel, ExchCXX::Spin);

std::pair<int,std::vector<double>> load_reference_density(ExchCXX::Spin);
std::pair<int,std::vector<double>> load_reference_sigma(ExchCXX::Spin);

double load_reference_exx( ExchCXX::Kernel );
