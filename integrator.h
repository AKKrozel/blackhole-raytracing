#pragma once

#include <vector>
#include <string>

class rk45_dormand_prince;

template <typename StopCondition>
rk45_dormand_prince integrate(double a0, double M0, const StopCondition& stop, double h, double t0,
	double tolerance_abs, double tolerance_rel, const std::vector<double>& y0);

void print_cartesian_2D(rk45_dormand_prince rk45, std::string filename, double phi_stop);

void print_cartesian_3D(rk45_dormand_prince rk45, std::string filename);

void print_radius_error(rk45_dormand_prince rk45, std::string filename, char orbit_type);

