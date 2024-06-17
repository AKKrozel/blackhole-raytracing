
#include "boyer_lindquist_metric.cpp"
//#include "integrator.h"

template <typename StopCondition>
rk45_dormand_prince integrate(double a0, double M0, const StopCondition& stop, double h, double t0,
	double tolerance_abs, double tolerance_rel, const std::vector<double>& y0) {

	boyer_lindquist_metric metric(a0, M0);

	auto dydx = [=, &metric](double x, const std::vector<double>& y) {
		metric.compute_metric(y[0], y[1]); // y[0] is r, y[1] is theta

		// compute the right hand sides of equations
		double ut = std::sqrt(metric.gamma11 * y[3] * y[3] + metric.gamma22 * y[4] * y[4] + metric.gamma33 * y[5] * y[5]) / metric.alpha;

		return std::vector<double>{
			metric.gamma11* (y[3] / ut),
				metric.gamma22* (y[4] / ut),
				metric.gamma33* (y[5] / ut) - metric.beta3,
				-metric.alpha * ut * metric.d_alpha_dr + y[5] * metric.d_beta3_dr - (1.0 / (2.0 * ut)) * (y[3] * y[3] * metric.d_gamma11_dr + y[4] * y[4] * metric.d_gamma22_dr + y[5] * y[5] * metric.d_gamma33_dr),
				-metric.alpha * ut * metric.d_alpha_dth + y[5] * metric.d_beta3_dth - (1.0 / (2.0 * ut)) * (y[3] * y[3] * metric.d_gamma11_dth + y[4] * y[4] * metric.d_gamma22_dth + y[5] * y[5] * metric.d_gamma33_dth),
				0.0
		};
	};

	// integrate with rk45
	rk45_dormand_prince rk45(6, tolerance_abs, tolerance_rel);
	rk45.integrate(dydx, stop, h, t0, y0);

	return rk45;

};

void print_cartesian_2D(rk45_dormand_prince rk45, std::string filename, double phi_stop) {

	// create csv file with position coordinates
	std::ofstream output_file_rk4(filename + ".csv");
	for (int i = 0; i < rk45.xs.size(); i++) {
		double r = rk45.result[i][0];
		double phi = rk45.result[i][2];

		if (phi <= phi_stop) {
			output_file_rk4 << r * std::cos(phi) << "," << r * std::sin(phi) << std::endl;
		}
	}
	output_file_rk4.close();

}

void print_cartesian_3D(rk45_dormand_prince rk45, std::string filename, double a) {
	
	// create csv file with position coordinates
	std::ofstream output_file_rk4(filename + ".csv");
	for (int i = 0; i < rk45.xs.size(); i++) {
		double r = rk45.result[i][0];
		double theta = rk45.result[i][1];
		double phi = rk45.result[i][2];

		output_file_rk4 << std::sqrt(r*r + a*a) * std::sin(theta) * std::cos(phi) << "," << std::sqrt(r*r + a*a) * std::sin(theta) * std::sin(phi) << "," << r * std::cos(theta) << std::endl;
	}
	output_file_rk4.close();

}

void print_radius_error(rk45_dormand_prince rk45, std::string filename, char orbit_type) {
	std::ofstream output_file_rk4(filename + ".csv");

	double r0 = 0.0;  // Initial radius
    
    // Determine r0 based on orbit type
    switch(orbit_type) {
        case 'A':
            r0 = 1 + std::sqrt(2);
            break;
        case 'B':
            r0 = 1 + std::sqrt(3);
            break;
        case 'C':
            r0 = 3;
            break;
        case 'D':
			r0 = 1 + 2 * std::sqrt(3);
			break;
		case 'E':
			r0 = 2;
			break;
	}

	for (int i = 0; i < rk45.xs.size(); i++) {
		double t = rk45.xs[i];
		double r = rk45.result[i][0];

		output_file_rk4 << t << "," << std::abs(r - r0)/r0 << std::endl;
	}
	output_file_rk4.close();
}

