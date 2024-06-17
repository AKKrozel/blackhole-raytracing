
#include "integrator.cpp"
#include <cfloat>
#include <iomanip>


using namespace std;

int main() {

	//Constants
	double M = 1.0; //Mass of Blackhole
	double a = 0.99; //Rotation Speed of Blackhole
	double I_em = 1.0; //Intensity of light emmited by accretion disk
	double aSQ = a * a; //a^2 is a frequently used constant
	double r_in = 5.0 * M; //Inner radius of ring
	double r_out = 20.0 * M; //Outter radius of ring
	double D = 1400.0; //Observer distance
	double theta_0 = 85.0; //Observer angle in degrees
	theta_0 = theta_0 * M_PI / 180.0; //Observer angle in radians
	int resolutionX = 150; //Number of pixels in x direction
	int resolutionY = 150; //Number of pixels in y direction
	double betaMax = 0.9; //Max angle to screen to side of observer positionin degrees
	double alphaMax = betaMax; /*9.0/16.0*/ //Max angle to screen above/below observer position in degrees
	betaMax *= M_PI / 180.0;
	alphaMax *= M_PI / 180.0;
	double distance_scalingY = 1.0 / resolutionY;
	double distance_scalingX = 1.0 / resolutionX;
	distance_scalingY += distance_scalingY * distance_scalingY; // correction to make the beginning and end pixels both be at the maximum angle
	distance_scalingX += distance_scalingX * distance_scalingX;
	double alphaIncrement = 2.0 * alphaMax * distance_scalingY;
	double betaIncrement = 2.0 * betaMax * distance_scalingX;

	//iterate through every pixel
	double r_stop_outter = (1.0 + 1e-5) * D * std::sqrt(1.0 + alphaMax * alphaMax + betaMax * betaMax); //A photon should never be able to get farther than the screen and still hit the screen if aimed at the black hole
	double r_H = M + std::sqrt(M * M - aSQ);
	double r_stop_inner = (1.0 + 1e-5) * r_H;
	auto stop = [&](double x, const std::vector<double>& y) {
		double radius = y[0];
		double theta = y[1];
		double collisionTolorance = 0.2; //degrees (angle around disk where detection is allowed)
		collisionTolorance = collisionTolorance * M_PI / 180.0; //radians
		return ((radius > r_in) && (radius < r_out) && (std::abs(theta - M_PI / 2.0) < collisionTolorance)) || (radius > r_stop_outter) || (radius < r_stop_inner); };
	
	
	std::ofstream output_file("pixelIntensities.csv");
	std::ofstream output_file_t("times.csv");


	double alpha = -1.0 * alphaMax;

	for (int i = 0; i < resolutionY; i++) {

		double beta = -1.0 * betaMax;

		for (int j = 0; j < resolutionX; j++) {

			//x,y position on screen
			double x_sc = D * beta;
			double y_sc = D * alpha;

			//initial generalized positions
			double r = std::sqrt(D * D + x_sc * x_sc + y_sc * y_sc);
			double theta = theta_0 - alpha;
			double phi = beta;

			//intermediary values
			double sin_theta = std::sin(theta);
			double cos_theta = std::cos(theta);
			double sinSQtheta = sin_theta * sin_theta;
			double cosSQtheta = cos_theta * cos_theta;
			double rSQ = r * r;
			double rhoSQ = rSQ + aSQ * cosSQtheta;
			double delta = rSQ + aSQ - 2.0 * M * r;
			double sigma = (rSQ + aSQ) * (rSQ + aSQ) - aSQ * delta * sinSQtheta;
			double g_thetatheta = rhoSQ;
			double g_rr = rhoSQ / delta;
			double g_phiphi = sigma * sinSQtheta / rhoSQ;

			//initial generalized momenta
			double u_r = -1.0 * std::sqrt(g_rr) * std::cos(beta) * std::cos(alpha);
			double u_theta = std::sqrt(g_thetatheta) * std::sin(alpha);
			double u_phi = std::sqrt(g_phiphi) * std::sin(beta) * std::cos(alpha);

			//raytrace by running integration alg
			vector<double> initialConditons{ r, theta, phi, u_r, u_theta, u_phi }; 
			double t_0 = 0.0;
			double h = 1e-30;
			double tolorance_abs = 1e-15;
			double tolorance_rel = 1e-15;

			rk45_dormand_prince rk45dp = integrate(a, M, stop, h, t_0, tolorance_abs, tolorance_rel, initialConditons);

			//final conditions
			vector<double> finalConditions = rk45dp.result[rk45dp.result.size() - 1];
			double r_f = finalConditions[0];

			//light intensity
			double I_obs = 0.0;

			if (r_f > r_stop_inner && r_f < r_stop_outter) {

				//final conditios
				double theta_f = finalConditions[1];
				double phi_f = finalConditions[2];
				double u_r_f = -1.0 * finalConditions[3]; //Need to reverse the generalized momenta directions since they were found by tracing backwards
				double u_theta_f = -1.0 * finalConditions[4];
				double u_phi_f = -1.0 * finalConditions[5];
				

				//intermediary values
				double r_fSQ = r_f * r_f;
				double sin_theta_f = std::sin(theta_f);
				double cos_theta_f = std::cos(theta_f);
				double sinSQtheta_f = sin_theta_f * sin_theta_f;
				double cosSQtheta_f = cos_theta_f * cos_theta_f;
				
				double rho_fSQ = r_fSQ + aSQ * cosSQtheta_f;
				double delta_f = r_fSQ + aSQ - 2.0 * M * r_f;
				double sigma_f = (r_fSQ + aSQ) * (r_fSQ + aSQ) - aSQ * delta_f * sinSQtheta_f;

				double g_tt_f = (2.0 * M * r_f / rho_fSQ) - 1.0;
				double g_tphi_f = -2.0 * M * a * r_f * sinSQtheta_f / rho_fSQ;
				double g_rr_f = rho_fSQ / delta_f;
				double g_thetatheta_f = rho_fSQ;
				double g_phiphi_f = sigma_f * sinSQtheta_f / rho_fSQ;

				double gamma11_f = 1.0 / g_rr_f;
				double gamma22_f = 1.0 / g_thetatheta_f;
				double gamma33_f = 1.0 / g_phiphi_f;

				double beta1 = 0.0;
				double beta2 = 0.0;
				double beta3 = -2.0 * M * a * r_f / sigma_f;

				double metricAlpha_f = std::sqrt(rho_fSQ * delta_f / sigma_f);
				double ut_f = std::sqrt(gamma11_f * u_r_f * u_r_f + gamma22_f * u_theta_f * u_theta_f + gamma33_f * u_phi_f * u_phi_f) / metricAlpha_f;
				double u_t_f = -metricAlpha_f * metricAlpha_f * ut_f + (u_r_f * beta1 + u_theta_f * beta2 + u_phi_f * beta3);

				double omega = 1.0 / (a + (std::pow(r_f, 1.5) / std::sqrt(M)));

				double onePlusZ = (1.0 + omega * u_phi_f / u_t_f) / std::sqrt(-1.0 * g_tt_f - omega * omega * g_phiphi_f - 2.0 * omega * g_tphi_f);
				
				//light intensity
				I_obs = I_em / std::pow(onePlusZ, 3.0);

			}

			//Output light path
			std::string filename = to_string(i) + "," +  to_string(j);
			print_cartesian_3D(rk45dp, filename, a);

			//Output time to screen
			std::vector<double> xs = rk45dp.get_xs();
			long double t = xs[xs.size() - 1];

			output_file << I_obs << ",";
			output_file_t << std::setprecision(60) << t << ",";

			beta += betaIncrement;
		}

		output_file << endl;
		output_file_t << endl;

		alpha += alphaIncrement;
	}

	output_file.close();
	output_file_t.close();

	
	return 0;

}
