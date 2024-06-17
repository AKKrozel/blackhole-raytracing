
#include "boyer_lindquist_metric.h"

boyer_lindquist_metric::boyer_lindquist_metric(double a0, double M0) {
	// Initialize the parameters a and M
	this->a = a0;
	this->M = M0;
}

void boyer_lindquist_metric::compute_metric(double r, double th) {
	// Compute all the required metric components in one go

	// useful constants
	double sin_th = std::sin(th);
	double cos_th = std::cos(th);
	double sinsq_th = sin_th * sin_th;
	double cossq_th = cos_th * cos_th;

	double rhosq = (r * r) + (this->a * this->a * cossq_th);
	double Delta = (r * r) + (this->a * this->a) - (2.0 * this->M * r);
	double Sigma = (r * r + this->a * this->a) * (r * r + this->a * this->a) - (this->a * this->a * Delta * sinsq_th);

	// compute non-derivative metrics
	this->alpha = std::sqrt(rhosq * Delta / Sigma);
	this->beta3 = -(2.0 * this->M * this->a * r / Sigma);

	this->g_00 = (2.0 * this->M * r / rhosq) - 1.0;
	this->g_03 = -(2.0 * this->M * this->a * r / rhosq) * sinsq_th;
	this->g_11 = rhosq / Delta;
	this->g_22 = rhosq;
	this->g_33 = (Sigma / rhosq) * sinsq_th;

	this->gamma11 = 1.0 / this->g_11;
	this->gamma22 = 1.0 / this->g_22;
	this->gamma33 = 1.0 / this->g_33;

	// compute r derivatives
	double d_rhosq_dr = 2.0 * r;
	double d_Delta_dr = 2.0 * r - 2.0 * M;
	double d_Sigma_dr = 2.0 * (r * r + this->a * this->a) * (2.0 * r) - (this->a * this->a) * d_Delta_dr * sinsq_th;

	this->d_alpha_dr = (1.0 / (2.0 * this->alpha)) * (Sigma * (rhosq * d_Delta_dr + Delta * d_rhosq_dr) - (rhosq * Delta * d_Sigma_dr)) / (Sigma * Sigma);
	this->d_beta3_dr = -(2.0 * this->M * this->a) * (Sigma - r * d_Sigma_dr) / (Sigma * Sigma);
	this->d_gamma11_dr = (rhosq * d_Delta_dr - Delta * d_rhosq_dr) / (rhosq * rhosq);
	this->d_gamma22_dr = -d_rhosq_dr / (rhosq * rhosq);
	this->d_gamma33_dr = (1.0 / sinsq_th) * (Sigma * d_rhosq_dr - rhosq * d_Sigma_dr) / (Sigma * Sigma);

	// compute theta derivatives (simplified because Delta is constant with respect to theta)
	double d_rhosq_dth = (this->a * this->a) * -2.0 * cos_th * sin_th;
	double d_Delta_dth = 0.0;
	double d_Sigma_dth = -(this->a * this->a) * Delta * (2.0 * sin_th * cos_th);

	this->d_alpha_dth = (1.0 / (2.0 * this->alpha)) * (Sigma * (d_rhosq_dth * Delta) - (rhosq * Delta * d_Sigma_dth)) / (Sigma * Sigma);
	this->d_beta3_dth = -(2.0 * this->M * this->a * r) * -d_Sigma_dth / (Sigma * Sigma);
	this->d_gamma11_dth = Delta * -d_rhosq_dth / (rhosq * rhosq);
	this->d_gamma22_dth = -d_rhosq_dth / (rhosq * rhosq);
	this->d_gamma33_dth = ((Sigma * sinsq_th * d_rhosq_dth) - rhosq * (Sigma * 2.0 * sin_th * cos_th + sinsq_th * d_Sigma_dth)) / (Sigma * Sigma * sinsq_th * sinsq_th);
	
}
