#include "Matrial.h"
#include "GreensFunction.h"

VecDoub EvaluationFunction(VecDoub &x) {
	// 1-d geometry info
	double k_ini = x[0];
	double phi_ini = x[1];
	int meshn = 300;
	double length = 30;
	vector<double> phi, phi_diff;
	for (int i = 0; i <= meshn; ++i) {
		phi.push_back(0);
		phi_diff.push_back(0);
	}
	// matrial cross section info
	SimpleMatrial core  (1.36, 0.0181, 0.0279);
	SimpleMatrial reflector (0.55, 0.0127,    0.0);

	double gf1, gf2;
	double b_core = sqrt((core.get_nusigf() / k_ini - core.get_siga()) / core.get_d());
	double b_refl = sqrt(reflector.get_siga() / reflector.get_d());
	// green's function
	double _x = length / meshn;
	gf1 = phi_star_core(_x, b_core);
	gf2 = phi_star_diff_core(_x, b_core);
	// in core region
	phi[0] = phi_ini, phi_diff[0] = 0.0;
	for (int i = 1; i <= meshn; ++i) {
		phi[i] = -2 * (gf1 * phi_diff[i - 1] + gf2 * phi[i - 1]);
		phi_diff[i] = (phi[i - 1] / 2 + gf2 * phi[i]) / gf1;
	}
	double phi_diff_core = phi_diff[meshn];
	// in the reflector region
	double phi_diff_reflector = -phi[meshn] * b_refl;

	VecDoub difference(2);
	difference[0] = abs(phi_diff_core * core.get_d() - phi_diff_reflector * reflector.get_d());
	difference[1] = abs(phi_ini-1.0);
	return difference;
}

















