#include "Overlaps.h"
#include <cmath>
#include <math.h>
#include "Parameters.h"
#include "gsl/gsl_sf_hermite.h"
#include "BasicFunctions.h"
#include "gsl/gsl_integration.h"

using namespace Parameters;


struct integration_params { int alk; int alb; int ik; int ib; };


double V_integration_function( double z, void* p )
{
	struct integration_params* params = ( struct integration_params* )p;
	int alk = (params->alk);
	int alb = (params->alb);
	int ik = (params->ik);
	int ib = (params->ib);

	double zk = (alk == 1) ? Overlaps::z1(z) : Overlaps::z2(z);
	double zb = (alb == 1) ? Overlaps::z1(z) : Overlaps::z2(z);
	double phiket = Overlaps::phiket(ik, zk);
	double phibra = Overlaps::phibra(ib, zb);
	double dv = Overlaps::DV(alk, z);
	return phiket * dv * phibra;
}


double S_integration_function(double z, void* p)
{
	struct integration_params* params = (struct integration_params*)p;
	int alk = (params->alk);
	int alb = (params->alb);
	int ik = (params->ik);
	int ib = (params->ib);

	double zk = (alk == 1) ? Overlaps::z1(z) : Overlaps::z2(z);
	double zb = (alb == 1) ? Overlaps::z1(z) : Overlaps::z2(z);
	double phiket = Overlaps::phiket(ik, zk);
	double phibra = Overlaps::phibra(ib, zb);
	return phiket * phibra;
}


double Overlaps::V(const double& z)
{
	using namespace std;
	return pow(omega, 2) * pow(pow(z, 2) - pow(distance / 2, 2), 2) / (2 * pow(distance, 2));
}


double Overlaps::z1(const double& z)
{
	return distance / 2 + z;
}


double Overlaps::z2(const double& z)
{
	return distance / 2 - z;
}


double Overlaps::Vho1(const double& z)
{
	using namespace std;
	return pow(omega, 2) * pow( z1(z), 2 ) / 2;
}


double Overlaps::Vho2(const double& z)
{
	using namespace std;
	return pow(omega, 2) * pow(z2(z), 2) / 2;
}


double Overlaps::DV1(const double& z)
{
	return V(z)-Vho1(z);
}


double Overlaps::DV2(const double& z)
{
	return V(z)-Vho2(z);
}


double Overlaps::DV(const int& alpha, const double& z)
{
	return (alpha == 1) ? DV1(z) : DV2(z);
}


double Overlaps::phiket(const int& nzk, const double& z_ket)
{
	using namespace std;
	return pow( omega / M_PI1, 1.0/4 ) * exp( -omega * pow( z_ket , 2 ) / 2) * 
		   gsl_sf_hermite( nzk, sqrt( omega ) * z_ket ) / 
		   sqrt( pow( 2, nzk ) * BasicFunctions::Factorial( nzk ) );
}


double Overlaps::phibra(const int& nzb, const double& z_bra)
{
	using namespace std;
	return pow(omega / M_PI1, 1.0 / 4) * exp(-omega * pow(z_bra, 2) / 2) *
		gsl_sf_hermite(nzb, sqrt(omega) * z_bra) /
		sqrt(pow(2, nzb) * BasicFunctions::Factorial(nzb));
}

int Overlaps::GetCount(const int& i, const int& j, const int& k, const int& z, const int& i_max, const int& j_max, const int& k_max, const int& z_max)
{
	return j_max * k_max * z_max * i + k_max * z_max * j + z_max * k + z;
}


double Overlaps::Eharmonic(const int& n_principal, const int& m_magnetic)
{
	return omega0 * (2.0 * n_principal + std::abs(m_magnetic) + 1.0);
}


double Overlaps::Eharmonic(const int& n_principal)
{
	double r = omega * (n_principal + 1.0 / 2);
	return r;
}


double Overlaps::VIntegrate(const int& alk, const int& alb, const int& ik, const int& ib)
{
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

	double result;
	double error;
	integration_params params = { alk, alb, ik, ib };

	double eps_abs = 1e-7;
	double eps_rel = 1e-7;

	gsl_function F;
	F.function = &V_integration_function;
	F.params = &params;

	gsl_integration_qagi(&F, eps_abs, eps_rel, 1000, w, &result, &error);

	gsl_integration_workspace_free(w);

	return result;
}


double Overlaps::SIntegrate(const int& alk, const int& alb, const int& ik, const int& ib)
{
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

	double result;
	double error;
	integration_params params = { alk, alb, ik, ib };

	double eps_abs = 1e-7;
	double eps_rel = 1e-7;

	gsl_function F;
	F.function = &S_integration_function;
	F.params = &params;

	gsl_integration_qagi(&F, eps_abs, eps_rel, 1000, w, &result, &error);

	gsl_integration_workspace_free(w);

	return result;
}

