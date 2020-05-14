#pragma once


class Overlaps
{

public:
	//Fuctions for the potential
	static double V(const double& z);
	static double z1(const double& z);
	static double z2(const double& z);
	static double Vho1(const double& z);
	static double Vho2(const double& z);
	static double DV1(const double& z);
	static double DV2(const double& z);
	static double DV(const int& alpha, const double& z);

	//Wavefunctions
	static double phiket(const int& nzk, const double& z_ket);
	static double phibra(const int& nzb, const double& z_bra);

	//Counter acquisition
	static int GetCount(const int& i, const int& j, const int& k, const int& z, const int& i_max, const int& j_max, const int& k_max, const int& z_max);
	static int GetCount(const int& i, const int& j, const int& i_max, const int& j_max);

public:
	//Harmonic Oscilator Energies
	static double Eharmonic(const int& n_principal, const int& m_magnetic);
	static double Eharmonic(const int& n_principal);

public:
	//integration
	static double VIntegrate(const int& alk, const int& alb, const int& ik, const int& ib);
	static double SIntegrate(const int& alk, const int& alb, const int& ik, const int& ib);
};

