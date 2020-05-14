// DoubleWellQuantumDot.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <cstdio>
#include <cmath>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>
#include <chrono>
#include "BasicFunctions.h"
#include "Overlaps.h"
#include "Parameters.h"
#include "FilesStrings.h"
#include <fstream>
#include <string>
#include <iostream>

using namespace Parameters;
using namespace std::chrono;

int main() {
	steady_clock::time_point startTime = steady_clock::now();

// INPUT HANDLING
/*
	int nzmax_in = 0;
	int nrmax_in = 0;

	//Handle input parameters
	std::string line;
	std::ifstream myfile("InputFile.txt");
	if (myfile.is_open()) {
		std::getline(myfile, line);
		nzmax_in = int(FilesStrings::extractDouble(line));
		std::cout << "nzmax =" << nzmax_in << '\n';
		std::getline(myfile, line);
		nrmax_in = int(FilesStrings::extractDouble(line));
		std::cout << "nrmax =" << nrmax_in << '\n';
		std::getline(myfile, line);
	}
	else
	{
		std::cout << "Input is missing!" << std::endl;
		std::cout << "Check code documentation about input file format." << std::endl;
		std::cout << "Terminating!" << std::endl;
		std::cin.get();
		return 1;
	}
*/
	double* Vmat = new double[nzmax * nzmax * 2 * 2]();
	double* Smat = new double[nzmax * nzmax * 2 * 2]();
	double* Ehoz = new double[nzmax]();
	double* Ehorm = new double[nrmax * mmax]();

	for (int ik = 0; ik < nzmax; ik++)
	{
		Ehoz[ik] = Overlaps::Eharmonic(ik);
	}

	for (int ik = 0; ik < nrmax; ik++) 
	{
		for (int mk = 0; mk < mmax; mk++)
		{
			int count = Overlaps::GetCount(ik,mk,nrmax,mmax);
			Ehorm[count] = Overlaps::Eharmonic(ik,mk);
		}
	}
	
	for (int ib = 0; ib < nzmax; ib++)
	{

		for (int ik = 0; ik < nzmax ; ik++)
		{

			for (int alb = 0; alb < 2; alb++)
			{

				for (int alk = 0; alk < 2; alk++)
				{
					int alk_in = alk + 1;
					int alb_in = alb + 1;
					int count = Overlaps::GetCount(ib, ik, alb, alk, nzmax, nzmax, 2, 2);
					Vmat[count] = Overlaps::VIntegrate(alk_in, alb_in, ik, ib);
					Smat[count] = Overlaps::SIntegrate(alk_in, alb_in, ik, ib);
				}
			}
		}
	}

	const int matsize = 4 * nzmax * nzmax * nrmax * nrmax;

	gsl_matrix* Hmat = gsl_matrix_calloc(matsize, matsize);
	gsl_matrix* Emat = gsl_matrix_calloc(matsize, matsize);

	for (int beb = 0; beb < 2; beb ++)
	{
		for (int bek = 0; bek < 2; bek++)
		{
			for (int alb = 0; alb < 2; alb++)
			{
				for (int alk = 0; alk < 2; alk++)
				{
					for (int jzb = 0; jzb < nzmax; jzb++)
					{
						for (int jzk = 0; jzk < nzmax; jzk++)
						{
							for (int izb = 0; izb < nzmax; izb++)
							{
								for (int izk = 0; izk < nzmax; izk++)
								{
									for (int jrb = 0; jrb < nrmax; jrb++)
									{
										for (int jrk = 0; jrk < nrmax; jrk++)
										{
											for (int irb = 0; irb < nrmax; irb++)
											{
												for (int irk = 0; irk < nrmax; irk++)
												{
													//counters for 2D matrices Hmat and Emat
													int i = irb * 4 * nrmax * nzmax * nzmax +
															jrb * 4 * nzmax * nzmax +
															2 * alb * nzmax * nzmax +
															beb * nzmax * nzmax +
															izb * nzmax +
															jzb;

													int j = irk * 4 * nrmax * nzmax * nzmax +
															jrk * 4 * nzmax * nzmax +
															2 * alk * nzmax * nzmax +
															bek * nzmax * nzmax +
															izk * nzmax +
															jzk;

													int count_a = Overlaps::GetCount(izb, izk, alb, alk, nzmax, nzmax, 2, 2);
													int count_b = Overlaps::GetCount(jzb, jzk, beb, bek, nzmax, nzmax, 2, 2);
													int count_e1 = Overlaps::GetCount(irk,m1,nrmax,mmax);
													int count_e2 = Overlaps::GetCount(jrk,m2,nrmax,mmax);

													double Hmat_el =	double(BasicFunctions::KroneckerDelta(irb, irk)) *
																		double(BasicFunctions::KroneckerDelta(jrb, jrk)) *
																	(	(Ehorm[count_e1] + Ehorm[count_e2] + Ehoz[izk] + Ehoz[jzk]) *
																		Smat[count_a] * Smat[count_b] +
																		Vmat[count_a] * Smat[count_b] +
																		Smat[count_a] * Vmat[count_b]);

													gsl_matrix_set(Hmat, i, j, Hmat_el);

													double Emat_el =	double(BasicFunctions::KroneckerDelta(irb, irk)) *
																		double(BasicFunctions::KroneckerDelta(jrb, jrk)) *
																		Smat[count_a] * Smat[count_b];

													gsl_matrix_set(Emat, i, j, Emat_el);

												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	delete[]Vmat;
	delete[]Smat;
	delete[]Ehorm;
	delete[]Ehoz;

	gsl_vector* eigenvalues =	gsl_vector_alloc(matsize);
	gsl_eigen_gensymm_workspace* w = gsl_eigen_gensymm_alloc(matsize);

	gsl_eigen_gensymm(Hmat, Emat, eigenvalues, w);
	gsl_sort_vector(eigenvalues);

	std::ofstream outputfile("EigenValuesSorted.txt");
	if (!outputfile.is_open())
	{
		std::cout << "Output file could not be opened! Terminating!" << std::endl;
		return 1;
	}

	for (int i = 0; i < matsize; i++)
	{
		double eigenvalue = gsl_vector_get(eigenvalues, i);
		outputfile << std::setprecision(10) << std::fixed;
		outputfile << eigenvalue << '\n';
		printf("eigenvalue =%g\n", eigenvalue);
	}

	gsl_eigen_gensymm_free(w);
	gsl_vector_free(eigenvalues);
	gsl_matrix_free(Hmat);
	gsl_matrix_free(Emat);

	steady_clock::time_point endTime = steady_clock::now();
	duration<float> elapsed = endTime - startTime;
	printf("Time elapsed since simulation start: %g\n", elapsed.count());

}