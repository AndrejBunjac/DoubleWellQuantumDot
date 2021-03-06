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

	int nzmax_in = 0;
	int nrmax_in = 0;

	//Handle input parameters
	/*PROVERITI CASTING IZ DOUBLE U INT*/
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

	int nzmax = nzmax_in;
	int nrmax = nrmax_in;

	double* Vmat = new double[nzmax * nzmax * 2 * 2]();
	double* Smat = new double[nzmax * nzmax * 2 * 2]();

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
													double er1 = Overlaps::Eharmonic(irk, m1);
													double er2 = Overlaps::Eharmonic(jrk, m2);
													double ez1 = Overlaps::Eharmonic(izk);
													double ez2 = Overlaps::Eharmonic(jzk);

													double Hmat_el =	double(BasicFunctions::KroneckerDelta(irb, irk)) *
																		double(BasicFunctions::KroneckerDelta(jrb, jrk)) *
																	(	(er1 + er2 + ez1 + ez2) *
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

	/*
	TESTING PURPOSES
	{
		FILE* f = fopen("test.dat", "wb");
		FilesStrings::print_matrix(f, Hmat);
		fclose(f);
	}
	*/

	delete[]Vmat;
	delete[]Smat;

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