#ifndef GOPPACODE__H
#define GOPPACODE__H

#include <iostream>
#include <cmath>
#include <vector>
#include <stack>
#include <fstream>

#include <NTL/ZZ.h>
#include <NTL/GF2X.h>
#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/GF2EXFactoring.h>
#include <NTL/ZZXFactoring.h>

#include <NTL/mat_GF2.h>

using namespace std;
using namespace NTL;

class GoppaCode{
public:
	GoppaCode(int n, int m, GF2EX g);
	GoppaCode();

	mat_GF2 Encode(mat_GF2 message);
	mat_GF2 Decode(mat_GF2 word, string mode="Patterson");
	mat_GF2 SyndromeDecode(GF2EX syndrome_poly, string mode="Patterson");
	GF2EX GetGoppaPolynomial(void);
	mat_GF2 GetParityCheckMatrix(void);
	mat_GF2 GetGeneratorMatrix(void);
	vector<GF2E> GetCodeLocators(void);
	vector<GF2EX> GetSyndromeCalculator(void);
	
	/*Assistant*/
 	GF2E _GetGf2eGenerator(int degree);
	void _split(GF2EX& p0, GF2EX& p1, const GF2EX p);
	ZZ _norm(GF2EX a, GF2EX b);
	void _lattice_basis_reduce(GF2EX& M, GF2EX& N, const GF2EX s);
	void _extended_euclidean(GF2EX& sigma, GF2EX& omega, GF2EX A, GF2EX B, const int degree);
	void _extended_euclidean_I(GF2EX& sigma, GF2EX& omega, GF2EX A, GF2EX B, const int degree);
	void _extended_euclidean_II(GF2EX& sigma, GF2EX& omega, GF2EX A, GF2EX B, const int degree);
private:
	int _n;
	int _m;
	GF2EX _g;	
	mat_GF2 _H_Goppa;
	mat_GF2 _G_Goppa;
	vector<GF2EX> _SyndromeCalculator;
	vector<GF2E> _codelocators;
};


void combination(vector<vector<int> >& out_comb_list, vector<int> in_comb_list, int n, int k);
void integer_factor(vector<int>& out_factor_list, int in_small_integer);

int GetRowVectorWeight(vector<GF2> input_vec);
int GetRowVectorWeight(mat_GF2 input_marix);
//void PrintResult(mat_GF2 data,string addr);
#endif
