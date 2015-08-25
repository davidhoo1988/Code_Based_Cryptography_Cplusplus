#ifndef MCELIECE__H
#define MCELIECE__H

#include "GoppaCode.h"
/*McElieceCryptosystem is setup according to http://en.wikipedia.org/wiki/McEliece_cryptosystem
**McEliece originally suggested security parameter sizes of n=1024, k=524, t=50, 
**resulting in a public key size of 524*(1024-524) = 262,000 bits*/

class McElieceCryptosystem{
public:
	McElieceCryptosystem(int n, int m, GF2EX g);
	mat_GF2 Encrypt(mat_GF2 message);
	mat_GF2 Decrypt(mat_GF2 message);
	mat_GF2 ScramblerMatrixGen(int k, string mode = "random");
	Mat<ZZ> LargeScramblerMatrixGen(int k, string mode = "random");
	mat_GF2 PermutationMatrixGen(int n, string mode = "random");
	Mat<ZZ> LargePermutationMatrixGen(int n, string mode = "random");

	mat_GF2 GetPublicKey(void);
	GoppaCode GetGoppaCode(void);
	mat_GF2 GetPermutationMatrix(void);
	mat_GF2 GetScramblerMatrix(void);
	Mat<ZZ> GetLargePermutationMatrix(void);
	Mat<ZZ> GetLargeScramblerMatrix(void);
private:
	GoppaCode _GoppaCode;
	GF2EX _g;
	int _m;
	int _n;
	mat_GF2 _P;//Permutation Matrix
	mat_GF2 _S;//Scrambler Matrix
	Mat<ZZ> _LargeP;
	Mat<ZZ> _LargeS;
	mat_GF2 _PublicKey;
};



#endif
