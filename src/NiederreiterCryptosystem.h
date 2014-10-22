#ifndef NIEDERREITER__H
#define NIEDERREITER__H

#include "GoppaCode.h"
#include "McElieceCryptosystem.h"

/*NiederreiterCryptosystem is setup according to http://en.wikipedia.org/wiki/Niederreiter_cryptosystem
/Courtois, Finiasz and Sendrier(CFS) showed how the Niederreiter cryptosystem can be used to derive a signature scheme.
/The choice of the code parameters is related to the probability that a random syndrome is decodable. 
/Courtois, Finiaz, and Sendrier suggest the parameter values n = 216 and t = 9.*/

class NiederreiterCryptosystem{
public:
	NiederreiterCryptosystem(int n, int m, GF2EX g, int delta);
	mat_GF2 Encrypt(mat_GF2 message);
	mat_GF2 Decrypt(mat_GF2 received_word);
	mat_GF2 Signature(mat_GF2 message, string mode);
	bool Verify(mat_GF2 signature, mat_GF2 message);

	mat_GF2 GetPublicKey(void);
	GoppaCode GetGoppaCode(void);	

private:
	GoppaCode _GoppaCode;
	GF2EX _g;
	int _m;
	int _n;
	int _delta;
	mat_GF2 _S; //Scrambler Matrix
	mat_GF2 _P; //Permutation Matrix
	mat_GF2 _PublicKey;
};

#endif
