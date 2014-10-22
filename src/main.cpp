#include "GoppaCode.h"
#include "McElieceCryptosystem.h"
#include "NiederreiterCryptosystem.h"
#include <time.h>

mat_GF2 GetRandomMessage(int message_length);
mat_GF2 GetRandomMessageWithWeight(int message_length, int message_weight);

const string OPTION = "McElieceCryptosystem";
int main()
{
	/*initiation*/
	int m=5; int t=5; int delta = 2;
        int n= int(pow(double(2),double(m)));
	GF2X gf2e_mod;
	BuildIrred(gf2e_mod,m); // generate an irreducible polynomial over GF(2^m)
	cout << "GF_2m field_poly: " << gf2e_mod << endl;    	
	GF2E::init(gf2e_mod);// construct GF(2^m)   
	GF2EX g;	
	BuildIrred(g,t);  //generate goppa polynomial with t degree  
	//cout << "goppa_poly: " << deg(g) << endl;
	GF2EX g2 = g*g;
	if (OPTION == "GoppaCode"){
	     	GoppaCode goppacode(n,m,g2);
	
		GF2EX g = goppacode.GetGoppaPolynomial();
		cout << "irr_poly: " << g << endl;
		/*Timer starts*/	
		long beginTime = clock();
		mat_GF2 message = GetRandomMessage(n-m*t);
		cout << "random message: " << message << endl;
		
		mat_GF2 codeword = goppacode.Encode(message);
		cout << "codeword: " << codeword << endl;
		
		mat_GF2 error = GetRandomMessageWithWeight(n,t);
		cout << "error: " << error << endl;
		
		mat_GF2 ciphertext = codeword + error;
		cout << "ciphertext: " << ciphertext << endl;
		
		/*Timer ends*/        
		long endTime = clock();
		cout << "Encryption Time Elpased: " << double(endTime-beginTime)/CLOCKS_PER_SEC << endl;        
		//Decrypt the ciphertext using algebraic syndrome decoding algorihtm
		/*Timer starts*/	
	 	beginTime = clock();        
		mat_GF2 recoveredtext1 = goppacode.Decode(ciphertext, "Euclidean");
		cout << "recovered message1: " << recoveredtext1 << endl;
		/*Timer ends*/        
		endTime = clock();
		cout << "Decryption Time Elpased(seconds): " << double(endTime-beginTime)/CLOCKS_PER_SEC << endl; 
		if (recoveredtext1 == codeword)
		        cout << "It works!" << endl;
		else
		        cout << "Something wrong!" << endl;
		
		
		//calculate sigma
		vector<GF2E> codelocators = goppacode.GetCodeLocators();
		GF2EX sigma = conv<GF2EX>(1);		
		GF2EX X; SetX(X);		
		for (int i = 0; i < error.NumCols(); i++)
			if (error[0][i] == 1)
				sigma *= X-codelocators[i];
		//calculate omega
		GF2EX omega = conv<GF2EX>(0);		
		for (int i = 0; i < error.NumCols(); i++)
			if (error[0][i] == 1)
				omega += sigma/(X-codelocators[i]);
	
		//Compute the syndrome polynomial.
		GF2EX syndrome_poly = conv<GF2EX>(0);
		for (int i = 0; i < error.NumCols(); i++)
			syndrome_poly += goppacode.GetSyndromeCalculator()[i]*error[0][i];	
		GF2EX tmp = sigma*syndrome_poly-omega;	
		GF2EX mod; rem(mod, tmp, goppacode.GetGoppaPolynomial());
		if (mod == 0)
			cout << "Key equation1 holds." << endl;
		else 
			cout << "Key equation1 does not hold." << endl;
		
		
				

	}

	else if (OPTION == "McElieceCryptosystem"){
		McElieceCryptosystem McEliece(n,m,g);
		mat_GF2 message = GetRandomMessage(n-m*t);
		cout << "random message: " << message << endl;
		
		mat_GF2 ciphertext = McEliece.Encrypt(message);
		cout << "ciphertext: " << ciphertext << endl;

		mat_GF2 recoveredtext1 = McEliece.Decrypt(ciphertext);
		cout << "recovered message1: " << recoveredtext1 << endl;
		if (recoveredtext1 == message)
		        cout << "It works!" << endl;
		else
		        cout << "Something wrong!" << endl;
	}

	else if (OPTION == "NiederreiterCryptosystem"){
		NiederreiterCryptosystem Niederreiter(n,m,g,delta); 
		mat_GF2 message = GetRandomMessageWithWeight(n,deg(g));
		cout << "random message: " << message << endl;

		mat_GF2 ciphertext = Niederreiter.Encrypt(message);
		cout << "ciphertext: " << ciphertext << endl; 

		mat_GF2 recoveredtext1 = Niederreiter.Decrypt(ciphertext);
		cout << "recovered message1: " << recoveredtext1 << endl;
		if (recoveredtext1 == message)
		        cout << "It works!" << endl;
		else
		        cout << "Something wrong!" << endl;

		mat_GF2 message_hash = GetRandomMessageWithWeight(m*t,rand()%(m*t));//length:m*t
		//mat_GF2 message_hash = ciphertext;
		//cout << "hashed message: " << message_hash << endl;
		/*Timer starts*/	
		long beginTime = clock();
		mat_GF2 signature = Niederreiter.Signature(message_hash,"random");
		//cout << "signature: " << signature << endl;

		
		bool verification = Niederreiter.Verify(signature, message_hash);
		//cout << "Siganture: " << verification << endl;
		/*Timer ends*/        
		long endTime = clock();
		cout << "Decryption Time Elpased(seconds): " << double(endTime-beginTime)/CLOCKS_PER_SEC << endl;

	
	}

	return 0;        
}




 /* Assistant Functions */

mat_GF2 GetRandomMessage(int message_length)
{
	mat_GF2 message;
	message.SetDims(1,message_length);
	for (int i = 0; i < message_length; i++){
		message[0][i] = random_GF2();
	}
	return message;
}

mat_GF2 GetRandomMessageWithWeight(int message_length, int message_weight)
{
	if(message_length < message_weight){
		cout <<	"mat_GF2 GetRandomMessageWithWeight(int message_length, int message_weight): Error! To much message_weight." << endl;
		exit(-1);	
	}
	mat_GF2 message;
	message.SetDims(1,message_length);
	vector<int> rng(message_length); 
	for(int i = 0; i < message_length; i++)
		rng[i] = i;
	for (int i = 0; i < message_weight; i++){
		int p = rng.size()*(rand()/double(RAND_MAX));
		message[0][rng[p]] = 1;
		for(int j = p+1; j < rng.size(); j++)
			rng[j-1] = rng[j];
		rng.resize(--message_length);
		
	}
	return message;
}




