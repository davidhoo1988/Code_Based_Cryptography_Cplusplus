#include "GoppaCode.h"
#include "McElieceCryptosystem.h"
#include "NiederreiterCryptosystem.h"
#include <time.h>
#include <fstream>

mat_GF2 GetRandomMessage(int message_length);
mat_GF2 GetRandomMessageWithWeight(int message_length, int message_weight);

const string OPTION = "NiederreiterCryptosystem";
int main()
{
	/*initiation*/
	int m = 16; int t = 9; int delta = 3;
    int n= int(pow(double(2),double(m)));
	GF2X gf2e_mod;

	//BuildIrred(gf2e_mod,m); // randomly generate an irreducible polynomial over GF(2^m)
	// fix the field polynomial to be f = x**16 + x**5 + x**3 + x**2 + 1
	SetCoeff(gf2e_mod,16,1);
	SetCoeff(gf2e_mod,5,1);
	SetCoeff(gf2e_mod,3,1);
	SetCoeff(gf2e_mod,2,1);
	SetCoeff(gf2e_mod,0,1);

	//cout << "GF_2m field_poly: " << gf2e_mod << endl;    	
	GF2E::init(gf2e_mod);// construct GF(2^m)   
	GF2EX g;	
	BuildIrred(g,t);  //generate goppa polynomial with t degree 
	GF2X g0, g1, g2, g3, g4, g5, g6, g7, g8, g9;
	
	
	//g0=[0 1 1 0 1 0 1 0 1 1 0 1 1 0 1 1]
	SetCoeff(g0,1,1);
	SetCoeff(g0,2,1);
	SetCoeff(g0,4,1);
	SetCoeff(g0,6,1);
	SetCoeff(g0,8,1);
	SetCoeff(g0,9,1);
	SetCoeff(g0,11,1);
	SetCoeff(g0,12,1);
	SetCoeff(g0,14,1);
	SetCoeff(g0,15,1);
	//g1=[0 1 1 1 0 0 1 1 1 1 0 0 1 0 1]
	SetCoeff(g1,1,1);
	SetCoeff(g1,2,1);
	SetCoeff(g1,3,1);
	SetCoeff(g1,6,1);
	SetCoeff(g1,7,1);
	SetCoeff(g1,8,1);
	SetCoeff(g1,9,1);
	SetCoeff(g1,12,1);
	SetCoeff(g1,14,1);
	//g2=[1 1 0 0 1 1 1 0 1 0 1 1 1 1 0 1]
	SetCoeff(g2,0,1);
	SetCoeff(g2,1,1);
	SetCoeff(g2,4,1);
	SetCoeff(g2,5,1);
	SetCoeff(g2,6,1);
	SetCoeff(g2,8,1);
	SetCoeff(g2,10,1);
	SetCoeff(g2,11,1);
	SetCoeff(g2,12,1);
	SetCoeff(g2,13,1);
	SetCoeff(g2,15,1);
	//g3=[0 0 0 1 1 1 1 1 1 1 1 1 1 1 1]
	SetCoeff(g3,3,1);
	SetCoeff(g3,4,1);
	SetCoeff(g3,5,1);
	SetCoeff(g3,6,1);
	SetCoeff(g3,7,1);
	SetCoeff(g3,8,1);
	SetCoeff(g3,9,1);
	SetCoeff(g3,10,1);
	SetCoeff(g3,11,1);
	SetCoeff(g3,12,1);
	SetCoeff(g3,13,1);
	SetCoeff(g3,14,1);	
	//g4=[1 1 1 1 1 0 1 0 0 1 0 1 0 0 1 1]
	SetCoeff(g4,0,1);
	SetCoeff(g4,1,1);
	SetCoeff(g4,2,1);
	SetCoeff(g4,3,1);
	SetCoeff(g4,4,1);
	SetCoeff(g4,6,1);
	SetCoeff(g4,9,1);
	SetCoeff(g4,11,1);
	SetCoeff(g4,14,1);
	SetCoeff(g4,15,1);
	//g5=[1 1 0 0 1 0 1 0 1 1 1 0 1 1]
	SetCoeff(g5,0,1);
	SetCoeff(g5,1,1);
	SetCoeff(g5,4,1);
	SetCoeff(g5,6,1);
	SetCoeff(g5,8,1);
	SetCoeff(g5,9,1);
	SetCoeff(g5,10,1);
	SetCoeff(g5,12,1);
	SetCoeff(g5,13,1);
	//g6=[1 0 1 1 0 0 0 0 1 0 1 1 1 1 0 1]
	SetCoeff(g6,0,1);
	SetCoeff(g6,2,1);
	SetCoeff(g6,3,1);
	SetCoeff(g6,8,1);
	SetCoeff(g6,10,1);
	SetCoeff(g6,11,1);
	SetCoeff(g6,12,1);
	SetCoeff(g6,13,1);
	SetCoeff(g6,15,1);
	//g7=[0 1 0 1 1 1 0 1 1 0 1 0 1 0 1 1]
	SetCoeff(g7,1,1);
	SetCoeff(g7,3,1);
	SetCoeff(g7,4,1);
	SetCoeff(g7,5,1);
	SetCoeff(g7,7,1);
	SetCoeff(g7,8,1);
	SetCoeff(g7,10,1);
	SetCoeff(g7,12,1);
	SetCoeff(g7,14,1);
	SetCoeff(g7,15,1);
	//g8=[0 1 0 0 1 1 0 0 0 0 1 1 1 0 1 1] 
	SetCoeff(g8,1,1);
	SetCoeff(g8,4,1);
	SetCoeff(g8,5,1);
	SetCoeff(g8,10,1);
	SetCoeff(g8,11,1);
	SetCoeff(g8,12,1);
	SetCoeff(g8,14,1);
	SetCoeff(g8,15,1);
	//g9=[1]
	SetCoeff(g9,0,1);

	SetCoeff(g,0,conv<GF2E>(g0));
	SetCoeff(g,1,conv<GF2E>(g1));
	SetCoeff(g,2,conv<GF2E>(g2));
	SetCoeff(g,3,conv<GF2E>(g3));
	SetCoeff(g,4,conv<GF2E>(g4));
	SetCoeff(g,5,conv<GF2E>(g5));
	SetCoeff(g,6,conv<GF2E>(g6));
	SetCoeff(g,7,conv<GF2E>(g7));
	SetCoeff(g,8,conv<GF2E>(g8));
	SetCoeff(g,9,conv<GF2E>(g9));
	
	cout << "goppa_poly: " << g << endl; 

	GF2EX g_square = g;
	if (OPTION == "GoppaCode"){
	    GoppaCode goppacode(n,m,g_square);
		GF2EX g = goppacode.GetGoppaPolynomial();
		
		/*Timer starts*/	
		long beginTime = clock();
		mat_GF2 message = GetRandomMessage(n-m*t);
		//cout << "random message: " << message << endl;
		
		mat_GF2 codeword = goppacode.Encode(message);
		//cout << "codeword: " << codeword << endl;
		
		mat_GF2 error = GetRandomMessageWithWeight(n,t);
		//cout << "error: " << error << endl;
		
		mat_GF2 ciphertext = codeword + error;
		//cout << "ciphertext: " << ciphertext << endl;
		
		/*Timer ends*/        
		long endTime = clock();
		cout << "Encryption Time Elpased: " << double(endTime-beginTime)/CLOCKS_PER_SEC << endl;        
		//Decrypt the ciphertext using algebraic syndrome decoding algorihtm
		/*Timer starts*/	
	 	beginTime = clock();        
		mat_GF2 recoveredtext1 = goppacode.Decode(ciphertext, "Euclidean");
		//cout << "recovered message1: " << recoveredtext1 << endl;
		/*Timer ends*/        
		endTime = clock();
		cout << "Decryption Time Elpased(seconds): " << double(endTime-beginTime)/CLOCKS_PER_SEC << endl; 
		if (recoveredtext1 == codeword)
		        cout << "It works!" << endl;
		else
		        cout << "Something wrong!" << endl;  
		
		//output the parity check matrix to the file
		const char *path="/home/david/Dropbox/GoppaCode/Code_Based_Cryptography_Cplusplus/src/Data/ParityCheckMatrix.txt";
		ofstream ParityCheckMatrix;
		ParityCheckMatrix.open(path);		
		if (ParityCheckMatrix.is_open()){
			ParityCheckMatrix << "Row: " << goppacode.GetGeneratorMatrix().NumRows()
			<< "\t" << "Col: " << goppacode.GetGeneratorMatrix().NumCols() << endl;
			/*for(int i = 0; i < t; i++){
				for(int j = 0; j < n; j++){
					ParityCheckMatrix << goppacode.GetParityCheckMatrixPoly()[i][j]<<"\t";
				}
				ParityCheckMatrix << endl;
			}*/
			ParityCheckMatrix << goppacode.GetParityCheckMatrixPoly();
			ParityCheckMatrix.close();
		} 
		else {
			cout << "file not saved." << endl;
		}
	}

	else if (OPTION == "McElieceCryptosystem"){
		McElieceCryptosystem McEliece(n,m,g);
		int k = McEliece.GetGoppaCode().GetGeneratorMatrix().NumRows();
		mat_GF2 message = GetRandomMessage(k);
		cout << "random message: "  << endl;
		
		mat_GF2 ciphertext = McEliece.Encrypt(message);
		cout << "ciphertext: "  << endl;

		mat_GF2 recoveredtext1 = McEliece.Decrypt(ciphertext);
		cout << "recovered message1: " << recoveredtext1 << endl;
		if (recoveredtext1 == message)
		        cout << "It works!" << endl;
		else
		        cout << "Something wrong!" << endl;
	}

	else if (OPTION == "NiederreiterCryptosystem"){
		NiederreiterCryptosystem Niederreiter(n,m,g_square,delta); 
		//output the parity check matrix to the file
		const char *path="/home/david/Dropbox/GoppaCode/Code_Based_Cryptography_Cplusplus/src/Data/NiederreiterPublicKey.txt";
		//const char *path= "/Users/David/Dropbox/GoppaCode/Code_Based_Cryptography_Cplusplus/src/Data/NiederreiterPublicKey.txt";
		ofstream PublicKeyMatrix;
		PublicKeyMatrix.open(path);	
		/*mat_GF2 message = GetRandomMessageWithWeight(n,deg(g_square)/1);
		getchar();
		cout << "random message: " << message << endl;

		vector<int> message_index;
		for (int i = 0; i < message.NumCols(); i++)
			if (message[0][i] == 1)
				message_index.push_back(i);
			
		if (PublicKeyMatrix.is_open()){	
			PublicKeyMatrix << "message_weight: " << message_index.size() << endl;
			PublicKeyMatrix << "message_index: " << endl;
			for (int i = 0; i < message_index.size(); i++)
				PublicKeyMatrix << message_index[i] << "\t";
			PublicKeyMatrix << endl;
			//PublicKeyMatrix << message << endl;

			PublicKeyMatrix << "ParityCheckMatrix Transposed: " << endl;
			PublicKeyMatrix << transpose(Niederreiter.GetGoppaCode().GetParityCheckMatrix()) << endl;
			PublicKeyMatrix.close();

		} 
		else {
			cout << "file not saved." << endl;
		}


		mat_GF2 ciphertext = Niederreiter.Encrypt(message);
		cout << "ciphertext: " << ciphertext << endl;
		getchar();

		mat_GF2 recoveredtext1 = Niederreiter.Decrypt(ciphertext,"Patterson");
		//cout << "recovered message1: " << recoveredtext1 << endl;
		if (recoveredtext1 == message)
		        cout << "It works!" << endl;
		else
		        cout << "Something wrong!" << endl;

		getchar();  
		*/  

		//mat_GF2 message_hash = GetRandomMessage(m*deg(g_square));//length:m*t
		mat_GF2 message_hash = GetRandomMessageWithWeight(n,10)*transpose(Niederreiter.GetGoppaCode().GetParityCheckMatrix());
		cout << "hashed message: " << message_hash << endl;
		if(PublicKeyMatrix.is_open()){
			PublicKeyMatrix << "hashed message:" << endl;
			PublicKeyMatrix << message_hash << endl;
		}
		PublicKeyMatrix.close();
		getchar();
		/*Timer starts*/	
		long beginTime = clock();
		mat_GF2 signature = Niederreiter.Signature(message_hash,"random");
		cout << "signature: " << signature << endl;

		vector<int> signature_index;
		for (int i = 0; i < signature.NumCols(); i++)
			if (signature[0][i] == 1)
				signature_index.push_back(i);
		PublicKeyMatrix.open(path,ios::app);	
		PublicKeyMatrix << "signature index: " << signature_index.size() << endl;
		for (int i = 0; i < signature_index.size(); i++)
			PublicKeyMatrix << signature_index[i] << endl;	
		
		bool verification = Niederreiter.Verify(signature, message_hash);
		cout << "Signature: " << verification << endl;
		/*Timer ends*/        
		long endTime = clock();
		PublicKeyMatrix.close();
		cout << "Decryption Time Elpased(seconds): " << double(endTime-beginTime)/CLOCKS_PER_SEC << endl;
	}

	return 0;        
}




 /* Assistant Functions */

mat_GF2 GetRandomMessage(int message_length)
{
	mat_GF2 message;
	srand(777);
	message.SetDims(1,message_length);

	message[0][34] = 1;
	message[0][35] = 1;
	message[0][36] = 1;
	message[0][31] = 1;
	message[0][91] = 1;
	message[0][3] = 1;
	message[0][7] = 1;
	message[0][5] = 1;
	message[0][19] = 1;
	message[0][20] = 1;
	/*int message_weight = rand()%message_length;
	
	vector<int> rng(message_length); 
	for(int i = 0; i < message_length; i++)
		rng[i] = i;

	for (int i = 0; i < message_weight; i++){
		int p = rng.size()*(rand()/double(RAND_MAX));
		message[0][rng[p]] = 1;

		for(int j = p+1; j < rng.size(); j++)
			rng[j-1] = rng[j];
		rng.resize(--message_length);
		
	}*/
	return message;
}

mat_GF2 GetRandomMessageWithWeight(int message_length, int message_weight)
{
	if(message_length < message_weight){
		cout <<	"mat_GF2 GetRandomMessageWithWeight(int message_length, int message_weight): Error! Too much message_weight." << endl;
		exit(-1);	
	}
	mat_GF2 message;
	message.SetDims(1,message_length);
	srand(73);
	message[0][9036] = 1;
	message[0][8431] = 1;
	message[0][36991] = 1;
	message[0][55007] = 1;
	message[0][10055] = 1;
	message[0][5619] = 1;
	message[0][20992] = 1;
	message[0][32792] = 1;
	message[0][49819] = 1;

	
	/*vector<int> rng(message_length); 
	for(int i = 0; i < message_length; i++)
		rng[i] = i;

	for (int i = 0; i < message_weight; i++){
		int p = rng.size()*(rand()/double(RAND_MAX));
		message[0][rng[p]] = 1;
		cout << "random message index " << rng[p] << endl;
		for(int j = p+1; j < rng.size(); j++)
			rng[j-1] = rng[j];
		rng.resize(--message_length);
		
	}*/
	return message;
}




