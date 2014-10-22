#include "NiederreiterCryptosystem.h"

NiederreiterCryptosystem::NiederreiterCryptosystem(int n, int m, GF2EX g, int delta)
{
	//Construct the Goppa code
	GoppaCode goppacode(n,m,g);
	int k = goppacode.GetGeneratorMatrix().NumRows();
	if(k != n-m*deg(g)){
		cout << "NiederreiterCrytosystem::NiederreiterCryptosystem(int n, int m, GF2EX g): Error! Incorrect Goppa code generator matrix." << endl;
		exit(-1);	
	}

	//Set up a random scrambler matrix([n-k]*[n-k])
	mat_GF2 S; S.SetDims(n-k,n-k);
	for(int i = 0; i < n-k; i++)
		for(int j = 0; j < n-k; j++)
			S[i][j] = (rand()/double(RAND_MAX))>0.5 ? 1 : 0;
        mat_GF2 S_saved = S;
	while(gauss(S_saved) < n-k){
		S[(n-k)*(rand()/double(RAND_MAX))][(n-k)*(rand()/double(RAND_MAX))] += 1;
                S_saved = S;
         }
	
	//Set up the permutation matrix
	mat_GF2 P; 
	int P_len = n; P.SetDims(P_len,P_len);
	vector<int> rng(n);
	for(int i = 0; i < P_len; i++)
		rng[i] = i;
	for (int i = 0; i < n; i++){
		int p = rng.size()*(rand()/double(RAND_MAX));
		P[i][rng[p]] = 1;
		for(int j = p+1; j < rng.size(); j++)
			rng[j-1] = rng[j];
		rng.resize(--P_len);
		
	}
	

	//Remember these values
	this->_GoppaCode = goppacode;
	this->_g = g;
	this->_m = m;
	this->_n = n;
	this->_delta = delta;
	this->_P = P;
	this->_S = S;
	this->_PublicKey = S*goppacode.GetParityCheckMatrix()*P;
	
}

mat_GF2 NiederreiterCryptosystem::Encrypt(mat_GF2 message)
{
	if(message.NumRows() !=1 || message.NumCols() != this->_PublicKey.NumCols()){
		cout << "mat_GF2 NiederreiterCrytosystem::Encrypt(mat_GF2 message): Error! Message is not of correct length." << endl;
		exit(-1);
	}
	//Note we output the ciphertext as a row vector whereas the original one outputs its transpose.	
	return transpose(this->_PublicKey*transpose(message));	
}

mat_GF2 NiederreiterCryptosystem::Decrypt(mat_GF2 received_word)
{
	received_word = transpose(received_word);
	if(received_word.NumRows() != this->_PublicKey.NumRows() || received_word.NumCols() != 1){
		cout << "mat_GF2 NiederreiterCryptosystem::Decrypt(mat_GF2 received_word): Error! Recieved word is not of correct dimension." << endl;
		exit(-1);	
	}	
	//Strip off the random scrambler matrix and decode the received word.
	received_word = inv(this->_S)*received_word;
	
	//Convert received_word([n-k]*[1]) into syndrome polynomial
	//To achieve this, received_word is first converted into a [t]*[m] matrix
	//Then convert this matrix to vec_GF2 -> GF2X -> GF2E.
	int t = deg(this->_g); int m = this->_m;
	mat_GF2 tmp; tmp.SetDims(t,m);
	for (int i = 0; i < t; i++)
		for (int j = 0; j < m; j++)
			tmp[i][j] = received_word[i*m+j][0];
	GF2EX syndrome_poly;
	for (int i = 0; i < t; i++)
		SetCoeff( syndrome_poly, i, conv<GF2E> (conv<GF2X>(tmp[i])) );	
	
	//Retrieve using syndrome decoding algorithm. 		
	received_word = this->_GoppaCode.SyndromeDecode(syndrome_poly);
		
	//Solve the system to determine the original message.
	received_word = received_word*this->_P;
	
	return received_word;					
	
}

mat_GF2 NiederreiterCryptosystem::Signature(mat_GF2 message, string mode)
{
	
	if(message.NumCols() != this->_m*deg(this->_g)){
		cout << "mat_GF2 NiederreiterCryptosystem::Signature(mat_GF2 message, string mode): Error! message length incorrect." << endl;
		exit(-1);
	}
	//Decrypt a message of length m*t using complete docoding algorithm, here we try to decode a 't+delta' code.
	mat_GF2 new_message = message;
	
	if (mode == "random"){
		//randomly select columns from PublicKey Matrix	and add these columns to message.
		//The indexes of these columns are stored in random_pos_list.
		int cnt = 0;
		while(1){
			
			cout << "signature " << ++cnt << endl;	
			new_message = message;
			vector<int> random_pos_list(this->_delta);
			vector<int> rng(this->_n);
			for(int i = 0; i < this->_n; i++)
				rng[i] = i;
			int rng_len = rng.size();
			
			for (int i = 0; i < this->_delta; i++){
				int p = rng.size()*(rand()/double(RAND_MAX));
				random_pos_list[i] = p;
				
				for(int j = 0; j < new_message.NumCols(); j++)
					new_message[0][j] += this->_PublicKey[j][p];
			
				for(int j = p+1; j < rng.size(); j++)
					rng[j-1] = rng[j];
			
				rng.resize(--rng_len);		
			}
			
			
			new_message = this->Decrypt(new_message);
		
			if(GetRowVectorWeight(new_message) == deg(this->_g)){
			//A decodable hash value is found. In this circumstance, the assumption of the above randomly selected
			//columns are all correct and we must strip off the indexes of these columns from new_message.		
				for (int i = 0; i < random_pos_list.size(); i++)
					new_message[0][random_pos_list[i]] += 1;
				break;
				
			}
			
			else
				cout << "not found, weight: " << GetRowVectorWeight(new_message) << endl;
			
					
		}
		return new_message;
	}

	else if (mode == "sequential"){
		vector<int> rng(this->_n);
		for(int i = 0; i < this->_n; i++)
			rng[i] = i;
		vector<vector<int> > cmb_list;
		combination(cmb_list, rng, this->_n, this->_delta);
		int dim1 = cmb_list.size(); int dim2 = cmb_list[0].size();
		int cnt = 0;
		for (int i = 0; i < dim1; i++){
			cout << "signature " << ++cnt << endl;			
			new_message = message;
			for (int j = 0; j < dim2; j++)
				for (int k = 0; k < new_message.NumCols(); k++)							
					new_message[0][k] += this->_PublicKey[k][cmb_list[i][j]];
								
			new_message = this->Decrypt(new_message);

			if(GetRowVectorWeight(new_message) == deg(this->_g)){
				for (int ii = 0; ii < dim2; ii++)
					new_message[0][cmb_list[i][ii]] += 1;
				cout << "index: ";
				for (int ii = 0; ii < dim2; ii++)
					cout << cmb_list[i][ii] << " ";
				cout << endl;
				return new_message;
				
			}
	
			else
				cout << "not found, weight: " << GetRowVectorWeight(new_message) << endl; 					
		}
		
		cout << "Signature failed!" << endl;
		exit(-1);	
	}	
}


bool NiederreiterCryptosystem::Verify(mat_GF2 signature, mat_GF2 message)
{
	//Encrypt signature and then compare it with the document's hash value.
	if (this->Encrypt(signature) == message)
			//This is a valid signature.
			return true;
	else 
		//This is an invalid signature.			
		return false;	
}
		
/*Accessors*/
mat_GF2 NiederreiterCryptosystem::GetPublicKey(void)
{
		return this->_PublicKey;
}

GoppaCode NiederreiterCryptosystem::GetGoppaCode(void)
{
		return this->_GoppaCode;
}



