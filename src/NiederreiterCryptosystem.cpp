#include "NiederreiterCryptosystem.h"
#define K_LIMIT 100000


NiederreiterCryptosystem::NiederreiterCryptosystem(int n, int m, GF2EX g, int delta)
{
	//Construct the Goppa code
	GoppaCode goppacode(n,m,g);
	//int k = goppacode.GetGeneratorMatrix().NumRows();
	int k = n - m*deg(g);
	if (k < n-m*deg(g)){
		cout << "NiederreiterCrytosystem::NiederreiterCryptosystem(int n, int m, GF2EX g): Error! Incorrect Goppa code generator matrix." << endl;
		exit(-1);	
	}
	if (k < K_LIMIT){
		//Set up a random scrambler matrix([mt]*[mt])
		mat_GF2 S = ScramblerMatrixGen(m*deg(g), "fixed");
	
		//Set up the permutation matrix([n]*[n])
		mat_GF2 P = PermutationMatrixGen(n, "fixed"); 
	
		//Remember these values
		this->_GoppaCode = goppacode;
		this->_g = g;
		this->_m = m;
		this->_n = n;
		this->_delta = delta;
		this->_P = P;
		this->_S = S;

		cout << "n=" << n << endl;
		cout << "k=" << k << endl;
		cout << "n-k=" << n-k << endl;
		cout << "H row" << goppacode.GetParityCheckMatrix().NumRows() << endl;
		cout << "H col" << goppacode.GetParityCheckMatrix().NumCols() << endl;
		this->_PublicKey = S*goppacode.GetParityCheckMatrix()*P;
		cout << "End of constructor1" << endl;	
	}
	else {
		//Set up a random scrambler matrix([mt]*[mt])
		Mat<ZZ> S = LargeScramblerMatrixGen(m*deg(g), "fixed");
	
		//Set up the permutation matrix([n]*[n])
		Mat<ZZ> P = LargePermutationMatrixGen(n,"fixed"); 
	
		//Remember these values
		this->_GoppaCode = goppacode;
		this->_g = g;
		this->_m = m;
		this->_n = n;
		this->_delta = delta;
		this->_LargeP = P;
		this->_LargeS = S;

		cout << "n=" << n << endl;
		cout << "k=" << k << endl;
		cout << "n-k=" << n-k << endl;
		cout << "H row" << goppacode.GetParityCheckMatrix().NumRows() << endl;
		cout << "H col" << goppacode.GetParityCheckMatrix().NumCols() << endl;
		// S*H*P
		int S_dim = S.NumRows();
		int P_dim = P.NumCols();
		//S*H
		this->_PublicKey = goppacode.GetParityCheckMatrix();

		for(int i = 0; i < S_dim; i++){
			int index1 = conv<int>(S[i][0]);
			int index2 = conv<int>(S[i][1]); 

			if (index1 != 0 || index2 != 0)
				this->_PublicKey[index2] += this->_PublicKey[index1];				
				
		}
		
		//(S*H)*P
		for(int i = 0; i < P_dim; i++){
			int index1 = conv<int>(P[0][i]);
			int index2 = conv<int>(P[1][i]);
			GF2 tmp;
			if (index1 != 0 || index2 != 0)
				for(int j = 0; j < S_dim; j++){
					tmp = this->_PublicKey[j][index1];
					this->_PublicKey[j][index1] = this->_PublicKey[j][index2];
					this->_PublicKey[j][index2] = tmp;
				}
		}
		cout << "End of constructor2" << endl;	
	}
	
	
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

mat_GF2 NiederreiterCryptosystem::Decrypt(mat_GF2 received_word, string mode)
{
	received_word = transpose(received_word);
	if(received_word.NumRows() != this->_PublicKey.NumRows() || received_word.NumCols() != 1){
		cout << "mat_GF2 NiederreiterCryptosystem::Decrypt(mat_GF2 received_word): Error! Recieved word is not of correct dimension." << endl;
		exit(-1);	
	}	
	//int k = this->_GoppaCode.GetGeneratorMatrix().NumRows();
	int k = K_LIMIT + 1;
	if (k < K_LIMIT) {
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
		received_word = this->_GoppaCode.SyndromeDecode(syndrome_poly,mode);
			
		//Solve the system to determine the original message.
		received_word = received_word*this->_P;
	}

	else {	
		//Strip off the random scrambler matrix and decode the received word.
		int S_dim = this->_LargeS.NumRows();
		for (int i = 0; i < S_dim; i++){
			int index1 = conv<int>(this->_LargeS[S_dim-1-i][0]);
			int index2 = conv<int>(this->_LargeS[S_dim-1-i][1]);
			if (index1 != 0 || index2 != 0)
				received_word[index2] += received_word[index1];	
		}		
		long endTime1 = clock();
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
		received_word = this->_GoppaCode.SyndromeDecode(syndrome_poly,mode);
			
		//Solve the system to determine the original message.
		//received_word = received_word*this->_P;
		int P_dim = this->_LargeP.NumCols();
		for(int i = 0; i < P_dim; i++){
			int index1 = conv<int>(this->_LargeP[0][i]);
			int index2 = conv<int>(this->_LargeP[1][i]);
			GF2 tmp;
			if (index1 != 0 || index2 != 0)
				for(int j = 0; j < 1; j++){
					tmp = received_word[j][index1];
					received_word[j][index1] = received_word[j][index2];
					received_word[j][index2] = tmp;
				}
		}

	}
	
	
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
	mat_GF2	paritycheckmatrix = this->GetGoppaCode().GetParityCheckMatrix();
	if (mode == "random"){
		//randomly select columns from PublicKey Matrix	and add these columns to message.
		//The indexes of these columns are stored in random_pos_list.
		srand(1);
		int cnt = 0;
		SetPrngSeed(100);
		while(1){			
			cout << "signature " << ++cnt << endl;	
			new_message = message;
			vector<int> random_pos_list(this->_delta);
			vector<int> rng(this->_n);
			for(int i = 0; i < this->_n; i++)
				rng[i] = i;
			int rng_len = rng.size();			

			for (int i = 0; i < this->_delta; i++){
				//int p = rand()%rng_len;
				int p = PrngGen("LSFR")%rng_len;
				random_pos_list[i] = p;
				cout << "p = " << p << endl;							

				for(int j = 0; j < new_message.NumCols(); j++){
					new_message[0][j] += paritycheckmatrix[j][p];
				}	
				/*for(int j = p+1; j < rng.size(); j++)
					rng[j-1] = rng[j];
			
				rng.resize(--rng_len);*/		
			}

			if (IterIrredTest(this->_g)){
				mat_GF2 origin_message = new_message;
				cout << "corrected syndrome: " << origin_message << endl;
				//convert origin_message into syndrome polynomial.
				int t = deg(this->_g); int m = this->_m;
				mat_GF2 tmp; tmp.SetDims(t,m);
				for (int i = 0; i < t; i++)
					for (int j = 0; j < m; j++)
						tmp[i][j] = origin_message[0][i*m+j];
				GF2EX syndrome_poly;
				for (int i = 0; i < t; i++)
					SetCoeff( syndrome_poly, i, conv<GF2E> (conv<GF2X>(tmp[i])) );	
				new_message = this->GetGoppaCode().SyndromeDecode(syndrome_poly,"Patterson");	
					
				//new_message = this->Decrypt(origin_message, "Patterson");		
				getchar();
				if(GetRowVectorWeight(new_message) <= deg(this->_g)){					
					cout << "message weight " << GetRowVectorWeight(new_message) << endl;
					vector<int> new_message_index;
					for (int i = 0; i < new_message.NumCols(); i++)
						if (new_message[0][i] == 1){
							new_message_index.push_back(i);
						}
					for (int i = 0; i < new_message_index.size(); i++)	
						cout << "message weight index " << new_message_index[i] << endl;
					cout << "c: " << new_message << endl;

				//A decodable hash value is found. In this circumstance, the assumption of the above randomly selected
				//columns are all correct and we must strip off the indexes of these columns from new_message.		
					for (int i = 0; i < random_pos_list.size(); i++){
						new_message[0][random_pos_list[i]] += 1;
						cout << "index corrected: " << random_pos_list[i] << endl;	
					}
					getchar();	
					return new_message;				
				}
			
				else
					cout << "not found, weight: " << GetRowVectorWeight(new_message) << endl;
			}
				
			else {
				mat_GF2 origin_message = new_message;

				new_message = this->Decrypt(origin_message, "Euclidean");	

				mat_GF2 Hc1;
				Hc1 = this->Encrypt(new_message);

				if (GetRowVectorWeight(new_message) != this->_n){
					cout << "WEIGHT: " << GetRowVectorWeight(new_message);
					getchar();
				}

				if(GetRowVectorWeight(new_message) <= deg(this->_g)/2){
					if (Hc1 == origin_message)
						cout << "no affect before or after the decryption." << endl;
					else
						cout << "affected before or after the decryption." << endl;

					cout << "message weight " << GetRowVectorWeight(new_message) << endl;
					vector<int> new_message_index;
					for (int i = 0; i < new_message.NumCols(); i++)
						if (new_message[0][i] == 1){
							new_message_index.push_back(i);
						}
					for (int i = 0; i < new_message_index.size(); i++)	
						cout << "message weight index " << new_message_index[i] << endl;
					cout << "c: " << new_message << endl;

					
					

					mat_GF2 Hdelta1;
					mat_GF2 delta; 
					delta.SetDims(1,this->_PublicKey.NumCols());
					for (int i = 0; i < random_pos_list.size(); i++)
						delta[0][random_pos_list[i]] = 1;
					Hdelta1 = this->Encrypt(delta);
					cout << "Hdelta" << Hdelta1 << endl;

					//mat_GF2 Hdelta = message-Hc1;
					mat_GF2 Hdelta = transpose(this->_PublicKey*transpose(delta));
					
					

					mat_GF2 delta1 = this->Decrypt(Hc1, "Euclidean");
					std::vector<int> delta1_index;
					for (int i = 0; i < delta1.NumCols(); i++)
						if (delta1[0][i] == 1){
							delta1_index.push_back(i);
						}

					cout << "delta: " <<  delta1_index.size() << endl;
					/*for (int i = 0; i < delta1_index.size(); i++)	
						cout << "delta index " << delta1_index[i] << endl;*/
					cout << "S-Hc-Hdelta: " << message-Hdelta1-Hc1 << endl;
					/*mat_GF2 Diff = this->Encrypt(this->Decrypt(message-Hdelta1-Hc1, "Euclidean"));
					std::vector<int> Diff_index;
					for (int i = 0; i < Diff.NumCols(); i++)
						if (Diff[0][i] == 1)
							Diff_index.push_back(i);
					cout << "Diff_index: " << Diff_index.size() << endl;
					for (int i = 0; i < Diff_index.size(); i++)
						cout << Diff_index[i] << endl;
					cout << "Diff: " << Diff << endl;*/
					


					//A decodable hash value is found. In this circumstance, the assumption of the above randomly selected
					//columns are all correct and we must strip off the indexes of these columns from new_message.		
					for (int i = 0; i < random_pos_list.size(); i++){
						new_message[0][random_pos_list[i]] += 1;
						cout << "index corrected: " << random_pos_list[i] << endl;	
					}
					vector<int> signature_index;
					for (int i = 0; i < new_message.NumCols(); i++)
						if (new_message[0][i] == 1){
							signature_index.push_back(i);
						}
					for (int i = 0; i < signature_index.size(); i++)	
						cout << "signature weight index " << signature_index[i] << endl;
					return new_message;	
				}
			
				else
					cout << "not found, weight: " << GetRowVectorWeight(new_message) << endl;
			}
			
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

		for (int i = 0; i < dim1; i++)
			cout << "cmb_list " << cmb_list[i][0] << endl;
		getchar();
		int cnt = 0;
		for (int i = 0; i < dim1; i++){
			cout << "signature " << ++cnt << endl;			
			new_message = message;
			for (int j = 0; j < dim2; j++)
				for (int k = 0; k < new_message.NumCols(); k++)							
					new_message[0][k] += this->_PublicKey[k][cmb_list[i][j]];
								
			if (IterIrredTest(this->_g)){
				new_message = this->Decrypt(new_message, "Patterson");

				if (GetRowVectorWeight(new_message) <= deg(this->_g)){
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

			else {
				new_message = this->Decrypt(new_message, "Euclidean");	
				
				if (GetRowVectorWeight(new_message) <= deg(this->_g)/2){
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
							
		}
		
		cout << "Signature failed!" << endl;
		exit(-1);	
	}	
}


bool NiederreiterCryptosystem::Verify(mat_GF2 signature, mat_GF2 message)
{
	//Encrypt signature and then compare it with the document's hash value.
	//cout << "Verify " << this->Encrypt(signature);
	cout << "Verfiy" << signature*transpose(this->GetGoppaCode().GetParityCheckMatrix()) << endl;
	cout << "Verify " << message << endl;
	if (this->Encrypt(signature) == message)
	//if (signature*transpose(this->GetGoppaCode().GetParityCheckMatrix()) == message)
			//This is a valid signature.
		return true;
	else 
		//This is an invalid signature.			
		return false;	
}

mat_GF2 NiederreiterCryptosystem::ScramblerMatrixGen(int k, string mode)
{
	mat_GF2 S; 
	if (mode == "random"){
		//generate a k*k scrambler matrix S
		// if S is small, a fully randomized matrix is constructed
		S.SetDims(k,k);
		for(int i = 0; i < k; i++)
		for(int j = 0; j < k; j++)
			S[i][j] = (rand()/double(RAND_MAX))>0.5 ? 1 : 0;
	    mat_GF2 S_saved = S;

	 	int rank = 0;
		while((rank = gauss(S_saved)) < k){
			S[k*(rand()/double(RAND_MAX))][k*(rand()/double(RAND_MAX))] += 1;
	        S_saved = S;
	    }
	}

	else if (mode == "fixed"){
		//generate a k*k scrambler matrix S
		// a fixed pattern is executed: r0->r1, r1->r2, r2->r3, r3->r4, r4->r5
		S.SetDims(k,k);
		for (int i = 0; i < k; i++)
			S[i][i] = 1;

		/*for (int i = 0; i < 5; i++)
			S[i+1] += S[i];*/ 

	}
	
	else {
		cout << "NiederreiterCryptosystem::ScramblerMatrixGen(int k, string mode=\"random\"): Error! Cannot find matched input mode." << endl;	
		exit(-1);	
	}
	
	return S;
}

Mat<ZZ> NiederreiterCryptosystem::LargeScramblerMatrixGen(int k, string mode)
{
	Mat<ZZ> S;
	//if S is too large, to avoid out of memory error and save time, a much simpler randomized matrix is constructed.
	S.SetDims(k,2);
	if (mode == "random"){
		for (int cnt = 0; cnt < k; cnt++){ // a number of k linear transformation are conducted
			S[cnt][0] = conv<ZZ>(k*(rand()/double(RAND_MAX)));
			S[cnt][1] = conv<ZZ>(k*(rand()/double(RAND_MAX)));
		}	
	}
	else if (mode == "fixed"){
			for (int i = 0; i < 5; i++){
				S[i][0] = i;
				S[i][1] = i+1;
			}
	}
	else {
		cout << "NiederreiterCryptosystem::LargeScramblerMatrixGen(int k, string mode=\"random\"): Error! Cannot find matched input mode." << endl;	
		exit(-1);
	}
	return S;
}

mat_GF2 NiederreiterCryptosystem::PermutationMatrixGen(int n, string mode)
{
	//generate a n*n permutation matrix P
	mat_GF2 P; 
	int P_len = n; 

	P.SetDims(P_len,P_len);
	if (mode == "random"){
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
	}

	else if (mode == "fixed"){
		//generate a n*n permutation matrix P
		//a fixed pattern is executed: c0<->c1, c1<->c2, c2<->c3, c3<->c4, c4<->c5
		for (int i = 0; i < n; i++)
			P[i][i] = 1;
		/*GF2 tmp;
		for (int i = 0; i < 5; i++){
			for (int j = 0; j < n; j++){
				tmp = P[j][i];
				P[j][i] = P[j][i+1];
				P[j][i+1] = tmp;
			}
		}*/
	}
	
	else {
		cout << "NiederreiterCryptosystem::PermutationMatrixGen(int n, string mode=\"random\"): Error! Cannot find matched input mode." << endl;	
		exit(-1);
	}
	return P;
}

Mat<ZZ> NiederreiterCryptosystem::LargePermutationMatrixGen(int n, string mode)
{
	Mat<ZZ> P; 
	int P_len = n; 
	P.SetDims(2,P_len);
	if (mode == "random"){
		for(int i = 0; i < P_len; i++){
		int rand1 = P_len*(rand()/double(RAND_MAX));
		int rand2 = P_len*(rand()/double(RAND_MAX));
		P[0][i] = conv<ZZ>(rand1);
		P[1][i] = conv<ZZ>(rand2);
		}
	}
	
	else if (mode == "fixed"){
			for (int i = 0; i < 5; i++){
				P[0][i] = i;
				P[1][i] = i+1;
			}	
	}

	else {
		cout << "NiederreiterCryptosystem::LargePermutationMatrixGen(int n, string mode=\"random\"): Error! Cannot find matched input mode." << endl;	
		exit(-1);
	}

	return P;
}

		
/*Accessors*/
mat_GF2 NiederreiterCryptosystem::GetPublicKey(void)
{
		return this->_PublicKey;
}

mat_GF2 NiederreiterCryptosystem::GetP(void)
{

	return this->_P;
}

mat_GF2 NiederreiterCryptosystem::GetS(void)
{

	return this->_S;
}

GoppaCode NiederreiterCryptosystem::GetGoppaCode(void)
{
		return this->_GoppaCode;
}

/*Helpers*/

void NiederreiterCryptosystem::SetPrngSeed(unsigned int seed){
	this->seed = seed;
}
int NiederreiterCryptosystem::PrngGen(string mode){
	if (mode == "LSFR"){
		//The generating polynomial is x^25 + x^3 + 1
		//this->seed = ((((this->seed >> 2) ^ (this->seed >> 24)) & 0x01) | (this->seed << 1)) & 0x7fffffff;
		//The generating  polynomial is x^31 + x^3 + 1
		this->seed = ((((this->seed >> 2) ^ (this->seed >> 30)) & 0x01) | (this->seed << 1)) & 0x7fffffff;
		return this->seed;
	}
	else if (mode == "LCG"){
		//% 2147483648u) & 0x7fffffff
		/*this->seed = 1103515245 * this->seed + 12345;
		return this->seed & 0x7fffffff;*/
		this->seed = this->seed * 1103515245 + 12345; 	 
		return (unsigned int)(this->seed >> 1);
		//return (unsigned int)(this->seed/65536) % 32768;		
	}
	else{
		cout << "Incorrect PRNG mode!" << endl;
		exit(-1);
	}
	
}

