#include "McElieceCryptosystem.h"
#define K_LIMIT 1

McElieceCryptosystem::McElieceCryptosystem(int n, int m, GF2EX g)
{
	//Construct Goppa code
	GoppaCode goppa_code(n,m,g);
	int k = goppa_code.GetGeneratorMatrix().NumRows();
	cout << "McElieceCryptosystem: code dimension k = " << k << endl;
	if(!k){
		cout << "McElieceCryptosystem::McElieceCryptosystem(int m, int n, GF2EX g): Error! Generator Matrix is empty." << endl;	
		exit(-1);	
	}

	if (k < K_LIMIT){
		//matrix way of doing it
		//Set up a random scrambler matrix
		mat_GF2 S = ScramblerMatrixGen(k,"fixed");
		cout << "McElieceCryptosystem: Scrambler S: " << endl;

		//Set up the permutation matrix
		mat_GF2 P = PermutationMatrixGen(n,"fixed"); 
		cout << "McElieceCryptosystem: Permutation P0 " << endl;

		this->_PublicKey = S*goppa_code.GetGeneratorMatrix()*P;
		
		//Remember these values
		this->_GoppaCode = goppa_code; 
		this->_g = g;
		this->_m = m;
		this->_n = n;
		this->_P = P;
		this->_S = S;
		cout << "End of constructor" << endl;
	}
		
	else {
		//linear transformation way of doing it
		//Set up a random scrambler matrix
		Mat<ZZ> S = LargeScramblerMatrixGen(k,"fixed");
		cout << "McElieceCryptosystem: Scrambler: " << endl;

		//Set up the permutation matrix
		Mat<ZZ> P = LargePermutationMatrixGen(n,"fixed"); 
		cout << "McElieceCryptosystem: Permutation: " << endl;

		int k = S.NumRows();
		int n = P.NumCols();
		//S*G
		this->_PublicKey = goppa_code.GetGeneratorMatrix();

		for(int i = 0; i < k; i++){
			int index1 = conv<int>(S[i][0]);
			int index2 = conv<int>(S[i][1]); 

			if (index1 != 0 || index2 != 0)
				this->_PublicKey[index2] += this->_PublicKey[index1];				
				
		}

		
		//(S*G)*P
		for(int i = 0; i < n; i++){
			int index1 = conv<int>(P[0][i]);
			int index2 = conv<int>(P[1][i]);
			GF2 tmp;
			if (index1 != 0 || index2 != 0)
				for(int j = 0; j < k; j++){
					tmp = this->_PublicKey[j][index1];
					this->_PublicKey[j][index1] = this->_PublicKey[j][index2];
					this->_PublicKey[j][index2] = tmp;
				}
		}
	

		//Remember these values
		this->_GoppaCode = goppa_code; 
		this->_g = g;
		this->_m = m;
		this->_n = n;
		this->_LargeP = P;
		this->_LargeS = S;
		cout << "End of constructor" << endl;
	}
	
}

mat_GF2 McElieceCryptosystem::ScramblerMatrixGen(int k, string mode)
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

		for (int i = 0; i < 5; i++)
			S[i+1] += S[i]; 
	}
	
	else {
		cout << "McElieceCryptosystem::ScramblerMatrixGen(int k, string mode=\"random\"): Error! Cannot find matched input mode." << endl;	
		exit(-1);	
	}
	
	return S;
}

Mat<ZZ> McElieceCryptosystem::LargeScramblerMatrixGen(int k, string mode)
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
		cout << "McElieceCryptosystem::LargeScramblerMatrixGen(int k, string mode=\"random\"): Error! Cannot find matched input mode." << endl;	
		exit(-1);
	}
	return S;
}

mat_GF2 McElieceCryptosystem::PermutationMatrixGen(int n, string mode)
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

			GF2 tmp;
			for (int i = 0; i < 5; i++){
				for (int j = 0; j < n; j++){
					tmp = P[j][i];
					P[j][i] = P[j][i+1];
					P[j][i+1] = tmp;
				}
			}
	}
	
	else {
		cout << "McElieceCryptosystem::PermutationMatrixGen(int n, string mode=\"random\"): Error! Cannot find matched input mode." << endl;	
		exit(-1);
	}
	return P;
}

Mat<ZZ> McElieceCryptosystem::LargePermutationMatrixGen(int n, string mode)
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
		cout << "McElieceCryptosystem::LargePermutationMatrixGen(int n, string mode=\"random\"): Error! Cannot find matched input mode." << endl;	
		exit(-1);
	}

	return P;
}

mat_GF2 McElieceCryptosystem::Encrypt(mat_GF2 message)
{
	//Ensure that the message is of the correct length(k).
	if (message.NumCols() != this->_PublicKey.NumRows()){
		cout << "mat_GF2 McElieceCryptosystem::Encrypt(mat_GF2 message): Error! message is not of the correct length." << endl;
		exit(-1);
	}
	//Get an error vector, ensuring that there are exactly t errors.
	int t = deg(this->_g), ncols =	this->_GoppaCode.GetGeneratorMatrix().NumCols();
	mat_GF2 err_vec; err_vec.SetDims(1,ncols);
	
	while (GetRowVectorWeight(err_vec) < t)
		err_vec[0][rand()%ncols] = 1;
	return message*this->_PublicKey + err_vec;
}

mat_GF2 McElieceCryptosystem::Decrypt(mat_GF2 message)
{
	//Ensure that the message is of the correct length(n).	
	if (message.NumCols() != this->_PublicKey.NumCols()){
		cout << "mat_GF2 McElieceCryptosystem::Decrypt(mat_GF2 message): Error! message is not of the correct length." << endl;
		exit(-1);		
	}
	int k = this->_GoppaCode.GetGeneratorMatrix().NumRows();
	int n = this->_GoppaCode.GetGeneratorMatrix().NumCols();
	
	if (this->_PublicKey.NumRows() < K_LIMIT){
		/*****STEP01*****/
		//Strip off the permutation and decode the received word.
		message *= inv(this->_P);
		
		/*****STEP02*****/
		message = this->_GoppaCode.Decode(message);
		
		/*****STEP03*****/
		//Solve the system to determine the original message.
		mat_GF2 G; 
		G = this->_GoppaCode.GetGeneratorMatrix();

	    mat_GF2 G_tran = transpose(G);
	    mat_GF2 G_tran_augment; G_tran_augment.SetDims(n, k+1);
	    for (int i = 0; i < G_tran.NumRows(); i++)
	            for (int j = 0; j < G_tran.NumCols(); j++)
	                    G_tran_augment[i][j] = G_tran[i][j];
	    for (int i = 0; i < message.NumCols(); i++)
	            G_tran_augment[i][G_tran.NumCols()] = message[0][i];
	    gauss(G_tran_augment);
	    

		mat_GF2 tmp1; tmp1.SetDims(k,k);
		mat_GF2 tmp2; tmp2.SetDims(k,1);
		for (int i = 0; i < k; i++)
			for (int j = 0; j < k+1; j++)
				if(j != k)
					tmp1[i][j] = G_tran_augment[i][j];
				else	
					tmp2[i][0] = G_tran_augment[i][j];	
		
		message = transpose(inv(tmp1)*tmp2);
		/*****STEP04*****/
		message = message * inv(this->_S);
		
	}
		
	else{
		//Strip off the permutation and decode the received word.
		cout << "STEP01" << endl;
		for(int i = 0; i < n; i++){
			GF2 tmp;
			int index1 = conv<int>(this->_LargeP[0][n-1-i]);
			int index2 = conv<int>(this->_LargeP[1][n-1-i]);
			if (index1 != 0 || index2 != 0){
				tmp = message[0][index1];
				message[0][index1] = message[0][index2];
				message[0][index2] = tmp;
			}
			
		}
		cout << "STEP02" << endl;
		message = this->_GoppaCode.Decode(message);

		//Solve the system to determine the original message.
		cout << "STEP03" << endl;
		mat_GF2 G; 
		G = this->_GoppaCode.GetGeneratorMatrix();
		
		int k = G.NumRows();
	    int n = G.NumCols();

	    mat_GF2 G_tran = transpose(G);
	    mat_GF2 G_tran_augment; G_tran_augment.SetDims(n, k+1);
	    for (int i = 0; i < G_tran.NumRows(); i++)
	            for (int j = 0; j < G_tran.NumCols(); j++)
	                    G_tran_augment[i][j] = G_tran[i][j];
	    for (int i = 0; i < message.NumCols(); i++)
	            G_tran_augment[i][G_tran.NumCols()] = message[0][i];
	    gauss(G_tran_augment);
	        
		mat_GF2 tmp1; tmp1.SetDims(k,k);
		mat_GF2 tmp2; tmp2.SetDims(k,1);
		for (int i = 0; i < k; i++)
			for (int j = 0; j < k+1; j++)
				if(j != k)
					tmp1[i][j] = G_tran_augment[i][j];
				else	
					tmp2[i][0] = G_tran_augment[i][j];	
		
		message = transpose(inv(tmp1)*tmp2);

		cout << "STEP04" << endl;
		for(int i = 0; i < k; i++){
			int index1 = conv<int>(this->_LargeS[i][1]);
			int index2 = conv<int>(this->_LargeS[i][0]);
			if (index1 != 0 || index2 != 0)
				message[0][index2] += message[0][index1];
		}
	}
	
    return message;    	
}



//Accessors
mat_GF2 McElieceCryptosystem::GetPublicKey(void)
{
	return this->_PublicKey;
}

GoppaCode McElieceCryptosystem::GetGoppaCode(void)
{
	return this->_GoppaCode;
}

mat_GF2 McElieceCryptosystem::GetScramblerMatrix(void)
{
	return this->_S;
}

mat_GF2 McElieceCryptosystem::GetPermutationMatrix(void)
{
	return this->_P;
}

Mat<ZZ> McElieceCryptosystem::GetLargeScramblerMatrix(void)
{
	return this->_LargeS;
}

Mat<ZZ> McElieceCryptosystem::GetLargePermutationMatrix(void)
{
	return this->_LargeP;
}



