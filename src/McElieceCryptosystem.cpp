#include "McElieceCryptosystem.h"

McElieceCryptosystem::McElieceCryptosystem(int n, int m, GF2EX g)
{
	//Construct Goppa code
	GoppaCode goppa_code(n,m,g);
	int k = goppa_code.GetGeneratorMatrix().NumRows();
	if(!k){

		cout << "McElieceCryptosystem::McElieceCryptosystem(int m, int n, GF2EX g): Error! Generator Matrix is empty." << endl;	
		exit(-1);	
	}
	//Set up a random scrambler matrix
	mat_GF2 S; S.SetDims(k,k);
	for(int i = 0; i < k; i++)
		for(int j = 0; j < k; j++)
			S[i][j] = (rand()/double(RAND_MAX))>0.5 ? 1 : 0;
        mat_GF2 S_saved = S;
	while(gauss(S_saved) < k){
		S[k*(rand()/double(RAND_MAX))][k*(rand()/double(RAND_MAX))] += 1;
                S_saved = S;
         }
	cout << "McElieceCryptosystem: Scrambler: " <<  endl;
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
	//cout << "McElieceCryptosystem: Permutation: " << P;

	//Remember these values
	this->_GoppaCode = goppa_code;
	this->_g = g;
	this->_m = m;
	this->_n = n;
	this->_P = P;
	this->_S = S;
	this->_PublicKey = S*goppa_code.GetGeneratorMatrix()*P;
	
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
	//Strip off the permutation and decode the received word.
	message *= inv(this->_P);
	message = this->_GoppaCode.Decode(message);
        
	//Solve the system to determine the original message.
	mat_GF2 SG; 
	SG.SetDims(this->_GoppaCode.GetGeneratorMatrix().NumRows(), this->_GoppaCode.GetGeneratorMatrix().NumCols());	
	SG = this->_S*this->_GoppaCode.GetGeneratorMatrix();
	
	int k = SG.NumRows();
        int n = SG.NumCols();
        mat_GF2 SG_tran = transpose(SG);
        mat_GF2 SG_tran_augment; SG_tran_augment.SetDims(n, k+1);
        for (int i = 0; i < SG_tran.NumRows(); i++)
                for (int j = 0; j < SG_tran.NumCols(); j++)
                        SG_tran_augment[i][j] = SG_tran[i][j];
        for (int i = 0; i < message.NumCols(); i++)
                SG_tran_augment[i][SG_tran.NumCols()] = message[0][i];
        gauss(SG_tran_augment);
        
	mat_GF2 tmp1; tmp1.SetDims(k,k);
	mat_GF2 tmp2; tmp2.SetDims(k,1);
	for(int i = 0; i < k; i++)
		for(int j = 0; j < k+1; j++)
			if(j != k)
				tmp1[i][j] = SG_tran_augment[i][j];
			else	
				tmp2[i][0] = SG_tran_augment[i][j];	
	tmp2 = inv(tmp1)*tmp2;
       

	
	return transpose(tmp2);
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






