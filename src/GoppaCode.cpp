#include "GoppaCode.h"

GoppaCode::GoppaCode()
{
	this->_n = 0;
	this->_m = 0;
	this->_g = conv<GF2EX>(0);
	this->_SyndromeCalculator = vector<GF2EX>(0);
	this->_codelocators = vector<GF2E>(0);
	
}

GoppaCode::GoppaCode(int n, int m, GF2EX g)
{	
	int t = deg(g);
	this->_n = n;
	this->_m = m;
	this->_g = g;

	/*find primitive root over GF(2^m), [0 0 0 1 0 0 0 1 1 0 1 1 1 0 0 1] by default*/
	GF2X primitive_root;
	SetCoeff(primitive_root,3,1);
	SetCoeff(primitive_root,7,1);
	SetCoeff(primitive_root,8,1);
	SetCoeff(primitive_root,10,1);
	SetCoeff(primitive_root,11,1);
	SetCoeff(primitive_root,12,1);
	SetCoeff(primitive_root,15,1);
	this->_Gf2eGenerator = conv<GF2E>(primitive_root);

	//this->_Gf2eGenerator = _GetGf2eGenerator(m);
	cout << "generator: " << this->_Gf2eGenerator << endl;
        /*Initialize the code locators*/
        _codelocators.resize(n);
        for (int i = 0; i < n-1; i++)
                _codelocators[i] = power(this->_Gf2eGenerator, conv<ZZ>(i+1));
        _codelocators[n-1] = GF2E::zero();
        /*for (int i = 0; i < _codelocators.size(); i++)
                cout << "cl: " << _codelocators[i] << endl;*/
	/*This is the same h to which Bernstein refers in his polynomial view of a Goppa code*/
	GF2EX h; h = 1; //h = 1
	GF2EX X; SetX(X); //X = 'X'
	/*Gamma is a list of the polynomials used to determine membership in the code*/
	vector<GF2EX> gamma(n); GF2EX inverse;
	vector<vector<GF2E> > H_check_poly(t, vector<GF2E>(n)); 
	for (int i = 0; i < n; i++){
		InvMod(inverse, X-_codelocators[i], g);
		gamma[i] = MulMod(h,inverse,g);
		//cout << "gamma: " << gamma[i] << endl;
		for (int j = 0; j < t; j++){
		/*Check to make sure the coefficient exists
		It may not if the polynomial is not of degree t*/
			if (j <= deg(gamma[i]))
				H_check_poly[j][i] = coeff(gamma[i],j);
			else
				H_check_poly[j][i] = 0;
			//cout << "H_check_poly: " << H_check_poly[j][i]; 		
		}		
	}

	_H_Goppa_Poly.SetDims(t, n);
	for (int i = 0; i < t; i++)
		for(int j = 0; j < n; j++){
			_H_Goppa_Poly[i][j] = H_check_poly[i][j];	
	}

	/*Construct the binary parity-check matrix for the Goppa code.
	Do so by converting each element of F_2^m to its binary representation*/
	//vector<vector<GF2> > vec_vec_H_Goppa(m*H_check_poly.size(), vector<GF2>(H_check_poly[0].size()));
	//mat_GF2 _H_Goppa;
	_H_Goppa.SetDims(m*H_check_poly.size(), H_check_poly[0].size()); 
	for (int i = 0; i < H_check_poly.size(); i++)
		for (int j = 0; j < H_check_poly[i].size(); j++){
			GF2X tmp = conv<GF2X>(H_check_poly[i][j]);
			int coeff_deg = deg(tmp);
			for (int k = i*m; k < (i+1)*m; k++){
				_H_Goppa[k][j] = coeff(tmp,k%m);
				//cout << _H_Goppa[k][j] << " ";			
			}
		}

	//cout << "__Goppa " << _H_Goppa << endl;
	/*Construct the generator matrix for our code by computinga basis for the null-space of H_Goppa. 
	The null-space is, by definition, the codewords of our code.*/
	mat_GF2 G_Goppa; 
	//kernel(_G_Goppa, transpose(_H_Goppa));
	//gauss(_G_Goppa);
	/*if (_G_Goppa.NumRows() != n-m*deg(g)){
		cout << "GoppaCode::GoppaCode(): G_Goppa("<< _G_Goppa.NumRows() <<") dimension error." << endl;
		exit(-1);
	}*/
	//cout << "_G_Goppa" << _G_Goppa << endl; PrintResult(_G_Goppa,"./matrix.txt");
	/*Construct the syndrome calculator. This will be used
	to simplify the calculation of syndromes for decoding.*/
	int codelocators_len = _codelocators.size();
	_SyndromeCalculator.resize(codelocators_len);
	for (int i = 0; i < codelocators_len; i++)
		_SyndromeCalculator[i] = InvMod(X - _codelocators[i], g);	

}
		

mat_GF2 GoppaCode::Encode(mat_GF2 message)
{      
	if(message.NumRows() != 1 || message.NumCols() < this->_n-this->_m*deg(this->_g)){
		cout << "mat_GF2 GoppaCode::Encode(mat_GF2 message): message dimension error." << endl;
		exit(-1);
	}	
	return message*this->_G_Goppa;
}


void GoppaCode::_split(GF2EX& p0, GF2EX& p1, const GF2EX p)
{
	/*split polynomial p over F into even part po and odd part p1 such that p(z) = p2 (z) + z p2 (z)*/
	int t = deg(p);
	GF2E tmp;
	for (int i = 0; i <= t; i++){
		tmp = coeff(p,i); //cout << "split: " << tmp << endl;
		for (int j = 0; j < this->_m-1; j++)//tmp_sqrt=tmp^{2^{m-1}}
			tmp = sqr(tmp);	
		
		if(i%2==0){
			SetCoeff(p0, i/2, tmp);
			
		}				
		else{
			SetCoeff(p1, i/2, tmp);	
		}
				
	}	
}

ZZ GoppaCode::_norm(GF2EX a, GF2EX b)
{
	/*This is the way in which Bernstein indicates the norm of a member of the lattice is
	to be defined->|a^2+X*b^2|*/
	GF2EX X; SetX(X);
	long degree = deg(a*a+X*b*b);
	return power(conv<ZZ>(2), degree);
}


void GoppaCode::_lattice_basis_reduce(GF2EX& M, GF2EX& N, const GF2EX s)
{
	int t = deg(this->_g);
	vector<GF2EX> a, b;
	a.push_back(conv<GF2EX>(0)); b.push_back(conv<GF2EX>(0));
	GF2EX q,r;
	DivRem(q, r, this->_g, s);
	a[0] = this->_g - q*s; b[0] = 0 - q;

	if (r == 0){
			M = s; N =1;			
			return;
		}
	/*If the norm is already small enough, we are done. 
	Otherwise, intialize the base case of the recursive prfailureocess.*/
	if (_norm(a[0],b[0]) > power(conv<ZZ>(2),t)){
		a.push_back(conv<GF2EX>(0)); b.push_back(conv<GF2EX>(0));
		DivRem(q,r,s,a[0]);
		a[1] = r; b[1] = 1 - q*b[0];
		
		if (a[1] == 0){
			M = s; N =1;			
			return;
		}
	}
	else{
		M = a[0]; N = b[0];
		return;
	}
	/*Continue subtracting integer multiples of the shorter vector from the longer until
	the produced vector has a small enough norm.*/
	int i = 1;
	while (_norm(a[i],b[i]) > power(conv<ZZ>(2),t)){
		a.push_back(conv<GF2EX>(0)); b.push_back(conv<GF2EX>(0));
		DivRem(q,r,a[i-1],a[i]);
		a[i+1] = r; b[i+1] = b[i-1] - q*b[i];
		i++;
		//cout <<"a["<<i<<"]: "<<a[i]<<endl;cout <<"b["<<i<<"]: "<<b[i]<<endl;
	}
	M = a[i]; N = b[i];	
	return;	
}

void GoppaCode::_lattice_basis_reduce_I(GF2EX& sigma, GF2EX& omega, GF2EX A, GF2EX B, const int degree)
{	
	// v*B = s mod A, deg(A) > deg(B)
	vector<GF2EX> s; vector<GF2EX> u; vector<GF2EX> v;
	s.push_back(A);	s.push_back(B);
	u.push_back(conv<GF2EX>(1));	u.push_back(conv<GF2EX>(0));
	v.push_back(conv<GF2EX>(0));	v.push_back(conv<GF2EX>(1));
	int i = 1;
	GF2EX q,r;
	while (deg(s[i]) > degree){
		i++;
		s.push_back(conv<GF2EX>(0)); u.push_back(conv<GF2EX>(0)); v.push_back(conv<GF2EX>(0));
		DivRem(q,r,s[i-2],s[i-1]);
		s[i] = s[i-2] - q*s[i-1];
		u[i] = u[i-2] - q*u[i-1];
		v[i] = v[i-2] - q*v[i-1];
		/*cout << "s" << i << s[i] << endl;
		cout << "q" << i << q << endl;
		cout <<"deg(s) " << deg(s[i]) << endl;*/
	}
	sigma = v[i];
	
	omega = s[i];
}

void GoppaCode::_xgcd(GF2EX& S, GF2EX& U, GF2EX& V, GF2EX A, GF2EX B)
{
	//Si = Ui*A + Vi*B
	vector<GF2EX> s; vector<GF2EX> u; vector<GF2EX> v;
	s.push_back(A);	s.push_back(B);
	u.push_back(conv<GF2EX>(1));	u.push_back(conv<GF2EX>(0));
	v.push_back(conv<GF2EX>(0));	v.push_back(conv<GF2EX>(1));
	int i = 1;
	GF2EX q,r;
	GF2EX X; SetX(X); 
	cout << "S[0]: " << s[0] << endl;
	cout << "S[1]: " << s[1] << endl;
	while (!IsZero(s[i])){
		i++;
		q = 0;
		int t = deg(s[i-2]) - deg(s[i-1]);
		s.push_back(conv<GF2EX>(0)); u.push_back(conv<GF2EX>(0)); v.push_back(conv<GF2EX>(0));
		DivRem(q,r,s[i-2],s[i-1]);
		s[i] = s[i-2] + q*s[i-1];
		u[i] = u[i-2] + q*u[i-1];
		v[i] = v[i-2] + q*v[i-1];
		/*GF2E coeff1, coeff2;
		GF2EX Q;
		while (t >= 0){	
			coeff1 = LeadCoeff(s[i-2]);
			coeff2 = LeadCoeff(s[i-1]);
			Q = coeff1 * inv(coeff2) * power(X, t);
			s[i] = s[i-2] + Q*s[i-1];
			u[i] = u[i-2] + Q*u[i-1];
			v[i] = v[i-2] + Q*v[i-1];
			q += Q;  
			s[i-2] = s[i];
			u[i-2] = u[i];
			v[i-2] = v[i];
			t = deg(s[i-2]) - deg(s[i-1]);	
		}*/		
		cout << "i-1: " << i-1 <<endl;
		cout << "S: " << s[i] << endl;
		cout << "q: " << q << endl;
		cout << "V before add: " << v[i-1] << endl; 
		cout << "V: " << q*v[i-1] << endl;
	}
	
		
	S = 1;
	U = u[i-1]/s[i-1];
	V = v[i-1]/s[i-1];
}

void GoppaCode::_extended_euclidean(GF2EX& sigma, GF2EX& omega, GF2EX A, GF2EX B, const int degree)
{	
	// v*B = s mod A, deg(A) > deg(B)
	vector<GF2EX> s; vector<GF2EX> u; vector<GF2EX> v;
	s.push_back(A);	s.push_back(B);
	u.push_back(conv<GF2EX>(1));	u.push_back(conv<GF2EX>(0));
	v.push_back(conv<GF2EX>(0));	v.push_back(conv<GF2EX>(1));
	int i = 1;
	GF2EX q,r;
	while (deg(s[i]) >= degree){
		i++;
		s.push_back(conv<GF2EX>(0)); u.push_back(conv<GF2EX>(0)); v.push_back(conv<GF2EX>(0));
		DivRem(q,r,s[i-2],s[i-1]);
		s[i] = s[i-2] - q*s[i-1];
		u[i] = u[i-2] - q*u[i-1];
		v[i] = v[i-2] - q*v[i-1];
		//cout << "degree sigma" << deg(v[i]) << endl;
	}
	sigma = v[i];
	
	omega = s[i];
}

void GoppaCode::_extended_euclidean_I(GF2EX& sigma, GF2EX& omega, GF2EX A, GF2EX B, const int degree) //output sigma & omega.
{	
	// v*B = s mod A, deg(A) > deg(B)
	vector<GF2EX> s; vector<GF2EX> u; vector<GF2EX> v;
	s.push_back(A);	s.push_back(B);
	u.push_back(conv<GF2EX>(1));	u.push_back(conv<GF2EX>(0));
	v.push_back(conv<GF2EX>(0));	v.push_back(conv<GF2EX>(1));
	int i = 1;
	GF2EX q;
	GF2EX X; SetX(X); 
	while (deg(s[i]) >= degree){
		i++; 
		q = 0;
		int t = deg(s[i-2]) - deg(s[i-1]);		
		s.push_back(conv<GF2EX>(0)); u.push_back(conv<GF2EX>(0)); v.push_back(conv<GF2EX>(0));
		GF2E coeff1, coeff2;
		GF2EX Q;
		//calculate s[i] and q.
		while (t >= 0){
			coeff1 = LeadCoeff(s[i-2]);
			coeff2 = LeadCoeff(s[i-1]);
			Q = coeff1 * inv(coeff2) * power(X, t);
			s[i] = s[i-2] + Q*s[i-1];
			q += Q;  
			s[i-2] = s[i];
			t = deg(s[i-2]) - deg(s[i-1]);								
		}
		u[i] = u[i-2] - q*u[i-1];
		v[i] = v[i-2] - q*v[i-1];
	}
	sigma = v[i];
	omega = s[i];
} 

void GoppaCode::_extended_euclidean_II(GF2EX& sigma, GF2EX& omega, GF2EX A, GF2EX B, const int degree)
{
	vector<GF2EX> s; vector<GF2EX> u; vector<GF2EX> v;
	s.push_back(A);	s.push_back(B);
	u.push_back(conv<GF2EX>(1));	u.push_back(conv<GF2EX>(0));
	v.push_back(conv<GF2EX>(0));	v.push_back(conv<GF2EX>(1));
	int i=1;	
	GF2E delta; GF2EX q;
	GF2EX X; SetX(X);
	while (deg(s[i]) >= degree){
		i++;
		delta = 1; q = 0;
		int t = deg(s[i-2])-deg(s[i-1]);
		s.push_back(conv<GF2EX>(0)); u.push_back(conv<GF2EX>(0)); v.push_back(conv<GF2EX>(0));
		GF2E coeff1, coeff2;
		GF2EX Q;
		
		while (t >= 0){
			coeff1 = LeadCoeff(s[i-2]);
			coeff2 = LeadCoeff(s[i-1]);
			Q = coeff1 * power(X,t);
			s[i] = coeff2*s[i-2] - Q*s[i-1];
			delta *= coeff2;				
			q = coeff2*q + Q;  
			s[i-2] = s[i];
			t = deg(s[i-2])-deg(s[i-1]);
		}
		
		u[i] = delta*u[i-2] - q*u[i-1];
		v[i] = delta*v[i-2] - q*v[i-1];
	}
	sigma = v[i];
	omega = s[i];
}

mat_GF2 GoppaCode::Decode(mat_GF2 word, string mode)
{
	
	if (word.NumRows()!=1 || word.NumCols()!= this->_SyndromeCalculator.size()){
		cout << "mat_GF2 GoppaCode::Decode(mat_GF2 word): Error! word dimension dismatches."<< endl;
		exit(-1);
	}	
	//Compute the syndrome necessary.
	GF2EX syndrome_poly = conv<GF2EX>(0);
	int tmp = word.NumCols();
	for (int i = 0; i < tmp; i++)
		syndrome_poly += this->_SyndromeCalculator[i]*word[0][i];
        cout << "Decode syndrome_poly: " << syndrome_poly << endl;
	mat_GF2 error = SyndromeDecode(syndrome_poly,mode);
	return word+error;
	
}
	
	
mat_GF2 GoppaCode::SyndromeDecode(GF2EX syndrome_poly, string mode)
{	
	//Unlike Decode(), SyndromeDecode uses the syndrome polynomial as its input, outputs the error pattern.	
	//cout << "Decode syndrome_poly: " << syndrome_poly << endl;
	if (mode == "Patterson"){
		//Compute the syndrome necessary for Patterson’s Algorithm.
		//Take the necessary square root
		GF2EX g0,g1; _split(g0,g1,this->_g);	
		GF2EX sqrt_X = (g0*InvMod(g1,this->_g))%this->_g;
		GF2EX T = InvMod(syndrome_poly,this->_g);
		cout << "INVERSE1: " << T << endl;

		_xgcd(g0, g1, T, this->_g, syndrome_poly);
		cout << "INVERSE2: " << g1*this->_g+T*syndrome_poly << endl;

		GF2EX T0, T1, X;
		SetX(X);  _split(T0,T1,T+X);
	 	GF2EX R; rem(R,T0+sqrt_X*T1,this->_g); //R = (T0+ sqrt_X*T1).mod(g);
	 	cout << "R(x): " << R << endl;
		//Perform lattice basis reduction.
		GF2EX alpha, beta;	
		//using Euclid's alg to calculate alpha, beta.
		_lattice_basis_reduce_I(beta, alpha, this->_g, R, deg(this->_g)/2);     
		//cout << "alpha" << alpha*alpha << endl;
		//cout << "beta" << beta*beta*X << endl;
		
		//Construct the error-locator polynomial.
		GF2EX sigma;
		sigma = (alpha*alpha) + (beta*beta)*X; 
		//cout << "sigma" << sigma << endl;
		sigma = sigma/LeadCoeff(sigma);
		cout << "norm sigma" << sigma << endl;
		//For every root of the error polynomial, correct the error induced at the corresponding index.
		int codelocators_len = _codelocators.size();
		int sigma_degree = deg(sigma);	
		mat_GF2 error; error.SetDims(1,codelocators_len);	

		GF2EX X_2m; 

		X_2m = X;
		for (int i = 0; i < this->_m; i++){
			X_2m = X_2m*X_2m%sigma;
		}
			
		
		if (!IsZero(X_2m-X)){
			cout << "Decode failure in Divisibility!" << endl;
			for (int i = 0; i < codelocators_len; i++)
				error[0][i] = 1;		
			return error; //Decode failure, 'error' = all ones.
		}

		else{
			for (int i = codelocators_len-1; i >= 0; i--){
				/*GF2E tmp = conv<GF2E>(0);
				for(int j = 0; j <= sigma_degree; j++)
					tmp += power(_codelocators[i], j)*coeff(sigma, j);
				cout << "tmp1 " << tmp << IsZero(tmp) << endl;*/

				GF2E tmp = coeff(sigma,sigma_degree);

				for(int j = 0; j <= sigma_degree-1; j++){
					tmp = _codelocators[i]*tmp + coeff(sigma,sigma_degree-j-1);
				}
				if (IsZero(tmp))
					error[0][i] = 1;									
			}	
			return error;
		}
	}

	else if (mode == "Euclidean"){
		GF2EX sigma, omega;
		_extended_euclidean_I(sigma,omega, this->_g,syndrome_poly,deg(this->_g)/2);
		
		sigma = sigma/LeadCoeff(sigma);
		//For every root of the error polynomial,
		//correct the error induced at the corresponding index.
		int codelocators_len = _codelocators.size();
		int sigma_degree = deg(sigma);	
		
		mat_GF2 error; error.SetDims(1,codelocators_len);
		GF2EX X; SetX(X); 
		GF2EX X_2m;
		X_2m = X;
		for (int i = 0; i < this->_m; i++)
			X_2m = X_2m*X_2m;

		if (!IsZero((X_2m-X)%sigma)){
			cout << "Decode failure!" << endl;		
			for (int i = 0; i < codelocators_len; i++)
				error[0][i] = 1;		
			return error; //Decode failure, 'error' = all ones.
	
		}

		for (int i = 0; i < codelocators_len; i++){
				GF2E tmp = conv<GF2E>(0);
				for(int j = 0; j <= sigma_degree; j++)
					tmp += power(_codelocators[i], j)*coeff(sigma, j);
				//cout << "tmp " << tmp << IsZero(tmp) << endl;
				if (IsZero(tmp)){
					error[0][i] = 1;
					//cout << "zero point" << endl;
				}				
		}	
		return error;
	}
	
}

GF2EX GoppaCode::GetGoppaPolynomial(void)
{
	return this->_g;
}

mat_GF2E GoppaCode::GetParityCheckMatrixPoly(void)
{
	return this->_H_Goppa_Poly;
}

mat_GF2 GoppaCode::GetParityCheckMatrix(void)
{
	return this->_H_Goppa;
}

mat_GF2 GoppaCode::GetGeneratorMatrix(void)
{
	return this->_G_Goppa;
}

vector<GF2E> GoppaCode::GetCodeLocators(void)
{
	return this->_codelocators;
}

vector<GF2EX> GoppaCode::GetSyndromeCalculator(void)
{
	return this->_SyndromeCalculator;
}

GF2E GoppaCode::GetGf2eGenerator()
{

	return this->_Gf2eGenerator;
}

GF2E GoppaCode::_GetGf2eGenerator(int degree)
{
	vector<int> final_factor_list;	vector<vector<int> > comb_list;
	integer_factor(final_factor_list, int (pow(double(2),double(degree)))-1);
	GF2E primitive_root;
        int factor_list_len = final_factor_list.size();
        /*for (int i = 0; i < final_factor_list.size(); i++)
                cout << final_factor_list[i] << " ";
        cout << endl;*/
	while (1){
		primitive_root = random_GF2E();
		bool output = true;
		if (IsZero(primitive_root))
			continue;
		
		for (int i = 1; i < factor_list_len; i++){
                        //cout << "i " << factor_list_len << endl;
			combination(comb_list,final_factor_list,factor_list_len,i);
                        int comb_list_nrows = comb_list.size();
                        int comb_list_ncols = comb_list[0].size();

			for (int j = 0; j < comb_list_nrows; j++){
				int exp = 1;				
				for (int k = 0; k < comb_list_ncols; k++)
					exp *= comb_list[j][k];
                                //cout << "exp= " << exp << endl;
				if (IsOne(power(primitive_root, conv<ZZ> (exp)))){
					output = false;
					break;	
				}				
				else{
					output = true;
					continue;
				}
			}
			if (output == false)
				break;
		}
		if (output == true)
			break;
	}	
	//cout << " found! primitiveroot " << primitive_root << endl;
	return primitive_root;
}

void combination(vector<vector<int> >& out_comb_list, vector<int> in_comb_list, int n, int k)
{
        if (n == 0 || k ==0){
        	cout << "combination: cannot choose zero elements from the set." << endl;
        	exit(-1);
        }


        int tmp1 = 1, tmp2 = 1;
        for (int i = n; i > n-k; i--)
                tmp1 *= i;
        for (int i = 1; i <=k; i++)
                tmp2 *= i;
        int nrows = tmp1/tmp2; int cnt = 0;
        vector<int> tmp(k);
        out_comb_list.resize(nrows);
        for(int i = 0; i < nrows; i++)
                out_comb_list[i].resize(k);
        stack<int> s;
        s.push(-1);
        while(!s.empty()){
                if (s.size() > k){
                    for(int i = 0; i < k; i++)
                        out_comb_list[cnt][i] = tmp[i];
                    s.pop();
                    cnt++;
                    continue;
                }

                int start = s.top() + 1;
                s.pop();
                for(int i = start; i < n; i++){
                    tmp[s.size()] = in_comb_list[i];
                    s.push(i);
                    s.push(i);
                    break;
                }
        }
        
}

void integer_factor(vector<int>& out_factor_list, int in_small_integer)
{/*factor small integer of type 'int'*/
	if (in_small_integer > 2){
		while (in_small_integer%2 == 0 && in_small_integer !=2){
			out_factor_list.push_back(2);
			in_small_integer /= 2;						
		}
		int val = sqrt(in_small_integer);		
		for (int i = 3; i < val; i += 2){
			while (val >= i){
				if (in_small_integer%i == 0){
					in_small_integer /= i;
					out_factor_list.push_back(i);

					val = sqrt(in_small_integer);
				}
				else
					break;					
			}			
		}
		out_factor_list.push_back(in_small_integer);

	}

	else if (in_small_integer < 0){
		cout << "GoppaCode.cpp: void integer_factor(factor_list& out_factor_list, int in_small_integer):\
                         Error, in_small_integer must be positive." << endl;
		exit(-1);	
	}	

	else{
		out_factor_list.push_back(in_small_integer);
	}
}


//This is a help function which will be useful for encryption.
int GetRowVectorWeight(mat_GF2 input_matrix)
{
	if (input_matrix.NumRows() != 1){
		cout << "GetRowVectorWeight(mat_GF2 input_marix): Error, input_matrix must be one-dimensional." << endl;
		exit(-1);
	}
	int weight = 0, ncols = input_matrix.NumCols();
	for (int i = 0; i < ncols; i++)
		if (input_matrix[0][i] == conv<GF2>(1))
			weight += 1;
	return weight;
		
}

int GetRowVectorWeight(vector<GF2> input_vec)
{
	int weight = 0; 
	int ncols = input_vec.size();
	for (int i =0; i < ncols; i++)
		if (input_vec[i] == conv<GF2>(1))
			weight += 1;
	return weight;
}

/*void PrintResult(mat_GF2 data,string addr)
{
	ofstream file(addr.c_str());
	if(file.is_open()){
		file << data;
	}
	file.close();
}*/


