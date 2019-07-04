//NOTE: the 2*N should be in power of 2

//IMPORTANT HEADERS AND MICROS
//------------------------------------------------------------------------------------------------------------------------
#include<bits/stdc++.h>
using namespace std;
#define ll long long int
#define fast ios_base::sync_with_stdio(NULL),cin.tie(NULL),cout.tie(NULL);
#define endl "\n"
const long double PI = 3.141592653589793238463;

//PROGRAM
//-----------------------------------------------------------------------------------------------------------------------
complex<double> find_w(int n)
{
	//formula: e^(i*Theta) = cos(Theta) + i*sin(Theta)
	return complex<double>(cos(2*PI/n),sin(2*PI/n));
}

complex<double> find_w_inverse(int n)
{
	//formula: e^(-i*Theta) = cos(-Theta) + i*sin(-Theta)
	return complex<double>(cos(-(2*PI/n)), sin(-(2*PI/n)));
}

vector<complex<double>> FFT(vector<complex<double>> A,int m,complex<double> w)
{
	if(m <= 1)
		return A;

	vector<complex<double>> A_odd;
	vector<complex<double>> A_even;

	for(int i = 0; i < m; i+=2)
	{
		A_even.push_back(A[i]);
		A_odd.push_back(A[i+1]);
	}

	vector<complex<double>> F_even = FFT(A_even, m/2, w*w);
	vector<complex<double>> F_odd = FFT(A_odd, m/2, w*w);

	vector<complex<double>> F(m);
	complex<double> temp = 1;
	for(int i = 0; i < m/2; i++)
	{
		F[i] = F_even[i] + temp*(F_odd[i]);
		F[i+m/2] = F_even[i] - temp*(F_odd[i]);
		temp = temp*w;
	}

	return F;
}

int main()
{
	// freopen("in.txt","r",stdin);
	// freopen("out.txt","w",stdout);

	fast

	int n;
	cin >> n; //input no of elements

	n = 2*n; 
	vector<complex<double>> a;
	vector<complex<double>> b;

	//init with 0
	for(int i = 0; i < n; i++)
	{
		a.push_back(0);
		b.push_back(0);
	}

	//input the polynomials A and B
	for(int i = 0; i < n/2; i++)	
		cin >> a[i];
	for(int i = 0; i < n/2; i++)
		cin >> b[i];

	//finding Omega and inverse of Omega 
	complex<double> w = find_w(n);
	complex<double> w_inverse = find_w_inverse(n);

	//finding FFT of equations A and B separately
	vector<complex<double>> F_A = FFT(a,n,w);
	vector<complex<double>> F_B = FFT(b,n,w);

	//combing the FFT of A and B
	vector<complex<double>> F_C(n);
	for(int i = 0; i < n; i++)
		F_C[i] = F_A[i] * F_B[i];

	//finding the inverse DFT and dividing it by 1/n
	vector<complex<double>> ans = FFT(F_C, n, w_inverse);
	double fac = 1.0/n;
	for(int i = 0; i < n; i++)
		ans[i] = ans[i]*fac;

	//output the real part of the result
	for(int i = 0; i < n; i++)
		cout << real(ans[i]) << ' ';
}
