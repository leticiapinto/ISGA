#include <iostream>
#include <algorithm>
#include <vector>
#include <random>
#include <ctime>
#include <cassert>

#include <fstream>
#include <string>

using namespace std;

#define CALIBRANDO

#define MAX_N 3000
int PERT_d = 1;
int GREED_p = 80;

typedef uint_fast64_t t_elapsed;      // Defined as tenth of a millisecond
static double timeMeasure = 10000.0;  // Defined as tenth of a millisecond
inline t_elapsed elapsed(bool reset = false)
{
	static std::clock_t start = std::clock();
	if ( reset ) start = std::clock();
	return (timeMeasure*double(std::clock()-start))/double(CLOCKS_PER_SEC);
}
t_elapsed tottime, bsttime;

string inst;
ofstream fout;
vector<int_fast16_t> P, W, I;
int_fast16_t C,N;
vector<double> F2;
vector<int_fast16_t> F1;
int_fast16_t maxFi;

struct t_sol {
	vector<int_fast8_t> X;
	int_fast32_t o;
	double w;
	void print(){
		fout << "Rslt " << inst  << " " << GREED_p << " " << PERT_d << " " << o << " " << tottime << " " << bsttime << endl;
		for (int i = 0; i < N; ++i)
			if (X[i]) fout << i << ", ";
		/* cerr << endl<< " *** "<<endl;
		for (int i = 0; i < N; ++i)
			if (X[i]==0) cerr << i << ", ";//*/
		fout << endl;
	};
};



void validate(t_sol S)
{
	int weight = 0;
	for (int i = 0; i < N; ++i)
	{
		if (S.X[i]) weight+=W[i];
	}
	//cerr <<"weight total  : "<< weight<<endl;
}


void readinst(string file)
{
	ifstream fin (file);

	if (fin.is_open())
	{
		string nombre;
		//fin >> nombre;
		fin >> N;
		//cerr << N << endl;

		P.resize(N*N);
		W.resize(N);
		F1.resize(N);
		F2.resize(N);

		for(auto p = P.begin(),_p = P.end(); p < _p; p+=N+1) fin >> *p;

		for(int i = 0; i < N; i++)
		{
			for (auto p = P.begin()+i*N+i+1, q = p+N-1, _p=P.begin()+(i+1)*N; p < _p; q+=N, ++p)
			{
				int y;
				fin >> y;
				*q = *p = y;
			}
		}
		fin >> C;
		fin >> C;

		for(auto &w: W) fin >> w;
		fin.close();
	}

	//cerr<<C<<endl;

	/*
	cerr<<"\nMatriz P :\n\n";
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			cerr<<P[i*N+j]<<"\t";
		}
		cerr<<"\n";
	}
	cerr<<"\n";//*/
}

// Solution representation:
// S saves indexes of the solution, but
// S[0] is the weight of the current solution S.

void Verify (t_sol S)
{
	int wS = 0;
	for (int i = 0; i < N; ++i)
	{
		if (S.X[i]) wS+=W[i];
	}
	assert(wS <= C);
	assert(wS == S.w);
	//cerr <<"weight total  : "<< wS <<endl;

	int oS = 0;
	for (int i = 0; i < N; ++i)
	{
		if (S.X[i]) oS+=P[i*N + i];
	}
	//cerr <<" sum : "<< oS <<endl;
	for (int i = 0; i < N; ++i)
	for (int j = i+1; j < N; ++j)
	{
		if (S.X[i] and S.X[j]) oS+=P[i*N + j];
	}
	//cerr <<" sum : " << oS << " " << S.o << endl;
	assert(oS == S.o);
}

void CalculateF1 (t_sol S)
{
	auto x = S.X.begin();
	auto f1 = F1.begin();
	auto w = W.begin();

	maxFi=0;
	for (int i = 0; i < N; ++i, ++x, ++f1, ++w)
	{
		if(*x or S.w + *w > C)
		{
			F1[i] = 0;
		}
		else
		{
			auto p = P.begin()+i*N;
			int_fast32_t OS2 = S.o + *(p + i);
			for (auto x = S.X.begin(), _x = S.X.end(); x<_x; ++x, ++p)
				if (*x) OS2 += (*p);
			F1[i] = OS2;
			if (F1[i]>F1[maxFi]) maxFi = i;
		}
	}
}

int_fast16_t RandomF1 ()
{
	int_fast16_t limit = GREED_p * F1[maxFi]/100.0;
	I.resize(MAX_N);
	auto i= I.begin();
	for (auto it = F1.begin(), _it = F1.end(); it<_it; ++it)
		if (*it >= limit) *(i++) = it - F1.begin();
	return *( I.begin() + rand()%(i-I.begin()) );
}

void CalculateF2 (t_sol S)
{
	auto x = S.X.begin();
	auto f2 = F2.begin();
	auto w = W.begin();

	maxFi=0;
	for (int i = 0; i < N; ++i, ++x, ++f2, ++w)
	{
		if(*x or S.w + *w > C)
		{
			F2[i] = 0.0;
		}
		else
		{
			auto p = P.begin()+i*N;
			int_fast32_t OS2 = S.o + *(p + i);

			for (auto x = S.X.begin(), _x = S.X.end(); x<_x; ++x, ++p)
				if (*x) OS2 += (*p);

			F2[i] = double(OS2)/double(S.w + W[i]);
			if (F2[i]>F2[maxFi]) maxFi = i;
		}
	}
}

int_fast16_t RandomF2 ()
{
	double limit = GREED_p * F2[maxFi]/100.0;
	I.resize(MAX_N);
	auto i= I.begin();
	for (auto it = F2.begin(), _it = F2.end(); it<_it; ++it)
		if (*it >= limit) *(i++) = it - F2.begin();
	return *( I.begin() + rand()%(i-I.begin()) );
}

void CreateGreedy (t_sol & S)
{
	S.X.assign(N,0);
	S.o = 0;
	S.w = 0;
	Verify(S);
	CalculateF2(S);
	while (F2[maxFi]>0.0)
	{
		S.X[maxFi] = 1;
		auto p = P.begin()+maxFi*N;
		int_fast32_t OS2 = S.o;// + *(p + maxFi);
		for (auto x = S.X.begin(), _x = S.X.end(); x<_x; ++x, ++p)
				if (*x) OS2 += (*p);
		S.o = OS2;
		S.w += W[maxFi];
		Verify(S);
		CalculateF2(S);
	}
}

void LocalSearch (t_sol & S)
{
//	cerr <<"LS 1 "<< S.o << endl;
	LocalSearchReset:
	auto x0 = S.X.begin();
//	cerr <<"LS 2"<<endl;
	long aaa = S.X.size(); aaa*=aaa;
	for (;x0<S.X.end();x0++)
	{
//		cerr <<"LS 3"<<endl;
		while (*x0) ++x0;
//		cerr <<"LS 4"<<endl;
		if (x0>=S.X.end()) break;
		auto i0 = x0 - S.X.begin();
//		cerr<<", w[i0]"<< W[i0]<<"i0 "<<i0<<endl;
		auto p = P.begin()+i0*N;
		auto w0 = S.w + W[i0];
		auto o0 = *(p+i0);
		for (auto x = S.X.begin(), _x = S.X.end(); x<_x; ++x, ++p)
				if (*x) o0 += (*p);
//		cerr <<"LS 5"<<endl;
		auto x1 = S.X.begin();

		for (;x1<S.X.end();++x1)
		{
			while (*x1==0) ++x1;
//			cerr <<"LS 6"<<endl;
			if(x1>=S.X.end()) break;
			auto i1 = x1-S.X.begin();
//			cerr<<"w0"<<w0<<", w[i1]"<< W[i1]<<"i1 "<<i1<<endl;
//			cerr <<"LS 7"<<endl;
			//if (!--aaa) { cerr << "dfkljadsja dñls"; exit;}
			if ( w0 - W[i1] < C )
			{
//				cerr <<"LS 8"<<endl;
				auto p = P.begin()+i1*N;
				auto o1 = *(p+i0);//equipara al reemplazo.
				for (auto x = S.X.begin(), _x = S.X.end(); x<_x; ++x, ++p)
					if (*x) o1 += (*p);
//				cerr <<"LS 9"<<endl;
				if (o0 > o1)
				{
//					cerr <<"LS 10"<<endl;
					*x0 = 1;
					*x1 = 0;
					S.w = S.w + W[i0] - W[i1];
					S.o = S.o + o0 - o1;
					Verify(S);
					//cerr << " Mejora " << o0 << "  " << o1 << " " << S.o << endl;
					goto LocalSearchReset;
				}
			}
		}
	}
	//cerr << " LS paso " << S.o << endl;
}

void IteratedGreedy ()
{
	t_sol bS, S;

	elapsed(true);
	CreateGreedy (S);

	LocalSearch (S);

	bS = S;
	//S.print();

	int nIts = 4*N;
	while (--nIts)
	{
		// randomly sorts the elements in X
		{
			I.resize(MAX_N);
			auto i= I.begin();
			for(auto it = S.X.begin(); it < S.X.end(); ++it )
				if (*it) *(i++) = it - S.X.begin();
			I.resize(i-I.begin());
			shuffle(I.begin(),I.end(), mt19937{random_device{}()});
		}

		// removes the PERT_d with their weight and benefit
		for (auto i = I.begin() + PERT_d - 1; i >= I.begin(); --i)
		{
			auto p = P.begin()+(*i)*N;
			int_fast32_t OS2 = S.o;
			for (auto x = S.X.begin(), _x = S.X.end(); x<_x; ++x, ++p)
					if (*x) OS2 -= (*p);
			S.o = OS2;
			S.w -= W[*i];
			S.X[*i] = 0;
			Verify(S);
		}
		// rebuilds solution
		CalculateF1(S);
		while (F1[maxFi]>0)
		{
			maxFi = RandomF1 ();
			S.X[maxFi] = 1;
			auto p = P.begin()+maxFi*N;
			int_fast32_t OS2 = S.o;
			for (auto x = S.X.begin(), _x = S.X.end(); x<_x; ++x, ++p)
					if (*x) OS2 += (*p);
			S.o = OS2;
			S.w += W[maxFi];
			Verify(S);
			CalculateF1(S);
		}
		//cerr<<"6"<<endl;
		LocalSearch (S);
		//cerr<<"7"<<endl;
		if (S.o > bS.o)
		{
			bS = S;
			//cerr <<" best : " << bS.o << endl;
			bsttime = elapsed();
		}
		else
		{
			S = bS;
		}
	}
	tottime = elapsed();
	bS.print();
}



main(int argc, char* argv[])
{
	string file = argv[1];
	cout<<endl<<"file : "<<file<<endl;
	cout<< "PERT_d ="<< PERT_d<< ",GREED_p = "<<GREED_p<<endl;
	W.reserve(MAX_N);
	P.reserve(MAX_N*MAX_N);
	
	//readinst("../1000_25/1000_25_3.dat");
	//readinst("../1000_50/1000_50_5.dat");
	//readinst("../1000_75/1000_75_4.dat");
	//readinst("../1000_100/1000_100_6.dat");
	readinst(file);


	//readinst("r_10_100_13.txt");

	IteratedGreedy();

}
