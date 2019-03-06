#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>

#include <cmath>


struct T
{
	T()
	= default;

	T(double v):m_element(v)
	{
	}

	double operator()() const
	{
		return  m_element;
	}

private:
	double m_element{};
};

T max(T A, T B)
{
	return std::max(A(), B());
}

T make_identity_add()
{
	return {- std::numeric_limits<double>::infinity()};
}

T make_identity_mul()
{
	return {};
}

T operator+(T A, T B)
{
	return max(A, B);
}

T operator*(T A, T B)
{
	return {A() + B()};
}

T pow(T A, unsigned long long p)
{
	T res = make_identity_mul();

	while(p>0)
	{
		--p;
		res = res * A;
	}

	return res;
}

std::ostream &operator<<( std::ostream &output, T A )
{
	output << A();
	return output;
}


struct Monome
{
	Monome()
	{
		m_powCoeff = std::unique_ptr<unsigned int[]>(new unsigned int[m_dimVar]);
	}
	Monome(std::initializer_list<unsigned int> l)
	{
		m_dimVar = l.size();

		if(m_dimVar == 0)
		{
			m_dimVar = 1;
			m_powCoeff = std::unique_ptr<unsigned int[]>(new unsigned int[m_dimVar]);
			m_powCoeff[0] = 0;
			return;
		}
		m_powCoeff = std::unique_ptr<unsigned int[]>(new unsigned int[m_dimVar]);
		unsigned long long i = 0;
		for (auto v:l)
		{
			m_powCoeff[i++] = v; //gcc-trunk -O1 -O2 same than (++i - 1)
		}
	}

	T operator()(std::initializer_list<T> l)
	{
		//compute product of each var with their power
		T res = make_identity_mul();
		unsigned long long i = 0;
		for (auto v:l)
		{
			if(i<m_dimVar)
			{
				res = res * pow(v, m_powCoeff[i++]);
			}
			else
			{
				return res;
			}

		}

		return res;
	}
private:
	unsigned long long m_dimVar{0};
	std::unique_ptr<unsigned int[]> m_powCoeff;

	friend std::ostream &operator<<( std::ostream &output, const Monome& M );
};

std::ostream &operator<<( std::ostream &output, const Monome& M )
{
	output << "[ ";
	for(unsigned long long i=0; i<M.m_dimVar-1; ++i)
		output << M.m_powCoeff[i] << ", ";
	output << M.m_powCoeff[M.m_dimVar-1];
	output << " ]";
	return output;
}


struct Polynome
{
	Polynome()
	= default;
	Polynome(std::initializer_list<std::pair<T, std::initializer_list<unsigned int> > > l)
	{
		m_nbMonome = l.size();
		m_coeffMonome = std::unique_ptr<std::pair<T, Monome>[]>(new std::pair<T, Monome>[m_nbMonome]);

		unsigned long long i = 0;
		for (auto m:l)
		{
			m_coeffMonome[i++] = m;
		}
	}

	T operator()(std::initializer_list<T> l) const
	{
		T res = make_identity_add();
		for(unsigned long long i=0; i<m_nbMonome; ++i)
		{
			res = res + m_coeffMonome[i].first*m_coeffMonome[i].second(l);
		}
		return res;
	}
private:
	unsigned long long m_nbMonome{0};
	std::unique_ptr<std::pair<T, Monome>[]> m_coeffMonome{};

	friend std::ostream &operator<<( std::ostream &output, const Polynome& P );
};

std::ostream &operator<<( std::ostream &output, const Polynome& P )
{


//	std::cout << P.m_nbMonome;
	for(unsigned long long i=0; i<P.m_nbMonome-1; ++i)
		output << P.m_coeffMonome[i].first;// << "*" << P.m_coeffMonome[i].second << " + ";
	output << P.m_coeffMonome[P.m_nbMonome - 1].first << "*" << P.m_coeffMonome[P.m_nbMonome - 1].second;
	return output;
}

struct Domain_1D
{
	Domain_1D()= default;

	Domain_1D(double inf, double sup, int N): m_inf(inf), m_sup(sup), m_N(N)
	{
		if(m_N<1)
			m_N = 1;
		if(m_N == 1)
			m_sup = m_inf; // Domain == {m_inf}

		m_point = std::unique_ptr<T[]>(new T[m_N]);

		m_point[0] = m_inf;
		m_point[m_N-1] = m_sup;

		double h = (m_sup - m_inf)/(m_N - 1);
		for (int i=1; i<m_N-1; ++i)
		{
			m_point[i] = m_inf + i*h;
		}

	}

	void operator()(const Polynome& P)
	{
		for(int i=0; i<m_N; ++i)
		{
			m_point[i] = P({m_point[i]});
		}
	}

	void saveToFile(const std::string& filename)
	{
		std::ofstream myfile;
		myfile.open(filename, std::ios::trunc);

		double x = m_inf;
		double h = (m_sup - m_inf)/(m_N - 1);

		for(int i= 0; i<m_N; ++i)
		{
			myfile << x << " " << m_point[i] << '\n';

			x += h;
		}

		myfile.close();
	}

private:
	double m_inf{-1.};
	double m_sup{1.};
	int m_N{3};

	std::unique_ptr<T[]> m_point{};

	friend std::ostream &operator<<( std::ostream &output, const Domain_1D& D );
};

std::ostream &operator<<( std::ostream &output, const Domain_1D& D )
{
	output << "[ ";
	for(int i=0; i<D.m_N-1; ++i)
		output << D.m_point[i] << ", ";
	output << D.m_point[D.m_N-1];
	output << " ]";
	return output;
}



int main()
{
	double h = .05;
	double binf = -50.;
	double bsup = -binf;

	Polynome P({{-1,{}}, {3,{1}}, {2,{2}}, {1,{3}} });
	Polynome Q({{-1,{}}, {2,{1}}, {-2,{2}}, {1,{3}} });

	Polynome T({{1,{}}, {3,{1}}, {1,{2}} });

	Domain_1D D_P(binf, bsup, static_cast<int>((bsup - binf)/h));
	Domain_1D D_Q(binf, bsup, static_cast<int>((bsup - binf)/h));

	D_P(P);
	D_Q(Q);

	D_P.saveToFile(R"(C:\Users\micro\AnacondaProjects\P.dat)");
	D_Q.saveToFile(R"(C:\Users\micro\AnacondaProjects\Q.dat)");

	return 0;
}
