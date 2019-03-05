#include <algorithm>
#include <iostream>
#include <memory>

#include <cmath>


struct T
{
	T()
	{
		m_element = 0.;
	}

	T(double v)
	{
		m_element = v;
	}

	double operator()() const
	{
		return  m_element;
	}

private:
	double m_element;
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
				res = res * pow(v, i++);
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
};

struct Polynome
{
	Polynome()
	{

	}
	Polynome(std::initializer_list<std::initializer_list<unsigned int> > l)
	{
		m_nbMonome = l.size();
		m_monome = std::unique_ptr<Monome[]>( new Monome[m_nbMonome]);
		unsigned long long i = 0;
		for (auto v:l)
		{
			m_monome[i++] = Monome(v);
		}
	}

private:
	unsigned long long m_nbMonome{0};
	std::unique_ptr<Monome[]> m_monome{};
};

int main()
{
	Monome M ({1,2,3,4,5});

	std::cout << "Dans T, 1 + -inf =";
	std::cout << T(1) + make_identity_add();
	std::cout << "\n";

	std::cout << "Dans T, 1 * 3 =";
	std::cout << T(1) * make_identity_add();
	std::cout << "\n";
	return 0;
}
