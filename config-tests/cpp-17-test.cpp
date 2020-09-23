#include <algorithm>

template<typename COMP>
class foo
{
  public:
	foo(int a, COMP&& b)
		: m_a(a), m_b(std::move(b)) {}
	
	int m_a;
	COMP m_b;
};

void bar(const int& b)
{
	int c = 1;
	auto f = new foo(c, [tag = c, b](const int& x)
		{ x < b; });
}
