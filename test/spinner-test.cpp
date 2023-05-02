#include <cif++/utilities.hpp>

#include <random>
#include <thread>

void test_one()
{
	std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(100, 1000);

	cif::progress_bar pb(10, "test");

	for (int i = 0; i < 10; ++i)
	{
		std::this_thread::sleep_for(std::chrono::milliseconds(distrib(gen)));

		pb.message("step " + std::to_string(i));
		pb.consumed(1);
	}
}

void test_two()
{
	cif::progress_bar pb(10, "test");


	for (int i = 0; i < 5; ++i)
		pb.consumed(1);
}

void test_three()
{
	using namespace std::literals;

	cif::progress_bar pb(10, "test");
	pb.consumed(10);

	std::this_thread::sleep_for(100ms);
}

int main()
{
	test_one();
	test_two();
	test_three();

	return 0;
}