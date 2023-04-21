#include <cif++/utilities.hpp>

#include <random>
#include <thread>

int main()
{
	cif::progress_bar pb(10, "test");

	std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(100, 1000);
 
	for (int i = 0; i < 10; ++i)
	{
		std::this_thread::sleep_for(std::chrono::milliseconds(distrib(gen)));

		pb.message("step " + std::to_string(i));
		pb.consumed(1);
	}

	return 0;
}