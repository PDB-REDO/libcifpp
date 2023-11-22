#include "cif++/utilities.hpp"

#include <catch2/catch_test_macros.hpp>

#include <random>
#include <thread>

TEST_CASE("test_one")
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

TEST_CASE("test_two")
{
	cif::progress_bar pb(10, "test");


	for (int i = 0; i < 5; ++i)
		pb.consumed(1);
}

TEST_CASE("test_three")
{
	using namespace std::literals;

	cif::progress_bar pb(10, "test");
	pb.consumed(10);

	std::this_thread::sleep_for(100ms);
}
