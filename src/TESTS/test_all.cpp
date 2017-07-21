#include <gtest/gtest.h>
#include <iostream>

int main(int argc, char **argv) {


	::testing::InitGoogleTest(&argc, argv);

	std::cout << 3333;
	return RUN_ALL_TESTS();
}
