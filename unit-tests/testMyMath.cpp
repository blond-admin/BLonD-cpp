#include <iostream>
#include <string>
#include <list>

#include <unistd.h>

#include <gtest/gtest.h>
#include "math_functions.h"
#include "utilities.h"
#include "constants.h"


TEST(testLinspace, test1) {
	ftype a[10];
	mymath::linspace(a, 0, -10, 10);
	ASSERT_DOUBLE_EQ(a[9], -10);
	ASSERT_DOUBLE_EQ(a[0], 0);
}

TEST(testLinspace, test2) {
	ftype a[10];
	mymath::linspace(a, 1, 1, 10);
	ASSERT_DOUBLE_EQ(a[0], 1);
	ASSERT_DOUBLE_EQ(a[9], 1);
}

TEST(testMean, test1) {
	ftype a[5] = { -2, -1, 0, 1, 2};
	ftype m = mymath::mean(a, 5);
	ASSERT_DOUBLE_EQ(m, 0);

	ftype b[5] = {1, 1, 1, 1, 1};
	m = mymath::mean(b, 5);
	ASSERT_DOUBLE_EQ(m, 1);

	ftype c[5] = {10, 12, 14, 16 , 18};
	m = mymath::mean(c, 5);
	ASSERT_DOUBLE_EQ(m, 14);

}

TEST(testSTD, test1) {
	ftype a[5] = {1, 1, 1, 1, 1};
	ftype b = mymath::standard_deviation(a, 5);
	ASSERT_DOUBLE_EQ(b, 0);

	ftype c[5] = {10, 12.4, 131, -22.3, 45.6};
	ftype m = mymath::mean(c, 5);
	ftype std1 = mymath::standard_deviation(c, 5, m);
	ftype std2 = mymath::standard_deviation(c, 5);
	ASSERT_DOUBLE_EQ(std1, std2);
}

TEST(testTrapezoid, test1) {
	ftype b[5] = {1, 2, 3, 4, 5};
	ftype trap = mymath::trapezoid(b, 1, 5);
	ASSERT_DOUBLE_EQ(trap, 12);

	ftype a[5] = {1.1, 1.2, 1.3, 1.4, 1.5};
	trap = mymath::trapezoid(a, 0.1, 5);
	ASSERT_DOUBLE_EQ(trap, 0.52);

	trap = mymath::trapezoid(a, 1, 5);
	ASSERT_DOUBLE_EQ(trap, 5.2);

	ftype c[10] = { -0.61, -0.51, 0.39, -0.54,
		0.67, 1.4, 1.1, 1.4, 0.16, 0.9};
	trap = mymath::trapezoid(c, 1, 10);
	ASSERT_NEAR(trap, 4.215, 1e-8);
}

TEST(testTrapezoid, test2) {
	ftype a[5] = {1.1, 1.2, 1.3, 1.4, 1.5};
	ftype trap = mymath::trapezoid(a, a, 5);
	ASSERT_NEAR(trap, 0.52, 1e-8);

	ftype c[10] = { -0.61, -0.51, 0.39, -0.54,
		0.67, 1.4, 1.1, 1.4, 0.16, 0.9};
	trap = mymath::trapezoid(c, c, 10);
	ASSERT_NEAR(trap, 0.21895, 1e-8);
}

TEST(testCumTrap, test1){
	std::string params = "../unit-tests/references/MyMath/CumTrap/";
	ftype epsilon = 1e-8;
	
	ftype a[10] =  { -0.61, -0.51, 0.39, -0.54,
		0.67, 1.4, 1.1, 1.4, 0.16, 0.9 };
	ftype *trap = mymath::cum_trapezoid(a, 1, 10);
	//dump(trap, 10, "trap\n");
	/*
	std::string result = util::exec(
		"python -c 'import scipy.integrate as ct;\
		 print ct.cumtrapz([-0.61, -0.51, 0.39, -0.54, 0.67, 1.4,\
		  1.1, 1.4, 0.16, 0.9], initial=0)'");
	std::replace(result.begin(), result.end(),'[', ' ');
	std::replace(result.begin(), result.end(),']', ' ');
	*/
	std::vector<ftype> v;
	util::read_vector_from_file(v, params + "cumtrap_reference.txt");
	//result.erase(0,1);
	//result.erase(result.end()-1,1);
	
	for (unsigned int i = 0; i < v.size(); ++i)
	{
		ftype ref = v[i];
		ftype real = trap[i];
		ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
	}
	/*
	std::istringstream ss(result);
	for (int i = 0; i < 10; ++i)
	{
		ftype b;
		ss >> b;
		ASSERT_NEAR(trap[i], b, 1e-8);
	}
	*/
}

/*
TEST(testSin, test1){
	ftype a = sin(0);
	std::string result = util::exec(
		"python -c 'import numpy as np;\
		 print np.sin(0)'");
	ftype num = std::stod(result);
	ASSERT_NEAR(a, num, 1e-8);

	a = sin(constant::pi / 2);
	result = util::exec(
		"python -c 'import numpy as np;\
		 print np.sin(np.pi/2)'");
	num = std::stod(result);
	ASSERT_NEAR(a, num, 1e-8);
}

TEST(testCos, test1){
	ftype a = cos(0);
	std::string result = util::exec(
		"python -c 'import numpy as np;\
		 print np.cos(0)'");
	ftype num = std::stod(result);
	ASSERT_NEAR(a, num, 1e-8);

	a = cos(constant::pi / 2);
	result = util::exec(
		"python -c 'import numpy as np;\
		 print np.cos(np.pi/2)'");
	num = std::stod(result);
	ASSERT_NEAR(a, num, 1e-8);

}
*/

int main(int ac, char* av[]) {
	::testing::InitGoogleTest(&ac, av);
	return RUN_ALL_TESTS();
}
