#include <gtest/gtest.h>

// A Unit Test deliberately designed to fail
class testFail : public ::testing::Test { };

TEST_F(testFail, thisShallFail)
{
   ASSERT_EQ(0, 1);
}

int main(int ac, char *av[])
{
   ::testing::InitGoogleTest(&ac, av);
   return RUN_ALL_TESTS();
}
