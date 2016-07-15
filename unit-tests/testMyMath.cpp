#include <iostream>
#include <string>
#include <list>

#include <gtest/gtest.h>
#include <blond/math_functions.h>
#include <blond/utilities.h>
#include <blond/constants.h>
#include <blond/configuration.h>
#include <blond/globals.h>


TEST(testLinspace, test1)
{
   ftype a[10];
   mymath::linspace(a, 0, -10, 10);
   ASSERT_DOUBLE_EQ(a[9], -10);
   ASSERT_DOUBLE_EQ(a[0], 0);
}

TEST(testLinspace, test2)
{
   ftype a[10];
   mymath::linspace(a, 1, 1, 10);
   ASSERT_DOUBLE_EQ(a[0], 1);
   ASSERT_DOUBLE_EQ(a[9], 1);
}

TEST(testMean, test1)
{
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

TEST(testSTD, test1)
{
   ftype a[5] = {1, 1, 1, 1, 1};
   ftype b = mymath::standard_deviation(a, 5);
   ASSERT_DOUBLE_EQ(b, 0);

   ftype c[5] = {10, 12.4, 131, -22.3, 45.6};
   ftype m = mymath::mean(c, 5);
   ftype std1 = mymath::standard_deviation(c, 5, m);
   ftype std2 = mymath::standard_deviation(c, 5);
   ASSERT_DOUBLE_EQ(std1, std2);
}


TEST(testMin, double_type)
{
   ftype b[10] = {0.8, 9.14 , 2.1, 2, 10.25, 1.3, 2.12, 3, 4, 5};
   int i = mymath::min<ftype>(b, 10);
   ASSERT_DOUBLE_EQ(b[i], 0.8);

   ftype a[5] = {10.1, 1.2, -1.3, 1.4, 1.5};
   i = mymath::min<ftype>(a, 5);
   ASSERT_DOUBLE_EQ(a[i], -1.3);
}

TEST(testMin, int_type)
{
   int b[10] = {1, 9 , 2, 2, 10, 0, 2, 3, 4, 5};
   int i = mymath::min<int>(b, 10);
   ASSERT_EQ(b[i], 0);

   int a[5] = {10, 1, -1, 1, 1};
   i = mymath::min<int>(a, 5);
   ASSERT_EQ(a[i], -1);
}

TEST(testMax, double_type)
{
   ftype b[10] = {0.8, 9.14 , 2.1, 2, 10.25, 1.3, 2.12, 3, 4, 5};
   int i = mymath::max<ftype>(b, 10);
   ASSERT_DOUBLE_EQ(b[i], 10.25);

   ftype a[5] = {10.1, 1.2, -1.3, 1.4, 1.5};
   i = mymath::max<ftype>(a, 5);
   ASSERT_DOUBLE_EQ(a[i], 10.1);
}

TEST(testMax, int_type)
{
   int b[10] = {1, 9 , 2, 2, 10, 0, 2, 3, 4, 5};
   int i = mymath::max<int>(b, 10);
   ASSERT_EQ(b[i], 10);

   int a[5] = {10, 1, -1, 100, 1};
   i = mymath::max<int>(a, 5);
   ASSERT_EQ(a[i], 100);
}

TEST(testTrapezoid, test1)
{
   ftype b[5] = {1, 2, 3, 4, 5};
   ftype trap = mymath::trapezoid<ftype>(b, 1, 5);
   ASSERT_DOUBLE_EQ(trap, 12);

   ftype a[5] = {1.1, 1.2, 1.3, 1.4, 1.5};
   trap = mymath::trapezoid<ftype>(a, 0.1, 5);
   ASSERT_DOUBLE_EQ(trap, 0.52);

   trap = mymath::trapezoid<ftype>(a, 1, 5);
   ASSERT_DOUBLE_EQ(trap, 5.2);

   ftype c[10] = { -0.61, -0.51, 0.39, -0.54,
                   0.67, 1.4, 1.1, 1.4, 0.16, 0.9
                 };
   trap = mymath::trapezoid<ftype>(c, 1, 10);
   ASSERT_NEAR(trap, 4.215, 1e-8);
}

TEST(testTrapezoid, test2)
{
   ftype a[5] = {1.1, 1.2, 1.3, 1.4, 1.5};
   ftype trap = mymath::trapezoid<ftype>(a, a, 5);
   ASSERT_NEAR(trap, 0.52, 1e-8);

   ftype c[10] = { -0.61, -0.51, 0.39, -0.54,
                   0.67, 1.4, 1.1, 1.4, 0.16, 0.9
                 };
   trap = mymath::trapezoid<ftype>(c, c, 10);
   ASSERT_NEAR(trap, 0.21895, 1e-8);
}

TEST(testCumTrap, test1)
{
   std::string params = "../unit-tests/references/MyMath/CumTrap/";
   ftype epsilon = 1e-8;

   ftype a[10] =  { -0.61, -0.51, 0.39, -0.54,
                    0.67, 1.4, 1.1, 1.4, 0.16, 0.9
                  };
   ftype *trap = mymath::cum_trapezoid<ftype>(a, 1, 10);

   std::vector<ftype> v;
   util::read_vector_from_file(v, params + "cumtrap_reference.txt");

   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = trap[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
   }

   delete[] trap;
}



TEST(testConvolution, test1)
{
   std::vector<ftype> c, a, b, v;
   a.resize(10);
   b.resize(20);
   mymath::linspace(a.data(), 0.f, 100.f, a.size());
   mymath::linspace(b.data(), 100.f, 1000.f, b.size());
   c.resize(a.size() + b.size() - 1);
   mymath::convolution(a.data(), a.size(),
                       b.data(), b.size(), c.data());


   std::string params = "../unit-tests/references/MyMath/convolution/";
   ftype epsilon = 1e-8;

   util::read_vector_from_file(v, params + "convolution1.txt");

   ASSERT_EQ(v.size(), c.size());

   epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = c[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of convolution1 failed on i "
            << i << std::endl;
   }

   v.clear();


}

TEST(testConvolution, test2)
{
   std::vector<ftype> c, a, b, v;
   a.resize(100, 0);
   b.resize(20);
   mymath::linspace(b.data(), 100.f, 1000.f, b.size());
   c.resize(a.size() + b.size() - 1);
   mymath::convolution(a.data(), a.size(),
                       b.data(), b.size(), c.data());


   std::string params = "../unit-tests/references/MyMath/convolution/";
   ftype epsilon = 1e-8;

   util::read_vector_from_file(v, params + "convolution2.txt");

   ASSERT_EQ(v.size(), c.size());

   epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = c[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of convolution2 failed on i "
            << i << std::endl;
   }

   v.clear();


}

TEST(testConvolution, test3)
{
   std::vector<ftype> c, a, b, v;
   a = {0, 0, 0, 0, 1, 1, 1, 0, 0, 0};
   b = {1, 1, 1, 1, 0, 0, 0, 1, 1, 1};
   c.resize(a.size() + b.size() - 1);
   mymath::convolution(a.data(), a.size(),
                       b.data(), b.size(), c.data());


   std::string params = "../unit-tests/references/MyMath/convolution/";
   ftype epsilon = 1e-8;

   util::read_vector_from_file(v, params + "convolution3.txt");

   ASSERT_EQ(v.size(), c.size());

   epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = c[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of convolution3 failed on i "
            << i << std::endl;
   }

   v.clear();


}

TEST(arange, test1)
{
   std::string params = "../unit-tests/references/MyMath/arange/";
   f_vector_t v;

   auto a = mymath::arange<ftype>(0, 100);

   ftype epsilon = 1e-8;

   util::read_vector_from_file(v, params + "arange1.txt");

   ASSERT_EQ(v.size(), a.size());

   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = a[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of a failed on i "
            << i << std::endl;
   }

}


TEST(arange, test2)
{
   std::string params = "../unit-tests/references/MyMath/arange/";
   f_vector_t v;
   auto a = mymath::arange<ftype>(0, 100, 2.5);

   ftype epsilon = 1e-8;

   util::read_vector_from_file(v, params + "arange2.txt");

   ASSERT_EQ(v.size(), a.size());

   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = a[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of a failed on i "
            << i << std::endl;
   }

}


TEST(arange, test3)
{
   std::string params = "../unit-tests/references/MyMath/arange/";
   f_vector_t v;

   auto a = mymath::arange<int>(10, 80, 2);
   util::read_vector_from_file(v, params + "arange3.txt");
   //util::dump(v, "integer vector");

   ASSERT_EQ(v.size(), a.size());

   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = a[i];
      ASSERT_EQ(ref, real)
            << "Testing of a failed on i "
            << i << std::endl;
   }

}




int main(int ac, char *av[])
{
   ::testing::InitGoogleTest(&ac, av);
   return RUN_ALL_TESTS();
}
