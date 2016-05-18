#include <iostream>
#include <string>
#include <list>

#include <unistd.h>

#include <gtest/gtest.h>
#include "math_functions.h"
#include "utilities.h"
#include "constants.h"
#include "configuration.h"


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

TEST(testTrapezoid, test1)
{
   ftype b[5] = {1, 2, 3, 4, 5};
   ftype trap = mymath::trapezoid(b, 1, 5);
   ASSERT_DOUBLE_EQ(trap, 12);

   ftype a[5] = {1.1, 1.2, 1.3, 1.4, 1.5};
   trap = mymath::trapezoid(a, 0.1, 5);
   ASSERT_DOUBLE_EQ(trap, 0.52);

   trap = mymath::trapezoid(a, 1, 5);
   ASSERT_DOUBLE_EQ(trap, 5.2);

   ftype c[10] = { -0.61, -0.51, 0.39, -0.54,
                   0.67, 1.4, 1.1, 1.4, 0.16, 0.9
                 };
   trap = mymath::trapezoid(c, 1, 10);
   ASSERT_NEAR(trap, 4.215, 1e-8);
}

TEST(testTrapezoid, test2)
{
   ftype a[5] = {1.1, 1.2, 1.3, 1.4, 1.5};
   ftype trap = mymath::trapezoid(a, a, 5);
   ASSERT_NEAR(trap, 0.52, 1e-8);

   ftype c[10] = { -0.61, -0.51, 0.39, -0.54,
                   0.67, 1.4, 1.1, 1.4, 0.16, 0.9
                 };
   trap = mymath::trapezoid(c, c, 10);
   ASSERT_NEAR(trap, 0.21895, 1e-8);
}

TEST(testCumTrap, test1)
{
   std::string params = "../unit-tests/references/MyMath/CumTrap/";
   ftype epsilon = 1e-8;

   ftype a[10] =  { -0.61, -0.51, 0.39, -0.54,
                    0.67, 1.4, 1.1, 1.4, 0.16, 0.9
                  };
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

   for (unsigned int i = 0; i < v.size(); ++i) {
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



TEST(testRFFT, rfft1)
{
   std::string params = "../unit-tests/references/MyMath/fft/";
   ftype epsilon = 1e-8;

   std::vector<ftype> v, in;
   std::vector<complex_t> out;
   in.resize(256);
   mymath::linspace(in.data(), 0.f, 100.f, in.size());
   mymath::rfft(in, 512, out);

   util::read_vector_from_file(v, params + "rfft1_real.txt");

   ASSERT_EQ(v.size(), out.size());

   epsilon = 1e-8;

   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = out[i].real();
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of out.real() failed on i "
            << i << std::endl;
   }
   v.clear();

   util::read_vector_from_file(v, params + "rfft1_imag.txt");

   ASSERT_EQ(v.size(), out.size());

   epsilon = 1e-8;

   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = out[i].imag();
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of out.imag() failed on i "
            << i << std::endl;
   }
   v.clear();

}


TEST(testRFFT, rfft2)
{
   std::string params = "../unit-tests/references/MyMath/fft/";
   ftype epsilon = 1e-8;



   std::vector<ftype> v, in;
   std::vector<complex_t> out;

   in.resize(256);
   mymath::linspace(in.data(), 0.f, 1000.f, in.size());
   mymath::rfft(in, 256, out);

   util::read_vector_from_file(v, params + "rfft2_real.txt");

   ASSERT_EQ(v.size(), out.size());

   epsilon = 1e-8;

   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = out[i].real();
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of out.real() failed on i "
            << i << std::endl;
   }
   v.clear();

   util::read_vector_from_file(v, params + "rfft2_imag.txt");

   ASSERT_EQ(v.size(), out.size());


   epsilon = 1e-8;

   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = out[i].imag();
      //printf("%+.8e\n",real);
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of out.imag() failed on i "
            << i << std::endl;
   }
   v.clear();

}

TEST(testRFFT, irfft)
{
   std::string params = "../unit-tests/references/MyMath/fft/";
   ftype epsilon = 1e-8;

   std::vector<ftype> v, a, b;
   //std::vector<complex_t> outA, outB;
   a.resize(256);
   b.resize(256);
   mymath::linspace(a.data(), 0.f, 100.f, a.size());
   mymath::linspace(b.data(), 0.f, 1000.f, b.size());

   std::vector<complex_t> az, bz;

   mymath::real_to_complex(a, az);
   mymath::real_to_complex(b, bz);
   //printf("ok here\n");

   mymath::fft(az, 512, az);
   //printf("ok here\n");

   mymath::fft(bz, 512, bz);

   std::transform(az.begin(), az.end(), bz.begin(),
                  az.begin(), std::multiplies<complex_t>());

   //printf("ok here\n");

   //for (unsigned int i = 0; i < outA.size(); ++i) {
   //   printf("outA * outB: %+.8e\n", std::abs(outA[i]));
   //}
   //az.clear();
   mymath::ifft(az, 512, az);
   //util::dump(a.data(), 10, "inverse complex fft");
   //printf("ok here\n");

   util::read_vector_from_file(v, params + "inverse_rfft.txt");

   ASSERT_EQ(v.size(), az.size());


   // WARNING!! Absolute difference is used on pusrpose!!
   epsilon = 1e-8;
   // !!!!
   int j = 0;
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = az[i].real();
      //ASSERT_NEAR(ref, real, epsilon /** std::max(fabs(ref), fabs(real))*/)
      if (std::max(fabs(ref), fabs(real)) < 1e-8) {
         /*
         ASSERT_DOUBLE_EQ(std::trunc(ref / epsilon), std::trunc(real/epsilon))
            << "Testing of az.real() failed on i "
            << i << std::endl;
         */
         j++;
      } else {
         //printf("%lf\n",real);
         ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
               << "Testing of az.real failed on i "
               << i << std::endl;
      }


   }
   //printf("%d\n", j);
   if (100.0 * j / v.size() > 10.0) {
      printf("Test leaves out more than 10 %% of data\n");
      printf("Maybe you should reconsider it?\n");

   }
   v.clear();


}

TEST(testConvolution, test1)
{
   std::vector<ftype> c, a, b, v;
   a.resize(10);
   b.resize(20);
   mymath::linspace(a.data(), 0.f, 100.f, a.size());
   mymath::linspace(b.data(), 100.f, 1000.f, b.size());
   c = mymath::convolution(a, b);


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
   c = mymath::convolution(a, b);


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
   a = {0,0,0,0,1,1,1,0,0,0};
   b = {1,1,1,1,0,0,0,1,1,1};
   c = mymath::convolution(a, b);


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


int main(int ac, char *av[])
{
   ::testing::InitGoogleTest(&ac, av);
   return RUN_ALL_TESTS();
}
