#include <iostream>

#include <blond/configuration.h>
#include <blond/fft.h>
#include <blond/math_functions.h>
#include <blond/utilities.h>
#include <gtest/gtest.h>

TEST(testIRFFT, rfft_even) {
    // std::string params = TEST_FILES"/MyMath/fft/";

    f_vector_t in, out4;
    complex_vector_t in2;
    complex_vector_t out1, out2, out3;

    // in.resize(256);
    // mymath::linspace(in.data(), 0.f, 100.f, in.size());

    in = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    // fft::rfft(in, 2 * in.size(), out1);
    fft::rfft(in, out1, in.size());

    // util::dump(out1.data(), out1.size(), "rfft\n");
    // in2.resize(in.size());
    fft::real_to_complex(in, in2);

    fft::fft(in2, out2, in2.size());
    // fft::fft(in2, 2 * in2.size(), out2);
    // util::dump(out2.data(), out2.size(), "fft\n");

    fft::ifft(out2, out3, out2.size());
    // util::dump(out3.data(), out3.size(), "ifft after fft\n");

    fft::irfft(out1, out4);
    // util::dump(out4.data(), out4.size(), "ifft after rfft\n");

    ASSERT_EQ(out4.size(), out3.size());

    auto epsilon = 1e-8;

    for (uint i = 0; i < out3.size(); ++i) {
        ftype ref = out3[i].real();
        ftype real = out4[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of out4 failed on i " << i << std::endl;
    }

    fft::destroy_plans();
}

TEST(testIRFFT, rfft_odd) {
    // std::string params = TEST_FILES"/MyMath/fft/";

    f_vector_t in, out4;
    complex_vector_t in2;
    complex_vector_t out1, out2, out3;

    // in.resize(256);
    // mymath::linspace(in.data(), 0.f, 100.f, in.size());

    in = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    // fft::rfft(in, 2 * in.size(), out1);
    fft::rfft(in, out1, in.size());

    // util::dump(out1.data(), out1.size(), "rfft\n");

    // util::dump(in.data(), in.size(), "in\n");
    // in2.resize(in.size());
    fft::real_to_complex(in, in2);

    // util::dump(in2.data(), in2.size(), "in2\n");

    fft::fft(in2, out2, in2.size());
    // fft::fft(in2, 2 * in2.size(), out2);
    // util::dump(out2.data(), out2.size(), "fft\n");

    fft::ifft(out2, out3, out2.size());
    // util::dump(out3.data(), out3.size(), "ifft after fft\n");

    fft::irfft(out1, out4, 15);
    // util::dump(out4.data(), out4.size(), "ifft after rfft\n");

    ASSERT_EQ(out4.size(), out3.size());

    auto epsilon = 1e-8;

    for (uint i = 0; i < out3.size(); ++i) {
        ftype ref = out3[i].real();
        ftype real = out4[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of out4 failed on i " << i << std::endl;
    }
    fft::destroy_plans();
}

TEST(testIRFFT, irfft_big) {
    std::string params = TEST_FILES"/MyMath/fft/irfft/";

    complex_vector_t in2;
    f_vector_t out;
    f_vector_t in1(10000);

    for (uint i = 0; i < in1.size(); ++i) {
        in1[i] = i;
    }

    fft::rfft(in1, in2, in1.size());

    fft::irfft(in2, out);

    std::vector<ftype> v;

    util::read_vector_from_file(v, params + "irfft1.txt");

    ASSERT_EQ(v.size(), out.size());

    ftype max = *max_element(out.begin(), out.end(), [](ftype i, ftype j) {
        return fabs(i) < fabs(j);
    });
    ftype epsilon = 1e-9 * max;

    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = out[i];

        ASSERT_NEAR(ref, real, epsilon) << "Testing of irfft failed on i " << i
                                        << std::endl;
    }
    fft::destroy_plans();
}

TEST(testIRFFT, irfft_big2) {

    std::string params = TEST_FILES"/MyMath/fft/irfft/";

    complex_vector_t in(101);
    f_vector_t out;

    for (uint i = 0; i < in.size(); ++i) {
        if (i < 100)
            in[i] = complex_t(std::sin(i), std::sqrt(i));
        else
            in[i] = complex_t(std::sin(i), 0.0);
    }

    fft::irfft(in, out);

    std::vector<ftype> v;

    util::read_vector_from_file(v, params + "irfft2.txt");

    ASSERT_EQ(v.size(), out.size());

    ftype max = *max_element(out.begin(), out.end(), [](ftype i, ftype j) {
        return fabs(i) < fabs(j);
    });
    ftype epsilon = 1e-9 * max;

    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = out[i];

        ASSERT_NEAR(ref, real, epsilon) << "Testing of irfft failed on i " << i
                                        << std::endl;
    }
    fft::destroy_plans();
}

TEST(testIRFFT, irfft_test) {

    std::string params = TEST_FILES"/MyMath/fft/irfft/";

    complex_vector_t in(101);
    f_vector_t out;

    for (uint i = 0; i < in.size(); ++i) {
        if (i < 100)
            in[i] = complex_t(std::sin(i), std::sqrt(i));
        else
            in[i] = complex_t(std::sin(i), 0.0);
    }

    fft::irfft(in, out);

    std::vector<ftype> v;

    util::read_vector_from_file(v, params + "irfft2.txt");

    ASSERT_EQ(v.size(), out.size());

    ftype max = *max_element(out.begin(), out.end(), [](ftype i, ftype j) {
        return fabs(i) < fabs(j);
    });
    ftype epsilon = 1e-9 * max;

    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = out[i];

        ASSERT_NEAR(ref, real, epsilon) << "Testing of irfft failed on i " << i
                                        << std::endl;
    }
    fft::destroy_plans();
}

TEST(testRFFT, rfft1) {
    std::string params = TEST_FILES"/MyMath/fft/";
    ftype epsilon = 1e-8;

    std::vector<ftype> v, in;
    std::vector<complex_t> out;
    in.resize(256);
    mymath::linspace(in.data(), 0.f, 100.f, in.size());
    fft::rfft(in, out, 512);

    util::read_vector_from_file(v, params + "rfft1_real.txt");

    ASSERT_EQ(v.size(), out.size());

    epsilon = 1e-8;

    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = out[i].real();
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of out.real() failed on i " << i << std::endl;
    }
    v.clear();

    util::read_vector_from_file(v, params + "rfft1_imag.txt");

    ASSERT_EQ(v.size(), out.size());

    epsilon = 1e-8;

    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = out[i].imag();
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of out.imag() failed on i " << i << std::endl;
    }
    v.clear();
    fft::destroy_plans();
}

TEST(testRFFT, rfft2) {
    std::string params = TEST_FILES"/MyMath/fft/";
    ftype epsilon = 1e-8;

    std::vector<ftype> v, in;
    std::vector<complex_t> out;

    in.resize(256);
    mymath::linspace(in.data(), 0.f, 1000.f, in.size());
    fft::rfft(in, out, 256);

    util::read_vector_from_file(v, params + "rfft2_real.txt");

    ASSERT_EQ(v.size(), out.size());

    epsilon = 1e-8;

    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = out[i].real();
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of out.real() failed on i " << i << std::endl;
    }
    v.clear();

    util::read_vector_from_file(v, params + "rfft2_imag.txt");

    ASSERT_EQ(v.size(), out.size());

    epsilon = 1e-8;

    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = out[i].imag();
        // printf("%+.8e\n",real);
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of out.imag() failed on i " << i << std::endl;
    }
    v.clear();
    fft::destroy_plans();
}

TEST(testRFFT, irfft) {
    std::string params = TEST_FILES"/MyMath/fft/";
    ftype epsilon = 1e-8;

    std::vector<ftype> v, a, b;
    // std::vector<complex_t> outA, outB;
    a.resize(256);
    b.resize(256);
    mymath::linspace(a.data(), 0.f, 100.f, a.size());
    mymath::linspace(b.data(), 0.f, 1000.f, b.size());

    std::vector<complex_t> az, bz;

    fft::real_to_complex(a, az);
    fft::real_to_complex(b, bz);
    // printf("ok here\n");

    fft::fft(az, az, 512);
    // printf("ok here\n");

    fft::fft(bz, bz, 512);

    std::transform(az.begin(), az.end(), bz.begin(), az.begin(),
                   std::multiplies<complex_t>());

    // printf("ok here\n");

    // for (unsigned int i = 0; i < outA.size(); ++i) {
    //   printf("outA * outB: %+.8e\n", std::abs(outA[i]));
    //}
    // az.clear();
    fft::ifft(az, az, 512);
    // util::dump(a.data(), 10, "inverse complex fft");
    // printf("ok here\n");

    util::read_vector_from_file(v, params + "inverse_rfft.txt");

    ASSERT_EQ(v.size(), az.size());

    // WARNING!! Absolute difference is used on pusrpose!!
    epsilon = 1e-8;
    // !!!!
    int j = 0;
    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = az[i].real();
        // ASSERT_NEAR(ref, real, epsilon /** std::max(fabs(ref), fabs(real))*/)
        if (std::max(fabs(ref), fabs(real)) < 1e-8) {
            /*
            ASSERT_DOUBLE_EQ(std::trunc(ref / epsilon),
            std::trunc(real/epsilon))
               << "Testing of az.real() failed on i "
               << i << std::endl;
            */
            j++;
        } else {
            // printf("%lf\n",real);
            ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
                << "Testing of az.real failed on i " << i << std::endl;
        }
    }
    // printf("%d\n", j);
    if (100.0 * j / v.size() > 10.0) {
        printf("Test leaves out more than 10 %% of data\n");
        printf("Maybe you should reconsider it?\n");
    }
    v.clear();
    fft::destroy_plans();
}

TEST(testRFFTFREQ, even_no_spacing) {
    ftype epsilon;

    std::vector<ftype> real, res;
    real = fft::rfftfreq(10);
    res = {0, 0.1, 0.2, 0.3, 0.4, 0.5};

    ASSERT_EQ(real.size(), res.size());
    epsilon = 1e-7;
    for (unsigned int i = 0; i < res.size(); ++i) {
        ftype ref = res[i];
        ftype real2 = real[i];
        ASSERT_NEAR(ref, real2, epsilon * std::max(fabs(ref), fabs(real2)))
            << "Testing of rfftfreq failed on i " << i << std::endl;
    }
}

TEST(testRFFTFREQ, odd_no_spacing) {
    ftype epsilon;

    std::vector<ftype> real, res;
    real = fft::rfftfreq(11);
    res = {0, 0.09090909, 0.18181818, 0.27272727, 0.36363636, 0.45454545};

    ASSERT_EQ(real.size(), res.size());
    epsilon = 1e-7;
    for (unsigned int i = 0; i < res.size(); ++i) {
        ftype ref = res[i];
        ftype real2 = real[i];
        ASSERT_NEAR(ref, real2, epsilon * std::max(fabs(ref), fabs(real2)))
            << "Testing of rfftfreq failed on i " << i << std::endl;
    }
}

TEST(testRFFTFREQ, even_spacing) {
    ftype epsilon;

    std::vector<ftype> real, res;
    real = fft::rfftfreq(10, 1.5);
    res = {0, 0.066666667, 0.133333333, 0.2, 0.26666666667, 0.3333333};

    ASSERT_EQ(real.size(), res.size());
    epsilon = 1e-7;
    for (unsigned int i = 0; i < res.size(); ++i) {
        ftype ref = res[i];
        ftype real2 = real[i];
        ASSERT_NEAR(ref, real2, epsilon * std::max(fabs(ref), fabs(real2)))
            << "Testing of rfftfreq failed on i " << i << std::endl;
    }
}

TEST(testRFFTFREQ, odd_spacing) {
    ftype epsilon;

    std::vector<ftype> real, res;
    real = fft::rfftfreq(11, 2.5);
    res = {0, 0.03636364, 0.07272727, 0.10909091, 0.14545455, 0.18181818};

    ASSERT_EQ(real.size(), res.size());
    epsilon = 1e-7;
    for (unsigned int i = 0; i < res.size(); ++i) {
        ftype ref = res[i];
        ftype real2 = real[i];
        ASSERT_NEAR(ref, real2, epsilon * std::max(fabs(ref), fabs(real2)))
            << "Testing of rfftfreq failed on i " << i << std::endl;
    }
}

int main(int ac, char* av[]) {
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
