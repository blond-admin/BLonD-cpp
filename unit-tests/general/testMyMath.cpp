#include <blond/blond.h>

#include <gtest/gtest.h>
#include <testing_utilities.h>

using namespace mymath;
using namespace std;


TEST(testLinspace, test1)
{
    double a[10];
    linspace(a, 0, -10, 10);
    ASSERT_DOUBLE_EQ(a[9], -10);
    ASSERT_DOUBLE_EQ(a[0], 0);
}

TEST(testLinspace, test2)
{
    double a[10];
    linspace(a, 1, 1, 10);
    ASSERT_DOUBLE_EQ(a[0], 1);
    ASSERT_DOUBLE_EQ(a[9], 1);
}

TEST(testLinspace, test3)
{
    auto a = linspace(-10.0, 10.0, 10);
    ASSERT_NEAR_LOOP({ -10.0, -7.777778, -5.555556, -3.3333333, -1.111111,
                       1.1111111, 3.33333333, 5.55555556, 7.777777778, 10.0
                     }, a, "linspace", 1e-5);
}


TEST(testMean, test1)
{
    double a[5] = { -2, -1, 0, 1, 2};
    double m = mean(a, 5);
    ASSERT_DOUBLE_EQ(m, 0);

    double b[5] = {1, 1, 1, 1, 1};
    m = mean(b, 5);
    ASSERT_DOUBLE_EQ(m, 1);

    double c[5] = {10, 12, 14, 16, 18};
    m = mean(c, 5);
    ASSERT_DOUBLE_EQ(m, 14);
}

TEST(testSTD, test1)
{
    double a[5] = {1, 1, 1, 1, 1};
    double b = standard_deviation(a, 5);
    ASSERT_DOUBLE_EQ(b, 0);

    double c[5] = {10, 12.4, 131, -22.3, 45.6};
    double m = mean(c, 5);
    double std1 = standard_deviation(c, 5, m);
    double std2 = standard_deviation(c, 5);
    ASSERT_DOUBLE_EQ(std1, std2);
}

TEST(testMin, double_type)
{
    double b[10] = {0.8, 9.14, 2.1, 2, 10.25, 1.3, 2.12, 3, 4, 5};
    int i = min<double>(b, 10);
    ASSERT_DOUBLE_EQ(b[i], 0.8);

    double a[5] = {10.1, 1.2, -1.3, 1.4, 1.5};
    i = min<double>(a, 5);
    ASSERT_DOUBLE_EQ(a[i], -1.3);
}

TEST(testMin, int_type)
{
    int b[10] = {1, 9, 2, 2, 10, 0, 2, 3, 4, 5};
    int i = min<int>(b, 10);
    ASSERT_EQ(b[i], 0);

    int a[5] = {10, 1, -1, 1, 1};
    i = min<int>(a, 5);
    ASSERT_EQ(a[i], -1);
}

TEST(testMax, double_type)
{
    double b[10] = {0.8, 9.14, 2.1, 2, 10.25, 1.3, 2.12, 3, 4, 5};
    int i = max<double>(b, 10);
    ASSERT_DOUBLE_EQ(b[i], 10.25);

    double a[5] = {10.1, 1.2, -1.3, 1.4, 1.5};
    i = max<double>(a, 5);
    ASSERT_DOUBLE_EQ(a[i], 10.1);
}

TEST(testMax, int_type)
{
    int b[10] = {1, 9, 2, 2, 10, 0, 2, 3, 4, 5};
    int i = max<int>(b, 10);
    ASSERT_EQ(b[i], 10);

    int a[5] = {10, 1, -1, 100, 1};
    i = max<int>(a, 5);
    ASSERT_EQ(a[i], 100);
}

TEST(testTrapezoid, test1)
{
    double b[5] = {1, 2, 3, 4, 5};
    double trap = trapezoid<double>(b, 1, 5);
    ASSERT_DOUBLE_EQ(trap, 12);

    double a[5] = {1.1, 1.2, 1.3, 1.4, 1.5};
    trap = trapezoid<double>(a, 0.1, 5);
    ASSERT_DOUBLE_EQ(trap, 0.52);

    trap = trapezoid<double>(a, 1, 5);
    ASSERT_DOUBLE_EQ(trap, 5.2);

    double c[10] = { -0.61, -0.51, 0.39, -0.54, 0.67, 1.4, 1.1, 1.4, 0.16, 0.9};
    trap = trapezoid<double>(c, 1, 10);
    ASSERT_NEAR(trap, 4.215, 1e-8);
}



TEST(testTrapezoid, test2)
{
    double a[5] = {1.1, 1.2, 1.3, 1.4, 1.5};
    double trap = trapezoid<double>(a, a, 5);
    ASSERT_NEAR(trap, 0.52, 1e-8);

    double c[10] = { -0.61, -0.51, 0.39, -0.54, 0.67, 1.4, 1.1, 1.4, 0.16, 0.9};
    trap = trapezoid<double>(c, c, 10);
    ASSERT_NEAR(trap, 0.21895, 1e-8);
}


TEST(testTrapezoid, test3)
{
    int b[5] = {1, 2, 3, 4, 5};
    double trap = trapezoid(b, 1.0, 5);
    ASSERT_DOUBLE_EQ(trap, 12);
}

TEST(testTrapezoid, test4)
{
    int_vector_t a = {0, 0, 0, 1, 0, 9, 3, 9,
                      50, 70, 90, 10, 5, 6, 2,
                      1, 1, 0, 0
                     };

    double trap = trapezoid(a.data(), 1.0, a.size());
    ASSERT_DOUBLE_EQ(trap, 257.0);

}


TEST(testTrapezoid, test5)
{
    int_vector_t a = {0, 0, 0, 1, 0, 9, 3, 9,
                      50, 70, 90, 10, 5, 6, 2,
                      1, 1, 0, 0
                     };

    double trap = trapezoid(a.data(), 0.12, a.size());
    ASSERT_NEAR(trap, 30.84, 1e-5);

}


TEST(testTrapezoid, test6)
{
    int_vector_t a = {0, 0, 0, 0, 1, 0, 1, 0,
                      1, 1, 0, 1, 0, 4, 9, 3,
                      9, 6, 12, 19, 23, 33, 31,
                      31, 37, 49, 53, 80, 85, 77,
                      110, 111, 142, 153, 168, 196,
                      208, 228, 247, 242, 257, 274,
                      282, 301, 265, 356, 310, 363,
                      334, 300, 353, 352, 312, 307,
                      259, 312, 246, 261, 245, 198,
                      211, 188, 171, 193, 139, 117,
                      109, 94, 86, 70, 82, 46, 37,
                      39, 34, 22, 24, 15, 10, 12,
                      9, 9, 6, 6, 2, 0, 3, 1, 3,
                      1, 2, 0, 0, 0, 0, 1, 0, 0,
                      0, 0
                     };

    double trap = trapezoid(a.data(), 2.10324675e-10, a.size());
    ASSERT_NEAR(trap, 2.10324675e-06, 1e-15);

}



TEST(testCumTrap, test1)
{
    std::string params = TEST_FILES "/MyMath/CumTrap/";
    double epsilon = 1e-8;

    double a[10] = { -0.61, -0.51, 0.39, -0.54, 0.67, 1.4, 1.1, 1.4, 0.16, 0.9};
    auto trap = cum_trapezoid<double>(a, 1, 0, 10);

    std::vector<double> v;
    util::read_vector_from_file(v, params + "cumtrap_reference.txt");

    for (unsigned int i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = trap[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)));
    }

    // delete[] trap;
}


TEST(testCumTrap, test2)
{
    // std::string params = TEST_FILES "/MyMath/CumTrap/";
    double epsilon = 1e-8;

    double a[10] = { -0.61, -0.51, 0.39, -0.54, 0.67, 1.4, 1.1, 1.4, 0.16, 0.9};
    auto trap = cum_trapezoid<double>(a, 0.51, 10);

    f_vector_t v{ -0.2856, -0.3162, -0.35445,
                  -0.3213, 0.20655, 0.84405,
                  1.48155, 1.87935, 2.14965};
    // util::read_vector_from_file(v, params + "cumtrap_reference.txt");

    for (unsigned int i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = trap[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)));
    }

    // delete[] trap;
}


TEST(testConvolution, test1)
{
    std::vector<double> c, a, b, v;
    a.resize(10);
    b.resize(20);
    linspace(a.data(), 0.f, 100.f, a.size());
    linspace(b.data(), 100.f, 1000.f, b.size());
    c.resize(a.size() + b.size() - 1);
    convolution(a.data(), a.size(), b.data(), b.size(), c.data());

    std::string params = TEST_FILES "/MyMath/convolution/";
    double epsilon = 1e-8;

    util::read_vector_from_file(v, params + "convolution1.txt");

    ASSERT_EQ(v.size(), c.size());

    epsilon = 1e-8;
    for (unsigned int i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = c[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of convolution1 failed on i " << i << std::endl;
    }

    v.clear();
}

TEST(testConvolution, test2)
{
    std::vector<double> c, a, b, v;
    a.resize(100, 0);
    b.resize(20);
    linspace(b.data(), 100.f, 1000.f, b.size());
    c.resize(a.size() + b.size() - 1);
    convolution(a.data(), a.size(), b.data(), b.size(), c.data());

    std::string params = TEST_FILES "/MyMath/convolution/";
    double epsilon = 1e-8;

    util::read_vector_from_file(v, params + "convolution2.txt");

    ASSERT_EQ(v.size(), c.size());

    epsilon = 1e-8;
    for (unsigned int i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = c[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of convolution2 failed on i " << i << std::endl;
    }

    v.clear();
}

TEST(testConvolution, test3)
{
    std::vector<double> c, a, b, v;
    a = {0, 0, 0, 0, 1, 1, 1, 0, 0, 0};
    b = {1, 1, 1, 1, 0, 0, 0, 1, 1, 1};
    c.resize(a.size() + b.size() - 1);
    convolution(a.data(), a.size(), b.data(), b.size(), c.data());

    std::string params = TEST_FILES "/MyMath/convolution/";
    double epsilon = 1e-8;

    util::read_vector_from_file(v, params + "convolution3.txt");

    ASSERT_EQ(v.size(), c.size());

    epsilon = 1e-8;
    for (unsigned int i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = c[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of convolution3 failed on i " << i << std::endl;
    }

    v.clear();
}

TEST(arange, test1)
{
    std::string params = TEST_FILES "/MyMath/arange/";
    f_vector_t v;

    auto a = arange<double>(0, 100);

    double epsilon = 1e-8;

    util::read_vector_from_file(v, params + "arange1.txt");

    ASSERT_EQ(v.size(), a.size());

    for (unsigned int i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = a[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of a failed on i " << i << std::endl;
    }
}

TEST(arange, test2)
{
    std::string params = TEST_FILES "/MyMath/arange/";
    f_vector_t v;
    auto a = arange<double>(0, 100, 2.5);

    double epsilon = 1e-8;

    util::read_vector_from_file(v, params + "arange2.txt");

    ASSERT_EQ(v.size(), a.size());

    for (unsigned int i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = a[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of a failed on i " << i << std::endl;
    }
}

TEST(arange, test3)
{
    std::string params = TEST_FILES "/MyMath/arange/";
    f_vector_t v;

    auto a = arange<int>(10, 80, 2);
    util::read_vector_from_file(v, params + "arange3.txt");
    // util::dump(v, "integer vector");

    ASSERT_EQ(v.size(), a.size());

    for (unsigned int i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = a[i];
        ASSERT_EQ(ref, real) << "Testing of a failed on i " << i << std::endl;
    }
}


TEST(interp, test1)
{
    auto epsilon = 1e-8;
    f_vector_t x{1.0, 4.0, 5.0, 7.0, 10.0};
    f_vector_t xp{0.0, 4.5, 9.0, 12.0};
    f_vector_t fp{10.0, 19.0, 28.0, 34.0};
    f_vector_t y;
    interp(x, xp, fp, y, fp.front(), fp.back());

    f_vector_t v{12.0, 18.0, 20., 24., 30.};
    ASSERT_EQ(v.size(), y.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = y[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of y failed on i " << i << std::endl;
    }
}


TEST(interp, test2)
{
    std::string params = TEST_FILES "/MyMath/lin_interp/test2/";

    auto epsilon = 1e-8;
    auto x = arange(0.0, 100.0, 0.3);
    f_vector_t xp{ 0, 20, 50, 80, 100};
    f_vector_t fp{ -1, 12, 0, -42, 42};
    f_vector_t y;
    interp(x, xp, fp, y, fp.front(), fp.back());

    double max = *max_element(y.begin(), y.end(), [](double i, double j) {
        return std::abs(i) < std::abs(j);
    });
    max = std::abs(max);

    f_vector_t v;
    util::read_vector_from_file(v, params + "interp.txt");
    ASSERT_EQ(v.size(), y.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = y[i];
        ASSERT_NEAR(ref, real, epsilon * max)
                << "Testing of y failed on i " << i << std::endl;
    }
}


TEST(interp, test3)
{
    std::string params = TEST_FILES "/MyMath/lin_interp/test3/";

    auto epsilon = 1e-8;
    auto x = arange(1.0, 99.0, 0.5);
    f_vector_t xp{10, 25, 55, 56, 70, 90};
    f_vector_t fp{ -3, 4, 27, 41, -12, 13};
    f_vector_t y;
    interp(x, xp, fp, y, fp.front(), fp.back());

    double max = *max_element(y.begin(), y.end(), [](double i, double j) {
        return std::abs(i) < std::abs(j);
    });
    max = std::abs(max);

    f_vector_t v;
    util::read_vector_from_file(v, params + "interp.txt");
    ASSERT_EQ(v.size(), y.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = y[i];
        EXPECT_NEAR(ref, real, epsilon * max)
                << "Testing of y failed on i " << i << std::endl;
    }
}


TEST(interp, test4)
{
    f_vector_t x = {2577.916649096185, 2577.916649096185,
                    2420.3176528693307, 1542.92854183888,
                    492.7387950151245, 31.438007915807102,
                    493.6019188781009, 1544.0287765342725,
                    2420.8570118212447, 2577.916649096186
                   };

    f_vector_t xp = {0.0, 149.26952997505427,
                     561.6140644052764, 1140.5209174462839,
                     1750.49167701547, 2248.757269331754,
                     2518.69415854002, 2546.478641180378,
                     2546.478641180378, 2546.478641180378
                    };

    f_vector_t yp = {8.860503368317807e-13, 3.336738187768568e-30,
                     2.8802773054147077e-59, 5.543815268967439e-100,
                     2.471804120154838e-152, 2.8346334050399317e-216,
                     1.5622782804689192e-291,
                     0.0, 0.0, 0.0
                    };

    f_vector_t ref = {0.0, 0.0, 1.033061209466945e-216,
                      1.886470228309715e-100, 5.573463993757814e-31,
                      6.99437184718444e-13, 5.50361904338525e-31,
                      1.8764705720834508e-100, 1.0273973488998069e-216, 0.0
                     };

    auto epsilon = 1e-8;
    auto real = interp(x, xp, yp);

    ASSERT_NEAR_LOOP(ref, real, "interp", epsilon);
}


TEST(flatten, test1)
{
    auto epsilon = 1e-8;
    f_vector_2d_t a(5, f_vector_t(5));
    f_vector_t b;
    for (uint i = 0; i < a.size(); i++) {
        for (uint j = 0; j < a[i].size(); j++) {
            a[i][j] = i * j;
            b.push_back(i * j);
        }
    }
    // cout << "flatten(a): " << flatten(a);
    // cout << "b: " << b;
    ASSERT_NEAR_LOOP(flatten(a), b, "flatten");
}


TEST(flatten, test2)
{
    auto epsilon = 1e-8;
    f_vector_2d_t a{{1, 2, 3, 4}, {5, 5, 5}, { -2, 2}};
    f_vector_t b {1, 2, 3, 4, 5, 5, 5, -2, 2};
    ASSERT_NEAR_LOOP(flatten(a), b, "flatten");
}


TEST(meshgrid, test1)
{
    int_vector_t a {1, 2, 3};
    int_vector_t b {10, 20};
    int_vector_2d_t c, d;
    meshgrid(a, b, c, d);
    ASSERT_EQ_LOOP(flatten(c), {1, 2, 3, 1, 2, 3}, "meshgrid 1st ret");
    ASSERT_EQ_LOOP(flatten(d), {10, 10, 10, 20, 20, 20}, "meshgrid 2nd ret");
}

TEST(meshgrid, test2)
{
    int_vector_t a {1, 2, 3};
    int_vector_t b {1, 2, 3};
    int_vector_2d_t c, d;
    meshgrid(a, b, c, d);
    ASSERT_EQ_LOOP(flatten(c), {1, 2, 3, 1, 2, 3, 1, 2, 3}, "meshgrid 1st ret");
    ASSERT_EQ_LOOP(flatten(d), {1, 1, 1, 2, 2, 2, 3, 3, 3}, "meshgrid 2nd ret");
}

TEST(meshgrid, test3)
{
    int_vector_t a {1, 2, 3};
    int_vector_t b {1, 2, 3};
    int_vector_2d_t c;
    meshgrid(a, b, c, c);
    ASSERT_EQ_LOOP(flatten(c), {1, 1, 1, 2, 2, 2, 3, 3, 3}, "meshgrid 1st-2nd ret");
}


TEST(random_choice, test1)
{
    auto epsilon = 1e-2;
    const int size = 2000000;
    int_vector_t elems {0,  1,  2,  3,  4,  5,  6,  7,  8,  9};
    f_vector_t weights {1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
    auto res = random_choice(elems, size, weights);
    f_vector_t freq(elems.size(), 0.);
    FOR(res, it) { freq[*it]++; }
    freq /= size;

    ASSERT_NEAR_LOOP(freq, weights / accumulate(ALL(weights), 0.0),
                     "frequency", epsilon);
}

TEST(random_choice, test2)
{
    auto epsilon = 1e-2;
    const int size = 2000000;
    int_vector_t elems {0,  1,  2,  3,  4,  5,  6,  7,  8,  9};
    f_vector_t weights {1., 0., 1., 0., 0., 1., 0., 1., 1., 0.};
    auto res = random_choice(elems, size, weights);
    f_vector_t freq(elems.size(), 0.);
    FOR(res, it) { freq[*it]++; }
    freq /= size;

    ASSERT_NEAR_LOOP(freq, weights / accumulate(ALL(weights), 0.0),
                     "frequency", epsilon);
}


TEST(random_choice, test3)
{
    auto epsilon = 1e-2;
    const int size = 2000000;
    int_vector_t elems {0,  1,  2,  3,  4,  5,  6,  7,  8,  9};
    f_vector_t weights {0., 0., 1., 0., 0., 0., 0., 0., 0., 0.};
    auto res = random_choice(elems, size, weights);
    f_vector_t freq(elems.size(), 0.);
    FOR(res, it) { freq[*it]++; }
    freq /= size;

    ASSERT_NEAR_LOOP(freq, weights / accumulate(ALL(weights), 0.0),
                     "frequency", epsilon);
}


TEST(random_choice, test4)
{
    auto epsilon = 1e-2;
    const int size = 2000000;
    int_vector_t elems {0,  1,  2,  3,  4};
    auto res = random_choice(elems, size);
    f_vector_t freq(elems.size(), 0.);
    FOR(res, it) { freq[*it]++; }
    // freq /= size;

    ASSERT_NEAR_LOOP(freq, f_vector_t(5, size / 5), "frequency", epsilon);
}





int main(int ac, char *av[])
{
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
