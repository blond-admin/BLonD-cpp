#include <iostream>
#include <blond/plots/plot_parameters.h>
#include <gtest/gtest.h>

using namespace std;


TEST(testPlotParameters, plot_voltage_programme1)
{
    f_vector_t time;
    f_vector_t voltage;
    for (int i = 0; i < 100; i++) {
        time.push_back(i);
        voltage.push_back(10.0*i + 1.0);
    }
    plot_voltage_programme(time, voltage);
}


int main(int ac, char *av[])
{
    python::initialize();
    ::testing::InitGoogleTest(&ac, av);
    auto ret = RUN_ALL_TESTS();
    python::finalize();
    return ret;
}
