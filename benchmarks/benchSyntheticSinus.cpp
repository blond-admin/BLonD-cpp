#ifdef __MINGW32__
#undef _GLIBCXX_DEBUG
#endif
#include <benchmark/benchmark.h>
#include <blond/math_functions.h>
#include <random>

const int max = 50000000;
const int min = 16;

static void CustomStep(benchmark::internal::Benchmark* b) {
    int i = min;
    while (i <= max) {
        i = i * 1.5;
        std::vector<int> args(1, i);
        b->Args(args);
    }
}

void fill(std::vector<double>& data) {
    auto size = data.size();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0, 2.0 * M_PI);
#pragma omp for
    for (int i = 0; i < size; ++i) {
        data[i] = dist(gen);
    }
}

static void sin_std(benchmark::State& state) {
    auto size = state.range(0);
    std::vector<double> data(size);
    fill(data);
    auto accumulator = 0.0;

    while (state.KeepRunning()) {
        state.PauseTiming();
        std::vector<double> v(size);
        state.ResumeTiming();
        for (int i = 0; i < size; ++i) {
            v[i] = sin(data[i]);
        }
        state.PauseTiming();
#pragma omp for
        for (int i = 0; i < size; ++i) {
            accumulator += v[i];
        }
        state.ResumeTiming();
    }
    std::cout << accumulator << std::endl;
}
BENCHMARK(sin_std)->Apply(CustomStep);

static void sin_std_parallel(benchmark::State& state) {
    auto size = state.range(0);
    std::vector<double> data(size);
    fill(data);
    auto accumulator = 0.0;

    while (state.KeepRunning()) {
        state.PauseTiming();
        std::vector<double> v(size);
        state.ResumeTiming();
#pragma omp for
        for (int i = 0; i < size; ++i) {
            v[i] = sin(data[i]);
        }
        state.PauseTiming();
#pragma omp for
        for (int i = 0; i < size; ++i) {
            accumulator += v[i];
        }
        state.ResumeTiming();
    }
    std::cout << accumulator << std::endl;
}
BENCHMARK(sin_std_parallel)->Apply(CustomStep);

static void sin_blond_parallel(benchmark::State& state) {
    auto size = state.range(0);
    std::vector<double> data(size);
    fill(data);
    auto accumulator = 0.0;
    while (state.KeepRunning()) {
        state.PauseTiming();
        std::vector<double> v(size);
        state.ResumeTiming();
#pragma omp for
        for (int i = 0; i < size; ++i) {
            v[i] = vdt::fast_sin(data[i]);
        }
        state.PauseTiming();
#pragma omp for
        for (int i = 0; i < size; ++i) {
            accumulator += v[i];
        }
        state.ResumeTiming();
    }
    std::cout << accumulator << std::endl;
}
BENCHMARK(sin_blond_parallel)->Apply(CustomStep);

static void sin_blond(benchmark::State& state) {
    auto size = state.range(0);
    std::vector<double> data(size);
    fill(data);
    auto accumulator = 0.0;

    while (state.KeepRunning()) {
        state.PauseTiming();
        std::vector<double> v(size);
        state.ResumeTiming();
        for (int i = 0; i < size; ++i) {
            v[i] = vdt::fast_sin(data[i]);
        }
        state.PauseTiming();
        for (int i = 0; i < size; ++i) {
            accumulator += v[i];
        }
        state.ResumeTiming();
    }
    std::cout << accumulator << std::endl;
}
BENCHMARK(sin_blond)->Apply(CustomStep);

#ifdef WITH_OPENCL
#include <boost/compute.hpp>

BOOST_COMPUTE_FUNCTION(double, gpu_sin, (double x), { return sin(x); });

static void sin_opencl(benchmark::State& state) {
    auto size = state.range(0);
    std::vector<double> data(size);
    fill(data);
    auto accumulator = 0.0;

    using namespace boost::compute;
    device gpu = system::default_device();
    context context(gpu);
    command_queue queue(context, gpu);

    while (state.KeepRunning()) {
        state.PauseTiming();
        std::vector<double> v(size);
        state.ResumeTiming();
        vector<double> device_out(v.size(), context);
        vector<double> device_in(data.size(), context);
        boost::compute::copy(data.begin(), data.end(), device_in.begin(),
                             queue);
        transform(device_in.begin(), device_in.end(), device_out.begin(),
                  gpu_sin, queue);
        // reduce(device_out.begin(), device_out.end(), &accumulator,
        // plus<double>());
        boost::compute::copy(device_out.begin(), device_out.end(), v.begin(),
                             queue);
        queue.finish();
        state.PauseTiming();
        for (int i = 0; i < size; ++i) {
            accumulator += v[i];
        }
        state.ResumeTiming();
    }
    std::cout << accumulator << std::endl;
}
BENCHMARK(sin_opencl)->Apply(CustomStep);

static void sin_opencl_short(benchmark::State& state) {
    auto size = state.range(0);
    std::vector<double> data(size);
    fill(data);
    auto accumulator = 0.0;

    using namespace boost::compute;
    device gpu = system::default_device();
    context context(gpu);
    command_queue queue(context, gpu);
    while (state.KeepRunning()) {
        state.PauseTiming();
        std::vector<double> v(size);
        state.ResumeTiming();
        vector<double> device_out(v.size(), context);
        vector<double> device_in(data.size(), context);
        boost::compute::copy(data.begin(), data.end(), device_in.begin(),
                             queue);
        transform(device_in.begin(), device_in.end(), device_out.begin(),
                  boost::compute::sin<double>(), queue);
        // reduce(device_out.begin(), device_out.end(), &accumulator,
        // plus<double>());
        boost::compute::copy(device_out.begin(), device_out.end(), v.begin(),
                             queue);
        queue.finish();
        state.PauseTiming();
        for (int i = 0; i < size; ++i) {
            accumulator += v[i];
        }
        state.ResumeTiming();
    }
    std::cout << accumulator << std::endl;
}
BENCHMARK(sin_opencl_short)->Apply(CustomStep);
#endif

int main(int argc, char** argv) {
    ::benchmark::Initialize(&argc, argv);
    omp_set_num_threads(omp_get_max_threads());
    ::benchmark::RunSpecifiedBenchmarks();
    return 0;
}