#undef _GLIBCXX_DEBUG
#include <benchmark/benchmark.h>
#include <blond/math_functions.h>
#include <boost/compute.hpp>
#include <random>

const int max = 32000000;
const int min = 16;

static void CustomStep(benchmark::internal::Benchmark* b) {
	int i = min;
	while ( i <= max) {
		i = i*1.5;
        std::vector<int> args(1, i);
        b->Args(args);
	}
}

void fill(std::vector<double> & data) {
    auto size = data.size();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0, 2.0*M_PI);
    for(int i = 0; i < size; ++i) {
        data[i] = dist(gen);
    }
}

static void sin_std(benchmark::State& state) {
    auto size = state.range(0);
	std::vector<double> data(size);
    fill(data);
	std::vector<double> v(size);
    auto accumulator = 0.0;

    while (state.KeepRunning()) {
        for(int i=0; i< size; ++i) {
            v[i] = sin(data[i]);
        }
        for(int i=0; i< size; ++i) {
            accumulator += v[i];
        }
    }
    std::cout << accumulator << std::endl;
}BENCHMARK(sin_std)->Apply(CustomStep);

static void sin_blond_par(benchmark::State& state) {
    auto size = state.range(0);
    std::vector<double> data(size);
    fill(data);
    std::vector<double> v(size);
    auto accumulator = 0.0;

    while (state.KeepRunning()) {
#pragma omp for
        for(int i=0; i< size; ++i) {
            v[i] = mymath::fast_sin(data[i]);
        }
        for(int i=0; i< size; ++i) {
            accumulator += v[i];
        }
    }
    std::cout << accumulator<< std::endl;
}BENCHMARK(sin_blond_par)->Apply(CustomStep);

static void sin_blond(benchmark::State& state) {
	auto size = state.range(0);
	std::vector<double> data(size);
	fill(data);
	std::vector<double> v(size);
	auto accumulator = 0.0;

	while (state.KeepRunning()) {
		for (int i = 0; i< size; ++i) {
			v[i] = mymath::fast_sin(data[i]);
		}
		for (int i = 0; i< size; ++i) {
			accumulator += v[i];
		}
	}
	std::cout << accumulator << std::endl;
}BENCHMARK(sin_blond)->Apply(CustomStep);

BOOST_COMPUTE_FUNCTION(double, gpu_sin, (double x), {
   return sin(x);
});

static void sin_opencl(benchmark::State& state) {
    auto size = state.range(0);
    std::vector<double> data(size);
    fill(data);
    std::vector<double> v(size);
    auto accumulator = 0.0;


    using namespace boost::compute;
    device gpu = system::default_device();
    context context(gpu);
    command_queue queue(
            context, gpu, command_queue::enable_profiling
    );
    vector<double> device_in(data.size(), context);
    boost::compute::copy(
            data.begin(), data.end(), device_in.begin(), queue
    );
    vector<double> device_out(v.size(), context);

    while (state.KeepRunning()) {

        transform(device_in.begin(), device_in.end(), device_out.begin(), gpu_sin, queue);
        //reduce(device_out.begin(), device_out.end(), &accumulator, plus<double>());
        boost::compute::copy(device_out.begin(), device_out.end(), v.begin(), queue);
        queue.finish();
        for(int i=0; i< size; ++i) {
            accumulator += v[i];
        }
    }
    std::cout << accumulator<< std::endl;
}BENCHMARK(sin_opencl)->Apply(CustomStep);

static void sin_opencl_short(benchmark::State& state) {
	auto size = state.range(0);
	std::vector<double> data(size);
	fill(data);
	std::vector<double> v(size);
	auto accumulator = 0.0;


	using namespace boost::compute;
	device gpu = system::default_device();
	context context(gpu);
	command_queue queue(
		context, gpu, command_queue::enable_profiling
	);
	vector<double> device_in(data.size(), context);
	boost::compute::copy(
		data.begin(), data.end(), device_in.begin(), queue
	);
	vector<double> device_out(v.size(), context);

	while (state.KeepRunning()) {

		transform(device_in.begin(), device_in.end(), device_out.begin(), boost::compute::sin<double>(), queue);
		//reduce(device_out.begin(), device_out.end(), &accumulator, plus<double>());
		boost::compute::copy(device_out.begin(), device_out.end(), v.begin(), queue);
		queue.finish();
		for (int i = 0; i< size; ++i) {
			accumulator += v[i];
		}
	}
	std::cout << accumulator << std::endl;
}BENCHMARK(sin_opencl_short)->Apply(CustomStep);

//BENCHMARK_MAIN();
int main(int argc, char** argv) {
		::benchmark::Initialize(&argc, argv); 
		::benchmark::RunSpecifiedBenchmarks();
		std::cin.get();
		return 0;
}