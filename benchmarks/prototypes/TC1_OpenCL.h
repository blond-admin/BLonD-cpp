/*
OpenCL performance implies that data stays in the device memory as long as possible.
Thus all algorithm-requiered data shall be provided to an external user only on demand. That implies 
protection of all vectors such as beam_dE (a user needs to ask for beam_dE and tall us if he changed it, 
each such operation is data transfer from PC to device.
In this sample we read all data at the constructor stage.
To get data out use destructor or write_data
We warm up all Device functions in constructor (text->binary->initiated, ready to run kernel)
We do not use BOOST_COMPUTE_CLOSURE+transform because:
a) BOOST_COMPUTE_STRINGIZE_SOURCE allows __constant/__local/read_only explicit setup
b) we can be sure when the Device functions are compiled (and that they are compiled only once)
Developers note: printf works in OpenCL kernels, reading OpenCL error codes is useful (see https://streamcomputing.eu/blog/2013-04-28/opencl-error-codes/)
*/

#ifdef WITH_OPENCL

//Set OpenCL compilation exceptions in debug mode
#ifndef NDEBUG
#define BOOST_COMPUTE_DEBUG_KERNEL_COMPILATION
#endif

#include <boost/compute.hpp>

//Set OpenCL extensions for float64 and atomic_inc
//#ifdef cl_khr_fp64
//#pragma OPENCL EXTENSION cl_khr_fp64 : enable
//#elif defined(cl_amd_fp64)
//#pragma OPENCL EXTENSION cl_amd_fp64 : enable
//#else
//#endif
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable


struct TC1_OpenCL {
	TC1_OpenCL() {
		using namespace boost::compute;
		//Device kernels source code
		const std::string source_kick = BOOST_COMPUTE_STRINGIZE_SOURCE(
			__kernel void kick(__global double* beam_dt,
								__global double* beam_dE, 
								double current_voltage,
								double current_omega_RF,
								double current_phi,
								double acc_kick_gpu
				) {
			int id = get_global_id(0);
			beam_dE[id] += acc_kick_gpu + current_voltage * sin(current_omega_RF * beam_dt[id] + current_phi);
		});

			const std::string source_drift = BOOST_COMPUTE_STRINGIZE_SOURCE(
			__kernel void drift(__global double* beam_dE,
								__global double* beam_dt,
								double T_x_coeff
				) {
				int id = get_global_id(0);
				beam_dt[id] += T_x_coeff * beam_dE[id];
		});

			const std::string source_histogram = BOOST_COMPUTE_STRINGIZE_SOURCE(
			__kernel void histogram(__global double* input,
				__global int* output,
				int cut_left,
				int cut_right
				) {
			int id = get_global_id(0);
			int key = input[id];
			if (key < cut_left || key > cut_right) {
				return;
			}
			atomic_inc(&output[key]);
		});

		//Setup device
		gpu = system::default_device();
		context= boost::compute::context(gpu);
		queue = command_queue(
			context, gpu
		);

		//data for Device kernels 
		auto Beam = Context::Beam;
		auto Slice = Context::Slice;

		n_macroparticles = vector<int>(Slice->n_macroparticles.size(), context);
		dE = vector<ftype>(Beam->dE.size(), context);
		dt = vector<ftype>(Beam->dt.size(), context);
		move_array(Slice->n_macroparticles, n_macroparticles);
		move_array(Beam->dE, dE);
		move_array(Beam->dt, dt);
		queue.finish();

		//Install programms
		build_kernel("kick", source_kick, kick_kernel);
		build_kernel("drift", source_drift, drift_kernel);
		build_kernel("histogram", source_histogram, histogram_kernel);

		queue.finish();

		//WarmUp CPU test data
		this->acceleration_kick.resize(Context::GP->n_turns + 1);
		for (uint i = 0; i < Context::RfP->E_increment.size(); ++i) {
			acceleration_kick[i] = -Context::RfP->E_increment[i];
		}
	}
	~TC1_OpenCL() {
		queue.flush();
		queue.finish();
	}
	//Saves data from device vectors onto host vectors
	void write_data() {
		auto GP = Context::GP;
		auto Beam = Context::Beam;
		move_array(n_macroparticles, Context::Slice->n_macroparticles);
		move_array(dE, Beam->dE);
		move_array(dt, Beam->dt);
		queue.finish();
	}

	void track() {
		RingAndRfSection_track();
		Slice_track();
		Context::RfP->counter++;
	}

	void RingAndRfSection_track() {
		kick(Context::RfP->counter);
		drift(Context::RfP->counter + 1);
	}

	void Slice_track() {
		auto Slice = Context::Slice;
		histogram(Slice->cut_left, Slice->cut_right);
	}
private:
	//GPU related members
	boost::compute::context context;
	boost::compute::command_queue queue;
	boost::compute::device gpu;
	//Kernels
	boost::compute::kernel kick_kernel;
	boost::compute::kernel drift_kernel;
	boost::compute::kernel histogram_kernel;
	//Input/Output
	boost::compute::vector<ftype> dE;
	boost::compute::vector<ftype> dt;
	boost::compute::vector<int> n_macroparticles;
	f_vector_t acceleration_kick;

	//Wrapper
	void kick(const uint index) {
		auto RfP = Context::RfP;
		kick(RfP->n_rf, RfP->voltage, RfP->omega_RF, RfP->phi_RF, index);
	}

	//Host implementation
	void kick(const int n_rf, 
		const f_vector_2d_t & __restrict voltage,
		const f_vector_2d_t  & __restrict omega_RF,
		const f_vector_2d_t & __restrict phi_RF, 
		const uint index) {
		using namespace boost::compute;
		auto acc_kick_gpu = acceleration_kick[index] / static_cast<ftype>(n_rf);
		
		for (int j = 0; j < n_rf; ++j) {
			auto current_voltage = voltage[j][index];
			auto current_omega_RF = omega_RF[j][index];
			auto current_phi = phi_RF[j][index];
			
			kick(dt, dE, current_voltage, current_omega_RF, current_phi, acc_kick_gpu);
		}
	}

	//wrapper
	inline void drift(const uint index) {
		auto GP = Context::GP;
		auto RfP = Context::RfP;
		auto Beam = Context::Beam;

		drift(GP->t_rev[index],
			RfP->length_ratio, 
			RfP->eta_0(index),
			RfP->beta(index),
			RfP->energy(index));
	}

	//Host implementation
	inline void drift(const ftype T0,
		const ftype length_ratio,
		const ftype eta_zero,
		const ftype beta, 
		const ftype energy) {
		const ftype T = T0 * length_ratio;
		const ftype T_x_coeff = T * eta_zero / (beta * beta * energy);

		drift(dE, dt, T_x_coeff);
	}

	//Host implementation
	inline void histogram(const ftype cut_left, const ftype cut_right) {
		histogram(dt, n_macroparticles, cut_left, cut_right);
	}

	//Device kernel call
	void kick(const boost::compute::vector<ftype> & beam_dt,
		boost::compute::vector<ftype> & beam_dE,
		double current_voltage,
		double current_omega_RF,
		double current_phi,
		double acc_kick_gpu) {

		using namespace boost::compute;
		kick_kernel.set_arg(0, beam_dt.get_buffer());
		kick_kernel.set_arg(1, beam_dE.get_buffer());
		kick_kernel.set_arg(2, current_voltage);
		kick_kernel.set_arg(3, current_omega_RF);
		kick_kernel.set_arg(4, current_phi);
		kick_kernel.set_arg(5, acc_kick_gpu);

		start_kernel(kick_kernel, beam_dE.size());
	}

	//Device kernel call
	void drift(const boost::compute::vector<ftype> & beam_dE,
		boost::compute::vector<ftype> & beam_dt,
		double T_x_coeff) {

		using namespace boost::compute;
		drift_kernel.set_arg(0, beam_dE.get_buffer());
		drift_kernel.set_arg(1, beam_dt.get_buffer());
		drift_kernel.set_arg(2, T_x_coeff);

		start_kernel(drift_kernel, beam_dE.size());
	}

	//Device kernel call
	void histogram(const boost::compute::vector<ftype> & input,
		boost::compute::vector<int> & output,
		int cut_left,
		int cut_right
	) {
		using namespace boost::compute;
		histogram_kernel.set_arg(0, input.get_buffer());
		histogram_kernel.set_arg(1, output.get_buffer());
		histogram_kernel.set_arg(2, cut_left);
		histogram_kernel.set_arg(3, cut_right);
		start_kernel(histogram_kernel, input.size());
	}

	//Utilety functions to read/write kernels, transfer 1d and 2d vectors
	void build_kernel(const std::string & function_name, const std::string & source, boost::compute::kernel & kernel) {
		using namespace boost::compute;
		program program =
			program::create_with_source(source, context);
		program.build();
		kernel = boost::compute::kernel(program, function_name);
	}

	void start_kernel(const boost::compute::kernel & kernel, int size) {
		const size_t origin = 0;
		int wg_size = gpu.max_work_group_size();
		auto rem = size % wg_size;
		size = wg_size * (size / wg_size);
		if (size  > 0) {
			queue.enqueue_1d_range_kernel(kernel, origin, size, wg_size);
		}
		if (rem > 0) {
			queue.enqueue_1d_range_kernel(kernel, size, size + rem, 0);
		}

		queue.finish();
	}

	template<class T1, class T2>
	void move_array(const T1 & from, T2 & to) {
		to.resize(from.size());
		boost::compute::copy(from.begin(), from.end(), to.begin(), queue);
	}

	template<class T1, class T2>
	void move_2d_from_host(const T1 & from, T2 & to) {
		auto counter = 0;
		int rowLength = from.begin()->size();
		to.reserve(from.size()*rowLength);
		for (auto & row : from) {
			boost::compute::copy(row.begin(), row.end(), to.begin() + counter*rowLength, queue);
			++counter;
		}
	}

	template<class T1, class T2>
	void move_2d_to_host(const T1 & from, T2 & to) {
		int rowLength = to.begin()->size();
		for (auto i = 0; i < from.size(); ++i) {
			boost::compute::copy(from.begin() + i*rowLength, from.end() + (i + 1)*rowLength, (to.begin() + i)->begin(), queue);
		}
	}
};
#endif