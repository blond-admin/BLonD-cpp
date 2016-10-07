/*
* Slices.cpp
*
*  Created on: Mar 22, 2016
*      Author: kiliakis
*/
#include <blond/beams/Slices.h>
#include <blond/globals.h>
#include <blond/math_functions.h>
#include <blond/fft.h>
#include <blond/python.h>

Slices::Slices(uint n_slices, int n_sigma, ftype cut_left, ftype cut_right,
               cuts_unit_t cuts_unit, fit_t fit_option,
               bool direct_slicing)
{
    this->n_slices = n_slices;
    this->cut_left = cut_left;
    this->cut_right = cut_right;
    this->cuts_unit = cuts_unit;
    this->fit_option = fit_option;
    this->n_sigma = n_sigma;
    this->n_macroparticles.resize(n_slices, 0);
    this->fFloatMacroparticles.resize(n_slices, 0.0);
    this->edges.resize(n_slices + 1, 0.0);
    this->bin_centers.resize(n_slices, 0.0);
    set_cuts();

    if (fit_option == gaussian) {
        bl_gauss = 0.0;
        bp_gauss = 0.0;
    }

    if (direct_slicing) track();
}

Slices::~Slices() { fft::destroy_plans(); }

void Slices::set_cuts()
{
    auto Beam = Context::Beam;
    /*
    *Method to set the self.cut_left and self.cut_right properties. This is
    done as a pre-processing if the mode is set to 'const_space', for
    'const_charge' this is calculated each turn.*
    *The frame is defined by :math:`n\sigma_{RMS}` or manually by the user.
    If not, a default frame consisting of taking the whole bunch +5% of the
    maximum distance between two particles in the bunch will be taken
    in each side of the frame.*
    */

    if (cut_left == 0 && cut_right == 0) {
        if (n_sigma == 0) {
            sort_particles();
            cut_left =
                Beam->dt.front() - 0.05 * (Beam->dt.back() - Beam->dt.front());
            cut_right =
                Beam->dt.back() + 0.05 * (Beam->dt.back() - Beam->dt.front());
        } else {
            ftype mean_coords = mymath::mean(Beam->dt.data(), Beam->dt.size());
            ftype sigma_coords = mymath::standard_deviation(
                                     Beam->dt.data(), Beam->dt.size(), mean_coords);
            cut_left = mean_coords - n_sigma * sigma_coords / 2;
            cut_right = mean_coords + n_sigma * sigma_coords / 2;
        }
    } else {
        cut_left = convert_coordinates(cut_left, cuts_unit);
        cut_right = convert_coordinates(cut_right, cuts_unit);
    }

    mymath::linspace(edges.data(), cut_left, cut_right, n_slices + 1);
    for (uint i = 0; i < bin_centers.size(); ++i)
        bin_centers[i] = (edges[i + 1] + edges[i]) / 2;
}

void Slices::sort_particles()
{
    /*
    *Sort the particles with respect to their position.*
    */

    auto Beam = Context::Beam;
    struct API particle {
        ftype dE;
        ftype dt;
        int id;
        bool operator<(const particle &other) const
        {
            return dt < other.dt;
        }
    };

    std::vector<particle> particles(Beam->dE.size());
    for (uint i = 0; i < particles.size(); i++)
        particles[i] = {Beam->dE[i], Beam->dt[i], Beam->id[i]};

    std::sort(particles.begin(), particles.end());

    for (uint i = 0; i < particles.size(); i++) {
        Beam->dE[i] = particles[i].dE;
        Beam->dt[i] = particles[i].dt;
        Beam->id[i] = particles[i].id;
    }
}

ftype Slices::convert_coordinates(const ftype cut,
                                  const cuts_unit_t type)
{
    auto RfP = Context::RfP;
    /*
    *Method to convert a value from one input_unit_type to 's'.*
    */
    if (type == s)
        return cut;
    else  /*if (type == rad)*/
        return cut / RfP->omega_RF[RfP->idx][RfP->counter];

}

void Slices::track()
{
    slice_constant_space_histogram();
    if (fit_option == fit_t::gaussian)
        gaussian_fit();
}

void Slices::slice_constant_space_histogram()
{
    /*
    *Constant space slicing with the built-in numpy histogram function,
    with a constant frame. This gives the same profile as the
    slice_constant_space method, but no compute statistics possibilities
    (the index of the particles is needed).*
    *This method is faster than the classic slice_constant_space method
    for high number of particles (~1e6).*
    */
    auto Beam = Context::Beam;

    histogram(Beam->dt.data(), n_macroparticles.data(), cut_left, cut_right,
              n_slices, Beam->n_macroparticles);
}

void Slices::histogram(const ftype *__restrict input,
                       int *__restrict output,
                       const ftype cut_left,
                       const ftype cut_right,
                       const int n_slices,
                       const int n_macroparticles)
{

    const ftype inv_bin_width = n_slices / (cut_right - cut_left);
    // const ftype range = cut_right - cut_left;
    // histogram is faster with ints
    typedef int hist_t;
    hist_t *h;
    #pragma omp parallel
    {
        const int threads = omp_get_num_threads();
        const int id = omp_get_thread_num();
        int tile = (n_macroparticles + threads - 1) / threads;
        int start = id * tile;
        int end = std::min(start + tile, n_macroparticles);
        const int row = id * n_slices;

        #pragma omp single
        h = (hist_t *)calloc(threads * n_slices, sizeof(hist_t));

        int *h_row = &h[row];

        for (int i = start; i < end; ++i) {
            const ftype a = input[i];
            if (a < cut_left || a > cut_right) continue;
            const int ffbin = (int)((a - cut_left) * inv_bin_width);
            h_row[ffbin]++;
        }
        #pragma omp barrier

        tile = (n_slices + threads - 1) / threads;
        start = id * tile;
        end = std::min(start + tile, n_slices);

        for (int i = start; i < end; i++)
            output[i] = 0;
        // memset(&output[start], 0, (end-start) * sizeof(ftype));

        for (int i = 0; i < threads; ++i) {
            const int r = i * n_slices;
            for (int j = start; j < end; ++j) {
                output[j] += h[r + j];
            }
        }
    }

    if (h) free(h);
}

void Slices::track_cuts()
{
    /*
    *Track the slice frame (limits and slice position) as the mean of the
    bunch moves.
    Requires Beam statistics!
    Method to be refined!*
    */
    auto Beam = Context::Beam;

    ftype delta = Beam->mean_dt - 0.5 * (cut_left + cut_right);
    cut_left += delta;
    cut_right += delta;
    for (int i = 0; i < n_slices + 1; ++i) {
        edges[i] += delta;
    }
    for (int i = 0; i < n_slices; ++i) {
        bin_centers[i] += delta;
    }
}

void Slices::smooth_histogram(const ftype *__restrict input,
                              ftype *__restrict output,
                              const ftype cut_left,
                              const ftype cut_right,
                              const int n_slices,
                              const int n_macroparticles)
{

    int i;
    ftype a;
    ftype fbin;
    ftype ratioffbin;
    ftype ratiofffbin;
    ftype distToCenter;
    int ffbin = 0;
    int fffbin = 0;
    const ftype inv_bin_width = n_slices / (cut_right - cut_left);
    const ftype bin_width = (cut_right - cut_left) / n_slices;

    for (i = 0; i < n_slices; i++) {
        output[i] = 0;
    }

    for (i = 0; i < n_macroparticles; i++) {
        a = input[i];
        if ((a < (cut_left + bin_width * 0.5)) ||
                (a > (cut_right - bin_width * 0.5)))
            continue;
        fbin = (a - cut_left) * inv_bin_width;
        ffbin = (int)(fbin);
        distToCenter = fbin - (ftype)(ffbin);
        if (distToCenter > 0.5)
            fffbin = (int)(fbin + 1.0);
        ratioffbin = 1.5 - distToCenter;
        ratiofffbin = 1 - ratioffbin;
        if (distToCenter < 0.5)
            fffbin = (int)(fbin - 1.0);
        ratioffbin = 0.5 - distToCenter;
        ratiofffbin = 1 - ratioffbin;
        output[ffbin] = output[ffbin] + ratioffbin;
        output[fffbin] = output[fffbin] + ratiofffbin;
    }
}

void Slices::slice_constant_space_histogram_smooth()
{
    /*
    At the moment 4x slower than slice_constant_space_histogram but smoother.
    */
    auto Beam = Context::Beam;

    smooth_histogram(Beam->dt.data(), fFloatMacroparticles.data(), cut_left,
                     cut_right, n_slices, Beam->n_macroparticles);
}

void Slices::rms()
{
    /*
    * Computation of the RMS bunch length and position from the line density
    (bunch length = 4sigma).*
    */
    auto Beam = Context::Beam;

    f_vector_t lineDenNormalized(n_slices); // = new ftype[n_slices];
    f_vector_t array(n_slices);             // = new ftype[n_slices];

    const auto timeResolution = bin_centers[1] - bin_centers[0];
    const auto trap = mymath::trapezoid(n_macroparticles.data(),
                                        timeResolution,
                                        n_slices);

    for (int i = 0; i < n_slices; ++i)
        lineDenNormalized[i] = n_macroparticles[i] / trap;
    for (int i = 0; i < n_slices; ++i)
        array[i] = bin_centers[i] * lineDenNormalized[i];

    bp_rms = mymath::trapezoid(array.data(), timeResolution, n_slices);

    for (int i = 0; i <  n_slices; ++i)
        array[i] = (bin_centers[i] - bp_rms) * (bin_centers[i] - bp_rms) *
                   lineDenNormalized[i];

    auto temp = mymath::trapezoid(array.data(), timeResolution, n_slices);
    bl_rms = 4 * std::sqrt(temp);
}


void Slices::fwhm(const ftype shift)
{

    /*
    * Computation of the bunch length and position from the FWHM
    assuming Gaussian line density.*
    */
    int max = *std::max_element(n_macroparticles.begin(),
                                n_macroparticles.end());
    ftype half_max = shift + 0.5 * (max - shift);
    ftype timeResolution = bin_centers[1] - bin_centers[0];

    // First aproximation for the half maximum values
    int taux1, taux2;

    int i = 0;
    while (i < n_slices && n_macroparticles[i] < half_max) i++;
    taux1 = i;
    // prev1 is one before taux1, if taux1 is 0, then prev1 is the last
    int prev1 = taux1 > 0 ? taux1 - 1 : n_slices - 1;

    i = n_slices - 1;
    while (i >= 0 && n_macroparticles[i] < half_max) i--;
    taux2 = i;

    // dprintf("taux1, taux2 = %d, %d\n", taux1, taux2);
    ftype t1, t2;

    if (taux1 < n_slices && taux2 < n_slices - 1 && taux2 >= 0) {
        t1 = bin_centers[taux1] -
             (n_macroparticles[taux1] - half_max) /
             (n_macroparticles[taux1] - n_macroparticles[prev1]) *
             timeResolution;

        t2 = bin_centers[taux2] +
             (n_macroparticles[taux2] - half_max) /
             (n_macroparticles[taux2] - n_macroparticles[taux2 + 1]) *
             timeResolution;
        bl_fwhm = 4 * (t2 - t1) / cfwhm;
        bp_fwhm = (t1 + t2) / 2;
    } else {
        bl_fwhm = nan("");
        bp_fwhm = nan("");
    }
}

// ftype Slices::fast_fwhm()
// {

//     /*
//     * Computation of the bunch length and position from the FWHM
//     assuming Gaussian line density.*
//     height = np.max(self.n_macroparticles)
//     index = np.where(self.n_macroparticles > height/2.)[0]
//     return cfwhm*(Slices.bin_centers[index[-1]] - Slices.bin_centers[index[0]])
//     */
//     auto Beam = Context::Beam;

//     uint max_i =
//         mymath::max(n_macroparticles.data(), Beam->n_macroparticles, 1);
//     ftype half_max = 0.5 * n_macroparticles[max_i];

//     int i = 0;
//     while (n_macroparticles[i] < half_max && i < (int)n_slices)
//         i++;
//     int taux1 = i;
//     i = n_slices - 1;
//     while (i >= 0 && n_macroparticles[i] < half_max)
//         i--;
//     int taux2 = i;
//     // update bp
//     return cfwhm * (bin_centers[taux2] - bin_centers[taux1]);
// }

void Slices::fwhm_multibunch() {}

void Slices::beam_spectrum_generation(uint n, bool onlyRFFT)
{

    fBeamSpectrumFreq = fft::rfftfreq(n, bin_centers[1] - bin_centers[0]);

    if (onlyRFFT == false) {
        f_vector_t v(n_macroparticles.begin(), n_macroparticles.end());
        fft::rfft(v, fBeamSpectrum, n, Context::n_threads);
    }
}

void Slices::beam_profile_derivative(f_vector_t &x,
                                     f_vector_t &derivative,
                                     std::string mode)
{
    /*
    The input is one of the two available methods for differentiating
    a function. The two outputs are the coordinate step and the discrete
    derivative of the Beam profile respectively.
    */

    x = bin_centers;
    const auto dist_centers = x[1] - x[0];
    if (mode == "filter1d") {

        f_vector_t temp(n_macroparticles.begin(), n_macroparticles.end());
        derivative = gaussian_filter1d(temp, 1, 1, "wrap");
        for (auto &d : derivative) d /= dist_centers;

    } else if (mode == "gradient") {

        f_vector_t temp(n_macroparticles.begin(), n_macroparticles.end());
        derivative = gradient(temp, dist_centers);

    } else if (mode == "diff") {

        derivative.resize(n_macroparticles.size() - 1);
        for (uint i = 1; i < n_macroparticles.size(); i++)
            derivative[i - 1] = (n_macroparticles[i] - n_macroparticles[i - 1])
                                / dist_centers;
        f_vector_t diffCenters(x.begin(), x.end() - 1);
        for (auto &d : diffCenters) d += dist_centers / 2;
        f_vector_t res;
        mymath::lin_interp(x, diffCenters, derivative,
                           res, derivative.front(), derivative.back());
        derivative = res;
    } else {
        std::cerr << "Optioin for derivative is not recognized.\n";
    }

}


f_vector_t Slices::gaussian_filter1d(f_vector_t &x, int sigma,
                                     int order, std::string mode)
{
    python::import();
    auto pFunc = python::import("utilities", "gaussian_filter1d");
    auto pX = python::convert_double_array(x.data(), x.size());
    auto pSigma = python::convert_int(sigma);
    auto pOrder = python::convert_int(order);
    auto pMode = python::convert_string(mode);
    auto ret = PyObject_CallFunctionObjArgs(pFunc, pX, pSigma, pOrder, pMode,
                                            NULL);
    assert(ret);
    auto npArray = (PyArrayObject *)(ret);
    int len = PyArray_SHAPE(npArray)[0];
    ftype *res = (ftype *) PyArray_DATA(npArray);

    return f_vector_t(&res[0], &res[len]);
}


f_vector_t Slices::gradient(f_vector_t &x, ftype dist)
{
    python::import();
    auto pFunc = python::import("utilities", "gradient");
    auto pX = python::convert_double_array(x.data(), x.size());
    auto pDist = python::convert_double(dist);
    auto ret = PyObject_CallFunctionObjArgs(pFunc, pX, pDist, NULL);
    assert(ret);

    auto npArray = (PyArrayObject *)(ret);
    int len = PyArray_SHAPE(npArray)[0];
    ftype *res = (ftype *) PyArray_DATA(npArray);
    return f_vector_t(&res[0], &res[len]);
}

// NOTE: if you specify transfer_function_plot == "true" then
// b, a = cheby2(nCoefficients, ...)
// this function must not be called with transfer_function_plot == "true" !!
void Slices::beam_profile_filter_chebyshev(std::map<std::string, std::string>
        filter_option, int &nCoefficients, f_vector_t &b, f_vector_t &a)
{
    python::import();
    auto pFunc = python::import("utilities", "beam_profile_filter_chebyshev");
    auto pNMacroparticles = python::convert_int_array(n_macroparticles.data(),
                            n_macroparticles.size());
    auto pResolution = python::convert_double(bin_centers[1] - bin_centers[0]);
    auto pNSlices = python::convert_int(n_slices);

    if (filter_option.find("transfer_function_plot") != filter_option.end()
            && filter_option["transfer_function_plot"] == "true") {
        std::cerr << "[beam_profile_filter_chebyshev] A complex vector must\n"
                  << "be pased as the last argument in order to call this\n"
                  << "function with transfer_function_plot ==true\n";
        exit(-1);
    }

    auto pFilterOption = python::convert_dictionary(filter_option);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pNMacroparticles,
                                            pResolution, pNSlices,
                                            pFilterOption, NULL);

    if (!ret) {
        std::cerr << "[beam_profile_filter_chebyshev] An error occured "
                  << "while running a python function\n";
        exit(-1);
    }

    auto npArray = (PyArrayObject *) ret;
    int len = PyArray_SHAPE(npArray)[0];
    ftype *res = (ftype *) PyArray_DATA(npArray);
    nCoefficients = res[0];

    b = f_vector_t(&res[1], &res[1 + 1 + nCoefficients]);
    a = f_vector_t(&res[2 + nCoefficients], &res[3 + 2 * nCoefficients]);
}

// NOTE
// this function must be called with transfer_function_plot == "true" !!
// if not then an error will happen
void Slices::beam_profile_filter_chebyshev(std::map<std::string, std::string>
        filter_option, int &nCoefficients,
        f_vector_t &transferFreq, complex_vector_t &transferGain)
{
    python::import();
    auto pFunc = python::import("utilities", "beam_profile_filter_chebyshev");
    auto pNMacroparticles = python::convert_int_array(n_macroparticles.data(),
                            n_macroparticles.size());
    auto pResolution = python::convert_double(bin_centers[1] - bin_centers[0]);
    auto pNSlices = python::convert_int(n_slices);

    if (filter_option.find("transfer_function_plot") == filter_option.end()
            || filter_option["transfer_function_plot"] != "true") {
        std::cerr << "[beam_profile_filter_chebyshev] A double vector must\n"
                  << "be pased as the last argument in order to call this\n"
                  << "function without transfer_function_plot == true\n";
        exit(-1);
    }
    auto pFilterOption = python::convert_dictionary(filter_option);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pNMacroparticles,
                                            pResolution, pNSlices,
                                            pFilterOption, NULL);

    if (!ret) {
        std::cerr << "[beam_profile_filter_chebyshev] An error occured "
                  << "while running a python function\n";
        exit(-1);
    }

    auto npArray = (PyArrayObject *) ret;
    int len = PyArray_SHAPE(npArray)[0];
    ftype *res = (ftype *) PyArray_DATA(npArray);

    if (len != 1 + 3 * n_slices) {
        std::cerr << "[beam_profile_filter_chebyshev] Unexpected result length\n";
        exit(-1);
    }

    nCoefficients = res[0];
    transferFreq = f_vector_t(&res[1], &res[1 + n_slices]);

    for (int i = 1 + n_slices; i < 1 + 2 * n_slices; i++) {
        transferGain.push_back({res[i], res[i + n_slices]});
    }
}


void Slices::gaussian_fit()
{
    auto Beam = Context::Beam;
    f_vector_t p0;
    double max = *std::max_element(n_macroparticles.begin(),
                                   n_macroparticles.end());
    // std::cout.precision(12);
    // std::cout << "max: " << max << "\n";
    if (bl_gauss == 0 && bp_gauss == 0) {
        p0 = {max,
              mymath::mean(Beam->dt.data(), Beam->dt.size()),
              mymath::standard_deviation(Beam->dt.data(), Beam->dt.size())
             };
    } else {
        p0 = {max,
              bp_gauss,
              bl_gauss / 4
             };
    }
    // util::dump(p0, "p0\n");

    python::import();
    auto pFunc = python::import("utilities", "curve_fit");
    auto pBinCenters = python::convert_double_array(bin_centers.data(),
                       bin_centers.size());
    auto pNMacroparticles = python::convert_int_array(n_macroparticles.data(),
                            n_macroparticles.size());
    auto pP0 = python::convert_double_array(p0.data(), p0.size());
    auto ret = PyObject_CallFunctionObjArgs(pFunc, pBinCenters,
                                            pNMacroparticles, pP0, NULL);
    if (!ret) {
        std::cerr << "[gaussian_fit] An error occured while running"
                  << "a python function\n";
        exit(-1);
    }

    auto npArray = (PyArrayObject *) ret;
    int len = PyArray_SHAPE(npArray)[0];
    ftype *pfit_gauss = (ftype *) PyArray_DATA(npArray);
    // std::cout << "pfit_gauss[1]: " <<pfit_gauss[1] << "\n";
    // std::cout << "pfit_gauss[2]: " <<pfit_gauss[2] << "\n";
    bl_gauss = 4 * std::abs(pfit_gauss[2]);
    bp_gauss = std::abs(pfit_gauss[1]);

}
