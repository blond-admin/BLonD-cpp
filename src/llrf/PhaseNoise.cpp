/*
* @Author: Konstantinos Iliakis
* @Date:   2016-06-21 11:51:39
* @Last Modified by:   Konstantinos Iliakis
* @Last Modified time: 2016-06-21 14:38:41
* Methods to generate RF phase noise from noise spectrum and feedback noise
* amplitude as a function of bunch length**
* :Authors: **Helga Timko**
*/

#include <blond/constants.h>
#include <blond/llrf/PhaseNoise.h>
#include <blond/math_functions.h>
#include <blond/fft.h>
#include <random>
using namespace blond;

void PhaseNoise::spectrum_to_phase_noise(f_vector_t &t, f_vector_t &dphi,
        const f_vector_t &freq_array,
        const f_vector_t &ReS,
        PhaseNoise::transform_t transform)
{

    // Resolution in time domain
    const auto ReSLen = ReS.size();
    const auto f_max = freq_array.back();
    uint nt = 0;
    uint dt = 0;

    if (transform == transform_t::transform_none ||
            transform == transform_t::r) {
        nt = 2 * (ReSLen - 1);
        dt = 1 / (2 * f_max);
    } else if (transform == transform_t::c) {
        nt = ReSLen;
        dt = 1 / f_max;
    } else {
        std::cerr << "ERROR: The choise of Fourier transform for the\n"
                  << "RF noise generation could not be recognized.\n"
                  << "Use 'r' or 'c'.\n";
        exit(-1);
    }

    // Generate white noise in time domain
    f_vector_t r1(nt);
    f_vector_t r2(nt);

    if (fSeed1 < 0 || fSeed2 < 0) {
        f_vector_t random;
        // std::string home = util::GETENV("TEST_FILES");
        util::read_vector_from_file(random, TEST_FILES "/normal_distribution.dat");
        for (uint i = 0; i < r1.size(); i++) {
            r1[i] = random[(2 * i) % random.size()];
            r2[i] = random[(2 * i + 1) % random.size()];
        }
    } else {
        std::default_random_engine gen(fSeed1);
        std::uniform_real_distribution<> dist(0, 1);
        for (auto &v : r1) v = dist(gen);
        gen.seed(fSeed2);
        for (auto &v : r2) v = dist(gen);
    }

    if (transform == transform_t::transform_none ||
            transform == transform_t::r) {
        f_vector_t Gt(nt);

        auto factor = 2 * constant::pi;
        auto f1 = [factor](double x) { return std::cos(factor * x); };
        std::transform(r1.begin(), r1.end(), r1.begin(), f1);

        auto f2 = [](double x) { return std::sqrt(-2 * std::log(x)); };
        std::transform(r2.begin(), r2.end(), r2.begin(), f2);

        std::transform(r1.begin(), r1.end(), r2.begin(), Gt.begin(),
                       std::multiplies<double>());

        // FFT to frequency domain
        complex_vector_t Gf;
        fft::rfft(Gt, Gf);

        // Multiply by desired noise probability density
        factor = 2 * f_max;
        auto f3 = [factor](double x) { return std::sqrt(factor * x); };

        r1.resize(ReS.size());
        std::transform(ReS.begin(), ReS.end(), r1.begin(), f3);

        auto f4 = [](double r, complex_t z) { return r * z; };
        std::transform(r1.begin(), r1.end(), Gf.begin(), Gf.begin(), f4);

        // fft back to time domain to get final phase shift
        Gt.clear();
        fft::irfft(Gf, Gt);

        // Use only real part of the phase shift and normalize
        t.resize(nt);
        mymath::linspace(t.data(), 0, nt * dt, nt);

        auto f6 = [](complex_t z) { return z.real(); };

        dphi.resize(Gt.size());

        std::transform(Gt.begin(), Gt.end(), dphi.begin(), f6);

    } else if (transform == transform_t::c) {

        complex_vector_t Gt(nt);

        const complex_t factor = 2 * constant::pi * complex_t(0, 1);
        auto f1 = [factor](double x) { return std::exp(factor * x); };
        std::transform(r1.begin(), r1.end(), Gt.begin(), f1);

        auto f2 = [](double x) { return std::sqrt(-2 * std::log(x)); };
        std::transform(r2.begin(), r2.end(), r2.begin(), f2);

        auto f3 = [](complex_t a, double b) { return a * b; };
        std::transform(Gt.begin(), Gt.end(), r2.begin(), Gt.begin(), f3);

        // FFT to frequency domain
        complex_vector_t Gf;
        fft::fft(Gt, Gf);

        // Multiply by desired noise probability density

        auto factor2 = f_max;
        auto f4 = [factor2](double x) { return std::sqrt(factor2 * x); };

        r1.resize(ReS.size());
        std::transform(ReS.begin(), ReS.end(), r1.begin(), f4);

        auto f5 = [](double r, complex_t z) { return r * z; };

        std::transform(r1.begin(), r1.end(), Gf.begin(), Gf.begin(), f5);
        // fft back to time domain to get final phase shift
        Gt.clear();
        fft::ifft(Gf, Gt);

        t.resize(nt);
        mymath::linspace(t.data(), 0, nt * dt, nt);

        auto f6 = [](complex_t z) { return z.real(); };

        dphi.resize(Gt.size());

        std::transform(Gt.begin(), Gt.end(), dphi.begin(), f6);
    }
}

f_vector_t PhaseNoise::spectrum_generation(int k, int nt, double df, double ampl,
        f_vector_t freq,
        predistortion_t predistortion)
{
    int nmin = std::floor(fFMin * fFs[k] / df);
    int nmax = std::ceil(fFMax * fFs[k] / df);
    int nf = nt / 2 + 1;

    f_vector_t spectrum;

    if (fPredistortion == predistortion_t::exponential) {
        spectrum.resize(nmin, 0);

        auto v1 = mymath::arange<double>(0, nmax - nmin + 1);

        auto f1 = [ampl, nmax, nmin](double x) {
            return ampl * std::exp(std::log(100.0) * x / (nmax - nmin));
        };

        std::transform(v1.begin(), v1.end(), v1.begin(), f1);

        spectrum.insert(spectrum.end(), v1.begin(), v1.end());

        spectrum.resize(nf, 0);
    } else if (fPredistortion == predistortion_t::linear) {
        spectrum.resize(nmin, 0);

        f_vector_t v1(nmax - nmin + 1);
        mymath::linspace(v1.data(), 0, ampl, nmax - nmin + 1);

        spectrum.insert(spectrum.end(), v1.begin(), v1.end());

        spectrum.resize(nf, 0);
    } else if (fPredistortion == predistortion_t::hyperbolic) {
        spectrum.resize(nmin, 0);

        auto v1 = mymath::arange<double>(nmin, nmax + 1);

        auto f1 = [ampl, nmax, nmin](double x) {
            return ampl * 1.0 / (1 + 0.99 * (nmin - x) / (nmax - nmin));
        };

        std::transform(v1.begin(), v1.end(), v1.begin(), f1);
        spectrum.insert(spectrum.end(), v1.begin(), v1.end());

        spectrum.resize(nf, 0);
    } else if (fPredistortion == predistortion_t::weightfunction) {
        spectrum.resize(nmin, 0);

        // frequency relative to fs0
        f_vector_t frel(&freq[nmin], &freq[nmax + 1]);

        std::transform(frel.begin(), frel.end(), frel.begin(),
                       std::bind2nd(std::divides<double>(), fFs[k]));

        // truncate center freqs
        auto f1 = [](double x) { return x > 0.999 ? 0.999 : x; };
        std::transform(frel.begin(), frel.end(), frel.begin(), f1);
        // rms bunch length in rad corresponding to 1.2 ns
        auto sigma = 0.754;
        auto gamma = 0.577216;
        auto pi = constant::pi;
        auto f2 = [sigma, gamma, pi](double x) {
            auto sigma2 = sigma * sigma;
            return std::pow(4 * pi * x / sigma2, 2) *
                   std::exp(-16 * (1 - x) / sigma2) +
                   0.25 * pow(1 +
                              8 * x / sigma2 *
                              std::exp(-8 * (1 - x) / sigma2) *
                              (gamma + std::log(8 * (1 - x) / sigma2) +
                               8 * (1 - x) / sigma2),
                              2);
        };

        const auto factor = ampl / f2(frel[0]);
        for (uint i = 0; i < frel.size(); ++i) {
            frel[i] = factor * f2(frel[i]);
            // frel[i] *= factor;
        }

        spectrum.insert(spectrum.end(), frel.begin(), frel.end());

        spectrum.resize(nf, 0);

    } else {
        spectrum.resize(nmin, 0);
        spectrum.resize(nmax + 1, ampl);
        spectrum.resize(nf, 0);
    }

    return spectrum;

    // std::cout << "mean dphi : " << mymath::mean(fDphi.data(), fDphi.size())
    // << "\n";

    // std::cout << "mean dphi : " << mymath::mean(fDphi.data(), fDphi.size())
    // << "\n";
}

PSBPhaseNoiseInjection::PSBPhaseNoiseInjection(
    double delta_f, uint corr_time, double fmin, double fmax,
    double initial_amplitude, int seed1, int seed2,
    predistortion_t predistortion, rescale_ampl_t rescale_amplitude)
{
    auto RfP = Context::RfP;
    auto GP = Context::GP;

    fDeltaF = delta_f;
    fCorr = corr_time;
    fFMin = fmin;
    fFMax = fmax;
    fAi = initial_amplitude;
    fPredistortion = predistortion;
    fSeed1 = seed1;
    fSeed2 = seed2;
    fNTurns = GP->n_turns;
    fDphi.resize(fNTurns + 1, 0);
    fRescaleAmpl = rescale_amplitude;

    fFs.resize(RfP->omega_s0.size());
    std::transform(RfP->omega_s0.begin(), RfP->omega_s0.end(), fFs.begin(),
                   std::bind2nd(std::divides<double>(), 2 * constant::pi));
}

PSBPhaseNoiseInjection::~PSBPhaseNoiseInjection() { fft::destroy_plans(); }

void PSBPhaseNoiseInjection::generate()
{
    auto GP = Context::GP;

    for (uint i = 0; i < fNTurns / fCorr; ++i) {

        // Scale amplitude to keep area (phase noise amplitude) constant
        uint k = i * fCorr; // Current time step
        double f_max = GP->f_rev[k] / 2;
        double f_s0 = fFs[k];

        int n_points_pos_f_incl_zero = (int)(f_max / fDeltaF) + 2;
        int nt = 2 * (n_points_pos_f_incl_zero - 1);
        nt = mymath::next_regular(nt);
        n_points_pos_f_incl_zero = nt / 2 + 1;
        double df = f_max / (n_points_pos_f_incl_zero - 1);

        double ampl = 0.0;
        if (fRescaleAmpl == rescale_ampl_t::with_sync_freq)
            ampl = fAi * fFs[0] / f_s0;
        else if (fRescaleAmpl == rescale_ampl_t::no_scaling)
            ampl = fAi;
        else {
            std::cerr << "ERROR: Illegal Value\n"
                      << "fRescaleAmpl should be one of "
                      << "with_sync_freq or no_scaling\n";
        }

        f_vector_t freq(n_points_pos_f_incl_zero);
        mymath::linspace(freq.data(), 0, f_max, n_points_pos_f_incl_zero);

        auto spectrum =
            spectrum_generation(k, nt, df, ampl, freq, fPredistortion);

        f_vector_t noise_t, noise_dphi;
        spectrum_to_phase_noise(noise_t, noise_dphi, freq, spectrum);

        fSeed1 += 239;
        fSeed2 += 158;

        // Fill phase noise array

        const uint kmax =
            i < fNTurns / fCorr - 1 ? (i + 1) * fCorr : fNTurns + 1;

        std::copy(noise_dphi.begin(), noise_dphi.begin() + kmax - k,
                  fDphi.begin() + k);

        auto rms_noise =
            mymath::standard_deviation(noise_dphi.data(), noise_dphi.size());
        // std::cout << "RF noise for time step " << noise_t[1]
        //           << " s (iter " << i << ") has r.m.s phase "
        //           << rms_noise << " rad (" << rms_noise * 180 / constant::pi
        //           << " deg)\n";
    }
}

LHCFlatSpectrum::LHCFlatSpectrum(uint time_points, uint corr_time, double fmin,
                                 double fmax, double initial_amplitude, int seed1,
                                 int seed2, predistortion_t predistortion)
{
    auto RfP = Context::RfP;
    auto GP = Context::GP;

    fNt = time_points;
    fCorr = corr_time;
    fFMin = fmin;
    fFMax = fmax;
    fAi = initial_amplitude;
    fPredistortion = predistortion;
    fSeed1 = seed1;
    fSeed2 = seed2;
    fNTurns = GP->n_turns;
    fDphi.resize(fNTurns + 1, 0);

    if (fPredistortion != predistortion_t::predistortion_none) {
        // Overwrite frequencies
        fFMin = 0.8571;
        fFMax = 1.001;
    }

    if (fNt < 2 * fCorr) {
        std::cerr << "ERROR: Need more time point in LHCFlatSpectrum.\n";
        exit(-1);
    }

    // Synchrotron frequency array
    auto phis =  RfP->calc_phi_s(RfP, RfParameters::acc_sys_t::as_single);

    fFs.resize(phis.size());

    for (uint i = 0; i < fFs.size(); ++i) {
        fFs[i] =
            constant::c / GP->ring_circumference *
            std::sqrt(RfP->harmonic[RfP->section_index][i]
                      * RfP->voltage[RfP->section_index][i] *
                      std::abs(RfP->eta_0[i] * std::cos(phis[i])) /
                      (2 * constant::pi * RfP->energy[i]));
    }
}

LHCFlatSpectrum::~LHCFlatSpectrum() { fft::destroy_plans(); }

void LHCFlatSpectrum::generate()
{

    auto GP = Context::GP;

    for (uint i = 0; i < fNTurns / fCorr; ++i) {

        // Scale amplitude to keep area (phase noise amplitude) constant
        int k = i * fCorr; // Current time step
        auto ampl = fAi * fFs[0] / fFs[k];

        // Calculate the frequency step
        int nf = fNt / 2 + 1; // #points in frequency domain
        double df = GP->f_rev[k] / fNt;

        f_vector_t freq(nf);
        mymath::linspace(freq.data(), 0, nf * df, nf);

        auto spectrum =
            spectrum_generation(k, fNt, df, ampl, freq, fPredistortion);

        f_vector_t noise_t, noise_dphi;
        spectrum_to_phase_noise(noise_t, noise_dphi, freq, spectrum);

        // spectrum_to_phase_noise();
        fSeed1 += 239;
        fSeed2 += 158;

        // Fill phase noise array

        const uint kmax =
            i < fNTurns / fCorr - 1 ? (i + 1) * fCorr : fNTurns + 1;

        for (uint i = 0; i < kmax - k; ++i)
            fDphi[i + k] = noise_dphi[i];

        auto rms_noise =
            mymath::standard_deviation(noise_dphi.data(), noise_dphi.size());

        // std::cout << "RF noise for time step " << noise_t[1]
        //           << " s (iter " << i << ") has r.m.s phase "
        //           << rms_noise << " rad (" << rms_noise * 180 / constant::pi
        //           << " deg)\n";
    }
}
