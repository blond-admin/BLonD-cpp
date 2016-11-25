/*
* Music.cpp
*
*  Created on: Nov 25, 2016
*      Author: kiliakis
*/
#include <blond/impedances/Music.h>
#include <blond/vector_math.h>
#include <blond/constants.h>
#include <blond/math_functions.h>
#ifdef PARALLEL
#include <parallel/algorithm>
#else
#include <algorithm>
#endif

using namespace std;

Music::Music(Beams *beam, const f_vector_t &resonator,
             int n_macroparticles, long int n_particles)
{
    Beam = beam;
    R_S = resonator[0];
    omega_R = resonator[1];
    Q = resonator[2];
    this->n_macroparticles = n_macroparticles;
    this->n_particles = n_particles;
    alpha = omega_R / (2 * Q);
    omega_bar = sqrt(omega_R * omega_R - alpha * alpha);
    constant = -constant::e * R_S * omega_R * n_particles
               / (n_macroparticles * Q);
    induced_voltage.resize(Beam->dt.size(), 0);
    induced_voltage[0] = constant / 2;
    coeff1 = -alpha / omega_bar;
    coeff2 = -R_S * omega_R / (Q * omega_bar);
    coeff3 = omega_R * Q / (R_S * omega_bar);
    coeff4 = alpha / omega_bar;
}

void Music::track()
{

    struct particle {
        double de;
        double dt;
        bool operator<(const particle &o) const
        {
            return dt < o.dt;
        }
    };

//     std::chrono::time_point<std::chrono::high_resolution_clock>start;
//     std::chrono::duration<double> duration(0.0);
//     start = std::chrono::system_clock::now();

    vector<particle> particles;
    particles.resize(n_macroparticles);
    for (int i = 0; i < n_macroparticles; i++)
        particles[i] = {Beam->dE[i], Beam->dt[i]};

#ifdef PARALLEL
    __gnu_parallel::sort(ALL(particles));
#else
    std::sort(ALL(particles));
#endif

    for (int i = 0; i < n_macroparticles; i++) {
        Beam->dE[i] = particles[i].de;
        Beam->dt[i] = particles[i].dt;
    }

//     duration = std::chrono::system_clock::now() - start;
//     std::cout << "sorting time: " << duration.count() << '\n';

//     start = std::chrono::system_clock::now();

    Beam->dE[0] += induced_voltage[0];
    double input_first_component = 1;
    double input_second_component = 0;
    for (int i = 0; i < n_macroparticles - 1; i++) {
        const double time_difference = Beam->dt[i + 1] - Beam->dt[i];
        const double exp_term = mymath::fast_exp(-alpha * time_difference);
        const double cos_term = mymath::fast_cos(omega_bar * time_difference);
        const double sin_term = mymath::fast_sin(omega_bar * time_difference);

        const double product_first_component =
            exp_term * ((cos_term + coeff1 * sin_term)
                        * input_first_component + coeff2 * sin_term
                        * input_second_component);

        const double product_second_component =
            exp_term * (coeff3 * sin_term * input_first_component
                        + (cos_term + coeff4 * sin_term)
                        * input_second_component);

        induced_voltage[i + 1] = constant * (0.5 + product_first_component);
        Beam->dE[i + 1] += induced_voltage[i + 1];
        input_first_component = product_first_component + 1;
        input_second_component = product_second_component;



//     duration = std::chrono::system_clock::now() - start;
//     std::cout << "tracking time: " << duration.count() << '\n';


    }
}