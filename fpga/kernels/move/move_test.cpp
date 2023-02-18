#include "move.hpp"

#include <cmath>
#include <iostream>
#include <iomanip>

#define TEST_TOLERANCE 1e-14

template<typename T>
bool is_near(T a, T b, T eps) {
  return std::abs(a - b) < eps;
}

constexpr int dims[] = {4, 4, 4};
constexpr int num = dims[0] * dims[1] * dims[2];
constexpr int size = 3 * num;
constexpr real_t box[] = {1, 1, 1};
constexpr real_t xc[] = {box[0]/2, box[1]/2, box[2]/2};
constexpr real_t dx[] = {
    box[0]/(dims[0]-1),
    box[1]/(dims[1]-1),
    box[2]/(dims[2]-1)
};

constexpr real_t dtc = 0.1/3;
constexpr real_t vel_mag = 1;

constexpr real_t mass = 1;
constexpr real_t temperature = 10;
constexpr real_t boltzmann = 1;

int main() {
    
  real_t v[size], a[size];
  real_t x[size], res[size];

  auto dx2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
  auto dx_mag = std::sqrt( dx2 );
  auto dt = dtc * dx_mag / vel_mag;

  auto dt_sqr = dt*dt;

  for (int i = 0, pos = 0; i < dims[0]; ++i)
    for (int j = 0; j < dims[1]; ++j)
        for (int k = 0; k < dims[2]; ++k, ++pos)
        {
            real_t xi[] = {dx[0]*i, dx[1]*j, dx[2]*k};

            real_t delta[] = {
                xi[0] - xc[0],
                xi[1] - xc[1],
                xi[2] - xc[2],
            };

            auto d2 = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
            auto dmag = std::sqrt( d2 );
            auto vfact = vel_mag / dmag;
            
            real_t vi[] = {delta[0]*vfact, delta[1]*vfact, delta[2]*vfact};
            real_t ai[] = {vi[0]*vel_mag, vi[1]*vel_mag, vi[2]*vel_mag};

            for (int d=0; d<3; ++d) {
                x  [3*pos + d] = xi[d];
                v  [3*pos + d] = vi[d];
                a  [3*pos + d] = ai[d];
                auto xnew = xi[d] + vi[d]*dt + 0.5*dt_sqr*ai[d];

                auto len = box[d];
                if (xnew > len)
                    xnew -= len;
                else if (xnew < 0)
                    xnew += len;
                res[3*pos + d] = xnew;
            }
            //std::cout << "moving pos=" << pos << " : ";
            //for (int d=0; d<3; ++d) std::cout << std::setw(8) << std::setprecision(4) << x[3*pos + d] << " ";
            //std::cout << " to ";
            //for (int d=0; d<3; ++d) std::cout << std::setw(8) << std::setprecision(4) << res[3*pos + d] << " ";
            //std::cout << std::endl;
        }

  move_kernel(x, v, a, dt, box, mass, temperature, boltzmann, num);

  for (int i = 0, pos = 0; i < dims[0]; ++i)
    for (int j = 0; j < dims[1]; ++j)
        for (int k = 0; k < dims[2]; ++k, ++pos)
        {
            for (int d=0; d<3; ++d) {
                if (!is_near(res[3*pos+d], x[3*pos+d], TEST_TOLERANCE)) {
                    std::cout << "*** Test failed. ***" << std::endl;
                    std::cout << "Error is " << std::abs(res[3*pos+d] - x[3*pos+d]) << std::endl;
                    std::cout << std::fixed << std::setprecision(4);
                    std::cout << "pos=" << pos << " : ";
                    for (int dd=0; dd<3; ++dd)
                        std::cout << res[3*pos+dd] << " ";
                    std::cout << " versus ";
                    for (int dd=0; dd<3; ++dd)
                        std::cout << x[3*pos+dd] << " ";
                    return EXIT_FAILURE;
                }
            }
        }

  std::cout << "Test passed." << std::endl;
  return EXIT_SUCCESS;
}
