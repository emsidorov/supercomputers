#include <iostream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <mpi.h>
#include <omp.h>


class Grid {
public:
    Grid(int Nx, int Ny, int Nz, double dx, double dy, double dz)
        : Nx(Nx), Ny(Ny), Nz(Nz), dx(dx), dy(dy), dz(dz) {
        array = new double[Nx * Ny * Nz]();
        for (int i = 0; i < Nx * Ny * Nz; ++i) {
            array[i] = 0;
        }
    }

    double get(int i, int j, int k) const {
        return array[index(i, j, k)];
    }

    void set(int i, int j, int k, double value) {
        array[index(i, j, k)] = value;
    }

    double* get_pointer(int i, int j, int k) {
        return &array[index(i, j, k)];
    }

    double get_diff(int i, int j, int k) {
        double diff_x, diff_y, diff_z;
        diff_x = (this->get(i - 1, j, k) - 2 * this->get(i, j, k) + this->get(i + 1, j, k)) / (dx * dx);
        diff_y = (this->get(i, j - 1, k) - 2 * this->get(i, j, k) + this->get(i, j + 1, k)) / (dy * dy);
        diff_z = (this->get(i, j, k - 1) - 2 * this->get(i, j, k) + this->get(i, j, k + 1)) / (dz * dz);
        return diff_x + diff_y + diff_z;
    }

    int size_x() const { return Nx; }
    int size_y() const { return Ny; }
    int size_z() const { return Nz; }

    void swap(Grid& other) {
        std::swap(this->array, other.array);
    }

    void print() {
        std::cout << "----------Grid Start----------\n";
        for (int i = 0; i < Nx; ++i) {
            for(int j = 0; j < Ny; ++j) {
                for(int k = 0; k < Nz; ++k) {
                    std::cout << i << " " << j << " " << k << " " << this->get(i, j, k) << std::endl;
                }
            }
        }
        std::cout << "---------Grid Finish----------\n";
    }

    void save(const std::string& filename) {
        std::ofstream file(filename.c_str());
        if (!file.is_open()) {
            throw std::runtime_error("Unable to open file for writing.");
        }

        file << "Nx: " << Nx << " Ny: " << Ny << " Nz: " << Nz << "\n";
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {
                for (int k = 0; k < Nz; ++k) {
                    file << this->get(i, j, k) << ",";
                }
            }
        }

        file.close();
    }

private:
    int Nx, Ny, Nz;
    double *array;
    double dx, dy, dz;

    int index(int i, int j, int k) const {
        return i * (Ny * Nz) + j * Nz + k;
    }
};


class Function {
public:
    Function(int grid_size, double dx, double dy, double dz, double Lx, double Ly, double Lz)
    : grid_size(grid_size + 1), dx(dx), dy(dy), dz(dz), Lx(Lx), Ly(Ly), Lz(Lz),
      a_t(0.5 * std::sqrt(4.0 / (Lx * Lx) + 1.0 / (Ly * Ly) + 1.0 / (Lz * Lz))) {}

    double u(double x, double y, double z, double t) const {
        return std::sin(x * 2.0 * M_PI / Lx) * std::sin(y * M_PI / Ly) * std::sin(z * M_PI / Lz) * std::cos(a_t * t + 2 * M_PI);
    }

    double u_from_idx(int i, int k, int j, double t) const {
        double x = i * dx;
        double y = j * dy;
        double z = k * dz;
        return this->u(x, y, z, t);
    }

    void fill_grid(Grid& grid, double t) {
        for (int i = 0; i < grid_size; ++i) {
            for (int j = 0; j < grid_size; ++j) {
                for (int k = 0; k < grid_size; ++k) {
                    double x = i * dx;
                    double y = j * dy;
                    double z = k * dz;
                    grid.set(i, j, k, u(x, y, z, t));
                }
            }
        }
    }

    void fill_block(Grid& grid, int* left, int* size, double t) {
        for (int i = 0; i < size[0]; ++i) {
            for (int j = 0; j < size[1]; ++j) {
                for (int k = 0; k < size[2]; ++k) {
                    double x = (left[0] + i - 1) * dx;
                    double y = (left[1] + j - 1) * dy;
                    double z = (left[2] + k - 1) * dz;
                    grid.set(i, j, k, u(x, y, z, t));
                }
            }
        }
    }

private:
    double Lx, Ly, Lz;
    int grid_size;
    double a_t;
    double dx, dy, dz;
};


double compute_error(const Grid& grid1, const Grid& grid2) {
    if (grid1.size_x() != grid2.size_x() || 
        grid1.size_y() != grid2.size_y() || 
        grid1.size_z() != grid2.size_z()) { return -1; }

    int size = grid1.size_x() * grid1.size_y() * grid1.size_z();
    double sum_diff = 0.0;

    for (int i = 0; i < grid1.size_x(); ++i) {
        for (int j = 0; j < grid1.size_y(); ++j) {
            for (int k = 0; k < grid1.size_z(); ++k) {
                double diff = grid1.get(i, j, k) - grid2.get(i, j, k);
                sum_diff += std::abs(diff);
            }
        }
    }

    return sum_diff / size;
}


int main(int argc, char *argv[]) {

// init MPI
    MPI_Init(&argc, &argv);
    int rank;
    int world;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);
    MPI_Request request = MPI_REQUEST_NULL;

//  parse args
    if (argc != 8) {
        if (!rank) std::cout << "Wrong number of arguments\n";
        MPI_Finalize();
        return 0;
    }

    double Lx = atof(argv[1]);
    double Ly = atof(argv[2]);
    double Lz = atof(argv[3]);
    int grid_size = atoi(argv[4]);
    double T = atof(argv[5]);
    int time_grid_size = atoi(argv[6]);
    bool save = atoi(argv[7]);

    double dx = Lx / grid_size;
    double dy = Ly / grid_size;
    double dz = Lz / grid_size;
    double dt = T / time_grid_size;

    if (!rank) std::cout << "Params: grid_size - " << grid_size << ", L - " << Lx << std::endl;


//  init cart
    int nproc[3] = {0, 0, 0};
    int periods[3] = {true, true, true};
    int coord[3];
    MPI_Dims_create(world, 3, nproc);
    MPI_Comm cart;
    MPI_Cart_create(MPI_COMM_WORLD, 3, nproc, periods, false, &cart);
    MPI_Cart_coords(cart, rank, 3, coord);

// init grid

    int left_borders[3];
    int right_borders[3];
    int size[3];
    int size_pad[3];

    for (int i = 0; i < 3; ++i) {
        int block_size = (grid_size + 1) / nproc[i];
        left_borders[i] = block_size * coord[i];
        right_borders[i] = block_size * (coord[i] + 1);
        if (coord[i] == nproc[i] - 1) {
            right_borders[i] += (grid_size + 1) % nproc[i];
        }
        size[i] = right_borders[i] - left_borders[i];
        size_pad[i] = size[i] + 2;
    }

    Grid grid(
        size_pad[0],
        size_pad[1],
        size_pad[2],
        dx, dy, dz
    );
    Grid grid_prev(
        size_pad[0],
        size_pad[1],
        size_pad[2],
        dx, dy, dz
    );

    Function function(grid_size, dx, dy, dz, Lx, Ly, Lz);

    function.fill_block(grid, left_borders, size_pad, 0);
    function.fill_block(grid_prev, left_borders, size_pad, dt);

    double max_err = 0;
    int i, j, k;
    double new_value, value_1, value_2, cur_err;

    double start_time = MPI_Wtime();

    #pragma omp parallel for
    for (i = 1; i < size_pad[0] - 1; ++i) {
        for (j = 1; j < size_pad[1] - 1; ++j) {
            for (k = 1; k < size_pad[2] - 1; ++k) {
                int grid_i = left_borders[0] + i - 1;
                int grid_j = left_borders[1] + j - 1;
                int grid_k = left_borders[2] + k - 1;

                new_value = grid.get(i, j, k) + (dt * dt) / (8 * M_PI * M_PI) * grid.get_diff(i, j, k);
                grid_prev.set(i, j, k, new_value);
                cur_err = std::fabs(function.u_from_idx(grid_i, grid_j, grid_k, dt) - new_value);
                max_err = std::max(cur_err, max_err);
            }
        }
    }

    // double global_max_err;
    // MPI_Reduce(&max_err, &global_max_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // if (!rank) std::cout << "Error at time " << dt << " : " << max_err << std::endl;

    int start[3] = {0, 0, 0};
    int x_size[3] = {1, size[1], size[2]};
    int y_size[3] = {size[0], 1, size[2]};
    int z_size[3] = {size[0], size[1], 1};
    int rank_prev, rank_next, is_left, is_right;
    MPI_Datatype x_type, y_type, z_type;
    MPI_Type_create_subarray(3, size_pad, x_size, start, MPI_ORDER_C, MPI_DOUBLE, &x_type);
    MPI_Type_create_subarray(3, size_pad, y_size, start, MPI_ORDER_C, MPI_DOUBLE, &y_type);
    MPI_Type_create_subarray(3, size_pad, z_size, start, MPI_ORDER_C, MPI_DOUBLE, &z_type);
    MPI_Type_commit(&x_type);
    MPI_Type_commit(&y_type);
    MPI_Type_commit(&z_type);


    for (double t = 2; t <= time_grid_size; ++t) {
        is_left = coord[0] == 0;
        is_right = coord[0] == nproc[0] - 1;
        MPI_Cart_shift(cart, 0, 1, &rank_prev, &rank_next);
        MPI_Isend(grid_prev.get_pointer(1 + is_left, 1, 1), 1, x_type, rank_prev, 1, MPI_COMM_WORLD, &request);
        MPI_Request_free(&request);
        MPI_Isend(grid_prev.get_pointer(size[0] - is_right, 1, 1), 1, x_type, rank_next, 2, MPI_COMM_WORLD, &request);
        MPI_Request_free(&request);

        is_left = coord[1] == 0;
        is_right = coord[1] == nproc[1] - 1;
        MPI_Cart_shift(cart, 1, 1, &rank_prev, &rank_next);
        if (!is_left) {
            MPI_Isend(grid_prev.get_pointer(1, 1, 1), 1, y_type, rank_prev, 3, MPI_COMM_WORLD, &request);
            MPI_Request_free(&request);
        } 
        if (!is_right) {
            MPI_Isend(grid_prev.get_pointer(1, size[1], 1), 1, y_type, rank_next, 4, MPI_COMM_WORLD, &request);
            MPI_Request_free(&request);
        }

        is_left = coord[2] == 0;
        is_right = coord[2] == nproc[2] - 1;
        MPI_Cart_shift(cart, 2, 1, &rank_prev, &rank_next);
        if (!is_left) {
            MPI_Isend(grid_prev.get_pointer(1, 1, 1), 1, z_type, rank_prev, 5, MPI_COMM_WORLD, &request);
            MPI_Request_free(&request);
        }
        if (!is_right) {
            MPI_Isend(grid_prev.get_pointer(1, 1, size[2]), 1, z_type, rank_next, 6, MPI_COMM_WORLD, &request);
            MPI_Request_free(&request);
        }

        MPI_Cart_shift(cart, 0, 1, &rank_prev, &rank_next);
        MPI_Recv(grid_prev.get_pointer(0, 1, 1), 1, x_type, rank_prev, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(grid_prev.get_pointer(size_pad[0] - 1, 1, 1), 1, x_type, rank_next, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        is_left = coord[1] == 0;
        is_right = coord[1] == nproc[1] - 1;
        MPI_Cart_shift(cart, 1, 1, &rank_prev, &rank_next);
        if (!is_left) {
            MPI_Recv(grid_prev.get_pointer(1, 0, 1), 1, y_type, rank_prev, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (!is_right) {
            MPI_Recv(grid_prev.get_pointer(1, size_pad[1] - 1, 1), 1, y_type, rank_next, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        is_left = coord[2] == 0;
        is_right = coord[2] == nproc[2] - 1;
        MPI_Cart_shift(cart, 2, 1, &rank_prev, &rank_next);
        if (!is_left) {
            MPI_Recv(grid_prev.get_pointer(1, 1, 0), 1, z_type, rank_prev, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (!is_right) {
            MPI_Recv(grid_prev.get_pointer(1, 1, size_pad[2] - 1), 1, z_type, rank_next, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        max_err = 0;
        #pragma omp parallel for
        for (int i = 1; i < size_pad[0] - 1; ++i) {
            for (int j = 1; j < size_pad[1] - 1; ++j) {
                for (int k = 1; k < size_pad[2] - 1; ++k) {
                    int grid_i = left_borders[0] + i - 1;
                    int grid_j = left_borders[1] + j - 1;
                    int grid_k = left_borders[2] + k - 1;

                    if (
                        grid_j == 0 || grid_k == 0 ||
                        grid_j == grid_size || grid_k == grid_size
                    ) {
                        new_value = 0;
                    } else {
                        value_1 = (dt * dt) / (4 * M_PI * M_PI) * grid_prev.get_diff(i, j, k);
                        value_2 = 2 * grid_prev.get(i, j, k) - grid.get(i, j, k);
                        new_value = value_1 + value_2;
                    }

                    grid.set(i, j, k, new_value);
                    cur_err = std::fabs(function.u_from_idx(grid_i, grid_j, grid_k, t * dt) - new_value);
                    max_err = std::max(cur_err, max_err);
                }
            }
        }
        grid_prev.swap(grid);

        // MPI_Reduce(&max_err, &global_max_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        // if (!rank) std::cout << "Error at time " << t * dt << " : " << global_max_err << std::endl;
    }

    double time = MPI_Wtime() - start_time;
    double max_time;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (!rank) std::cout << "Parallel job time:" << max_time << std::endl;
    MPI_Finalize();
    return 0;
}
