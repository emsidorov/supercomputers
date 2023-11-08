#include <iostream>
#include <cstdlib> // Для использования функции atof и atoi
#include <cmath>
#include <algorithm>


class Grid {
public:
    Grid(int grid_size, double dx, double dy, double dz)
        : Nx(grid_size + 1), Ny(grid_size + 1), Nz(grid_size + 1), dx(dx), dy(dy), dz(dz) {
        array = new double[Nx * Ny * Nz]();
        for (int i = 0; i < Nx * Ny * Nz; ++i) {
            array[i] = 0;
        }
    }

    ~Grid() {
        delete[] array;
    }

    double get(int i, int j, int k) const {
        return array[index(i, j, k)];
    }

    void set(int i, int j, int k, double value) {
        array[index(i, j, k)] = value;
    }

    double get_diff(int i, int j, int k) {
        double diff_x = (this->get(i - 1, j, k) - 2 * this->get(i, j, k) + this->get(i + 1, j, k)) / (dx * dx);
        double diff_y = (this->get(i, j - 1, k) - 2 * this->get(i, j, k) + this->get(i, j + 1, k)) / (dy * dy);
        double diff_z = (this->get(i, j, k - 1) - 2 * this->get(i, j, k) + this->get(i, j, k + 1)) / (dz * dz);
        return diff_x + diff_y + diff_z;
    }

    int size_x() const { return Nx; }
    int size_y() const { return Ny; }
    int size_z() const { return Nz; }

    void swap(Grid& other) {
        if (this->Nx != other.Nx || this->Ny != other.Ny || this->Nz != other.Nz) {
            throw std::invalid_argument("Grid sizes do not match and cannot be swapped.");
        }
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

    void set_border(Grid& grid, int i, int j, int k, double t) {
        if (j == 0 && k == 0) {
            grid.set(i, j, k, 0);
        } else {
            grid.set(i, j, k, this->u_from_idx(i, j, k, t));
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
        grid1.size_z() != grid2.size_z()) {
        throw std::invalid_argument("Grid sizes do not match.");
    }

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
    if (argc != 7) {
        std::cerr << "Usage: " << argv[0] << " Lx Ly Lz grid_size T time_size" << std::endl;
        return 1;
    }

    double Lx = atof(argv[1]);
    double Ly = atof(argv[2]);
    double Lz = atof(argv[3]);
    int grid_size = atoi(argv[4]);
    double T = atof(argv[5]);
    int time_grid_size = atoi(argv[6]);

    double dx = Lx / grid_size;
    double dy = Ly / grid_size;
    double dz = Lz / grid_size;
    double dt = T / time_grid_size;

    Grid grid(grid_size, dx, dy, dz);
    Grid grid_prev(grid_size, dx, dy, dz);

    Function function(grid_size, dx, dy, dz, Lx, Ly, Lz);

    function.fill_grid(grid, 0);
    function.fill_grid(grid_prev, dt);

    double max_err = 0;
    for (int i = 1; i < grid_size; ++i) {
        for (int j = 1; j < grid_size; ++j) {
            for (int k = 1; k < grid_size; ++k) {
                double new_value = grid.get(i, j, k) + (dt * dt) / (4 * M_PI * M_PI) * grid.get_diff(i, j, k);
                grid_prev.set(i, j, k, new_value);
                double cur_err = std::fabs(function.u_from_idx(i, j, k, dt) - new_value);
                max_err = std::max(cur_err, max_err);
            }
        }
    }
    std::cout << "Error at time " << dt << " : " << max_err << std::endl;

    for (double t = 2; t <= 20; ++t) {
        max_err = 0;
        for (int i = 1; i < grid_size; ++i) {
            for (int j = 1; j < grid_size; ++j) {
                for (int k = 1; k < grid_size; ++k) {
                    double value_1 = (dt * dt) / (2 * M_PI * M_PI) * grid_prev.get_diff(i, j, k);
                    double value_2 = 2 * grid_prev.get(i, j, k) - grid.get(i, j, k);
                    double new_value = value_1 + value_2;
                    grid.set(i, j, k, new_value);

                    double cur_err = std::fabs(function.u_from_idx(i, j, k, t * dt) - new_value);
                    max_err = std::max(cur_err, max_err);
                }
            }
        }
        grid_prev.swap(grid);
        // function.fill_grid(grid_prev, t * dt);
        std::cout << "Error at time " << t * dt << " : " << max_err << std::endl;
    }
    return 0;
}
