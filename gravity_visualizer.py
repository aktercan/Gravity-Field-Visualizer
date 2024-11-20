import numpy as np
import matplotlib.pyplot as plt


def read_potential_coefficients(filename):
    """Reads spherical harmonic coefficients and header data from a GFC file."""
    with open(filename, 'r') as file:
        has_errors = True

        max_degree = None
        GM = None
        R = None

        while True:
            line = file.readline()
            if not line:
                break
            line = line.strip()
            if not line:
                continue

            keyword = line.split()[0]

            if keyword == 'max_degree':
                max_degree = int(line.split()[1])
            elif keyword == 'radius':
                R = float(line.split()[1])
            elif keyword == 'earth_gravity_constant':
                GM = float(line.split()[1])
            elif keyword == 'errors':
                if line.split()[1] == 'no':
                    has_errors = False
            elif keyword == 'end_of_head':
                break

        if max_degree is None or GM is None or R is None:
            raise ValueError("Header does not contain necessary information")

        cnm = np.zeros((max_degree + 1, max_degree + 1))
        snm = np.zeros((max_degree + 1, max_degree + 1))
        cnm_error = np.zeros((max_degree + 1, max_degree + 1))
        snm_error = np.zeros((max_degree + 1, max_degree + 1))

        while True:
            line = file.readline()
            if not line:
                break
            line = line.strip()
            if not line or not line.startswith('gfc'):
                continue

            cells = line.split()
            n = int(cells[1])
            m = int(cells[2])
            cnm[n, m] = float(cells[3])
            snm[n, m] = float(cells[4])

            if has_errors:
                cnm_error[n, m] = float(cells[5])
                snm_error[n, m] = float(cells[6])

    return cnm, snm, GM, R


def filter_coefficients_gaussian(radius, max_degree):
    """Applies a Gaussian filter to coefficients up to max_degree."""
    R = 6378.137
    b = np.log(2) / (1 - np.cos(radius / R))

    wn = np.zeros(max_degree + 1)
    wn[0] = 1
    wn[1] = (1 + np.exp(-2 * b)) / (1 - np.exp(-2 * b)) - 1 / b

    for n in range(2, max_degree + 1):
        wn[n] = -(2 * n - 1) / b * wn[n - 1] + wn[n - 2]

    return wn


def legendre_functions(theta, max_degree):
    """Computes normalized associated Legendre polynomials."""
    global legendre_factor1, legendre_factor2

    if 'legendre_factor1' not in globals() or legendre_factor1.shape[0] < max_degree + 1:
        legendre_factor1 = np.zeros((max_degree + 1, max_degree + 1))
        legendre_factor2 = np.zeros((max_degree + 1, max_degree + 1))
        legendre_factor1[1, 1] = np.sqrt(3)
        for n in range(2, max_degree + 1):
            legendre_factor1[n, n] = np.sqrt((2 * n + 1) / (2 * n))
        for m in range(max_degree):
            for n in range(m + 1, max_degree + 1):
                f = (2 * n + 1) / ((n + m) * (n - m))
                legendre_factor1[n, m] = np.sqrt(f * (2 * n - 1))
                legendre_factor2[n, m] = np.sqrt(f * (n - m - 1) * (n + m - 1) / (2 * n - 3))

    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    Pnm = np.zeros((max_degree + 1, max_degree + 1))
    Pnm[0, 0] = 1

    for n in range(1, max_degree + 1):
        Pnm[n, n] = legendre_factor1[n, n] * sin_theta * Pnm[n - 1, n - 1]

    for m in range(max_degree):
        Pnm[m + 1, m] = legendre_factor1[m + 1, m] * cos_theta * Pnm[m, m]

    for m in range(max_degree):
        for n in range(m + 2, max_degree + 1):
            Pnm[n, m] = (legendre_factor1[n, m] * cos_theta * Pnm[n - 1, m]
                         - legendre_factor2[n, m] * Pnm[n - 2, m])

    return Pnm


def gravity_disturbance(Lambda, Theta, cnm, snm, GM, R, max_degree=30):
    """Calculates gravity disturbance on a grid of points."""
    L = len(Lambda)
    T = len(Theta)
    gd = np.zeros((T, L))

    for i in range(L):
        for j in range(T):
            sum_nm = 0
            Pnm = legendre_functions(Theta[j], max_degree)

            for n in range(max_degree + 1):
                for m in range(n + 1):
                    Pnm_value = Pnm[n, m]
                    cosm_L = np.cos(m * Lambda[i])
                    sinm_L = np.sin(m * Lambda[i])
                    cnm_term = cnm[n, m] * cosm_L
                    snm_term = snm[n, m] * sinm_L
                    sum_nm += (n + 1) * (cnm_term + snm_term) * Pnm_value

            gd[j, i] = GM / (R ** 2) * sum_nm

    return gd


def show_grid(lambda_vals, theta_vals, grid, titletext):
    """Plots the computed gravity grid."""
    plt.figure(figsize=(12, 8))
    plt.contourf(np.degrees(lambda_vals), 90 - np.degrees(theta_vals), grid, 70, cmap='jet')
    plt.colorbar()

    try:
        coast = np.loadtxt('coast.dat')
        plt.plot(coast[:, 0], coast[:, 1], 'k', linewidth=1)
    except OSError:
        print("Coastline data not found. Skipping coastline plot.")

    plt.title(titletext, fontsize=16, fontweight='bold')
    plt.xlabel('Longitude', fontsize=12)
    plt.ylabel('Latitude', fontsize=12)
    plt.axis('tight')
    plt.contour(np.degrees(lambda_vals), 90 - np.degrees(theta_vals), grid, 70, colors='k', linestyles='solid',
                linewidths=0.5)
    plt.show()


def calculate_gravity_anomalies(file1, file2, max_degree, filter_radius):
    """High-level function to compute gravity anomalies between two datasets."""
    cnm1, snm1, GM, R = read_potential_coefficients(file1)
    cnm2, snm2, GM, R = read_potential_coefficients(file2)

    cnm_diff = cnm2[:max_degree + 1, :max_degree + 1] - cnm1[:max_degree + 1, :max_degree + 1]
    snm_diff = snm2[:max_degree + 1, :max_degree + 1] - snm1[:max_degree + 1, :max_degree + 1]

    wn = filter_coefficients_gaussian(filter_radius, max_degree)
    for n in range(max_degree + 1):
        cnm_diff[n, :] *= wn[n]
        snm_diff[n, :] *= wn[n]

    Theta = np.arange(0, 181) * np.pi / 180
    Lambda = np.arange(-180, 181) * np.pi / 180
    gd = gravity_disturbance(Lambda, Theta, cnm_diff, snm_diff, GM, R, max_degree)
    return Lambda, Theta, gd


def main():
    first_file = 'EGSIEM_COMB_90_NEQ_2007_03.gfc'
    second_file = 'EGSIEM_COMB_90_NEQ_2007_09.gfc'
    max_degree = 30
    filter_radius = 500

    Lambda, Theta, gd = calculate_gravity_anomalies(first_file, second_file, max_degree, filter_radius)
    show_grid(Lambda, Theta, gd, 'Gravity Anomalies')


if __name__ == '__main__':
    main()