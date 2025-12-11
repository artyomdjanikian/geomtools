#include <Eigen/Dense>
#include <vector>
#include <iostream>

struct Point3D
{
    double x, y, z;
};

Eigen::Vector3d computeGradient(
    const Point3D &p0,
    double f0,
    const std::vector<Point3D> &neighbors,
    const std::vector<double> &values)
{
    size_t N = neighbors.size();
    Eigen::MatrixXd A(N, 3);
    Eigen::VectorXd b(N);

    for (size_t i = 0; i < N; ++i)
    {
        A(i, 0) = neighbors[i].x - p0.x;
        A(i, 1) = neighbors[i].y - p0.y;
        A(i, 2) = neighbors[i].z - p0.z;

        b(i) = values[i] - f0;
    }

    // Solve normal equations: g = (A^T A)^{-1} A^T b
    Eigen::Vector3d gradient = (A.transpose() * A).ldlt().solve(A.transpose() * b);

    return gradient;
}

int main()
{
    Point3D p0{0.0, 0.0, 0.0};
    double f0 = 1.0;

    std::vector<Point3D> neighbors = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {1.0, 1.0, 0.0}};

    std::vector<double> values = {2.0, 3.0, 4.0, 5.0};

    Eigen::Vector3d grad = computeGradient(p0, f0, neighbors, values);

    std::cout << "Gradient = " << grad.transpose() << std::endl;
    return 0;
}
#include <Eigen/Dense>
#include <vector>
#include <cmath>

struct Point3D
{
    double x, y, z;
};

struct Vec3
{
    double x, y, z;
};

//
// Weight kernel: compact quartic kernel.
// r = distance / h
//
double kernelWeight(double r)
{
    if (r >= 1.0)
        return 0.0;
    double s = 1.0 - r * r;
    return s * s;
}

//
// Compute weighted MLS derivatives for a scalar function.
// Produces gradient and full Hessian.
//
void mlsScalarDerivatives(
    const Point3D &p0,
    double f0,
    const std::vector<Point3D> &nbrs,
    const std::vector<double> &fvals,
    double h,
    Eigen::Vector3d &gradient,
    Eigen::Matrix3d &hessian)
{
    const int N = nbrs.size();
    Eigen::MatrixXd A(N, 9);
    Eigen::VectorXd b(N);
    Eigen::VectorXd w(N);

    for (int i = 0; i < N; ++i)
    {
        double dx = nbrs[i].x - p0.x;
        double dy = nbrs[i].y - p0.y;
        double dz = nbrs[i].z - p0.z;

        double r = std::sqrt(dx * dx + dy * dy + dz * dz) / h;
        w(i) = kernelWeight(r);

        A(i, 0) = dx;
        A(i, 1) = dy;
        A(i, 2) = dz;
        A(i, 3) = 0.5 * dx * dx;
        A(i, 4) = 0.5 * dy * dy;
        A(i, 5) = 0.5 * dz * dz;
        A(i, 6) = dx * dy;
        A(i, 7) = dx * dz;
        A(i, 8) = dy * dz;

        b(i) = fvals[i] - f0;
    }

    // Weighted normal-equation solve: (A^T W A)c = A^T W b
    Eigen::MatrixXd AW = A.transpose() * w.asDiagonal();
    Eigen::VectorXd rhs = AW * b;
    Eigen::MatrixXd ATA = AW * A;

    Eigen::VectorXd c = ATA.ldlt().solve(rhs);

    // Gradient
    gradient(0) = c(0);
    gradient(1) = c(1);
    gradient(2) = c(2);

    // Hessian
    hessian << c(3), c(6), c(7),
        c(6), c(4), c(8),
        c(7), c(8), c(5);
}

//
// Compute Laplacian of scalar field f = trace(Hessian)
//
double mlsLaplacian(const Eigen::Matrix3d &H)
{
    return H(0, 0) + H(1, 1) + H(2, 2);
}

//
// Compute divergence and curl for vector field u = (u,v,w)
//
void mlsVectorDerivatives(
    const Point3D &p0,
    const Vec3 &u0,
    const std::vector<Point3D> &nbrs,
    const std::vector<Vec3> &uvals,
    double h,
    double &divergence,
    Vec3 &curl)
{
    const int N = nbrs.size();

    // Separate components
    std::vector<double> ux(N), uy(N), uz(N);
    for (int i = 0; i < N; ++i)
    {
        ux[i] = uvals[i].x;
        uy[i] = uvals[i].y;
        uz[i] = uvals[i].z;
    }

    Eigen::Vector3d gx, gy, gz;
    Eigen::Matrix3d Htmp;

    mlsScalarDerivatives(p0, u0.x, nbrs, ux, h, gx, Htmp);
    mlsScalarDerivatives(p0, u0.y, nbrs, uy, h, gy, Htmp);
    mlsScalarDerivatives(p0, u0.z, nbrs, uz, h, gz, Htmp);

    // Divergence = du/dx + dv/dy + dw/dz
    divergence = gx(0) + gy(1) + gz(2);

    // Curl = ∇×u = (dw/dy - dv/dz, du/dz - dw/dx, dv/dx - du/dy)
    curl.x = gz(1) - gy(2);
    curl.y = gx(2) - gz(0);
    curl.z = gy(0) - gx(1);
}

✔ Example Usage

    Point3D p0{0, 0, 0};
double f0 = 1.0;

std::vector<Point3D> nbrs = ...;
std::vector<double> fvals = ...;

double h = 0.2; // MLS smoothing radius

Eigen::Vector3d grad;
Eigen::Matrix3d H;
mlsScalarDerivatives(p0, f0, nbrs, fvals, h, grad, H);

double lap = mlsLaplacian(H);
std::cout << "Gradient = " << grad.transpose() << "\n";
std::cout << "Laplacian = " << lap << "\n";

Vector field example :

    Vec3 u0{1, 2, 3};
std::vector<Vec3> uvals = ...;

double div;
Vec3 curl;

mlsVectorDerivatives(p0, u0, nbrs, uvals, h, div, curl);

std::cout << "Divergence = " << div << "\n";
std::cout << "Curl = " << curl.x << " " << curl.y << " " << curl.z << "\n";

👍 Benefits of MLS over simple least squares

        Handles irregular point distributions robustly

            Smooth
            derivatives(C² fit)

                Local radius h controls smoothing

    Noise
    - resistant

    Works on point clouds and particle systems(SPH - like)

        If you want,
    I can also provide :

✅ Automatic selection of optimal smoothing radius
✅ Gaussian or Wendland kernels
✅ Fast neighbor precomputation(KD - tree)
✅ GPU / CUDA MLS version
✅ Vector Laplacian(important in Navier–Stokes)

                        Just tell me !Get smarter responses,
    upload files and images, and more.

                                 ChatGPT can make mistakes.Check important info.