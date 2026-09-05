#ifndef LMC_ANSYS_INCLUDE_CLUSTER_GEOMETRY_H
#define LMC_ANSYS_INCLUDE_CLUSTER_GEOMETRY_H

#include <Eigen/Dense>

#include <vector>

using namespace std;
using namespace Eigen;


/**
 * @brief Calculates the position-based gyration tensor and associated
 * geometric properties of a cluster.
 *
 * The ClusterGeometry class performs a purely geometric analysis of a
 * collection of Cartesian positions. No atomic masses, chemical species,
 * lattice information, energetic quantities, or cluster-identification
 * information are used in the calculation.
 *
 * For a cluster containing N Cartesian positions r_i, the geometric center
 * is defined as
 *
 *     r_c = (1 / N) * sum_i r_i .
 *
 * The gyration tensor is then calculated as
 *
 *     S = (1 / N) * sum_i
 *         (r_i - r_c) (r_i - r_c)^T .
 *
 * The gyration tensor is a real symmetric 3 x 3 matrix. In the original
 * Cartesian coordinate system it has the form
 *
 *          | S_xx  S_xy  S_xz |
 *     S =  | S_xy  S_yy  S_yz |
 *          | S_xz  S_yz  S_zz | .
 *
 * After diagonalization, the tensor is represented in its principal-axis
 * coordinate system as
 *
 *           | lambda_x^2       0             0       |
 *     S' =  |     0        lambda_y^2        0       |
 *           |     0            0         lambda_z^2  | ,
 *
 * where
 *
 *     lambda_x^2 <= lambda_y^2 <= lambda_z^2 .
 * 
 * @see https://en.wikipedia.org/wiki/Gyration_tensor
 *
 * IMPORTANT:
 *
 * lambda_x^2, lambda_y^2, and lambda_z^2 are the actual eigenvalues of the
 * gyration tensor.
 *
 * The squared notation is used because every component of the gyration
 * tensor is constructed from a product of Cartesian displacements and
 * therefore has units of length squared. Consequently, the eigenvalues also
 * have units of length squared.
 *
 * The corresponding characteristic lengths are
 *
 *     lambda_x = sqrt(lambda_x^2),
 *
 *     lambda_y = sqrt(lambda_y^2),
 *
 *     lambda_z = sqrt(lambda_z^2),
 *
 * with
 *
 *     lambda_x <= lambda_y <= lambda_z .
 *
 * The squared radius of gyration is
 *
 *     R_g^2 =
 *         lambda_x^2
 *         + lambda_y^2
 *         + lambda_z^2 ,
 *
 * which is equivalent to
 *
 *     R_g^2 = Tr(S) .
 *
 * Therefore,
 *
 *     R_g =
 *         sqrt(
 *             lambda_x^2
 *             + lambda_y^2
 *             + lambda_z^2
 *         ) .
 *
 * Additional shape descriptors are calculated from the principal moments.
 *
 * The asphericity is
 *
 *     b =
 *         lambda_z^2
 *         - 0.5 * (lambda_x^2 + lambda_y^2) .
 *
 * The acylindricity is
 *
 *     c =
 *         lambda_y^2
 *         - lambda_x^2 .
 *
 * The relative shape anisotropy is
 *
 *     kappa^2 =
 *         [ b^2 + (3 / 4)c^2 ] / R_g^4 .
 *
 * The relative shape anisotropy is dimensionless and lies between zero and
 * one:
 *
 *     kappa^2 = 0
 *
 * corresponds to an isotropic spatial distribution, while
 *
 *     kappa^2 = 1
 *
 * corresponds to an ideal one-dimensional distribution of points.
 *
 * @note All calculations are position-based and unweighted. This class does
 *       not calculate a mass-weighted gyration tensor.
 *
 * @note The Cartesian positions supplied to this class must already be
 *       unwrapped with respect to periodic boundary conditions. A cluster
 *       crossing a periodic simulation-cell boundary must therefore be
 *       represented as one spatially continuous object before constructing
 *       ClusterGeometry.
 *
 * @note The x, y, and z labels associated with lambda_x, lambda_y, and
 *       lambda_z refer to the principal-axis coordinate system of the
 *       gyration tensor. They do not necessarily correspond to the original
 *       Cartesian x, y, and z axes of the simulation cell.
 *
 * @note The eigenvectors defining the principal axes are determined only up
 *       to sign. Therefore, a principal-axis vector v and the vector -v
 *       represent the same physical axis.
 */
class ClusterGeometry
{
public:

  /**
   * @brief Constructs the geometric description of a cluster.
   *
   * All geometric quantities are calculated during construction from the
   * supplied unwrapped Cartesian positions.
   *
   * The constructor calculates:
   *
   * - cluster size,
   * - geometric center,
   * - gyration tensor,
   * - principal moments,
   * - principal axes,
   * - squared radius of gyration,
   * - radius of gyration,
   * - asphericity,
   * - acylindricity,
   * - relative shape anisotropy.
   *
   * @param unwrappedCartesianPositions
   *        Cartesian positions of all spatial points belonging to the
   *        cluster. The coordinates must already be unwrapped across
   *        periodic boundaries.
   *
   * @throws invalid_argument
   *         If no Cartesian positions are supplied or if any position
   *         contains a non-finite coordinate.
   */
  explicit ClusterGeometry(
      const vector<Vector3d> &unwrappedCartesianPositions);


  /**
   * @brief Returns the number of spatial points defining the cluster.
   *
   * @return Number of Cartesian positions used in the calculation.
   */
  size_t GetClusterSize() const;


  /**
   * @brief Returns the geometric center of the cluster.
   *
   * The geometric center is the arithmetic mean of all supplied positions:
   *
   *     r_c = (1 / N) * sum_i r_i .
   *
   * Every spatial point contributes equally.
   *
   * @return Geometric center in Cartesian coordinates.
   */
  const Vector3d &GetGeometricCenter() const;


  /**
   * @brief Returns the position-based gyration tensor.
   *
   * The tensor is calculated as
   *
   *     S = (1 / N) * sum_i
   *         (r_i - r_c) (r_i - r_c)^T .
   *
   * The returned tensor is represented in the original Cartesian coordinate
   * system.
   *
   * @return Symmetric 3 x 3 gyration tensor.
   */
  const Matrix3d &GetGyrationTensor() const;


  /**
   * @brief Returns the principal moments of the gyration tensor.
   *
   * The returned vector contains
   *
   *     [lambda_x^2, lambda_y^2, lambda_z^2]
   *
   * with
   *
   *     lambda_x^2 <= lambda_y^2 <= lambda_z^2 .
   *
   * These quantities are the actual eigenvalues of the gyration tensor and
   * therefore have units of length squared.
   *
   * @return Eigenvalues of the gyration tensor ordered from smallest to
   *         largest.
   */
  const Vector3d &GetPrincipalMoments() const;


  /**
   * @brief Returns lambda_x^2.
   *
   * lambda_x^2 is the smallest eigenvalue of the gyration tensor.
   *
   * @return Smallest principal moment, with units of length squared.
   */
  double GetLambdaXSquared() const;


  /**
   * @brief Returns lambda_y^2.
   *
   * lambda_y^2 is the intermediate eigenvalue of the gyration tensor.
   *
   * @return Intermediate principal moment, with units of length squared.
   */
  double GetLambdaYSquared() const;


  /**
   * @brief Returns lambda_z^2.
   *
   * lambda_z^2 is the largest eigenvalue of the gyration tensor and
   * corresponds to the principal direction having the largest spatial
   * variance.
   *
   * @return Largest principal moment, with units of length squared.
   */
  double GetLambdaZSquared() const;


  /**
   * @brief Returns lambda_x.
   *
   * lambda_x is not itself an eigenvalue of the gyration tensor. It is the
   * characteristic length obtained from the smallest eigenvalue:
   *
   *     lambda_x = sqrt(lambda_x^2).
   *
   * @return Smallest principal characteristic length.
   */
  double GetLambdaX() const;


  /**
   * @brief Returns lambda_y.
   *
   * lambda_y is obtained from
   *
   *     lambda_y = sqrt(lambda_y^2).
   *
   * @return Intermediate principal characteristic length.
   */
  double GetLambdaY() const;


  /**
   * @brief Returns lambda_z.
   *
   * lambda_z is obtained from
   *
   *     lambda_z = sqrt(lambda_z^2).
   *
   * It corresponds to the principal direction having the largest spatial
   * variance.
   *
   * @return Largest principal characteristic length.
   */
  double GetLambdaZ() const;


  /**
   * @brief Returns the principal axes of the gyration tensor.
   *
   * The normalized eigenvectors are stored as columns:
   *
   *     column(0) -> principal x-axis -> lambda_x^2
   *
   *     column(1) -> principal y-axis -> lambda_y^2
   *
   *     column(2) -> principal z-axis -> lambda_z^2 .
   *
   * Consequently, column(2) represents the direction of maximum spatial
   * variance of the cluster.
   *
   * @return Matrix containing the three normalized principal-axis vectors.
   */
  const Matrix3d &GetPrincipalAxes() const;


  /**
   * @brief Returns the squared radius of gyration.
   *
   * The squared radius of gyration is
   *
   *     R_g^2 =
   *         lambda_x^2
   *         + lambda_y^2
   *         + lambda_z^2 .
   *
   * This quantity is also equal to the trace of the gyration tensor.
   *
   * @return Squared radius of gyration, with units of length squared.
   */
  double GetRadiusOfGyrationSquared() const;


  /**
   * @brief Returns the radius of gyration.
   *
   * The radius of gyration is
   *
   *     R_g =
   *         sqrt(
   *             lambda_x^2
   *             + lambda_y^2
   *             + lambda_z^2
   *         ) .
   *
   * It provides a rotationally invariant measure of the overall spatial
   * extent of the cluster.
   *
   * @return Radius of gyration.
   */
  double GetRadiusOfGyration() const;


  /**
   * @brief Returns the asphericity of the cluster.
   *
   * The asphericity is defined as
   *
   *     b =
   *         lambda_z^2
   *         - 0.5 * (lambda_x^2 + lambda_y^2) .
   *
   * For an isotropic distribution where all three principal moments are
   * equal, the asphericity is zero.
   *
   * @return Asphericity b.
   */
  double GetAsphericity() const;


  /**
   * @brief Returns the acylindricity of the cluster.
   *
   * The acylindricity is defined as
   *
   *     c =
   *         lambda_y^2
   *         - lambda_x^2 .
   *
   * Since the principal moments are stored in ascending order, the
   * acylindricity is non-negative.
   *
   * @return Acylindricity c.
   */
  double GetAcylindricity() const;


  /**
   * @brief Returns the relative shape anisotropy.
   *
   * The relative shape anisotropy is calculated as
   *
   *     kappa^2 =
   *         [ b^2 + (3 / 4)c^2 ] / R_g^4 .
   *
   * It is dimensionless and provides a normalized measure of cluster-shape
   * anisotropy.
   *
   *     kappa^2 = 0
   *
   * corresponds to an isotropic distribution, whereas
   *
   *     kappa^2 = 1
   *
   * corresponds to an ideal linear distribution.
   *
   * @return Relative shape anisotropy kappa^2.
   */
  double GetRelativeShapeAnisotropy() const;


private:

  /**
   * @brief Validates the supplied Cartesian coordinates.
   *
   * @throws invalid_argument
   *         If the position vector is empty or contains non-finite values.
   */
  void ValidateCartesianPositions() const;


  /**
   * @brief Calculates the geometric center of the cluster.
   *
   * The result is stored in geometricCenter_.
   */
  void CalculateGeometricCenter();


  /**
   * @brief Calculates the position-based gyration tensor.
   *
   * The tensor is calculated relative to geometricCenter_ and stored in
   * gyrationTensor_.
   */
  void CalculateGyrationTensor();


  /**
   * @brief Diagonalizes the gyration tensor.
   *
   * The gyration tensor is real and symmetric, so a self-adjoint
   * eigensolver is used.
   *
   * The resulting eigenvalues are stored as
   *
   *     [lambda_x^2, lambda_y^2, lambda_z^2]
   *
   * with
   *
   *     lambda_x^2 <= lambda_y^2 <= lambda_z^2 .
   *
   * The corresponding eigenvectors are stored as columns of principalAxes_
   * using the same ordering.
   *
   * @throws runtime_error
   *         If diagonalization of the gyration tensor fails.
   */
  void CalculatePrincipalMoments();


  /**
   * @brief Calculates R_g^2 and R_g from the principal moments.
   */
  void CalculateRadiusOfGyration();


  /**
   * @brief Calculates the shape descriptors from the principal moments.
   *
   * The calculated quantities are:
   *
   * - asphericity,
   * - acylindricity,
   * - relative shape anisotropy.
   */
  void CalculateShapeDescriptors();


  /**
   * @brief Unwrapped Cartesian positions defining the cluster.
   */
  vector<Vector3d> unwrappedCartesianPositions_;


  /**
   * @brief Number of spatial points defining the cluster.
   */
  size_t clusterSize_ = 0;


  /**
   * @brief Geometric center of the cluster in Cartesian coordinates.
   */
  Vector3d geometricCenter_ = Vector3d::Zero();


  /**
   * @brief Position-based gyration tensor.
   *
   * Tensor components have units of length squared.
   */
  Matrix3d gyrationTensor_ = Matrix3d::Zero();


  /**
   * @brief Eigenvalues of the gyration tensor.
   *
   * Stored as
   *
   *     [lambda_x^2, lambda_y^2, lambda_z^2]
   *
   * with
   *
   *     lambda_x^2 <= lambda_y^2 <= lambda_z^2 .
   *
   * These quantities are the actual eigenvalues of the gyration tensor.
   */
  Vector3d principalMoments_ = Vector3d::Zero();


  /**
   * @brief Principal-axis vectors corresponding to the principal moments.
   *
   * Columns correspond respectively to
   *
   *     lambda_x^2,
   *     lambda_y^2,
   *     lambda_z^2.
   */
  Matrix3d principalAxes_ = Matrix3d::Identity();


  /**
   * @brief Squared radius of gyration R_g^2.
   */
  double radiusOfGyrationSquared_ = 0.0;


  /**
   * @brief Radius of gyration R_g.
   */
  double radiusOfGyration_ = 0.0;


  /**
   * @brief Asphericity b.
   */
  double asphericity_ = 0.0;


  /**
   * @brief Acylindricity c.
   */
  double acylindricity_ = 0.0;


  /**
   * @brief Relative shape anisotropy kappa^2.
   */
  double relativeShapeAnisotropy_ = 0.0;
};





#endif // LMC_ANSYS_INCLUDE_CLUSTER_GEOMETRY_H