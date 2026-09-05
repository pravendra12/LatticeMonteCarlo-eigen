#include "ClusterGeometry.h"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <stdexcept>


ClusterGeometry::ClusterGeometry(
    const vector<Vector3d> &unwrappedCartesianPositions)
    : unwrappedCartesianPositions_(unwrappedCartesianPositions)
{
  ValidateCartesianPositions();

  clusterSize_ = unwrappedCartesianPositions_.size();

  CalculateGeometricCenter();

  CalculateGyrationTensor();

  CalculatePrincipalMoments();

  CalculateRadiusOfGyration();

  CalculateShapeDescriptors();
}


void ClusterGeometry::ValidateCartesianPositions() const
{
  if (unwrappedCartesianPositions_.empty())
  {
    throw invalid_argument(
        "ClusterGeometry: no Cartesian positions were supplied.");
  }


  for (const auto &cartesianPosition :
       unwrappedCartesianPositions_)
  {
    if (!cartesianPosition.allFinite())
    {
      throw invalid_argument(
          "ClusterGeometry: a Cartesian position contains a non-finite value.");
    }
  }
}


void ClusterGeometry::CalculateGeometricCenter()
{
  geometricCenter_.setZero();


  for (const auto &cartesianPosition :
       unwrappedCartesianPositions_)
  {
    geometricCenter_ += cartesianPosition;
  }


  geometricCenter_ /=
      static_cast<double>(clusterSize_);
}


void ClusterGeometry::CalculateGyrationTensor()
{
  gyrationTensor_.setZero();


  for (const auto &cartesianPosition :
       unwrappedCartesianPositions_)
  {
    const Vector3d positionFromGeometricCenter =
        cartesianPosition - geometricCenter_;


    gyrationTensor_ +=
        positionFromGeometricCenter
        * positionFromGeometricCenter.transpose();
  }


  gyrationTensor_ /=
      static_cast<double>(clusterSize_);
}


void ClusterGeometry::CalculatePrincipalMoments()
{
  SelfAdjointEigenSolver<Matrix3d> eigenSolver(
      gyrationTensor_);


  if (eigenSolver.info() != Success)
  {
    throw runtime_error(
        "ClusterGeometry: failed to diagonalize the gyration tensor.");
  }


  /*
   * Eigen::SelfAdjointEigenSolver returns eigenvalues in ascending order.
   *
   * Therefore the returned ordering already matches the gyration-tensor
   * convention:
   *
   *     lambda_x^2 <= lambda_y^2 <= lambda_z^2.
   *
   * No reordering is necessary.
   */
  principalMoments_ =
      eigenSolver.eigenvalues();


  /*
   * Eigenvectors are returned as columns and correspond directly to the
   * eigenvalues above:
   *
   *     column(0) -> lambda_x^2
   *     column(1) -> lambda_y^2
   *     column(2) -> lambda_z^2.
   */
  principalAxes_ =
      eigenSolver.eigenvectors();
}


void ClusterGeometry::CalculateRadiusOfGyration()
{
  /*
   * By definition:
   *
   *     R_g^2 =
   *         lambda_x^2
   *         + lambda_y^2
   *         + lambda_z^2
   *
   * which is also equal to Tr(S).
   */
  radiusOfGyrationSquared_ =
      max(0.0, principalMoments_.sum());


  radiusOfGyration_ =
      sqrt(radiusOfGyrationSquared_);
}


void ClusterGeometry::CalculateShapeDescriptors()
{
  const double lambdaXSquared =
      principalMoments_(0);

  const double lambdaYSquared =
      principalMoments_(1);

  const double lambdaZSquared =
      principalMoments_(2);


  /*
   * Asphericity:
   *
   *     b =
   *         lambda_z^2
   *         - 1/2 (lambda_x^2 + lambda_y^2)
   */
  asphericity_ =
      lambdaZSquared
      - 0.5 * (lambdaXSquared + lambdaYSquared);


  /*
   * Acylindricity:
   *
   *     c =
   *         lambda_y^2 - lambda_x^2
   */
  acylindricity_ =
      lambdaYSquared - lambdaXSquared;


  /*
   * Relative shape anisotropy:
   *
   *     kappa^2 =
   *         [b^2 + (3/4)c^2] / R_g^4.
   *
   * A cluster whose points all occupy the same position has R_g = 0.
   * In that degenerate case there is no directional anisotropy, so the
   * quantity is set to zero to avoid division by zero.
   */
  if (radiusOfGyrationSquared_ == 0.0)
  {
    relativeShapeAnisotropy_ = 0.0;

    return;
  }


  const double radiusOfGyrationFourth =
      radiusOfGyrationSquared_
      * radiusOfGyrationSquared_;


  relativeShapeAnisotropy_ =
      (
          asphericity_ * asphericity_
          + 0.75 * acylindricity_ * acylindricity_
      )
      / radiusOfGyrationFourth;
}


size_t ClusterGeometry::GetClusterSize() const
{
  return clusterSize_;
}


const Vector3d &ClusterGeometry::GetGeometricCenter() const
{
  return geometricCenter_;
}


const Matrix3d &ClusterGeometry::GetGyrationTensor() const
{
  return gyrationTensor_;
}


const Vector3d &ClusterGeometry::GetPrincipalMoments() const
{
  return principalMoments_;
}


double ClusterGeometry::GetLambdaXSquared() const
{
  return principalMoments_(0);
}


double ClusterGeometry::GetLambdaYSquared() const
{
  return principalMoments_(1);
}


double ClusterGeometry::GetLambdaZSquared() const
{
  return principalMoments_(2);
}


double ClusterGeometry::GetLambdaX() const
{
  return sqrt(max(0.0, GetLambdaXSquared()));
}


double ClusterGeometry::GetLambdaY() const
{
  return sqrt(max(0.0, GetLambdaYSquared()));
}


double ClusterGeometry::GetLambdaZ() const
{
  return sqrt(max(0.0, GetLambdaZSquared()));
}


const Matrix3d &ClusterGeometry::GetPrincipalAxes() const
{
  return principalAxes_;
}


double ClusterGeometry::GetRadiusOfGyrationSquared() const
{
  return radiusOfGyrationSquared_;
}


double ClusterGeometry::GetRadiusOfGyration() const
{
  return radiusOfGyration_;
}


double ClusterGeometry::GetAsphericity() const
{
  return asphericity_;
}


double ClusterGeometry::GetAcylindricity() const
{
  return acylindricity_;
}


double ClusterGeometry::GetRelativeShapeAnisotropy() const
{
  return relativeShapeAnisotropy_;
}