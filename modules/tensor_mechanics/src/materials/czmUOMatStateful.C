//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "czmUOMatStateful.h"
// #include "MooseMesh.h"
#include "RotationMatrix.h"

registerMooseObject("TensorMechanicsApp", czmUOMatStateful);

template <>
InputParameters
validParams<czmUOMatStateful>()
{
  InputParameters params = validParams<Material>();

  params.addRequiredParam<UserObjectName>(
      "displacement_jump_UO", "the name of the computing material property across an interface");

  params.addRequiredParam<Real>("DeltaU0",
                                "a scalar representing the effective displacement"
                                "at maximum stress");
  params.addRequiredParam<Real>("MaxAllowableTraction",
                                "a scalar representing the maximum allowable stress");
  params.addRequiredParam<Real>("Beta",
                                "a coeffcient for waiting the effect o shear traction in mix mode");

  params.addClassDescription("this material class is used when defining a "
                             "cohesive zone model to store stafeul properties");
  return params;
}

czmUOMatStateful::czmUOMatStateful(const InputParameters & parameters)
  : Material(parameters),
    _deltaU0(getParam<Real>("DeltaU0")),
    _maxAllowableTraction(getParam<Real>("MaxAllowableTraction")),
    _beta(getParam<Real>("Beta")),
    _displacement_jump_UO(getUserObject<DispJumpUO_QP>("displacement_jump_UO")),
    _displacement_jump(declareProperty<std::vector<Real>>("displacement_jump")),
    _displacement_jump_local(declareProperty<std::vector<Real>>("displacement_jump_local")),
    _max_effective_jump(declareProperty<Real>("max_effective_jump")),
    _max_effective_jump_old(getMaterialPropertyOld<Real>("max_effective_jump")),
    _effective_jump(declareProperty<Real>("effective_jump")),
    _effective_jump_old(getMaterialPropertyOld<Real>("effective_jump")),
    _traction(declareProperty<std::vector<Real>>("traction")),
    _traction_local(declareProperty<std::vector<Real>>("traction_local")),
    _effective_traction(declareProperty<Real>("effective_traction")),
    _traction_spatial_derivatives(
        declareProperty<std::vector<std::vector<Real>>>("traction_spatial_derivatives")),
    _traction_spatial_derivatives_local(
        declareProperty<std::vector<std::vector<Real>>>("traction_spatial_derivatives_local")),
    _residual(declareProperty<std::vector<Real>>("czmResidual")),
    _jacobian(declareProperty<std::vector<std::vector<Real>>>("czmJacobian"))

{
}

void
czmUOMatStateful::initQpStatefulProperties()
{
  _max_effective_jump[_qp] = 0.0;
  _effective_jump[_qp] = 0.0;
}

Real
czmUOMatStateful::computeEffectiveJump()
{
  Real effective_jump = 0;
  for (unsigned int i = 0; i < 3; i++)
  {
    Real temp = _displacement_jump_local[_qp][i] * _displacement_jump_local[_qp][i];
    if (i > 0)
      temp *= _beta * _beta;
    effective_jump += temp;
  }
  return std::sqrt(effective_jump);
}

Real
czmUOMatStateful::computeEffectiveTractionNonLinear(Real effective_jump)
{
  // compute non linear traction divide by effective_jump
  Real effective_traction_nl = 0;
  if (effective_jump != 0)
  {
    Real d_norm = effective_jump / _deltaU0;
    effective_traction_nl = std::exp(1) * _maxAllowableTraction / _deltaU0 * std::exp(-1 * d_norm);
  }
  return effective_traction_nl;
}

Real
czmUOMatStateful::computeEffectiveTractionLinear()
{
  // Real effective_traction_l = 0;
  return computeEffectiveTractionNonLinear(_max_effective_jump_old[_qp]);
}

bool
czmUOMatStateful::checkLoadUnload()
{
  bool is_unloading = false;
  if (_effective_jump[_qp] < _effective_jump_old[_qp]) /*unloading*/
    is_unloading = true;

  return is_unloading;
}

Real
czmUOMatStateful::computeEffectiveTraction()
{
  if (!checkLoadUnload()) /* loading */
  {
    if (_max_effective_jump_old[_qp] < _deltaU0) /* not damaged*/
      return computeEffectiveTractionNonLinear(_effective_jump[_qp]);
    else /*damaged*/
    {
      if (_effective_jump[_qp] <= _max_effective_jump_old[_qp]) /* linear region*/
        return computeEffectiveTractionLinear();
      else
        return computeEffectiveTractionNonLinear(_effective_jump[_qp]); /* NON linear region*/
    }
  }
  else /* unload */
  {
    return computeEffectiveTractionLinear();
  }
}

void
czmUOMatStateful::computeQpProperties()
{

  _displacement_jump[_qp].resize(3, 0);
  _displacement_jump_local[_qp].resize(3, 0);
  _traction[_qp].resize(3, 0);
  _traction_local[_qp].resize(3, 0);
  _traction_spatial_derivatives[_qp].resize(3, std::vector<Real>(3, 0));
  _traction_spatial_derivatives_local[_qp].resize(3, std::vector<Real>(3, 0));

  _residual[_qp].resize(3, 0);
  _jacobian[_qp].resize(3, std::vector<Real>(3, 0));

  _displacement_jump[_qp] =
      _displacement_jump_UO.getDisplacementJump(_current_elem->id(), _current_side, _qp);

  RealTensorValue RotationGlobal2Local =
      RotationMatrix::rotVec1ToVec2(_normals[_qp], RealVectorValue(1, 0, 0));
  _displacement_jump_local[_qp] = rotateVector(_displacement_jump[_qp], RotationGlobal2Local);

  _traction_local[_qp] = computeTractionLocal();
  _traction_spatial_derivatives_local[_qp] = computeTractionSpatialDerivativeLocal();
  _traction[_qp] = rotateVector(_traction_local[_qp], RotationGlobal2Local, /*inverse =*/true);
  _traction_spatial_derivatives[_qp] = rotateTensor2(
      _traction_spatial_derivatives_local[_qp], RotationGlobal2Local, /*inverse =*/true);

  _residual[_qp] = _traction[_qp];
  _jacobian[_qp] = _traction_spatial_derivatives[_qp];
  if (_effective_jump[_qp] > _max_effective_jump_old[_qp])
  {
    std::cout << "update _max_effective_jump" << std::endl;
    _max_effective_jump[_qp] = _effective_jump[_qp];
  }

  if (_current_elem->id() == 1234 && _qp == 0)
  {
    std::cout << "ELEM: " << _current_elem->id() << " SIDE:" << _current_side << " QP: " << _qp
              << std::endl;
    std::cout << "    _normals: " << _normals[_qp](0) << " " << _normals[_qp](1) << " "
              << _normals[_qp](2) << " " << std::endl;
    std::cout << "UNLOAD: " << checkLoadUnload() << std::endl;

    std::cout << "_max_effective_jump_old: " << _max_effective_jump_old[_qp] << std::endl;
    std::cout << "_max_effective_jump: " << _max_effective_jump[_qp] << std::endl;
    std::cout << "_effective_jump: " << _effective_jump[_qp] << std::endl;
    // std::cout << "    _traction_local: " << _traction_local[_qp][0] << " "
    //           << _traction_local[_qp][1] << " " << _traction_local[_qp][2] << " " << std::endl;
    //
    // std::cout << "    _traction_derivatives[0][:]: "
    //           << _traction_spatial_derivatives_local[_qp][0][0] << " "
    //           << _traction_spatial_derivatives_local[_qp][0][1] << " "
    //           << _traction_spatial_derivatives_local[_qp][0][2] << " " << std::endl;
    // std::cout << "    _traction_derivatives[1][:]: "
    //           << _traction_spatial_derivatives_local[_qp][1][0] << " "
    //           << _traction_spatial_derivatives_local[_qp][1][1] << " "
    //           << _traction_spatial_derivatives_local[_qp][1][2] << " " << std::endl;
    // std::cout << "    _traction_derivatives[2][:]: "
    //           << _traction_spatial_derivatives_local[_qp][2][0] << " "
    //           << _traction_spatial_derivatives_local[_qp][2][1] << " "
    //           << _traction_spatial_derivatives_local[_qp][2][2] << " " << std::endl;
    // std::cout << "    RotationGlobal2Local: " << RotationGlobal2Local(0, 0) << " "
    //           << RotationGlobal2Local(0, 1) << " " << RotationGlobal2Local(0, 2) << " "
    //           << std::endl;
    // std::cout << "                          " << RotationGlobal2Local(1, 0) << " "
    //           << RotationGlobal2Local(1, 1) << " " << RotationGlobal2Local(1, 2) << " "
    //           << std::endl;
    // std::cout << "                          " << RotationGlobal2Local(2, 0) << " "
    //           << RotationGlobal2Local(2, 1) << " " << RotationGlobal2Local(2, 2) << " "
    //           << std::endl;
    // std::cout << "    _displacement_jump: " << _displacement_jump[_qp][0] << " "
    //           << _displacement_jump[_qp][1] << " " << _displacement_jump[_qp][2] << " "
    //           << std::endl;
    // std::cout << "    _displacement_jump_local: " << _displacement_jump_local[_qp][0] << " "
    //           << _displacement_jump_local[_qp][1] << " " << _displacement_jump_local[_qp][2] << "
    //           "
    //           << std::endl;
    // // std::cout << "    _traction_spatial_derivatives_local: " <<
    // // _traction_spatial_derivatives_local << std::endl; std::cout << "
    // // _traction_spatial_derivatives: " <<  _traction_spatial_derivatives << std::endl;
    // std::cout << "    _residual: " << _residual[_qp][0] << " " << _residual[_qp][1] << " "
    //           << _residual[_qp][2] << " " << std::endl;
    // // std::cout << "    _jacobian: " << _jacobian << std::endl;
  }
}

std::vector<Real>
czmUOMatStateful::computeTractionLocal()
{

  std::vector<Real> Tlocal(3, 0);
  _effective_jump[_qp] = computeEffectiveJump();
  _effective_traction[_qp] = computeEffectiveTraction();

  for (unsigned int i = 0; i < 3; i++)
  {
    Real d_jump = _displacement_jump_local[_qp][i];
    if (i > 0)
      d_jump *= _beta * _beta;
    Tlocal[i] = _effective_traction[_qp] * d_jump;
  }

  return Tlocal;
}

std::vector<std::vector<Real>>
czmUOMatStateful::computeTractionSpatialDerivativeLocalNonLinear()
{
  std::vector<std::vector<Real>> TractionDerivativeLocalNonLinear(3, std::vector<Real>(3, 0));
  // Non linear
  // T_i/Delat_i = T_eff/Delta*(dij*Delta_i[B^2] - (Delta_i[B^2])^2/Delta_eff
  // T_eff/Delta is returned by computeEffectiveTractionNonLinear
  // dij Kronecker delta

  // non linear case
  for (unsigned int i = 0; i < 3; i++) /*i is jump component*/
  {
    for (unsigned int j = 0; j < 3; j++) /*j is partial derivatives component*/
    {
      Real temp_diag = 0;

      // diagonal term only present when  i == j
      if (i == j)
      {
        temp_diag = 1;
        if (j > 0)
          temp_diag *= _beta * _beta;
      }
      // off diagonal always present
      Real temp_off_diag = 0;
      if (_effective_jump[_qp] != 0)
      {
        temp_off_diag -= _displacement_jump_local[_qp][i] / (_deltaU0 * _effective_jump[_qp]);
        if (i > 0)
          temp_off_diag *= _beta * _beta;

        temp_off_diag *= _displacement_jump_local[_qp][j];
        if (j > 0)
          temp_off_diag *= _beta * _beta;
      }
      TractionDerivativeLocalNonLinear[i][j] =
          _effective_traction[_qp] * (temp_diag + temp_off_diag);
    }
  }

  return TractionDerivativeLocalNonLinear;
}

std::vector<std::vector<Real>>
czmUOMatStateful::computeTractionSpatialDerivativeLocalLinear()
{
  std::vector<std::vector<Real>> TractionDerivativeLocaLinear(3, std::vector<Real>(3, 0));
  // Linear
  // T_i/Delat_i = T_eff/Delta*(dij*Delta_i[B^2] - (Delta_i[B^2])^2/Delta_eff
  // T_eff/Delta is returned by computeEffectiveTractionNonLinear
  // dij Kronecker delta

  Real max_traction = computeEffectiveTractionNonLinear(_max_effective_jump_old[_qp]);
  // non linear case
  for (unsigned int i = 0; i < 3; i++) /*i is jump component*/
  {
    TractionDerivativeLocaLinear[i][i] = max_traction;
    if (i > 0)
      TractionDerivativeLocaLinear[i][i] *= _beta * _beta;
  }

  return TractionDerivativeLocaLinear;
}

std::vector<std::vector<Real>>
czmUOMatStateful::computeTractionSpatialDerivativeLocal()
{
  std::vector<std::vector<Real>> TractionDerivativeLocal(3, std::vector<Real>(3, 0));

  // non linear case
  if (!checkLoadUnload()) /* loading */
  {
    if (_max_effective_jump_old[_qp] < _deltaU0) /* not damaged*/
      return computeTractionSpatialDerivativeLocalNonLinear();
    else /*damaged*/
    {
      if (_effective_jump[_qp] <= _max_effective_jump_old[_qp]) /* linear region*/
        return computeTractionSpatialDerivativeLocalLinear();
      else
        return computeTractionSpatialDerivativeLocalNonLinear(); /* NON linear region*/
    }
  }
  else /* unload */
  {
    return computeTractionSpatialDerivativeLocalLinear();
  }
}

std::vector<Real>
czmUOMatStateful::rotateVector(const std::vector<Real> v,
                               const RealTensorValue R,
                               const bool inverse /*= false*/)
{
  RealTensorValue R_loc = R;
  if (inverse)
    R_loc = R_loc.transpose();

  std::vector<Real> vrot(3, 0);

  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      vrot[i] += v[j] * R_loc(i, j);
  return vrot;
}

std::vector<std::vector<Real>>
czmUOMatStateful::rotateTensor2(const std::vector<std::vector<Real>> T,
                                const RealTensorValue R,
                                const bool inverse /*= false*/)
{
  RealTensorValue R_loc = R;
  if (inverse)
    R_loc = R_loc.transpose();

  std::vector<std::vector<Real>> trot(3, std::vector<Real>(3, 0));
  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      for (unsigned int k = 0; k < 3; k++)
        for (unsigned int l = 0; l < 3; l++)
          trot[i][j] += R_loc(i, k) * R_loc(j, l) * T[k][l];
  return trot;
}
