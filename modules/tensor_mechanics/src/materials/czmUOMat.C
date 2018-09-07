//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "czmUOMat.h"
// #include "MooseMesh.h"
#include "RotationMatrix.h"

registerMooseObject("TensorMechanicsApp", czmUOMat);

template <>
InputParameters
validParams<czmUOMat>()
{
  InputParameters params = validParams<Material>();

  params.addRequiredParam<UserObjectName>(
      "displacement_jump_UO", "the name of the computing material property across an interface");

  params.addRequiredParam<std::vector<Real>>(
      "DeltaU0",
      "a vector containing the displacement value at which maximum"
      " traction occurs for the normal(1st) and tangential(2nd) "
      " direction.");
  params.addRequiredParam<std::vector<Real>>("MaxAllowableTraction",
                                             "a vector containing the maximum allowed traction"
                                             "for the normal(1st) and tangential(2nd) direction.");
  params.addClassDescription("this material class is used when defining a "
                             "cohesive zone model to store stafeul properties");
  return params;
}

czmUOMat::czmUOMat(const InputParameters & parameters)
  : Material(parameters),
    _deltaU0(getParam<std::vector<Real>>("DeltaU0")),
    _maxAllowableTraction(getParam<std::vector<Real>>("MaxAllowableTraction")),
    _displacement_jump_UO(getUserObject<DispJumpUO_QP>("displacement_jump_UO")),
    _displacement_jump(declareProperty<std::vector<Real>>("displacement_jump")),
    _displacement_jump_local(declareProperty<std::vector<Real>>("displacement_jump_local")),
    _traction(declareProperty<std::vector<Real>>("traction")),
    _traction_local(declareProperty<std::vector<Real>>("traction_local")),
    _traction_spatial_derivatives(
        declareProperty<std::vector<std::vector<Real>>>("traction_spatial_derivatives")),
    _traction_spatial_derivatives_local(
        declareProperty<std::vector<std::vector<Real>>>("traction_spatial_derivatives_local"))

{
}

void
czmUOMat::computeQpProperties()
{

  _displacement_jump[_qp].resize(3, 0);
  _displacement_jump_local[_qp].resize(3, 0);
  _traction[_qp].resize(3, 0);
  _traction_local[_qp].resize(3, 0);
  _traction_spatial_derivatives[_qp].resize(3, std::vector<Real>(3, 0));
  _traction_spatial_derivatives_local[_qp].resize(3, std::vector<Real>(3, 0));

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
}

std::vector<Real>
czmUOMat::rotateVector(const std::vector<Real> v,
                       const RealTensorValue R,
                       const bool inverse /*= false*/)
{
  RealTensorValue R_loc = R;
  if (inverse)
    R_loc = R_loc.transpose();

  std::vector<Real> vrot(3, 0);

  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      vrot[i] += v[i] * R_loc(i, j);
  return vrot;
}

std::vector<std::vector<Real>>
czmUOMat::rotateTensor2(const std::vector<std::vector<Real>> T,
                        const RealTensorValue R,
                        const bool inverse /*= false*/)
{
  RealTensorValue R_loc = R;
  if (inverse)
    R_loc = R_loc.transpose();

  std::vector<std::vector<Real>> trot(3, std::vector<Real>(3, 0));
  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int k = 0; k < 3; k++)
      for (unsigned int j = 0; j < 3; j++)
        for (unsigned int l = 0; l < 3; l++)
          trot[i][k] += T[k][l] * R_loc(i, j) * R_loc(j, l);
  return trot;
}

std::vector<Real>
czmUOMat::computeTractionLocal()
{
  std::vector<Real> Tlocal(3, 0);
  // convention N, T, S
  Real temp, X, expX, A_i, B_i;
  X = 0;
  for (unsigned int i = 0; i < 3; i++)
  {
    temp = _displacement_jump_local[_qp][i] / _deltaU0[i];
    if (i > 0)
    {
      temp *= temp; // square for shear component
    };
    X += temp;
  }
  expX = std::exp(-X);
  for (unsigned int i = 0; i < 3; i++)
  {
    if (i == 0)
    {
      temp = std::exp(1);
    }
    else
    {
      temp = std::sqrt(2 * std::exp(1));
    }
    A_i = _maxAllowableTraction[i] * temp;
    B_i = _displacement_jump_local[_qp][i] / _deltaU0[i];
    Tlocal[i] = A_i * B_i * expX;
  }
  return Tlocal;
}

std::vector<std::vector<Real>>
czmUOMat::computeTractionSpatialDerivativeLocal()
{
  std::vector<std::vector<Real>> TractionDerivativeLocal(3, std::vector<Real>(3, 0));
  // this function compute partial derivates of Tn[0][:], Tt[1][:], Ts[2][:]
  // w.r.t. dun, dut, dus
  // T_i = A_i*B_i*exp(-X) with:
  // A_i = \sigma_i,max * (\alpha_i*e)^{1/\alpha_i} with \alpha_i = 1 for i==n
  // \alpha_i = 2 for i!=n
  // B_i = \delta_u,i / \delta_0,i
  // X = sum_i=1^3{(\delta_u,i / \delta_0,i)^\alpha_i}  with \alpha_i = 1 for i==n
  // \alpha_i = 2 for i!=n
  // dTi_duj = A_i * ( dBi_duj * exp(-X) + B_i * exp(-X) * dX_duj  )
  //         = A_i * ( exp(-X) * (dBi_duj + B_i * dX_duj ) )
  // convention N, T, S
  unsigned int i, j;
  Real expX, temp, X;
  // compute X and the exponential term
  temp = 0;
  X = 0;
  for (i = 0; i < 3; i++)
  {
    temp = _displacement_jump_local[_qp][i] / _deltaU0[i];
    if (i > 0)
      temp *= temp;
    X += temp;
  }
  expX = std::exp(-X);
  // compute partial derivatives in local coordaintes w.r.t. the master surface siplacement
  //            | dTn/dun dTn/dut dTn/dus |
  // dTi_duj  = | dTt/dun dTt/dut dTt/dus | = _TractionDerivativeLocal[i][j]
  //            | dTs/dun dTs/dut dTs/dus |
  Real A_i, B_i;
  Real dBi_dui, dX_duj;
  for (i = 0; i < 3; i++)
  {
    // compute A_i
    if (i == 0) // alpha = 1
      A_i = std::exp(1);
    else // alpha = 2
      A_i = std::sqrt(2 * std::exp(1));
    A_i *= _maxAllowableTraction[i];
    // compute B_i
    B_i = _displacement_jump_local[_qp][i] / _deltaU0[i];
    for (j = 0; j < 3; j++)
    {
      // add term for diagonal entry dBi_dui
      dBi_dui = 0;
      if (i == j)
      {
        dBi_dui = 1 / _deltaU0[j];
      }
      // compute the derivative of the argument of exponential
      if (j == 0) // alpha = 1
        dX_duj = 1. / _deltaU0[j];
      else // alpha = 2
        dX_duj = 2. * _displacement_jump_local[_qp][i] / (_deltaU0[j] * _deltaU0[j]);
      TractionDerivativeLocal[i][j] =
          A_i * expX * (dBi_dui - B_i * dX_duj); // the minus sign is due to exp(-X)
    }
  }
  return TractionDerivativeLocal;
}
