//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CZMUOBasedMaterial.h"
#include "RotationMatrix.h"

registerMooseObject("TensorMechanicsApp", CZMUOBasedMaterial);

template <>
InputParameters
validParams<CZMUOBasedMaterial>()
{
  InputParameters params = validParams<Material>();

  params.addRequiredParam<UserObjectName>(
      "traction_separation_UO",
      "the name of the user object including the traction separation law");
  params.addRequiredParam<UserObjectName>(
      "displacement_jump_UO",
      "the name of the UO collecting bulk material property across an interface");
  params.addClassDescription("this material class is used when defining a "
                             "cohesive zone model");
  return params;
}

CZMUOBasedMaterial::CZMUOBasedMaterial(const InputParameters & parameters)
  : Material(parameters),
    _displacement_jump_UO(getUserObject<DispJumpUO_QP>("displacement_jump_UO")),
    _traction_separation_UO(getUserObject<CZMTractionSeparationUOBase>("traction_separation_UO")),
    _displacement_jump(declareProperty<std::vector<Real>>("displacement_jump")),
    _displacement_jump_local(declareProperty<std::vector<Real>>("displacement_jump_local")),
    _displacement_jump_local_old(
        getMaterialPropertyOld<std::vector<Real>>("displacement_jump_local")),
    _traction(declareProperty<std::vector<Real>>("traction")),
    _traction_local(declareProperty<std::vector<Real>>("traction_local")),
    _traction_spatial_derivatives(
        declareProperty<std::vector<std::vector<Real>>>("traction_spatial_derivatives")),
    _traction_spatial_derivatives_local(
        declareProperty<std::vector<std::vector<Real>>>("traction_spatial_derivatives_local")),
    _czm_residual(declareProperty<std::vector<Real>>("czm_residual")),
    _czm_jacobian(declareProperty<std::vector<std::vector<Real>>>("czm_jacobian")),
    _n_uo_czm_properties(_traction_separation_UO.getNumberStatefulMaterialProperies())

{

  if (_n_uo_czm_properties > 0)
  {
    // initialize the userobject material property container
    _uo_czm_properties.resize(_n_uo_czm_properties);
    _uo_czm_properties_old.resize(_n_uo_czm_properties);
    for (unsigned int mp_index = 0; mp_index < _n_uo_czm_properties; mp_index++)
    {
      // declare a material property
      _uo_czm_properties[mp_index] = &declareProperty<std::vector<Real>>(
          _traction_separation_UO.getStatefulMaterialProperyName(mp_index));
      // for a stateful material property get the old value
      _uo_czm_properties_old[mp_index] = &getMaterialPropertyOld<std::vector<Real>>(
          _traction_separation_UO.getStatefulMaterialProperyName(mp_index));
    }
  }
}

void
CZMUOBasedMaterial::computeQpProperties()
{
  _displacement_jump[_qp].resize(3, 0);
  _displacement_jump_local[_qp].resize(3, 0);
  _traction[_qp].resize(3, 0);
  _traction_local[_qp].resize(3, 0);
  _traction_spatial_derivatives[_qp].resize(3, std::vector<Real>(3, 0));
  _traction_spatial_derivatives_local[_qp].resize(3, std::vector<Real>(3, 0));

  _czm_residual[_qp].resize(3, 0);
  _czm_jacobian[_qp].resize(3, std::vector<Real>(3, 0));

  RealTensorValue RotationGlobal2Local =
      RotationMatrix::rotVec1ToVec2(_normals[_qp], RealVectorValue(1, 0, 0));

  _displacement_jump[_qp] =
      _displacement_jump_UO.getDisplacementJump(_current_elem->id(), _current_side, _qp);

  _displacement_jump_local[_qp] = rotateVector(_displacement_jump[_qp], RotationGlobal2Local);

  _traction_local[_qp] = _traction_separation_UO.computeTractionLocal(_qp);
  _traction_spatial_derivatives_local[_qp] =
      _traction_separation_UO.computeTractionSpatialDerivativeLocal(_qp);

  if (_n_uo_czm_properties > 0)
    for (unsigned int mp_index = 0; mp_index < _n_uo_czm_properties; mp_index++)
      (*_uo_czm_properties[mp_index])[_qp] =
          _traction_separation_UO.getNewStatefulMaterialProperty(mp_index, _qp);

  _traction[_qp] = rotateVector(_traction_local[_qp], RotationGlobal2Local, /*inverse =*/true);
  _traction_spatial_derivatives[_qp] = rotateTensor2(
      _traction_spatial_derivatives_local[_qp], RotationGlobal2Local, /*inverse =*/true);

  _czm_residual[_qp] = _traction[_qp];
  _czm_jacobian[_qp] = _traction_spatial_derivatives[_qp];
}

void
CZMUOBasedMaterial::initQpStatefulProperties()
{
  RealTensorValue RotationGlobal2Local =
      RotationMatrix::rotVec1ToVec2(_normals[_qp], RealVectorValue(1, 0, 0));

  _displacement_jump[_qp].resize(3, 0);
  _displacement_jump_local[_qp].resize(3, 0);

  // _displacement_jump[_qp] =
  //     _displacement_jump_UO.getDisplacementJump(_current_elem->id(), _current_side, _qp);
  // /// is there a special procedure for reload?
  // _displacement_jump_local[_qp] = rotateVector(_displacement_jump[_qp], RotationGlobal2Local);

  /// is there a special procedure for reload?
  if (_n_uo_czm_properties > 0)
  {
    for (unsigned int mp_index = 0; mp_index < _n_uo_czm_properties; mp_index++)
      (*_uo_czm_properties[mp_index])[_qp] =
          _traction_separation_UO.getStatefulMaterialProperysIntialValues(mp_index);
  }
}

std::vector<Real>
CZMUOBasedMaterial::rotateVector(const std::vector<Real> v,
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
CZMUOBasedMaterial::rotateTensor2(const std::vector<std::vector<Real>> T,
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
