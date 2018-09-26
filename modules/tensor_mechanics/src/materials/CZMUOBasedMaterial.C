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
  params.addParam<UserObjectName>(
      "unload_traction_separation_UO",
      "the name of the user object including the traction separation law");
  params.addParam<UserObjectName>(
      "coopenetration_penalty_UO",
      "the name of the user object including the traction separation law");
  params.addRequiredParam<UserObjectName>(
      "displacement_jump_UO",
      "the name of the UO collecting bulk material property across an interface");
  params.addParam<Real>("coopenetration_penalty", 1, "copenation penalty factor");
  params.addClassDescription("this material class is used when defining a "
                             "cohesive zone model");
  return params;
}

CZMUOBasedMaterial::CZMUOBasedMaterial(const InputParameters & parameters)
  : Material(parameters),
    _displacement_jump_UO(getUserObject<DispJumpAndNormalsUO_QP>("displacement_jump_UO")),
    _traction_separation_UO(getUserObject<CZMTractionSeparationUOBase>("traction_separation_UO")),
    _unload_traction_separation_UO(
        parameters.isParamSetByUser("unload_traction_separation_UO")
            ? getUserObject<CZMTractionSeparationUOBase>("unload_traction_separation_UO")
            : _traction_separation_UO),
    _coopenetration_penalty_UO(
        parameters.isParamSetByUser("coopenetration_penalty_UO")
            ? getUserObject<CZMTractionSeparationUOBase>("coopenetration_penalty_UO")
            : _traction_separation_UO),
    _coopenetration_penalty(getParam<Real>("coopenetration_penalty")),
    _selected_CZM_UO(&_unload_traction_separation_UO),
    _displacement_jump(declareProperty<RealVectorValue>("displacement_jump")),
    _displacement_jump_local(declareProperty<RealVectorValue>("displacement_jump_local")),
    _displacement_jump_local_old(
        getMaterialPropertyOld<RealVectorValue>("displacement_jump_local")),
    _traction(declareProperty<RealVectorValue>("traction")),
    _traction_local(declareProperty<RealVectorValue>("traction_local")),
    _traction_spatial_derivatives(declareProperty<RankTwoTensor>("traction_spatial_derivatives")),
    _traction_spatial_derivatives_local(
        declareProperty<RankTwoTensor>("traction_spatial_derivatives_local")),
    _czm_residual(declareProperty<RealVectorValue>("czm_residual")),
    _czm_jacobian(declareProperty<std::vector<std::vector<Real>>>("czm_jacobian")),
    _normals_MP(declareProperty<RealVectorValue>("normals_MP")),
    _normals_neighbor_MP(declareProperty<RealVectorValue>("normals_neighbor_MP")),
    _normals_average(declareProperty<RealVectorValue>("normals_average")),

    _uo_id(0),
    _n_uo_czm_properties(_traction_separation_UO.getNumberStatefulMaterialProperties()),
    _n_non_stateful_uo_czm_properties(
        _traction_separation_UO.getNumberNonStatefulMaterialProperties())

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
          _traction_separation_UO.getStatefulMaterialPropertyName(mp_index));
      // for a stateful material property get the old value
      _uo_czm_properties_old[mp_index] = &getMaterialPropertyOld<std::vector<Real>>(
          _traction_separation_UO.getStatefulMaterialPropertyName(mp_index));
    }
  }

  if (_n_non_stateful_uo_czm_properties > 0)
  {
    // initialize the userobject material property container
    _uo_non_stateful_czm_properties.resize(_n_non_stateful_uo_czm_properties);
    for (unsigned int mp_index = 0; mp_index < _n_non_stateful_uo_czm_properties; mp_index++)
    {
      // declare a material property
      _uo_non_stateful_czm_properties[mp_index] = &declareProperty<std::vector<Real>>(
          _traction_separation_UO.getNonStatefulMaterialPropertyName(mp_index));
    }
  }
}

void
CZMUOBasedMaterial::computeQpProperties()
{
  // resize non stateful mp
  if (_n_non_stateful_uo_czm_properties > 0)
    for (unsigned int mp_index = 0; mp_index < _n_non_stateful_uo_czm_properties; mp_index++)
      (*_uo_non_stateful_czm_properties[mp_index])[_qp].resize(
          _traction_separation_UO.getNonStatefulMaterialPropertySize(mp_index));

  _czm_jacobian[_qp].resize(3, std::vector<Real>(3, 0));
  _normals_MP[_qp] = _displacement_jump_UO.getNormalMaster(_current_elem->id(), _current_side, _qp);
  _normals_neighbor_MP[_qp] =
      _displacement_jump_UO.getNormalSlave(_current_elem->id(), _current_side, _qp);

  for (unsigned int i = 0; i < 3; i++)
  {
    _normals_average[_qp](i) = (_normals_MP[_qp](i) + _normals_neighbor_MP[_qp](i)) / 2;
    // if (_normals_MP[_qp](i) != _normals[_qp](i))
    // {
    //   std::cout << " normals are wrong: UO " << _normals_MP[_qp](i) << " MOOSE" <<
    //   _normals[_qp](i)
    //             << std::endl;
    // }
    // else
    //   std::cout << " normals_average " << _normals_average[_qp](i) << std::endl;
  }
  RealTensorValue RotationGlobal2Local =
      RotationMatrix::rotVec1ToVec2(_normals_average[_qp], RealVectorValue(1, 0, 0));

  _displacement_jump[_qp] =
      _displacement_jump_UO.getDisplacementJump(_current_elem->id(), _current_side, _qp);

  _displacement_jump_local[_qp] = rotateVector(_displacement_jump[_qp], RotationGlobal2Local);

  if (_n_uo_czm_properties > 0)
    for (unsigned int mp_index = 0; mp_index < _n_uo_czm_properties; mp_index++)
      (*_uo_czm_properties[mp_index])[_qp] =
          _traction_separation_UO.getNewStatefulMaterialProperty(_qp, mp_index);
  if (_n_non_stateful_uo_czm_properties > 0)
    for (unsigned int mp_index = 0; mp_index < _n_non_stateful_uo_czm_properties; mp_index++)
      (*_uo_non_stateful_czm_properties[mp_index])[_qp] =
          _traction_separation_UO.getNewNonStatefulMaterialProperty(_qp, mp_index);

  selectCzmUO();

  _traction_local[_qp] = _selected_CZM_UO->computeTractionLocal(_qp);
  //
  _traction_spatial_derivatives_local[_qp] =
      _selected_CZM_UO->computeTractionSpatialDerivativeLocal(_qp);

  if (_displacement_jump_local[_qp](0) < 0)
  {
    _traction_local[_qp](0) *= _coopenetration_penalty;
    for (unsigned int i = 0; i < 3; i++)
      _traction_spatial_derivatives_local[_qp](i, 0) *= _coopenetration_penalty;
  }

  _traction[_qp] = rotateVector(_traction_local[_qp], RotationGlobal2Local, /*inverse =*/true);
  _traction_spatial_derivatives[_qp] = rotateTensor2(
      _traction_spatial_derivatives_local[_qp], RotationGlobal2Local, /*inverse =*/true);

  _czm_residual[_qp] = _traction[_qp];
  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      _czm_jacobian[_qp][i][j] = _traction_spatial_derivatives[_qp](i, j);
}

void
CZMUOBasedMaterial::initQpStatefulProperties()
{
  RealTensorValue RotationGlobal2Local =
      RotationMatrix::rotVec1ToVec2(_normals[_qp], RealVectorValue(1, 0, 0));

  for (unsigned int i = 0; i < 3; i++)
  {
    _displacement_jump[_qp](i) = 0;
    _displacement_jump_local[_qp](i) = 0;
  }

  /// is there a special procedure for reload?
  if (_n_uo_czm_properties > 0)
  {
    for (unsigned int mp_index = 0; mp_index < _n_uo_czm_properties; mp_index++)
    {
      (*_uo_czm_properties[mp_index])[_qp].resize(
          _traction_separation_UO.getStatefulMaterialPropertySize(mp_index));

      (*_uo_czm_properties[mp_index])[_qp] =
          _traction_separation_UO.getStatefulMaterialPropertysIntialValues(mp_index);
    }
  }
}

RealVectorValue
CZMUOBasedMaterial::rotateVector(const RealVectorValue v,
                                 const RealTensorValue R,
                                 const bool inverse /*= false*/)
{
  RealTensorValue R_loc = R;
  if (inverse)
    R_loc = R_loc.transpose();

  RealVectorValue vrot;

  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      vrot(i) += v(j) * R_loc(i, j);
  return vrot;
}

RankTwoTensor
CZMUOBasedMaterial::rotateTensor2(const RankTwoTensor T,
                                  const RealTensorValue R,
                                  const bool inverse /*= false*/)
{
  RealTensorValue R_loc = R;
  if (inverse)
    R_loc = R_loc.transpose();

  RankTwoTensor trot = T;
  trot.rotate(R_loc);
  return trot;
}

void
CZMUOBasedMaterial::selectCzmUO()
{
  if (_n_uo_czm_properties > 0)
  {
    unsigned int _uo_id = _traction_separation_UO.checkLoadUnload(_qp);

    switch (_uo_id)
    {
      case (unsigned int)0:
        _selected_CZM_UO = &_traction_separation_UO;
        break;
      case (unsigned int)1:
        _selected_CZM_UO = &_unload_traction_separation_UO;
        break;
      case (unsigned int)2:
        _selected_CZM_UO = &_coopenetration_penalty_UO;
        break;
      default:
        mooseError("CZMUOBasedMaterial:: something is wrong, selecting wrong UO");
    }
  }
  else
    _selected_CZM_UO = &_traction_separation_UO;
}
