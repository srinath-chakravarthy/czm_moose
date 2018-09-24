//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CZMUnloadLinear.h"
#include "Material.h"
#include "MooseError.h"

registerMooseObject("TensorMechanicsApp", CZMUnloadLinear);
template <>
InputParameters
validParams<CZMUnloadLinear>()
{
  InputParameters params = validParams<CZMTractionSeparationUOBase>();
  params.addParam<unsigned int>("n_stateful_mp", 0, "number of stateful material properties");
  params.addParam<unsigned int>(
      "n_non_stateful_mp", 0, "number of NON-stateful material properties");
  params.addParam<std::string>("max_effective_jump_mp_name",
                               "max_effective_jump",
                               "name fo the mp representing the effective jump");
  params.addParam<std::string>("max_linear_traction_mp_name",
                               "max_effective_traction",
                               "name fo the mp representing the maximum linear traction_peak");
  params.addParam<std::string>("weighted_displacement_jump_mp_name",
                               "weighted_displacement_jump",
                               "name fo the mp representing the weighted dispalcement jump");
  params.addParam<std::string>("displacement_jump_weights_mp_name",
                               "displacement_jump_weights",
                               "name fo the mp representing the weighted applied to each component "
                               "of the displacement jump (come from non linear cohesive law)");
  params.addClassDescription("Linear unloading for damged material");
  return params;
}

CZMUnloadLinear::CZMUnloadLinear(const InputParameters & parameters)
  : CZMTractionSeparationUOBase(parameters),
    _max_effective_jump(getMaterialPropertyOldByName<std::vector<Real>>(
        getParam<std::string>("max_effective_jump_mp_name"))),
    _max_effective_traction(getMaterialPropertyOldByName<std::vector<Real>>(
        getParam<std::string>("max_linear_traction_mp_name"))),
    _weighted_displacement_jump(getMaterialPropertyByName<std::vector<Real>>(
        getParam<std::string>("weighted_displacement_jump_mp_name"))),
    _displacement_jump_weights(getMaterialPropertyByName<std::vector<Real>>(
        getParam<std::string>("displacement_jump_weights_mp_name")))

{
}

RealVectorValue
CZMUnloadLinear::computeTractionLocal(unsigned int qp) const
{
  RealVectorValue TractionLocal;

  Real T = 0;
  if (_max_effective_jump[qp][0] > 0)
  {
    T = _max_effective_traction[qp][0] / _max_effective_jump[qp][0];
    for (unsigned int i = 0; i < 3; i++)
      TractionLocal(i) = T * _weighted_displacement_jump[qp][i];
  }
  else
    mooseError("CZMUnloadLinear:: should not be called if call _max_effective_jump ==0 ");
  return TractionLocal;
}

RankTwoTensor
CZMUnloadLinear::computeTractionSpatialDerivativeLocal(unsigned int qp) const
{
  RankTwoTensor TractionDerivativeLocal;

  if (_max_effective_jump[qp][0] > 0)
  {
    Real T = _max_effective_traction[qp][0] / _max_effective_jump[qp][0];
    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int j = 0; j < 3; j++)
        if (i == j)
          TractionDerivativeLocal(i, j) = T * _displacement_jump_weights[qp][i];
        else
          TractionDerivativeLocal(i, j) = 0;
  }
  else
    mooseError("CZMUnloadLinear:: should not be called if call _max_effective_jump ==0 ");

  return TractionDerivativeLocal;
}
