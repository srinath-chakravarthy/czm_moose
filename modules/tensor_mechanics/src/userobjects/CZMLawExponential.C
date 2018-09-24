//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CZMLawExponential.h"
#include "Material.h"
#include "MooseError.h"

registerMooseObject("TensorMechanicsApp", CZMLawExponential);
template <>
InputParameters
validParams<CZMLawExponential>()
{
  InputParameters params = validParams<CZMTractionSeparationUOBase>();
  params.addClassDescription("Simple Exponential cohseive law model, with damage");
  params.addParam<unsigned int>("n_stateful_mp", 3, "number of stateful material properties");
  params.addParam<std::vector<std::string>>(
      "stateful_mp_names",
      std::vector<std::string>{"effective_jump", "max_effective_jump", "max_effective_traction"},
      "name of stateful material properties");
  params.addParam<std::vector<unsigned int>>("stateful_mp_sizes",
                                             std::vector<unsigned int>{1, 1, 1},
                                             "size of each stateful material properties");
  params.addParam<std::vector<Real>>("stateful_mp_initial_values",
                                     std::vector<Real>{0, 0, 0},
                                     "intial material proeprties values");
  params.addParam<unsigned int>(
      "n_non_stateful_mp", 3, "number of NON-stateful material properties");
  params.addParam<std::vector<std::string>>("non_stateful_mp_names",
                                            std::vector<std::string>{"weighted_displacement_jump",
                                                                     "effective_traction_nl",
                                                                     "displacement_jump_weights"},
                                            "name of NON stateful material properties");
  params.addParam<std::vector<unsigned int>>("non_stateful_mp_sizes",
                                             std::vector<unsigned int>{3, 1, 3},
                                             "size of each stateful material properties");

  params.addRequiredParam<Real>(
      "displacement_jump_peak",
      "the value of effective displacement jump at wich peak traction occurs");
  params.addRequiredParam<Real>("traction_peak", "the value of peak effective traction");
  params.addParam<Real>("beta", 0.5, "coefficinet weighting the effect of shear displacement jump");
  return params;
}

CZMLawExponential::CZMLawExponential(const InputParameters & parameters)
  : CZMTractionSeparationUOBase(parameters),
    _displacement_jump_peak(getParam<Real>("displacement_jump_peak")),
    _traction_peak(getParam<Real>("traction_peak")),
    _beta(getParam<Real>("beta")),
    _effective_jump(getMaterialPropertyByName<std::vector<Real>>(_stateful_mp_names[0])),
    _effective_jump_old(getMaterialPropertyOldByName<std::vector<Real>>(_stateful_mp_names[0])),
    _max_effective_jump_old(getMaterialPropertyOldByName<std::vector<Real>>(_stateful_mp_names[1])),
    _max_effective_traction_old(
        getMaterialPropertyOldByName<std::vector<Real>>(_stateful_mp_names[2])),
    _weighted_displacement_jump(
        getMaterialPropertyByName<std::vector<Real>>(_non_stateful_mp_names[0])),
    _effective_traction(getMaterialPropertyByName<std::vector<Real>>(_non_stateful_mp_names[1]))
{
}

///
unsigned int
CZMLawExponential::checkLoadUnload(const unsigned int qp) const
{
  // first step _max_effective_jump_old ==0
  if (_max_effective_jump_old[qp][0] == 0)
  {
    if (_displacement_jump[qp](0) >= 0)
      return 0;
    else
      return 2;
  }
  // no copenatration
  if (_displacement_jump[qp](0) >= 0)
  {
    // loading
    if (_effective_jump[qp][0] >= _effective_jump_old[qp][0])
    {
      if (_effective_jump[qp][0] >= _max_effective_jump_old[qp][0]) // no damage
        return 0;                                                   // select non linear CZM
      else
        return 1; // linear
                  // check for damage
    }
    // unloading
    else
      return 1; // linear
  }
  else // copenetration
    return 2;
}

RealVectorValue
CZMLawExponential::computeTractionLocal(unsigned int qp) const
{
  RealVectorValue TractionLocal;

  for (unsigned int i = 0; i < 3; i++)
    TractionLocal(i) = _effective_traction[qp][0] * _weighted_displacement_jump[qp][i];

  return TractionLocal;
}

RankTwoTensor
CZMLawExponential::computeTractionSpatialDerivativeLocal(unsigned int qp) const
{

  RankTwoTensor TractionDerivativeLocal;
  Real T_eff = _effective_traction[qp][0];
  Real D_eff = _effective_jump[qp][0];
  Real beta2 = std::pow(_beta, 2);

  if (D_eff > 0)
  {
    Real D_eff_D_p = D_eff * _displacement_jump_peak;

    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int j = 0; j < 3; j++)
      {
        Real diag_term = 0;
        if (i == j)
        {
          diag_term += 1;
          if (i > 0)
            diag_term *= beta2;
        }

        Real offdiag_term =
            _weighted_displacement_jump[qp][i] * _weighted_displacement_jump[qp][j] / D_eff_D_p;

        TractionDerivativeLocal(i, j) = T_eff * (diag_term - offdiag_term);
      }
  }
  else if (D_eff == 0)
  {
    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int j = 0; j < 3; j++)
      {
        if (i == j)
        {
          TractionDerivativeLocal(i, j) = T_eff;
          if (i > 0)
            TractionDerivativeLocal(i, j) *= beta2;
        }
        else
        {
          TractionDerivativeLocal(i, j) = 0;
        }
      }
  }
  return TractionDerivativeLocal;
}

Real
CZMLawExponential::getEffectiveJump(unsigned int qp) const
{
  Real effective_jump = 0;
  for (unsigned int i = 0; i < 3; i++)
  {
    Real temp = _displacement_jump[qp](i) * _displacement_jump[qp](i);
    if (i > 0)
      temp *= _beta * _beta;
    effective_jump += temp;
  }
  return std::sqrt(effective_jump);
}

Real
CZMLawExponential::getEffectiveTraction(unsigned int qp) const
{
  /// T_eff =e*Tp/Dp*exp(-D_eff/D_p)
  Real effective_traction_nl = 0;

  if (_effective_jump[qp][0] != 0)
  {
    Real d_norm = _effective_jump[qp][0] / _displacement_jump_peak;
    effective_traction_nl =
        std::exp(1) * _traction_peak / _displacement_jump_peak * std::exp(-1 * d_norm);
  }
  return effective_traction_nl;
}

std::vector<Real>
CZMLawExponential::getNewNonStatefulMaterialProperty(unsigned int qp, unsigned int mp_index) const
{
  std::vector<Real> temp(getNonStatefulMaterialPropertySize(mp_index), 0);
  if (mp_index == 0) /*WEIGHTED DISPLACEMENT JUMP*/
  {
    for (unsigned int i = 0; i < 3; i++)
    {
      temp[i] = _displacement_jump[qp](i);
      if (i > 0)
        temp[i] *= _beta * _beta;
    }
  }
  else if (mp_index == 1) /*effective traction*/
    temp[0] = getEffectiveTraction(qp);
  else if (mp_index == 2) /*DISPALCEMENT JUMP WEIGTHS*/
    for (unsigned int i = 0; i < 3; i++)
    {
      temp[i] = 1;
      if (i > 0)
        temp[i] *= _beta * _beta;
    }

  return temp;
}

std::vector<Real>
CZMLawExponential::getNewStatefulMaterialProperty(unsigned int qp, unsigned int mp_index) const
{
  std::vector<Real> temp(getStatefulMaterialPropertySize(mp_index), 0);
  if (mp_index == 0) /*EFFECTIVE JUMP*/
    temp[0] = getEffectiveJump(qp);
  else if (mp_index == 1) /*MAX EFFECTIVE JUMP*/
  {
    temp[0] = _max_effective_jump_old[qp][0];
    if (_effective_jump[qp][0] > _max_effective_jump_old[qp][0] && _displacement_jump[qp](0) > 0)
      temp[0] = _effective_jump[qp][0]; // new maximum effective jump
  }
  else if (mp_index == 2) /*MAX EFFECTIVE TRACTION*/
  {
    temp[0] = _max_effective_traction_old[qp][0];
    if (_effective_jump[qp][0] > _max_effective_jump_old[qp][0] && _displacement_jump[qp](0) > 0)
      temp[0] = getEffectiveTraction(qp) * _effective_jump[qp][0];
  }

  return temp;
}
