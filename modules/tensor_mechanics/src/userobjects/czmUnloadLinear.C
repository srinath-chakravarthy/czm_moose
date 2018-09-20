//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "czmUnloadLinear.h"
#include "Material.h"
#include "MooseError.h"

registerMooseObject("TensorMechanicsApp", czmUnloadLinear);
template <>
InputParameters
validParams<czmUnloadLinear>()
{
  InputParameters params = validParams<CZMTractionSeparationUOBase>();
  params.addClassDescription("Simple Exponential cohseive law model, with damage");
  return params;
}

czmUnloadLinear::czmUnloadLinear(const InputParameters & parameters)
  : CZMTractionSeparationUOBase(parameters)
{
}

std::vector<Real>
czmUnloadLinear::computeTractionLocal(unsigned int qp) const
{
  std::vector<Real> TractionLocal(3, 0);

  Real T = _max_effective_traction_old[qp][0] / _max_effective_jump_old[qp][0];

  for (unsigned int i = 0; i < 3; i++)
  {
    Real d_jump = _displacement_jump[qp][i];
    if (i > 0)
      d_jump *= _beta * _beta;
    TractionLocal[i] = T * d_jump;
  }

  return TractionLocal;
}


std::vector<std::vector<Real>>
CohesiveLaw_Exponential::computeTractionSpatialDerivativeLocal(unsigned int qp) const
{
  std::vector<std::vector<Real>> TractionDerivativeLocal(3, std::vector<Real>(3, 0));
  TractionDerivativeLocal = getTractionSpatialDerivativeNonLinear(qp);
  return TractionDerivativeLocal;
}








std::vector<Real>
czmUnloadLinear::getTractionLinear(unsigned int qp) const
{

}

Real
czmUnloadLinear::getEffectiveTractionLinear(unsigned int qp) const
{
  Real effective_traction_l = 0;
  Real current_effective_jump = getEffectiveJump(qp);
  if (current_effective_jump != 0)
  {
    effective_traction_l = _max_effective_traction_old[qp][0] * current_effective_jump;
  }
  return effective_traction_l;
}
