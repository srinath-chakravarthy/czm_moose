//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef COHESIVELAW_EXPONENTIAL_H
#define COHESIVELAW_EXPONENTIAL_H

#include "CZMTractionSeparationUOBase.h"

class CohesiveLaw_Exponential;

template <>
InputParameters validParams<CohesiveLaw_Exponential>();

/**
Traction sepration law basic user object
 */
class CohesiveLaw_Exponential : public CZMTractionSeparationUOBase
{
public:
  CohesiveLaw_Exponential(const InputParameters & parameters);

  std::vector<Real> computeTractionLocal(unsigned int qp) const override;
  std::vector<std::vector<Real>>
  computeTractionSpatialDerivativeLocal(unsigned int qp) const override;

  std::vector<Real> getNewStatefulMaterialProperty(unsigned int qp,
                                                   unsigned int mp_index) const override;

protected:
  // cohesive law parameters
  const Real _displacement_jump_peak;
  const Real _traction_peak;
  const Real _beta;
  const MaterialProperty<std::vector<Real>> & _max_effective_jump_old;
  const MaterialProperty<std::vector<Real>> & _max_effective_traction_old;

  Real getEffectiveJump(unsigned int /*qp*/) const;
  Real getEffectiveTractionNonLinear(unsigned int /*qp*/) const;
  Real getEffectiveTractionLinear(unsigned int /*qp*/) const;
  std::vector<Real> getTractionNonLinear(unsigned int /*qp*/) const;
  std::vector<Real> getTractionLinear(unsigned int /*qp*/) const;
};

#endif // COHESIVELAW_EXPONENTIAL_H
