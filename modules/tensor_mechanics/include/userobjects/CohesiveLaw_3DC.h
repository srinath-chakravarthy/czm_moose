//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef COHESIVELAW_3DC_H
#define COHESIVELAW_3DC_H

#include "CZMTractionSeparationUOBase.h"

class CohesiveLaw_3DC;

template <>
InputParameters validParams<CohesiveLaw_3DC>();

/**
Traction sepration law basic user object
 */
class CohesiveLaw_3DC : public CZMTractionSeparationUOBase
{
public:
  CohesiveLaw_3DC(const InputParameters & parameters);

  std::vector<Real> computeTractionLocal(unsigned int qp) const override;
  std::vector<std::vector<Real>>
  computeTractionSpatialDerivativeLocal(unsigned int qp) const override;

protected:
  // cohesive law parameters
  const std::vector<Real> _deltaU0;
  const std::vector<Real> _maxAllowableTraction;
};

#endif // COHESIVELAW_3DC_H
