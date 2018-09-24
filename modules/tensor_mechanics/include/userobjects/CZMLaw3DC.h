//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef CZMLAW3DC_H
#define CZMLAW3DC_H

#include "CZMTractionSeparationUOBase.h"

class CZMLaw3DC;

template <>
InputParameters validParams<CZMLaw3DC>();

/**
Traction sepration law basic user object
 */
class CZMLaw3DC : public CZMTractionSeparationUOBase
{
public:
  CZMLaw3DC(const InputParameters & parameters);

  RealVectorValue computeTractionLocal(unsigned int qp) const override;
  RankTwoTensor computeTractionSpatialDerivativeLocal(unsigned int qp) const override;

protected:
  // cohesive law parameters
  const std::vector<Real> _deltaU0;
  const std::vector<Real> _maxAllowableTraction;
};

#endif // CZMLAW3DC_H
