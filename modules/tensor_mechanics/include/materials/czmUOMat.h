//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef INTERFACEUOMATERIAL_H
#define INTERFACEUOMATERIAL_H

#include "Material.h"
#include "DispJumpUO_QP.h"

class czmUOMat;

template <>
InputParameters validParams<czmUOMat>();

/**
 *
 */
class czmUOMat : public Material
{
public:
  czmUOMat(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  // cohesive law paramters
  const std::vector<Real> _deltaU0;
  const std::vector<Real> _maxAllowableTraction;

  /// User objects computing the displacement jump
  const DispJumpUO_QP & _displacement_jump_UO;

  MaterialProperty<std::vector<Real>> & _displacement_jump;
  MaterialProperty<std::vector<Real>> & _displacement_jump_local;
  MaterialProperty<std::vector<Real>> & _traction;
  MaterialProperty<std::vector<Real>> & _traction_local;
  MaterialProperty<std::vector<std::vector<Real>>> & _traction_spatial_derivatives;
  MaterialProperty<std::vector<std::vector<Real>>> & _traction_spatial_derivatives_local;
  MaterialProperty<std::vector<Real>> & _residual;
  MaterialProperty<std::vector<std::vector<Real>>> & _jacobian;

  std::vector<Real> rotateVector(const std::vector<Real> /*v*/,
                                 const RealTensorValue /*R*/,
                                 const bool inverse = false);

  std::vector<std::vector<Real>> rotateTensor2(const std::vector<std::vector<Real>> /*T*/,
                                               const RealTensorValue /*R*/,
                                               const bool inverse = false);

  std::vector<Real> computeTractionLocal();
  std::vector<std::vector<Real>> computeTractionSpatialDerivativeLocal();
};

#endif // INTERFACEUOMATERIAL_H
