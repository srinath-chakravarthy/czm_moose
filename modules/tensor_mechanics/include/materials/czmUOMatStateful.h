//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef CZMUOMATSTATEFUL_H
#define CZMUOMATSTATEFUL_H

#include "Material.h"
#include "DispJumpUO_QP.h"

class czmUOMatStateful;

template <>
InputParameters validParams<czmUOMatStateful>();

/**
 *
 */
class czmUOMatStateful : public Material
{
public:
  czmUOMatStateful(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  virtual void initQpStatefulProperties() override;

  // cohesive law paramters
  const Real _deltaU0;
  const Real _maxAllowableTraction;
  const Real _beta;

  /// User objects computing the displacement jump
  const DispJumpUO_QP & _displacement_jump_UO;

  MaterialProperty<std::vector<Real>> & _displacement_jump;
  MaterialProperty<std::vector<Real>> & _displacement_jump_local;
  MaterialProperty<Real> & _max_effective_jump;
  const MaterialProperty<Real> & _max_effective_jump_old;
  MaterialProperty<Real> & _effective_jump;
  const MaterialProperty<Real> & _effective_jump_old;
  MaterialProperty<std::vector<Real>> & _traction;
  MaterialProperty<std::vector<Real>> & _traction_local;
  MaterialProperty<Real> & _effective_traction;
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
  Real computeEffectiveJump();
  Real computeEffectiveTraction();
  Real computeEffectiveTractionNonLinear(Real /*effective_jump*/);
  Real computeEffectiveTractionLinear();
  bool checkLoadUnload();
  std::vector<std::vector<Real>> computeTractionSpatialDerivativeLocalNonLinear();
  std::vector<std::vector<Real>> computeTractionSpatialDerivativeLocalLinear();
};

#endif // CZMUOMATSTATEFUL_H
