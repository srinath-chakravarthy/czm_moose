//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef CZMUOBASEDMATERIAL_H
#define CZMUOBASEDMATERIAL_H

#include "Material.h"
#include "DispJumpUO_QP.h"
#include "CZMTractionSeparationUOBase.h"
class CZMUOBasedMaterial;
template <>
InputParameters validParams<CZMUOBasedMaterial>();
/**
 *
 */
class CZMUOBasedMaterial : public Material
{
public:
  CZMUOBasedMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  virtual void initQpStatefulProperties() override;

  /// User objects computing the displacement jump
  const DispJumpUO_QP & _displacement_jump_UO;

  /// User objectets defining the traction separation law
  /// non linear TS law
  const CZMTractionSeparationUOBase & _traction_separation_UO;
  /// unloading behavior for stateful laws
  const CZMTractionSeparationUOBase & _unload_traction_separation_UO;
  /// penalty for copentration behavior
  const CZMTractionSeparationUOBase & _coopenetration_penalty_UO;

  /// the disaplcement jump in global coordiantes
  MaterialProperty<std::vector<Real>> & _displacement_jump;

  /// the disaplcement jump in natural element coordiantes
  MaterialProperty<std::vector<Real>> & _displacement_jump_local;
  /// the disaplcement jump in natural element coordiantes at the previous time step
  const MaterialProperty<std::vector<Real>> & _displacement_jump_local_old;

  /// the value of the Traction in global coordiantes
  MaterialProperty<std::vector<Real>> & _traction;

  /// the value of the Traction in natural element coordiantes
  MaterialProperty<std::vector<Real>> & _traction_local;

  /// the value of the traction derivatives in global coordiantes
  MaterialProperty<std::vector<std::vector<Real>>> & _traction_spatial_derivatives;

  /// the value of the traction derivatives in natural element coordiantes
  MaterialProperty<std::vector<std::vector<Real>>> & _traction_spatial_derivatives_local;

  /// the material property in which the residual is stored
  MaterialProperty<std::vector<Real>> & _czm_residual;

  /// the material property in which the jacobian is stored
  MaterialProperty<std::vector<std::vector<Real>>> & _czm_jacobian;

  const unsigned int _n_uo_czm_properties;
  std::vector<MaterialProperty<std::vector<Real>> *> _uo_czm_properties;
  std::vector<const MaterialProperty<std::vector<Real>> *> _uo_czm_properties_old;

  // /// rotation matrix rotating a vector V to (0, 0, 1) i.e. _VLocal = R*_V
  // RealTensorValue _RotationGlobal2Local;

  /// Rotate a vector "T" via the rotation matrix "R".
  /// inverse rotation is achieved by setting "inverse" = true
  std::vector<Real> rotateVector(const std::vector<Real> /*V*/,
                                 const RealTensorValue /*R*/,
                                 const bool inverse = false);

  /// Rotate a rank2 tensor "T" via the rotation matrix "R".
  /// inverse rotation is achieved by setting "inverse" = true
  std::vector<std::vector<Real>> rotateTensor2(const std::vector<std::vector<Real>> /*T*/,
                                               const RealTensorValue /*R*/,
                                               const bool inverse = false);
};

#endif // CZMUOBASEDMATERIAL_H
