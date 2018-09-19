//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef CZMTRACTIONSEPARATIONUOBASE_H
#define CZMTRACTIONSEPARATIONUOBASE_H

#include "SideUserObject.h"
#include "RankTwoTensor.h"

class CZMTractionSeparationUOBase;

template <>
InputParameters validParams<CZMTractionSeparationUOBase>();

/**
Traction sepration law base user object
 */
class CZMTractionSeparationUOBase : public SideUserObject
{
public:
  CZMTractionSeparationUOBase(const InputParameters & parameters);

  /// @{ Block all methods that are not used in explicitly called UOs
  virtual void initialize() override;
  virtual void execute() override final;
  virtual void finalize() override final;
  virtual void threadJoin(const UserObject &) override final;

  /// return the number of stateful material properties
  unsigned int getNumberStatefulMaterialProperties() const;

  /// return the history variable name
  std::string getStatefulMaterialPropertyName(const unsigned int /*mp_index*/) const;

  /// return the size of a history variable
  unsigned int getStatefulMaterialPropertySize(const unsigned int /*mp_index*/) const;

  /// return the intial values of the history variables
  std::vector<Real> getStatefulMaterialPropertysIntialValues(const unsigned int /*mp_index*/) const;

  /// retrun the initial values of a given material property
  virtual std::vector<Real> getNewStatefulMaterialProperty(const unsigned int /*qp*/,
                                                           const unsigned int /*mp_index*/) const;

  /// method returning if we are loading or unloading the material.
  /// this must bd overdden as different cohesive laws check load unload, differently
  virtual bool checkLoadUnload(const unsigned int /*qp*/) const;

  ///method computing the effective jump according to the give traction sepration law
  virtual Real getEffectiveJump(const unsigned int /*qp*/) const;

  /// method returning the traction value in local coordinates
  virtual std::vector<Real> computeTractionLocal(const unsigned int /*qp*/) const;

  /// method returning the traction derivates in local coordinates
  virtual std::vector<std::vector<Real>>
  computeTractionSpatialDerivativeLocal(const unsigned int /*qp*/) const;

protected:
  /// number of history variables present in the model
  const unsigned int _n_stateful_mp;
  const std::vector<std::string> _stateful_mp_names;
  const std::vector<unsigned int> _stateful_mp_sizes;
  const std::vector<std::vector<Real>> _stateful_mp_initial_values;

  const std::string _displacement_jump_mp_name;
  const MaterialProperty<std::vector<Real>> & _displacement_jump;
  const MaterialProperty<std::vector<Real>> & _displacement_jump_old;

  std::vector<std::vector<Real>> ResizeInitialValues() const;
};

#endif // CZMTRACTIONSEPARATIONUOBASE_H
