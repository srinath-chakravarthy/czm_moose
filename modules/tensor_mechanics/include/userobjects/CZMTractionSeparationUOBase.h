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

#include "DiscreteElementUserObject.h"
#include "RankTwoTensor.h"

class CZMTractionSeparationUOBase;

template <>
InputParameters validParams<CZMTractionSeparationUOBase>();

/**
Traction sepration law basic user object
 */
class CZMTractionSeparationUOBase : public DiscreteElementUserObject
{
public:
  CZMTractionSeparationUOBase(const InputParameters & parameters);

  /// @{ Block all methods that are not used in explicitly called UOs
  // virtual void initialize() override;
  // virtual void execute() override final;
  // virtual void finalize() override final;
  // virtual void threadJoin(const UserObject &) override final;

  /// return the number of stateful material properties
  unsigned int getNumberHistoryVariables() const;

  /// return the intiali values of the history variables
  std::vector<Real> getHistoryVariablesIntialValues() const;

  /// method returning if we are loading or unloading the material.
  /// this must bd overdden as different cohesive laws check load unload, differently
  virtual bool checkLoadUnload(const std::vector<Real> /*current_jump*/,
                               const std::vector<Real> /*old_jump*/) const;

  virtual std::vector<Real> getNewStatefulMaterialProperty() const;

  // /// initialization of stateful material properties
  // virtual void initStatefulMaterialProperty(unsigned int /*materialPropertyID*/,
  //                                           std::vector<Real> & /*statefulePropertyValue*/)
  //                                           const;

  /// method updating stateful material properties
  // virtual void
  // updateStatefulMaterialProperty(unsigned int /*qp*/,
  //                                unsigned int /*materialPropertyID*/,
  //                                std::vector<Real> & /*statefulePropertyValue*/,
  //                                const std::vector<Real> & /*statefulePropertyValue_old*/) const;

  /// method returning the traction value in local coordinates
  virtual std::vector<Real>
  computeTractionLocal(const std::vector<Real> /*current_stateful_MP*/,
                       const std::vector<Real> /*old_statedul_MP*/,
                       const std::vector<Real> /*current_displacement_jump*/,
                       const std::vector<Real> /*old_displacement_jump*/,
                       const std::vector<Real> other_current_required_mp = std::vector<Real>(),
                       const std::vector<Real> other_old_required_mp = std::vector<Real>());

  /// method returning the traction derivates in local coordinates
  virtual std::vector<std::vector<Real>> computeTractionSpatialDerivativeLocal(
      const std::vector<Real> /*current_stateful_MP*/,
      const std::vector<Real> /*old_statedul_MP*/,
      const std::vector<Real> /*current_displacement_jump*/,
      const std::vector<Real> /*old_displacement_jump*/,
      const std::vector<Real> other_current_required_mp = std::vector<Real>(),
      const std::vector<Real> other_old_required_mp = std::vector<Real>());

protected:
  /// number of history variables present in the model
  const unsigned int _n_history_variables;
  const std::vector<Real> _history_variables_initial_values;
  /// The dispalcement jump accross the interface
  // const MaterialProperty<RealVectorValue> & _JumpLocal;
  // // the variable containing the list of the stateful
  // // material properties variables
  // const std::vector<std::string> _cohesive_law_stateful_properties_names;
};

#endif // CZMTRACTIONSEPARATIONUOBASE_H
