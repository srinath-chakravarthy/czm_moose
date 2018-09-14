//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CZMTractionSeparationUOBase.h"
// #include "Material.h"
#include "MooseError.h"

template <>
InputParameters
validParams<CZMTractionSeparationUOBase>()
{
  InputParameters params = validParams<DiscreteElementUserObject>();
  params.addClassDescription(
      "User Object implementing basic functions for traction separations law");
  params.set<ExecFlagEnum>("execute_on") = EXEC_CUSTOM;
  params.suppressParameter<ExecFlagEnum>("execute_on");
  params.addParam<unsigned int>("n_history_variables", 0, "number of stateful cohesive properties");
  params.addParam<std::vector<Real>>("history_variables_initial_values",
                                     std::vector<Real>(0),
                                     "number of stateful cohesive properties");

  return params;
}

CZMTractionSeparationUOBase::CZMTractionSeparationUOBase(const InputParameters & parameters)
  : DiscreteElementUserObject(parameters),
    _n_history_variables(getParam<unsigned int>("n_history_variables")),
    _history_variables_initial_values(
        getParam<std::vector<Real>>("history_variables_initial_values"))
// _cohesive_law_stateful_properties_names(
//     getParam<std::vector<std::string>>("cohesive_law_stateful_properties_names"))
{
  if (_n_history_variables != _history_variables_initial_values.size())
    mooseError("CZMTractionSeparationUOBase:: length of history_variables_initial_values does not "
               "match n_history_variables");
}

// void
// CZMTractionSeparationUOBase::statefulMaterialPropertyNames(
//     std::vector<std::string> & materialPropertyNames) const
// {
//   unsigned int num_prop = _cohesive_law_stateful_properties_names.size();
//   materialPropertyNames.resize(num_prop);
//
//   for (unsigned int i = 0; i < num_prop; i++)
//     materialPropertyNames[i] = _cohesive_law_stateful_properties_names[i];
// }

unsigned int
CZMTractionSeparationUOBase::getNumberHistoryVariables() const
{
  return _n_history_variables;
}

std::vector<Real>
CZMTractionSeparationUOBase::getHistoryVariablesIntialValues() const
{
  if (_n_history_variables == 0)
    mooseError("CZMTractionSeparationUOBase:: there are no history variables to initialize. This "
               "method should not be called");
  return _history_variables_initial_values;
}

bool
CZMTractionSeparationUOBase::checkLoadUnload(const std::vector<Real> /*current_jump*/,
                                             const std::vector<Real> /*old_jump*/) const
{
  mooseError("CZMTractionSeparationUOBase::checkLoadUnload should never "
             "be called directly but always subclassed");
}

std::vector<Real>
CZMTractionSeparationUOBase::getNewStatefulMaterialProperty() const
{
  mooseError("CZMTractionSeparationUOBase::getNewStatefulMaterialProperty should never "
             "be called directly but always subclassed");
}

std::vector<Real>
CZMTractionSeparationUOBase::computeTractionLocal(
    const std::vector<Real> /*current_stateful_MP*/,
    const std::vector<Real> /*old_statedul_MP*/,
    const std::vector<Real> /*current_displacement_jump*/,
    const std::vector<Real> /*old_displacement_jump*/,
    const std::vector<Real> /*other_current_required_mp = std::vector<Real>()*/,
    const std::vector<Real> /*other_old_required_mp = std::vector<Real>()*/)
{
  mooseError("CZMTractionSeparationUOBase::computeTractionLocal should never "
             "be called directly but always subclassed");
}

std::vector<std::vector<Real>>
CZMTractionSeparationUOBase::computeTractionSpatialDerivativeLocal(
    const std::vector<Real> /*current_stateful_MP*/,
    const std::vector<Real> /*old_statedul_MP*/,
    const std::vector<Real> /*current_displacement_jump*/,
    const std::vector<Real> /*old_displacement_jump*/,
    const std::vector<Real> /*other_current_required_mp = std::vector<Real>()*/,
    const std::vector<Real> /*other_old_required_mp = std::vector<Real>()*/)
{
  mooseError("CZMTractionSeparationUOBase::computeTractionLocal should never "
             "be called directly but always subclassed");
}

// ovverride standard UO functions
// void
// CZMTractionSeparationUOBase::initialize()
// {
// }
//
// void
// CZMTractionSeparationUOBase::execute()
// {
//   mooseError("execute CZMTractionSeparationUOBase must be called explicitly from Materials");
// }
//
// void
// CZMTractionSeparationUOBase::finalize()
// {
//   mooseError("finalize CZMTractionSeparationUOBase must be called explicitly from
//   Materials");
// }
//
// void
// CZMTractionSeparationUOBase::threadJoin(const UserObject &)
// {
//   mooseError("threadJoin CZMTractionSeparationUOBase must be called explicitly from
//   Materials");
// }
