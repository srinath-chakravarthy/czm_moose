//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CZMTractionSeparationUOBase.h"
#include "Material.h"
#include "MooseError.h"

template <>
InputParameters
validParams<CZMTractionSeparationUOBase>()
{
  InputParameters params = validParams<SideUserObject>();
  params.addClassDescription(
      "User Object implementing basic functions for traction separations law");
  params.set<ExecFlagEnum>("execute_on") = EXEC_CUSTOM;
  params.suppressParameter<ExecFlagEnum>("execute_on");
  params.addRequiredParam<unsigned int>("n_stateful_mp", "number of stateful material properties");
  params.addParam<std::vector<std::string>>(
      "stateful_mp_names", std::vector<std::string>(0), "name of stateful material properties");
  params.addParam<std::vector<unsigned int>>("stateful_mp_sizes",
                                             std::vector<unsigned int>(0),
                                             "size of each stateful material properties");
  params.addParam<std::vector<Real>>(
      "stateful_mp_initial_values", std::vector<Real>(0), "intial material proeprties values");
  params.addRequiredParam<unsigned int>("n_non_stateful_mp",
                                        "number of NON-stateful material properties");
  params.addParam<std::vector<std::string>>("non_stateful_mp_names",
                                            std::vector<std::string>(0),
                                            "name of NON stateful material properties");
  params.addParam<std::vector<unsigned int>>("non_stateful_mp_sizes",
                                             std::vector<unsigned int>(0),
                                             "size of each NON stateful material properties");
  params.addRequiredParam<std::string>("displacement_jump_mp_name",
                                       "the name of the material property hosting the displacement "
                                       "jump computed in natural element coordinate system");

  return params;
}

CZMTractionSeparationUOBase::CZMTractionSeparationUOBase(const InputParameters & parameters)
  : SideUserObject(parameters),
    _n_stateful_mp(getParam<unsigned int>("n_stateful_mp")),
    _stateful_mp_names(getParam<std::vector<std::string>>("stateful_mp_names")),
    _stateful_mp_sizes(getParam<std::vector<unsigned int>>("stateful_mp_sizes")),
    _stateful_mp_initial_values(ResizeInitialValues()),
    _n_non_stateful_mp(getParam<unsigned int>("n_non_stateful_mp")),
    _non_stateful_mp_names(getParam<std::vector<std::string>>("non_stateful_mp_names")),
    _non_stateful_mp_sizes(getParam<std::vector<unsigned int>>("non_stateful_mp_sizes")),
    _displacement_jump(getMaterialPropertyByName<std::vector<Real>>(
        getParam<std::string>("displacement_jump_mp_name"))),
    _displacement_jump_old(getMaterialPropertyOldByName<std::vector<Real>>(
        getParam<std::string>("displacement_jump_mp_name")))

{
  if (_n_stateful_mp != _stateful_mp_names.size())
  {
    std::cout << "_n_stateful_mp: " << _n_stateful_mp
              << " stateful_mp_names: " << _stateful_mp_names.size() << std::endl;

    for (unsigned int i = 0; i < _stateful_mp_names.size(); i++)
      std::cout << _stateful_mp_names[i] << std::endl;
    std::cout << std::endl;
    mooseError("CZMTractionSeparationUOBase:: n_stateful_mp does match with the number of supplied "
               "material properies names  stateful_mp_names");
  }
  if (_n_stateful_mp != _stateful_mp_sizes.size())

    mooseError("CZMTractionSeparationUOBase:: n_stateful_mp does match with the number of supplied "
               "material properies sizes stateful_mp_sizes");

  unsigned int total_size_stateful_mp_initial_values = 0;
  for (unsigned int i = 0; i < _stateful_mp_initial_values.size(); i++)
    total_size_stateful_mp_initial_values += _stateful_mp_initial_values[i].size();

  if (std::accumulate(_stateful_mp_sizes.begin(), _stateful_mp_sizes.end(), 0) !=
      (int)total_size_stateful_mp_initial_values)
  {
    std::cout << "_stateful_mp_sizes: "
              << std::accumulate(_stateful_mp_sizes.begin(), _stateful_mp_sizes.end(), 0)
              << " _stateful_mp_initial_values: " << (int)_stateful_mp_initial_values.size()
              << std::endl;

    mooseError("CZMTractionSeparationUOBase:: the number of supplied initial values does not match "
               "the material properies size");
  }
}

unsigned int
CZMTractionSeparationUOBase::getNumberStatefulMaterialProperties() const
{
  return _n_stateful_mp;
}

std::string
CZMTractionSeparationUOBase::getStatefulMaterialPropertyName(unsigned int mp_index) const
{
  return _stateful_mp_names[mp_index];
}

unsigned int
CZMTractionSeparationUOBase::getStatefulMaterialPropertySize(unsigned int mp_index) const
{
  return _stateful_mp_sizes[mp_index];
}

std::vector<Real>
CZMTractionSeparationUOBase::getStatefulMaterialPropertysIntialValues(unsigned int mp_index) const
{
  return _stateful_mp_initial_values[mp_index];
}

std::vector<Real>
CZMTractionSeparationUOBase::getNewStatefulMaterialProperty(unsigned int /*qp*/,
                                                            unsigned int /*mp_index*/) const
{
  mooseError("CZMTractionSeparationUOBase::getNewStatefulMaterialProperty should never "
             "be called directly but always subclassed");
}

unsigned int
CZMTractionSeparationUOBase::getNumberNonStatefulMaterialProperties() const
{
  return _n_non_stateful_mp;
}

std::string
CZMTractionSeparationUOBase::getNonStatefulMaterialPropertyName(unsigned int mp_index) const
{
  return _non_stateful_mp_names[mp_index];
}

unsigned int
CZMTractionSeparationUOBase::getNonStatefulMaterialPropertySize(unsigned int mp_index) const
{
  return _non_stateful_mp_sizes[mp_index];
}

std::vector<Real>
CZMTractionSeparationUOBase::getNewNonStatefulMaterialProperty(unsigned int /*qp*/,
                                                               unsigned int /*mp_index*/) const
{
  mooseError("CZMTractionSeparationUOBase::getNewStatefulMaterialProperty should never "
             "be called directly but always subclassed");
}

bool
CZMTractionSeparationUOBase::checkLoadUnload(const unsigned int /*qp*/) const
{
  mooseError("CZMTractionSeparationUOBase::checkLoadUnload should never "
             "be called directly but always subclassed");
}

std::vector<Real>
CZMTractionSeparationUOBase::computeTractionLocal(unsigned int /*qp*/) const
{
  mooseError("CZMTractionSeparationUOBase::computeTractionLocal should never "
             "be called directly but always subclassed");
}

std::vector<std::vector<Real>>
CZMTractionSeparationUOBase::computeTractionSpatialDerivativeLocal(unsigned int /*qp*/) const
{
  mooseError("CZMTractionSeparationUOBase::computeTractionLocal should never "
             "be called directly but always subclassed");
}

std::vector<std::vector<Real>>
CZMTractionSeparationUOBase::ResizeInitialValues() const
{
  std::vector<std::vector<Real>> temp;
  std::vector<Real> temp_init_values = getParam<std::vector<Real>>("stateful_mp_initial_values");
  unsigned int n_subvector = _stateful_mp_sizes.size();
  temp.resize(n_subvector);
  unsigned int c = 0;
  for (unsigned int i = 0; i < n_subvector; i++)
  {
    unsigned int subvector_size = _stateful_mp_sizes[i];
    temp[i].resize(subvector_size);
    for (unsigned int j = 0; j < subvector_size; j++)
    {
      temp[i][j] = temp_init_values[c];
      c += 1;
    }
  }
  return temp;
}

Real
CZMTractionSeparationUOBase::getEffectiveJump(const unsigned int /*qp*/) const
{
  mooseError("getEffectiveJump CZMTractionSeparationUOBase must be overridden");
}

// ovverride standard UO functions
void
CZMTractionSeparationUOBase::initialize()
{
}

void
CZMTractionSeparationUOBase::execute()
{
  mooseError("execute CZMTractionSeparationUOBase must be called explicitly from Materials");
}

void
CZMTractionSeparationUOBase::finalize()
{
  mooseError("finalize CZMTractionSeparationUOBase must be called explicitly from Materials");
}

void
CZMTractionSeparationUOBase::threadJoin(const UserObject &)
{
  mooseError("threadJoin CZMTractionSeparationUOBase must be called explicitly from Materials");
}
