//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CZMNormals.h"

registerMooseObject("TensorMechanicsApp", CZMNormals);

template <>
InputParameters
validParams<CZMNormals>()
{
  InputParameters params = validParams<Material>();

  params.addClassDescription("this material class compute normals at each QP");
  return params;
}

CZMNormals::CZMNormals(const InputParameters & parameters)
  : Material(parameters), _normals_MP(declareProperty<RealVectorValue>("normals_MP"))
{
}

void
CZMNormals::computeQpProperties()
{
  _normals_MP[_qp] = _normals[_qp];
}
