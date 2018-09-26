//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef CZMNORMALS_H
#define CZMNORMALS_H

#include "Material.h"

// Forward declaration
class CZMNormals;

template <>
InputParameters validParams<CZMNormals>();

/**
 * CZMNormals is material returning normals as MP at each QP.
 */
class CZMNormals : public Material
{
public:
  CZMNormals(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  MaterialProperty<RealVectorValue> & _normals_MP;
  MaterialProperty<unsigned int> & __elem;
  MaterialProperty<unsigned int> & __qp;
};

#endif // CZMNORMALS_H
