//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef DISPJUMPANDNORMALSUO_QP_H
#define DISPJUMPANDNORMALSUO_QP_H

#include "InterfaceUserObject.h"

class DispJumpAndNormalsUO_QP;

template <>
InputParameters validParams<DispJumpAndNormalsUO_QP>();

/**
 *
 */
class DispJumpAndNormalsUO_QP : public InterfaceUserObject
{
public:
  DispJumpAndNormalsUO_QP(const InputParameters & parameters);
  virtual ~DispJumpAndNormalsUO_QP();

  virtual void initialize();
  virtual void execute();
  virtual void finalize() { return; };
  virtual void threadJoin(const UserObject & /*uo*/) { return; };

  RealVectorValue getDisplacementJump(dof_id_type elem, unsigned int side, unsigned int qp) const;

  RealVectorValue getNormalMaster(dof_id_type elem, unsigned int side, unsigned int qp) const;
  RealVectorValue getNormalSlave(dof_id_type elem, unsigned int side, unsigned int qp) const;

protected:
  /// this map is used for storing data at QPs.
  /// keys<element_id, side_id>, values<vector (1 elem_per_QP)<vector<Real (_mean_mat_prop, _var_jump)>>>
  std::map<std::pair<dof_id_type, unsigned int>, std::vector<std::vector<RealVectorValue>>>
      _map_values;
  const VariableValue & _ux;
  const VariableValue & _ux_neighbor;
  const VariableValue & _uy;
  const VariableValue & _uy_neighbor;
  const VariableValue & _uz;
  const VariableValue & _uz_neighbor;
  const MaterialProperty<RealVectorValue> & _normals_MP;
  const MaterialProperty<RealVectorValue> & _normals_MP_neighbor;
};

#endif // DISPJUMPANDNORMALSUO_QP_H
