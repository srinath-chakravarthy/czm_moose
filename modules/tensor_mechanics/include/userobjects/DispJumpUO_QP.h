//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef DISPJUMPUO_QP_H
#define DISPJUMPUO_QP_H

#include "InterfaceUserObject.h"

class DispJumpUO_QP;

template <>
InputParameters validParams<DispJumpUO_QP>();

/**
 *
 */
class DispJumpUO_QP : public InterfaceUserObject
{
public:
  DispJumpUO_QP(const InputParameters & parameters);
  virtual ~DispJumpUO_QP();

  virtual void initialize();
  virtual void execute();
  virtual void finalize() { return; };
  virtual void threadJoin(const UserObject & /*uo*/) { return; };

  std::vector<Real> getDisplacementJump(dof_id_type elem, unsigned int side, unsigned int qp) const;

protected:
  /// this map is used for storing data at QPs.
  /// keys<element_id, side_id>, values<vector (1 elem_per_QP)<vector<Real (_mean_mat_prop, _var_jump)>>>
  std::map<std::pair<dof_id_type, unsigned int>, std::vector<std::vector<Real>>> _map_values;
  const VariableValue & _ux;
  const VariableValue & _ux_neighbor;
  const VariableValue & _uy;
  const VariableValue & _uy_neighbor;
  const VariableValue & _uz;
  const VariableValue & _uz_neighbor;
};

#endif // DISPJUMPUO_QP_H
