//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DispJumpUO_QP.h"
#include "MooseMesh.h"
registerMooseObject("TensorMechanicsApp", DispJumpUO_QP);

template <>
InputParameters
validParams<DispJumpUO_QP>()
{
  InputParameters params = validParams<InterfaceUserObject>();
  params.addParam<MaterialPropertyName>("diffusivity",
                                        0.0,
                                        "The name of the diffusivity material property that "
                                        "will be used in the flux computation.");
  params.addParam<bool>(
      "use_old_prop",
      false,
      "A Boolean to indicate whether the current or old value of a material prop should be used.");
  params.addRequiredCoupledVar("disp_x", "displacement in X");
  params.addCoupledVar("disp_y", "displacement in Y");
  params.addCoupledVar("disp_z", "displacement in Z");
  return params;
}

DispJumpUO_QP::DispJumpUO_QP(const InputParameters & parameters)
  : InterfaceUserObject(parameters),
    _ux(coupledValue("disp_x")),
    _ux_neighbor(coupledNeighborValue("disp_x")),
    _uy(_mesh.dimension() >= 2 ? coupledValue("disp_y") : _zero),
    _uy_neighbor(_mesh.dimension() >= 2 ? coupledNeighborValue("disp_y") : _zero),
    _uz(_mesh.dimension() >= 3 ? coupledValue("disp_z") : _zero),
    _uz_neighbor(_mesh.dimension() >= 3 ? coupledNeighborValue("disp_z") : _zero)
{
}

DispJumpUO_QP::~DispJumpUO_QP() {}

void
DispJumpUO_QP::initialize()
{

  // define the boundary map nad retrieve element side and boundary_ID
  std::vector<std::tuple<dof_id_type, unsigned short int, boundary_id_type>> elem_side_bid =
      _mesh.buildSideList();

  // retrieve on which boudnary this UO operates
  std::set<BoundaryID> boundaryList = boundaryIDs();

  // clear map values
  _map_values.clear();

  // initialize the map_values looping over all the element and sides
  for (unsigned int i = 0; i < elem_side_bid.size(); i++)
  {
    // check if this boundary
    // if this element side is part of the boundary then add elements to the map
    if (boundaryList.find(std::get<2>(elem_side_bid[i])) != boundaryList.end())
    {

      // make pair
      std::pair<dof_id_type, unsigned int> elem_side_pair =
          std::make_pair(std::get<0>(elem_side_bid[i]), std::get<1>(elem_side_bid[i]));
      // initialize map elemenet
      std::vector<std::vector<Real>> var_values(0, std::vector<Real>(3, 0));

      // add entry to the value map
      _map_values[elem_side_pair] = var_values;
    }
  }
}

void
DispJumpUO_QP::execute()
{

  // find the entry on the map
  auto it = _map_values.find(std::make_pair(_current_elem->id(), _current_side));
  if (it != _map_values.end())
  {

    // insert two vector value for each qp
    auto & vec = _map_values[std::make_pair(_current_elem->id(), _current_side)];
    vec.resize(_qrule->n_points());

    // loop over qps and do stuff
    for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
    {
      vec[qp].resize(3, 0);

      // compute displacement jump
      vec[qp][0] = _ux_neighbor[qp] - _ux[qp];
      vec[qp][1] = _uy_neighbor[qp] - _uy[qp];
      vec[qp][2] = _uz_neighbor[qp] - _uz[qp];
    }
  }
  else
    mooseError("DispJumpUO_QP::execute cannot fine the required element and side");
}

std::vector<Real>
DispJumpUO_QP::getDisplacementJump(dof_id_type elem, unsigned int side, unsigned int qp) const
{
  auto dispJump = _map_values.find(std::make_pair(elem, side));
  if (dispJump != _map_values.end())
  {
    return dispJump->second[qp];
  }
  else
    mooseError("DispJumpUO_QP::getDisplacementJump can't find the given qp");
}
