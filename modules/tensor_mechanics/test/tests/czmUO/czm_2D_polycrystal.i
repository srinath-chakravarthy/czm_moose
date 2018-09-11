[Mesh]
  file = poly2d.msh
  parallel_type = REPLICATED
[]

[MeshModifiers]
  [./breakmesh]
    type = BreakMeshByBlock
    # split_interface = false
  [../]

  [./add_side_sets]
     type = SideSetsFromNormals
     normals = '0  -1  0
                0  1  0
                -1 0  0
                1  0  0'
     fixed_normal = true
     new_boundary = 'bottom top left right'
     depends_on = breakmesh
   [../]

  # [./bottom]
  #   type = SideSetsAroundSubdomain
  #   depends_on = 'breakmesh'
  #   new_boundary = 'bottom'
  #   normal = '0 -1 0'
  # [../]
  # [./top]
  #   type = SideSetsAroundSubdomain
  #   depends_on = 'breakmesh'
  #   new_boundary = 'top'
  #   normal = '0 1 0'
  # [../]
  # [./left]
  #   type = SideSetsAroundSubdomain
  #   depends_on = 'breakmesh'
  #   new_boundary = 'left'
  #   normal = '-1 0 0'
  # [../]
  # [./right]
  #   type = SideSetsAroundSubdomain
  #   depends_on = 'breakmesh'
  #   new_boundary = 'right'
  #   normal = '1 0 0'
  # [../]

[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Modules/TensorMechanics/Master]
  strain = SMALL
  add_variables = true
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]

[]
[AuxVariables]
  [./sxx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./syy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./sxy]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]
[AuxKernels]
  [./sxx]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    variable = sxx
    # block = '1 2 3'
  []
  [./syy]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    variable = syy
    # block = '1 2 3'
  []
  [./sxy]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    variable = sxy
    # block = '1 2 3'
  []
[]
[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
  [../]
[]
[BCs]
  [./bottom_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = bottom
    function = -0.001*t
  [../]
  [./top_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = 0.001*t
  [../]
  [./left_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = left
    function = -0.001*t
  [../]
  [./right]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = 0.001*t
  [../]

  # [./bottom_z]
  #   type = DirichletBC
  #   variable = disp_z
  #   boundary = bottom_1
  #   value = 0.0
  # [../]


[]
[InterfaceKernels]
  [./interface_x]
    type = czmInterfaceKernel
    variable = disp_x
    neighbor_var = disp_x
    disp_1 = disp_y
    disp_1_neighbor = disp_y
    disp_index = 0
    boundary = 'interface'
  [../]
  [./interface_y]
    type = czmInterfaceKernel
    variable = disp_y
    neighbor_var = disp_y
    disp_1 = disp_x
    disp_1_neighbor = disp_x
    disp_index = 1
    boundary = 'interface'
  [../]
[]
[UserObjects]
  [./displacement_jump_uo]
    type = DispJumpUO_QP
    disp_x = disp_x
    disp_y = disp_y
    boundary = 'interface'
    execute_on = 'initial LINEAR timestep_end'
  [../]
[]

[Materials]
  [./Elasticity_tensor]
    type = ComputeElasticityTensor
    # block = 'ANY'
    fill_method = symmetric_isotropic
    C_ijkl = '0.3 0.5e8'
  [../]
  [./strain]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    # block = 'ANY'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    # block = 'ANY'
  [../]
  [./gap]
    type = czmUOMat
    is_interface_material = true
    boundary = 'interface'
    displacement_jump_UO = 'displacement_jump_uo'
    DeltaU0 = '0.01 0.005 0.005'
    MaxAllowableTraction = '100 70 70'
  [../]
[]
 [Preconditioning]
   [./SMP]
     type = SMP
     full = true
   [../]
 []
[Executioner]
  # Preconditisoned JFNK (default)
  type = Transient
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  # petsc_options_value = 'hypre     boomerang'
  solve_type = NEWTON
  nl_abs_tol = 1e-6
  nl_rel_tol = 1e-6
  nl_max_its = 5
  l_tol = 1e-10
  l_max_its = 50
  start_time = 0.0
  dt = 50.0
  dtmin = 1e-5
  end_time = 1000
  line_search = none
[]
[Outputs]
  [./out]
    type = Exodus
  [../]
[]
[Postprocessors]
  [./sxx_3G]
    type = ElementAverageValue
    variable = sxx
    execute_on = 'initial timestep_end'
  [../]
  [./syy_3G]
    type = ElementAverageValue
    variable = syy
    execute_on = 'initial timestep_end'
  [../]
  [./sxy_3G]
    type = ElementAverageValue
    variable = sxy
    execute_on = 'initial timestep_end'
  [../]
[]
