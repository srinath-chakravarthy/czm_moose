[Mesh]
  file = poly.msh
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
                1  0  0
                0 0 -1
                0 0  1'
     fixed_normal = true
     new_boundary = 'bottom top left right rear front'
     depends_on = breakmesh
   [../]

[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
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
  [./disp_z]
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
  [./szz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./syz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./sxz]
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
  [./szz]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
    variable = szz
    # block = '1 2 3'
  []
  [./syz]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 2
    variable = syz
    # block = '1 2 3'
  []
  [./sxz]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 2
    variable = sxz
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
    displacements = 'disp_x disp_y disp_z'
  [../]
[]
[BCs]
  # [./bottom_x]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = bottom
  #   value = 0.0
  # [../]
  # [./bottom_y]
  #   type = FunctionDirichletBC
  #   variable = disp_y
  #   boundary = bottom
  #   function = loadUnloadFunction_NEG
  # [../]
  # [./bottom_z]
  #   type = DirichletBC
  #   variable = disp_z
  #   boundary = bottom
  #   value = 0.0
  # [../]
  # [./top_x]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = top
  #   value = 0.0
  # [../]
  # [./top_y]
  #   type = FunctionDirichletBC
  #   variable = disp_y
  #   boundary = top
  #   function = loadUnloadFunction
  # [../]
  # [./top_z]
  #   type = DirichletBC
  #   variable = disp_z
  #   boundary = top
  #   value = 0.0
  # [../]
  [./left_x]
    type = DirichletBC
    variable = disp_y
    boundary = left
    value = 0.0
  [../]
  [./left_y]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = left
    function = loadUnloadFunction_NEG
  [../]
  [./left_z]
    type = DirichletBC
    variable = disp_z
    boundary = left
    value = 0.0
  [../]
  [./right_x]
    type = DirichletBC
    variable = disp_y
    boundary = right
    value = 0.0
  [../]
  [./right_y]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = loadUnloadFunction
  [../]
  [./right_z]
    type = DirichletBC
    variable = disp_z
    boundary = right
    value = 0.0
  [../]
  # [./left_x]
  #   type = FunctionDirichletBC
  #   variable = disp_x
  #   boundary = left
  #   function = loadUnloadFunction_NEG
  # [../]
  # [./right_x]
  #   type = FunctionDirichletBC
  #   variable = disp_x
  #   boundary = right
  #   function = loadUnloadFunction
  # [../]
  # [./rear_x]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = rear
  #   value = 0.0
  # [../]
  # [./rear_y]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = rear
  #   value = 0.0
  # [../]
  # [./rear_z]
  #   type = FunctionDirichletBC
  #   variable = disp_z
  #   boundary = rear
  #   function = loadUnloadFunction_NEG
  # [../]
  # [./front_x]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = front
  #   value = 0.0
  # [../]
  # [./front_y]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = front
  #   value = 0.0
  # [../]
  # [./front_z]
  #   type = FunctionDirichletBC
  #   variable = disp_z
  #   boundary = front
  #   function = loadUnloadFunction
  # [../]
[]
[Functions]
  [./loadUnloadFunction]
    type = PiecewiseLinear
    x = '0 4    10 16   20  30 40   65  70 120'
    y = '0 0.8  0  1.2  0   2   0   5   0   10'
  [../]
  [./loadUnloadFunction_NEG]
    type = PiecewiseLinear
    x = '0 4    10 16   20  30 40   65  70 120'
    y = '0 -0.8  0  -1.2  0   -2   0   -5   0   -10'
  [../]
[]
[InterfaceKernels]
  [./interface_x]
    type = czmInterfaceKernel
    variable = disp_x
    neighbor_var = disp_x
    disp_1 = disp_y
    disp_1_neighbor = disp_y
    disp_2 = disp_z
    disp_2_neighbor = disp_z
    disp_index = 0
    boundary = 'interface'
  [../]
  [./interface_y]
    type = czmInterfaceKernel
    variable = disp_y
    neighbor_var = disp_y
    disp_1 = disp_x
    disp_1_neighbor = disp_x
    disp_2 = disp_z
    disp_2_neighbor = disp_z
    disp_index = 1
    boundary = 'interface'
  [../]
  [./interface_z]
    type = czmInterfaceKernel
    variable = disp_z
    neighbor_var = disp_z
    disp_1 = disp_x
    disp_1_neighbor = disp_x
    disp_2 = disp_y
    disp_2_neighbor = disp_y
    disp_index = 2
    boundary = 'interface'
  [../]
[]
[UserObjects]
  [./displacement_jump_uo]
    type = DispJumpUO_QP
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    boundary = 'interface'
    execute_on = 'initial LINEAR NONLINEAR timestep_end'
  [../]
[]



[Materials]
  [./Elasticity_tensor]
    type = ComputeElasticityTensor
    # block = '1 2 3'
    fill_method = symmetric_isotropic
    C_ijkl = '0.3 0.5e8'
  [../]
  [./strain]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y disp_z'
    # block = '1 2 3'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    # block = '1 2 3'
  [../]
  [./gap]
    type = czmUOMatStateful
    is_interface_material = true
    boundary = 'interface'
    displacement_jump_UO = 'displacement_jump_uo'
    DeltaU0 = '2'
    MaxAllowableTraction = '100'
    Beta = 0.5
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
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu        superlu_dist'
  # petsc_options_value = 'hypre     boomerang'
  solve_type = PJFNK
  nl_abs_tol = 1e-6
  nl_rel_tol = 1e-6
  nl_max_its = 5
  l_tol = 1e-10
  l_max_its = 50
  start_time = 0.0
  dt = 1
  # dtmin = 0.1
  end_time = 120
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
    # block = 3
  [../]
  [./syy_3G]
    type = ElementAverageValue
    variable = syy
    execute_on = 'initial timestep_end'
    # block = 3
  [../]
  [./szz_3G]
    type = ElementAverageValue
    variable = szz
    execute_on = 'initial timestep_end'
    # block = 3
  [../]
  [./syz_3G]
    type = ElementAverageValue
    variable = syz
    execute_on = 'initial timestep_end'
    # block = 3
  [../]
  [./sxz_3G]
    type = ElementAverageValue
    variable = sxz
    execute_on = 'initial timestep_end'
    # block = 3
  [../]
  [./sxy_3G]
    type = ElementAverageValue
    variable = sxy
    execute_on = 'initial timestep_end'
    # block = 3
  [../]
[]
