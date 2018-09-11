[Mesh]
  file = poly_3d.e
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
  [./right_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = 0.001*t
  [../]
  [./rear_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = rear
    function = -0.001*t
  [../]
  [./front_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = front
    function = 0.001*t
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
    execute_on = 'initial LINEAR timestep_end'
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
    type = czmUOMat
    is_interface_material = true
    boundary = 'interface'
    displacement_jump_UO = 'displacement_jump_uo'
    DeltaU0 = '1.0 0.5 0.5'
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
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu        superlu_dist'
  # petsc_options_value = 'hypre     boomerang'
  solve_type = newton
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-6
  nl_max_its = 5
  l_tol = 1e-10
  l_max_its = 50
  start_time = 0.0
  dt = 10.0
  end_time = 1000
  # line_search = none
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
