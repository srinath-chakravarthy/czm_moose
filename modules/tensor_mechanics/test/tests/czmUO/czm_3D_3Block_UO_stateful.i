[Mesh]
  file = coh3D_3Blocks.e
  parallel_type = REPLICATED
[]

[MeshModifiers]
  [./breakmesh]
    type = BreakMeshByBlock
    # split_interface = false
  [../]

  [./bottom_block_1]
    type = SideSetsAroundSubdomain
    depends_on = 'breakmesh'
    block = '1'
    new_boundary = 'bottom_1'
    normal = '0 0 -1'
  [../]
  [./top_block_2]
    type = SideSetsAroundSubdomain
    depends_on = 'breakmesh'
    block = '2'
    new_boundary = 'top_2'
    normal = '0 0 1'
  [../]
  [./top_block_3]
    type = SideSetsAroundSubdomain
    depends_on = 'breakmesh'
    block = '3'
    new_boundary = 'top_3'
    normal = '0 0 1'
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
    block = '1 2 3'
  []
  [./syy]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    variable = syy
    block = '1 2 3'
  []
  [./szz]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
    variable = szz
    block = '1 2 3'
  []
  [./syz]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 2
    variable = syz
    block = '1 2 3'
  []
  [./sxz]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 2
    variable = sxz
    block = '1 2 3'
  []
  [./sxy]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    variable = sxy
    block = '1 2 3'
  []
[]
[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y disp_z'
  [../]
[]
[BCs]
  [./bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = bottom_1
    value = 0.0
  [../]
  [./bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom_1
    value = 0.0
  [../]
  [./bottom_z]
    type = DirichletBC
    variable = disp_z
    boundary = bottom_1
    value = 0.0
  [../]
  [./top2_x]
    type = DirichletBC
    variable = disp_x
    boundary = top_2
    value = 0.0
  [../]
  [./top2_y]
    type = DirichletBC
    variable = disp_y
    boundary = top_2
    value = 0.0
  [../]
  [./top2_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = top_2
    function = 0.02*t
  [../]
  [./top3_x]
    type = DirichletBC
    variable = disp_x
    boundary = top_3
    value = 0.0
  [../]
  [./top3_y]
    type = DirichletBC
    variable = disp_y
    boundary = top_3
    value = 0.0
  [../]
  [./top3_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = top_3
    function = 0.02*t
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
    execute_on = 'initial NONLINEAR LINEAR timestep_end'
  [../]
  [./cohesive_law_exponential]
    type = CohesiveLaw_Exponential
    stateful_mp_names = 'max_effective_jump max_effective_traction'
    stateful_mp_sizes ='1 1'
    stateful_mp_initial_values ='0 0'
    displacement_jump_peak = 1
    traction_peak = 150
    displacement_jump_mp_name = 'displacement_jump_local'
    boundary = 'interface'
  [../]
[]

[Materials]
  [./Elasticity_tensor]
    type = ComputeElasticityTensor
    block = '1 2 3'
    fill_method = symmetric_isotropic
    C_ijkl = '0.3 0.5e8'
  [../]
  [./strain]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y disp_z'
    block = '1 2 3'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = '1 2 3'
  [../]
  [./gap]
    type = CZMUOBasedMaterial
    is_interface_material = true
    boundary = 'interface'
    displacement_jump_UO = 'displacement_jump_uo'
    traction_separation_UO = 'cohesive_law_exponential'
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
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-6
  nl_max_its = 5
  l_tol = 1e-10
  l_max_its = 50
  start_time = 0.0
  dt = 10
  end_time = 500
  dtmin = 1
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
    block = 3
  [../]
  [./syy_3G]
    type = ElementAverageValue
    variable = syy
    execute_on = 'initial timestep_end'
    block = 3
  [../]
  [./szz_3G]
    type = ElementAverageValue
    variable = szz
    execute_on = 'initial timestep_end'
    block = 3
  [../]
  [./syz_3G]
    type = ElementAverageValue
    variable = syz
    execute_on = 'initial timestep_end'
    block = 3
  [../]
  [./sxz_3G]
    type = ElementAverageValue
    variable = sxz
    execute_on = 'initial timestep_end'
    block = 3
  [../]
  [./sxy_3G]
    type = ElementAverageValue
    variable = sxy
    execute_on = 'initial timestep_end'
    block = 3
  [../]
[]
