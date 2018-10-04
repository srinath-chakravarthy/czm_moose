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
  [./all]
    strain = SMALL
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz stress_yz stress_xz stress_xy'
    use_finite_deform_jacobian = FALSE
  [../]
[]



[Functions]
  [./loadUnloadFunction]
    type = PiecewiseLinear
    x = '0 4     8 14     21      32    42   67   92 142'
    y = '0 0.08  0  0.12  -0.02   0.2   0    0.5   0   1'
    # x = '0 1    2 '
    # y = '0 -0.2 0 '
  [../]
  [./loadUnloadPressure]
    type = PiecewiseLinear
    x = '0 4     8 14     21      32    42   67   92 142'
    y = '0 0.08  0  0.12  -0.02   0.2   0    0.5   0   1'
    # x = '0 1    2 '
    # y = '0 -0.2 0 '
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
  # [./top2_z]
  #   type = FunctionDirichletBC
  #   variable = disp_z
  #   boundary = top_2
  #   function = loadUnloadFunction
  # [../]
  [./top2_z]
    type = DirichletBC
    variable = disp_z
    boundary = top_2
    value = 0.0
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
    function = loadUnloadFunction
  [../]
[]
[InterfaceKernels]
  [./interface_x]
    type = CZMInterfaceKernel
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
    type = CZMInterfaceKernel
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
    type = CZMInterfaceKernel
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
  [./cohesive_law_exponential]
    type = CZMLawExponential
    displacement_jump_peak = 0.1
    traction_peak = 150
    displacement_jump_mp_name = 'displacement_jump_local'
    boundary = 'interface'
  [../]
  [./cohesive_law_unload_linear]
    type = CZMUnloadLinear
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
    unload_traction_separation_UO  = 'cohesive_law_unload_linear'
    coopenetration_penalty = 1e3
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
  dt = 1
  end_time = 142
  # dtmin = 1
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
    variable = stress_xx
    execute_on = 'initial timestep_end'
    block = 3
  [../]
  [./syy_3G]
    type = ElementAverageValue
    variable = stress_yy
    execute_on = 'initial timestep_end'
    block = 3
  [../]
  [./szz_3G]
    type = ElementAverageValue
    variable = stress_zz
    execute_on = 'initial timestep_end'
    block = 3
  [../]
  [./syz_3G]
    type = ElementAverageValue
    variable = stress_yz
    execute_on = 'initial timestep_end'
    block = 3
  [../]
  [./sxz_3G]
    type = ElementAverageValue
    variable = stress_xz
    execute_on = 'initial timestep_end'
    block = 3
  [../]
  [./sxy_3G]
    type = ElementAverageValue
    variable = stress_xy
    execute_on = 'initial timestep_end'
    block = 3
  [../]
  [./disp_3Z]
    type = ElementAverageValue
    variable = disp_z
    execute_on = 'initial timestep_end'
    block = 3
  [../]
  [./sxx_2G]
    type = ElementAverageValue
    variable = stress_xx
    execute_on = 'initial timestep_end'
    block = 2
  [../]
  [./syy_2G]
    type = ElementAverageValue
    variable = stress_yy
    execute_on = 'initial timestep_end'
    block = 2
  [../]
  [./szz_2G]
    type = ElementAverageValue
    variable = stress_zz
    execute_on = 'initial timestep_end'
    block = 2
  [../]
  [./syz_2G]
    type = ElementAverageValue
    variable = stress_yz
    execute_on = 'initial timestep_end'
    block = 2
  [../]
  [./sxz_2G]
    type = ElementAverageValue
    variable = stress_xz
    execute_on = 'initial timestep_end'
    block = 2
  [../]
  [./sxy_2G]
    type = ElementAverageValue
    variable = stress_xy
    execute_on = 'initial timestep_end'
    block = 2
  [../]
  [./disp_top3_z]
    type = SideAverageValue
    variable = disp_z
    execute_on = 'initial timestep_end'
    boundary = 'top_3'
  [../]
[]
