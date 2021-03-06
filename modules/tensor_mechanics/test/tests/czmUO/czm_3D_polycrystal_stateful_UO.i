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
   [./left_n1]
     type = AddExtraNodeset
     new_boundary = 'left_n1'
     nodes = '2 69 114 191 246'
     depends_on = breakmesh
   [../]
   [./left_n2]
     type = AddExtraNodeset
     new_boundary = 'left_n2'
     nodes = '7 37 25 2 45'
     depends_on = breakmesh
   [../]
   [./right_n1]
     type = AddExtraNodeset
     new_boundary = 'right_n1'
     nodes = '225 288 310 383'
     depends_on = breakmesh
   [../]
   [./right_n2]
     type = AddExtraNodeset
     new_boundary = 'right_n2'
     nodes = '228 276 296 533'
     depends_on = breakmesh
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
  [../]
[]


[BCs]
  [./left_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'left_n1 left_n2'
    value = 0.0
  [../]
  [./left_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./left_z]
    type = DirichletBC
    variable = disp_z
    boundary = left_n1
    value = 0.0
  [../]
  [./right_x]
    type = DirichletBC
    variable = disp_y
    boundary = 'right_n1'
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
    boundary = right_n1
    value = 0.0
  [../]
[]

[Functions]
  [./loadUnloadFunction]
    type = PiecewiseLinear
    x = '0 4      8  14     21      32    42   67   92 142'
    y = '0 0.002  0  0.012  -0.002   0.02   0    0.05   0   0.1'
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
    displacement_jump_peak = 0.01
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
    fill_method = symmetric_isotropic
    C_ijkl = '0.3 0.5e8'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
  [./gap]
    type = CZMUOBasedMaterial
    is_interface_material = true
    boundary = 'interface'
    displacement_jump_UO = 'displacement_jump_uo'
    traction_separation_UO = 'cohesive_law_exponential'
    unload_traction_separation_UO = 'cohesive_law_unload_linear'
    coopenetration_penalty = 100
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
  solve_type = NEWTON
  nl_abs_tol = 1e-4
  nl_rel_tol = 1e-6
  nl_max_its = 50
  l_tol = 1e-10
  l_max_its = 5
  start_time = 0.0
  dt = 1
  dtmin = 0.1
  end_time = 142
  line_search = none
[]
[Outputs]
  [./out]
    type = Exodus
  [../]
[]
[Postprocessors]
  [./sxx]
    type = ElementAverageValue
    variable = stress_xx
    execute_on = 'initial timestep_end'
  [../]
  [./syy]
    type = ElementAverageValue
    variable = stress_yy
    execute_on = 'initial timestep_end'
  [../]
  [./szz]
    type = ElementAverageValue
    variable = stress_zz
    execute_on = 'initial timestep_end'
  [../]
  [./syz]
    type = ElementAverageValue
    variable = stress_yz
    execute_on = 'initial timestep_end'
  [../]
  [./sxz]
    type = ElementAverageValue
    variable = stress_xz
    execute_on = 'initial timestep_end'
  [../]
  [./sxy]
    type = ElementAverageValue
    variable = stress_xy
    execute_on = 'initial timestep_end'
  [../]
  [./disp_rightX]
    type = SideAverageValue
    variable = disp_x
    execute_on = 'initial timestep_end'
    boundary = 'right'
  [../]
[]
