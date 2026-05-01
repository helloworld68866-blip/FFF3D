下面给你一份**“按论文可复现的 DEC3D 数值实现手册”**。它的目标不是还原作者的每一行源码，而是让一个只会 C++ 的程序员，能把论文里描述的那套 \**C++ / MPI / HYPRE\** 程序骨架重新搭出来，并且搭得**尽量贴近论文架构**。先说边界：这份手册能忠实复现论文明确写出的架构——单流体双温、球坐标 3D 欧拉、径向动网格、宏分区、HLLC+PPM 水动力、全隐式热/辐射/alpha、HYPRE 解 3D 稀疏非对称系统——但不能保证 1:1 恢复原作者源码，因为论文没有公开全部 limiter 常数、opacity/group 表、输入 deck、I/O 格式和 HYPRE 调参。你上传的论文就是这里的依据。

## 1. 先把这套程序想成什么

DEC3D 不是“一个大方程一次性全求解”的程序，而是一个**按物理模块分步推进**的 3D 球坐标并行程序。论文把它写成五个算子：

- 水动力 `H`
- 热传导 `T`
- 电子-离子能量交换 `E`
- 辐射输运 `R`
- alpha 粒子输运 `A`

整体时间推进采用 **Lie splitting**。更关键的是：**水动力只在粗网格上显式推进，用它决定放宽后的时间步长；热传导、辐射、alpha 都回到细网格上做全隐式求解。** 这就是 DEC3D 最重要的程序组织原则。

论文对 DEC3D 的总体定义也很明确：它是**单流体、双温度、全球坐标、3D 并行欧拉代码**；水动力采用 **HLLC 近似 Riemann 求解器**，再用 **PPM** 把一阶上风结果提高到三阶；热、辐射和 alpha 的隐式扩散由 **HYPRE 的 AMG 预条件 GMRES** 来解，而且 **不做方向分裂**。

------

## 2. 你应该复现的“最小真实架构”

如果你照着论文搭，最终程序应当长成这样：

### 2.1 网格层

- **细网格**：真实物理状态所在网格，尺寸 `Nr × Ntheta × Nphi`
- **粗网格**：只给显式水动力用，按每个径向层单独在 `theta/phi` 上做宏分区
- 网格是**球坐标单元中心**格式：
  - `r_i = (i - 1/2) Δr`
  - `θ_j = (j - 1/2) Δθ`
  - `φ_k = (k - 1/2) Δφ`
- `Δr = R(t)/Nr`, `Δθ = π/Ntheta`, `Δφ = 2π/Nphi`
- 外边界半径 `R(t)` 会随时间更新，满足
  `dR/dt = β v_CM^shell(t^n)`，其中 `β` 在论文里取 `1~2` 范围。

### 2.2 并行层

论文展示的球坐标版 DEC3D 是**径向子域分解**：每个 rank 持有一段 `r`，而本地持有完整的 `θ × φ` 面。这样做的直接好处是：

1. MPI 只需和左右两个径向邻居交换 halo。
2. 极点和方位角周期映射可以在 rank 内部处理。
3. 宏分区算法按“每个径向层”做，正好和这种分解方式一致。

论文还明确说 **PPM 三阶显式水动力需要 3 层 ghost cell**，而这些 radial ghost 是靠 MPI 边界交换补齐的。

### 2.3 物理状态层

从工程上，建议你把状态拆成下面几组：

```cpp
struct Grid {
    int nr_local, ntheta, nphi;
    int ng = 3;                 // PPM 需要 3 层 ghost
    int i_global_begin;         // 本 rank 对应的全局 r 起点
    double R, dr, dtheta, dphi;
    double beta;
};

struct HydroState {
    Array3D rho;                // 密度
    Array3D mom_r, mom_t, mom_p;// 三个动量分量
    Array3D E;                  // 总能量
    Array3D chi_e;              // chi_e = Pe^(3/5)，论文就是这样平流电子压强
};

struct ThermoState {
    Array3D Pe, Pi;
    Array3D Te, Ti;
};

struct RadiationState {
    std::vector<Array3D> Ug;    // 每个能群一个辐射能量密度场
};

struct AlphaState {
    Array3D eps_alpha;          // alpha 能量密度
};

struct MacroMap {
    std::vector<int> Mtheta;    // 每个 r 层的 theta 合并倍率
    std::vector<int> Mphi;      // 每个 r 层的 phi 合并倍率
};
```

这里最容易漏掉的是 `chi_e = P_e^(3/5)`。论文不是把 `Pe` 直接放进 HLLC 平流，而是把
`∂t Pe^(3/5) + ∇·(v Pe^(3/5)) = 0`
当作**标量平流方程**来处理，然后再恢复 `Pe`。这点必须照着做，否则你的程序骨架就和论文不一样了。

------

## 3. 工程目录建议

按“只会 C++ 也能维护”的方式，我建议目录直接这样分：

```text
dec3d/
  CMakeLists.txt
  app/
    main.cpp
  core/
    array3d.hpp
    grid.hpp
    indexer.hpp
    constants.hpp
  mpi/
    cart_radial_comm.hpp
    halo_exchange.hpp
  mesh/
    spherical_geometry.hpp
    moving_radius.hpp
    ghost_remap_pole_phi.hpp
    macro_zoning.hpp
    prolong_restrict.hpp
  hydro/
    conserved.hpp
    primitive.hpp
    eos.hpp
    hllc.hpp
    ppm.hpp
    moving_mesh_flux.hpp
    hydro_driver.hpp
  diffusion/
    stencil7_builder.hpp
    hypre_solver.hpp
    diffusion_problem.hpp
  physics/
    thermal_diffusion.hpp
    ei_equilibration.hpp
    radiation_group.hpp
    radiation_bc.hpp
    alpha_diffusion.hpp
  io/
    input.hpp
    vtk_output.hpp
    checkpoint.hpp
  tests/
    shock_tube.cpp
    thermal_slab.cpp
    radiation_slab.cpp
```

核心思想是：**先写一个通用 `DiffusionProblem` + `HypreSolver`，再让热传导、辐射、alpha 都只是“换参数喂给它”。** 因为论文的 3D 隐式离散正是这么设计的：统一写成一个一般扩散方程。

------

## 4. 球坐标网格、边界和 ghost 怎么做

### 4.1 单元中心

球坐标一定要用**单元中心**，这是论文明确强调的，因为这样可以避开原点和极点奇点。

### 4.2 外边界

在外径方向，论文给出的水动力边界是“零流入”类型：外层 ghost 令密度和总压强外推，径向速度按零流入处理。热扩散/alpha 扩散在外边界可先做**零通量 Neumann**。辐射边界则不要和热扩散共用，单独封装，因为论文把辐射边界单独成节处理。

### 4.3 原点与极点

这里不要硬写 “`vr=0`”“`Tghost=Tcell`” 之类拍脑袋条件。论文说得很清楚：**原点和极点不是简单边界，而是通过 `theta/phi` 的几何映射把 ghost cell 映射到内部单元。** 也就是说，你需要写一个 `ghost_remap_pole_phi()`，先做索引映射，再进入数值通量或隐式系数计算。

------

## 5. 宏分区（macro-zoning）到底怎么落地

这是论文里最“像代码设计”而不是“像物理公式”的部分。

### 5.1 为什么要做

球坐标靠近原点、靠近极点时，单元的弧长

- `ΔSθ = r_i Δθ`
- `ΔSφ = r_i sinθ_j Δφ`
  会非常小，显式 CFL 时间步被这些小角单元卡死。论文专门说：**这个问题只影响像 PPM 水动力这种显式双曲系统，不影响隐式扩散系统。** 所以才有“粗网格做水动力，细网格做扩散”的分工。

### 5.2 论文里的判据

论文给了非常明确的粗化判据：

- 当 `ΔSθ < Δr/2` 或 `ΔSφ < Δr/2` 时，该单元被视为“细单元”，需要粗化；
- 粗化只发生在角向；
- `theta` 和 `phi` 的粗化倍率分别取二的幂：
  - `Mtheta = 2^(Ltheta)`
  - `Mphi = 2^(Lphi)`
- 目标是让粗网格弧长满足：
  - `r_i Δθ Mtheta > Δr/2`
  - `r_i sinθ_j Δφ Mphi > Δr/2`

这意味着你的 `Ntheta`、`Nphi` 最好本来就取 2 的幂。

### 5.3 代码里怎么写

直接按论文的三步走：

1. **detect**：按每个径向层扫描，找出需要粗化的角向区块。
2. **generate coarse mesh**：给该径向层生成 `(Mtheta, Mphi)`。
3. **restrict / prolong**：
   - 显式水动力前：细网格 primitive → 粗网格 primitive
   - 水动力后：粗网格结果 → 回填/插值到细网格

注意：**宏分区只服务于显式 hydro 子步，不要把它套进 HYPRE 的隐式求解。** 论文的原意就是这样。

------

## 6. 水动力模块怎么搭

### 6.1 你需要实现什么

论文的水动力部分，工程上可以拆成四层：

1. **守恒变量/原始变量互转**
2. **HLLC 近似 Riemann 求解器**
3. **PPM 重构**
4. **动网格通量修正**

### 6.2 为什么是 HLLC

论文选 HLLC 的原因不是“它常见”，而是因为 3D 下单元界面上存在**切向速度**，必须把**接触间断**包含进来；HLLC 比 HLL 更完整，适合 3D Euler。

### 6.3 电子压强的处理

论文非常有辨识度的一点是：

- 单流体总压 `P` 走正常的 Euler 水动力；
- 电子压强不直接作为独立能量方程求，而是平流 `Pe^(3/5)`；
- 然后恢复 `Pe`，再用 `Pi = P - Pe` 算离子压强。

所以你的 hydro 更新后要做：

```cpp
Pe = pow(chi_e, 5.0 / 3.0);
Pi = P_total - Pe;
Te = Pe / ne;
Ti = Pi / ni;
```

### 6.4 PPM 和降阶

论文的 PPM 不是“全域都三阶硬上”。它的单调性处理思想很明确：**如果局部检测到非单调/噪声，就退回更耗散的一阶上风 HLLC；平滑区域才用高阶重构。** 你在工程里把这件事做成 `LimiterDecision` 就够了。

### 6.5 动网格

DEC3D 的网格只沿 **径向** 运动，不在 `theta/phi` 动。论文使用的是**保守型 moving-mesh 有限体积做法**，把单元面速度带进数值通量里，而且作者明确说：**把 moving-mesh 通量和 HLLC 水动力同时更新，比先单独做网格平移再做 hydro 更不耗散。**

------

## 7. HYPRE 隐式求解器怎么搭

这里是整个项目最值得一次写好的地方。

### 7.1 论文给出的统一形式

DEC3D 的 3D 隐式扩散统一写成：

# [ A^n \frac{T^{n+1} - T^n}{\Delta t}

\nabla \cdot D^n \nabla T^{n+1}
+
C^n T^{n+1}
+
B^n
]

论文明确说：热传导、辐射、alpha 都可以从这个通式读出各自的 `A,B,C,D`。

### 7.2 你真正该写的类

把 HYPRE 封成两个类就够了：

```cpp
class DiffusionProblem {
public:
    virtual void build_coefficients(...) = 0;  // 填 A,B,C,D
    virtual void apply_boundary(...) = 0;
};

class HypreSolver {
public:
    void assemble_matrix_7pt(...);
    void assemble_rhs(...);
    void solve_gmres_boomeramg(...);
};
```

### 7.3 7 点模板怎么装

论文给出的离散形式本质上就是一个**球坐标单元中心 7 点模板**：中心点 + 六个邻居 `r±, θ±, φ±`。如果你把式子整理成标准线性系统：

[
a_P T_P^{n+1}

- a_E T_E^{n+1}
- a_W T_W^{n+1}
- a_N T_N^{n+1}
- a_S T_S^{n+1}
- a_T T_T^{n+1}
- a_B T_B^{n+1}
  = rhs
  ]

那么按论文的记号，可以直接写成：

```cpp
aP = 1 + Zpr + Zmr + Zpt + Zmt + Zpp + Zmp - dt / A * C;
aE = -Zpr;  aW = -Zmr;
aN = -Zpt;  aS = -Zmt;
aT = -Zpp;  aB = -Zmp;
rhs = Told + dt / A * B;
```

这就是你组装 HYPRE 稀疏矩阵时每一行该填的东西。

### 7.4 HYPRE 该怎么选

论文的结论很清楚：

- 这些 3D 隐式系统是**大规模、稀疏、非对称**；
- 尤其辐射，因为 opacity 空间变化很强，更非对称；
- 作者最后用的是 **AMG 预条件的 GMRES**；
- 而且 **不做方向分裂**。

所以你的第一版直接用：

- `HYPRE_IJMatrix / HYPRE_IJVector`
- 转成 `ParCSR`
- 线性求解器：`ParCSRGMRES`
- 预条件器：`BoomerAMG`

就已经和论文架构非常接近了。工程上别先折腾 BiCGStab、SOR、多重网格自己写版本。

------

## 8. 各物理模块怎么往这个通用隐式框架里塞

### 8.1 电子/离子热传导

直接复用通用扩散框架：

- 电子：`T = Te`, `A = 3/2 ne`, `D = κe`, `C=0`, `B=0`
- 离子：`T = Ti`, `A = 3/2 ni`, `D = κi`, `C=0`, `B=0`

论文的热传导是 Spitzer-Härm 型，并在 3D 里采用**全隐式、空间二阶**。

### 8.2 电子-离子能量交换

这个模块反而最简单。它不是空间耦合方程，所以**不要上 HYPRE**。论文给了一个点值的指数衰减形式：

[
(Te - Ti)^{n+1} = (Te - Ti)^n \exp(-2\Delta t/\tau_{ei})
]

工程上就逐单元 local update。

### 8.3 多群辐射

辐射模块要分两半：

1. **流体平流 + (P dV)**：这部分在 hydro 子步里处理；
2. **扩散 + 发射/吸收刚性源项**：这部分在隐式步里处理。

论文特别强调：辐射在热点区会很快向外扩散，而且 opacity 空间变化强，因此必须全隐式，还需要 flux limiter 处理光学薄热点里的“扩散流不能大于真实粒子流”的问题。

把每个能群都看成一个 `DiffusionProblem` 即可：

- `T = Ug`
- `A = 1`
- `D = Dg`
- `C = - c κg^P`
- `B = c κg^P Bg`

然后再把各群与物质的净交换累加回电子能量。

### 8.4 alpha 粒子

论文的 alpha 模型同样可以直接塞进通用扩散器：

- `T = ε_alpha`
- `A = 1`
- `D = D_alpha`
- `C = -1 / τ_alpha,e`
- `B = n_D n_T <σv>_DT E_alpha0`

而且论文明确把它和热传导、辐射一起归入“全隐式扩散型方程”。

------

## 9. 主时间循环，按这个写就对了

这是最接近“源码骨架”的版本：

```cpp
while (t < t_end) {
    // 0. 更新外边界半径 R(t)
    update_outer_radius(grid, shell_cm_velocity, dt_guess);

    // 1. 先把角向/极点 ghost 映射补齐，再做径向 MPI halo
    remap_pole_and_phi_ghosts_local(grid, state);
    exchange_radial_halo_mpi(grid, hydro_state, /*ng=*/3);

    // 2. 宏分区：只给 hydro 用
    detect_fine_cells_and_build_macro_map(grid, macro_map);
    restrict_fine_to_coarse(grid, macro_map, fine_primitive, coarse_primitive);

    // 3. 用粗网格做显式水动力，拿到 relaxed dt
    dt = compute_cfl_on_coarse_mesh(grid, macro_map, coarse_primitive);
    hydro_hllc_ppm_moving_mesh_step(grid, macro_map, coarse_state, dt);

    // 4. 把 coarse 结果回填到 fine
    prolongate_coarse_to_fine(grid, macro_map, coarse_state, fine_state);

    // 5. 恢复热力学量
    recover_pressures_and_temperatures(fine_state);

    // 6. 细网格上做隐式热传导
    solve_electron_thermal_diffusion_hypre(...);
    solve_ion_thermal_diffusion_hypre(...);

    // 7. 电子-离子交换（逐点）
    electron_ion_relaxation_local(...);

    // 8. 多群辐射：逐群隐式
    for (int g = 0; g < ngroups; ++g) {
        solve_radiation_group_hypre(g, ...);
    }
    update_electron_energy_from_radiation(...);

    // 9. alpha 粒子隐式扩散
    solve_alpha_diffusion_hypre(...);
    deposit_alpha_heating_to_electrons(...);

    // 10. 诊断与输出
    compute_diagnostics(...);
    write_output(...);

    t += dt;
}
```

这段主循环背后的论文依据只有两句：**hydro 在 coarse states 上显式推进以放宽 `dt`；扩散型模块在 fine states 上全隐式推进，并使用这个 relaxed `dt`。**

------

## 10. 真正的实现顺序：别一上来就全做

最稳的复现路线是这样的：

### 第一步：串行、球坐标、静网格、纯 hydro

先做：

- cell-centered spherical grid
- conserved/primitive conversion
- HLLC
- 一阶版本先跑通

目标：shock tube 和 1D 球壳测试不炸。

### 第二步：加 PPM 和单调性降阶

先不加 MPI，不加辐射。把 HLLC + PPM 在球坐标上跑稳。论文专门用 shock tube 验证了这件事。

### 第三步：加动网格

只让 `r` 方向动，通量里带面速度。到这一步，hydro 骨架才和 DEC3D 接近。

### 第四步：加径向 MPI 分解和 3 层 halo

这是论文球坐标版最接近原架构的并行方式。先只交换 hydro 变量。

### 第五步：加宏分区

先写 `detect -> coarse map -> restrict -> prolongate`，只对 hydro 生效。

### 第六步：写通用 HYPRE 扩散器

先用一个最简单的纯扩散 benchmark 跑通 7 点模板和 GMRES+BoomerAMG。然后再挂热传导。

### 第七步：加电子-离子交换

这一步最简单，逐点更新。

### 第八步：加多群辐射

最后再做 flux limiter、发射吸收、群循环和电子能量回写。辐射是最难的，不要提前。

### 第九步：加 alpha

复用同一个 HYPRE 框架。

------

## 11. 你至少要做的三个 benchmark

论文的验证思路，翻译成你的工程任务，就是这三个：

1. **shock tube**
   - 验 HLLC
   - 验 PPM
   - 验降阶是否抑制噪声
2. **纯热扩散 / 热交换 slab**
   - 验 7 点隐式装配
   - 验 HYPRE 求解器
   - 验 Te/Ti 交换是否稳定
3. **clean 1D/3D stagnation profile**
   - 验整体 hydro + thermal 的耦合
   - 看密度、压强、速度、电子温度轮廓是否合理
     论文用 LILAC 做过这类对比。

------

## 12. 最容易写错的地方

### 12.1 把所有模块都放在细网格上显式推进

这样你会立刻被极点和原点的小角单元卡死，完全背离 DEC3D 的宏分区思想。

### 12.2 直接在球坐标原点写死边界条件

论文的处理不是“原点零通量”这种简单 BC，而是 ghost cell 的几何映射。

### 12.3 给热、辐射、alpha 各写一套不同的线性求解器

这会把工程复杂度炸掉。论文其实已经给你统一成了一个一般扩散方程。

### 12.4 先做辐射再做 hydro

顺序别反。DEC3D 的骨架是：**coarse hydro 决定 dt，fine implicit modules 用这个 dt。**

------

## 13. 一句话总结成“程序员语言”

你要复现的不是一个“复杂物理软件”，而是一个非常清晰的程序模式：

**球坐标细网格保存真实状态 → 每步按径向层做宏分区 → 粗网格上用 HLLC+PPM+径向动网格做显式 hydro → 细网格上把热传导/辐射/alpha 全都改写成统一的 7 点隐式线性系统 → 交给 HYPRE 的 AMG+GMRES 去解 → 电子-离子交换逐点更新。** 这就是论文里的 DEC3D 骨架。

你要是愿意，我下一条可以直接给你补一版**“可开工的 CMake + 类定义 + main.cpp 骨架”**。
