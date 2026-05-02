附件是我现有项目代码，我希望增加一个纯二维模拟mode，专门用于模拟只有轴对称扰动的情形，应该怎么做？给我详细的实现方案

我会先快速梳理压缩包里的目录和主时间推进/网格/输出结构，再给出一个“最小侵入”的二维轴对称 mode 方案：数据结构、输入参数、循环边界、算子改动、验证路径和分阶段落地顺序。

下面是我建议的实现方案。核心判断是： **不要另写一套 `Array2D + 2D Hydro + 2D Diffusion` 分叉代码** 。你现在的代码已经有球坐标 `Array3D(r,theta,phi)` 、m=0 Legendre 扰动初始化、方向开关式 hydro sweep、HYPRE 分布式扩散组装。最小侵入、风险最低的“纯二维轴对称 mode”应当定义为：

> **内部仍使用 `Array3D` ，但令 `phi_cells = 1` ，该唯一 φ 单元代表完整 $0\\sim2\\pi$ 方位角；所有 φ 方向通量、φ sweep、φ HYPRE coupling 全部关闭，只在 `(r,\theta)` 上演化。**

这样得到的是数学上的轴对称二维球坐标模拟，而不是“3D 里放一个 m=0 扰动然后复制很多 φ 面”。

论文/项目依据也支持这个方向：DEC2D/DEC3D 都有 spherical mesh、conservative moving mesh、radiation/alpha 模块，而 DEC3D 的 3D 化主要是并行和 HYPRE 扩散求解扩展；球坐标采用 cell-centered 网格以避开奇点，并且 moving mesh 只沿径向运动。 论文里的 single-mode database 也明确区分了 `m=0` 的 2-D modes 和 `m≠0` 的 3-D modes，所以你这个需求在物理语义上就是“DEC3D 代码里的 DEC2D/axisymmetric runtime path”。

---

## 1\. 目标行为定义

新增一个输入参数，例如：

```markdown
[mesh]
dimensionality = axisymmetric_2d
radial_cells = 256
theta_cells = 128
phi_cells = 1
geometry = spherical
moving_mesh = false
macro_zoning = true

[perturbation]
enabled = true
type = single_mode_radial_velocity
ell = 20
m = 0
target = radial_velocity_cm_s
```

这个 mode 的强约束是：

1. `phi_cells` 必须等于 `1` 。
2. `perturbation.m` 必须等于 `0` 。
3. 只允许轴对称扰动，例如 $P\_\\ell(\\cos\\theta)$ 或 $Y\_\\ell^0$，不能允许 `m != 0` 。
4. `mom_phi` / `v_phi` 必须初始化为 0，并在每步诊断中检查它是否保持在 roundoff 级别。
5. Hydro 只做 `r` 和 `theta` sweep，不做 `phi` sweep。
6. Diffusion matrix 只含 `r±` 和 `theta±` coupling，不含 `phi±` coupling。
7. 唯一 φ 单元的几何权重必须是完整 `2π` ，即体积、面积、总质量、总能量、产额都代表完整旋转体，而不是一个单位弧度 wedge。

这比现在简单设置 `m=0, phi_cells=32` 更纯粹：现在的做法仍然存了 32 个完全相同的 φ 面，还会浪费 hydro/diffusion/output 成本；新 mode 才是真正的二维自由度。

---

## 2\. 总体架构：保留 Array3D，不要新建 Array2D

你的 `CanonicalState` 已经把所有权威状态存在 `Array3D` 里：

```markdown
rho(r, theta, phi)
mom_r(r, theta, phi)
mom_theta(r, theta, phi)
mom_phi(r, theta, phi)
e_fluid_total(r, theta, phi)
e_electron(r, theta, phi)
radiation_groups[g](r, theta, phi)
alpha_state(...)
```

建议第一版不要动这个设计，只让 `phi` extent 变成 1。原因：

- 绝大多数 kernel 可复用。
- HYPRE row ownership 会自然变成 `Nr * Ntheta * 1` 。
- 输出、checkpoint、diagnostics 不需要全部重写。
- 后续若 profile 显示 `Array3D(phi=1)` 仍有瓶颈，再引入真正 `Array2D` 做内存布局优化。

换句话说，第一版目标不是“代码类型系统二维化”，而是“PDE 自由度二维化”。

---

## 3\. 输入 deck 和配置层改动

### 3.1 新增枚举

在 `src/io/input_deck.hpp` 里给 `MeshConfig` 加一个字段：

```markdown
enum class Dimensionality {
    Full3D,
    Axisymmetric2D
};

struct MeshConfig {
    Dimensionality dimensionality = Dimensionality::Full3D;
    ...
};
```

也可以先用字符串，但 enum 更安全。

### 3.2 修改 parser 校验

你现在 parser 会要求 `phi_cells` 是偶数，所以 `phi_cells=1` 会失败。需要把逻辑改成：

```markdown
if (cfg.mesh.dimensionality == Dimensionality::Axisymmetric2D) {
    if (cfg.mesh.geometry != "spherical") {
        return error("axisymmetric_2d requires spherical geometry");
    }
    if (cfg.mesh.phi_cells != 1u) {
        return error("axisymmetric_2d requires phi_cells = 1");
    }
} else {
    if (cfg.mesh.phi_cells == 0u || cfg.mesh.phi_cells % 2u != 0u) {
        return error("full_3d requires positive even phi_cells");
    }
}
```

扰动校验也要加硬约束：

```markdown
if (cfg.mesh.dimensionality == Dimensionality::Axisymmetric2D &&
    cfg.perturbation.enabled &&
    cfg.perturbation.m != 0) {
    return error("axisymmetric_2d only supports m = 0 perturbations");
}
```

你现有 `PerturbationConfig` 已经基本只支持 `m=0` single-mode radial velocity，所以这里是把已有事实升级成显式 mode 合约。

### 3.3 runtime report 必须写清楚

输出 report 里加：

```markdown
mesh_dimensionality=axisymmetric_2d
active_hydro_directions=r,theta
phi_cells=1
azimuthal_weight=2pi
phi_sweep_executed=false
```

这很重要。否则以后看结果时很容易把 “2D axisym” 和 “3D m=0 replicated” 混掉。

---

## 4\. 球坐标几何层

你的 `BuildSphericalGeometry` 已经固定使用：

```markdown
theta ∈ [0, π]
phi ∈ [0, 2π]
```

所以当 `phi_cells=1` 时，唯一 φ face 区间自然就是 `[0,2π]` 。这正是二维轴对称所需的体积权重。这里不要改成 `Δφ = 1` 或 `Δφ = 2π / Nphi_eff` 之类，否则总质量、总能量、产额都会错。

需要加一个测试：

```markdown
// axisymmetric geometry volume test
sum(cell_volume over r,theta,phi=0)
≈ 4/3 * π * R^3
```

这样可以证明 `phi_cells=1` 代表完整球体，而不是一个 wedge。

---

## 5\. 初始化层

你的初始化里已经有轴对称单模：

```markdown
AxisymmetricSingleMode(theta) = LegendreP(ell, cos(theta))
```

这很好。第一版只需要加三件事：

### 5.1 强制 m=0

axisymmetric mode 下不能接受 `m != 0` 。

### 5.2 强制 mom\_phi = 0

初始化 profile 如果给了 `v_phi` 或 `mom_phi` 非零，建议直接报错，而不是偷偷清零：

```markdown
if (axisymmetric_2d && max_abs(mom_phi) > tolerance) {
    return error("axisymmetric_2d requires zero azimuthal velocity");
}
```

如果只是 roundoff，则可以清零并写 diagnostic。

### 5.3 初始化报告写扰动公式

例如：

```markdown
initial_perturbation_formula=vr += A*abs(vr_r0)*f(r)*P_l(cos(theta))
initial_perturbation_m=0
```

这样后续看 log 能确认不是误跑了 3D spherical harmonic。

---

## 6\. 边界 remap：这是最容易踩坑的地方

现有球坐标边界里，pole/origin remap 通常需要做半圈 φ 映射：

```markdown
phi -> phi + π
```

所以当前代码要求 `phi_cells` 为偶数。axisymmetric 2D 下 `phi_cells=1` 时，半圈映射仍然合法，因为所有 φ 方向值相同，映射回唯一 φ cell 即可。

### 6.1 hydro pole/origin 映射

修改类似 `MapPhiAcrossPole(phi, phi_cells)` 的函数：

```markdown
std::size_t MapPhiAcrossPole(std::size_t phi, std::size_t phi_cells) {
    if (phi_cells == 1u) {
        return 0u;
    }
    if (phi_cells == 0u || phi_cells % 2u != 0u) {
        throw std::invalid_argument("phi_cells must be one or even");
    }
    return (phi + phi_cells / 2u) % phi_cells;
}
```

现有向量符号变换原则可以保留：

- 穿过 `theta` pole： `mom_theta` 变号， `mom_phi` 变号。
- 穿过 `r=0` origin： `mom_r` 变号， `mom_phi` 变号。
- `mom_phi=0` 时这些变换不会产生非轴对称流。

### 6.2 scalar remap

`src/mesh/boundary/spherical_scalar_remap.*` 也要做同样处理。否则 hydro 可以跑，但热传导/辐射/alpha 的 HYPRE assembly 在 origin/pole remap 校验时仍会因为 `phi_cells=1` 失败。

关键原则：

```markdown
if (phi_cells == 1u) {
    half_turn_phi = 0u;
}
```

而不是关闭 pole/origin remap。轴对称二维仍然需要正确的球坐标 pole/origin parity。

---

## 7\. Hydro 改动

你的 hydro options 已经有方向开关，这是实现二维 mode 的最好入口：

```markdown
StaticGridHydroOptions {
    apply_radial_sweep;
    apply_theta_sweep;
    apply_phi_sweep;
    use_ppm_reconstruction;
    use_macro_zoning;
    ...
}
```

在 runtime options 构造处加：

```markdown
if (config.mesh.dimensionality == Dimensionality::Axisymmetric2D) {
    options.apply_radial_sweep = true;
    options.apply_theta_sweep = true;
    options.apply_phi_sweep = false;
    options.macro_zoning_use_phi_ppm = false;
}
```

### 7.1 CFL 估计

CFL 也必须只看 active directions。即：

```markdown
dt = min(dt_r, dt_theta);
```

不要让 `dt_phi` 参与。 `phi_cells=1` 时 `Δφ=2π` 很大，通常不会限制时间步，但为了语义清晰，还是应该明确排除。

### 7.2 macro-zoning

axisymmetric 2D 下仍然可以保留 θ 向 macro-zoning，尤其靠近 pole 的小角单元仍可能限制显式 hydro。

但 `Mphi` 应固定为 1：

```markdown
if (axisymmetric_2d) {
    macro_map.Mphi = 1;
}
```

宏分区逻辑可以继续扫描 `r,theta` ，只是不再做 φ 粗化，也不再做 φ coarse sweep。

### 7.3 moving mesh / ALE

论文里 DEC2D/DEC3D spherical moving mesh 本来就只沿径向动，这和轴对称二维完全兼容。

但你代码里有一个 direct macro-ALE-HLLC MPI path，可能内部默认三方向都存在。建议分两步落地：

第一步先支持：

```markdown
axisymmetric_2d + moving_mesh = false
```

跑通所有 hydro/diffusion。

第二步再打开：

```markdown
axisymmetric_2d + moving_mesh = true
```

并确保 ALE path 同样跳过 φ sweep，只保留径向 moving-face flux 和 θ flux。

---

## 8\. HYPRE 扩散层

这里你的现有代码很接近可用状态。分布式 diffusion assembly 如果发现：

```markdown
global_phi_cells <= 1
```

就应该不装配 φ 邻居 coupling。这样每个 scalar equation 的未知数自然变成：

```markdown
Nr * Ntheta * 1
```

也就是二维轴对称矩阵。

需要重点检查两点：

1. scalar remap 不再拒绝 `phi_cells=1` 。
2. matrix stencil 没有偷偷给唯一 φ cell 加 periodic self-neighbor。

目标矩阵结构应是：

```markdown
center
r_minus, r_plus
theta_minus, theta_plus
```

没有：

```markdown
phi_minus, phi_plus
```

热传导、辐射、alpha 都走同一个 diffusion 框架。论文里这些 diffusion-type equations 使用 hydro 得到的 relaxed timestep，再隐式推进，这个结构在二维 mode 下不需要改变。

---

## 9\. 辐射、热传导、alpha 模块

这三类物理模块不应该各自写二维版本。

只要 diffusion operator 正确降为 `(r,theta)` ，它们自然成为二维轴对称版本：

- 电子热传导： `Te(r,theta)`
- 离子热传导： `Ti(r,theta)`
- 多群辐射： `Ug(r,theta)`
- alpha 能量密度： `eps_alpha(r,theta)`

辐射的物质耦合、Planck/Rosseland opacity、group source 等都仍然逐 cell 计算，不涉及 φ 导数，因此无需特殊改动。论文中的多群辐射方程本身就是守恒/扩散形式，加入流体 advection 和 work terms 后仍可按方向项关闭 φ 项。

---

## 10\. 输出和 checkpoint

第一版建议保持现有 3D 输出格式，但 `phi_cells=1` 。例如：

```markdown
output_00010.csv3d
r,theta,phi,rho,...
...
```

其中 `phi` 只有一个值。这样不会破坏现有 post-processing。

同时新增一个可选 2D 输出：

```markdown
[output]
axisymmetric_2d_csv = true
```

输出列：

```markdown
r, theta, rho, vr, vtheta, Te, Ti, P, ...
```

checkpoint/restart 里必须写：

```markdown
dimensionality=axisymmetric_2d
phi_cells=1
```

并禁止把 axisymmetric checkpoint 用 full3d deck 直接读入，除非你专门写了 expand-to-3D 工具。

---

## 11\. 每步必须加的 invariant diagnostics

axisymmetric mode 下建议每个 major stage 后记录：

```markdown
max_abs_mom_phi
max_abs_v_phi
phi_sweep_executed
hypre_phi_couplings
global_scalar_unknowns
total_mass
total_energy
```

硬性判断：

```markdown
if (max_abs_mom_phi > 1e-12 * characteristic_momentum) {
    fail_or_warn("axisymmetric invariant violated: nonzero mom_phi");
}
```

也可以允许 roundoff 自动清零：

```markdown
if (max_abs_mom_phi < tiny) {
    mom_phi = 0;
}
```

但不要默默清掉大值，因为那会掩盖真正的数值错误。

---

## 12\. 推荐测试清单

### 12.1 input deck 测试

新增：

```markdown
axisymmetric_2d accepts phi_cells=1
axisymmetric_2d rejects phi_cells=2
axisymmetric_2d rejects perturbation m=1
full_3d rejects phi_cells=1
```

### 12.2 几何测试

```markdown
sum cell volumes for Nr,Ntheta,phi=1 equals 4/3*pi*R^3
```

### 12.3 hydro ghost 测试

覆盖：

```markdown
theta pole remap works with phi_cells=1
origin radial remap works with phi_cells=1
mom_theta/mom_r signs are correct
mom_phi remains zero
```

### 12.4 hydro sweep 测试

跑一个小 case：

```markdown
Nr=16, Ntheta=16, Nphi=1
H only
```

断言：

```markdown
apply_phi_sweep=false
phi_sweep_executed=false
state extents are 16x16x1
```

### 12.5 diffusion/HYPRE 测试

断言矩阵行数：

```markdown
global_rows = Nr * Ntheta
```

断言 stencil 没有 φ 邻居。

### 12.6 与 3D replicated m=0 对比

跑两个小 case：

```markdown
A: full_3d, phi_cells=8, m=0
B: axisymmetric_2d, phi_cells=1, m=0
```

比较同一时间的 φ 平均：

```markdown
rho_2d(r,theta) ≈ mean_phi rho_3d(r,theta,phi)
vr_2d(r,theta)  ≈ mean_phi vr_3d(r,theta,phi)
Te_2d(r,theta)  ≈ mean_phi Te_3d(r,theta,phi)
```

如果 full3d 里 φ sweep 对完全 φ-invariant 状态只产生 roundoff，二者应非常接近。这个测试是最重要的物理回归。

---

## 13\. 分阶段落地顺序

### Phase A：配置 + 几何 + 初始化

目标：程序能读入：

```markdown
dimensionality = axisymmetric_2d
phi_cells = 1
m = 0
```

并完成初始化、输出初始场。

验收：

```markdown
total_mass 正确
mom_phi = 0
volume sum 正确
```

### Phase B：边界 remap

目标：pole/origin remap 允许 `phi_cells=1` 。

验收：

```markdown
theta pole ghost test pass
origin ghost test pass
scalar remap test pass
```

### Phase C：hydro 二维推进

目标：hydro 只跑 `r,theta` 。

验收：

```markdown
H-only axisymmetric case 跑若干步不炸
phi_sweep_executed=false
mom_phi remains zero
```

### Phase D：隐式模块

目标： `T/E/R/A` 全部支持 `phi_cells=1` 。

验收：

```markdown
H+T+E smoke pass
H+T+E+R smoke pass
H+T+E+R+A smoke pass
HYPRE rows = Nr*Ntheta
```

### Phase E：对比验证

目标：axisymmetric 2D 和 full3D m=0 replicated case 对比。

验收：

```markdown
L1/Linf difference within tolerance
mass/energy histories agree
RT growth history agree for P_l(cos theta)
```

---

## 14\. 你最终应该避免的两个错误实现

第一，不要把二维 mode 做成：

```markdown
phi_cells = 2 或 4
```

然后说这是二维。那仍然是 3D 存储和 3D φ 拓扑，只是很粗。

第二，不要把唯一 φ 单元当成一个小 wedge。axisymmetric 2D 的唯一 φ 单元必须代表完整 `2π` ，否则总质量、总能量、yield、burn average 都会少一个 `2π` 因子。

---

## 15\. 最小代码改动总表

| 区域                           | 改动                                                         |
| ------------------------------ | ------------------------------------------------------------ |
| `input_deck.hpp/cpp`           | 加 `dimensionality=axisymmetric_2d` ；允许该 mode 下 `phi_cells=1` ；强制 `m=0` |
| `spherical_mesh.cpp`           | 基本不用改；加 volume test，确认 `phi_cells=1` 覆盖 `2π`     |
| `profile_initializer.cpp`      | 保持 `P_l(cosθ)` ；axisym 下强制 `mom_phi=0`                 |
| `theta_pole_boundary.cpp`      | `phi_cells==1` 时 half-turn map 返回 0                       |
| `radial_boundary_contract.cpp` | origin half-turn 复用 `phi_cells==1 -> 0`                    |
| `spherical_scalar_remap.cpp`   | scalar half-turn validation 允许 `phi_cells==1`              |
| `dec3d_app.cpp`                | axisym 下设置 `apply_phi_sweep=false` 、 `macro_zoning_use_phi_ppm=false` |
| `static_grid_hydro.cpp`        | 确认所有 φ sweep 都受 `apply_phi_sweep` 控制                 |
| ALE/MPI hydro path             | moving mesh path 也要跳过 φ sweep                            |
| HYPRE diffusion                | 确认 `global_phi_cells<=1` 不装 φ coupling                   |
| output/checkpoint              | 记录 dimensionality；可选 2D CSV 输出                        |
| tests                          | parser、geometry、ghost、hydro、diffusion、3D parity regression |

---

一句话总结： **把 axisymmetric 2D 做成 `spherical Array3D(r,theta,phi=1)` 的一等 runtime mode；唯一 φ cell 覆盖完整 `2π` ，关闭所有 φ 方向 PDE 操作，保留 r/theta hydro、radial moving mesh、HYPRE 隐式扩散和现有状态结构。** 这条路线最符合你现在代码结构，也最容易和论文里的 DEC2D/DEC3D 物理设定对齐。