可以，下面我直接按你当前 `62.zip` 的现状，给一份**能落地到代码目录、阶段顺序、测试与验收都明确**的 `phase-p2-thermal` 实现方案。

先说总判断：

**你现在的仓库已经把 P2 的外壳准备好了，但真正的 P2 物理与数值主体几乎还是空的。**  
我看了 `62.zip`，当前：

- `src/thermal/` 为空
- `src/diffusion/common/`、`src/diffusion/matrix/`、`src/diffusion/hypre/` 都还没文件
- `src/state/thermo_state/` 也还没文件
- 但 `runtime` 已经有 `StageId::thermal` / `StageId::equilibration`，并且 `p2` 的 required stage mask 已经要求 `H/T/E`
- `canonical_state` 已经有 `e_electron`，cached-field 元数据里也已经预留了 `Te/Ti/Pe/Pi/conductivity` 等缓存位
- `hydro` 已经按 `chi_e -> Pe -> E_electron` 这条链把电子能量通道打通了

所以 P2 最正确的做法不是再补“更多 hydro”，而是：

## **在不动 P1 主干的前提下，补齐 `thermo recovery + diffusion infrastructure + thermal operator + equilibration operator + P2 benchmarks`**

---

# 一、P2 的目标和边界

按 Woo 论文，DEC3D 的主时间推进是  
`H -> T -> E -> R -> A`，其中 hydro 在 coarse states 上显式推进，后续热传导/电子离子交换/辐射/alpha 都在 **fine states 上**推进。热传导和电子离子交换对应方程分别是：  
\[
\partial_t\left(\frac{P_e}{\gamma-1}\right)=\nabla\cdot \kappa_e \nabla T_e,\qquad
\partial_t\left(\frac{P_i}{\gamma-1}\right)=\nabla\cdot \kappa_i \nabla T_i
\]
以及  
\[
\partial_t T_e=-\frac{1}{\tau_{ei}}(T_e-T_i),\qquad
\partial_t T_i=-\frac{1}{\tau_{ei}}(T_i-T_e)
\]
。fileciteturn27file4 fileciteturn26file4

并且论文明确说，hydro 之后的这些传输过程是**在 fine states 上、用 hydro 给出的 relaxed `Δt`、隐式求解**。fileciteturn26file4 fileciteturn28file15

所以 P2 的边界应当非常硬：

1. **只做 `T` 和 `E`**
2. `R/A` 继续不注册
3. **macro-zoning 不参与 P2**
4. **ALE 不参与 P2**
5. **thermal/equilibration 都只在 fine mesh authoritative state 上工作**
6. benchmark 先按论文那套“流体速度不更新、只看温度与能量传输”来做。论文在热传导与热传导+equilibration benchmark 里就是这样干的。fileciteturn26file11 fileciteturn27file18

---

# 二、P2 要守住的 authoritative 规则

这个部分你仓库里的 acceptance 其实已经定了，我按你当前 state 结构翻译成实现规则：

## T stage
`T` 要写：
- `E_electron`
- `E_fluid_total`

原因很直接：  
电子热传导改电子能量，离子热传导改离子能量；而 canonical state 并不单独存 `E_ion`，它把离子能量隐含在

\[
E_{ion}=E_{fluid\_total}-E_{kinetic}-E_{electron}
\]

里。所以只要离子能量变了，`E_fluid_total` 就必须跟着写。

## E stage
`E` 只写：
- `E_electron`

并且必须保持：
- `E_fluid_total` 不变（在浮点容差内）

这条很关键。  
因为 electron-ion equilibration 是**电子与离子热库之间的内部分账**，不是往系统里加能量，也不是往系统里抽能量。

---

# 三、我建议的实现顺序

不要一口气上全套。  
最稳的是按下面 6 个 workstream 来：

---

## P2-0：先把热学恢复层补出来
先做一个**不带任何 solver** 的 `thermo view`，这是所有后续模块的入口。

### 新文件
- `src/state/thermo_state/thermo_view.hpp`
- `src/state/thermo_state/thermo_view.cpp`

### 它负责恢复什么
从 canonical authoritative state 恢复：

- `rho`
- `v_r / v_theta / v_phi`
- `E_kinetic`
- `E_ion`
- `P_e`
- `P_i`
- `T_e`
- `T_i`
- `n_e`
- `n_i`

### 恢复逻辑
沿用你现在 hydro 的 `γ=5/3` 口径：

\[
P_e=(\gamma-1)E_e,\qquad
E_i=E_{fluid\_total}-E_{kinetic}-E_e,\qquad
P_i=(\gamma-1)E_i
\]

然后再由 number density closure 给出

\[
T_e = \frac{P_e}{n_e},\qquad
T_i = \frac{P_i}{n_i}
\]

### 这里最重要的工程决定
由于你当前代码还不是完整物理单位闭合系统，我建议在 P2 里把 **number-density closure 单独封装**，不要散在 thermal/equilibration 里。

比如加一个小结构：

```cpp
struct ThermoClosureOptions {
  double zbar = 1.0;
  double ion_mass_code = 1.0;
  bool use_code_unit_density_as_number_density = true;
};
```

然后把“`n_i = rho / mbar_i`, `n_e = Zbar * n_i`”或者更简化的 code-unit 版本，全都放这里。

### 红测
- `test_thermo_view_recovery.cpp`
- 重点查：
  - `E_i >= 0`
  - `P_e/P_i/T_e/T_i` 恢复正确
  - 改 `E_electron` 后缓存失效
  - 改 `E_fluid_total` 后缓存失效

---

## P2-1：补 diffusion infrastructure，但只做“数值骨架”
这层不要夹带 physics 语义。

### 新文件建议
- `src/diffusion/common/diffusion_boundary.hpp/.cpp`
- `src/diffusion/common/diffusion_diagnostics.hpp/.cpp`
- `src/diffusion/matrix/scalar_diffusion_7pt.hpp/.cpp`
- `src/diffusion/hypre/hypre_scalar_solver.hpp/.cpp`

### 这一层应该提供什么
只提供：

1. **7 点模板装配辅助**
2. **边界条件处理**
3. **Hypre bridge**
4. **matrix / rhs / solve diagnostics**

不要让它替 physics 决定：
- 这是 thermal 还是 radiation
- `A/B/C/D` 是什么
- 写哪个 authoritative field

这些必须留在 operator 自己手里。

### 3D 离散口径
Woo 论文把 3D 扩散统一写成：

\[
A^n\frac{T^{n+1}-T^n}{\Delta t}
=
\nabla\cdot D^n \nabla T^{n+1}
+
C^n T^{n+1}
+
B^n
\]

并说 thermal/radiation/alpha 都从这个通式里读 `A/B/C/D`。fileciteturn27file8

对 P2 thermal 来说，直接代入就是：

- `A = 3/2 n_s^n`
- `D = κ_s^n`
- `B = 0`
- `C = 0`

其中 `s=e/i`。这和论文 Eq. (5.11)(5.12) 是一致的。fileciteturn27file4

### Solver 选型
直接按 thesis 口径走：

- `HYPRE AMG-preconditioned GMRES`

Woo 明确写了 DEC3D 用 HYPRE 的 AMG 预条件 GMRES 来做 3D thermal/radiation/alpha 隐式扩散，而且**不做方向分裂**。fileciteturn28file14

### 红测
- `test_diffusion_7pt_row_assembly.cpp`
- `test_hypre_scalar_solver_smoke.cpp`
- `test_diffusion_boundary_zero_flux.cpp`

---

## P2-2：先上 `ThermalOperator`，把 T 路径单独跑通
这个阶段只做 thermal diffusion，不做 equilibration。

### 新文件
- `src/thermal/thermal_coefficients.hpp/.cpp`
- `src/thermal/thermal_operator.hpp/.cpp`

### operator contract
和 hydro 一样：
- `bind(...)`
- `estimate_dt(...)`
- `advance(...)`

### `estimate_dt`
P2 thermal 是隐式的，不该再给 hard cap。  
建议：
- `hard_cap_dt = inf`
- `soft_advice_dt` 可选给一个 diffusive advisory
- 或者先也返回 `inf`，只把 diagnostics 做全

### `advance()`
拆成两次独立 solve：

1. **electron solve**
2. **ion solve**

它们共享：
- same mesh
- same matrix builder
- same Hypre bridge

但系数不同：
- `A_e, κ_e, T_e`
- `A_i, κ_i, T_i`

### conductivity closure
这里我的建议很明确：

#### 第一刀
先做：
- `Spitzer`
- `LeeMore`

都留接口，但 benchmark 先用 **Lee-More + degeneracy** 这条路径作为主验证。  
因为 Woo 的 benchmark 结果明确说：
- 用 Lee-More Coulomb log + degeneracy correction，和 LILAC 对得更好
- 而把整个 `κ` 做错误的 space-lagging 是明显错的
- “simple mean approximation” 的 interface coefficient 是被验证过的。fileciteturn27file16 fileciteturn27file6

### face coefficient
这里建议你直接照 thesis benchmark 的结论来：

## **cell-interface diffusion coefficient 用 simple mean**
也就是：

\[
(\text{geom}\cdot \kappa)_{face}
=
\frac{(\text{geom}\cdot \kappa)_{L}+(\text{geom}\cdot \kappa)_{R}}{2}
\]

不要一开始就搞“把整个 κ 拉出去做 space-lagging”那条坏路。论文已经说那条是错的。fileciteturn27file16

### 边界条件
这点也别偷懒。

Woo 对 thermal diffusion 的 3D 口径是：

- **外边界**：zero heat flux
- **原点**：不是 1D/2D 那种简单零梯度，而是 full 3D condition，允许热流通过 origin
- 也就是 3D thermal 需要和你现有球坐标 ghost/remap 体系对齐。fileciteturn27file1 fileciteturn27file8

所以 thermal 的 scalar ghost 规则应是：

- `phi`：周期
- `theta` pole：标量 remap，**不做符号翻转**
- `r=0`：full 3D across-origin scalar remap
- outer edge：Neumann zero-flux

### 写回 canonical state
solve 出 `T_e^{n+1}`、`T_i^{n+1}` 后：

\[
E_e^{n+1} = \frac{3}{2} n_e^n T_e^{n+1}
\]
\[
E_i^{n+1} = \frac{3}{2} n_i^n T_i^{n+1}
\]
\[
E_{fluid\_total}^{n+1}
=
E_{kinetic}^{n}
+
E_e^{n+1}
+
E_i^{n+1}
\]

注意：
- `mom_*` 不动
- `E_kinetic` 用 stage entry 的 current momentum 算
- `updated_fields = {e_electron, e_fluid_total}`

### diagnostics
每次 solve 都必须报：
- matrix built
- rhs built
- solver invoked
- converge / fail
- iteration count
- residual 或等效证据

不然你 acceptance 里“thermal solve path lacks matrix or solve evidence”会直接红。

### 红测
- `test_thermal_write_set.cpp`
- `test_thermal_stage_matrix_evidence.cpp`
- `test_thermal_pure_slab_spitzer.cpp`
- `test_thermal_pure_slab_leemore.cpp`

---

## P2-3：再上 `EquilibrationOperator`
这个阶段不要上 HYPRE。  
电子离子交换是逐 cell ODE，本地更新就够了。

### 新文件
- `src/equilibration/equilibration_coefficients.hpp/.cpp`
- `src/equilibration/equilibration_operator.hpp/.cpp`

### 关键实现
不要直接照抄一个“更新 `T_e/T_i` 完事”的版本。  
因为你 acceptance 要求：

- `E` stage 只写 `E_electron`
- 且必须保持 `E_fluid_total` 不变

所以最稳的做法是：

### 先让温差按 thesis 衰减
论文给了：

\[
(T_e-T_i)^{n+1}=(T_e-T_i)^n \exp(-2\Delta t/\tau_{ei})
\]

fileciteturn26file1

### 再结合总热能守恒重构新温度
设：

\[
C_e=\frac32 n_e,\qquad C_i=\frac32 n_i
\]
\[
S = C_e T_e^n + C_i T_i^n
\]
\[
\Delta^{n+1} = (T_e^n-T_i^n)e^{-2\Delta t/\tau_{ei}}
\]

则：

\[
T_i^{n+1}=\frac{S-C_e\Delta^{n+1}}{C_e+C_i}
\]
\[
T_e^{n+1}=T_i^{n+1}+\Delta^{n+1}
\]

这样同时满足：

- thesis 给出的温差衰减
- 以及你项目 acceptance 要求的总流体能量守恒

然后：

\[
E_e^{n+1}=C_e T_e^{n+1}
\]

只写回 `e_electron`，  
`e_fluid_total` 保持不动。

### `tau_ei`
第一刀建议做成独立 closure helper：
- `ComputeElectronIonRelaxationTime(...)`

并把所有 thesis 外的 regularization / floor / unit choice 写进 assumption ledger。

### diagnostics
最少给：
- `p2.equilibration.executed`
- `p2.equilibration.energy_conserved`
- `p2.equilibration.max_abs_deltaTeTi_before`
- `p2.equilibration.max_abs_deltaTeTi_after`

### 红测
- `test_equilibration_preserves_e_fluid_total.cpp`
- `test_equilibration_reduces_deltaT.cpp`
- `test_equilibration_zero_dt_no_change.cpp`

---

## P2-4：runtime wiring
这一块你外壳已经有了，但还没把 operator 真挂上去。

### 要做的事
1. 在 runtime 初始化时，`p2` 真实注册：
   - `hydro`
   - `thermal`
   - `equilibration`
2. `radiation` / `alpha` 继续不注册
3. stage 顺序固定：
   - `H -> T -> E`
4. hydro 给出的 `dt` 继续作为全局步长；`T/E` 只给 advisory，不抢最终 dt

这和 Woo 的 H/T/E/R/A Lie splitting 顺序是对上的。fileciteturn26file4

### 红测
- `test_p2_runtime_required_stages.cpp`
- `test_p2_forbidden_stage_registration.cpp`
- `test_p2_stage_order.cpp`

---

## P2-5：benchmark 与 acceptance
P2 不要一上来就做复杂 3D 热斑。  
先照 thesis 把 1D benchmark 全打完。

### Benchmark 1：纯 thermal diffusion，Spitzer
论文热传导 benchmark 的 slab 初值是：

- `R = 100 μm`
- `r0 = 50 μm`
- `Te(r<r0)=5 keV`
- `Te(r>=r0)=0.5 keV`
- `rho(r<r0)=50 g/cm^3`
- `rho(r>=r0)=100 g/cm^3`

并且**不更新 fluid momentum**。fileciteturn26file0 fileciteturn26file11

### Benchmark 2：纯 thermal diffusion，Lee-More
同样 slab，但用 Lee-More Coulomb log + degeneracy。

### Benchmark 3：thermal + equilibration
同样 slab，加上 E stage。  
论文图 5.34 就是这条 benchmark。fileciteturn26file6

### Benchmark 4：1D NIF replay（manual acceptance）
最后再做 thesis 图 5.36 那种综合 1D 对照，用来验：
- hydro
- thermal
- equilibration  
一起工作。Woo 说这一组和 LILAC 对得不错。fileciteturn27file16

---

# 四、我建议的文件清单

按你当前仓库结构，P2 至少要新建这些文件：

```text
src/state/thermo_state/
  thermo_view.hpp
  thermo_view.cpp

src/diffusion/common/
  diffusion_boundary.hpp
  diffusion_boundary.cpp
  diffusion_diagnostics.hpp
  diffusion_diagnostics.cpp

src/diffusion/matrix/
  scalar_diffusion_7pt.hpp
  scalar_diffusion_7pt.cpp

src/diffusion/hypre/
  hypre_scalar_solver.hpp
  hypre_scalar_solver.cpp

src/thermal/
  thermal_coefficients.hpp
  thermal_coefficients.cpp
  thermal_operator.hpp
  thermal_operator.cpp

src/equilibration/
  equilibration_coefficients.hpp
  equilibration_coefficients.cpp
  equilibration_operator.hpp
  equilibration_operator.cpp
```

测试方面建议最少新增：

```text
tests/contract/
  test_thermo_view_recovery.cpp
  test_diffusion_7pt_row_assembly.cpp
  test_hypre_scalar_solver_smoke.cpp
  test_thermal_write_set.cpp
  test_thermal_stage_matrix_evidence.cpp
  test_equilibration_preserves_e_fluid_total.cpp
  test_equilibration_reduces_deltaT.cpp
  test_p2_runtime_required_stages.cpp

tests/acceptance/
  test_p2_thermal_slab_spitzer.cpp
  test_p2_thermal_slab_leemore.cpp
  test_p2_thermal_equilibration_slab.cpp

tests/manual/
  benchmark_p2_nif_1d_implosion.cpp
```

---

# 五、最容易踩坑的地方

## 1. 不要让 thermal/equilibration 去碰宏分区
P2 一律 fine-grid authoritative state。Woo 也是这么组织的。fileciteturn28file15

## 2. 不要把 `Te/Ti` 变成新的 authoritative truth
authoritative 仍然只能是：
- `e_fluid_total`
- `e_electron`

`Te/Ti/Pe/Pi/κ` 只是 recovered/operator-local。

## 3. E stage 绝不能偷偷改 `e_fluid_total`
这条一定要单独立红测。

## 4. thermal benchmark 里不要更新 momentum
论文 benchmark 为了看热传导与交换本身，就是把流体速度固定住。fileciteturn27file18

## 5. face coefficient 不要一开始就搞复杂 space-lagging
先用 thesis 认可的 **simple mean approximation**。fileciteturn27file16

---

# 六、我建议的提交顺序

最稳的是 5 刀，不要一刀全做：

### 刀 1
`thermo_view + recovery + cache invalidation tests`

### 刀 2
`diffusion 7pt builder + hypre smoke + diagnostics`

### 刀 3
`thermal operator (Spitzer/LeeMore) + pure slab benchmark`

### 刀 4
`equilibration operator + local energy-conserving relaxation`

### 刀 5
`runtime p2 wiring + 1D integrated NIF replay`

---

# 七、一句话版执行策略

**P2 不该从“大而全的 3D 物理集成”开始，而应该从：**
1. 把 `Te/Ti` 从 canonical state 正确恢复出来  
2. 把 3D fine-grid 隐式 diffusion 数值骨架写好  
3. 先跑纯 thermal slab  
4. 再跑 thermal+equilibration slab  
5. 最后再挂进 `H->T->E` runtime

这条路线最稳，也最符合 Woo 论文和你当前仓库状态。
