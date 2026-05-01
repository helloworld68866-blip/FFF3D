**先把完整 PPM 做完，再碰 MPI**，这是对的。

而且你现在这版代码里，差的东西其实非常明确，不是泛泛地“再高阶一点”：

- 现在的 `ReconstructPpmLine(...)` **没有 `dt`、没有 `dx`、没有特征分解**。
- 它在成功时直接写的是  
  `result.interfaces.push_back({face_state, face_state, false});`  
  也就是**左右界面态完全相同**。
- `BuildPpmInterfaceFluxes(...)` 也**没有把时间步和线几何传进 reconstruction**，所以它从 API 上就不可能变成 traced-interface PPM。

所以你现在要做的，不是“继续往现有 face-value 壳上补 patch”，而是把它改成真正的：

**direction-local characteristic + parabolic cell profile + traced left/right interface states + troubled-cell first-order fallback**

这 4 件事一旦做完，PPM 就算真收住了；MPI 可以后置。

---

## 先定义：对你这套代码来说，什么叫“完整 PPM”

你现在不用把目标定成“把论文里所有花活一次做满”，否则很容易失控。

对你当前这个仓库，**full PPM closeout** 的最小正确含义是：

1. **在 direction-local primitive 上建立局部 1D 特征基**
   变量还是你现在这组：
   \[
   W_{\text{dir}} = (\rho, v_n, v_{t1}, v_{t2}, P, \chi_e)
   \]

2. **在特征变量上做抛物重构**
   不是 face 上直接算一个 `7/12` 值就完事，而是对每个 cell 建：
   - 左边值 `q_L`
   - 右边值 `q_R`
   - 抛物曲率 `q_6`

3. **按特征速度做 tracing**
   真正得到：
   - `W_{i+1/2}^L`
   - `W_{i+1/2}^R`
   
   这两个态一般**不应该相同**。

4. **troubled cell/interface 直接降到 first-order HLLC**
   不能 clip 一下继续冒充高阶。

你先把这 4 个做完，就已经足够把当前 ledger 里那句  
“not full characteristic / traced-interface PPM closeout”  
真正收掉了。

---

## 最重要的实现判断：先别加新 option，直接把现有 PPM 实现升级

不要再造一个 `PPM2`、`full_ppm`、`characteristic_ppm` 之类的新用户选项。

现在最稳的做法是：

- **保留** `use_ppm_reconstruction`
- **保留** `reconstruction_ghost_layers`
- **升级** `ppm_reconstruction.*` 的内部实现
- `static_grid_hydro.cpp` 只负责把 `dt` 和 line metric 喂进去

这样以后做 MPI，不会再为“选项语义分叉”付二次代价。

---

# 具体怎么做

## 第一步：先把 reconstruction API 改对

这是第一刀，而且必须先动。

你现在的签名大概是：

```cpp
PpmLineReconstructionResult ReconstructPpmLine(
    const std::vector<HydroConservativeState>& ghosted_line,
    HydroDirection direction,
    std::size_t interior_cells,
    std::size_t ghost_layers) noexcept;
```

这不够，因为 full PPM 一定需要 `dt/dx`。

我建议直接改成：

```cpp
PpmLineReconstructionResult ReconstructPpmLine(
    const std::vector<HydroConservativeState>& ghosted_line,
    HydroDirection direction,
    std::size_t interior_cells,
    std::size_t ghost_layers,
    double dt_s,
    const std::vector<double>& effective_cell_widths) noexcept;
```

要求：

- `effective_cell_widths.size() == ghosted_line.size()`
- 每个线单元都能拿到自己的 `dx_eff`
- reconstruction 内部自己算 `sigma = |lambda| dt / dx_eff`

这一改非常关键。  
**不把 `dt` 和 `dx` 放进 reconstruction，后面所有“traced-interface PPM”都只是名字像。**

---

## 第二步：把当前 face-level 实现拆成“cell profile -> traced interface”

你现在是“对 face 直接出一个 face_value”。  
完整 PPM 必须改成“先对 cell 建 profile，再从 profile trace 到 face”。

我建议在 `ppm_reconstruction.cpp` 内部引入 3 个内部结构体，放 `.cpp` 里，不一定暴露到头文件：

```cpp
struct CharacteristicBasis {
  std::array<double, 6> lambda{};
  std::array<std::array<double, 6>, 6> R{};
  std::array<std::array<double, 6>, 6> L{};
  bool success{false};
};

struct PpmCellProfile {
  CharacteristicBasis basis;
  std::array<double, 6> c_center{};
  std::array<double, 6> c_left{};
  std::array<double, 6> c_right{};
  std::array<double, 6> c6{};
  DirectionalPrimitiveState w_left_edge{};
  DirectionalPrimitiveState w_right_edge{};
  bool troubled{false};
};

struct PpmExecutionStats {
  std::size_t troubled_cell_count{0};
  std::size_t downgraded_interface_count{0};
  std::size_t characteristic_basis_failures{0};
};
```

思路是：

1. `ghosted_line` 先全转成 `DirectionalPrimitiveState`
2. 对每个能建 profile 的 cell，做：
   - 特征基
   - 特征投影
   - 抛物 profile
   - edge monotonicity reset
   - primitive 物理性检查
3. 最后再按 face 装配 `left/right traced states`

---

## 第三步：实现 6 变量的局部特征基

你的方向局部变量是：

\[
W = [\rho, u, v, w, p, \chi]^T
\]

其中：

- `u = v_n`
- `v,w` 是两切向速度
- `chi = chi_e`

对应的特征速度就是：

\[
\lambda = \{u-c,\; u,\; u,\; u,\; u,\; u+c\}
\]

这里 4 个 `u` 分别对应：

- 接触波
- 两个切向剪切波
- `chi_e` 标量波

右特征向量你可以先用这一组，足够干净：

- acoustic `u-c`：
  \[
  r_- = [1,\; -c/\rho,\; 0,\; 0,\; c^2,\; 0]^T
  \]
- contact：
  \[
  r_c = [1,\; 0,\; 0,\; 0,\; 0,\; 0]^T
  \]
- shear-1：
  \[
  r_{t1} = [0,\; 0,\; 1,\; 0,\; 0,\; 0]^T
  \]
- shear-2：
  \[
  r_{t2} = [0,\; 0,\; 0,\; 1,\; 0,\; 0]^T
  \]
- scalar：
  \[
  r_{\chi} = [0,\; 0,\; 0,\; 0,\; 0,\; 1]^T
  \]
- acoustic `u+c`：
  \[
  r_+ = [1,\; c/\rho,\; 0,\; 0,\; c^2,\; 0]^T
  \]

### 这里有个非常实用的建议

**左特征矩阵 `L` 先不要手推硬编码，先数值求逆。**

也就是：

- 先拼 `R`
- 再做一个小的 `6x6` Gauss-Jordan inverse
- 求 `L = R^{-1}`

原因很简单：

- 你现在最缺的是“先做对”，不是“先做到极致快”
- 6x6 数值逆在这个阶段完全够用
- 还能顺手做失败保护：一旦矩阵病态或非有限，直接标 troubled / downgrade

后面真要优化，再把 `L` 的闭式写死。

---

## 第四步：在特征变量上做真正的 cell-wise PPM profile

这一段是 full PPM 的核心。

### 4.1 先做特征投影

对每个 cell `i`，用它自己的 `basis_i`，把邻域 primitive 投到 characteristic space：

\[
C_k = L_i W_k
\]

这里建议 **以 cell-centered basis 为主**，不要一上来就做 Roe/face-average basis。  
对你这个仓库，cell-centered 已经够做 full closeout，而且实现简单很多。

### 4.2 用受限斜率构造 face edge value

不要再直接只用现在那个 `7/12` face stencil。  
完整做法建议改成：

先对每个 characteristic component 做 TVD slope：

\[
\delta_i^{TVD} =
\begin{cases}
0, & (C_i - C_{i-1})(C_{i+1} - C_i) \le 0 \\
\mathrm{sign}(\cdot)\min\left(\frac{|(C_{i+1}-C_{i-1})|}{2}, 2|C_i-C_{i-1}|, 2|C_{i+1}-C_i|\right), & \text{otherwise}
\end{cases}
\]

然后 face interpolation 用这条：

\[
C_{i+1/2} =
\frac{1}{2}(C_i + C_{i+1})
-\frac{1}{6}(\delta_{i+1}^{TVD} - \delta_i^{TVD})
\]

这条比你现在的裸 `7/12` 更稳，而且本质上就是 PPM 常用的受限形式。

### 4.3 组 cell profile

对 cell `i`：

- `C_L,i = C_{i-1/2}`
- `C_R,i = C_{i+1/2}`
- `ΔC_i = C_R,i - C_L,i`
- `C6_i = 6*C_i - 3*(C_L,i + C_R,i)`

### 4.4 做 monotone reset

这一段不能省。

先做第一道：

如果

\[
(C_R - C_i)(C_i - C_L) \le 0
\]

那就说明这个 cell 在这个特征分量上不该重构，直接设：

- `C_L = C_i`
- `C_R = C_i`

然后再做 overshoot reset。  
你可以直接用 thesis 里那种不对称 edge reset 规则：

- 如果某一边离中心过远、另一边过近，就把远的一边按 `2x` 规则拉回。

这一步做完以后，再把 `c_left/c_right` 变回 primitive edge state，检查：

- `rho > 0`
- `p > 0`
- `chi_e >= 0`
- finite

只要任一 edge 非物理，就把**整个 cell** 标成 troubled，不要继续“修一修”。

---

## 第五步：做真正的 characteristic tracing

这是你现在完全没有的东西。

对每个 characteristic component `m`，定义：

\[
\sigma_m = \frac{|\lambda_m| \Delta t}{\Delta x_{\text{eff}}}
\]

建议：

- 如果 `sigma > 1 + eps`，直接记 troubled
- 算术上再 clamp 到 `[0,1]`

### 5.1 从 cell profile 积分到界面

对于 cell `i` 的右界面 `i+1/2`，右行波用：

\[
I^+_{m,i}
=
C_{R,i}
-\frac{\sigma_m}{2}
\left[
\Delta C_i
-
\left(1-\frac{2}{3}\sigma_m\right) C6_i
\right]
\]

对于 cell `i` 的左界面 `i-1/2`，左行波用：

\[
I^-_{m,i}
=
C_{L,i}
+\frac{\sigma_m}{2}
\left[
\Delta C_i
+
\left(1-\frac{2}{3}\sigma_m\right) C6_i
\right]
\]

### 5.2 face 左右态怎么组

对界面 `i+1/2`：

- 左态来自左 cell `i`
- 右态来自右 cell `i+1`

我建议你第一版先用这个最清楚的装配法：

对于左态 `W_{i+1/2}^L`：

- 若 `lambda_m >= 0`，用 `I^+_{m,i}`
- 若 `lambda_m < 0`，用 `C_{R,i}[m]`

对于右态 `W_{i+1/2}^R`：

- 若 `lambda_m <= 0`，用 `I^-_{m,i+1}`
- 若 `lambda_m > 0`，用 `C_{L,i+1}[m]`

然后分别乘各自 cell 的 `R_i` / `R_{i+1}` 回到 primitive space。

这版实现没有那么“教科书炫技”，但已经是**真正的 characteristic traced-interface PPM** 了，而且非常适合你现在的代码架构。

### 5.3 fallback 规则

界面 `i+1/2` 只要出现下面任一项，就直接退回 first-order：

- 左 cell 或右 cell 被标 troubled
- traced left/right state 非物理
- traced left/right state 非有限
- 后续 HLLC solve 失败

退化语义就保持你现在这套：

- `left = primitive_i`
- `right = primitive_{i+1}`
- `downgraded_to_first_order = true`

不要 clip 后继续走高阶。

---

## 第六步：把 `static_grid_hydro.cpp` 的 PPM 接口改成“传 dt 和 line metric”

这个改动只要做对一次，后面 MPI 基本不用返工。

### 6.1 改 `BuildPpmInterfaceFluxes(...)`

它现在缺的就是 `dt` 和 `dx_eff`。  
建议改成：

```cpp
bool BuildPpmInterfaceFluxes(
    const std::vector<HydroConservativeState>& ghosted_line,
    SweepDirection direction,
    std::size_t interior_cells,
    std::size_t ghost_layers,
    double dt_s,
    const std::vector<double>& effective_cell_widths,
    std::vector<HydroConservativeState>& interface_fluxes,
    std::size_t& downgraded_interface_count,
    std::string& failure_reason)
```

### 6.2 在 `static_grid_hydro.cpp` 里新增 line metric helper

加一个跟 `ExtractGhostedLine(...)` 并列的 helper，专门抽取这条 line 的 `dx_eff`。

对你现在的 mesh，这一步其实很好做，因为当前几何是规则的：

- **radial line**  
  \[
  \Delta x_{\text{eff}} = \Delta r
  \]

- **theta line**（固定某个 radial shell）  
  \[
  \Delta x_{\text{eff}} = r_c \Delta \theta
  \]

- **phi line**（固定 radial 和 theta）  
  \[
  \Delta x_{\text{eff}} = r_c \sin(\theta_c)\Delta \phi
  \]

也就是说，在你当前 mesh builder 下，**一条 line 上的 `dx_eff` 实际上是常数**。  
这会把实现难度降很多。

### 6.3 这一步很关键

不要为了这次 full PPM 去改 ghost substrate。  
你现在的方向 ghost 准备已经够用：

- radial: override / prepared ghosts
- theta: pole remap
- phi: periodic

full PPM 这次只需要**消费它们**，不用重写它们。

---

## 第七步：保留 HLLC，不要动 Riemann solver

你这次不需要改 `hllc_solver.*`。

架构上最稳的是：

- `ppm_reconstruction`：负责 characteristic reconstruction + tracing
- `hllc_solver`：负责 Riemann solve
- `static_grid_hydro`：负责 flux divergence
- `HydroOperator`：负责 runtime execution / writeback

也就是说，PPM 做出来以后，`BuildPpmInterfaceFluxes(...)` 仍然是：

1. `ReconstructPpmLine(...)` 得到每个 face 的 `left/right`
2. `MakeDirectionalConservativeState(...)`
3. `SolveHllcRiemann(...)`
4. `FromDirectionalFlux(...)`

这条责任边界别打散。

---

## 第八步：测试怎么补，才算真收住

这一段很重要。  
你现在已经有两条测试，但它们还不够证明“完整 PPM”。

### 8.1 先把 `test_ppm_reconstruction.cpp` 升级成 full-PPM contract

至少加 3 个断言：

#### A. traced state 必须依赖 `dt`
同一条平滑单调 line，分别用 `dt_small` 和 `dt_large` 调 reconstruction：

- 至少一个界面的 `left` 或 `right` 要发生变化

这个断言非常值钱，因为它直接排除了“没有 tracing 的假 PPM”。

#### B. smooth case 下 `left != right`
对一个平滑梯度 line，非 troubled 界面上：

- `interface.left` 和 `interface.right` 不应完全相同

这条会直接打掉你当前实现，因为你现在成功分支就是 `{face_state, face_state, false}`。

#### C. troubled case 继续保留 downgrade
保留你现在那种振荡/低压坏 stencil 测试，确认：

- `downgraded_interface_count > 0`
- downgrade 后左右态仍 physical

### 8.2 加一条 basis/unit test

建议新增一个 `test_ppm_characteristic_basis.cpp`：

- 检查 `L * R ≈ I`
- 检查 `A * R ≈ R * Λ`

这里 `A` 你可以直接按 primitive Jacobian 写出来。  
这条测试会让你的 characteristic path 非常硬。

### 8.3 升级 `test_hydro_operator_ppm_reconstruction.cpp`

除了现在已有的 diagnostics，再加两条：

- `result.execution_evidence->implementation_id == "p1.hydro.operator.ppm_static_grid"`
- diagnostics 里不要再只剩“first-order static-grid HLLC update”那种误导性语义

### 8.4 再补一条 single-rank 数值验证

在做 MPI 之前，我建议至少补一条更像“closeout”的单 rank case：

#### 方案一：径向 shock tube / radial wave with PPM
复用 `test_hydro_spherical_radial_wave_sanity.cpp` 的 scaffold，但打开：

- `use_ppm_reconstruction = true`
- 最好先只开 radial sweep

检查：

- 路径执行成功
- budget 仍然完整
- 相比 first-order，接触/跃迁更陡，或者至少 downgrade 比例可审计

#### 方案二：smooth line convergence
更干净一点的做法，是给 `ReconstructPpmLine(...)` 做一个解析平滑函数线测试：

- 用光滑函数构造 line
- 检查 face traced state 对解析解的误差随分辨率下降

这条更能证明“高阶”，而不只是“能跑”。

---

## 第九步：顺手把 observability 修正掉

这两个地方现在就该一起改：

### 9.1 `HydroOperator` 成功诊断语义
现在 PPM 成功路径里还会追加：

- `p1.hydro.stage.hllc_static_grid`
- `"first-order static-grid HLLC update"`

PPM 做满之后，这句会明显不对。  
建议改成类似：

- `p1.hydro.stage.hllc_static_grid`
- `"hydro stage executed static-grid HLLC flux assembly"`
- 如果 PPM 打开，再额外加  
  `"hydro stage executed direction-local characteristic PPM traced-interface reconstruction before HLLC flux assembly"`

### 9.2 失败时的 implementation_id
`FailedHydroStageResult(...)` 现在是硬编码 first-order implementation id。  
你应该让它吃一个可选参数；PPM 失败时也要能看出是 PPM 失败。

这个在你后面做 MPI 调试时会特别值钱。

---

# 我建议的施工顺序

按这个顺序最稳：

### commit 1
只改 API 和 tests，让 full-PPM failing-first 先立起来

- `ReconstructPpmLine(... dt_s, dx_eff ...)`
- `test_ppm_reconstruction.cpp` 新增：
  - dt-sensitive
  - left/right not identical
- 先让它红

### commit 2
加 characteristic basis + 6x6 inverse + unit test

- `BuildCharacteristicBasis(...)`
- `Invert6x6(...)`
- `test_ppm_characteristic_basis.cpp`

### commit 3
加 cell profile build

- `PpmCellProfile`
- slope limiter
- interface edge values
- monotone reset
- primitive positivity checks

### commit 4
加 traced-interface assembly

- `I^+ / I^-`
- `W_face^L / W_face^R`
- interface-level fallback

### commit 5
接入 `static_grid_hydro.cpp`

- `ExtractLineEffectiveWidths(...)`
- `BuildPpmInterfaceFluxes(... dt, dx_eff ...)`
- diagnostics 升级

### commit 6
补 single-rank closeout tests

- `test_hydro_operator_ppm_reconstruction.cpp`
- `test_hydro_spherical_radial_wave_sanity.cpp` 的 PPM 版本
- 最好再来一条 theta/phi PPM 执行 case

### commit 7
更新 `assumption-ledger.md`

把现在的 A-P1-008 从：

- temporary narrow primitive PPM slice

改成：

- single-rank characteristic/traced-interface PPM 已闭合
- remaining gap is real MPI ng=3 high-order closeout

---

# 这次先不要做的事

为了不把 full PPM 做成无底洞，这几件事先别绑在这一轮里：

1. **不要先做 MPI**
   先把 single-rank traced-interface 路径收干净。

2. **不要先做 contact steepening / flattening**
   这些是 tuning，不是 first closeout blocker。

3. **不要把 geometric source 强行塞进 characteristic tracing**
   你当前 hydro path 已经接受 old-time explicit source step，这一轮先把 homogeneous directional PPM 做实。

4. **不要动 HLLC solver**
   这轮只换 interface states，不换 Riemann solver。

---

# 你什么时候可以说“PPM 已完成，MPI 还没做”

满足下面 5 条，就可以这么说：

1. `ReconstructPpmLine(...)` 显式依赖 `dt/dx`
2. smooth 非 troubled case 下，face `left/right` 明显不相同
3. reconstruction 在 characteristic space 上完成，不再是 primitive face-value shell
4. troubled interface 会退回 first-order HLLC，而不是 clip
5. `HydroOperator` 在 single-rank 的 `r/theta/phi` case 上都能真实执行这条路径，并有清楚 diagnostics

做到这一步，再去做 real MPI，就很顺了，因为那时 MPI 只是在运这条已经闭合的高阶路径，而不是一边修数值核一边修并行外壳。

