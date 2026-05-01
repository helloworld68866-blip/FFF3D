## 把 `γ` 做成**运行时配置参数**，不是编译期常量

## 现在代码里最直接的问题是：

## 它在当前源码里被多处直接调用，包括：

- `src/state/hydro_state/hydro_view.hpp`
- `src/state/hydro_state/hydro_view.cpp`
- `src/hydro/riemann/hllc_solver.cpp`
- `src/hydro/reconstruction/ppm_reconstruction.cpp`
- `src/hydro/driver/hydro_operator.cpp`

所以如果你想让 `γ` 可配置，**不要再让这些地方去调用一个全局 constexpr**，而要改成：

> **由 runtime / case config 提供一个 `gamma`，在 `bind(...)` 时注入 hydro operator，再显式传给所有数值核。**

这和你当前系统设计的方向是一致的：operator 不该自己发明物理常量，runtime/context 才是统一配置入口。`P1` 的 hydro operator 也已经被设计成 runtime-facing contract。fileciteturn59file10turn59file0

---

# 我建议的最小改造方案

下面按“最小改动集”来列。

---

## 1. 增加一个正式的 `HydroEquationOfState` / `HydroParameters`
最小版可以非常简单：

```cpp id="skkcd7"
struct HydroParameters {
    double gamma = 5.0 / 3.0;
};
```

更稳一点可以叫：

```cpp id="2xcsi3"
struct HydroEquationOfState {
    double gamma = 5.0 / 3.0;
};
```

### 放哪里
建议放到：

- `src/hydro/driver/`
或
- `src/state/hydro_state/`

如果你问我更推荐哪边：  
**放在 `hydro` 层，但由 runtime/config 注入。**

原因是：
- 当前 `P1` 只需要 hydro 用到 `γ`
- 以后 `P2-P4` 如果也要用，再把它提升成更通用的 EOS config object

---

## 2. 把 `HydroIdealGasGamma()` 去掉，改成显式传参
当前这些函数应该改：

### `hydro_view.hpp/.cpp`
把：

- `ElectronPressureFromElectronEnergyDensity(...)`
- `ElectronEnergyDensityFromPressure(...)`

改成：

```cpp id="rjfo5u"
double ElectronPressureFromElectronEnergyDensity(
    double electron_energy_density,
    double gamma) noexcept;

double ElectronEnergyDensityFromPressure(
    double electron_pressure,
    double gamma) noexcept;
```

### 原因
因为现在这两条都依赖：
\[
P_e = (\gamma - 1) E_{electron}
\qquad\text{和}\qquad
E_{electron}=\frac{P_e}{\gamma-1}.
\]

只要 `γ` 可配置，这两条必须跟着变。

---

## 3. `RecoverPrimitiveState(...)` / `MakeConservativeState(...)` / `SoundSpeed(...)` 全部显式带 `gamma`
当前 `hllc_solver.cpp` 和 `hydro_operator.cpp` 里多处直接用：

\[
c_s = \sqrt{\gamma P / \rho}
\]

以及
\[
E = \frac{P}{\gamma-1} + \frac12 \rho v^2
\]

所以这些函数都要改成：

```cpp id="0wgt0e"
HydroPrimitiveState RecoverPrimitiveState(
    const HydroConservativeState& conservative,
    double gamma) noexcept;

HydroConservativeState MakeConservativeState(
    const HydroPrimitiveState& primitive,
    double gamma) noexcept;

double SoundSpeed(
    const HydroPrimitiveState& primitive,
    double gamma) noexcept;
```

### 受影响的位置
至少包括：

- `src/hydro/riemann/hllc_solver.cpp`
- `src/hydro/reconstruction/ppm_reconstruction.cpp`
- `src/hydro/driver/hydro_operator.cpp`

---

## 4. `HydroOperator` 在 `bind(...)` 时持有 `gamma`
最小改法：

```cpp id="ezisx8"
class HydroOperator {
public:
    bool bind(..., double gamma, ...);
private:
    double gamma_ = 5.0 / 3.0;
};
```

然后在：

- `estimate_dt()`
- `advance()`

里都用 `gamma_`，再往下传给：
- `RecoverPrimitiveState`
- `MakeConservativeState`
- `HLLC`
- `PPM`

### 为什么不要继续用全局函数
因为你以后：
- Sedov 用 `1.4`
- Noh 用 `5/3`
- 其它 hydro-only regression 可能还要不同 `γ`

这一定是 **case/runtime 级别** 的选择，不该是编译期常量。

---

## 5. 在 case/config 里加一个正式字段
你之前文档里其实已经有过类似配置草案。最简单就是在 case deck 里统一放：

```json id="52p1ho"
"hydro": {
  "gamma": 1.4
}
```

或者更通用一点：

```json id="oef8c1"
"physics": {
  "gamma": 1.4
}
```

### 我建议哪种
如果你只想先做 `P1`，放在 `hydro.gamma` 也行。  
如果你考虑后面 `P2-P4` 也统一复用，我更推荐：

```json id="5s69kk"
"physics": {
  "gamma": 1.6666666667
}
```

因为 Woo 论文当前总模型就是一个统一 ideal-gas `γ` 假设，不是 hydro 一家独占。fileciteturn59file2turn59file5

---

## 6. 你必须改的测试
这一步很重要，不然 Codex 很容易只做半套。

### A. 新增 `gamma` contract test
例如：

- `gamma = 5/3` 时，当前回归结果不变
- `gamma = 1.4` 时，`SoundSpeed`、`RecoverPrimitiveState`、`E ↔ P` 转换结果正确变化

### B. 更新 `HLLC` 测试
当前 `HLLC` 波速、星态、声速都依赖 `γ`。  
所以要加一条：

- 同一左右态，`γ=1.4` 和 `γ=5/3` 的波速/星态不同，且方向正确

### C. 更新 `estimate_dt` 测试
因为
\[
c_s = \sqrt{\gamma P/\rho}
\]
不同 `γ` 会直接影响 CFL。  
所以要加一条：

- 在同一几何、同一状态下，`γ` 增大时 `dt` 应变小

### D. 更新 `hydro_view` writeback 测试
因为
\[
P_e = (\gamma-1)E_e
\]
不同 `γ` 会改变 `chi_e \leftrightarrow P_e \leftrightarrow E_e` 这条链。

---

# 你现在最容易漏掉的地方

## 1. `chi_e` 路径
你现在的 hydro-local 电子通道语义是：

\[
E_e \rightarrow P_e \rightarrow \chi_e
\]
和
\[
\chi_e \rightarrow P_e \rightarrow E_e
\]

这里的：
\[
P_e = (\gamma-1) E_e
\]
只要 `γ` 改了，这条链就全部变。  
这一步最容易漏。

## 2. `PPM` 的 characteristic basis
当前 `ppm_reconstruction.cpp` 里构造特征基、声速都要用 `γ`。  
如果这里只改了一半，你会得到：
- primitive 恢复用一个 `γ`
- characteristic basis 用另一个（旧的 5/3）

这会直接把 PPM 搞坏。

## 3. `HydroBudgetResidualSummary`
如果它里面有任何从 `E` 推 `P` 或 `c_s` 的地方，也必须跟着带 `γ`。  
否则数值工件会和真实求解不一致。

---

# 我对“是否应该做成可配置”的建议

## 如果你的目标是：
### “尽量保持论文原模型”
那默认值必须还是：
\[
\gamma = 5/3
\]
而且文档里要明确：

> configurable gamma is a thesis-external engineering extension;  
> thesis-default remains `5/3`.

这条建议我很强。

## 如果你的目标是：
### “马上做 Sedov / Noh”
那也应该先把 `γ` 配置化，因为：
- Sedov 你大概率会想跑 `1.4`
- Noh 经典解又是 `5/3`

否则你的验证 case 会互相打架。

---

# 我建议你怎么记到文档里

在 `assumption-ledger.md` 里加一条类似：

- `assumption_id`: `p1-hydro-gamma-configurable`
- `phase`: `P1`
- `module`: `hydro`
- `thesis_gap`: thesis assumes ideal-gas gamma = 5/3; verification cases may require alternate gamma values
- `decision`: make gamma runtime-configurable while preserving 5/3 as the repository default
- `status`: `temporary`
- `impact_on_gates`: Sedov/Noh verification and hydro dt/HLLC/PPM paths must use case gamma consistently
- `impact_on_result_interpretation`: results with gamma != 5/3 are verification-mode runs, not direct thesis-default DEC3D reproductions

---

# 最后一行，最直接的实现清单

如果你问我“**要改什么**”，最小清单就是：

1. 去掉 `HydroIdealGasGamma()` 这个全局 constexpr  
2. 新增 `HydroParameters{ double gamma; }`  
3. `HydroOperator::bind(...)` 注入 `gamma`  
4. 所有依赖 `γ` 的函数改成显式传参：
   - `ElectronPressureFromElectronEnergyDensity`
   - `ElectronEnergyDensityFromPressure`
   - `RecoverPrimitiveState`
   - `MakeConservativeState`
   - `SoundSpeed`
   - `HLLC`
   - `PPM characteristic basis`
   - `estimate_dt`
5. case/config 增加 `gamma`
6. 加 `gamma` 相关 contract tests
7. assumption ledger 记一笔：默认仍是 thesis 的 `5/3`

