给你一版**能直接落到你当前 DEC3D 仓库里的 `Sedov / Noh` case 草案**，目标是回答两件事：

1. **如何证明你正确实现了球坐标系下的 hydro 方程**
2. **如何证明你能正确解析激波 / 冲击波**

先把验证逻辑说清楚：

- `Case 2/3/4` 这类你现在已有的 stronger checks，主要证明的是  
  **实现路径稳定、true-3D 路径存在、MPI/seam/macro/ALE 不会马上把解搞坏**。  
  这在你当前 `P1` 文档里也已经是正式证据口径的一部分。fileciteturn57file0turn57file1
- 但它们**不能替代**标准 PDE / shock 验证。
- 真要证明“球坐标 hydro 方程”和“激波解析能力”对，你需要加：
  - **Sod**：验证 HLLC / PPM 对 shock / contact / rarefaction 的解析能力
  - **Sedov**：验证球坐标下强球对称激波
  - **Noh**：验证收缩流、反射 shock、强压缩激波

所以，Sedov / Noh 不是锦上添花，而是你当前 `P1` 从“实现路径验证”走向“物理 / 方程验证”的关键一步。

---

# 一、什么时候开始做这两个 case

按你当前代码状态，**已经有资格开始做**，但建议第一轮先用最干净的配置：

- `single-rank`
- `ALE = off`
- `macro_zoning = off`

原因很简单：  
先验证 **PDE 本体** 和 **shock capture 本体**；等这两条站住，再把：

- `ALE on/off`
- `macro on/off`
- `single-rank / multi-rank MPI`

一层层叠回去。

这和你当前 `P1` 文档的阶段语义是一致的：`P1` 里 hydro、PPM、球坐标源项、ALE、macro-zoning、MPI 都是正式能力，但数值 closeout 不该一口气全叠在一个新 case 上。fileciteturn57file12turn57file5

---

# 二、Sedov spherical blast case 草案

## 1. 用途
Sedov–Taylor 爆轰最适合证明：

- 球坐标几何项写对了
- 强球对称激波传播对了
- 你的 3D backend 没把一个球对称问题搞坏

FLASH 文档明确把 Sedov 列成 hydrodynamics test problem，用来检查**强激波和非平面几何**。fileciteturn57file3

---

## 2. 几何与运行模式

### 第一轮基线
- `single-rank`
- `ALE = off`
- `macro_zoning = off`

### 第二轮实现路径回归
- `ALE = on`
- `macro_zoning = off`

### 第三轮实现路径回归
- `ALE = on`
- `macro_zoning = on`
- 再加 `multi-rank MPI`

---

## 3. 建议的 deck 草案

建议文件名：

- `cases/hydro_only/p1-sedov-spherical.json`

草案如下：

```json id="7qkq9j"
{
  "case": {
    "type": "hydro_only",
    "name": "p1_sedov_spherical",
    "case_contract_version": "p1-sedov-v1"
  },
  "mesh": {
    "nr": 256,
    "ntheta": 16,
    "nphi": 32,
    "r_min": 0.0,
    "r_max": 1.0,
    "nghost": 3
  },
  "runtime": {
    "t_end_s": 0.05,
    "max_steps": 100000,
    "cfl": 0.5,
    "dt_min_s": 1.0e-12,
    "dt_max_s": 1.0e-2
  },
  "hydro": {
    "gamma": 1.4,
    "ppm_enabled": true,
    "ale_enabled": false,
    "macro_zoning_enabled": false
  },
  "mpi": {
    "radial_decomposition": true
  },
  "init": {
    "kind": "sedov_spherical",
    "rho0": 1.0,
    "p0": 1.0e-5,
    "v_r0": 0.0,
    "v_theta0": 0.0,
    "v_phi0": 0.0,
    "blast_energy": 1.0,
    "blast_radius": 0.01,
    "electron_energy_fraction": 0.5
  },
  "output": {
    "write_profiles": true,
    "write_shell_maps": true,
    "write_shock_tracking": true,
    "write_budget_residual": true
  }
}
```

---

## 4. 初值说明

推荐：

\[
\rho_0 = 1,\qquad p_0 = 10^{-5},\qquad \gamma = 1.4
\]

速度全零。  
在中心半径 \(\delta r\) 内沉积总能量 \(E_{blast}\)。

为了和你当前双温状态模型一致，初始化时还要选一个电子能量比例，例如：
\[
E_{electron} = f_e \, E_{internal},\qquad f_e = 0.5
\]
这只是 hydro-only 初始化约定，记入 assumption ledger 即可。

---

## 5. 解析 / 半解析对比量

Sedov 不是简单一条闭式 profile，但至少有非常明确的 shock radius 标度：

\[
R_s(t)=C_3(\gamma)\left(\frac{E\,t^2}{\rho_0}\right)^{1/5}.
\]

这条足够用来做**第一层强激波验证**。  
FLASH 文档也明确用 Sedov 数值结果和解析 blast 解对比 shock 结构与 radial bins。fileciteturn57file3

建议你至少对比：

### A. shock radius
从数值解里提取 shock front 位置 \(R_s^{num}(t)\)，对比参考 \(R_s^{ref}(t)\)

### B. jump quantities
shock 后方刚过 front 的：
- `rho`
- `v_r`
- `p`

### C. radial-binned profiles
把 3D spherical 解做 radial bin，得到：
- `rho(r,t)`
- `v_r(r,t)`
- `p(r,t)`
- `Te(r,t)`

然后和 Sedov 参考 profile 对照  
（如果你现在还没接完整解析 profile evaluator，先做 shock radius + jump + radial bin 也够第一轮）

---

## 6. 通过阈值初稿

### baseline（single-rank, ALE off, macro off）
建议先定：

- shock radius 相对误差  
  \[
  \le 5\%
  \]
- shock 后 `rho` plateau / jump 相对误差  
  \[
  \le 10\%
  \]
- shell 上非球对称 mode power ratio  
  \[
  \le 10^{-6}
  \]
- `m_theta`、`m_phi` 相对于 `m_r` 的比例  
  \[
  \frac{L1(|m_\theta|+|m_\phi|)}{L1(|m_r|)} \le 10^{-8}
  \]

### stronger regression
分辨率翻倍后要求：
- `R_s(t)` 误差下降
- shock 厚度以 cell 数计不恶化

---

## 7. 必须产出的工件

建议目录：

- `analysis/output/case_sedov_spherical`

至少要有：

- `summary.txt`
- `case_manifest.txt`
- `shock_radius_vs_time.txt`
- `radial_profile_t*.txt`
- `budget_residual.txt`
- `shell_map_rho_t*.txt`
- `shell_map_te_t*.txt`
- `shell_map_velocity_t*.txt`
- 对应 png

---

# 三、Noh spherical implosion case 草案

## 1. 用途
Noh 更适合证明：

- 汇聚流和 outward shock 处理是否对
- 强压缩 shock capture 是否稳
- 动能 → 内能转换是否合理
- 原点附近是否有 wall-heating / density dip 等问题

这是 shock/implosion 验证里非常经典的标准题。

---

## 2. 运行模式

第一轮仍建议：

- `single-rank`
- `ALE = off`
- `macro_zoning = off`

然后再做：
- `ALE on/off`
- `macro on/off`
- `multi-rank MPI`

---

## 3. 建议的 deck 草案

建议文件名：

- `cases/hydro_only/p1-noh-spherical.json`

草案如下：

```json id="ydmw9a"
{
  "case": {
    "type": "hydro_only",
    "name": "p1_noh_spherical",
    "case_contract_version": "p1-noh-v1"
  },
  "mesh": {
    "nr": 256,
    "ntheta": 16,
    "nphi": 32,
    "r_min": 0.0,
    "r_max": 1.0,
    "nghost": 3
  },
  "runtime": {
    "t_end_s": 0.6,
    "max_steps": 100000,
    "cfl": 0.5,
    "dt_min_s": 1.0e-12,
    "dt_max_s": 1.0e-2
  },
  "hydro": {
    "gamma": 1.6666666667,
    "ppm_enabled": true,
    "ale_enabled": false,
    "macro_zoning_enabled": false
  },
  "mpi": {
    "radial_decomposition": true
  },
  "init": {
    "kind": "noh_spherical",
    "rho0": 1.0,
    "p0": 1.0e-12,
    "u0": -1.0,
    "v_theta0": 0.0,
    "v_phi0": 0.0,
    "electron_energy_fraction": 0.5
  },
  "output": {
    "write_profiles": true,
    "write_shell_maps": true,
    "write_shock_tracking": true,
    "write_budget_residual": true,
    "write_wall_heating_metrics": true
  }
}
```

注意：理论上经典 Noh 常写 `p0 = 0`，但代码里通常用极小正值替代，避免完全零压导致某些 EOS/恢复链数值病。这个要在 assumption ledger 里记一笔。

---

## 4. 解析解（球对称、\(\gamma=5/3\)、\(\rho_0=1\)、\(u_0=-1\)）

经典 spherical Noh 下：

### shock 位置
\[
r_s(t)=\frac{t}{3}
\]

### shocked region 密度
\[
\rho_2 = 64
\]

### shocked region 速度
\[
u_2 = 0
\]

这三条是最重要、最稳的比较量。  
压力 plateau 你也可以做，但第一轮先抓这三条已经足够有说服力。

---

## 5. 你应该对比什么

### A. shock front position
从数值解提取 shock 位置 \(r_s^{num}(t)\)，对比：
\[
r_s^{ref}(t)=t/3
\]

### B. plateau density
shock 后方密度 plateau 应接近：
\[
\rho_2 = 64
\]

### C. post-shock velocity
shock 后方应接近静止：
\[
u_2 = 0
\]

### D. origin artifact
看原点附近是否出现明显：
- wall heating
- density dip
- 非物理尖峰

---

## 6. 通过阈值初稿

### baseline（single-rank, ALE off, macro off）
建议先定：

- shock position 相对误差  
  \[
  \le 5\%
  \]
- shocked density plateau 相对误差  
  \[
  \le 10\%
  \]
- shocked velocity plateau  
  \[
  |u_2^{num}| \le 0.05
  \]
- `m_theta`、`m_phi` 相对于 `m_r` 的 leakage  
  \[
  \frac{L1(|m_\theta|+|m_\phi|)}{L1(|m_r|)} \le 10^{-8}
  \]

### stronger regression
- 分辨率翻倍后，shock position 和 plateau 误差下降
- origin artifact 不恶化

---

## 7. 必须产出的工件

建议目录：

- `analysis/output/case_noh_spherical`

至少要有：

- `summary.txt`
- `case_manifest.txt`
- `shock_position_vs_time.txt`
- `radial_profile_t*.txt`
- `budget_residual.txt`
- `wall_heating_diagnostics.txt`
- `shell_map_rho_t*.txt`
- `shell_map_velocity_t*.txt`
- 对应 png

---

# 四、这两个 case 如何接进你当前 `P1` 验收链

我建议它们不要替代你现在已有的 `Case 2/3/4`，而是补在前面：

## baseline analytic / self-similar verification
1. `Sod`
2. `Sedov spherical`
3. `Noh spherical`

## implementation-path stronger checks
4. `case2`
5. `case3`
6. `case4`

也就是说：

- **Sedov / Noh** 负责证明  
  “球坐标 hydro 方程本体和 shock capture 对”
- **case2/3/4** 负责证明  
  “你的 true-3D + MPI + seam + macro/ALE 路径稳”

这两类验证不要混为一谈。

---

# 五、我建议你现在先怎么落地

最实用的顺序是：

## 第一阶段
- `p1-sedov-spherical.json`
- `p1-noh-spherical.json`
- `single-rank`
- `ALE off`
- `macro off`

先把：
- shock position
- jump/plateau
- radial profile
- shell 对称性

做出来。

## 第二阶段
在这两个 case 上加：
- `ALE on/off`
- `macro on/off`

## 第三阶段
在这两个 case 上加：
- `single-rank / multi-rank MPI parity`

> **Sedov 用来证明：你在球坐标下能正确传播强球对称激波；**
> **Noh 用来证明：你能正确处理汇聚流和 outward shock。**
>
> **两者第一轮都应该用：single-rank + ALE off + macro off。**
>
> **通过后，再把 ALE / macro / MPI 一层层叠回去。**