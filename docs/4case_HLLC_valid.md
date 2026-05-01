可以。基于你**当前这版 ghost substrate + 一阶 full directional static-grid HLLC**，我建议在进入 `PPM` 之前，先做下面 **4 个数值验收 case**。这 4 个 case 的目标不是“把 `P1` 做完”，而是先证明：

- 当前低阶 3D 球坐标 hydro 核**方程写全了**
- `theta/phi` 路径不是假执行
- 极点/周期边界不会偷偷污染解
- 这套地基值得继续往 `PPM` 叠

这也符合你当前 `P1` 文档的边界：`P1` 要的是 thesis-level explicit hydro core，而 `P1a` 不能被误判成整个 `P1`；同时 `P1` 结束时必须有 true-3D hydro-only case 和 `rho / hydro-only Te / velocity` 的时间序列输出。`PPM` 也已经被你写死为需要三层 hydro ghost 的真高阶路径。fileciteturn19file0turn16file9turn16file2

---

## Case 1：Uniform-state invariance（静止均匀态不变性）

### 目的
这是你现在最重要的 blocker case。  
它用来检查：

- 几何源项和通量离散是否在离散层面基本平衡
- `theta/pole` 和 `phi/periodic` 边界是否会凭空制造结构
- 你的事务性 hydro commit 是否真的不污染真源

### 初值
全域统一：

- `rho = rho0`
- `mom_r = mom_theta = mom_phi = 0`
- `E_fluid_total = E0`
- `E_electron = Ee0`
- 因而 `chi_e` 也是常量恢复出的 operator-local 常量

建议做 **两个子版本**：

#### 1A. source-off
只跑三方向 HLLC sweep，不加几何源项。  
这条应该是最容易做到**机器精度级不变**的。

#### 1B. source-on
打开当前球坐标几何源项步。  
这条如果还漂，就说明你现在的 source 离散和通量离散还没真正配平。

### 必看指标
逐 cell 比较一步后的：

- `rho`
- `mom_r`
- `mom_theta`
- `mom_phi`
- `E_fluid_total`
- `chi_e`

是否与初始相同。  
这里不要只看“壳层内角向一致”，而要看**是否真的没漂**。

### 通过标准建议
- `source-off`：必须近似机器精度不变  
- `source-on`：如果不过，不要先调宽阈值，先把它当 **blocker** 来查  
  因为这条正好直接暴露你现在的离散平衡问题

### 建议输出
- `rho` slice
- `hydro-only Te` slice
- `|v|` slice  
虽然理论上应该全平，但这三张图非常适合肉眼发现“哪里悄悄长出结构了”。

---

## Case 2：3D backend 上的 spherical shock / wave case（球对称波传播）

### 目的
这条不是为了证明 3D，而是为了证明：

- 你现在的核已经能在 3D backend 上做**可信的径向波传播**
- `theta/phi` 路径不会把一个球对称问题搞坏
- 径向 boundary contract 不再是会污染最基本数值结果的 placeholder

### 初值
做一个**角向完全均匀**、只在 `r` 上有不连续的 case，例如：

- radial shock tube
- radial pressure pulse
- 简化球对称 blast / acoustic pulse

核心要求是：

- 所有 `(theta,phi)` 上初值完全相同
- 仍然走 full directional 3D hydro path，而不是旁路 1D solver

### 必看指标
- 径向剖面的 `rho / v_r / P / hydro-only Te`
- `mom_theta`、`mom_phi` 是否仍接近 0
- 不同 `(theta,phi)` 壳线之间是否几乎重合
- 单 rank / 多 rank MPI 的 parity

### 通过标准建议
- 不要求你一上来就和精确解逐点严丝合缝
- 但至少要满足：
  - 波前位置合理
  - 没有明显伪振荡
  - 不凭空长出横向动量
  - 多 rank 不在径向 seam 处长出额外伪波

### 建议输出
- 若干时刻的径向 profile
- 同时输出一张 2D slice，让你确认球对称初值没有被 3D 路径搞坏

---

## Case 3：True-3D low-mode static-grid sanity case（真 3D 低模静态网格检验）

### 目的
这是你进入 `PPM` 前最重要的 **true 3D case**。  
它用来证明：

- 当前 full directional 路径已经不是“3D 数组上的 radial-only 更新”
- `theta/phi` 通量路径真的在起作用
- 你后面 `P1` 要交付的 true-3D plot artifact 是有根的，不是最后临时拼出来的

### 初值
建议用一个**温和**的角向结构 case，不要一上来就太复杂。  
例如：

- 在某个壳层附近给 `v_r` 一个低模扰动  
  \[
  \delta v_r \propto f(r) Y_2^1(\theta,\phi)
  \]
- 或者在 `P` 上给一个低模角向结构
- 但先保持：
  - static grid
  - no ALE
  - no macro-zoning
  - no PPM

幅值建议先小一点，像 1%–3% 这种量级，目的是看方向路径和边界 contract，不是马上测强非线性。

### 必看指标
- `rho(t)` 的 3D slice / shell map
- `hydro-only Te(t)` 的 3D slice / shell map
- `|v|(t)` 或 velocity vector field
- mode projection：
  - 幅值 `A_lm(t)`
  - 相位 `phi_lm(t)`

### 必看对照
要和 **radial-only 路径**做一次对照。  
如果 radial-only 和 full directional 的结果几乎一样，那就说明 `theta/phi` 路径虽然执行了，但贡献还不够可信。

### 通过标准建议
- full directional 与 radial-only **必须明显不同**
- 不允许出现：
  - 明显的 pole 裂缝
  - `phi` seam 裂缝
  - 无由来的 NaN / 爆大动量
- 至少能稳定输出一小段时间序列图，而不是只跑一步

---

## Case 4：Pole-adjacent seam stress case（极区与缝合面压力测试）

### 目的
这是专门给你现在刚补上的 ghost substrate 准备的。  
它的目标不是漂亮物理解，而是：

- 把 `theta` 极点 remap
- `phi` 周期 ghost
- full directional 3D sweep

一起压到极区附近，看它们会不会互相打架。

### 初值
建议造一个**靠近极区集中**的角向结构，例如：

- 在 `theta ≈ 0` 和 `theta ≈ π` 附近给一个局部扰动
- 同时 `phi` 上要有非平凡结构（不能是 `m=0`）
- `r` 上尽量简单，避免把问题和复杂径向结构混在一起

一个很实用的形式是：

- 角向包络集中在 polar cap
- `phi` 方向再乘上 `cos(m\phi)` 或 `sin(m\phi)`

这样就能同时压到：
- pole remap
- `phi` periodic seam
- `theta/phi` directional flux

### 必看指标
- 极区附近是否出现 NaN
- `mom_theta / mom_phi` 是否异常放大
- polar cap 两侧的 shell map 是否连续
- `phi` seam 是否出现裂缝
- 单 rank / 多 rank MPI 下极区和径向 rank 边界是否额外长伪波

### 通过标准建议
这条不要求物理上多漂亮，但至少要满足：
- 不炸
- 不裂
- 不在极区附近长出明显非物理尖峰
- 单 rank / 多 rank 的 seam 行为一致

### 建议输出
- 北极 / 南极附近的局部切片
- shell map
- `|v|` 图  
这条图像往往比标量更有用。

---

# 这 4 个 case 的执行顺序

我建议严格按这个顺序跑：

1. **Uniform-state invariance**  
2. **Spherical shock / wave case**  
3. **True-3D low-mode sanity case**  
4. **Pole-adjacent seam stress case**

原因很简单：

- 1 和 2 是“底层核有没有写对”
- 3 和 4 才是“3D 路径和边界系统有没有真站住”

如果 1 或 2 还不稳，就别急着解释 3 和 4 的现象。

---

# 每个 case 都建议统一产出的工件

为了后面进入 `PPM` 前能真正做 closeout，我建议这 4 个 case 都产出同一套最小工件：

- `rho` 时间序列图
- `hydro-only Te` 时间序列图
- `velocity magnitude` 或 `velocity vector` 时间序列图
- `field-update summary`
- `execution evidence`
- `numerical check summary`

其中 true-3D 的 case（Case 3 和 4）再额外加：

- shell map
- mode projection（至少 Case 3）
- seam / polar diagnostics（至少 Case 4）

---

# 一句话总结

如果你只想记一句最实用的话：

> **进入 `PPM` 之前，先用这 4 个 case 把当前的一阶 full-directional spherical hydro 核“数值上站住”：**
>
> 1. 静止均匀态真不变  
> 2. 球对称波传播不被 3D 路径搞坏  
> 3. 真 3D 低模 case 能明显区别于 radial-only  
> 4. 极区和周期缝合面在 3D 路径下不炸不裂
