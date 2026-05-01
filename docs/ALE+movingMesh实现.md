先给最核心的两句话：

1. **网格当然要改。**径向 face 位置 \(r_{i+1/2}\) 必须从 \(t^n\) 更新到 \(t^{n+1}\)。  
2. **但不能只改网格。**论文的 moving-mesh 不是“先改几何，再拿旧的 cell-average 顶着用”，而是要对 **cell-integrated conserved quantities \((VQ)\)** 做更新：
\[
\frac{(VQ)^{n+1}-(VQ)^n}{\Delta t}
=
\big[A_r(F_r^*-\dot r\,Q^*)\big]_{i-1/2}^{i+1/2}
+\big[A_\theta F_\theta^*\big]_{j-1/2}^{j+1/2}
+\big[A_\phi F_\phi^*\big]_{k-1/2}^{k+1/2}
+\big(VS\big)^n.
\]
这就是 thesis 里 3D spherical conservative moving mesh 的正式更新式。mesh 只沿径向运动，moving boundary 对应的项是 \(A_r \dot r Q^*\)。fileciteturn44file0turn43file9

---

# 一、论文锁死的 ALE / moving mesh 骨架

Woo 论文明确锁死了 4 件事：

### 1. 只有径向 face 在动
论文明确说：
> In DEC2D and DEC3D, the mesh moves along the radial direction only. fileciteturn43file9turn44file0

所以：
- `theta` 面不动
- `phi` 面不动
- 只有 `r_{i+1/2}` 是时变的

### 2. moving-mesh 更新用 Leibniz rule 推导
论文先用：
\[
\int_D \frac{\partial Q}{\partial t}\,dV
=
\frac{d}{dt}\int_D Q\,dV
-
\int_{\partial D} Q\,v_{\partial D}\cdot dA
\]
然后把 moving boundary 的 contribution 写成：
\[
\int_{\partial D} Q\,v_{\partial D}\cdot dA
=
\big[A_r \dot r\,Q\big]_{i-1/2}^{i+1/2}.
\]
这一步非常关键：**moving mesh 本质上就是一个“网格边界输运项”**。fileciteturn44file0turn43file9

### 3. 3D spherical 的最终 moving-mesh 更新式就是 Eq. (5.112)
也就是上面那条：
\[
(VQ)^{n+1}
=
(VQ)^n
-\Delta t\,\Delta\big[A_r(F_r^*-\dot rQ^*)\big]
-\Delta t\,\Delta\big[A_\theta F_\theta^*\big]
-\Delta t\,\Delta\big[A_\phi F_\phi^*\big]
+\Delta t\,(VS)^n.
\]
这里：
- `F^*` 是 Riemann solver 的 upwind flux
- `Q^*` 是 Riemann solver 给出的 upwind interface state  
论文也是这么写的。fileciteturn43file9turn44file0

### 4. 外边界半径按壳层 CM 速度推进
你的设计文档已经把这一点写成正式 hydro design：
\[
\frac{dR}{dt}=\beta\,v_{CM}^{shell}(t^n)
\]
也就是最外层 radial face 的速度来自 shell center-of-mass velocity。这个是你现在项目里的 thesis-consistent 约束。fileciteturn43file4turn43file1

---

# 二、你现在代码里应该推进的变量是什么

按你已经锁死的项目规范，hydro 的局部包是：

\[
Q =
[\rho,\; mom_r,\; mom_\theta,\; mom_\phi,\; E_{fluid,total},\; \chi_e].
\]

这里：

- `chi_e` 是 operator-local
- `E_electron` 不直接当 hydro 守恒变量推进
- hydro 结束后再走  
  \[
  \chi_e \rightarrow P_e \rightarrow E_{electron}
  \]
  写回 canonical truth。fileciteturn43file4turn43file1

所以 ALE 更新时，真正要更新的是：

\[
U^{ext}_{ijk} \equiv V_{ijk} Q_{ijk}
\]

不是直接更新 cell-average `Q`。

---

# 三、`w_face` 怎么求

你问得最直接，这里给你最明确的口径。

## 1. 外边界 face 速度
最外层 face：

\[
w_{N_r+1/2} = \dot R^n = \beta\,v_{CM}^{shell}(t^n).
\]

这一步是你当前 thesis-consistent 设计里锁死的。fileciteturn43file4turn43file1

## 2. 内部 face 速度
论文片段没有给唯一的内部 predictor，所以这是**工程补全**，必须继续记 assumption ledger。  
当前最稳的做法是：同伦/比例缩放

如果：
- old outer radius = \(R^n\)
- new outer radius = \(R^{n+1}\)

那内部 face 直接按比例缩放：
\[
r_{i+1/2}^{\,n+1}
=
r_{i+1/2}^{\,n}\frac{R^{n+1}}{R^n}.
\]

于是：
\[
w_{i+1/2}
=
\frac{r_{i+1/2}^{\,n+1}-r_{i+1/2}^{\,n}}{\Delta t}.
\]

## 3. 统一定义
不管你内部 predictor 最后用哪一种，**代码里统一必须写成**：

\[
w_{i+1/2}
=
\frac{r_{i+1/2}^{\,n+1}-r_{i+1/2}^{\,n}}{\Delta t}.
\]

这样：
- proposal
- ALE flux
- runtime commit

三者语义一致。

---

# 四、`U_*`（更准确是 `Q^*`）怎么求

这是你当前 ALE 实现最容易做错的地方。

## 错误做法
不要做成：

- 先求 stationary-face 的 `F_HLLC(0)`
- 然后只按 `S_*` 的符号选左星态或右星态
- 再做 `F - w_face U_*`

这在 moving interface 下**不对**。

## 正确做法：按 `w_face` 在波扇中的位置选 upwind branch
你先照正常 HLLC 求：

- `S_L`
- `S_R`
- `S_*`
- `Q_L^*`
- `Q_R^*`
- `F_L^*`
- `F_R^*`

然后用 `w_face` 和波速比较，选真正的 moving-interface upwind 状态：

\[
Q_{upwind}(w_{face})=
\begin{cases}
Q_L, & w_{face}\le S_L,\\[4pt]
Q_L^*, & S_L < w_{face}\le S_*,\\[4pt]
Q_R^*, & S_* < w_{face}\le S_R,\\[4pt]
Q_R, & S_R < w_{face}.
\end{cases}
\]

对应的 upwind flux：

\[
F_{upwind}(w_{face})=
\begin{cases}
F_L, & w_{face}\le S_L,\\[4pt]
F_L^*, & S_L < w_{face}\le S_*,\\[4pt]
F_R^*, & S_* < w_{face}\le S_R,\\[4pt]
F_R, & S_R < w_{face}.
\end{cases}
\]

然后才有：

\[
F_r^{ALE}
=
F_{upwind}(w_{face})
-
w_{face}\,Q_{upwind}(w_{face}).
\]

### 这就是论文里的 `F^* - rdot Q^*`
论文的 `F^*` 和 `Q^*` 写成星号，是因为它们是 **Riemann solver 的 upwind solution**，不是“固定取左星态”或者“固定取 \(x/t=0\) 的界面值”。fileciteturn43file9turn44file0

---

# 五、真正 conservative 的 moving-mesh 更新要怎么做

这是另一个关键点：**不能只改 `Q`，必须改 `(VQ)`。**

## 正确更新量
定义：
\[
U^{ext,n}_{ijk}
=
V^n_{ijk} Q^n_{ijk}.
\]

然后更新：

\[
U^{ext,n+1}_{ijk}
=
U^{ext,n}_{ijk}
-\Delta t
\Big[
A_{r,i+1/2}^n F_{r,i+1/2}^{ALE}
-
A_{r,i-1/2}^n F_{r,i-1/2}^{ALE}
\Big]
\]
\[
-\Delta t
\Big[
A_{\theta,j+1/2}^n F_{\theta,j+1/2}^*
-
A_{\theta,j-1/2}^n F_{\theta,j-1/2}^*
\Big]
-\Delta t
\Big[
A_{\phi,k+1/2}^n F_{\phi,k+1/2}^*
-
A_{\phi,k-1/2}^n F_{\phi,k-1/2}^*
\Big]
+\Delta t\,V^n_{ijk} S^n_{ijk}.
\]

最后再恢复 cell-average：

\[
Q^{n+1}_{ijk}
=
\frac{U^{ext,n+1}_{ijk}}{V^{n+1}_{ijk}}.
\]

## 为什么不能只更新 `Q`
如果你只是：
- 用 old volume 做 flux divergence
- 直接更新 `Q`
- runtime 再 commit 新 mesh

那 `Q` 和新 `V^{n+1}` 就不一致，守恒性会坏。  
论文写 `((VQ)^{n+1}-(VQ)^n)/Δt` 不是装饰，是真正的保守结构。fileciteturn43file9turn44file0

---

# 六、和你当前项目最一致的实现顺序

你现在这套工程边界已经定得很清楚了，所以最稳的代码结构应该是：

## Step 1：old-time canonical state + old-time mesh
从 canonical truth 读：
- `rho`
- `mom_r`
- `mom_theta`
- `mom_phi`
- `E_fluid_total`
- `chi_e`

和 old mesh：
- `r_face^n`
- `V^n`
- `A_r^n`
- `A_theta^n`
- `A_phi^n`

## Step 2：用 old-time state 生成 radial ALE proposal
输出：
- `r_face^{n+1}`
- `w_face`
- `R_new`

hydro 自己只能 **propose**，不能 commit。

## Step 3：构造 extensive work state
\[
U^{ext,n}=V^n Q^n
\]

## Step 4：求 radial ALE flux
对每个 radial face：

1. reconstruction（当前一阶或 PPM）
2. HLLC：求 `S_L,S_R,S_*`
3. 求 `Q_L^*,Q_R^*,F_L^*,F_R^*`
4. 用 `w_face` 选 `Q_upwind, F_upwind`
5. 得到
   \[
   F_r^{ALE}=F_{upwind}-w_{face}Q_{upwind}
   \]

## Step 5：求 `theta/phi` 静态 flux
这两个方向仍然用普通 HLLC / PPM flux，不加 ALE。

## Step 6：对 `U^{ext}` 做 full update
\[
U^{ext,n+1}
=
U^{ext,n}
-\Delta t\,\Delta(A_r F_r^{ALE})
-\Delta t\,\Delta(A_\theta F_\theta)
-\Delta t\,\Delta(A_\phi F_\phi)
+\Delta t\,V^n S^n.
\]

source 仍然是你当前 phase-one 的 old-time explicit step，**不要**再额外把 moving-mesh 修正塞进 source。

## Step 7：除以 `V^{n+1}` 恢复 `Q^{n+1}`
\[
Q^{n+1} = U^{ext,n+1} / V^{n+1}
\]

## Step 8：authorized writeback
走你已经有的：
\[
\chi_e \rightarrow P_e \rightarrow E_{electron}
\]
再做 thermodynamic feasibility check，最后一次性 commit 到 canonical state。

## Step 9：runtime commit canonical mesh
只有 runtime 可以把：
- `r_face^n -> r_face^{n+1}`
- `R^n -> R^{n+1}`  
写进 canonical mesh。

---

# 七、你接下来必须补的两条 failing-first tests

如果你真要把 `Workstream 6b` 做成论文一致，我建议最少先补这两条红测。

## 测试 1：moving-interface branch 必须按 `w_face` 选，不是按 `0`
构造一个 radial Riemann 问题，让：

- `0` 落在左星区
- `w_face` 落在右星区

然后要求：
- `Q_upwind(w_face)` 必须和 stationary-face 选出来的不一样
- 当前如果还按 `S_*` 相对 0 去选，就应该先失败

## 测试 2：ALE update 必须对 `(VQ)` 保守
构造一个非零 mesh motion 的最小 case，检查：

\[
(VQ)^{n+1}
-
(VQ)^n
+
\Delta t\,\Delta\big(A(F-wQ)\big)
-
\Delta t\,VS
\approx 0
\]

至少先对：
- `rho`
- `mom_r`
- `E_fluid_total`

做。  
如果你现在还是只更新 `Q`，这条测试一定会暴露问题。

---

# 八、一句话总结

如果你只想记一句：

> **和论文一致的 ALE + moving mesh 做法，不是“提个 proposal，再把静止界面的 HLLC flux 减一个 `w_face U_*`”，而是：**
> 1. **按 `dR/dt = β v_CM^{shell}` 及其内部 predictor 生成 `r_face^{n+1}`，统一定义 `w_face=(r^{n+1}-r^n)/Δt`；**
> 2. **在 moving face 上按 `w_face` 在 HLLC 波扇中的位置选真正的 upwind `Q^*` 和 `F^*`；**
> 3. **用**
>    \[
>    (VQ)^{n+1}=(VQ)^n-\Delta t\,\Delta\big(A(F-\dot r Q)\big)+\Delta t\,(VS)^n
>    \]
>    **更新 extensive conserved quantities；**
> 4. **最后再除以 `V^{n+1}` 恢复新时刻 cell-average，并由 runtime commit 新网格。**

这就是论文锁死的 conservative moving mesh 结构；  
内部 face speed predictor、CFL prefactor、lower/origin boundary 这些仍然是你当前项目的工程补全，应该继续留在 assumption ledger 里。fileciteturn43file9turn44file0turn42file10

如果你愿意，我下一条可以直接按你当前代码结构，把  
**“Workstream 6b 需要改哪些函数、改成什么接口、加哪两条红测”**  
写成可执行清单。