好，直接给你**能写代码的 3D hydro 离散方法**，不再谈 plan。下面这版是按 **Woo 论文的 DEC3D hydro 主线** 和你现在已经锁定的 `P1` 规范来写的：**单流体 Euler 主系统 + `\chi_e = P_e^{3/5}` 标量平流 + HLLC + PPM + 球坐标几何源项 + 径向 ALE + macro-zoning + radial MPI**。论文明确说，DEC3D 是 **full spherical, 3-D parallel Eulerian code**，电子压强方程 `∂t P_e^{3/5} + ∇·(v P_e^{3/5}) = 0` 作为 **scalar advection equation** 在单流体 HLLC 里求，再用 **PPM** 把一阶结果提升到三阶；macro-zoning 用来放宽显式 hydro 的 Courant 限制；moving mesh 的径向面通量与 HLLC 集成。fileciteturn23file1turn23file7turn22file3turn22file6

---

# 1. 先解什么方程

论文在球坐标下的 hydro 方程写成：  
\[
\frac{\partial \mathbf Q}{\partial t}
+\frac{1}{r^2}\frac{\partial (r^2 \mathbf F_r)}{\partial r}
+\frac{1}{r\sin\theta}\frac{\partial (\sin\theta\,\mathbf F_\theta)}{\partial \theta}
+\frac{1}{r\sin\theta}\frac{\partial \mathbf F_\phi}{\partial \phi}
=\mathbf S ,
\]
其中状态向量和几何源项是：  
\[
\mathbf Q =
\begin{bmatrix}
\rho\\
\rho v_r\\
\rho v_\theta\\
\rho v_\phi\\
\varepsilon
\end{bmatrix},
\qquad
\mathbf S =
\begin{bmatrix}
0\\
2P/r+\rho(v_\theta^2+v_\phi^2)/r\\
-\rho v_\theta v_r/r+\cot\theta\,(\rho v_\phi^2+P)/r\\
-\rho v_r v_\phi/r-\cot\theta\,\rho v_\theta v_\phi/r\\
0
\end{bmatrix},
\]
而三个方向的物理通量是  
\[
\mathbf F_r =
\begin{bmatrix}
\rho v_r\\
\rho v_r^2+P\\
\rho v_r v_\theta\\
\rho v_r v_\phi\\
h v_r
\end{bmatrix},
\quad
\mathbf F_\theta =
\begin{bmatrix}
\rho v_\theta\\
\rho v_r v_\theta\\
\rho v_\theta^2+P\\
\rho v_\theta v_\phi\\
h v_\theta
\end{bmatrix},
\quad
\mathbf F_\phi =
\begin{bmatrix}
\rho v_\phi\\
\rho v_r v_\phi\\
\rho v_\theta v_\phi\\
\rho v_\phi^2+P\\
h v_\phi
\end{bmatrix},
\]
其中  
\[
\varepsilon=\frac{P}{\gamma-1}+\frac12\rho(v_r^2+v_\theta^2+v_\phi^2),
\qquad
h=\varepsilon+P.
\]
这些就是你要在 `P1` 里离散的**主方程**。fileciteturn23file3turn23file5

如果按你当前项目规范来实现，全局 authoritative truth 不是论文里的 `Q=[ρ,ρv_r,ρv_θ,ρv_φ,\varepsilon]` 直接照搬，而是：
\[
[\rho,\; mom_r,\; mom_\theta,\; mom_\phi,\; E_{fluid,total}]
\]
再加一个 **operator-local**
\[
\chi_e = P_e^{3/5}.
\]
也就是说，hydro 真正推进的局部守恒包应当是：
\[
[\rho,\; mom_r,\; mom_\theta,\; mom_\phi,\; E_{fluid,total},\; \chi_e].
\]
这正是你已经写进 `dec3d-hydro-design.md` 和 `phase-p1-implementation-plan.md` 的约束。fileciteturn22file5turn22file1

---

# 2. 网格与几何量怎么离散

DEC3D 的球坐标必须用 **cell-centered** 网格；论文明确强调球坐标里取 `α = 0.5`，就是单元中心格式，以避开球坐标奇点。fileciteturn22file2turn22file4

设：

- 径向面：\(r_{i-1/2}, r_{i+1/2}\)
- 极角面：\(\theta_{j-1/2}, \theta_{j+1/2}\)
- 方位角面：\(\phi_{k-1/2}, \phi_{k+1/2}\)

单元中心取：
\[
r_i=\frac{r_{i-1/2}+r_{i+1/2}}{2},\quad
\theta_j=\frac{\theta_{j-1/2}+\theta_{j+1/2}}{2},\quad
\phi_k=\frac{\phi_{k-1/2}+\phi_{k+1/2}}{2}.
\]

真正写有限体积更新时，你最需要的是**体积和三个面面积**：

单元体积：
\[
V_{ijk}
=
\frac{1}{3}\left(r_{i+1/2}^3-r_{i-1/2}^3\right)
\left(\cos\theta_{j-1/2}-\cos\theta_{j+1/2}\right)
\left(\phi_{k+1/2}-\phi_{k-1/2}\right).
\]

径向面面积：
\[
A^{(r)}_{i\pm1/2,j,k}
=
r_{i\pm1/2}^2
\left(\cos\theta_{j-1/2}-\cos\theta_{j+1/2}\right)
\Delta\phi_k.
\]

极角面面积：
\[
A^{(\theta)}_{i,j\pm1/2,k}
=
\frac12\left(r_{i+1/2}^2-r_{i-1/2}^2\right)
\sin\theta_{j\pm1/2}
\Delta\phi_k.
\]

方位角面面积：
\[
A^{(\phi)}_{i,j,k\pm1/2}
=
\frac12\left(r_{i+1/2}^2-r_{i-1/2}^2\right)
\Delta\theta_j.
\]

这些公式是球坐标几何直接积分出来的；真正代码里要把它们都缓存起来，不要每次 sweep 临时重算。论文把 DEC3D 定义为全球球坐标 3D 代码，且后面的 macro-zoning 判据也是基于这些弧长与几何尺度来的。fileciteturn23file7turn22file6

---

# 3. 静态网格上的 3D 有限体积更新

先说**不带 ALE、不带 macro-zoning** 的最基本 3D hydro 更新。

对任意守恒变量向量 \(\mathbf U\)（这里就是 \([\rho,mom_r,mom_\theta,mom_\phi,E_{fluid,total},\chi_e]\)），一时间步更新写成：

\[
(V\mathbf U)^{n+1}_{ijk}
=
(V\mathbf U)^n_{ijk}
-\Delta t
\Big[
A^{(r)}_{i+1/2}\mathbf F^{*}_{r,i+1/2}
-
A^{(r)}_{i-1/2}\mathbf F^{*}_{r,i-1/2}
\Big]
\]
\[
-\Delta t
\Big[
A^{(\theta)}_{j+1/2}\mathbf F^{*}_{\theta,j+1/2}
-
A^{(\theta)}_{j-1/2}\mathbf F^{*}_{\theta,j-1/2}
\Big]
-\Delta t
\Big[
A^{(\phi)}_{k+1/2}\mathbf F^{*}_{\phi,k+1/2}
-
A^{(\phi)}_{k-1/2}\mathbf F^{*}_{\phi,k-1/2}
\Big]
+\Delta t\,V_{ijk}\,\mathbf S_{ijk}.
\]

如果你做的是 phase-one 的 split hydro，那么建议按你之前定好的制度，用：

- source half-step
- `r` sweep
- `theta` sweep
- `phi` sweep
- source half-step

也就是：
\[
U^n \xrightarrow{S/2}
U^{(a)} \xrightarrow{r}
U^{(b)} \xrightarrow{\theta}
U^{(c)} \xrightarrow{\phi}
U^{(d)} \xrightarrow{S/2}
U^{n+1}.
\]

这个顺序和你之前锁定的 `hydro3d` 设计是一致的，默认 sweep 次序是 `r -> theta -> phi`。fileciteturn22file3

---

# 4. 每个方向都要“局部一维化”，这才是真 3D hydro

这是你现在最需要钉死的地方。

**真正的 3D hydro** 不是“在 3D 数组上做径向一维更新”，而是对三个方向都构造局部一维 Riemann 问题。

对每个方向，你都要定义 direction-local primitive：

\[
W_{\text{dir}}=(\rho, v_n, v_{t1}, v_{t2}, P, \chi_e).
\]

三个方向的映射是：

### r 方向
\[
v_n=v_r,\quad v_{t1}=v_\theta,\quad v_{t2}=v_\phi.
\]

### \(\theta\) 方向
\[
v_n=v_\theta,\quad v_{t1}=v_r,\quad v_{t2}=v_\phi.
\]

### \(\phi\) 方向
\[
v_n=v_\phi,\quad v_{t1}=v_r,\quad v_{t2}=v_\theta.
\]

这就是为什么 Woo 论文特别强调：HLLC 在 3D 里必须保留接触间断，因为**切向速度在界面两侧是要被带着走的**，这正是 3D shear flow 的来源。论文明说：HLL 的 two-wave model 对 1D 才完整，而 3D 必须用 Toro 的 three-wave HLLC，因为接触间断要把切向速度结构保留下来。fileciteturn23file8turn22file8

---

# 5. HLLC 怎么离散

你要的 HLLC，最少要长成下面这样。

设左、右界面态是 \(W_L, W_R\)，法向速度分别是 \(v_{n,L}, v_{n,R}\)，声速分别是 \(c_L,c_R\)。先用 Davis 波速估计给出：

\[
S_L = \min(v_{n,L}-c_L,\; v_{n,R}-c_R),
\]
\[
S_R = \max(v_{n,L}+c_L,\; v_{n,R}+c_R).
\]

然后接触波速度：
\[
S_*=
\frac{
P_R-P_L+\rho_L v_{n,L}(S_L-v_{n,L})-\rho_R v_{n,R}(S_R-v_{n,R})
}{
\rho_L(S_L-v_{n,L})-\rho_R(S_R-v_{n,R})
}.
\]

接触波两边压强连续：
\[
P_L^*=P_R^*=P^*,
\]
法向速度连续：
\[
v_{n,L}^*=v_{n,R}^*=S_*.
\]

Woo 论文还明确给出：3D HLLC 里切向速度在各自星区保持与本侧一致，也就是：
\[
v_{t1,K}^*=v_{t1,K},\qquad
v_{t2,K}^*=v_{t2,K},\qquad K=L,R.
\]
这正是 3D HLLC 与 1D HLL 的差别。fileciteturn23file8

然后星区密度：
\[
\rho_K^*
=
\rho_K \frac{S_K-v_{n,K}}{S_K-S_*},
\qquad K=L,R.
\]

星区总能量：
\[
E_K^*
=
\frac{
(S_K-v_{n,K})E_K - P_K v_{n,K} + P^* S_*
}{
S_K-S_*
}.
\]

构造星区守恒向量 \(U_K^*\) 之后，HLLC 通量就是标准四分支：

\[
F_{\text{HLLC}}=
\begin{cases}
F_L, & 0\le S_L,\\[4pt]
F_L + S_L(U_L^*-U_L), & S_L\le 0 \le S_*,\\[4pt]
F_R + S_R(U_R^*-U_R), & S_*\le 0 \le S_R,\\[4pt]
F_R, & S_R\le 0.
\end{cases}
\]

这条式子对 \(r,\theta,\phi\) 三个方向都成立，只是 `normal/tangential` 的映射不同。

---

# 6. `chi_e` 的离散：必须跟 HLLC 接触波走

这一条一定别让 Codex 自由发挥。

论文结论里已经写得很清楚：
\[
\partial_t P_e^{3/5} + \nabla\cdot(\vec v P_e^{3/5}) = 0
\]
在单流体 HLLC 里作为 **scalar advection equation** 来解。fileciteturn23file1turn21file3

所以在 hydro 内部，你必须先从 canonical 电子通道恢复：
\[
P_e=(\gamma-1)E_{electron},
\qquad
\chi_e=P_e^{3/5}.
\]

然后把 \(\chi_e\) 当成随流标量，**用和 HLLC 接触/迎风一致的选择**来出通量。  
不要单独写一个 ad hoc 的 upwind flux。

最稳的实现方式是把它当作一个和密度同型的守恒变量来做：

\[
F_{\chi} =
\begin{cases}
v_{n,L}\chi_L, & 0\le S_L,\\[4pt]
F_{\chi,L}+S_L(\chi_L^*-\chi_L), & S_L\le 0\le S_*,\\[4pt]
F_{\chi,R}+S_R(\chi_R^*-\chi_R), & S_*\le 0\le S_R,\\[4pt]
v_{n,R}\chi_R, & S_R\le 0,
\end{cases}
\]
其中
\[
\chi_K^*=\chi_K \frac{S_K-v_{n,K}}{S_K-S_*}.
\]

hydro 结束后再走回写链：
\[
P_e = \chi_e^{5/3},
\qquad
E_{electron} = \frac{P_e}{\gamma-1}.
\]

这就是你当前项目规范里“`chi_e` operator-local、通过 `chi_e -> P_e -> E_electron` authorized writeback”的数学实现。fileciteturn22file5turn22file1

---

# 7. PPM 怎么做，才算“真 PPM”

你的 `P1 Numerical Non-Negotiables` 已经写得很对：**必须是真 PPM，不是 PLM，也不是只有 admissibility 壳。**fileciteturn23file14

我建议你在代码里按下面这套步骤做 direction-local characteristic PPM。

## 7.1 沿 sweep 方向抽线
对每个 sweep 方向，固定另外两维，只在一条 line 上重构：

- `r` sweep：固定 `(j,k)`，沿 `i`
- `theta` sweep：固定 `(i,k)`，沿 `j`
- `phi` sweep：固定 `(i,j)`，沿 `k`

用的变量是：
\[
W_{\text{dir}}=(\rho, v_n, v_{t1}, v_{t2}, P, \chi_e).
\]

## 7.2 特征分解
对局部 1D Euler + tangential + scalar 系统，特征速度集合是：
\[
\lambda = \{v_n-c,\; v_n,\; v_n,\; v_n,\; v_n,\; v_n+c\}.
\]

其中：

- 左、右声波：\(v_n\pm c\)
- 接触波：\(v_n\)
- 两个切向剪切波：\(v_n\)
- 标量 \(\chi_e\) 也跟 \(v_n\)

所以做 PPM 时，不要直接在物理变量上粗暴插值，最好先投影到 characteristic variables，再做抛物重构和限制。

## 7.3 抛物线重构
对每个单元 \(i\)，构造左右边界值 \(q_{L,i}, q_{R,i}\) 和抛物系数 \(q_{6,i}\)：

\[
q_i(x) = q_{L,i}
+ \xi \Big[(q_{R,i}-q_{L,i}) + q_{6,i}(1-\xi)\Big],
\qquad 0\le \xi \le 1,
\]
其中
\[
q_{6,i}=6q_i - 3(q_{L,i}+q_{R,i}).
\]

边界值的生成可以用 Colella–Woodward 风格的三点/四点 stencil，但一定要做：

- monotonicity constraint
- positivity/admissibility 检查

## 7.4 特征追踪到界面
在时间步 \(\Delta t\) 内，不同波只从抛物线的一部分区域传播到界面。所以要按特征速度做 trace，得到界面左右态：

\[
W_{i+1/2}^{L,\;traced},\qquad W_{i+1/2}^{R,\;traced}.
\]

如果你先不想一次把 characteristic tracing 做到最满，也至少要保证：

- 左右界面态来自**抛物线**而不是线性重构
- troubled-cell 时直接降到 first-order HLLC
- 不是“把坏掉的高阶界面态裁一下继续冒充 PPM”

这条边界在你前面的 design 里也已经写死了。fileciteturn22file3turn23file14

## 7.5 什么时候降阶
一旦任一候选界面态不满足：

- \(\rho \ge \rho_{min}\)
- \(P \ge P_{min}\)
- \(\chi_e \ge \chi_{min}\)
- 从 \(\chi_e\) 恢复出的 \(P_e\) 与总压 \(P\) 仍可形成可行热库

或者出现明显非单调/troubled-cell，就直接：

- 用当前 sweep 方向的 **cell-centered state**
- 重新构造 **piecewise constant left/right states**
- 走 **first-order HLLC**

不要对失败的高阶界面态做 clip 后继续冒充高阶。这个 fallback 语义你前面也锁过。fileciteturn22file3

---

# 8. 球坐标几何源项怎么离散

论文已经把球坐标几何源项 \(\mathbf S\) 给出来了。对 `P1`，最稳的做法就是按你前面定好的规则：

- **source half-step**
- `r` sweep
- `theta` sweep
- `phi` sweep
- **source half-step**

也就是 Strang-like 居中源项更新。

在 phase-one，你完全可以把 source step 限定成：

- 只更新真正进入 \(\mathbf S\) 的动量分量
- 不改 `rho`
- 不改 `E_fluid_total`
- 不改 `chi_e`

这和你之前锁定的 hydro3d contract 一致。fileciteturn23file3turn22file3

单元级离散就是：

\[
(m_r)^{n+1/2} = (m_r)^n + \frac12 \Delta t\, S_r(U^n),
\]
\[
(m_\theta)^{n+1/2} = (m_\theta)^n + \frac12 \Delta t\, S_\theta(U^n),
\]
\[
(m_\phi)^{n+1/2} = (m_\phi)^n + \frac12 \Delta t\, S_\phi(U^n).
\]

别在 source step 里偷偷加第二次 ALE 修正。你前面的 design 已经把这条写死了。fileciteturn22file3

---

# 9. 真正的 3D advance 必须怎么接

这是你刚才一直在骂我、也是最关键的点：**不要把 3D 数组上的径向更新当成 3D hydro。**

真正的 full directional static-grid hydro 至少要这样做：

## 9.1 `r` sweep
- ghost ready
- 抽 `(j,k)` 方向线
- PPM / first-order 构造界面态
- HLLC 求 \(\mathbf F_r^*\)
- 径向面更新

## 9.2 `theta` sweep
- refresh needed ghosts
- 抽 `(i,k)` 方向线
- 用 \((\rho,v_\theta,v_r,v_\phi,P,\chi_e)\) 做局部 1D 问题
- HLLC 求 \(\mathbf F_\theta^*\)
- 极角面更新

## 9.3 `phi` sweep
- refresh needed ghosts
- 抽 `(i,j)` 方向线
- 用 \((\rho,v_\phi,v_r,v_\theta,P,\chi_e)\) 做局部 1D 问题
- HLLC 求 \(\mathbf F_\phi^*\)
- 方位面更新

然后把方向局部通量映射回 canonical 动量分量：

- `r` 方向：`mom_n -> mom_r`
- `theta` 方向：`mom_n -> mom_theta`
- `phi` 方向：`mom_n -> mom_phi`

这一步在你之前的 hydro3d design 里其实已经写得很清楚了，只是你现在的新项目 plan 还没把“什么时候必须落地真正的 `theta/phi` advance”写得够硬。fileciteturn22file3turn22file5

---

# 10. 径向 ALE 怎么接进去

论文对 moving mesh 的核心结论是：**把 moving-mesh 通量直接并到 HLLC 界面通量里，比把网格平移和 hydro 分开做更不耗散。**fileciteturn23file1turn23file9

所以，径向面通量必须改成：

\[
F_r^{ALE}=F_r^{HLLC}-w_{face}U_*.
\]

其中：

- \(w_{face}\)：径向面速度
- \(U_*\)：HLLC 接触波选中的界面星区守恒状态

实现上最稳的是：

1. 先用 shell-CM velocity 给出外边界速度  
   论文和 notes 都强调外边界更新是  
   \[
   \frac{dR}{dt} = \beta v_{CM}^{shell}(t^n)
   \]
   的思路。fileciteturn22file2
2. 再按同伦缩放或你已经定好的 radial-only mesh proposal 给出各径向面速度 \(w_{r,i+1/2}\)
3. 只在 **`r` 向 Riemann 通量**里使用 ALE 修正
4. canonical mesh commit 仍由 runtime 做，不让 hydro 私自改 mesh

---

# 11. macro-zoning 怎么并到 hydro 里

这是论文式 `P1` 不能少的部分。Woo 论文明确说，macro-zoning 是为了解决球坐标近极点/近原点处显式 PPM hydro 的小时间步问题；判断依据是弧长：
\[
\Delta S_\theta = r_i\Delta\theta,\qquad
\Delta S_\phi = r_i\sin\theta_j\Delta\phi,
\]
局部显式时间步受：
\[
\Delta t=\min[\Delta r/v_{max},\ \Delta S_\theta/v_{max},\ \Delta S_\phi/v_{max}]
\]
约束。fileciteturn22file0turn22file6

因此，你的 macro-zoning 必须这样并：

### 11.1 每个时间步先做 fine mesh detection
对每个 radial layer 计算：
- \(\Delta S_\theta\)
- \(\Delta S_\phi\)
- \(v_{max}=\max(v_r\pm c_s,\ v_\theta\pm c_s,\ v_\phi\pm c_s)\)

### 11.2 生成 coarse map
按 thesis-style 判据，把弧长小于 \(\Delta r/2\) 的角向细网格合并成粗网格。  
而且 `phi` coarse factor 必须允许依赖 `theta-band`，不能只按 shell 一个数。这个你前面的 design 也写死了。fileciteturn22file3turn22file6

### 11.3 守恒 restrict 到 coarse hydro package
如果你按你当前项目规范做，我建议：
- 对 \([\rho,mom_r,mom_\theta,mom_\phi,E_{fluid,total},\chi_e]\)  
  做 **cell integral** 的守恒 transfer
- 这是你当前工程化重建的正确做法  
  （论文在文字上说的是 prolongation/restriction of primitive variables，但为了保证有限体积守恒，你现在这套设计用 conservative transfer 更稳——这应记入 assumption ledger）

### 11.4 在 coarse hydro workspace 上执行完整显式 hydro
包括：
- source half-step
- `r/\theta/\phi` sweeps
- source half-step

### 11.5 保守 prolong 回 fine canonical state
然后再做 fine-grid 上的最终 thermodynamic recovery / writeback。

---

# 12. ghost / MPI 怎么服务 3D hydro

如果你按当前论文式架构做的是 radial MPI 分解，那么 ghost contract 应该固定成：

1. `exchange_radial_halo(...)`
2. `fill_phi_periodic_and_theta_pole(...)`

这一步对 fine 和 coarse workspace 都成立。  
而且在 `fine -> coarse` 之后，你还必须做一次 **coarse ghost bootstrap / refresh**，保证进入任一方向 sweep 前，该方向 reconstruction 所需的 ghost 都是当前态的。这条你前面的 `hydro3d` design 已经写死了。fileciteturn22file3

极点 remap 时要特别小心：

- 标量：不变
- 向量分量：按局部基向量做符号和分量变换

不然 `theta/phi` 方向的 3D advance 一接上，极点附近马上会出伪波。

---

# 14. 一句话收束

如果你要一句最能指导代码的话，那就是：

> **真正的 3D hydro 离散，不是“在 3D 数组上做径向 HLLC”，而是：**
> **在球坐标有限体积框架下，对 `r/\theta/\phi` 三个方向分别构造局部一维 Riemann 问题，用 direction-local HLLC + PPM 出界面通量，再加球坐标几何源项；然后只在径向通量中并入 ALE 修正，并在 thesis-style macro-zoning 的 coarse hydro workspace 上执行这整条显式更新。**
