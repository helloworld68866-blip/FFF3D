下面给你一个**可直接用于 DEC3D 减速阶段起算的合成 1D 剖面方案**。它的目标不是复现某个具体 shot 的 LILAC 输出，而是构造一个物理上自洽、数值上稳定、能正确触发壳层减速和内界面 RT 增长的初始态。

核心判断是：你现在不该从“激光烧蚀—加速阶段”开始补全，而应构造一个**壳层已经达到最大内爆速度、内壳面刚进入减速阶段**的剖面。Woo 论文的 DEC3D 数据库正是这样做的：OMEGA 77068 的基准起点为 \(t_0=2.28\ \mathrm{ns}\)，此时壳层达到最大内爆速度，内壳面在 \(r_0=68\ \mu\mathrm{m}\)，并在这里加载径向速度扰动。fileciteturn3file16

---

## 1. 物理图像：你要初始化的不是靶丸，而是“减速阶段初态”

ICF 中，外层材料被激光或 X 射线烧蚀，产生类似火箭的反冲，将剩余壳层向中心加速；NIF 文献中给出的典型内爆速度可接近 \(400\ \mathrm{km/s}\)，而停滞阶段则是壳层和 DT 燃料的动能转化为中心低密度热斑及其外侧高密度燃料层的内能。citeturn843636view0turn843636view2 减速阶段的本质就是：中心热斑压力开始反推高速内爆壳层，内壳面发生 deceleration RT；非对称会留下 residual kinetic energy，从而降低热斑压力和产额。PRL 的 OMEGA 高性能直接驱动实验也指出，相对一维对称模拟，实验推断热斑压力约低 40%，三维低模畸变是主要退化机制之一。citeturn595722view1

因此你的初始剖面应满足四个条件：

1. 中心热斑已经被 shock 加热，有 \(P_{\rm hs}\sim\) 数 Gbar；
2. 内壳面附近有冷而密的 DT 壳层，仍以 \(v_0\sim 250\)–\(350\ \mathrm{km/s}\) 向内运动；
3. 热斑速度场应在 \(r=0\) 正则，即 \(v_r(0)=0\)，不能整块平移；
4. 压力、密度、速度都要用若干网格宽度平滑，否则 PPM 会在第一步产生非物理 ringing。

---

## 2. 推荐的 OMEGA-like 合成起始剖面

全部用 cgs 单位实现最方便：

\[
1\ \mu{\rm m}=10^{-4}\ {\rm cm},\qquad 1\ {\rm Gbar}=10^{15}\ {\rm dyn/cm^2}.
\]

建议先采用下面这个参数组：

| 量                                |                                               推荐值 | 说明                                                 |
| --------------------------------- | ---------------------------------------------------: | ---------------------------------------------------- |
| 内壳面 / 热斑边界 \(R_{\rm if}\)  |                                   \(68\ \mu{\rm m}\) | 对齐 Woo 77068 起点                                  |
| 压缩壳层外边界 \(R_{\rm sh,out}\) |                            \(84\)–\(88\ \mu{\rm m}\) | 压缩后壳层厚度约 \(16\)–\(20\ \mu{\rm m}\)           |
| 计算域外边界 \(R_{\rm dom}\)      |                          \(105\)–\(120\ \mu{\rm m}\) | 给外侧低密度缓冲层                                   |
| 最大内爆速度 \(v_0\)              |                       \(3.0\times 10^7\ {\rm cm/s}\) | \(300\ \mathrm{km/s}\)                               |
| 热斑压力 \(P_{\rm hs,0}\)         |                          \(2.5\)–\(4.0\ {\rm Gbar}\) | 初始减速压力，不要直接设到停滞压力                   |
| 热斑中心密度 \(\rho_{\rm hc}\)    |                        \(1.0\)–\(2.0\ {\rm g/cm^3}\) | 边缘可升到 \(2\)–\(4\ {\rm g/cm^3}\)                 |
| 壳层峰值密度 \(\rho_{\rm sh}\)    |                          \(30\)–\(50\ {\rm g/cm^3}\) | 对 OMEGA 压缩燃料+残余壳层合理                       |
| 壳层 adiabat \(\alpha\)           |                                      \(3.0\)–\(3.5\) | Woo 77068 给出 \(\alpha=3.2\) fileciteturn3file16 |
| 外侧缓冲密度 \(\rho_{\rm out}\)   | \(10^{-3}\rho_{\rm sh}\) 到 \(10^{-2}\rho_{\rm sh}\) | 防止真空数值问题                                     |

对于 OMEGA 目标尺度，这个选择也与现代 OMEGA 直接驱动冷冻 DT 靶的数量级一致：OMEGA 是 60 束、30 kJ、351 nm 直接驱动平台；典型冷冻靶约 1 mm 直径，有薄 CH 外壳和几十微米 DT 冰层。citeturn843636view3 Woo 的 77068 例子更具体：外半径 \(430\ \mu{\rm m}\)，CD ablator 厚 \(8\ \mu{\rm m}\)，DT ice 厚 \(50\ \mu{\rm m}\)，收敛比约 20。fileciteturn3file16

---

## 3. 具体函数形式

定义平滑阶跃函数：

\[
H(r;R,w)=\frac12\left[1+\tanh\left(\frac{r-R}{w}\right)\right].
\]

取

\[
w_{\rm if}=2\Delta r\text{ 到 }4\Delta r,\qquad
w_{\rm out}=2\Delta r\text{ 到 }4\Delta r.
\]

然后定义三个区域权重：

\[
C(r)=1-H(r;R_{\rm if},w_{\rm if}),
\]

\[
S(r)=H(r;R_{\rm if},w_{\rm if})\left[1-H(r;R_{\rm sh,out},w_{\rm out})\right],
\]

\[
O(r)=H(r;R_{\rm sh,out},w_{\rm out}).
\]

其中 \(C\) 是中心热斑，\(S\) 是压缩壳层，\(O\) 是外侧低密度缓冲区。

### 3.1 密度剖面

推荐：

\[
\rho_{\rm hs}(r)=\rho_{\rm hc}\left[1+a_\rho\left(\frac{r}{R_{\rm if}}\right)^2\right],
\qquad a_\rho=0.5\text{ 到 }1.0.
\]

壳层密度可先用弱梯度：

\[
\rho_{\rm sh}(r)=\rho_{\rm sh,0}
\left[1-q_\rho \frac{r-R_{\rm if}}{R_{\rm sh,out}-R_{\rm if}}\right],
\qquad q_\rho=0.2\text{ 到 }0.4.
\]

最后：

\[
\rho(r)=C(r)\rho_{\rm hs}(r)+S(r)\rho_{\rm sh}(r)+O(r)\rho_{\rm out}.
\]

一个稳健默认值是：

\[
\rho_{\rm hc}=1.5\ {\rm g/cm^3},\quad
a_\rho=1,\quad
\rho_{\rm sh,0}=35\ {\rm g/cm^3},\quad
q_\rho=0.25,\quad
\rho_{\rm out}=0.03\ {\rm g/cm^3}.
\]

这会给出中心热斑 \(\rho\sim1.5\ {\rm g/cm^3}\)，热斑边缘 \(\rho\sim3\ {\rm g/cm^3}\)，壳层 \(\rho\sim25\)–\(35\ {\rm g/cm^3}\)。

### 3.2 壳层压力：用 adiabat，而不是随便设温度

冷壳层应该低熵、可压缩。用电子简并 Fermi 压力定义 adiabat：

\[
\alpha=\frac{P_{\rm sh}}{P_F},
\]

\[
P_F=\frac{(3\pi^2)^{2/3}}{5}
\frac{\hbar^2}{m_e}n_e^{5/3}.
\]

对等摩尔 DT、完全电离近似：

\[
\bar m_i=\frac{m_D+m_T}{2}\simeq 2.5m_p,\qquad
n_i=\frac{\rho}{\bar m_i},\qquad
n_e=n_i.
\]

壳层压力：

\[
P_{\rm sh}(r)=\alpha P_F[\rho_{\rm sh}(r)].
\]

取 \(\alpha=3.2\) 时，\(\rho_{\rm sh}=30\)–\(40\ {\rm g/cm^3}\) 会给出 \(P_{\rm sh}\sim2\)–\(3.2\ {\rm Gbar}\)，对应壳层温度约几十到一百多 eV，比较适合作为冷压缩燃料壳层。

### 3.3 热斑压力

热斑压力可先设为近似均匀：

\[
P_{\rm hs}(r)=P_{\rm hs,0}.
\]

默认：

\[
P_{\rm hs,0}=3.0\ {\rm Gbar}.
\]

更稳健的范围是：

\[
P_{\rm hs,0}=1.1\text{ 到 }1.5\times P_{\rm sh}(R_{\rm if}^+).
\]

这样热斑已经能反推内壳面，但不会一开始就把壳层炸开。总压力剖面：

\[
P(r)=C(r)P_{\rm hs,0}+S(r)P_{\rm sh}(r)+O(r)P_{\rm out},
\]

其中

\[
P_{\rm out}=10^{-3}P_{\rm hs,0}\text{ 到 }10^{-2}P_{\rm hs,0}.
\]

### 3.4 电子/离子压力与温度

第一版建议热平衡初始化：

\[
P_e=P_i=\frac12P.
\]

然后：

\[
T_e=\frac{P_e}{n_e},\qquad
T_i=\frac{P_i}{n_i}.
\]

如果你的温度单位是 keV：

\[
T_{\rm keV}=\frac{T_{\rm erg}}{1.602176634\times10^{-9}}.
\]

在上面的默认参数下，热斑 \(T_e\simeq T_i\sim1\)–\(3\ {\rm keV}\)，壳层 \(T_e\simeq T_i\sim0.05\)–\(0.15\ {\rm keV}\)。这正好形成“热斑—冷壳”的导热驱动结构，但不会让第一步热传导爆掉。

如果你想更接近真实 ICF，可在内界面附近给电子一个轻微热预热层：

\[
P_e/P=0.5+0.05\exp\left[-\left(\frac{r-R_{\rm if}}{4w_{\rm if}}\right)^2\right],
\]

但第一轮调试不建议这么做。

### 3.5 径向速度剖面

热斑内必须正则：

\[
v_{\rm hs}(r)=-\epsilon_h v_0\frac{r}{R_{\rm if}},
\qquad \epsilon_h=0.2\text{ 到 }0.4.
\]

壳层近似整体向内运动：

\[
v_{\rm sh}(r)=-v_0\left[1-\epsilon_s\frac{r-R_{\rm if}}{R_{\rm sh,out}-R_{\rm if}}\right],
\qquad \epsilon_s=0\text{ 到 }0.15.
\]

总速度：

\[
v_r(r)=C(r)v_{\rm hs}(r)+S(r)v_{\rm sh}(r)+O(r)v_{\rm out}(r),
\]

外侧缓冲区可取

\[
v_{\rm out}(r)=-0.5v_0
\]

或逐渐衰减到 0。第一版建议直接取 \(-0.5v_0\)，避免外侧缓冲层静止造成强剪切。

默认：

\[
v_0=3.0\times10^7\ {\rm cm/s}.
\]

并设：

\[
v_\theta=v_\phi=0.
\]

Woo 的 DEC3D 单模数据库也是在一维基态上加载径向速度扰动，初始横向速度为零。其扰动幅度范围为 \(1\%\)–\(14\%\)，作用在内壳面附近。

---

## 4. 三维扰动加载方式

建议完全沿用 Woo 的思想：只扰动 \(v_r\)，不要一开始扰动 \(\rho\)、\(P\)、\(T_e\)。写成：

\[
v_r(r,\theta,\phi)=v_{r0}(r)
+
A_v v_0 f_\ell(r)Y_\ell^m(\theta,\phi).
\]

其中：

\[
A_v=0.01\text{ 到 }0.14.
\]

径向包络函数可直接用 Woo 论文中的形式：

\[
f_\ell(r)=
\frac12\left(\frac{r}{r_0}\right)^\ell
\left[1-\tanh\left(\frac{r-r_0}{0.015r_0}\right)\right]
+
\frac12\left(\frac{r_0}{r}\right)^\ell
\left[1+\tanh\left(\frac{r-r_0}{0.015r_0}\right)\right],
\]

其中

\[
r_0=R_{\rm if}=68\ \mu{\rm m}.
\]

这个函数在内壳面附近为 1，向内外衰减，避免把太多初始扰动能量放到非界面区域。Woo 论文就是用这个思想在内壳面加载球谐径向速度扰动。

---

## 5. 辐射与 alpha 初值

如果你打开多群辐射，建议不要把 \(U_g\) 全部设为 0。否则第一步会产生很强的非物理辐射弛豫。DEC3D 的辐射方程包含物质–辐射交换项 \(c\kappa_g^P(B_g-\langle U_g\rangle)\)，且 Rosseland/Planck opacity 分别用于扩散与发射吸收。

推荐：

\[
U_g(r)=B_g[T_e(r)]
\]

作为 LTE 初值，至少在热斑和壳层内这样设。外侧缓冲区可设为

\[
U_g=10^{-6}B_g[T_e(R_{\rm sh,out})].
\]

alpha 能量密度初值建议：

\[
\varepsilon_\alpha(r,t_0)=0.
\]

因为减速阶段起点还没有显著燃烧；后续 alpha 源项由

\[
n_Dn_T\langle\sigma v\rangle_{DT}E_{\alpha0}
\]

自然产生。Woo 论文中的 alpha 模型本身就是 one-group diffusion 加源项和电子沉积损失项。

---

## 6. 一组可直接试跑的默认参数

\[
R_{\rm if}=68\ \mu{\rm m},
\quad
R_{\rm sh,out}=86\ \mu{\rm m},
\quad
R_{\rm dom}=115\ \mu{\rm m}.
\]

\[
\rho_{\rm hc}=1.5\ {\rm g/cm^3},
\quad
a_\rho=1.0,
\quad
\rho_{\rm sh,0}=35\ {\rm g/cm^3},
\quad
q_\rho=0.25,
\quad
\rho_{\rm out}=0.03\ {\rm g/cm^3}.
\]

\[
P_{\rm hs,0}=3.0\times10^{15}\ {\rm dyn/cm^2},
\quad
\alpha_{\rm sh}=3.2,
\quad
P_{\rm out}=3.0\times10^{12}\ {\rm dyn/cm^2}.
\]

\[
v_0=3.0\times10^7\ {\rm cm/s},
\quad
\epsilon_h=0.3,
\quad
\epsilon_s=0.1.
\]

