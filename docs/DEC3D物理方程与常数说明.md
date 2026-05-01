下面是一份可直接保存的文档。它只讲 **DEC3D 涉及的物理方程**，以及这些方程中各个常数、系数和物理量**怎么定义、怎么计算、从哪里来**。我刻意不谈任何离散、求解器、网格、时间推进和并行实现。DEC3D 在论文中的物理定位是：**单流体、双温度、三维辐射流体模型**，把系统拆成水动力、电子/离子热传导、电子-离子温度弛豫、多群辐射输运和 alpha 粒子输运五类物理过程。

---

# DEC3D 物理方程与常数说明

## 1. 基本假设与状态变量

论文对等离子体作了三条核心假设。第一，等离子体视为**单流体**，所以质量密度和速度场只有一套；第二，热学上是**双温度**，电子温度 $T_e$ 和离子温度 $T_i$ 分开演化；第三，采用**理想气体、完全电离**近似，且比热比取 $\gamma=5/3$，平均 DT 离子电荷取 $\bar Z=1$。论文写成电子和离子各自满足理想气体关系 $P_{e/i}=n_{e/i}T_{e/i}$；这里的温度常常直接用能量单位表示，所以式子里不总是显式写出 $k_B$。
单流体质量密度定义为
$$
\rho=n_i\bar m_i,\qquad \bar m_i=\frac{m_D+m_T}{2}.
$$
也就是说，离子平均质量由氘和氚的离子质量平均得到。这组假设的物理意思是：流体力学层面只跟踪一套 $\rho,\vec v,P$，但能量交换层面必须额外跟踪电子和离子两个热库，因为电子热传导、辐射耦合、alpha 电子阻曳和电子-离子碰撞都会让 $T_e\neq T_i$。

---

## 2. 水动力主方程

DEC3D 的主欧拉方程组是无粘单流体守恒律：

$$
\partial_t \rho + \nabla\cdot(\rho \vec v)=0, \tag{5.1}
$$

$$
\partial_t(\rho \vec v)+\nabla\cdot(\rho \vec v\otimes \vec v+\hat{\mathbb I}P)=0, \tag{5.2}
$$

$$
\partial_t\left(\frac{P}{\gamma-1}+\frac12\rho v^2\right) +\nabla\cdot\left[ \vec v\left(\frac{\gamma P}{\gamma-1}+\frac12\rho v^2\right)+\vec Q \right] =S. \tag{5.3}
$$

其中 $\hat{\mathbb I}$ 是单位张量，$\otimes$ 是外积；$\vec Q$ 和 $S$ 不是“数值通量”的意思，而是物理上的总热流和总热源项，包含热传导、辐射和 alpha 加热的贡献。

这三式中真正需要外部常数或外部物理输入的部分有四类：

一是 $\gamma=5/3$。这是理想单原子/完全电离等离子体近似下的比热比。

二是 $\bar m_i=(m_D+m_T)/2$。这里需要氘、氚离子质量 $m_D,m_T$。论文把它们当作基本粒子常数使用。

三是总压强

$$
P=P_e+P_i.
$$

如果把温度写成热力学温标，则应理解为

$$
P_e=n_e k_B T_e,\qquad P_i=n_i k_B T_i.
$$

如果把 $T_e,T_i$ 直接写成能量单位，论文就写成 $P=nT$。这是同一物理关系的两种单位写法。前者是我为避免单位混淆补充的说明。论文明确给出了 $P_{e/i}=n_{e/i}T_{e/i}$。

四是能量方程中的 $\vec Q,S$。这两个量本身并不独立定义，而是由后续三类热传输模块共同决定：热传导、辐射、alpha。

---

## 3. 电子和离子热传导

论文把热传导单独写成两条扩散方程：

$$
\partial_t\left(\frac{P_e}{\gamma-1}\right)=\nabla\cdot \kappa_e \nabla T_e, \tag{5.11}
$$

$$
\partial_t\left(\frac{P_i}{\gamma-1}\right)=\nabla\cdot \kappa_i \nabla T_i. \tag{5.12}
$$

这两式的核心不是形式，而是热导率 $\kappa_e,\kappa_i$ 怎么算。

### 3.1 电子热导率

论文采用 Spitzer–Härm 型电子热导率，并加上 Lee–More 简并修正：

$$
\kappa_e = \underbrace{ \frac{20(2/\pi)^{3/2}(k_B T_e)^{5/2}k_B} {m_e^{1/2}\bar Z e^4 \ln\Lambda_{ei}^{LM}} }_{\kappa_L} \underbrace{ \frac{0.095(\bar Z+0.24)}{1+0.24\bar Z} }_{\delta} f_{LM}. \tag{5.232}
$$

论文同时给出

$$
\kappa_i=\kappa_e\sqrt{\frac{m_e}{m_i}}.
$$

也就是说，离子热导率不是单独另造一套模型，而是从电子热导率乘上一个质量比因子得到。

这条式子里涉及的常数和输入如下。

首先是基础常数：

- $k_B$：玻尔兹曼常数。论文在 5.232 中显式出现。

- $e$：元电荷。论文显式出现。

- $m_e$：电子质量。论文显式出现。

- $m_i$：离子质量。论文通过 $\kappa_i=\kappa_e\sqrt{m_e/m_i}$ 使用。

然后是等离子体状态量：

- $T_e$：电子温度。

- $\bar Z$：有效离子电荷。对 DT 完全电离近似，前文取 $\bar Z=1$。

再是两个修正因子：

- $\delta$：Spitzer prefactor，论文把它解释为由于热电效应导致的 Lorentz 气体热导率修正。

- $$f_{LM}$$：Lee–More 简并修正因子。

### 3.2 Lee–More 简并修正因子

论文给出

$$
f_{LM}(n_e,T_e) = 1+\frac{3\pi^5}{51200}\left(\frac{T_F}{T_e}\right)^3\delta^{-2}, \tag{5.233}
$$

其中电子费米温度为

$$
T_F=\frac{\hbar^2}{2m_e k_B}(3\pi^2 n_e)^{2/3}.
$$

因此，$f_{LM}$ 的计算需要电子数密度 $n_e$、电子温度 $T_e$、电子质量 $m_e$、玻尔兹曼常数 $k_B$、约化普朗克常数 $\hbar$，以及前一式中的 $\delta$。

### 3.3 电子-离子 Coulomb 对数

论文显式给出了 Lee–More 版 Coulomb 对数的结构：

$$
\ln\Lambda_{ei}^{LM} = \max\left[ \frac12 \ln\!\left(1+(b_{\max}/b_{\min})^2\right),\ 2 \right].
$$

因此它需要两个截断尺度 $b_{\max}$ 和 $b_{\min}$。论文在当前可见文本里没有继续展开它们的具体表达式，但从这条式子本身已经能看出：$\ln\Lambda_{ei}^{LM}$ 不是常数，而是由微观碰撞几何和等离子体状态决定的。

论文还提到另一套 Spitzer Coulomb 对数 $\ln\Lambda_{ei}^{Spitzer}$，并在 benchmark 中比较了 Spitzer 与 Lee–More 两种选择。可惜当前解析文本没有把 Eq. (5.236) 的完整代数式展开出来，所以这里不强行补写，以免误引。能确定的是：DEC3D 物理讨论中两种 Coulomb 对数都被使用过，而 Lee–More 版本配合 $f_{LM}$ 与 LILAC 的一致性更好。

---

## 4. 电子–离子温度弛豫

电子与离子的能量交换在论文中写成一对局部弛豫方程：

$$
\partial_t T_e=-\frac{1}{\tau_{ei}}(T_e-T_i), \tag{5.13}
$$

$$
\partial_t T_i=-\frac{1}{\tau_{ei}}(T_i-T_e). \tag{5.14}
$$

这说明在这个模型里，电子和离子通过同一个弛豫时间 $\tau_{ei}$ 向共同温度靠近。

论文还给出这套模型对应的温差衰减规律：

$$
(T_e-T_i)^{n+1}=(T_e-T_i)^n\exp\!\left(-\frac{2\Delta t}{\tau_{ei}}\right). \tag{5.243}
$$

如果只看连续物理含义，它表达的是

$$
T_e-T_i \propto e^{-2t/\tau_{ei}},
$$

也就是电子与离子温差会指数衰减到零。

需要坦率说明一点：论文当前可访问的解析文本没有把 Eq. (5.241)——也就是 $\tau_{ei}$ 的完整代数表达式——完整展现出来。我们能可靠读到的是：$\tau_{ei}$ 来自电子-离子碰撞率，且其倒数 $1/\tau_{ei}$ 足够大时，会快速促成电子和离子接近热平衡。所以在这份文档里，我不凭记忆去补一条“标准 Spitzer $\tau_{ei}$”公式，避免把论文外公式误当成论文原式。

---

## 5. 多群辐射输运

论文从最一般的 Boltzmann 输运方程和频率依赖辐射输运方程出发：

$$
\frac{1}{v}\frac{\partial f}{\partial t}+\vec\Omega\cdot\nabla f=C(f), \tag{5.176}
$$

$$
\frac{1}{c}\frac{\partial I_\nu}{\partial t} +\vec\Omega\cdot\nabla I_\nu = \varepsilon_\nu-\kappa_\nu I_\nu. \tag{5.177}
$$

其中 $I_\nu$ 是频率 $\nu$ 上的谱辐射强度，$\varepsilon_\nu$ 是发射率，$\kappa_\nu$ 是不透明度。

DEC3D 实际采用的是**多群、通量限制的辐射扩散近似**。每个频率群 $g$ 的辐射能量满足：

$$
\partial_t \langle U_g\rangle + \nabla\cdot(\vec v \langle U_g\rangle) + \langle P_g\rangle \nabla\cdot \vec v = \nabla\cdot \bar D_g(\kappa_g^R)\nabla \langle U_g\rangle + c\kappa_g^P(B_g-\langle U_g\rangle). \tag{5.15}
$$

这条方程中最关键的物理系数有五个：$c,\kappa_g^R,\kappa_g^P,\bar D_g,B_g$。

### 5.1 光速 $c$

这是最直接的基础常数。它在辐射–物质交换项 $c\kappa_g^P(B_g-\langle U_g\rangle)$ 中出现。

### 5.2 Rosseland 群平均不透明度

论文给出

$$
\kappa_g^R(\rho,T_e) = \frac{\int_{\nu_g}^{\nu_{g+1}} \frac{\partial B_\nu(\nu,T_e)}{\partial T_e}\,d\nu} {\int_{\nu_g}^{\nu_{g+1}} \frac{1}{\kappa_\nu(\rho,T_e)} \frac{\partial B_\nu(\nu,T_e)}{\partial T_e}\,d\nu}. \tag{5.217}
$$

因此，$\kappa_g^R$ 的输入有：

- 质量密度 $\rho$；

- 电子温度 $T_e$；

- 单频不透明度 $\kappa_\nu(\rho,T_e)$；

- 群边界 $\nu_g,\nu_{g+1}$；

- 黑体谱函数 $B_\nu$。

### 5.3 Planck 群平均不透明度

论文给出

$$
\kappa_g^P(\rho,T_e) = \frac{\int_{\nu_g}^{\nu_{g+1}} \kappa_\nu(\rho,T_e)B_\nu(\nu,T_e)\,d\nu} {\int_{\nu_g}^{\nu_{g+1}}B_\nu(\nu,T_e)\,d\nu}. \tag{5.218}
$$

输入与 $\kappa_g^R$ 类似，但加权方式不同。论文还特别说明：Rosseland 平均更适合光学厚介质，Planck 平均更适合光学薄介质。

### 5.4 单频不透明度 $\kappa_\nu$

这是辐射方程里最重要的“外部物理输入”之一。论文说明 DEC2D/DEC3D 的 LTE tabular opacities 来自 **Astrophysical Opacity Library**。也就是说，$\kappa_\nu$ 不是在代码内部从零推导出来，而是由外部原子物理库随 $(\rho,T_e)$ 查表得到。

### 5.5 群黑体源项 $B_g$

论文把每个频率群的自发射能量密度写成

$$
B_g=a\,(k_B T_e)^4\,b(u_g,u_{g+1}), \tag{5.220}
$$

其中

$$
u_g=\frac{h\nu_g}{k_B T_e}.
$$

这里 $b(u_g,u_{g+1})$ 是该群在整体黑体谱中所占的权重因子，取值在 0 到 1 之间；$a$ 是黑体常数，论文写为

$$
a=\frac{8\pi^5}{15h^3c^3}.
$$

因此 $B_g$ 的计算需要：

- 普朗克常数 $h$；

- 玻尔兹曼常数 $k_B$；

- 光速 $c$；

- 电子温度 $T_e$；

- 群边界频率 $\nu_g,\nu_{g+1}$；

- 归一化 Planck 分布积分得到的群权重 $b(u_g,u_{g+1})$。

---

## 6. alpha 粒子输运

论文将 alpha 粒子输运建模为 **Atzeni one-group diffusion model**，描述的是 alpha 粒子在电子阻曳下减速并沉积能量的过程。它给出的主方程是

$$
\frac{\partial \varepsilon_\alpha}{\partial t} = \nabla\cdot D_\alpha \nabla \varepsilon_\alpha + n_D n_T\langle \sigma v\rangle_{DT} E_{\alpha0} - \frac{\varepsilon_\alpha}{\tau_{\alpha e}}. \tag{5.263}
$$

这里 $\varepsilon_\alpha$ 是 alpha 粒子能量密度，三项右端分别对应：空间扩散、DT 反应产额源项、以及对电子的能量沉积损失。

### 6.1 alpha 出生能量和初速度

论文明说 alpha 粒子出生能量为

$$
E_{\alpha0}=3.5~\text{MeV}.
$$

相应初速度 $v_{\alpha0}$ 由动能关系确定。若按非相对论近似，

$$
v_{\alpha0}=\sqrt{\frac{2E_{\alpha0}}{m_\alpha}}.
$$

这条速度公式是物理补充说明，不是论文显式写出的等式；论文只说 alpha 粒子以出生能量 $E_{\alpha0}=3.5$ MeV 和初速度 $v_{\alpha0}$ 出生。

### 6.2 alpha–电子弛豫时间

论文给出

$$
\tau_{\alpha e} = \frac{3 m_\alpha (k_B T_e)^{3/2}} {8\sqrt{2\pi m_e}\,n_e\,\bar Z_\alpha^2 e^4 \ln\Lambda_{\alpha e}^{Spitzer}}. \tag{5.262}
$$

这条式子告诉你 $\tau_{\alpha e}$ 的输入包括：

- alpha 质量 $m_\alpha$；

- 电子质量 $m_e$；

- 电子温度 $T_e$；

- 电子数密度 $n_e$；

- alpha 电荷数 $\bar Z_\alpha=2$；

- 元电荷 $e$；

- Spitzer 型 alpha–electron Coulomb 对数 $\ln\Lambda_{\alpha e}^{Spitzer}$。

论文还说明，$\ln\Lambda_{\alpha e}^{Spitzer}$ 的形式与电子–离子 Spitzer Coulomb 对数相同，只是把电荷数替换成 $Z_\alpha=2$。

### 6.3 alpha 扩散系数

论文进一步定义电子阻曳导致的扩散平均自由程

$$
\lambda_{\text{drag}}=\frac{v_{\alpha0}\tau_{\alpha e}}{9},
$$

然后定义扩散系数

$$
D_\alpha=v_{\alpha0}\lambda_{\text{drag}}.
$$

所以 $D_\alpha$ 实际上由 $v_{\alpha0}$ 和 $\tau_{\alpha e}$ 决定。

### 6.4 DT 反应源项

源项写成

$$
n_D n_T\langle \sigma v\rangle_{DT} E_{\alpha0}.
$$

其中：

- $n_D,n_T$ 分别是氘、氚离子数密度；

- $\langle \sigma v\rangle_{DT}$ 是 DT 反应率；

- $E_{\alpha0}=3.5$ MeV 是每个 alpha 粒子携带的出生能量。

论文特别说明，$\langle \sigma v\rangle_{DT}$ 来自 **Bosch–Hale model**。因此它也是外部反应率模型输入，而不是由 DEC3D 自己从核反应微观理论现推。

### 6.5 低能后期行为

论文还明确提到：当 alpha 粒子动能降到

$$
E_\alpha<0.5~\text{MeV}
$$

时，其运动将由与离子的横向散射主导。这是一条重要的物理注释，因为它说明 one-group 电子阻曳扩散模型主要覆盖的是高能到中能的主沉积阶段，低能尾部的离子散射物理并没有被单独展开成新的控制方程。

---

## 7. 各类常数的来源总结

为了实际查阅方便，可以把 DEC3D 物理方程中出现的“常数/系数”分成三类。

### 7.1 真正常数

这类量不随时空变化：

$$
\gamma=\frac53,\quad \bar Z=1,\quad \bar Z_\alpha=2,
$$

$$
k_B,\ h,\ \hbar,\ c,\ e,\ m_e,\ m_D,\ m_T,\ m_\alpha,
$$

$$
E_{\alpha0}=3.5\ \text{MeV}.
$$

其中 $\gamma,\bar Z,E_{\alpha0}$ 是论文模型设定；其余是基本物理常数和粒子质量。

### 7.2 由局部等离子体状态计算的量

这类量每个时空点都不同：

$$
\rho=n_i\bar m_i,\quad P_e=n_eT_e,\quad P_i=n_iT_i,
$$

$$
\kappa_e,\ \kappa_i,\ f_{LM},\ T_F,\ \ln\Lambda_{ei}^{LM},
$$

$$
\kappa_g^R,\ \kappa_g^P,\ B_g,
$$

$$
\tau_{\alpha e},\ D_\alpha.
$$

这些量不是独立输入，而是由 $\rho,n_e,n_i,T_e,T_i$ 等局部状态再加一些基本常数计算得到。

### 7.3 外部库或外部物理模型提供的量

这类量论文没有在 DEC3D 内部从头闭合：

- 单频不透明度 $\kappa_\nu(\rho,T_e)$：来自 Astrophysical Opacity Library。

- DT 反应率 $\langle \sigma v\rangle_{DT}$：来自 Bosch–Hale model。

- Spitzer/Lee–More Coulomb 对数的具体截断模型：来自外部等离子体输运理论。

---

## 8. 最终总纲

如果把这份文档再压缩成一句话，DEC3D 的物理闭合就是：
一套单流体欧拉方程负责
$$
\rho,\vec v,P
$$
的整体动力学；两条热扩散方程负责 $T_e,T_i$ 的空间传热；一对局部弛豫方程负责电子–离子之间的能量交换；一组多群辐射扩散方程负责辐射能量密度 $\langle U_g\rangle$ 与物质的吸收/发射耦合；一条 alpha 能量扩散方程负责 DT 反应产生的 alpha 粒子加热。所有这些模块共同决定能量方程中的 $\vec Q$ 和 $S$。