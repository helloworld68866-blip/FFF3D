现在代码里看到的

> `macro_zoning on && use_ppm_reconstruction=true -> fail`

不是因为论文不允许，而是因为你当前实现只做到 **first-order static-grid angular macro coupling**，还没有把 **coarse hydro work view 上的 PPM reconstruction 路径** 做出来。论文结论明确写的是：DEC3D 里 **HLLC 给 first-order hydrodynamic solution，随后用 PPM 提升到 third-order**；同时 **macro-zoning 在 every time step 把 fine mesh 映到 coarser mesh**，目的是放宽显式 hydro 的时间步。也就是说，论文里这两者不是互斥的，而是同一条显式 hydro 路径的一部分。

下面我给你一版**详细、能落代码**的方案，分成两层说：

---

# 一、先把原则说死：PPM 必须在 coarse hydro work view 上做，不是在 fine 上做完再塞回 coarse

这是最关键的一条。

论文锁死的是：

1. 每个时间步先做 macro-zoning detection / coarse generation  
2. **fine 物理量通过 volume-averaging 映到 coarse mesh**  
3. 在 coarse mesh 上做显式 Riemann hydro update  
4. 再把 coarse 上的新物理量映回 fine mesh  
5. macro-zoning **只作用于显式 hydrodynamics**，thermal / alpha / radiation 仍在 fine mesh 上解。fileciteturn64file0turn64file5turn64file6

同时，你自己的 hydro 设计文档已经把 coarse/fine transfer domain 锁成了：

- `rho`
- `m_r`
- `m_theta`
- `m_phi`
- `E`
- `chi_e`

并明确 coarse 表示只是 **short-lived hydro work view**，不是第二套真源。fileciteturn63file18

所以正确口径就是：

> **fine -> coarse 先做 restrict，PPM 在 coarse hydro package 上做，HLLC 也在 coarse 上做，更新后再 prolong 回 fine。**

不能搞成：
- fine 上先 PPM 重构
- 再把高阶界面态投到 coarse

那样 stencil、ghost、limiter、admissibility 全都会错位。

---

# 二、为什么你现在会 fail

因为你当前实现其实是：

- `macro_zoning`：已经有 detect / restrict / prolong / active hydro coupling
- 但这个 coupling 只支持 **first-order static-grid angular updates**
- 一旦 `use_ppm_reconstruction=true`，说明你还没有：
  - coarse `ng=3` ghost bootstrap 真接到 reconstruction 前
  - coarse `theta/phi` line 上的 PPM
  - coarse PPM 界面态喂 HLLC
  - coarse path 上的 troubled-cell / fallback

所以现在 fail 只是说明：

> **你还没把 coarse-grid 高阶路径做出来。**

不是理论上不能。

---

# 三、最小改造方案：先让 `macro + PPM` 真能一起开（不一次性拉满所有东西）

如果你想最快把当前代码从
- “macro on + PPM fail”
推进到
- “macro on + PPM true”

我建议按下面这个最小可行路径。

---

## Step 1：给 coarse hydro package 做完整的 `ng=3` current-state ghost substrate

你现在 fine 路径已经有了：

- `FillHydroRadialGhosts(...)`
- `FillHydroThetaPoleGhosts(..., ng=3)`
- `FillHydroPhiPeriodicGhosts(..., ng=3)`
- `PrepareDirectionalGhosts(...)`

你要做的是让 **coarse package 也能走同样的 ghost contract**。  
这和你之前 hydro3d 设计里写的要求一致：**fine -> coarse restrict 之后，进入任一方向 sweep 前，都必须做 coarse ghost bootstrap / refresh**；ghost 顺序固定是：

1. `exchange_radial_halo(...)`
2. `fill_phi_periodic_and_theta_pole(...)`

而且 PPM 需要 3 层 ghost。Woo 论文也明确写了 spherical PPM hydro 需要三层 ghost，极点和 `phi` 周期 ghost 用 3D 映射关系得到。fileciteturn63file18 fileciteturn63file14

### 你现在代码里应该加什么
对 coarse hydro work view，加这层接口：

```cpp id="j5jwi2"
PrepareDirectionalGhosts(
    coarse_package,
    coarse_geometry,
    Direction dir,
    ghost_layers = 3)
```

这里的 `coarse_geometry` 不能偷用 fine geometry，必须是 coarse map 对应的：
- coarse `theta` band centers / faces
- coarse `phi` band centers / faces
- coarse cell volumes
- coarse face areas

---

## Step 2：PPM reconstruction 改成能吃 coarse hydro package

你现在已经有 `ppm_reconstruction.*`，下一步不是重写它，而是让它接受 coarse package 的 line stencil。

### reconstruction 的输入仍然应该是 direction-local primitive
也就是你以前锁死的：

\[
W_{\text{dir}}=(\rho, v_n, v_{t1}, v_{t2}, P, \chi_e)
\]

并且：
- `E` 不是独立重构的变量
- 界面态总能量从 reconstructed primitive 一致恢复  
这条也在你 hydro3d design 里写死了。fileciteturn63file18

### 你现在该做的
对 coarse hydro package 上每条 line：
- 先用 coarse ghost 准备好的 current-state line
- 再做 PPM reconstruction
- 再做 troubled-cell 检测
- 非法态直接降回 first-order HLLC

不要让 macro path 一上来就绕过 fallback contract。

---

## Step 3：在 coarse hydro path 上把 PPM 界面态喂给 HLLC
当前一阶 macro path 的核心是：

- coarse cell-centered state
- HLLC
- update

要变成：

- coarse PPM left/right interface states
- HLLC(interface states)
- update

也就是现在这条 guard：

```cpp id="5g1f4r"
if (use_ppm_reconstruction) fail;
```

应该变成：

```cpp id="p5x6fz"
if (use_ppm_reconstruction) {
    // coarse-grid PPM path
} else {
    // first-order coarse path
}
```

---

## Step 4：prolong 先保持“守恒优先”，不要一上来做高阶 prolong
这是最重要的策略建议之一。

如果你当前 prolong 还是：
- constant fill
或者
- limited linear angular prolong

我建议你在第一轮 **macro + PPM** 联动时，先保留它。  
因为你现在真正想证明的是：

> **coarse hydro 上的 PPM + HLLC 真工作了。**

不是：
> **coarse hydro + high-order prolong 也同时完美。**

所以第一轮完全可以：

- restrict：守恒
- coarse update：高阶 PPM + HLLC
- prolong：先守恒、可 admissible、但不必高阶

等这一条走通，再升级 prolong。

---

## Step 5：先只开 single-rank + static-grid + angular macro + PPM
如果你是“最快想把当前 fail 改掉”，那这一步的范围应该收得很窄：

### 先允许
- `single-rank`
- `static-grid`
- `macro_zoning = on`
- `use_ppm_reconstruction = true`
- 只对 `theta/phi` 的 coarse angular path 开 PPM

### 先继续禁止
- `macro + PPM + ALE`
- `macro + PPM + real MPI`
- `macro + PPM + advanced prolong`

这样最稳。

---

# 四、和论文最一致的最终方案：`macro + PPM + real MPI + ALE` 一起开

如果你不是只想修当前 fail，而是想一次性把终局路线写清楚，那按论文和你项目设计，完整方案应该是这样：

---

## Phase A：每个时间步先基于 `geom_n` 做 macro map
这条你设计里已经锁死了：

- macro map 读 `geom_n`
- 整步冻结
- 不混 `geom_np1`

并且论文路线是：
- 每个时间步先 detect fine cells
- 再 coarse mesh generation
- 再 prolong / restrict。fileciteturn63file18 fileciteturn64file9

---

## Phase B：fine -> coarse restrict
对 hydro package：

\[
[\rho,\ m_r,\ m_\theta,\ m_\phi,\ E,\ \chi_e]
\]

做**守恒 restrict**。  
这一点在你设计里已经锁死。fileciteturn63file18

---

## Phase C：coarse ghost bootstrap / refresh
在 coarse package 上做：

1. radial halo exchange（real MPI）
2. `phi` periodic
3. `theta` pole remap

而且：
- 进入任一方向 sweep 前
- 当前方向所需 coarse ghost 必须是 current state

这条是你之前 hydro3d 设计写死的，而且论文的 spherical PPM 也要求 3 层 ghost。fileciteturn63file18 fileciteturn63file14

---

## Phase D：coarse hydro update
这一步才是真正的论文式显式 hydro：

1. source step（或你当前 phase-one old-time explicit source）
2. `r` sweep
3. `theta` sweep
4. `phi` sweep

每个 sweep 都是：

- coarse PPM reconstruction
- coarse HLLC
- current coarse geometry / face area / volume
- if `r` sweep and ALE enabled:
  \[
  F_r^{ALE} = F_{upwind}(w_{face}) - w_{face} Q_{upwind}(w_{face})
  \]

也就是：
- **PPM 和 macro coarse update 一起工作**
- **real MPI 通过 coarse halo exchange 提供 ghost**
- **ALE 只进 radial flux**

---

## Phase E：coarse -> fine prolong
更新后的 coarse hydro package 再 prolong 回 fine。  
第一轮可以保持：
- constant fill
或
- limited linear angular prolong

后面再升级。

---

## Phase F：fine-grid recovery / writeback / runtime commit
最后在 fine canonical state 上做：

- final recovery
- thermodynamic feasibility
- writeback
- runtime mesh commit（如果 ALE 开启）

这和你一直锁的 canonical truth 语义是一致的。

---

# 五、如果你真想一次性把 `macro + PPM + real MPI + ALE` 做完，我建议的实现顺序

虽然你说想一次性做完，但我还是建议你在代码里按下面顺序递进，不然非常容易缠成一团：

### 7b.5
**single-rank + static-grid + angular macro + PPM**  
目标：把当前 `use_ppm_reconstruction && macro_zoning` 的 fail 去掉

### 7c.5
**real MPI + angular macro + PPM**  
目标：让 coarse ghost bootstrap 真消费 `ng=3`

### 7c.6
**real MPI + angular macro + PPM + ALE(radial only)**  
目标：径向 ALE correction 和 coarse hydro path 真联动

也就是说，虽然最终目标可以是一口气描述成：

> `macro + PPM + real MPI + ALE`

但代码落地时，最好仍是分 3 个窄纵切。  
否则一旦错了，你根本不知道：
- 是 coarse ghost 错
- PPM 错
- HLLC 错
- ALE branch 错
- 还是 prolong 错

---

# 六、你现在应该先补哪些 failing-first tests

如果你真的要推进“两者一起开”，我建议最少先写这 6 条红测：

## 1. coarse PPM 必须真消费 `ng=3`
`test_macro_zoning_ppm_requires_ng3_ghosts.cpp`

## 2. macro + PPM 在 single-rank static-grid 下能跑
`test_macro_zoning_ppm_single_rank.cpp`

## 3. macro + PPM 不是偷偷 fallback 成一阶
`test_macro_zoning_ppm_not_first_order_fallback.cpp`

## 4. macro + PPM + real MPI 真消费 coarse halo
`test_macro_zoning_ppm_mpi_ghost_bootstrap.cpp`

## 5. macro + PPM + ALE 的 radial branch 真走 `w_face` 选 branch
`test_macro_zoning_ppm_ale_branch_selection.cpp`

## 6. macro + PPM + ALE + real MPI 失败时不污染 canonical truth
`test_macro_zoning_ppm_ale_mpi_transaction.cpp`

---

# 七、一句话总结

如果你只想记一句：

> **不是“PPM 和 macro-zoning 不能一起开”，而是你当前代码还没把“coarse hydro work view 上的 PPM reconstruction + coarse ghost bootstrap + HLLC +（可选）ALE”这条路径做出来。**
>
> **要让它们一起开，正确方法是：**
> 1. **fine -> coarse 守恒 restrict**
> 2. **coarse 当前态 ghost 准备（含 real MPI halo）**
> 3. **coarse PPM reconstruction**
> 4. **coarse HLLC**
> 5. **如果开 ALE，只在 radial coarse flux 上加 moving-face correction**
> 6. **coarse -> fine 守恒 prolong**
> 7. **fine canonical state recovery / writeback / runtime commit**

这条方案和论文路线不冲突，而且是你当前代码最合理的实现路径。  
论文锁死的是：`HLLC -> PPM` 的高分辨显式 hydro、macro-zoning 每步 coarse 化、real radial MPI 架构、以及 moving mesh 径向通量修正；coarse PPM 的具体代码组织、prolong 公式、以及内部 `w_face` predictor 仍然是你当前项目的工程补全，需要继续留在 assumption ledger。