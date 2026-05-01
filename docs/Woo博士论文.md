# The Influence of Three Dimensional Effects on the Core Dynamics and Nuclear Measurements of Inertial Fusion Implosions

by 

Ka Ming Woo 

Submitted in Partial Fulfillment of the 

Requirements for the Degree 

Doctor of Philosophy 

Supervised by Professor Riccardo Betti 

Department of Physics and Astronomy 

Arts, Sciences and Engineering 

Schools of Arts and Sciences 

University of Rochester 

Rochester, New York 

2019 

# Dedication

To my mother, and people who inspired me. 

# Table of Contents

Biographical Sketch vi 

Acknowledgments xi 

Abstract xii 

Contributors and Funding Sources xiv 

List of Tables xvi 

List of Figures xvii 

1 Research Objective 1 

1.1 Introduction to ICF . 1 

1.1.1 Rayleigh-Taylor instability . . . . 6 

1.1.2 Nuclear Diagnostics . . . . 7 

1.2 Understanding 3-D implosion asymmetries . . . 9 

2 Impact of Residual Kinetic Energy 12 

2.1 Motivation . 12 

2.2 DEC3D simulation database . 12 

2.3 Property of low-mode hot-spot volumes . . 15 

2.4 A 3-D hot-spot model 21 

2.4.1 Dynamics of 3-D hot-spot volumes . . . 21 

iv 

2.4.2 An adiabatic invariant 3-D model for low modes . . . . . . 27 

2.4.3 A non-adiabatic invariant 3-D model for high modes . . . . 34 

2.5 Conclusion 36 

# 3 Impact of Hot-Spot Flow Anisotropy 38

3.1 Motivation . 38 

3.2 Ion temperature measurement 39 

3.2.1 Properties of non-stagnating hot-spot flows . . . . . . 39 

3.2.2 Neutron energy spectrum model . . . . . . 44 

3.3 Neutron energy spectra . 48 

3.3.1 Single-mode ion-temperature characteristics . . . . . . 48 

3.3.2 Skewness and kurtosis for non-Gaussian spectra . . . . . . . 53 

3.4 A velocity variance decomposition technique . . 56 

3.4.1 Physical origin of velocity variance . . . . 56 

3.4.2 Single-mode velocity variance . . . . 62 

3.4.3 Multi-mode velocity variance . . . 75 

3.5 Hot-spot residual kinetic energy 78 

3.6 Conclusion . 82 

# 4 Impact of Isotropic Flows within the Hot Spot 84

4.1 Motivation . 84 

4.2 Physical meaning of variance and covariance 85 

4.3 Properties of isotropic velocity variance . . . 91 

4.4 Effects on DD/DT ion-temperature ratio . . . . . 97 

4.5 Multi-mode DD/DT ion-temperature ratio . . . . 102 

4.6 Diagnosing for hot spot flow isotropy . . . . . 107 

4.6.1 The first approximate closure . . 107 

4.6.2 Error propagation analysis for the first closure . . . . . . . . 111 

4.6.3 The second approximate closure . . . . . . 113 

4.6.4 The third approximate closure . . . . . . 114 

4.6.5 Integrated performance in multi-mode simulation . . . . . . 116 

4.6.6 Error propagation analysis in the integration test . . . . . . 119 

4.6.7 $T _ { \mathrm { D D } } / T _ { \mathrm { D T } }$ analysis . . . . . 120 

4.7 Effects of hot-spot flow anisotropy on yield degradation . . . . . . . 123 

4.8 Conclusion . 130 

# 5 DEC3D Computer Code 131

5.1 Motivation . 131 

5.2 Physical models 134 

5.2.1 Governing equations . . . . . 134 

5.2.2 Cartesian & spherical mesh discretization . . . . . . . . 138 

5.2.3 Macro-zoning . . . . 144 

5.3 Hydrodynamics 147 

5.3.1 Conservative moving mesh . . . . . 147 

5.3.2 Approximate Riemann solvers . . . . . 173 

5.4 Radiation transport . 188 

5.4.1 Multi-group flux-limited radiation diffusion . . . . . . . . . . 189 

5.4.2 Implementation in DEC2D and DEC3D . . . . . . . . . . . 201 

5.4.3 Boundary condition . . . . . . . 204 

5.4.4 Parallel multi-group radiation transport in DEC2D . . . . . 207 

5.4.5 Radiative ablative stabilization of RT instabilities . . . . . . 208 

5.5 Electron and ion heat conductions . 210 

5.5.1 Spitzer H¨arm thermal diffusion . . . . . . 210 

5.5.2 Electron and ion equilibration . . . . . 212 

5.5.3 Fully implicit second-order in space discretization . . . . . . 213 

5.5.4 Benchmark tests . 218 

5.6 Alpha Particle Transport . . 222 

vi 

5.7 Conclusion . 223 

6 Conclusion 225 

Bibliography 227 

# Biographical Sketch

The author was born in China. He moved to Hong Kong in 1994. He attended The Chinese University of Hong Kong from 2007 to 2012. He graduated with a Bachelor of Science degree in the theoretical stream of Physics and a minor in Mathematics in 2010. He received a Master of Philosophy degree in Physics in 2012. During the Master study, he was advised by Professor Simon Yu on stability issues of transportating heavy-ion beams through ion induction accelerators. He developed a method to diagnose the ion beam energy published in Nuclear Instruments and Methods, and a technique to stabilize the ion beam transportation published in Physical Review Accelerators and Beams. Inspired by the heavy-ion inertial fusion, he came to the University of Rochester in 2012 to conduct research on hydrodynamic instabilities of inertial confinement fusion implosions advised by Professor Riccardo Betti at the Laboratory of Laser Energetics as a Horton Fellow. He received a Master of Arts degree in Physics in 2014. 

# Presentations and Publications

# First-Author Publications



1. “The Influence of Three Dimensional Effects on the Core Dynamics and Nuclear Measurements of Inertial Fusion Implosions”, K. M. Woo, et al., manuscript in final preparation 





2. “Impact of Three-Dimensional Hot-Spot Flow Asymmetry on Ion-Temperature Measurements in Inertial Confinement Fusion Experiments”, K. M. Woo, R. Betti, D. Shvarts, O. M. Mannion, D. Patel, V. N. Goncharov, K. S. Anderson, P. B. Radha, J. P. Knauer, A. Bose, V. Gopalaswamy, A. R. Christopherson, E. M. Campbell, J. Sanz and H. Aluie, Physics of Plasmas, volume 25, 102710 (2018) 





3. “Effects of Residual Kinetic Energy on Yield Degradation and Ion Temperature Asymmetries in Inertial Confinement Fusion Implosions”, K. M. Woo, R. Betti, D. Shvarts, A. Bose, D. Patel, R. Yan, P.-Y. Chang, O. M. Mannion, R. Epstein, J. A. Delettrez, M. Charissis, K. S. Anderson, P. B. Radha, A. Shvydky, I. V. Igumenshchev, V. Gopalaswamy, A. R. Christopherson, J. Sanz and H. Aluie, Physics of Plasmas, volume 25, 052704 (2018) 



# Co-Author Publications



1. “Tripled yield in direct-drive laser fusion through statistical modelling”, V. Gopalaswamy, R. Betti, J. P. Knauer, N. Luciani, D. Patel, K. M. Woo, A. Bose, I. V. Igumenshchev, E. M. Campbell, K. S. Anderson, K. A. Bauer, M. J. Bonino, D. Cao, A. R. Christopherson, G. W. Collins, T. J. B. Collins, J. R. Davies, J. A. Delettrez, D. H. Edgell, R. Epstein, C. J. Forrest, D. H. Froula, V. Y. Glebov, V. N. Goncharov, D. R. Harding, S. X. Hu, D. W. Jacobs-Perkins, R. T. Janezic, J. H. Kelly, O. M. Mannion, A. Maximov, F. J. Marshall, D. T. Michel, S. Miller, S. F. B. Morse, J. Palastro, J. Peebles, P. B. Radha, S. P. Regan, S. Sampat, T. C. Sangster, A. B. Sefkow, W. Seka, R. C. Shah, W. T. Shmyada, A. Shvydky, C. Stoeckl, A. A. Solodov, W. Theobald, J. D. Zuegel, M. Gatu Johnson, R. D. Petrasso, C. K. Li,J. A. Frenje, Nature, volume 565, 581-586 (2019) 





2. “Theory of alpha heating in inertial fusion: Alpha-heating metrics and the onset of the burning-plasma regime”, A. R. Christopherson, R. Betti, J. Howard, K. M. Woo, A. Bose, E. M. Campbell and V. Gopalaswamy, Physics of Plasmas, volume 25, 072704 (2018) 





3. “A comprehensive alpha-heating model for inertial confinement fusion”, A. R. Christopherson, R. Betti, A. Bose, J. Howard, K. M. Woo, E. M. Campbell, J. Sanz, B. K. Spears, Physics of Plasmas, volume 25, 012703 (2018) 





4. “Analysis of trends in experimental observables: Reconstruction of the implosion dynamics and implications for fusion yield extrapolation for direct-drive cryogenic targets on OMEGA”,A. Bose, R. Betti, D. Mangino, K. M. Woo, D. Patel, A. R. Christopherson, V. Gopalaswamy, O. M. Mannion, S. P. Regan, V. N. Goncharov, D. H. Edgell, C. J. Forrest, J. A. Frenje, M. Gatu Johnson, V. Yu Glebov, I. V. Igumenshchev, J. P. Knauer, F. J. Marshall, P. B. Radha, R. Shah, C. Stoeckl, W. Theobald, T. C. Sangster, D. Shvarts and E. M. Campbell, Physics of Plasmas, volume 25, 062701 (2018) 





5. “The National Direct-Drive Program: OMEGA to the National Ignition Facility”, S. P. Regan and V. N. Goncharov and T. C. Sangster and E. M. Campbell and R. Betti and K. S. Anderson and T. Bernat and A. Bose and T. R. Boehly and M. J. Bonino and D. Cao and R. Chapman and T. J. B. Collins and R. S. Craxton and A. K. Davis and J. A. Delettrez and D. H. Edgell and R. Epstein and M. Farrell and C. J. Forrest and J. A. Frenje and D. H. Froula and M. Gatu Johnson and C. Gibson and V. Yu. Glebov and A. Greenwood and D. R. Harding and M. Hohenberger and S. X. Hu and H. Huang and J. Hund and I. V. Igumenshchev and D. W. Jacobs-Perkins and R. T. Janezic and M. Karasik and R. L. Keck and J. H. Kelly and T. J. Kessler and J. P. Knauer and T. Z. Kosc and S. J. Loucks and J. A. Marozas and F. J. Marshall and R. L. McCrory and P. W. McKenty and D. D. Meyerhofer and D. T. Michel and J. F. Myatt and S. P. Obenschain and R. D. Petrasso and N. Petta and P. B. Radha and M. J. Rosenberg and A. J. Schmitt and M. J. Schmitt and M. Schoff and W. Seka and W. T. Shmayda and M. J. Shoup III and A. Shvydky and A. A. Solodov and C. Stoeckl and W. Sweet and C. Taylor and R. Taylor and W. Theobald and J. Ulreich and M. D. Wittman and K. M. Woo and J. D. Zuegel, Fusion Science and Technology, volume 73, 89-97 (2018) 





6. “Electron Shock Ignition of Inertial Fusion Targets”, W. L. Shang, R. Betti, S. X. Hu, K. M. Woo, L. Hao, C. Ren, A. R. Christopherson, A. Bose and W. Theobald, Phys. Rev. Lett., volume 119, 195001 (2017) 





7. “The physics of long- and intermediate-wavelength asymmetries of the hot spot: Compression hydrodynamics and energetics”, A. Bose, R. Betti, D. Shvarts and K. M. Woo, Physics of Plasmas, volume 24, 102704 (2017) 





8. “Demonstration of Fuel Hot-Spot Pressure in Excess of 50 Gbar for Direct-Drive, Layered Deuterium-Tritium Implosions on OMEGA”, S. P. Regan, V. N. Goncharov, I. V. Igumenshchev, T. C. Sangster, R. Betti, A. Bose, T. R. Boehly, M. J. Bonino, E. M. Campbell, D. Cao, T. J. B. Collins, R. S. Craxton, A. K. Davis, J. A. Delettrez, D. H. Edgell, R. Epstein, C. J. Forrest, J. A. Frenje, D. H. Froula, M. Gatu Johnson, V. Yu. Glebov, D. R. Harding, M. Hohenberger, S. X. Hu, D. Jacobs-Perkins, R. Janezic, M. Karasik, R. L. Keck, J. H. Kelly, T. J. Kessler, J. P. Knauer, T. Z. Kosc, S. J. Loucks, J. A. Marozas, F. J. Marshall, R. L. McCrory, P. W. McKenty, D. D. Meyerhofer, D. T. Michel, J. F. Myatt, S. P. Obenschain, R. D. Petrasso, P. B. Radha, B. Rice, M. J. Rosenberg, A. J. Schmitt, M. J. Schmitt, W. Seka, W. T. Shmayda, M. J. Shoup, A. Shvydky, S. Skupsky, A. A. Solodov, C. Stoeckl, W. Theobald, J. Ulreich, M. D. Wittman, K. M. Woo, B. Yaakobi and J. D. Zuegel, Phys. Rev. Lett., volume 117, 025001 (2016) 





9. “Core conditions for alpha heating attained in direct-drive inertial confinement fusion”, A. Bose, K. M. Woo, R. Betti, E. M. Campbell, D. Mangino, A. R. Christopherson, R. L. McCrory, R. Nora, S. P. Regan, V. N. Goncharov, T. C. Sangster, C. J. Forrest, J. Frenje, M. Gatu Johnson, V. Yu Glebov, J. P. Knauer, F. J. Marshall, C. Stoeckl and W. Theobald, Phys. Rev. E, volume 94, 011201 (2016) 





10. “Alpha Heating and Burning Plasmas in Inertial Confinement Fusion”, R. Betti, A. R. Christopherson, B. K. Spears, R. Nora, A. Bose, J. Howard, Howard, K. M. Woo, Howard, M. J. Edwards, and J. Sanz, Phys. Rev. Lett., volume 114, 255003 (2015) 





11. “Hydrodynamic scaling of the deceleration-phase Rayleigh-Taylor instability”, A. Bose, K. M. Woo, R. Nora, and R. Betti, Physics of Plasmas, volume 22, 072702 (2015) 





12. “Theory of hydro-equivalent ignition for inertial fusion and its applications to OMEGA and the National Ignition Facility”, R. Nora, R. Betti, K. S. Anderson, A. Shvydky, A. Bose, K. M. Woo, A. R. Christopherson, J. A. Marozas, T. J. B. Collins, P. B. Radha, S. X. Hu, R. Epstein, F. J. Marshall, R. L. McCrory, T. C. Sangster, and D. D. Meyerhofer, Physics of Plasmas, volume 21, 056316 (2014) 



# First-Author Presentations



1. “Impact of Three-Dimensional Hot-Spot Flow Asymmetry on Ion-Temperature Measurements in Inertial Confinement Fusion Experiments”, K. M. Woo, R. Betti, D. Shvarts, O. M. Mannion, D. Patel, V. N. Goncharov, K. S. Anderson, P. B. Radha, J. P. Knauer, A. Bose, V. Gopalaswamy, A. R. Christopherson, E. M. Campbell, J. Sanz and H. Aluie, 60th Annual Meeting of the APS Division of Plasma Physics (2018) 





2. “Three-Dimensional Studies of the Effect of Residual Kinetic Energy on Yield Degradation”, K. M. Woo, R. Betti, A Bose, D. Patel and V. Gopalaswamy, 59th Annual Meeting of the APS Division of Plasma Physics (2017) 





3. “Study of Yield and Pressure Degradation in Inertial Confinement Fusion”, K. M. Woo, R. Betti, R. Yan, H. Aluie, A. Bose, D. X. Zhao and V. Gopalaswamy, 58th Annual Meeting of the APS Division of Plasma Physics (2016) 





4. “Three-Dimensional Simulations of the Deceleration Phase of Inertial Fusion Implosions”, K. M. Woo, R. Betti, A. Bose, R. Epstein, J. A. Delettrez, K. S. Anderson, R. Yan, P.-Y. Chang, D. Jonathan and M. Charissis, 57th Annual Meeting of the APS Division of Plasma Physics (2015) 





5. “The Three-Dimensional Hydrocode DEC3D with Multigroup Radiation Transport”, K. M. Woo, R. Epstein, J. A. Delettrez, A. Bose, R. Betti and K. S. Anderson, 56th Annual Meeting of the APS Division of Plasma Physics (2014) 



# Acknowledgments

I would like to thank my thesis advisor, Professor Riccardo Betti, for his support and insightful advice on my research studies, as well as his patience waiting for the final version of DEC3D code. I learned an important thing from Riccardo to use a strong physics intuition to solve problems. 

I would also like to thank my colleagues and collaborators from the Laboratory of Laser Energetics (LLE) including Dr. Po-Yu Chang, Dr. Rui Yan, Dr. Micheal Charissis, Dr. Reuben Epstein, Professor Dov Shvarts, Professor Javier Sanz Recio, Professor Hussein Aluie, and Owen Mannon for their helpful discussion in coding, parallel computation, numerical methods of radiation transport, theories of hydrodynamics instabilities, and experimental methods. 

I would also like to thank my group-mates Ryan, Dan, Arijit, Alison, Varchas, Dhrumir and Aarne for creating a fun research environment; all LLE colleagues in the illustration group for their huge effort to polish my APS presentation slides; all teachers from the University of Rochester including Professor Hagen for his elegant solution of symmetry to solve eigenvalues for a big matrix, Professor Das for his excellent memory during teaching, Laura for her coordination of my graduate study, and my friends in the University Badminton Club. 

I would also like to thank the dissertation committee for reviewing my thesis. 

Finally, I would also like to thank all my family members. 

# Abstract

Hydrodynamic instability is one of the primary sources of degrading the fusion yields in inertial confinement fusion (ICF) experiments. The presence of nonuniformities during the hot spot formation leads to dominant experimental signatures of implosion asymmetries. The physical mechanism of how hydrodynamic instabilities manifest themselves in experimental observable plays an important role to interpret three-dimensional (3-D) effects on ICF experimental data. 

In the first part of the thesis, we describe the development of a 3-D radiationhydrodynamic Eulerian spherical moving-mesh parallel code DEC3D to model the deceleration-phase Rayleigh-Taylor instability. The new code implements advanced modern numerical methods including the high-resolution shock-capturing technique the piecewise parabolic method for hydrodynamics, the macro-zoning technique to treat small time-step problems of the spherical mesh, and the integration of HYPRE to solve the implicit multi-group radiation diffusion. A singlemode and multi-mode simulation database was established to study the relations between 3-D hydrodynamic effects and implosion asymmetries. 

In the second part of the thesis, two comprehensive physical models were developed: (1) to explain the effects of the residual kinetic energy on the degradation of fusion yields and hot-spot pressures, and the property of larger hot-spot volumes for low modes, and (2) to explain the effects of 3-D hot-spot flow asymmetries on the variations of ion-temperature measurements. An analytical method of velocity variance decomposition was developed to infer the minimum ion temperatures and 

explain the physical mechanism of larger apparent ion temperatures than the true thermal ion temperatures. 

# Contributors and Funding Sources

This work was supervised by a dissertation committee consisting of Professor Riccardo Betti (advisor) of the Department of Mechanical Engineering, Professor Sarada G. Rajeev of the Department of Physics and Astronomy, Professor Adam Sefkow of the Department of Mechanical Engineering, and Professor Pierre-Alexandre Gourdain of the Department of Physics and Astronomy. Graduate study was supported by the Department of Physics and Astronomy at the University of Rochester, and a dissertation research Horton Fellowship. 

This material is based upon work supported by the Department of Energy National Nuclear Security Administration under Award Number DE-NA0001944 and DE-NA0003856, the University of Rochester, and the New York State Energy Research and Development Authority. Partial supported was provided by DOE Office of Fusion Energy Sciences grant DE-SC0014318. 

This report was prepared as an account of work sponsored by an agency of the U.S. Government. Neither the U.S. Government nor any agency thereof, nor any of their employees, makes any warranty, express or implied, or assumes any legal liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or process disclosed, or represents that its use would not infringe privately owned rights. Reference herein to any specific commercial product, process, or service by trade name, trademark, manufacturer, or otherwise does not necessarily constitute or imply its endorsement, recommendation, or favoring by the U.S. Government or any agency thereof. The views and 

opinions of authors expressed herein do not necessarily state or reflect those of the 

U.S. Government or any agency thereof. 

# List of Tables

4.1 DEC3D multi-mode perturbation (unit for $\sqrt { \sigma _ { i j } }$ is km/s) . . . . . . 118 

4.2 Performance of extrapolation thermal ion temperatures by inferring DD and DT along the same LOS at six different locations. Each LOS output DD and DT ion temperatures according to Brysk ion temperature formula in Eq. (4.1): $T _ { \mathrm { i } , X } ^ { \mathrm { i n f e r r e d } } = T _ { \mathrm { i } } ^ { \mathrm { t h e r m a l } } + ( m _ { \mathrm { n } } +$ $m _ { X } ) g _ { i } g _ { j } \sigma _ { i j } ^ { X }$ , where “X”denotes $\mathrm { H e ^ { 3 } }$ for DD reactions and $\alpha$ for DT reactions. Values of six hot-spot flow parameters $\sigma _ { i j } ^ { X }$ are shown in Table (4.1). In the method 3, thermal ion temperatures can be extrapolated directly from DD and DT ion temperature measured along the same LOS at one location given by $\hat { M } _ { 2 } ^ { - 1 } \cdot ( T _ { \mathrm { L O S } } ^ { \mathrm { D T } } , T _ { \mathrm { L O S } } ^ { \mathrm { D D } } ) ^ { T } \ .$ · $( 1 , 0 ) ^ { T }$ , where $T$ denotes for the transpose of a row vector into a column vector. In the method 2, the extrapolated thermal ion temperature is the averaged of six extrapolated thermal ion temperature from each LOS given by $\begin{array} { r } { \frac { 1 } { 6 } \hat { M } _ { 2 } ^ { - 1 } \cdot ( \vec { T } _ { 6 } ^ { \mathrm { D T } } \cdot \hat { e } _ { 6 } , \vec { T } _ { 6 } ^ { \mathrm { D D } } \cdot \hat { e } _ { 6 } ) ^ { T } \cdot ( 1 , 0 ) ^ { T } } \end{array}$ . In method 3, the extrapolated thermal ion temperature is shown closer to the burn or neutron-averaged thermal ion temperature. 

5.1 The summary for DEC2D and DEC3D. . . . . 132 

5.2 The summary for Cartesian-mesh DEC2D and DEC3D.[WBB+15] 142 

# List of Figures

1.1 DT fusion reactions in the laboratory scale and inside the core of the sun. The pressure and temperature of an ICF hot spot is comparable to the core of the sun. 2 

1.2 Laser energy deposition in direct-drive: (1) the hot electron preheat and (2) the cross-beam-energy-transfer degrade the implosion performance. 3 

1.3 The four stages in direct-drive ICF implosion: (1) shock launching and propagation, (2) acceleration phase, (3) deceleration phase and hot spot formation, and (4) disassembly phase and thermonuclear burn propagation. 4 

1.4 The adiabat shaping techniques [AB04] in (a) method of relaxed mass density profile, (b) method of decaying pressure. 5 

1.5 Comparison of classical and ablative RT instabilities. The mass ablation caused the thermal heat flux produces a dynamic pressure ${ \scriptstyle \frac { 1 } { 2 } } \rho _ { 1 } v _ { \mathrm { b } } ^ { 2 }$ , which reduces the pressure difference between the heavy and light fluids and the effective force acting on the perturbed heavy fluid mass. 6 

1.6 When non-stagnating hot-spot flow velocity is large, the neutroninferred ion-temperature measurements along different line-of-sight vary for the same imploding target. 8 

2.1 The initial mass-density profile for OMEGA shot 77068 at the beginning of deceleration phase $t _ { 0 } = 2 . 2 8$ ns. The initial radial velocity perturbation $\delta v _ { r } ( t _ { 0 } )$ is applied on the inner shell surface $r _ { 0 } = 6 8 \mu \mathrm { m }$ . The blue line shows the shape function $f ( r )$ that is unity at $r _ { 0 }$ . . . 13 

2.2 The summary of DEC3D single-mode simulation database for the shot 77068. (1) HLLC approximate Riemann solver with PPM highresolution method; (2) HYPRE implicit thermal diffusion; and (3) high angular resolution $1 2 8 \times 2 5 6$ for $\theta$ and $\phi$ zones. . . . 15 

2.3 3-D profiles for mass densities and vorticity magnitudes $| \omega |$ at stagnation: for 2-D single-mode $\ell = 8 , m = 0$ in (a) and (b); for 3-D single-mode $\ell = 8 , m = 4$ in (c) and (d), where $\vec { \omega } = \vec { \nabla } \times \vec { v }$ . . . . . . 19 

2.4 The 3-D hot-spot models for low and high modes. The red part is the burn volume while the blue part is the shell. As the shell implodes, the burn volume ${ \cal V } _ { \mathrm { b } } ( t )$ is shrinking with the boundary velocity $\vec { v } _ { \mathrm { b } } ( t )$ . For high modes, the cold bubbles shown by the white part do not contribute DT fusion reactions, and the actual burn volume is smaller than the total hot-spot volume bounded by the perturbed inner shell surface. . . . 22 

2.5 The 3-D hot-spot shapes for the electron-temperature contour surface (1, 1.5, 2, 2.5 keV) for high mode $Y _ { \ell = 1 2 } ^ { m = 6 }$ at stagnation. The length scale is the same in all figures. The left-most hot-spot shape includes the cold bubbles as shown in Fig. (3.1. It corresponds to the total perturbed hot-spot volume while the right-most figure is the hot core for DT fusion reactions. . . . 23 

2.6 Summary of single-mode energetics at stagnation for the same 14% initial velocity perturbations. Plot (a) is the burn-averaged hot spot and (b) is the $T _ { \mathrm { e } } \geq 1$ -keV hot spot. In (a) and (b), the blue dots are $\mathrm { K E _ { t o t a l } ^ { 3 D } } / \mathrm { K E _ { t o t a l } ^ { 1 D } }$ , black dots are $\mathrm { I E } _ { \mathrm { s h e l l } } ^ { 3 D } / \mathrm { I E } _ { \mathrm { s h e l l } } ^ { 1 D }$ and the red dots are $\mathrm { I E _ { h o t s p o t } ^ { 3 D } / I E _ { h o t s p o t } ^ { 1 D } }$ . The 3-D Legendre modes $\ell = 4 , 6 , 8 , 1 0 , 1 2$ with $m = \ell / 2$ are denoted by `.5 on the $x$ axis. (c) Measurement of $P _ { \mathrm { 1 k e V } } V _ { \mathrm { 1 k e V } } ^ { 5 / 3 }$ 5/3 with respect to 1-D at stagnation using the hot-spot definition with $T _ { \mathrm { e } } ~ \geq ~ 1$ keV shows the conservation of adiabatic parameter $P V ^ { 5 / 3 }$ for all modes. . . . . 

2.7 The blue curves indicate the analytic model relations for (a) $\hat { P } =$ $I \mathrm { \hat { E } } _ { \mathrm { H S } } ^ { 2 . 5 }$ , (b) $\hat { V } = \hat { I E } _ { \mathrm { H S } } ^ { - 1 . 5 }$ , (c) $\hat { P } = ( 1 - \mathrm { R K E } ) ^ { 2 . 5 }$ , and (d) $\hat { V } = ( 1 -$ RKE)−1.5. The shorthand notations are $\hat { Q } = Q _ { \mathrm { 3 D } } / Q _ { \mathrm { 1 D } }$ . . . . . . . . 28 

2.8 (a) The 1-D scaling relation for the hot-spot mass at stagnation $\hat { M } = \hat { P } ^ { 5 / 7 } \hat { V } ^ { 1 7 / 2 1 } \hat { \tau } ^ { 2 / 7 }$ for low modes $\ell = 1$ to 5 using the hot-spot definition $T _ { \mathrm { e } } \geq 1$ keV. The shorthand notations are $\hat { Q } _ { \mathrm { 3 D } } = Q _ { \mathrm { 3 D } } / Q _ { \mathrm { 1 D } }$ . (b) Comparison of simulated YOC against RKE and the analytic models for low modes $\ell = 1$ to 6 and (c) for high modes $\ell = 7$ to 12. The solid blue curve is the yield degradation model YOC = $( 1 - \mathrm { R K E } ) ^ { 4 . 4 }$ , while the dashed blue curve is yield degradation model $\mathrm { Y O C } = ( 1 - \mathrm { R K E } ) ^ { 5 . 5 }$ . . 33 

2.9 (a) For low mode $\ell = 1$ , the degradation of hot-spot pressure and the increasing hot-spot volume are strong functions of residual kinetic energies. (b) For the high mode $Y _ { \ell = 1 2 } ^ { m = 6 }$ , the burn volume is reduced and the core pressure is increased due to the growth of converging RT spikes. The burn-averaged quantities are used in (a) and (b). (c) The simulated YOC is compared against the hot-spot volume for high modes $\ell = 7$ to 12. The blue curve is the yield degradation model $\mathrm { Y O C } = \hat { V }$ , where $\hat { Q } = Q _ { \mathrm { 3 D } } / Q _ { \mathrm { 1 D } }$ . . . . . . 34 

3.1 The mass density profile for (a) low mode $\ell = 1$ and (b) high mode $Y _ { \ell = 1 2 } ^ { m = 6 }$ at stagnation with electron temperature contours at 1, 2, 3 and 4 keV. Black arrows indicate the fluid velocity field. Since the cores at temperatures $T _ { \mathrm { e } } ~ \geq ~ 4$ keV are approximately spherical in shape, the core compression is contributed mainly by the radial velocity component. For high mode $Y _ { \ell = 1 2 } ^ { m = 6 }$ , the bubbles are cold, characterized by low neutron production rates. . . 

3.2 The kinetic energy density profiles for (a) low mode $\ell = 1$ and (b) high mode $Y _ { \ell = 1 2 } ^ { m = 6 }$ at stagnation. For low mode $\ell = 1$ , both the RT spike and the jet contribute to residual kinetic energy. For high mode $Y _ { \ell = 1 2 } ^ { m = 6 }$ , the residual kinetic energy is dominated mainly by the nonstagnating RT spikes. The blue spherical outline is the residual kinetic energy of the unshocked part of the shell. . . . . . 

3.3 Velocity magnitude for low mode $\ell = 1$ at stagnation. A jet at 500 km/s is flowing along the negative $z$ direction. The yield-over-clean is 0.74 and the burn-averaged hot-spot velocity magnitude is 166 km/s in this simulation. The red curve indicates the hot core at a 2.3-keV electron temperature, showing that the hot core is advected by the jet with the maximum flow velocity. . . 

3.4 (a) The inferred DT ion temperatures at stagnation by IRIS3D using six detectors at different LOS’s at $- X , + X , - Y , + Y , - Z , + Z$ for low mode $\ell = 1$ (red curve) and high mode $Y _ { \ell = 1 2 } ^ { m = 6 }$ (blue curve). (b) The inferred DT ion temperatures at stagnation by IRIS3D using 16 detectors at different LOS’s including $- X , + X , - Y , + Y , - Z , + Z$ and other 10 typical nTOF (neutron time of flight) diagnostics on OMEGA. . . . 

3.5 (a) DEC3D mass density profile for the single-mode $Y _ { \ell = 1 } ^ { m = 0 }$ at stagnation, simulated by 7% initial velocity perturbation. The electron temperature of 2.4 to 2.45 keV is shown by the red hemisphere located at the origin. The arrows indicate the fluid velocity vectors. A jet with neutron-averaged velocity $\left. v _ { z } \right. = - 1 7 6 ~ \mathrm { k m / s }$ is shown flowing through the central part of the hot spot toward the negative $z$ direction. $\theta$ is the angle measured in the clockwise direction between the positive $z$ axis and LOS. (b) Comparison of the inferred ion temperature measurements from the neutron energy spectrum model in Eq. (3.8) with IRIS3D[WRF18] by post-processing DEC3D mode $\ell = 1$ hydrodynamic data in (a). . . 

3.6 Comparison of neutron-inferred ion temperatures between the synthetic neutron energy spectrum model by Eq. (3.8) and IRIS3D by Monte-Carlo simulations of neutron transport. . . . . 48 

3.7 DEC3D kinetic energy density profile in (a) and DEC3D mass density profile in (b) at stagnation for high mode $\ell = 4 0 , m = 2 0$ . The green contour surface is $T _ { \mathrm { e } } = 1$ keV. The vortices of high modes are localized within the cold bubbles, $T _ { \mathrm { e } } < 1$ keV, that do not produce significant amount of neutrons. The Doppler shift term resulting from vortices of large- $\ell$ single-mode perturbations is negligible because of low contribution to the burn distribution. 49 

3.8 Result of neutron energy spectrum model by post-processing DEC3D hydrodynamic data to compare the Doppler shift of the mean neutron energy and the bulk velocity broadening in (a) for low mode $\ell = 1$ and (b) for high mode $\ell = 4 0 , m = 2 0$ . The unnormalized Doppler-shifted neutron energy spectra sampling inside the cold bubble $f ( E _ { \mathrm { n } } ) _ { \mathrm { b u b b l e } }$ and the high-temperature hot core $f ( E _ { \mathrm { n } } ) _ { \mathrm { c o r e } }$ are compared with the normalized neutron energy spectrum $f _ { \mathrm { L O S } } ( E _ { \mathrm { n } } ) e$ observed at $+ z$ . 

3.9 Measurement of (a) $R _ { T }$ by the neutron energy spectrum model and IRIS3D. (b) ratios of maximum and minimum neutron-inferred ion temperatures to thermal ion temperatures by the neutron energy spectrum model. 3-D modes $Y _ { \ell = 4 } ^ { m = 2 }$ , $Y _ { \ell = 6 } ^ { m = 3 }$ Y , $Y _ { \ell = 8 } ^ { m = 4 }$ , $Y _ { \ell = 1 0 } ^ { m = 5 }$ and $Y _ { \ell = 1 2 } ^ { m = 6 }$ are denoted by $\ell = 4 . 5$ , 6.5, 8.5, 10.5 and 12.5 on the $x$ -axis. . . . 

3.10 DEC3D low mode $\ell = 2$ simulation by 7% initial velocity perturbation. The velocity fields are plotted on the top of mass density profiles on the $x - y$ and $z - x$ planes. The middle is the 3-D electron temperature contour surface to visualize the 3-D configuration of the high-temperature hot core. Non-stagnating hot-spot fluid velocity disturbance driven by the pair of RT spikes along the $z$ -axis and the expanding radial flow structure inside the large donut-shape warm bubble contribute a significant amount of non-translational residual kinetic energies. . . . 

3.11 (a)The synthetic neutron energy spectrum for mode $\ell = 2$ simulated by 14% initial velocity perturbation for a LOS located on the equator along the positive $x$ -axis. Positive excess kurtosis are pbserved around the tails of the neutron energy spectrum. (b)The synthetic neutron energy spectrum for mode $\ell = 1$ simulated by 7% initial velocity perturbation for a LOS located on the north pole. The downward-flowing jet leads to the formation of negative skewness and negative kurtosis relative to an observer located at the north pole. . . . . 

3.12 Measurement of kurtosis $M _ { 4 }$ for the single-mode spectrum $\ell = 1 - 1 2$ simulated by 1% − 14% initial velocity perturbations. For a normal Gaussian distribution, the kurtosis is 3, while the excess kurtosis is zero which is defined by subtracting 3 from kurtosis. . . . . . 

3.13 Comparison of inferred ion temperature ratio $R _ { T }$ between the Brysk temperatures by Eq. (3.30) through measuring the velocity variance and the synthetic neutron energy spectrum through (a) inferring ion temperatures from the FWHM of the neutron energy spectrum $f _ { \mathrm { L O S } } ( E _ { \mathrm { n } } )$ and (b) inferring ion temperatures from the velocity variance of the neutron velocity spectrum $f _ { L O S } ( v _ { \mathrm { n } } )$ . . . . . 

3.14 Comparison of maximum to minimum neutron-inferred ion temperatures between IRIS3D and the Brysk ion temperatures to validate Eq. (3.44) using DEC3D single-mode database $\ell = 1 - 1 2$ with different initial velocity perturbations $\delta v / v _ { 0 } = 0 . 0 1 - 0 . 1 4$ . . . . . 

3.15 A numerical test to show the term $\hat { \delta } _ { \mathrm { p } , 0 } \cdot \vec { T } _ { \mathrm { t h } }$ has vanishing contribution in $4 \pi$ solid angles by taking the thermal ion temperature as 3 keV. Six LOS’s are indicated by red dots. The coordinates $( \theta , \phi )$ in unit of radians are $( 1 . 0 7 0 0 , 0 . 8 3 1 4 ) _ { 1 }$ , $( 1 . 0 8 4 7 , 3 . 5 8 8 8 ) _ { 2 }$ , (2.0345, 2.8274)3, (0.6705, 4.3563)4, $( 1 . 5 3 3 4 , 2 . 8 1 4 2 ) _ { 5 }$ , and (1.4831, 5.4412)6 

54 

55 

3.16 (a) Prediction of neutron-inferred ion temperatures $T _ { \mathrm { p } }$ in the full map by Eq. (3.55) using six neutron-inferred ion temperatures, shown by the red dots, at six nTOF locations in OMEGA. The neutron-averaged thermal ion temperature $\langle T _ { \mathrm { i } } \rangle _ { \mathrm { b } } \equiv T _ { \mathrm { i } } ^ { \mathrm { t h e r m a l } }$ is 3.551 keV, and the minimum of neutron-inferred ion temperature $T _ { \mathrm { i , p r e d } } ^ { \mathrm { t h e r m a l } } \equiv$ T appro.infi,min is 3.525 keV. (b) Comparison of neutron-inferred ion tem- $T _ { \mathrm { i , m i n } } ^ { \mathrm { a p p r o . i n f } }$ peratures simulated by IRIS3D and the neutron energy spectrum model using 16 LOS’s with that of simulated by the Brysk ion temperatures by Eq. (3.55) using six LOS’s. . . . 70 

3.17 DEC3D single-mode simulations with $\delta v / v _ { 0 } = 0 . 0 7$ initial velocity perturbations. (a) Comparison of the directional and total velocity variances. Mode $\ell = 2$ exhibits the largest total velocity variance. (b) Comparison of the minimum inferred ion-temperature formula in Eq. (3.56) shown by open and solid red circles with the improved approximation formula in Eq. (3.66) shown by open and solid blue squares with isotropic variance separation. Significant improved prediction of thermal ion temperatures is observed. . . . . . 71 

3.18 The full map of ion-temperature measurement variation for low mode $\ell = 4 , m = 2$ by Eq. (3.45) through direct computation of variance and covariance from DEC3D hydrodynamic data. Two dim regions are resulted from unequal directional variance $\sigma _ { x x } \neq \sigma _ { y y }$ in Eq. (3.59). . . . . 72 

3.19 Investigation of 3-D effects of covariance terms for the multimode perturbation $\ell = 1 0 , m = 5$ , and $\ell = 1$ by post-processing DEC3D hydrodynamic data to obtain directional variance and covariance: (a) use Eq. (3.44) to obtain $\begin{array} { r } { \hat { T } _ { i } ^ { \mathrm { { i n f e r r e d } } } = \hat { T } _ { i } ^ { \mathrm { { t h e r m a l } } } + \sum _ { i = 1 } ^ { 3 } \sigma _ { i i } g _ { i } g _ { i } + } \end{array}$ $\textstyle \sum _ { i \neq j } \sigma _ { i j } g _ { i } g _ { j }$ ; (b) use the same formula but neglects directional variance terms to obtain $\begin{array} { r } { \hat { T } _ { i } ^ { \mathrm { i n f e r r e d } } = \hat { T } _ { i } ^ { \mathrm { t h e r m a l } } + \sum _ { i \neq j } \sigma _ { i j } g _ { i } g _ { j } } \end{array}$ . . . . 75 

3.20 Plot (a) compares $R _ { T }$ in DEC3D two-mode simulations (1) with a dominant low mode $\ell = 1$ in the blue curve and (2) with a dominant high mode $\ell = 1 0 , m = 5$ in the red curve. Plot (b) compares the corresponding yield-over-clean. . . . 

3.21 (a) The neutron-inferred ion-temperature sky map simulated by IRIS3D at high resolution by post-processing a strongly distorted DEC3D multimode simulation with an initial perturbation spectrum $\begin{array} { r l } { \Delta v / v _ { 0 } \sum _ { \ell = 1 } ^ { 1 2 } Y _ { \ell } ^ { m = \ell / 2 } } \end{array}$ and YOC = 0.36. (b) DEC3D mass density and velocity field profiles on the $x - y$ plane. A developed jet structure is observed in the $x$ direction. The shape of the distorted hot spot is indicated by $T _ { \mathrm { e } } = 0 . 5$ -keV contour surface. . 

3.22 Reconstruction of 3-D neutron-inferred ion-temperature profiles by Brysk ion temperature model for the strongly perturbed multi-mode simulation shown in Fig. (3.21) (a) using Eq. (3.44) by including three directional variance terms denoted by $\mathrm { v a r } = \sigma _ { x x } g _ { x } g _ { x } +$ $\sigma _ { y y } g _ { y } g _ { y } + \sigma _ { z z } g _ { z } g _ { z }$ and three covariance terms denoted by $\begin{array} { r l } { \mathrm { c o v } } & { { } = } \end{array}$ $2 \sigma _ { x y } g _ { x } g _ { y } + 2 \sigma _ { y z } g _ { y } g _ { z } + 2 \sigma _ { z x } g _ { z } g _ { x }$ , (b) including only the three directional variance terms, and (c) including only the three covariance terms. The neutron-averaged thermal ion temperature is $\langle T _ { \mathrm { i } } \rangle _ { \mathrm { b } } ^ { I R I S ^ { g } D }$ = 2.7 keV. . . . 

3.23 DEC3D 77068 single-mode database. Plot (a) validates Eq. (3.73) 

by comparing the sum of three neutron-inferred ion temperatures 

at $\hat { x }$ , $\hat { y }$ , and $\hat { z }$ to the thermal ion temperature and total variance. 

The red dots are low mode $\ell = 1$ . Plot (b) compares the minimum 

neutron-inferred ion temperatures defined by Eq. (3.56) with the 

neutron-averaged thermal ion temperatures for all single modes. 

The red dots are low mode $\ell = 1$ and the blue dots are low mode $\ell =$ 

2. Plot (c) compares the average neutron-inferred ion temperatures 

over three orthogonal directions with hot-spot and shell residual 

kinetic energies to validate Eq. (3.76). 82 

4.1 Comparison of single-mode burn-averaged linear velocities in (a) 

and bilinear velocities in (b) at stagnations simulated by DEC3D 

at 7% initial velocity perturbation for OMEGA shot 77068. Mode 

$\ell = 1$ has large burn-averaged $z$ -velocity due to the jet. Burn-

averaged bilinear velocities in (b) are negligible. The combined 

results in (a) and (b) lead to vanishing single-mode covariance. . . 87 

4.2 A level diagram to represent the transition of levels of the absolute 

magnitudes of covariance from single modes to multi-modes, and 

finally to a fully-developed turbulence by random mixing of multi-

modes. 90 

4.3 Comparison of time-integrated burn-averaged directional-variance 

$\langle \sigma _ { x x } \rangle$ and $\langle \sigma _ { z z } \rangle$ against different levels of initial velocity perturba-

tions. The transition of $\left. \sigma _ { z z } \right. \mathrm { ~ < ~ } \left. \sigma _ { x x } \right.$ is observed in large mode 

$\ell = 2$ perturbations with $\triangle v / v _ { 0 } > 1 0 \%$ in $D E C 3 D$ simulations. . . 93 

4.4 Comparison of the ratio $\sigma _ { z z } ^ { \mathrm { D T } } / \sigma _ { x x } ^ { \mathrm { D T } }$ for other modes $\ell = 3 - 1 2$ 

investigate the transition phenomenon. . . . 95 

4.5 (a) Comparison of 1-D profiles of normalized total ion number density and normalized thermal ion temperature for OMEGA shot 77068 at stagnation. The blue solid line indicates the range of influence caused short-wavelength mode perturbations from within the region $\mathcal { R } _ { \mathrm { e d g e } } : r _ { c } < r < R _ { \mathrm { h s } }$ due to a higher burn-weight factor $Y _ { \mathrm { c e l l } } ^ { \mathrm { D D } } / Y _ { \mathrm { t o t a l } } ^ { \mathrm { D D } } > Y _ { \mathrm { c e l l } } ^ { \mathrm { D  T } } / Y _ { \mathrm { t o t a l } } ^ { \mathrm { D ^ { \prime } I } }$ √ . The blue dashed line indicates the range of influence caused by long-wavelength mode perturbations within the region $\mathcal { R } _ { \mathrm { c o r e } } : r < r _ { c }$ due to a higher burn-weight factor $Y _ { \mathrm { c e l l } } ^ { \mathrm { D 1 T } } / Y _ { \mathrm { t o t a l } } ^ { \mathrm { D 1 } } > Y _ { \mathrm { c e l l } } ^ { \mathrm { D 1 D } } / Y _ { \mathrm { t o t a l } } ^ { \mathrm { D 1 D } }$ . (b) Comparison of DD and DT burnweight factor defined by $Y _ { \mathrm { c e l l } } ( r ) / Y _ { \mathrm { t o t a l } }$ as function of radius. Arbitrary unit is used in the $y$ -axis to show the difference of burn-weight factors essentially. . . 

4.6 Comparison of fluid velocity field vectors on $x - z$ planes for low modes $\ell = 1 - 2$ and a high mode $\ell = 1 2 , m = 6$ . Contours of electron temperatures at are shown to outline the region of high burn-weights. All plots are in the same spatial scale shown by the rulers of 10µm in the corner of mode $\ell = 1$ plot. The color contour shows the spatial profile of $n _ { \mathrm { i } } ^ { 2 } T _ { \mathrm { i } } ^ { \beta } r ^ { 2 }$ , where $\beta = 3 . 9$ is taken to study the spatial distribution of DT burn-weight factor. The increasing burn weights are shown by bright regions closer to hot spot interface within warm bubbles only. . . 

4.7 Comparison of the ratio $T _ { \mathrm { m i n } } ^ { \mathrm { L D D } } / T _ { \mathrm { m i n } } ^ { \mathrm { L D } }$ against the fraction of residual kinetic energy $f _ { \mathrm { r k e } }$ for the mode spectrum $\ell = 1 - 1 2$ over all levels of 1%-14% initial velocity perturbations. The back solid line is Eq. (4.39) by taking the reaction product mass ratio $R _ { m } ^ { \mathrm { { D D / D T } } } = 0 . 8$ and isotropic velocity variance ratio $R _ { \sigma } ^ { \mathrm { D D / D T } } = 1$ . The blue arrows show the trend of decreasing DD/DT minimum inferred ion-temperature ratio as a result of increasing isotropic velocity variance. . . . 

99 

4.8 The full map of three different covariance terms to understand the distributions of maximum and minimum values at all LOS angles $\theta$ and $\phi$ . (a) shows $x$ - $y$ covariance term $\mathrm { c o v _ { 1 2 } } / \sigma _ { 1 2 } = \sin ^ { 2 } \theta \sin 2 \phi$ , (b) shows $y$ - $z$ covariance term $\mathrm { c o v } _ { 2 3 } / \sigma _ { 2 3 } = \sin 2 \theta \sin \phi$ , (c) shows $z$ - $x$ covariance term $\cos _ { 3 1 } / \sigma _ { 3 1 } = \sin 2 \theta \cos \phi$ . (d) overlaps all contour lines with values $\pm 0 . 7 5$ to examine the influence of superposition effect. . . . . 

4.9 Comparison of DD/DT minimum inferred ion-temperature ratio with the fraction of residual kinetic energy in multi-mode perturbations. The initial spectrum $\begin{array} { r } { A _ { \ell } ^ { m } = \sum _ { i = 1 } ^ { N } ( \triangle v / v _ { 0 } ) Y _ { \ell } ^ { m } ( \theta + \theta _ { i } , \phi + \phi _ { i } ) } \end{array}$ is obtained by superposition of a given single mode with $N$ -set of random phases, where $N = 2 0$ was used. The black solid-dashed line is the analytic curve by Eq. (4.39) with product mass ratio $R _ { m } ^ { \mathrm { { D D / D T } } } ~ = ~ 0 . 8$ and DD/DT isotropic velocity variance ratio $R _ { \sigma } ^ { \mathrm { D D / D T } } = 1$ . The blue solid-dashed line is the same analytic curve by substituting Murphy’s definition of fraction of residual kinetic energy $f _ { \mathrm { r e k } } = 4 f _ { \mathrm { r e k } } ^ { \mathrm { M } }$ Single mode perturbations only provide a small range of $f _ { \mathrm { r e k } } ^ { \mathrm { s i n g l e - m o d e } } : 0 - 0 . 3$ corresponding to a weak degradation of DD/DT minimum ion-temperature ratios to below unity. The mechanism to produce large isotropic source require multi-mode perturbations to fill in the hot spot with numerous isotropic flows while significantly degrade the thermal ion temperature to push $f _ { \mathrm { r e k } } = ( m _ { \mathrm { n } } + m _ { \alpha } ) \sigma _ { \mathrm { i s o } } ^ { \mathrm { l y T } } / T _ { \mathrm { i } } ^ { \mathrm { t h e r m a l } }$ to transit from single-mode regime f single−moderek into multi-mode regime f mulrek Jrek $f _ { \mathrm { r e k } } ^ { \mathrm { s i n g l e - m o d e } }$ ti−mode : 0.3 − 1. The legend $f _ { \mathrm { r e k } } ^ { \mathrm { m u l t i - m o d e } } : 0 . 3 - 1$ of different color points is the same as in Fig. (4.7). The purple and green points that result in the least ratio $T _ { \mathrm { m i n } } ^ { \mathrm { D D } } / T _ { \mathrm { m i n } } ^ { \mathrm { D T } } \sim 0 . 9$ are single-mode random-phase simulations for mode $Y _ { \ell = 4 } ^ { m = 2 }$ and $Y _ { \ell = 6 } ^ { m = 3 }$ respectively, meaning that their flow structure are highly isotropic under random phase mixings. Data lies on the black curve implying the accuracy of analytic formula by Eq. (4.39). 

4.10 Comparison of isotropic velocity variance formula in Eq. (4.24) in the single-mode random-phase simulations between DD and DT. $\sigma _ { \mathrm { i s o } } ^ { \mathrm { D D } }$ and $\sigma _ { \mathrm { i s o } } ^ { \mathrm { D T } }$ isotropic velocity variance are approximately equal at the leading order. Small deviations are caused by the shifts in positions of peak burn-weight between DD and DT. 

106 

4.11 The search for the $6 t h$ LOS angles $\theta _ { 6 }$ and $\phi _ { 6 }$ in NIF is obtained by discretizing $\theta$ and $\phi$ angles into 16 and 32 uniform mesh. Current NIF five LOS are NITOF at ( $\theta _ { 1 } = 9 0 , \phi _ { 1 } = 3 1 5$ ) for DT, Spec-A at ( $\theta _ { 2 } = 1 1 6 , \phi _ { 2 } = 3 1 6$ ) for DT and DD, Spec-SP at ( $\theta _ { 3 } = 1 6 1 , \phi _ { 3 } =$ 56) for DT and DD, MRS at ( $\theta _ { 4 } = 7 3 , \phi _ { 4 } = 3 2 4$ ) for DT, and Spec-E at ( $\theta _ { 5 } = 9 0 , \phi _ { 5 } = 1 7 4$ ) for DT and DD. The peak-to-valley are defined by $\hat { \delta } _ { \mathrm { m a x } } = \mathrm { M a x } [ \hat { \delta } ( \theta , \phi ) ] _ { 4 \pi }$ and $\hat { \delta } _ { \mathrm { m i n } } = \mathrm { M i n } [ \hat { \delta } ( \theta , \phi ) ] _ { 4 \pi }$ , where the departure matrix is given by $\hat { \delta } = \hat { I } - \hat { M } _ { \mathrm { n e w } } \cdot \hat { M } _ { \mathrm { L O S } } ^ { - 1 }$ . . . . . . . . 110 

4.12 (a) The map of standard deviation $\mathrm { v a r } [ \vec { T } _ { \mathrm { n e w } } ( \theta , \phi ) ] ^ { 1 / 2 }$ for OMEGA six LOS. (b) The positions of re-allocated six LOS to minimize error propagations. . . . . 111 

4.13 (a) Implementation of 6-LOS method in Eq. (4.61) to extrapolate the true minimum DT inferred ion temperatures in OMEGA experiments according to method 1. The red data are minimum DT ion temperatures over the purple regions in Fig. (4.12-a) that has error propagation less than 0.2 keV. The gray circles are the minimum of 6 nTOF experimental DT ion temperatures. The red data has less spread and shows increasing temperature dependence of experimental yield. (b) Comparison of fitting exponent $b _ { \mathrm { { f i t } } }$ against different tolerance of error propagation. The first data at error of 0.2 keV gives the most robust fitting exponent $b _ { \mathrm { f i t } } = 3 . 4 8$ . . . . . . . 113 

4.14 Investigation of the performance of Eq. (4.72) in extrapolating the minimum DD ion temperature in OMEGA experiments. Gray circles are DD ion temperatures measured by 13.4 m nTOF while red circles are extrapolated DD minimum ion temperatures by Eq. (4.72). The minimum DT ion temperatures are taken as: (a) the minimum among all available nTOF’s DT measurements i.e. $T _ { \mathrm { m i n } } ^ { \mathrm { D T } } = \mathrm { M i n } [ T ^ { \mathrm { D T } } ]$ , and (b) the extrapolated minimum DT ion temperature from the 6-LOS method by searching the minimum over the purple region in Fig. (4.12-a). The sample size is the same as in Fig. (4.13). . . . . . . 116 

4.15 DEC3D multi-mode simulation for a strongly distorted hot-spot. Left is the 3-D electron temperature contour surface at 0.8 keV. DD and DT ion temperatures are inferred along the same LOS at six different locations. Right is the 2-D $x$ - $z$ plane for the hot spot electron temperature at stagnation. Black arrows are hot spot fluid velocity vectors. The size of arrow heads increase with the magnitude of fluid velocities. The red contour line is the electron temperature at 0.55 keV. . . . . . . . 117 

4.16 The same DEC3D multi-mode simulation as described in Fig. (4.15) for Brysk ion temperature in (a) and the neutron-inferred hot-spot flow velocities in (b). A strong correlation is observed between the flow velocity vector of the jet and positions of maximum neutroninferred ion temperatures. . . . . . 117 

4.17 (a) Comparison of relative changes between the second and the third term in Eq. (4.76). (b) Comparison of relative changes between the first and the second term in Eq. (4.77). . . . . . 121 



4.18 Comparison of DD to DT ion-temperature ratios $T _ { \mathrm { D D } } ^ { 1 3 . 4 \mathrm { m } } / T _ { \mathrm { D T } } ^ { \mathrm { P e t a l } }$ measured along the same LOS in OMEGA experiments, indicated by black circles, with the same multi-mode simulation described in Figs. (4.17-a, b and c), indicated by red circles. . . . . 122 





4.19 (c) The full map variation of DD/DT neutron-inferred ion temperature ratios given by Eq. (4.76). General features are strongly correlated with DT neutron-inferred ion temperature and neutroninferred hot-spot flow velocity asymmetries. . 123 





4.20 Investigation the fluid properties of mode $\ell = 1$ : (a) examine the small isotropic velocity variance in Eq. (4.85) (b) examine the fraction of the non-translational hot-spot residual kinetic energy with respect to the total hot-spot residual kinetic energy in Eq. (4.89) . 125 





4.21 Examination of Eq. (4.92) where $\begin{array} { r } { \xi = \frac { 1 } { 4 } ( R _ { T } - 1 ) } \end{array}$ and $R _ { T } = T _ { \mathrm { m a x } } ^ { \mathrm { I R I S 3 D } } / T _ { \mathrm { m i n } } ^ { \mathrm { I R I S 3 D } }$ . For mode $\ell = 1$ , $T _ { \mathrm { m a x } } ^ { \mathrm { l R l S 3 D } }$ are obsvered at the north and south poles while $T _ { \mathrm { m i n } } ^ { \mathrm { l K I S 3 D } }$ is observed at the equator. . . . . 128 





4.22 A correlation study between OMEGA experimental yields and DD/DT neutron-inferred ion-temperature measurement asymmetries. The impact of hot-spot flow isotropy represented by the term $T _ { \mathrm { m i n } } ^ { \mathrm { L D D } } / T _ { \mathrm { m i n } } ^ { \mathrm { D T } }$ is a function of non-translational hot-spot residual kinetic energies given by Eq. (4.39). The impact of hot-spot flow anisotropy represented by the term $1 - \xi$ is a function of the shell residual kinetic energies for low modes. The impact of different fusion reactivities between DD and DT is represented by the term $1 - T _ { \mathrm { D D } } ^ { \mathrm { L O S ( 1 3 . 4 m ) } } / T _ { \mathrm { D T } } ^ { \mathrm { L O S ( p e t a l ) } }$ T LOS(13.4m)DD /T LOS(petal)DT . GLILAC10 is the mode 10 growth factor $G _ { 1 0 } ^ { \mathrm { L I L A C } }$ obtained from 1-D $L I L A C$ simulations. 129 





5.1 The 3-D electron temperature contour at 1 keV at stagnation for the single-mode $Y _ { \ell = 1 2 } ^ { m = 6 }$ . The resolution is $1 2 8 \times 1 2 8 \times 2 5 6$ for $r , \theta , \phi$ zones, respectively. . . . . . . . . 132 





5.2 The benchmark test for a full 3-D, deceleration-phase, clean simulation for the shot 77068. Comparison of solutions between LILAC (black circles) and DEC3D (solid red curves) for (a) density, (b) total pressure, (c) radial velocity, and (d) electron temperature profiles at stagnation using resolution 128 $\times$ 64×128 in $r , \theta , \phi$ directions. The 1-D burn radius at stagnation for shot 77068 is ∼ 19 µm, and the return shock is located at ∼ 28 $\mu$ m. . . . 133 





5.3 The benchmark test for neutron productions without alpha heating. Comparison of solutions between $L I L A C$ (black circles) and DEC3D (solid red curves) for the temporal history of neutron rate. The fusion reactivity $< \sigma v >$ scales with ion temperatures $T _ { \mathrm { i } }$ in a power law $\sim T _ { \mathrm { i } } ^ { 3 . 8 5 }$ for the temperature range $0 . 2 \ : < \ : T _ { \mathrm { i } } \ : < \ : 5$ keV in the BUCKY [HMS05] fusion reactivity model. . . . . . . . . 134 





5.4 The flow chart for executing different physical module in one timestep in DEC3D. 136 





5.5 The Cartesian topology domain decomposition [CLS94] for the Cartesian mesh version DEC3D [WBB $^ +$ 15] through the message-passinginterface (MPI) for a 3-D mass density profile. . . . 140 





5.6 Exchanging boundary information between neighboring sub-domains through MPI send and receive functions in one time-step. . . . . . 141 





5.7 Three ghost cells are required at each boundary of the radial subdomain for the third-order explicit PPM hydro update. . . . 143 





5.8 The radial domain decomposition of $D E C 3 D$ into six sub-domains for a 3-D mass density profile for a single-mode $\ell = 6$ simulation. . 144 





5.9 Example of coarse mesh generation of azimuthal mesh along the poles in DEC3D. 145 





5.10 Illustration of 3-D coarse mesh generation near the origin and along the poles in different radial subdomains. 146 





5.11 The principle of 1-D finite-volume moving-mesh method. The cell interfaces are moving to the left in the next time step level $t ^ { n + 1 }$ so that the cell-averaged quantity $Q _ { i } ^ { n + 1 }$ is defined by the linear interpolation to account for the contribution of $Q _ { i - 1 } ^ { n }$ and $Q _ { i } ^ { n }$ in the previous time-step level $t ^ { n }$ . . . . 150 





5.12 van Leer’s geometrical interpolation for flux limiters. By defining the slopes $\triangle Q _ { i + 1 / 2 } = Q _ { i + 1 } - Q _ { i }$ and $\triangle Q _ { i - 1 / 2 } = Q _ { i } - Q _ { i - 1 }$ , a cell is considered as a numerical noise if the product of slopes is negative $\triangle Q _ { i + 1 / 2 } \triangle Q _ { i - 1 / 2 } < 0$ , and no cell reconstruction is needed $\triangle Q _ { i } = 0$ . Otherwise, the cell profile is monotonic and is reconstructed by addling a linear slope $\triangle Q _ { i } = \left( \triangle Q _ { i + 1 / 2 } + \triangle Q _ { i - 1 / 2 } \right) / 2$ . . . . . . . . 157 





5.13 Significant numerical noises are observed in MacCormack scheme across the sharp shell interface, whereas the solution of MUSCL scheme HLLC approximate Riemann solver is numerically stable. . 159 





5.14 Comparison of 1-D LeVeque wave-propagation algorithm between a slow moving-mesh at $\dot { x } = - 1$ and a fast moving-mesh at ${ \dot { x } } =$ $- 2$ that updates a operator-split moving mesh using Eq. (5.73), followed by a hydro update using the second-order MUSCL-HLLC Riemann solver in a shock tube problem. . . . 161 





5.15 The upwind solution in the first-order Lax-Wendroff update is equivalent to a linear interpolation to weigh all possible cell contents leaving and entering into a cell. . . . . . . . . . 162 





5.16 Performance of 3-D parallel simulations for modes $\ell = 1 0$ and $\ell = 2 0$ using the Cartesian-mesh version DEC3D implemented with the Cartesian moving-mesh algorithm described by Eq. (5.97). . . . . . 168 





5.17 The implementation of moving-mesh in Cartesian and spherical geometries in DEC2D and DEC3D. Leibniz integral rule states that the total rate of change of a cell content $\begin{array} { r } { \frac { d } { d t } \int _ { \mathcal { D } ( t ) } Q ( \vec { x } , t ) d V } \end{array}$ is the sum of the rate of change of the cell content $\begin{array} { r } { \int _ { D ( t ) } \frac { \partial } { \partial t } Q ( \vec { x } , t ) d V } \end{array}$ within the cell volume $\mathcal { D } ( t )$ and the rate of change of the cell content $\begin{array} { r } { \int _ { \partial D ( t ) } Q ( \vec { x } , t ) \vec { v } _ { \partial D ( t ) } \cdot \hat { n } d S } \end{array}$ due to the moving cell boundary $\partial \mathcal { D } ( t )$ , where $\hat { n }$ is the unit vector normal to the cell surface $d S$ . . . . . 171 





5.18 The wave propagation diagram on the $x - t$ plane to explain the two-wave model of HLL approximate Riemann solver. . . . 173 





5.19 The wave propagation diagram on the $x - t$ plane to illustrate the main features for the three-wave model of HLLC approximate Riemann solver. . . . . . 177 





5.20 Good agreements were obtained in the benchmark tests between first and second order HLLC approximate Riemann solvers and the exact solutions in a shock tube problem. The first-order upwind HLLC solution is numerically diffusive, and is applied to update any cell with non-monotonic profile to damp the numerical noises across . . . . . 183 





5.21 Comparison of numerical noise damping capabilities on a mass density profile between MUSCL scheme for HLLC and MacCormack scheme with artificial numerical viscosities for a mode $\ell = 2 0$ for a NIF implosion simulation. . . 184 





5.22 Performance of second-order MUSCL scheme for HLLC approximate Riemann solver for a 2-D mode $\ell = 2 0$ on resolving the fluid velocity on the left and the mass density on the right. . . . . . . . 184 





5.23 Benchmark test of PPM high resolution method to boost HLLC upwind solution to the third order in the spherical-mesh version DEC3D with LILAC’s result in a shock tube problem. . . . . . . . 188 





5.24 The plot of Kramers free-free diffusion mean free path $\ell _ { \nu } = 1 / \kappa _ { \nu }$ in Eq. (5.186) using the mass density and electron profiles for a 1-D NIF implosion at stagnation. The diffusion coefficient is changed significantly in space from the hot spot, to the cold shell and finally to the outer vacuum. 203 





5.25 Benchmark tests of the 12-group radiation transport in DEC2D with $L I L A C$ in a 1-D spherical slab. The initial electron temperature and mass density profiles are $T _ { \mathrm { e } } ( r < r _ { 0 } ) = 5 ~ \mathrm { k e V }$ and $T _ { \mathrm { e } } ( r \geq$ $r _ { 0 } ) = 0 . 5$ keV; $\rho ( r < r _ { 0 } ) = 5 0 ~ \mathrm { g / c m ^ { 3 } }$ and $\rho ( r \geq r _ { 0 } ) = 1 0 0 ~ \mathrm { g / c m ^ { 3 } }$ . . 205 





5.26 The benchmark test of multi-group radiation transport implemented in DEC2D. NIF* standards for National Ignition Facility. 206 





5.27 The flow chart of parallel multi-group radiation transport implemented in DEC2D. The diffusion equation for each group is solved independently by multiple core in parallel for each time step, followed by summing over all groups to obtained the net energy exchange with the plasma through an explicit update on the electron temperature. . . . . . 207 





5.28 The performance of parallel multi-group radiation transport implemented in DEC2D. The parallel code for 4-group radiation transport is about $\sim 2 \times$ $\times$ faster than the serial code. The total amount of real computational time saved by the parallel code is increased with the resolution from $N = 2 0 0$ to $N = 8 0 0$ in the test. . . . . . 207 





5.29 Comparison of effects of radiation ablative stabilization on NIF implosions between with and without radiation transport. 208 





5.30 Effect of radiation ablative stabilization for a 2-D low mode $\ell = 6$ (above), and a high mode $\ell = 2 0$ (below). . . . . 209 





5.31 The boundary condition of zero temperature gradient is applied at the origin by replacing $B _ { 1 }$ with $B _ { 1 } ^ { * }$ and at the edge of the simulation domain by replacing $B _ { N }$ with $B _ { N } ^ { * }$ for a 1-D thermal diffusion in spherical geometry. Temperatures are discretized at the center of a cell i.e., $r _ { i } = ( i - 1 / 2 ) \triangle r$ . . . . . . . 215 





5.32 A 1-D benchmark test with $L I L A C$ for electron and ion thermal diffusions within a spherical slab using Spitzer Coulomb logarithm $\ln \Lambda _ { \mathrm { e i } } ^ { \mathrm { S p i t z e r } }$ in Eq. (5.236), and without Lee-More degeneracy correction factor $f _ { \mathrm { L M } }$ in Eq. (5.233). . . . . 218 





5.33 A 1-D benchmark test with $L I L A C$ for electron and ion thermal diffusions within a spherical slab using Lee-More Coulomb logarithm $\mathrm { l n } \Lambda _ { \mathrm { e i } } ^ { \mathrm { L M } }$ in Eq. (5.234), and with Lee-More degeneracy correction factor $f _ { \mathrm { L M } }$ in Eq. (5.233). . . . . . 219 





5.34 A 1-D benchmark test with $L I L A C$ for electron and ion thermal diffusions and equilibration using Lee-More Coulomb logarithm and degeneracy correction factor. Ions are quasi static in the plasma due to heavier mass than electrons. However, there are rapid heat exchanges between electrons and ions due to fast electron and ion collisions. . . . . . . 219 





5.35 (a) Effects of space-lagging approximations in computing the cell interface diffusion coefficient $( r ^ { 2 } \kappa ) _ { i \pm 1 / 2 } ^ { n }$ on the mass density profile at stagnation for a NIF implosion. The resolution is 200 cells in the 1-D code. (b) A benchmark test to validate the simple mean approximation to compute the cell interface coefficient at high resolution with 1000 cells in the 1-D code. . . . . 220 



5.36 The DEC2D benchmark test with $L I L A C$ for a 1-D NIF implosion 

to validate hydrodynamics, electron and ion thermal diffusions and 

equilibration. . 221 

# 1 Research Objective

# 1.1 Introduction to ICF

The nuclear fusion reaction in a Deuterium (D) and Tritium (T) plasma, 

$$
\mathrm {D} + \mathrm {T} \rightarrow {} _ {2} ^ {4} \mathrm {H e} (3. 5 \mathrm {M e V}) + \mathrm {n} (1 4. 1 \mathrm {M e V}), \tag {1.1}
$$

produces a 3.5-MeV alpha particle and a 14.1-MeV neutron. This reaction requires plasma temperatures high enough to overcome the Coulomb repulsion and brings D and T nuclei into contact. In nature, fusion reactions occur at the center of the sun, where the temperature is around $\sim 1 0 ^ { 7 }$ K. In laboratory experiments, magnetic confinement fusion (MCF) and inertial confinement fusion (ICF) [VRMV92, AtV04] are two approaches to confine high-temperature fusion plasmas. Figure (1.1) shows the typical diameter of an ICF capsule in laboratory scales. The main goal of today’s fusion experiments is to demonstrate the viability of fusion as an energy source. 

This thesis focuses on inertial confinement fusion via lasers. Laser fusion scheme uses the direct and indirect drive approach [Lin95]. Two main laser facilities are currently used for laser-driven inertial confinement fusion experiments. The OMEGA laser at the University of Rochester can deliver about 30 kJ of laser energy in direct-drive experiments. The National Ignition Facility (NIF) deliveries about 2 MJ of ultra-violet laser light and is designed to explore the feasibility of fusion ignition. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/53ec7a2c30924c7f954190a428f3d97087a1b1acbb04577b4f80c6f95800cb99.jpg)



Figure 1.1: DT fusion reactions in the laboratory scale and inside the core of the sun. The pressure and temperature of an ICF hot spot is comparable to the core of the sun.


In a direct-drive ICF implosion, laser beams irradiate uniformly to compress a spherical capsule. In an direct-drive ICF implosion, laser beams irradiate the inner walls of a cylindrical gold can called a hohlraum. A fraction of the laser energy is converted into X rays, which provide a more uniform compression on the spherical capsule than a direct-drive implosion. Figure (1.1) shows the difference of laser irradiation schemes between direct-drive and indirect-drive. 

This thesis focuses on direct-drive ICF experiments and its applications on the OMEGA laser. The 60-beam, 30-kJ, 351-nm OMEGA laser system [BBC+97], high power ( $\sim$ 30TW) laser beams are used to symmetrically compress a spherical capsule to form a hot and dense plasma core, called the hot spot, at multi-keV temperatures and tens of g/cc density. A nominal OMEGA cryogenic target [CAB $^ +$ 15] is made of a thin Carbon-Hydrogen (CH) plastic ablator ( $\sim 8$ µm thick and $\sim 1$ mm diameter) shell enclosing an equal-molar Deuterium-Tritium (DT) ice layer (∼ 50 µm thick), filled with a low-density DT gas in thermodynamic equilibrium with the ice. 

The implosion is driven by the ablation pressure produced by the laser irradiation of the outer shell surface. The laser energy is transferred to electrons in 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/1ed08794ee7e45acdd485ba509761d65accd84d7da2b9bbc7d7d511e38169278.jpg)



Figure 1.2: Laser energy deposition in direct-drive: (1) the hot electron preheat and (2) the cross-beam-energy-transfer degrade the implosion performance.


a low-density coronal plasma by means of inverse bremsstrahlung radiation absorption [ZR02]. As laser rays propagate deeper into the coronal plasma, their trajectories are refracted and experience a total reflection [Kai00] near the critical density $n _ { \mathrm { c } }$ . Parametric instabilities [Kru88] from laser-plasma interactions (LPI) occurring below $n _ { \mathrm { c } }$ degrade the ICF implosion performance in two ways: (1) the generation of hot electrons that preheat and increase the adiabat of the cold shell thus reducing the compression of the hot spot [ZB07]; (2) the acoustic wave transferred energy from an incoming ray to an outgoing ray, a process known as cross-beam energy transfer (CBET) [MFS+17] which reduces the total amount of laser energy deposited on the ablation surface. The laser propagation is sketched in Fig. (1.2) for illustration purposes. 

The absorbed laser energy is conducted from the critical surface to an ablation surface, where the CH ablator on the outer shell surface is being ablated away in the coronal plasma leading to the so-called rocket effect [AtV04] that accelerates the cold shell inward. The heat transport across the thermal conduction zone is non-local [SNB00], defined by a longer electron thermal diffusion mean-free-path 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/7f857597d863149960e919d9893c05fefd68cc58f7894b13007fc8cd224d0193.jpg)



Figure 1.3: The four stages in direct-drive ICF implosion: (1) shock launching and propagation, (2) acceleration phase, (3) deceleration phase and hot spot formation, and (4) disassembly phase and thermonuclear burn propagation.


than the local plasma temperature gradient scale length, due to the steep electron temperature profile. [SDMV81] 

Figure (1.3) summarizes the four stages for a direct-drive ICF implosion including: (1) shock launching and propagation, (2) acceleration phase, (3) deceleration phase and hot spot formation, and (4) disassembly phase and thermonuclear burn propagation. 

In the shock launching and propagation phase, the foot of the laser pulse drives a shock wave propagating radially inward, while compressing the DT ice mass density by about 4 times. Before the shock arrives at the shell inner surface, the laser pulse is ramped up rapidly to drive a sequence of compressive waves that accurately catch up the shock at the shell inner surface, because the shock wave is subsonic relative to the post-shocked fluid velocity. As the shock breaks out onto the inner shell surface, a rarefaction wave is sent radially outward and the whole shell is set in motion. Soon after the rarefaction wave arrives at the outer shell surface, the whole shell accelerates. The shock dynamics and the corresponding 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/1d244bed75fd61efcc32b1b9c24b9ea8e7d749825417a4daa74714b3952df191.jpg)



(a)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/30b534e5a6c8b615db78c18cbbfe57877fac0a83006db03c3e59d971f30540eb.jpg)



(b)



Figure 1.4: The adiabat shaping techniques [AB04] in (a) method of relaxed mass density profile, (b) method of decaying pressure.


laser pulse shape are shown together in Fig. (1.3). 

The spatial profile of the cold shell adiabat, defined by $\alpha \sim P / \rho ^ { 5 / 3 }$ where $P$ is the shell pressure and $\rho$ is the shell mass density, can be shaped by launching a prepulse. There are two types of adiabat-shaping techniques: (1) a relaxed shell mass density profile shown in Fig. (1.4)-(a) produced by introducing a prepulse followed by the laser turned off, a shock launched by a foot in the laser pulse shape and an adiabatic compression; and (2) a decaying shell pressure shown in Fig. (1.4)-(b) produced by introducing a prepulse immediately followed by an adiabatic compression. The resulting shell adiabat in both cases decreases monotonically in radius from the outer shell surface to the inner shell surface, so that a larger outer shell adiabat $\alpha _ { \mathrm { o u t } } > \alpha _ { \mathrm { i n n } }$ mitigates the RT growth during the acceleration phase. 

In the acceleration phase, Rayleigh-Taylor (RT) [Ray83] instability grows on the outer shell surface. The feedout of RT instabilities from the outer shell surface into the inner shell surface seeds an initial perturbation for the deceleration-phase RT growth. At the end of the acceleration phase, the shell achieves its peak implosion velocity, while the main shock rebounding from the center breaks out on the inner shell surface. The shell is first impulsively decelerated by the shock, 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/e2af5f93084fcc3d2fe561e55df662ad4afac8fde4059b0814d17221a1fde0dd.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/53ae9b37ab1b6e55174be9362475c3f2a2be2492eb274fa3edd1a2a86b371a11.jpg)



Figure 1.5: Comparison of classical and ablative RT instabilities. The mass ablation caused the thermal heat flux produces a dynamic pressure ${ \scriptstyle \frac { 1 } { 2 } } \rho _ { \mathrm { l } } v _ { \mathrm { b } } ^ { 2 }$ , which reduces the pressure difference between the heavy and light fluids and the effective force acting on the perturbed heavy fluid mass.


and subsequently slowed down by the hot spot pressure (continuous deceleration). 

In the compression phase, the deceleration RT leads to an incomplete conversion of the shell kinetic energy into hot-spot internal energy through PdV work, resulting in the reduction of fusion yields, distortion of hot-spot shapes, and generation of residual kinetic energy in non-stagnating RT spikes and non-stagnating hot-spot flows. 

In the disassembly phase or the burn propagation phase, ignition occurs when the rate of alpha particle energy deposition exceeds than the rates of all energy loss mechanisms including (1) the expansion losses, (2) the radiative cooling, (3) the loss of heat and (4) loss of alpha particles from the hot spot. 

# 1.1.1 Rayleigh-Taylor instability

During the implosion, nonuniformities seed the Rayleigh-Taylor (RT) hydrodynamic instabilities [Ray83], in which the tip of RT spike grows in the direction of a lower effective gravitational potential energy. Here the “effective gravity” is the inertial acceleration in the moving frame of reference of the imploding shell. The “effective gravity” points from the heavy dense shell to the light ablated plasma 

during the acceleration phase. During the deceleration phase, the “effective gravity” points from the dense shell towards the lighter hot spot plasma. The 3-D structure of bubble-and-spike grows in a way to reduce the total effective gravitational potential energy as a function of time, thus converting the reduced total potential energy into fluid kinetic energy of RT spikes and bubbles. Figure (1.5) illustrates the interchange instability between the denser and the lighter DT. 

In the presence of thermal heat flux, the shell material is ablated off the perturbed interface at the blow-off velocity $\vec { v } _ { \mathrm { b } }$ and produces a dynamic pressure ${ \scriptstyle \frac { 1 } { 2 } } \rho _ { 1 } v _ { \mathrm { b } } ^ { 2 }$ in the light fluid region, where $\rho _ { \mathrm { l } }$ is the light fluid mass density. The resulting pressure difference between the heavy $P _ { \mathrm { h } }$ and light $P _ { \mathrm { l } }$ fluids is reduced: $\begin{array} { r } { ( P _ { \mathrm { h } } - P _ { \mathrm { l } } ) S = \triangle m \frac { d ^ { 2 } } { d t ^ { 2 } } \xi } \end{array}$ , where $S$ is the perturbed surface area, $\triangle m$ is the perturbed heavy fluid mass and $\xi$ is the vertical displacement of perturbation. After the thermal ablation of the cold shell material on the perturbed interface, the effective vertical displacement $\xi$ is reduced, a process known as ablative stabilization. In the acceleration phase, the source of heat flux is caused by the electron heat conduction, which transfers the absorbed laser energy from the critical surface to the ablation surface. In the deceleration phase, the heat flux flows from the high-temperature hot spot into the cold shell. The growth of RT instabilities in the four stages of ICF implosions is summarized in Fig. (1.5). 

# 1.1.2 Nuclear Diagnostics

Implosion asymmetries are measured by (1) nuclear diagnostics measuring hot-spot ion temperatures and neutron yields, and (2) X-ray diagnostics used to measure hot-spot sizes and shapes. In this thesis, we focus on in the neutron diagnostics including the time-of-flight (nTOF) detector and the magnetic recoil neutron spectrometer (MRS) [JFC+12], which capture and measure the energy spectra of the fusion neutrons escaping from the hot spot. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/e8cac9eb85060ee1a4ebd00291bcb6326d08eafa9f0e7fba79660350972cc98a.jpg)



Figure 1.6: When non-stagnating hot-spot flow velocity is large, the neutroninferred ion-temperature measurements along different line-of-sight vary for the same imploding target.


RT instabilities introduce a significant amount of non-stagnating hot-spot flows at the time of peak compression. As a result, the neutron velocity is Dopplershifted [Mur14] and the width of the neutron energy spectrum is broadened. The magnitudes of the Doppler velocity broadening vary with the line-of-sight (LOS) locations, depending on the spatial variation of the hot-spot flow velocity. Since the widths of the neutron energy spectra are used to infer the hot-spot ion temperatures, the flow effects can lead to significant variations in ion-temperature measurements along different lines of sights. 

The basic neutron kinematics during DT nuclear reactions is shown in Fig. (1.6), in which the hot spot is perturbed by mode $\ell = 1$ . The 14.1-MeV neutron escapes from the DT center-of-mass (CM) frame, and is captured by a nTOF or MRS detector along a given LOS. The measured neutron velocity $v _ { n }$ is a sum of the 14.1-MeV initial velocity $v _ { 0 }$ , a small component of positive velocity $v _ { \mathrm { r e l } } > 0$ due to the relative kinetic energy of DT ion-pair in their CM frame, a thermal velocity component $v _ { \mathrm { t h } } ^ { \mathrm { C M } }$ due to the motion of DT center-of-mass frame, and a flow velocity component $v _ { \mathrm { H o w } }$ that boosts the CM frame velocity. Only when $v _ { \mathrm { r e l } }$ , $v _ { \mathrm { t h } } ^ { \mathrm { C M } }$ and $v _ { \mathrm { H o w } }$ are all ideally zero, the shape function of the neutron energy spectrum is a delta function sharply peaked at 14.1-MeV. Otherwise, the neutron energy 

spectrum is approximately Gaussian with a width uniquely proportional to the square root of the neutron velocity variance of $v _ { n } - v _ { 0 }$ , 

$$
\operatorname {v a r} \left[ v _ {n} - v _ {0} \right] = \operatorname {v a r} \left[ v _ {\text {t h}} ^ {\mathrm {C M}} \right] + \operatorname {v a r} \left[ v _ {\text {f l o w}} \right] + \operatorname {v a r} \left[ v _ {\text {r e l}} \right]. \tag {1.2}
$$

The first term in Eq. (1.2) is the thermal velocity broadening contributed by the DT center-of-mass randomness motion at a finite temperature, the second term is the Doppler velocity broadening caused by the non-stagnating hot-spot flow velocity. The last term is the velocity broadening caused by the relative DT kinetic energy, and is typically negligible. Detailed derivation for Eq. (1.2) using relativistic neutron kinematics is discussed in Chapter 3. Effects of non-stagnating hot-spot flow motions captured by the term var[vflow] is the subject in this thesis. 

# 1.2 Understanding 3-D implosion asymmetries

The impact of implosion asymmetry on degrading ICF implosion performance is reflected on the Lawson’s performance parameter $P \tau$ , where $P$ and $\tau$ are the hot-spot pressure and the energy confinement time. Since the value of the hotspot pressure is inferred from the experimental yield and the neutron-inferred ion temperature, 3-D effects such as non-stagnating hot-spot fluid motions resulting in larger apparent ion temperatures than the true thermal ion temperatures, can affect the assessment of implosion performance. 

Developing a practical physical model for 3-D implosion asymmetry enable us to understand the mechanism and the formation of experimental signatures in low-mode and high-mode distorted implosions. This work requires modeling 3-D implosions, followed by post-processing the hydrodynamic data with neutrontransport codes to reconstruct the neutron energy spectra, or radiation-transport codes to reconstruct the X-ray images for distorted hot spots. Qualitatively, the 

comparison of post-processed data with that of experiments provides a clue to interpret the sources of perturbations in ICF implosions. Hydrodynamic instabilities play an important role driving these implosion asymmetries. A useful physical model must be able to quantitatively explain the signatures of implosion asymmetries in experimental observable, and provides the prediction capability to improve the experimental designs. 

In this thesis, a three-dimensional (3-D) radiation-hydrodynamics computer code DEC3D was developed to study the effects of 3-D non-uniformities. DEC3D provides an efficient simulation platform to study 3-D effects of RT instabilities in the deceleration phase of ICF implosions. A deceleration-phase single-mode and multi-mode simulation database was generated to study 3-D hydrodynamic scaling relations for RT instabilities and experimental signature of implosion asymmetries. The application of DEC3D simulation database led to an improved understanding of effects of residual kinetic energy (RKE) on yield degradations [WBS+18a] and the impact of three-dimensional hot-spot flow asymmetry on ion-temperature measurements [WBS+18b]. A comprehensive code description is presented in Chapter 5. This chapter includes an overview on approximate Riemann solver techniques for advanced computational fluid dynamics simulations, multi-group radiation diffusion models, alpha-particle transport models. 

The main physics results in this thesis can be summarized as follows. First, an analytical 3-D hot-spot model was developed [WBS $^ +$ 18a] to explain the degradation of fusion yields, hot-spot pressures, and the property of larger hot-spot volumes for low modes, defined by Legendre mode number $\ell \leq 6$ for spherical harmonic perturbations [AtV04]. The 3-D hydrodynamic effect on degrading the hot spot formation is explained in terms of growing residual kinetic energy (RKE), defined by the total amount of non-stagnating fluid kinetic energy at the time of peak compression. During the evolution of low modes in the deceleration phase, a large fraction of the shell kinetic energy is converted into fluid kinetic energy 

for non-stagnating RT spikes and bubbles. The physical mechanism of degrading the hydrodynamic efficiency in converting the shell’s kinetic energy into hot-spot internal energy through PdV work is captured by the 3-D hot-spot model. Low modes are shown to exhibit a large variation in ion-temperature measurements among different lines of sights. These results are presented in Chapter 2. 

Second, a generalized method of velocity variance decomposition was developed [WBS+18b] for the non-relativistic Brysk ion temperatures. [Bry73, Mur14, Mun16]. The velocity variance $\mathrm { v a r } [ v _ { \mathrm { H o w } } ]$ in Eq. (1.2) is decomposed into six hotspot flow parameters that uniquely determine the variation of ion-temperature measurements along different lines of sights. This technique is valid for hot spot distorted by single modes, multi modes, or turbulent flows. A well-behaved relation between the fluid properties of isotropic and anisotropic flow structures and the location of LOS for ion-temperature measurements is formed. This property is used to predict the true minimum ion temperature using ion-temperature measurements at six LOS. Low modes $\ell = 2$ are shown to exhibit a large minimum ion temperature for modest amplitude perturbations, defined by the pair of RT spikes not reaching the center. These results are presented in Chapter 3. 

Third, the method of velocity variance decomposition is applied to analyze OMEGA experimental data to extrapolate the minimum DD neutron-inferred ion temperatures T DDmin. $T _ { \mathrm { m i n } } ^ { \mathrm { L D D } }$ An improved correlation between the experimental yield and $T _ { \mathrm { m i n } } ^ { \mathrm { L D D } }$ was observed. An approximate closure was developed to diagnose the hotspot flow isotropy using DD and DT ion-temperature measurements. The physical mechanism for larger apparent ion temperatures than the true thermal ion temperatures due to the isotropic flow was studied through the analysis of the isotropic part of the velocity variance. These results are presented in Chapter 4. 

Chapter 6 is our conclusion. 

# 2 Impact of Residual Kinetic Energy

# 2.1 Motivation

In inertial confinement fusion (ICF) experiments, the nonuniformities in the target and in the laser illumination seed Rayleigh–Taylor (RT) instabilities in the acceleration and deceleration phases of ICF implosions thereby degrading the implosion performance. Qualitative explanations for the decreasing hydrodynamic efficiency to convert shell’s kinetic energy into hot spot’s internal energy were reported by Scott et al. [SCB $^ +$ 13], Spears et al. [SEH $^ +$ 14], Kritcher et al. [KTB+14] and Gu et al. [GDF $^ { 1 + }$ 14]. The complexity of nonlinear hydrodynamics is the main challenge to the derivation of analytic models to describe 3-D effects. Sanz et al. [SB05a, SGC+05] reported the first semi-analytic 3-D hot-spot model for the deceleration phase through linearized perturbation analysis and numerical solution for the 3-D Poisson equation for the perturbed velocity field. 

In this chapter, a new technique to derive an analytic 3-D hot-spot model [WBS+18a] for the deceleration phase of ICF implosions is presented, based on the idea to connect 1-D and 3-D hot-spot parameters through the conservation of the hot-spot adiabatic parameter [WBY+16]. 

# 2.2 DEC3D simulation database

DEC3D is used to develop a simulation database of different OMEGA implosions with single-mode and multi-mode perturbations. This provides a platform 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/ac6ed19f862f11dbc8939a3bcb951dd3a450104ee22037d0921564fed3f94e74.jpg)



Figure 2.1: The initial mass-density profile for OMEGA shot 77068 at the beginning of deceleration phase $t _ { 0 } = 2 . 2 8$ ns. The initial radial velocity perturbation $\delta v _ { r } ( t _ { 0 } )$ is applied on the inner shell surface r0 = 68 µm. The blue line shows the shape function $f ( r )$ that is unity at $r _ { 0 }$ .


for systematic studies of the degradation of implosion performance by 3-D nonuniformities in the deceleration phase of ICF implosions. 

Implosion 77068 is considered here as a test case because it is one of the bestperforming implosions to date on OMEGA, with a hot-spot pressure $> 5 0$ GBar [RGI $^ +$ 16]. The target has an 8-µm thickness CD ablator, an initial outside radius of 430 $\mu$ m and a DT-ice thickness of 50 µm. The implosion has adiabat $\alpha = 3 . 2$ , defined by the mass-averaged ratio of the fuel pressure to the Fermi-degenerate pressure $P / P _ { \mathrm { F e r m i } }$ , and a convergence ratio of about 20. It was shown that the core conditions for a hydrodynamic-equivalent 77068 implosion at a NIF energy scale could lead to significant alpha heating [BWB $^ +$ 16]. 

The DEC3D single-mode database is generated by introducing initial radial velocity perturbations $\delta v _ { r } ( t _ { 0 } )$ on the inner shell surface located at $r _ { 0 } = 6 8$ µm for OMEGA shot 77068 at the beginning of the deceleration phase. The time $t _ { 0 } = 2 . 2 8$ ns is when the shell has reached the maximum implosion velocity. Figure (2.1) shows the initial one-dimensional (1-D) profiles for $\rho , v _ { r } , P , P _ { \mathrm { e } }$ at $t _ { 0 }$ obtained from the 1-D $L I L A C$ simulation of 77068; $L I L A C$ is a 1-D Lagrangian radiation-

hydrodynamic code [DER $^ +$ 87] routinely used for target designs at the Laboratory of Laser Energetics. Its main capabilities include refractive ray tracing, crossbeam energy transfer, nonlocal electron thermal transport [CMD15], and firstprinciple equation of state [HMGS11]. Initial transverse velocities $v _ { \theta } , v _ { \phi }$ are zero. The initial perturbation is applied on the radial velocity component in forms of spherical-harmonic modes, 

$$
v _ {r} (\theta , \phi , r, t _ {0}) = v _ {r} ^ {L I L A C} (r, t _ {0}) + \frac {\triangle v}{v _ {0}} f (r) Y _ {\ell} ^ {m} (\theta , \phi). \tag {2.1}
$$

Legendre modes $\ell = 1$ to 12 are studied, and the initial radial velocity perturbation is $\begin{array} { r } { \delta v _ { r } ( t _ { 0 } ) = \frac { \triangle v } { v _ { 0 } } f ( r ) Y _ { \ell } ^ { m } ( \theta , \phi ) } \end{array}$ , where $v _ { 0 } ( r _ { 0 } , t _ { 0 } )$ is the fluid velocity. The sphericalharmonic functions are normalized and orthogonal. The single-mode database contains 2-D $m = 0$ modes with $\ell = 1 - 1 2$ , and 3-D $m \neq 0$ modes with $m = \ell _ { \mathrm { e v e n } } / 2$ only for even $\ell$ -modes. The perturbation levels $\triangle v / v _ { 0 }$ in the simulation database are in the range of $1 \%$ to $1 4 \%$ , large enough to degrade implosions with yield-overclean $\mathrm { ( Y O C } \equiv Y _ { \mathrm { 3 D } } / Y _ { \mathrm { 1 D } }$ , defined by the ratio of 3-D to 1-D yields) $\sim 0 . 2$ in the deceleration phase simulations. Figure (2.2) shows the 3-D electron temperature contour surface at 1 keV for all single modes in the DEC3D single-mode simulation database. The burn surface defined at 1 keV has the property of vanishing enthalpy and heat flux exchanged leading to the conservation of the hot spot adiabatic parameter $P _ { \mathrm { H S } } V _ { \mathrm { H S } } ^ { 5 / 3 }$ between 3-D and 1-D hot spots, where $P _ { \mathrm { H S } }$ and $V _ { \mathrm { H S } }$ are hotspot pressure and volume respectively. $f ( r )$ is a shape function that is unity at $r _ { 0 }$ and decays radially away from $r _ { 0 }$ to avoid loading too much initial energy perturbation on the sharp interface at $r _ { 0 }$ , 

$$
f (r) = \frac {1}{2} \left(\frac {r}{r _ {0}}\right) ^ {\ell} \left[ 1 - \tanh \left(\frac {r - r _ {0}}{0 . 0 1 5 r _ {0}}\right) \right] + \frac {1}{2} \left(\frac {r _ {0}}{r}\right) ^ {\ell} \left[ 1 + \tanh \left(\frac {r - r _ {0}}{0 . 0 1 5 r _ {0}}\right) \right]. \tag {2.2}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/c85fb0e58fa1ce4523956b24be83fa77abf5bc4219eb705ddb1a90552ca0f826.jpg)



Figure 2.2: The summary of DEC3D single-mode simulation database for the shot 77068. (1) HLLC approximate Riemann solver with PPM high-resolution method; (2) HYPRE implicit thermal diffusion; and (3) high angular resolution $1 2 8 \times 2 5 6$ for $\theta$ and $\phi$ zones.


# 2.3 Property of low-mode hot-spot volumes

We began with the investigation of 3-D deceleration-phase hydrodynamic instabilities for low modes. We expect larger 3-D burn-averaged hot-spot volumes and lower pressures than 1-D’s, 

$$
V _ {\mathrm {H S}} ^ {\mathrm {3 D}} (t) / V _ {\mathrm {H S}} ^ {\mathrm {1 D}} (t) \geq 1. \tag {2.3}
$$

This can be shown by using the subsonic flow approximation of the hot-spot fluid energy equation. 

$$
\frac {\partial}{\partial t} \left(\frac {P}{\gamma - 1}\right) + \vec {\nabla} \cdot \vec {v} \left(\frac {\gamma P}{\gamma - 1}\right) = \vec {\nabla} \cdot \kappa_ {\mathrm {e}} \vec {\nabla} T _ {\mathrm {e}} + \dot {\varepsilon} _ {\alpha} - \vec {\nabla} \cdot \vec {F} _ {\mathrm {r}}, (2. 4)
$$

where $P$ is the hot-spot fluid pressure, $\vec { v }$ is the hot-spot fluid velocity, $\gamma = 5 / 3$ is the ratio of ideal gas specific heats, $\kappa _ { \mathrm { e } }$ is the electron thermal conductivity and $T _ { \mathrm { e } }$ is the electron temperature. $\dot { \varepsilon } _ { \alpha } = \frac { 3 } { 2 } P ^ { 2 } / S _ { \alpha }$ is the rate of alpha particle energy deposition. The proportionality constant is $\begin{array} { r } { S _ { \alpha } = \frac { 2 4 T _ { \mathrm { i } } ^ { 2 } } { \langle \sigma v \rangle _ { \mathrm { D T } } E _ { \alpha } } } \end{array}$ hσvi E , where Ti is the ion 24T 2i $\mathrm { \Delta } T _ { \mathrm { i } }$ thermal temperature, $\langle \sigma v \rangle _ { \mathrm { D T } }$ is the DT fusion reactivity and $E _ { \alpha } = 3 . 5$ -MeV is the 

alpha particle birth kinetic energy. For simplicity $7 < T _ { \mathrm { i } } ^ { \prime } < 2 3$ KeV are considered here so that $S _ { \alpha }$ is about constant. $\vec { F _ { \mathrm { r } } } = \vec { F _ { \mathrm { r } } } ^ { \mathrm { f f } } + \vec { F _ { \mathrm { r } } } ^ { \mathrm { a b s } }$ is the total radiation energy flux leaving the hot spot through the free-free X-ray emission $\vec { F } _ { \mathrm { r } } ^ { \mathrm { f f } }$ and the absorption of low-energy photons $\vec { F } _ { \mathrm { r } } ^ { \mathrm { a b s } }$ . The rate of bremsstrahlung emission loss within the hot spot $\begin{array} { r } { \dot { \varepsilon } _ { \mathrm { r } } ^ { \mathrm { f f } } = \int _ { V _ { \mathrm { H S } } ^ { 3 \mathrm { D } } ( t ) } \vec { \nabla } \cdot \vec { F } _ { \mathrm { r } } ^ { \mathrm { f f } } d V = C _ { \mathrm { r } } P ^ { 2 } T _ { \mathrm { e } } ^ { - 3 / 2 } \mathrm { W } / \mathrm { m } ^ { 3 } \ \vert } \end{array}$ [AtV04] is obtained by an volume integration, with a proportionality constant $C _ { \mathrm { r } } = 3 . 8 8 \times 1 0 ^ { - 2 9 } Z ^ { 3 } / ( 1 + Z ) ^ { 2 }$ , where $Z$ is the atomic number, the hot-spot pressure $P$ in Pacal and the electron temperature $T _ { \mathrm { e } }$ in Joule. 

Both $\vec { \nabla } \cdot \kappa _ { \mathrm { e } } \vec { \nabla } T _ { \mathrm { e } }$ and $\vec { \nabla } { \cdot } \vec { F } _ { \mathrm { r } } ^ { \mathrm { a b s } }$ are not considered as loss terms for low modes after integrating Eq. (2.4) over the 3-D hot-spot volume $V _ { \mathrm { H S } } ^ { \mathrm { 3 D } }$ , because the rate of the increasing electron internal energy $\frac { 3 } { 2 } \dot { N _ { \mathrm { e } } } T _ { \mathrm { e } }$ due to the inflow of ablated cold shell material (which is assumed being raised to the same hot-spot electron temperature rapidly) from the inner shell surface compensates the rate of the decreasing electron internal energy $\textstyle { \frac { 3 } { 2 } } N _ { \mathrm { e } } { \dot { T } } _ { \mathrm { e } }$ due to the drop of electron temperature. 

$$
\int_ {V _ {\mathrm {H S}} ^ {\mathrm {3 D}} (t)} \left(\vec {\nabla} \cdot \kappa_ {\mathrm {e}} \vec {\nabla} T _ {\mathrm {e}} - \vec {\nabla} \cdot \vec {F} _ {\mathrm {r}} ^ {\mathrm {a b s}}\right) d V = \frac {3}{2} \left(\dot {N} _ {\mathrm {e}} T _ {\mathrm {e}} + N _ {\mathrm {e}} \dot {T} _ {\mathrm {e}}\right) = 0 \tag {2.5}
$$

The bracket $\left( { \dot { N } } _ { \mathrm { e } } T _ { \mathrm { e } } + N _ { \mathrm { e } } { \dot { T } } _ { \mathrm { e } } \right) < 0$ in Eq. (2.5), however, is negative for high modes, because some ablated cold shell material is trapped within the cold bubbles, but not recycled back to the burn-averaged hot-spot volume. 

The consequence of the subsonic flow approximation for low modes results in a spatially uniform hot-spot pressure, 

$$
P _ {\mathrm {H S}} ^ {\mathrm {3 D}} (t) = 2 n _ {\mathrm {H S}, \mathrm {i}} ^ {\mathrm {3 D}} (\vec {x}, t) T _ {\mathrm {H S}, \mathrm {i}} ^ {\mathrm {3 D}} (\vec {x}, t), \tag {2.6}
$$

where n HS,i $n _ { \mathrm { H S , i } } ^ { \mathrm { 3 D } }$ 3D is the hot-spot ion number density and the thermal equilibrium between ion and electron thermal temperatures $T _ { \mathrm { H S , i } } ^ { \mathrm { 3 D } } = T _ { \mathrm { H S , e } } ^ { \mathrm { 3 D } }$ is assumed. Equation (2.6) is not valid for mid and high modes, because flows within the bubbles become 

close to sonic. Substituting the following vector identity, 

$$
\vec {\nabla} \cdot \frac {\vec {r}}{3} = 1, \tag {2.7}
$$

into Eq. (2.4), treating the hot-spot pressure $P$ and the alpha heating source term $\dot { \varepsilon } _ { \alpha }$ as uniform in space and neglecting the radiative cooling term $\dot { \varepsilon } _ { r } ^ { \mathrm { H } }$ for a alphaheating dominated plasma leads to the following simplified form of the energy equation [SB05b], 

$$
\vec {\nabla} \cdot \vec {v} = \vec {\nabla} \cdot \left[ \frac {\gamma - 1}{\gamma P} \left(\kappa_ {\mathrm {e}} \vec {\nabla} T _ {\mathrm {e}} - \vec {F} _ {\mathrm {r}}\right) + \frac {(\gamma - 1) \dot {\varepsilon} _ {\alpha} - \dot {P}}{3 \gamma P} \vec {r} + \vec {v} _ {\mathrm {r o t}} + \vec {\nabla} \phi \right]. (2. 8)
$$

$\vec { v } _ { \mathrm { r o t } }$ is the rotational flow that satisfies the zero divergence property $\vec { \nabla } \cdot \vec { v } _ { \mathrm { r o t } } = 0$ , and $\vec { \nabla } \phi$ is the potential flow induced by RT instabilities that satisfies the Laplacian equation $\nabla ^ { 2 } \phi$ . Therefore, the 3-D burn-averaged hot-spot fluid velocity $\vec { v }$ for low modes is a sum of the potential flow $\vec { v } _ { \mathrm { p o t } }$ and the rotational flow $\vec { v } _ { \mathrm { r o t } }$ . The spatial dependence of $\vec { v } ( \vec { x } , t )$ is caused by the spatial variation of the hot-spot ion number density and the hot-spot ion thermal temperature. 

$$
\vec {v} (\vec {x}, t) = \underbrace {\frac {\gamma - 1}{\gamma P} \left(\kappa_ {\mathrm {e}} \vec {\nabla} T _ {\mathrm {e}} - \vec {F} _ {\mathrm {r}}\right) + \frac {(\gamma - 1) \dot {\varepsilon} _ {\alpha} - \dot {P}}{3 \gamma P} \vec {r} + \vec {\nabla} \phi} _ {\vec {v} _ {\mathrm {p o t}}} + \vec {v} _ {\mathrm {r o t}}. (2. 9)
$$

Equation (2.9) contains information about the flow structure for low modes including the total perturbed hot-spot volume to be analyzed in this chapter, as well as the velocity variance to be analyzed in Chapter 3-4. Since the low-mode hot-spot pressure is treated as a variable being uniform in space in Eq. (2.9), the anisotropic part of the hot-spot fluid velocity variance is determined by the spatial variations of the electron thermal ion temperature and the vorticity, whereas the isotropy part is determined by the radially convergent flow which is proportional to the vector $\bar { r }$ . By treating the hot-spot pressure as uniform in space in Eq. 

(2.9), and assuming a zero correlation between each term on the right-hand-side of Eq. (2.9), leads to the conclusion that the hot-spot fluid velocity variance for low modes is, 

$$
\begin{array}{l} \mathrm {v a r} [ \vec {v} \cdot \hat {d} ] = \mathrm {v a r} \left[ \frac {\gamma - 1}{\gamma P} \left(\kappa_ {\mathrm {e}} \vec {\nabla} T _ {\mathrm {e}} - \vec {F} _ {\mathrm {r}}\right) \cdot \hat {d} \right] + \mathrm {v a r} \left[ \frac {(\gamma - 1) \dot {\varepsilon} _ {\alpha} - \dot {P}}{3 \gamma P} \vec {r} \cdot \hat {d} \right] \\ + \operatorname {v a r} [ \vec {\nabla} \phi \cdot \hat {d} ] + \operatorname {v a r} [ \vec {v} _ {\mathrm {r o t}} \cdot \hat {d} ], \tag {2.10} \\ \end{array}
$$

where $\hat { d }$ is the line-of-sight unit vector. Figure (2.3) compares the 3-D profiles of vorticity magnitude $| \vec { \omega } |$ between a 2-D single-mode $\ell = 8 , m = 0$ and a 3-D single-mode $\ell = 8 , m = 4$ , simulated by the Cartesian-mesh version DEC3D. The boundaries of Cartesian-mesh domain decomposition are shown by the vertical and horizontal gaps in the plots. The topology of 2-D vorticity appears to be ring-structure in-shape rotating around the $z$ -axis, whereas the topology of 3-D vorticity exhibits the shape of a spherical-structure. As a result, the velocity variance for 2-D single-modes with $m = 0$ are more anisotropic than that for 3-D single-modes with $m \neq 0$ . The impact of this property on neutron-inferred ion-temperature measurement asymmetry will be discussed in Chapter 4. 

To understand the formation of vortex rings [Dav15, Saf92] within the hot spot as shown in Fig. (2.3), the general evolution equation of vorticity is obtained by taking the curl of the inviscid fluid momentum equation, 

$$
\partial_ {t} \vec {\omega} + \vec {v} \cdot \vec {\nabla} \vec {\omega} = (\vec {\omega} \cdot \vec {\nabla}) \vec {v} - \vec {\omega} (\vec {\nabla} \cdot \vec {v}) + \frac {\vec {\nabla} \rho \times \vec {\nabla} P}{\rho^ {2}}, \tag {2.11}
$$

where the identity $\vec { \nabla } ( v ^ { 2 } / 2 ) = ( \vec { v } \cdot \vec { \nabla } ) \vec { v } + \vec { v } \times \vec { \omega }$ is used. The divergence term can be absorbed using the mass density equation $\begin{array} { r } { \vec { \nabla } \cdot \vec { v } = - ( 1 / \rho ) \frac { D \rho } { D t } } \end{array}$ , where $\begin{array} { r } { \frac { D } { D t } = \partial _ { t } + \vec { v } \cdot \vec { \nabla } } \end{array}$ is the material time derivative. 

$$
\frac {D}{D t} \left(\frac {\vec {\omega}}{\rho}\right) = \left(\frac {\vec {\omega}}{\rho} \cdot \vec {\nabla}\right) \vec {v} + \frac {\vec {\nabla} \rho \times \vec {\nabla} P}{\rho^ {3}}. (2. 1 2)
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/1932528bd221311ce92819b37905462ae27e8bdbf95c034f6f18907f434160dc.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/02ecc56813880d7455962b820839a372c918bbf1fedc84ea1962cca014ed1dfc.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/55550bd594e936631748c42fde1a52b4ca76f9a0c91653d87d0770df8ac09d1a.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/702b2ac4e9248d4abe271ba99d81b589fb8e9a7013a0f56218ce1a28980d7de9.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/df7e819b27af3265f60b262e37986baebfd1b1050a0f3be3f7b6d1da9fd4e14a.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/ff5684add6c9720a624a08eaf43faeb070e5f3d733dbad54339f3d41396f1618.jpg)



Figure 2.3: 3-D profiles for mass densities and vorticity magnitudes $| \omega |$ at stagnation: for 2-D single-mode $\ell = 8 , m = 0$ in (a) and (b); for 3-D single-mode $\ell = 8 , m = 4$ in (c) and (d), where $\vec { \omega } = \vec { \nabla } \times \vec { v }$ .


The dynamics of vorticity within the hot spot in Eq. (2.12) is the following. First, the baroclinic term $\vec { \nabla } \rho \times \vec { \nabla } P$ is non-zero whenever the gradients of density and pressure are not parallel to each other across the perturbed interface of the cold shell. The vorticity generated by the baroclinic term is originally localized on the perturbed interface, and is advected into the central hot spot. In the nonlinear stage of RT instabilities, the vorticity is deposited on the surface of the mushroom structure of RT spikes, [SBRR04] and is ablated into the hot spot driven by the incoming heat flux due to the temperature difference between the high-temperature hot spot and the cold shell. Second, the vorticity stretching term $( \vec { \omega } \cdot \vec { \nabla } ) \vec { v }$ is non-zero only for 3-D non-axisymmetric flows, which holds for all 3-D $m \neq 0$ modes such as $\ell = 8 , m = 4$ in (c) and (d) of Fig. (2.3). The presence of vorticity stretching can intensifies the transport of fluid mass from the head of bubbles to the tip of spikes, [DG90, DGDH93] leading to a thinner shell and a higher risk of losing confinement. 

Integrating Eq. (2.9) over the total hot-spot volume $V _ { \mathrm { H S } } ^ { \mathrm { 3 D } }$ bounded by the perturbed cold shell interface, where the total heat flux is zero $\kappa _ { \mathrm { e } } \vec { \nabla } T _ { \mathrm { e } } - \vec { F _ { \mathrm { r } } } ^ { \mathrm { a b s } }  0$ as explained by Eq. (2.5), leads to an evolution equation of the hot spot volume. The rate of change of the total perturbed burn-averaged hot-spot volume, which is defined by $\begin{array} { r } { d V / d t = \int \vec { \nabla } \cdot \vec { v } d ^ { 3 } x = \int \vec { v } \cdot d \vec { S } } \end{array}$ , for low modes is, 

$$
\frac {d}{d t} V (t) = \left[ \frac {(\gamma - 1) \dot {\varepsilon} _ {\alpha} - \dot {P}}{\gamma P} \right] \int_ {V _ {\mathrm {H S}} ^ {\mathrm {3 D}}} \frac {\vec {r}}{3} \cdot d \vec {S} - \left(\frac {\gamma - 1}{\gamma P}\right) \int_ {V _ {\mathrm {H S}} ^ {\mathrm {3 D}}} \vec {\nabla} \cdot \vec {F} _ {\mathrm {r}} ^ {\mathrm {f f}} d V, \qquad (2. 1 3)
$$

where the surface integral $\begin{array} { r } { \int _ { V _ { \mathrm { H S } } ^ { 3 \mathrm { D } } } \frac { \vec { r } } { 3 } \cdot d \vec { S } = \int _ { V _ { \mathrm { H S } } ^ { 3 \mathrm { D } } } \vec { \nabla } \cdot \frac { \vec { r } } { 3 } d ^ { 3 } x } \end{array}$ is equal to the total hot-spot volume by using the vector identity in Eq. (2.7), whereas the volume integral vanish because of the zero divergence for rotational flows. Therefore, the time 

evolution equation is obtained for low modes’ 3-D total hot-spot volume: 

$$
\frac {d V _ {\mathrm {H S}} ^ {\mathrm {3 D}}}{d t} = \left[ \frac {(\gamma - 1) (\dot {\varepsilon} _ {\alpha} - \dot {\varepsilon} _ {r} ^ {\mathrm {f f}}) - \dot {P} _ {\mathrm {H S}} ^ {\mathrm {3 D}}}{\gamma P _ {\mathrm {H S}} ^ {\mathrm {3 D}}} \right] V _ {\mathrm {H S}} ^ {\mathrm {3 D}}, \tag {2.14}
$$

and equivalently, for low modes’ 3-D hot-spot adiabatic parameter: 

$$
\frac {d}{d t} \ln (P V ^ {\gamma}) _ {\mathrm {H S}} ^ {\mathrm {3 D}} = (\gamma - 1) \frac {\dot {\varepsilon} _ {\alpha} - \dot {\varepsilon} _ {r} ^ {\mathrm {f f}}}{P _ {\mathrm {H S}} ^ {\mathrm {3 D}}}. (2. 1 5)
$$

For OMEGA implosions, alpha heating and radiation loss are weak enough to be neglected compared with PdV work. Taking $( \dot { \varepsilon } _ { \alpha } + \dot { \varepsilon } _ { r } ) / ( \frac { 3 } { 2 } P _ { \mathrm { H S } } ^ { \mathrm { 3 D } } ) \ll 1$ in Eq. (2.15), yields $\begin{array} { r } { \frac { d } { d t } \ln { \left( P V ^ { \gamma } \right) } _ { \mathrm { H S } } ^ { \mathrm { 3 D } } = 0 } \end{array}$ for the conservation of the 3-D hot-spot adiabatic parameter $P V ^ { \gamma }$ for low modes relating 1-D and 3-D hydrodynamic quantities through the invariant $P V ^ { \gamma }$ : 

$$
\frac {V _ {\mathrm {H S}} ^ {\mathrm {3 D}}}{V _ {\mathrm {H S}} ^ {\mathrm {1 D}}} = \left(\frac {P _ {\mathrm {H S}} ^ {\mathrm {3 D}}}{P _ {\mathrm {H S}} ^ {\mathrm {1 D}}}\right) ^ {- 1 / \gamma}. \tag {2.16}
$$

Equation (2.16) states that low modes have larger 3-D hot-spot volumes than 1- D’s as long as the burn-averaged hot-spot pressures degrade compared with 1-D’s $P _ { \mathrm { H S } } ^ { \mathrm { 3 D } } / P _ { \mathrm { H S } } ^ { \mathrm { 1 D } } < 1$ . 

# 2.4 A 3-D hot-spot model

# 2.4.1 Dynamics of 3-D hot-spot volumes

Sketches of burn volumes for low and high modes on a 2-D plane are shown in Fig. (2.4). For low modes, bubbles are hot enough to sustain DT fusion reactions. The volume of hot bubbles is a part of the burn volume $V _ { b } ( t )$ that contributes to neutron productions. As the shell implodes and RT spikes grow, the rate of change in the burn volume is determined by the boundary velocity $\vec { v } _ { b } ( t )$ of the 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/b4f933e43bb7f6944d3a8b9372ab880961482a3f05285c95f22fa692078efb60.jpg)



TC14056J1


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/040deb370943269306482710f49fd19a13421136740e76d8661fdd1494a079fd.jpg)



Figure 2.4: The 3-D hot-spot models for low and high modes. The red part is the burn volume while the blue part is the shell. As the shell implodes, the burn volume $V _ { \mathrm { b } } ( t )$ is shrinking with the boundary velocity ${ \vec { v } } _ { \mathrm { b } } ( t )$ . For high modes, the cold bubbles shown by the white part do not contribute DT fusion reactions, and the actual burn volume is smaller than the total hot-spot volume bounded by the perturbed inner shell surface.


burn surface. For high modes, the boundary velocity is not equal to the fluid velocity of the perturbed inner shell surface but depends on the spatial profile of the neutron rate at different times. Inside the burn volume, physical properties are approximately uniform in space. For high modes, bubbles are too cold for fusion reactions, the burn volume is reduced to the clean volume surrounded by RT spikes, and the boundary velocity is not the same as the fluid velocity on the perturbed inner shell surface as shown in Fig. (2.4). 

A simple 3-D hot-spot model is shown in Fig. (2.4). For high modes, the cold bubbles do not contribute to DT fusion reactions and the actual burn volume is smaller than the perturbed hot-spot volume [KS01, BBSW17]. The volumeaveraged quantities with hot-spot electron temperatures greater than 1 keV are denoted by $Q \mathrm { { _ { 1 k e V } } }$ . The 1-keV contour is assumed to enclose the burn volume. Volumes of cold bubbles of the distorted hot spot are assumed to not contribute to fusion yields. The neutron-averaged or burn-averaged quantities are denoted by $Q _ { \mathrm { b } }$ . In this section, the analytical yield degradation model is first derived for low modes using the volume-averaged definition. For low modes, $Q _ { \mathrm { { 1 k e V } } } ~ \simeq ~ Q _ { \mathrm { { b } } }$ 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/8aa75795e0e9b735d78b9a950affa34b1b2c88b811c2c857265625c76308810f.jpg)



TC14058J1



Figure 2.5: The 3-D hot-spot shapes for the electron-temperature contour surface (1, 1.5, 2, 2.5 keV) for high mode $Y _ { \ell = 1 2 } ^ { m = 6 }$ at stagnation. The length scale is the same in all figures. The left-most hot-spot shape includes the cold bubbles as shown in Fig. (3.1. It corresponds to the total perturbed hot-spot volume while the right-most figure is the hot core for DT fusion reactions.


because neutrons produced within the warm bubbles are of the same order of magnitude of neutrons produced within the hot core [BBSW17]. The neutronaveraged quantities are used to characterize the yield degradation model for high modes because the burn volume for high modes is much smaller than the total perturbed hot-spot volume due to the formation of cold bubbles [BBSW17] as shown in Fig. (2.5). 

For low modes, the yield degradation is mainly caused by the incomplete conversion of shell’s kinetic energy into hot spot’s internal energy. To quantify yield degradations, an accurate 3-D hot-spot volume tracking is required in order to distinguish the useful part of clean volumes’ internal energy that produce fusion yields and the useless part of bubbles’ internal energy. 

The survey of deceleration-phase single-mode energetics at stagnation for the shot 77068 is shown in Fig. (2.6) using the burn-averaged and volume-averaged definitions. In Fig. (2.6)-(c) the red dots are Legendre modes with $m \ = \ 0$ , while blue dots are 3-D Legendre modes $\ell = 4 , 6 , 8 , 1 0 , 1 2$ with $m = \ell / 2$ . The 1-D energies at stagnation are $\mathrm { K E _ { t o t } ^ { \mathrm { 1 D } } = 2 0 8 \ J }$ , $\mathrm { I E _ { s h e l l , b } ^ { 1 D } = 8 3 0 \ J }$ , $\mathrm { I E _ { h o t s p o t , b } ^ { 1 D } = 5 0 4 }$ 

J, $\mathrm { I E _ { s h e l l ,1 k e V } ^ { 1 D } = 5 7 0 ~ J }$ , IE1Dhotspot,1keV = 764 J. Following behaviors are observed for initial velocity perturbations equal to $1 4 \ \%$ of the implosion velocity. The 3-D single modes with $m = \ell / 2$ are shown to have more total residual kinetic energy than 2-D single modes. The burn-averaged shell internal energy is shown to increase with Legendre mode number because of inclusion of cold bubbles for high modes, whereas the burn-averaged hot-spot internal energy is shown to decrease with mode numbers because of reduction in burn volumes. For the hot spot defined by $T _ { \mathrm { e } } \geq 1$ keV, internal energies for both hot spot and shell do not show significant variation among different modes. A trend indicates low modes to have more total residual kinetic energies at stagnation than high modes. The low mode $\ell = 1$ is shown to have the largest hot-spot kinetic energy within the burn volume at stagnation. 

A remarkable feature of $\mathrm { K E _ { H S } / I E _ { H S } } ( t _ { \mathrm { s t a g } } ) \leq 0 . 1 4$ is observed in the simulation database, which suggests that the subsonic flow approximation for the energy equation of hot spot is an appropriate assumption in both 1-D and 3-D implosions. Therefore, the exact hot-spot fluid energy equation 

$$
\frac {\partial}{\partial t} \left(\frac {P}{\gamma - 1} + \frac {1}{2} \rho v ^ {2}\right) + \vec {\nabla} \cdot \vec {v} \left(\frac {\gamma P}{\gamma - 1} + \frac {1}{2} \rho v ^ {2}\right) = Q, \tag {2.17}
$$

is simplified into 

$$
\frac {\partial}{\partial t} \frac {3}{2} P + \vec {\nabla} \cdot \frac {5}{2} P \vec {v} = Q \tag {2.18}
$$

by neglecting the fluid kinetic energy density term. Here the ratio of specific heats for a monatomic ideal gas is $\gamma = ( D + 2 ) / D = 5 / 3$ with three degrees of freedom $ \textit { D } = 3$ . The total heat exchange rate $Q = Q _ { \alpha } + Q _ { \mathrm { t h e r m a l } } + Q _ { \mathrm { r a d } }$ between the hot spot and the surrounding plasma is given by the summation of alpha heating rates $Q _ { \alpha }$ , heat-conduction loss rates $\mathrm {  { Q _ { t h e r m a l } } }$ , and radiation loss rates $Q _ { \mathrm { r a d } }$ . A direct volume integration for any scalar quantity $F$ over a 3-D control volume ${ \cal V } _ { \mathrm { b } } ( t )$ with a 3-D control surface $A _ { \mathrm { b } } ( t )$ moving at a control velocity $\vec { v } _ { \mathrm { b } }$ is given by Leibniz integral 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/ecb401e867c1b31e4e52c68f039bdd63ac544ceffdf3d8bcb44cdc493b57e59f.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/aa2d885d9f31c33b9ce427f80eabf0e06148ab3af2735a8b470d90f81427f961.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/2ce7ba704d4cf124f22be36a19d7028f6d8a4c5cc2d110a835dd2c2bce11b0de.jpg)



Figure 2.6: Summary of single-mode energetics at stagnation for the same $1 4 \%$ initial velocity perturbations. Plot (a) is the burn-averaged hot spot and (b) is the e are $T _ { \mathrm { e } } \geq 1$ $\mathrm { I E } _ { \mathrm { s h e l l } } ^ { 3 D } / \mathrm { I E } _ { \mathrm { s h e l l } } ^ { 1 D }$ -keV hot spot. In (a) and (b), the blue dots are and the red dots are $\mathrm { I E _ { h o t s p o t } ^ { 3 D } / I E _ { h o t s p o t } ^ { 1 D } }$ total total. The 3-D Legendre modes $\mathrm { K E _ { t o t a l } ^ { 3 D } } / \mathrm { K E _ { t o t a l } ^ { 1 D } }$ , black dots $\ell = 4 , 6 , 8 , 1 0 , 1 2$ with respect to 1-D at stagnation using thehows the conservation of adiabatic paramet with $m = \ell / 2$ are denoted by `.5 on the $x$ axis. (c) Measurement h t definition withfor all modes. $P _ { \mathrm { 1 k e V } } V _ { \mathrm { 1 k e V } } ^ { 5 / 3 }$ $T _ { \mathrm { e } } \geq 1$ $P V ^ { 5 / 3 }$


rule, 

$$
\frac {\mathrm {d}}{\mathrm {d} t} \int_ {V _ {\mathrm {b}} (t)} F d V = \int_ {V _ {\mathrm {b}} (t)} \frac {\partial F}{\partial t} d V + \int_ {A _ {\mathrm {b}} (t)} F \vec {v} _ {\mathrm {b}} \cdot d \vec {A}, \tag {2.19}
$$

which means that the total rate of change for any scalar quantity within a control volume is a sum of the rate of change by itself and the rate of change due to the moving enclosed surface. Substitute $F = { \frac { 3 } { 2 } } P$ into Eq. (2.19), 

$$
\int_ {V _ {\mathrm {b}} (t)} \frac {\partial}{\partial t} \frac {3}{2} P d V = \frac {\mathrm {d}}{\mathrm {d} t} \int_ {V _ {\mathrm {b}} (t)} \frac {3}{2} P d V - \int_ {A _ {\mathrm {b}} (t)} \frac {3}{2} P \vec {v} _ {\mathrm {b}} \cdot d \vec {A}, \tag {2.20}
$$

and integrate Eq. (2.18) over a burn volume $V _ { \mathrm { b } } ( t )$ enclosed by a burn surface $A _ { \mathrm { b } } ( t )$ moving at a user-defined control velocity ${ \vec { v } } _ { \mathrm { b } } ( t )$ , 

$$
\int_ {V _ {\mathrm {b}} (t)} \frac {\partial}{\partial t} \frac {3}{2} P d V + \int_ {A _ {\mathrm {b}} (t)} \frac {5}{2} P \vec {v} \cdot d \vec {A} = \int_ {V _ {\mathrm {b}} (t)} Q d V. \tag {2.21}
$$

Define the volume-averaged hot spot pressure $\begin{array} { r } { P _ { \mathrm { { b } } } ( t ) \ = \ \int _ { V _ { \mathrm { { b } } } ( t ) } P d V / V _ { \mathrm { { b } } } } \end{array}$ and the volume-averaged total heat exchange rate $\begin{array} { r } { Q _ { \mathrm { b } } ( t ) = \int _ { V _ { \mathrm { b } } ( t ) } Q d V / V _ { \mathrm { b } } } \end{array}$ such that Eq. 

(2.21) becomes, 

$$
\frac {3}{2} \frac {\mathrm {d}}{\mathrm {d} t} \left(P _ {\mathrm {b}} V _ {\mathrm {b}}\right) - \int_ {A _ {\mathrm {b}} (t)} \frac {3}{2} P \vec {v} _ {\mathrm {b}} \cdot d \vec {A} + \int_ {A _ {\mathrm {b}} (t)} \frac {5}{2} P \vec {v} \cdot d \vec {A} = Q _ {\mathrm {b}} V _ {\mathrm {b}}. \tag {2.22}
$$

The $P d V$ work done on the user-defined burn volume 

$$
P _ {\mathrm {b}} \frac {\mathrm {d} V _ {\mathrm {b}}}{\mathrm {d} t} = P _ {\mathrm {b}} \int_ {A _ {\mathrm {b}} (t)} \vec {v} _ {\mathrm {b}} \cdot d \vec {A} \tag {2.23}
$$

is obtained by substituting $\triangle P ( \vec { x } , t ) = P ( \vec { x } , t ) - P _ { \mathrm { b } } ( t )$ into the second term in Eq. (2.22). The hot spot pressure is assumed to be sufficiently uniform in the userdefined burn volume so that the contribution from the $\begin{array} { r } { \int _ { A _ { \mathrm { b } } ( t ) } \triangle P ( \vec { x } , t ) \vec { v } _ { \mathrm { b } } \cdot d \vec { A } } \end{array}$ term in Eq. (2.24) is neglected for the burn-averaged definition for hot-spot pressures. 

$$
- \int_ {A _ {\mathrm {b}} (t)} \frac {3}{2} P \vec {v} _ {\mathrm {b}} \cdot d \vec {A} = P _ {\mathrm {b}} \frac {\mathrm {d} V _ {\mathrm {b}}}{\mathrm {d} t} + \int_ {A _ {\mathrm {b}} (t)} \triangle P (\vec {x}, t) \vec {v} _ {\mathrm {b}} \cdot d \vec {A} - \int_ {A _ {\mathrm {b}} (t)} \frac {5}{2} P \vec {v} _ {\mathrm {b}} \cdot d \vec {A}. \tag {2.24}
$$

The subsonic hot-spot energy equation in Eq. (2.22) becomes, 

$$
\frac {3}{2} \frac {\mathrm {d}}{\mathrm {d} t} \left(P _ {\mathrm {b}} V _ {\mathrm {b}}\right) + P _ {\mathrm {b}} \frac {\mathrm {d} V _ {\mathrm {b}}}{\mathrm {d} t} = Q _ {\mathrm {b}} V _ {\mathrm {b}} - \int_ {A _ {\mathrm {b}} (t)} \frac {5}{2} P \delta \vec {v} \cdot d \vec {A}, \tag {2.25}
$$

where $\delta \vec { v } = \vec { v } - \vec { v } _ { \mathrm { b } }$ . The surface integral measures the flux of enthalpy energy density $h = \gamma P / ( \gamma - 1 ) = \textstyle { \frac { 5 } { 2 } } P$ moving across the burn surface $A _ { \mathrm { b } } ( t )$ , because the control velocity $\vec { v } _ { \mathrm { b } }$ is not the same as the local fluid velocity $\vec { v }$ . Equation (2.25) is the subsonic flow approximation of the hot-spot energy equation derived for burn-averaged quantities $Q _ { \mathrm { b } }$ , which is valid for both low and high modes. In Fig. (2.4) the burn volume $V _ { \mathrm { b } }$ for low and high modes is shown in red. $V _ { \mathrm { b } }$ includes the region of warm bubbles for low modes because of $V _ { \mathrm { 1 k e V } } \simeq V _ { \mathrm { b } }$ but excludes the cold bubbles for high modes. In this work, the volume-averaged quantities $Q \mathrm { { _ { 1 k e V } } }$ are adopted to describe low modes $\ell \leq 6$ , whereas the burn-averaged quantities $Q _ { \mathrm { b } }$ are adopted to describe high modes $\ell \geq 7$ . 

Substitute $P ( \vec { x } , t ) = P _ { \mathrm { b } } + \triangle P ( \vec { x } , t )$ into Eq. (2.25), the hot-spot energy equation can be rewritten as 

$$
\frac {\mathrm {d}}{\mathrm {d} t} \ln \left(P _ {\mathrm {b}} V _ {\mathrm {b}} ^ {5 / 3}\right) = \frac {2 Q _ {\mathrm {b}}}{3 P _ {\mathrm {b}}} - \frac {5}{3 V _ {\mathrm {b}}} \int_ {A _ {\mathrm {b}} (t)} \delta \vec {v} \cdot d \vec {A} - \underbrace {\frac {5}{3 P _ {\mathrm {b}} V _ {\mathrm {b}}} \int_ {A _ {\mathrm {b}} (t)} \triangle P (\vec {x} , t) \delta \vec {v} \cdot d \vec {A}} _ {\text {v a n i s h e s f o r b u r n - a v e r a g e d} P _ {\mathrm {b}}}, \tag {2.26}
$$

and define the hot-spot adiabatic parameter PbV 5/3b , which varies as a function $P _ { \mathrm { b } } V _ { \mathrm { b } } ^ { 5 / 3 }$ of time, depending on the total heat exchange rate $Q _ { \mathrm { b } } ( t )$ , and the enthalpy flux moving across the burn surface $A _ { \mathrm { b } } ( t )$ . The surface integral for $\triangle P$ vanishes in the burn-averaged definition for hot-spot pressures. Simultaneously, the hot-spot entropy increases with time as a result of heat-transferring processes within the hot spot such as radiation, electron and ion heat conduction, and equilibration. 

# 2.4.2 An adiabatic invariant 3-D model for low modes

When the hot-spot surface is defined by $T _ { \mathrm { e } } ~ = ~ 1$ keV, the velocity of the burn surface $\vec { v } _ { \mathrm { b } }$ is about the same as the fluid velocity $\vec { v }$ on the perturbed inner shell surface, i.e., $\delta \vec { v } \simeq 0$ . It leads to the first property of vanishing enthalpy flux in Eqs. (2.25) and (2.26) for both low and high modes. The second property is vanishing $\mathrm {  { Q _ { t h e r m a l } } }$ on the perturbed inner shell surface for both low and high modes because the heat leaving the hot spot is recycled in the form of internal and kinetic energies of the plasma ablated off the hot-spot’s inner shell surface [BUL $^ +$ 01]. Therefore, in the absence of alpha heating and radiation loss $Q _ { \alpha } = Q _ { \mathrm { r a d } } = 0$ , the hotspot adiabatic parameter $P _ { \mathrm { b } } V _ { \mathrm { b } } ^ { 5 / 3 }$ is approximately conserved for both low and high modes [WBY+16] for volume-averaged quantities with Te ≥ 1 keV, which is validated in Fig. (2.6) 

$$
(P _ {\mathrm {3 D}} V _ {\mathrm {3 D}} ^ {5 / 3}) _ {\mathrm {1 k e V}} (t _ {\mathrm {s t a g}}) = (P _ {\mathrm {1 D}} V _ {\mathrm {1 D}} ^ {5 / 3}) _ {\mathrm {1 k e V}} (t _ {\mathrm {s t a g}}), \qquad \qquad (2. 2 7)
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/dd2f83bd4d874de1938422d0990dfd65b458c445c0b1156d8e155be38c08b9f4.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/42d7f335823db15d56dc9debc05f1272dff9e7dc4d0099511c599a1fccba34c5.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/3a7d8bb5172303a61246eda7024a6981e66159ee167fdf09c9f4c6f393067d8b.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/33cb85f5451ae9954271c49ec9ef1aacb193c44915f958ce7bc0c928137facac.jpg)



Figure 2.7: The blue curves indicate the analytic model relations for (a) $\hat { P } =$ $\hat { I E } _ { \mathrm { H S } } ^ { 2 . 5 }$ , (b) $\hat { V } = \hat { I E } _ { \mathrm { H S } } ^ { - 1 . 5 }$ , (c) $\hat { P } = ( 1 - \mathrm { R K E } ) ^ { 2 . 5 }$ , and (d) $\hat { V } = ( 1 - \mathrm { R K E } ) ^ { - 1 . 5 }$ . The shorthand notations are $\hat { Q } = Q _ { \mathrm { 3 D } } / Q _ { \mathrm { 1 D } }$ .


where $t _ { \mathrm { s t a g } }$ is the bang time. In the burn-average definition, however, the hot-spot adiabatic parameter decreases with Legendre mode numbers because of reduction in burn volumes. The conservation of hot-spot adiabatic parameter for low modes leads to the first general 3-D relation that connects the hot-spot pressure and the hot-spot volume in terms of hot-spot internal energetics [WBY+16] 

$$
\frac {P _ {\mathrm {1 k e V}} ^ {\mathrm {3 D}}}{P _ {\mathrm {1 k e V}} ^ {\mathrm {1 D}}} = \left(\frac {\mathrm {I E} _ {\mathrm {1 k e V , H S}} ^ {\mathrm {3 D}}}{\mathrm {I E} _ {\mathrm {1 k e V , H S}} ^ {\mathrm {1 D}}}\right) ^ {5 / 2} \tag {2.28}
$$

and 

$$
\frac {V _ {1 \mathrm {k e V}} ^ {\mathrm {3 D}}}{V _ {1 \mathrm {k e V}} ^ {\mathrm {1 D}}} = \left(\frac {\mathrm {I E} _ {1 \mathrm {k e V} , \mathrm {H S}} ^ {\mathrm {3 D}}}{\mathrm {I E} _ {1 \mathrm {k e V} , \mathrm {H S}} ^ {\mathrm {1 D}}}\right) ^ {- 3 / 2}. \tag {2.29}
$$

Equations (2.28) and (2.29) are validated in Fig. (2.7). When $Q _ { \alpha }$ and $Q _ { \mathrm { r a d } }$ are important in the presence of strong alpha heating or dominant radiation losses, the hot spot adiabatic parameter is not constant in time. A time-dependent integrating factor modifies the adiabatic parameter, 

$$
(P _ {\mathrm {b}} V _ {\mathrm {b}} ^ {5 / 3}) (t _ {\mathrm {s t a g}}) = (P _ {\mathrm {b}} V _ {\mathrm {b}} ^ {5 / 3}) (t _ {0}) \mathrm {E x p} \left(\int_ {t _ {0}} ^ {t _ {\mathrm {s t a g}}} \frac {2 Q _ {\mathrm {b}}}{3 P _ {\mathrm {b}}} d t\right). \qquad (2. 3 0)
$$

The exponent factor appears as a multiplier modifying the expression of the hot spot pressure and the hot spot volume in Eqs. (2.28) and (2.29) only. The following analysis of the total energy conservation and the mass ablation rate remain unchanged. The conservation of total energy provides the second general 3-D relation that connects the hot-spot’s internal energies to the shell and the hot-spot’s kinetic energies at stagnation. The shell’s initial kinetic energy is converted into the hot-spot’s internal energy through $P d V$ work 

$$
\mathrm {I E} _ {\mathrm {1 k e V , H S}} ^ {\mathrm {s t a g}} \simeq \mathrm {K E} _ {\mathrm {t o t}} ^ {\mathrm {m a x}} (t _ {0}) - \mathrm {K E} _ {t o t} ^ {\mathrm {s t a g}} - \mathrm {I E} _ {\mathrm {S H}} ^ {\mathrm {s t a g}}. \tag {2.31}
$$

$\mathrm { K E } _ { \mathrm { t o t } } ^ { \mathrm { m a x } }$ and KE $\begin{array} { l } { \mathrm { s t a g } } \\ { \mathrm { t o t } } \end{array}$ are the total kinetic energy measured in the simulation domain at the beginning of the deceleration phase at time $t _ { 0 }$ and at stagnation $t _ { \mathrm { s t a g } }$ respectively. For shot 77068 the initial total kinetic energy $\mathrm { K E } _ { \mathrm { t o t } } ^ { \mathrm { m a x } } ( t _ { 0 } ) = 1 3 3 0 \ J$ $J$ is large enough to neglect the initial shell and hot-spot internal energies in Eq. (2.31). We define the normalized residual kinetic energy as 

$$
\mathrm {R K E} = \frac {\mathrm {K E} _ {\mathrm {t o t}} ^ {\mathrm {3 D}} (t _ {\mathrm {s t a g}}) - \mathrm {K E} _ {\mathrm {t o t}} ^ {\mathrm {1 D}} (t _ {\mathrm {s t a g}})}{\mathrm {K E} _ {\mathrm {t o t}} ^ {\mathrm {m a x}} (t _ {0})}, \tag {2.32}
$$

Apart from $\mathrm { K E } _ { \mathrm { t o t } } ^ { \mathrm { m a x } } ( t _ { 0 } )$ is the initial value; other variables are measured at stagnation. For simplicity, the subscripts for $t _ { \mathrm { s t a g } }$ and $t _ { 0 }$ are dropped, and divide both sides of Eq. (2.31) with respect to 1-D values: 

$$
\frac {\mathrm {I E} _ {\mathrm {1 k e V , H S}} ^ {\mathrm {3 D}}}{\mathrm {I E} _ {\mathrm {1 k e V , H S}} ^ {\mathrm {1 D}}} = \frac {1 - \hat {K} _ {\mathrm {t o t}} ^ {\mathrm {3 D}} - \hat {I} _ {\mathrm {S H}} ^ {\mathrm {3 D}}}{1 - \hat {K} _ {\mathrm {t o t}} ^ {\mathrm {1 D}} - \hat {I} _ {\mathrm {S H}} ^ {\mathrm {1 D}}}, \tag {2.33}
$$

where the simplified labels are defined by $\hat { K } _ { \mathrm { t o t } } ^ { \mathrm { 3 D / 1 D } } \equiv \mathrm { K E _ { t o t } ^ { \mathrm { 3 D / 1 D } } / K E _ { t o t } ^ { \mathrm { m a x } } }$ ≡ KE3D/tot and ISH ˆ3D/1D $\hat { I } _ { \mathrm { S H } } ^ { \mathrm { 3 D / 1 D } } \equiv$ $\mathrm { I E _ { S H } ^ { 3 D / 1 D } / K E _ { t o t } ^ { m a x } }$ IE3DSH . The right-hand side of Eq. (2.33) can be expanded into 

$$
\frac {1 - \mathrm {R K E} - \mathrm {R I E} _ {\mathrm {S H}} - \left(\hat {K} _ {\text {t o t}} ^ {\mathrm {3 D}} + \hat {I} _ {\mathrm {S H}} ^ {\mathrm {3 D}}\right) \left(\hat {K} _ {\text {t o t}} ^ {\mathrm {1 D}} + \hat {I} _ {\mathrm {S H}} ^ {\mathrm {1 D}}\right)}{1 - \left(\hat {K} _ {\text {t o t}} ^ {\mathrm {1 D}} + \hat {I} _ {\mathrm {S H}} ^ {\mathrm {1 D}}\right) ^ {2}}, \tag {2.34}
$$

where the normalized residual shell internal energy is defined as $\mathrm { R I E _ { S H } } \equiv \hat { I } _ { \mathrm { S H } } ^ { \mathrm { 3 D } } - \hat { I } _ { \mathrm { S H } } ^ { \mathrm { 1 D } }$ . Equation (2.34) can be simplified by retaining only the leading term RKE RIESH and neglecting the quadratic terms leading to 

$$
\frac {\mathrm {I E} _ {1 \mathrm {k e V} , \mathrm {H S}} ^ {3 \mathrm {D}}}{\mathrm {I E} _ {1 \mathrm {k e V} , \mathrm {H S}} ^ {1 \mathrm {D}}} \simeq 1 - \mathrm {R K E}. \tag {2.35}
$$

From the survey of single-mode energetics shown in Fig. (2.6), the change of volume-averaged shell internal energies measured at stagnation for various single modes is significantly less than the change in the total residual kinetic energies that justifies the approximation RKE $\gg$ RIESH. The changes in hot-spot pressure and hot-spot volume can be expressed as a unique function of the total residual kinetic energy by rewriting Eqs. (2.28) and (2.29): 

$$
P _ {1 \mathrm {k e V}} ^ {\mathrm {3 D}} / P _ {1 \mathrm {k e V}} ^ {\mathrm {1 D}} \simeq (1 - \mathrm {R K E}) ^ {5 / 2} \tag {2.36}
$$

and 

$$
V _ {1 \mathrm {k e V}} ^ {\mathrm {3 D}} / V _ {1 \mathrm {k e V}} ^ {\mathrm {1 D}} \simeq (1 - \mathrm {R K E}) ^ {- 3 / 2}. \tag {2.37}
$$

Equations (2.36) and (2.37) provide the fundamental explanation for 3-D hydrodynamic behavior for low modes. The larger hot-spot volume and lower hot-spot pressure observed in low modes are caused by the increasing total residual kinetic energy. Equations (2.36) and (2.37) are validated by $D E C 3 D$ simulations in Fig. (2.7) using the hot-spot definition with $T _ { \mathrm { e } } ~ \geq ~ 1$ keV. Observe that the scalings are not affected significantly by neglecting the quadratic terms $( \hat { K } _ { \mathrm { t o t } } ^ { \mathrm { 1 D } } + \hat { I } _ { \mathrm { S H } } ^ { \mathrm { 1 D } } ) ^ { 2 }$ and $( \hat { K } _ { \mathrm { t o t } } ^ { \mathrm { 3 D } } + \hat { I } _ { \mathrm { S H } } ^ { \mathrm { 3 D } } ) ( \hat { K } _ { \mathrm { t o t } } ^ { \mathrm { 1 D } } + \hat { I } _ { \mathrm { S H } } ^ { \mathrm { 1 D } } )$ in Eq. (2.34). 

The fusion reactivity scales with ion temperatures as a power law $< \sigma v > \sim T _ { \mathrm { i } } ^ { 4 }$ for ion temperatures $T _ { \mathrm { i } } = 1$ to 5 keV, as shown in Fig. (5.3), which leads to a simple estimation of neutron yield $Y \sim n ^ { 2 } < \sigma v > V \tau$ or $Y \sim n ^ { 2 } T _ { \mathrm { i } } ^ { 4 } V \tau$ , where $n$ is the hot-spot ion number density and $\tau$ is the burn width. The yield can be 

expressed in terms of the hot-spot pressure and hot-spot volume $Y \sim P ^ { 4 } V ^ { 3 } M ^ { - 2 } \tau$ by substituting the ideal gas relation $P \sim n T _ { \mathrm { i } }$ and the hot-spot mass $M = \rho V$ . Because the hot-spot pressure $P$ and the hot-spot volume $V$ are related to the hotspot internal energies and residual kinetic energies through adiabatic conditions and energy conservation, the yield scaling $Y \sim P ^ { 4 } V ^ { 3 } M ^ { - 2 } \tau$ describes the general 3-D hot-spot conditions, regardless of linear or nonlinear RT instabilities. The yield’s dependence on $P ^ { 4 } V ^ { 3 } M ^ { - 2 }$ comes from the fusion reactivity scaling $T _ { \mathrm { i } } ^ { 4 }$ , which changes to $T _ { \mathrm { i } } ^ { 2 }$ for fusion plasma [BCBW16] ranging from 6 to 20 keV. The yieldover-clean is approximated as 

$$
\mathrm {Y O C} \simeq \left(\frac {P _ {\mathrm {3 D}}}{P _ {\mathrm {1 D}}}\right) ^ {4} \left(\frac {V _ {\mathrm {3 D}}}{V _ {\mathrm {1 D}}}\right) ^ {3} \left(\frac {M _ {\mathrm {3 D}}}{M _ {\mathrm {1 D}}}\right) ^ {- 2} \left(\frac {\tau_ {\mathrm {3 D}}}{\tau_ {\mathrm {1 D}}}\right). \tag {2.38}
$$

The yield degradation can be shown to be a strong function of the residual kinetic energy of the compressed shell by substituting Eqs. (2.28), (2.29), and (2.35) into $P$ and $V$ terms. The yield-over-clean is further simplified into 

$$
\mathrm {Y O C} \simeq (1 - \mathrm {R K E}) ^ {5. 5} \left(\frac {M _ {\mathrm {3 D}}}{M _ {\mathrm {1 D}}}\right) ^ {- 2} \left(\frac {\tau_ {\mathrm {3 D}}}{\tau_ {\mathrm {1 D}}}\right). \tag {2.39}
$$

The rate of change in the hot-spot mass is given by $\dot { M } = \rho v _ { a } S$ , where $v _ { a }$ is the mass ablation velocity on the inner shell surface. We assume that for low modes, the ablation rate scaling is similar to the predictions of the 1-D theory. Onedimensional approximations [BUL+01] for the mass ablation velocity $\rho v _ { a } \sim T _ { \mathrm { i } } ^ { 5 / 2 } / R$ and the perturbed hot-spot surface area $S \sim R ^ { 2 }$ are used. The total gain of hot-spot mass caused by thermal mass ablation on the inner shell surface over a characteristic time $\tau$ is 

$$
\frac {M}{\tau} \sim T _ {\mathrm {i}} ^ {5 / 2} R, \tag {2.40}
$$

where $\tau$ is the burn duration. A simple scaling for the perturbed hot-spot mass is obtained by substituting the ion temperature $T _ { \mathrm { i } } = m _ { \mathrm { D T } } P V / M$ in terms of $P$ and 

$V$ , where $m _ { \mathrm { D T } }$ is the DT ion mass. 

$$
\frac {M _ {\mathrm {3 D}}}{M _ {\mathrm {1 D}}} \simeq \left(\frac {P _ {\mathrm {3 D}}}{P _ {\mathrm {1 D}}}\right) ^ {5 / 7} \left(\frac {V _ {\mathrm {3 D}}}{V _ {\mathrm {1 D}}}\right) ^ {1 7 / 2 1} \left(\frac {\tau_ {\mathrm {3 D}}}{\tau_ {\mathrm {1 D}}}\right) ^ {2 / 7}. \tag {2.41}
$$

As shown in Fig. (2.8), 1-D approximations for the mass ablation rate and hotspot surface area provide a reasonable estimation for the perturbed hot-spot mass for low modes. Equation (2.41) can be validated using the hot-spot mass scaling relation in Ref. [ZB07]: 

$$
M _ {h} = \left\{\int_ {0} ^ {\tau} P _ {h} \left(t ^ {\prime}\right) ^ {\beta} \left[ R _ {h} \left(t ^ {\prime}\right) ^ {3 \gamma} P _ {h} \left(t ^ {\prime}\right) \right] ^ {(\nu + 1 / 3) / \gamma} d t ^ {\prime} \right\} ^ {1 / (\nu + 1)}, \tag {2.42}
$$

where $t = 0$ is the beginning time of the deceleration-phase, $\gamma = 5 / 3$ for the ideal gas, $\nu = 5 / 2$ for Spitzer thermal conductivity, and $\beta = 4 / 5$ for the self-similar flow solution for the 1-D deceleration phase model in the absence of alpha and radiation transport. The subscript $h$ denotes the volume-averaged hot-spot quantities. By approximating the time integration over a characteristic burn time $\tau$ , the hot-spot mass at stagnation resulting from the mass ablation off the inner shell surface is 

$$
M _ {h} \simeq \left\{\tau P _ {h} ^ {4 / 5} \left[ V _ {h} ^ {\gamma} P _ {h} \right] ^ {1 7 / 1 0} \right\} ^ {2 / 7} = P _ {h} ^ {5 / 7} V _ {h} ^ {1 7 / 2 1} \tau^ {2 / 7}. \tag {2.43}
$$

By eliminating the hot-spot mass term in Eq. (2.39), 

$$
\mathrm {Y O C} \simeq (1 - \mathrm {R K E}) ^ {6 1 / 1 4} \left(\frac {\tau_ {\mathrm {3 D}}}{\tau_ {\mathrm {1 D}}}\right) ^ {3 / 7}. \tag {2.44}
$$

The effect of ablation off the inner shell surface relaxes the dependence of yield degradation on residual kinetic energies from $( 1 - \mathrm { R K E } ) ^ { 5 . 5 }$ to $( 1 - \mathrm { R K E } ) ^ { 4 . 4 }$ , where $6 1 / 1 4 \sim 4 . 4 $ . The burn truncation, about $0 . 9 \le \tau _ { \mathrm { { 3 D } } } / \tau _ { 1 \mathrm { { D } } } \le 1 . 0 5$ , is typically a small effect and can be neglected. Therefore, the YOC is reduced to a strong function of the residual kinetic energy of the compressed shell, and hydrodynamic 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/46b9c0e64aa0d3b092b6180a5c4a860c4a7575912a05ba6400c25da485d655af.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/eaa7a166289d6ba52d3e779b4b92868ba6a7a07e718f54b4dcaced901eae8584.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/b0d89882881f689416383c43ca2b7d93f5c2277cd665797a1e1b61b77f98cf96.jpg)



Figure 2.8: (a) The 1-D scaling relation for the hot-spot mass at stagnation $\hat { M } =$ $\hat { P } ^ { 5 / 7 } \hat { V } ^ { 1 7 / 2 1 } \hat { \tau } ^ { 2 / 7 }$ for low modes $\ell = 1$ to 5 using the hot-spot definition $T _ { \mathrm { e } } ~ \geq ~ 1$ keV. The shorthand notations are $\hat { Q } _ { \mathrm { 3 D } } = Q _ { \mathrm { 3 D } } / Q _ { \mathrm { 1 D } }$ . (b) Comparison of simulated YOC against RKE and the analytic models for low modes $\ell = 1$ to 6 and (c) for high modes $\ell = 7$ to 12. The solid blue curve is the yield degradation model $\mathrm { Y O C } = ( 1 - \mathrm { R K E } ) ^ { 4 . 4 }$ , while the dashed blue curve is yield degradation model $\mathrm { Y O C } = ( 1 - \mathrm { R K E } ) ^ { 5 . 5 } .$ .


instabilities play an important role in causing the yield degradation to scale as 

$$
\mathrm {Y O C} \simeq (1 - \mathrm {R K E}) ^ {4. 4}. \tag {2.45}
$$

$\mathrm { Y O C } ~ = ~ ( 1 - \mathrm { R K E } ) ^ { \mathrm { s . } \mathrm { { o } A } }$ from pure hydrodynamic instabilities without thermal transport and YOC = (1 − RKE)4.4 that includes ablation driven by thermal losses are compared in Fig. (2.8) for low modes $\ell = 1$ to 6 and high modes $\ell \geq 7$ . The yield is shown to decrease monotonically with residual kinetic energies for both low and high modes when the hot spot is defined by $T _ { \mathrm { e } } ~ \geq ~ 1$ keV. Longwavelength perturbations are shown to be well approximated by the adiabatic implosion model. The energetic behaviors of pressure degradation, increasing hotspot volumes and yield degradation between low modes $\ell = 1$ and $\ell = 2$ are indistinguishable because the conservations of energy and the hot-spot adiabatic parameter $P V ^ { 5 / 3 }$ are general properties for arbitrarily distorted hot spot in 3-D. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/4980b269da394cea5a9b1eda8abe7ad9be892d9aa3abe9c8ea9672bdbb83d58e.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/d99f0e332d5f3454e2dc19ec6d78fdb553812cc5d6a5b4c9d2fa19861d782a04.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/ab9127068b455168e70b6d97e567be57a3ea20dc2117103de2afdde318364bbf.jpg)



Figure 2.9: (a) For low mode $\ell = 1$ , the degradation of hot-spot pressure and the increasing hot-spot volume are strong functions of residual kinetic energies. (b) For the high mode $Y _ { \ell = 1 2 } ^ { m = 6 }$ , the burn volume is reduced and the core pressure is increased due to the growth of converging RT spikes. The burn-averaged quantities are used in (a) and (b). (c) The simulated YOC is compared against the hot-spot volume for high modes $\ell = 7$ to 12. The blue curve is the yield degradation model $\mathrm { Y O C } = \hat { V }$ , where $\hat { Q } = Q _ { \mathrm { 3 D } } / Q _ { \mathrm { 1 D } }$ .


# 2.4.3 A non-adiabatic invariant 3-D model for high modes

The hot spot model for high modes requires to use the burn-averaged hot-spot pressure $P _ { \mathrm { b } }$ and volume $V _ { \mathrm { b } }$ . For high modes $\ell \geq 7$ the hot-spot adiabatic parameter PbV 5/3b i $P _ { \mathrm { b } } V _ { \mathrm { b } } ^ { 5 / 3 }$ s not conserved because of reduction in burn volumes. Figure (2.5) shows the decreasing hot-spot volume with electron temperatures. The cold bubbles shown by the 1-keV contour surface do not contribute to a significant fraction of neutron productions compared with the core part shown by the 2.5-keV contour surface. DEC3D single-mode simulation database to validate the yield degradation model for high modes in terms of reduction in burn volume derived in Refs. [KS01] and [CBS $^ +$ 10]: 

$$
\mathrm {Y O C} \simeq \frac {V _ {\mathrm {b}} ^ {\mathrm {3 D}}}{V _ {\mathrm {b}} ^ {\mathrm {1 D}}}. \tag {2.46}
$$

The hot-spot pressure and ion temperature between 3-D and 1-D are assumed equivalent in Ref. [CBS+10]. The 3-D effect in Ref. [CBS+10] is introduced in terms of reduction in burn volumes as a result of the nonlinear growth of RT spikes. On average about 10% variation in hot-spot pressures $P _ { \mathrm { b } } ^ { \mathrm { 3 D } } / P _ { \mathrm { b } } ^ { \mathrm { 1 D } }$ for Legendre modes $\ell \geq 7$ is observed in our simulation database. The slight increase 

of $ { P _ { \mathrm { 3 D } } } /  { P _ { \mathrm { 1 D } } }$ shown in Fig. (2.9) is due to the perfect symmetry of the RT spikes converging to the center. This property is lost for low modes since the pressure degradation is a strong function of RKE. Figure (2.9) validates Eq. (2.46) where the yield degradation for high modes is dominated by the reduction in burn volume in agreement with the 2-D results of Refs. [BBSW17, KS01]. 

To derive the yield degradation model for high modes in Eq. (2.46), the unbounded growth of the hot-spot adiabatic parameter must be included to account for all kinds of heat transfers across the fluid elements and across the burn surface of the 3-D hot spot. The flow correction parameter reported in Ref. [BBSW17] is applied to describe the growth of the hot-spot adiabatic parameter from the beginning of the deceleration phase to stagnations in 1D and 3D are $P _ { \mathrm { b , 1 D / 3 D } } V _ { \mathrm { b , 1 D / 3 D } } ^ { \gamma } ( t _ { \mathrm { s t a g } } ^ { \mathrm { 1 D / 3 D } } ) = P _ { \mathrm { b } } V _ { \mathrm { b } } ^ { \gamma } ( t _ { 0 } ) e ^ { F _ { \mathrm { 1 D / 3 D } } }$ . Introduce the normalized variable $\hat { Q } _ { \mathrm { b } } = Q _ { \mathrm { b } } ^ { \mathrm { 3 D } } / Q _ { \mathrm { b } } ^ { \mathrm { 1 D } }$ to denote the ratio of 3-D to 1-D burn-averaged quantities, 

$$
\hat {P} _ {\mathrm {b}} \hat {V} _ {\mathrm {b}} ^ {\gamma} = e ^ {\triangle F}, \tag {2.47}
$$

where $\triangle F = F _ { \mathrm { 3 D } } - F _ { \mathrm { 1 D } }$ which can be positive nor negative. Therefore, the yield degradation model of $\mathrm { Y O C } ~ = ~ \hat { P } _ { b } ^ { 4 } \hat { V } _ { b } ^ { 3 } \hat { M } _ { b } ^ { - 2 } \hat { \tau }$ in Eq. (2.38) is a function of the normalized hot-spot internal energy $\hat { I } _ { \mathrm { b } } = \hat { P } _ { b } \hat { V } _ { b }$ , which is related to residual kinetic energies $\hat { I } _ { \mathrm { b } } = 1 - \mathrm { R K E }$ through Eq. (2.35) such that 

$$
\hat {P} _ {\mathrm {b}} = \hat {I} _ {\mathrm {b}} ^ {5 / 2} e ^ {- \frac {3}{2} \triangle F}, \tag {2.48}
$$

$$
\hat {V} _ {\mathrm {b}} = \hat {I} _ {\mathrm {b}} ^ {- 3 / 2} e ^ {\frac {3}{2} \triangle F}. \tag {2.49}
$$

The condition of $( e ^ { \triangle F } / \hat { I } _ { \mathrm { b } } ) _ { \mathrm { h i g h } } \ll 1$ is observed in Eq. (2.49) for high modes in order to describe a significant reduction in burn volume for a hot spot filled with cold bubbles, and correspondingly the burn-averaged hot-spot pressure $\hat { P } _ { \mathrm { b } } =$ $( \hat { I } _ { \mathrm { b } } / e ^ { \triangle F } ) ^ { 3 / 2 } \hat { I } _ { \mathrm { b } }$ in Eq. (2.48) starts to rise such as in the secondary-piston effect. 

Substitute Eqs. (2.48)–(2.49) in Eq. (2.38), a generalized YOC model that is valid for low and high modes in the presence of strong alpha is, 

$$
\mathrm {Y O C} = (1 - \mathrm {R K E}) ^ {1 1 / 2} e ^ {- \frac {3}{2} \triangle F} \hat {M} _ {\mathrm {b}} ^ {- 2} \hat {\tau}. \tag {2.50}
$$

Substitute the 3-D hot-spot mass approximation $\hat { M } _ { \mathrm { b } } = \hat { P } _ { \mathrm { b } } ^ { 5 / 7 } \hat { V } _ { \mathrm { b } } ^ { 1 7 / 2 1 } \hat { \tau } ^ { 2 / 7 }$ in above expression, 

$$
\mathrm {Y O C} \simeq \hat {P} _ {\mathrm {b}} ^ {1 8 / 7} \hat {V} _ {\mathrm {b}} ^ {2 9 / 2 1} \hat {\tau} ^ {3 / 7}. \tag {2.51}
$$

For high modes, the burn-averaged hot-spot pressures are approximately unity according to Fig. (2.9)-(b), so that the resulting yield degradation model YOC $\simeq$ $\hat { V } _ { \mathrm { b } } ^ { 2 9 / 2 1 }$ is reduced to be a unique function of the burn-averaged hot-spot volume $\hat { V } _ { \mathrm { b } }$ by neglecting the small variation of 3-D energy confinement times with respect to 1-D’s, in agreement with Eq. (2.46). 

# 2.5 Conclusion

The 3-D radiation-hydrodynamic code DEC3D is used to study the yield degradation caused by Rayleigh-Taylor instabilities in the deceleration phase of inertial confinement fusion. A systematic investigation using DEC3D’s synthetic singlemode database indicates that the yield degradation caused by low- and mid-mode nonuniformities is a strong function of the residual kinetic energy. This result agrees with a simple YOC model assuming the hot spot satisfying adiabatic implosion model and subsonic flow approximation. The dependence of the YOC on residual kinetic energy (RKE) is also in agreement with 2-D HYDRA simulation results from Kritcher et al. The simulated YOC is well approximated by the pure hydrodynamic instability curve $( 1 - \mathrm { R K E } ) ^ { 5 . 5 }$ and the thermal-driven curve (1−RKE)4.4. The degradation of hot-spot pressure and increasing hot-spot volume 

for low modes can be explained in terms of increasing residual kinetic energies. For high modes the yield degradation is dominated by the reduction in burn volumes $\mathrm { Y O C } = V _ { \mathrm { b } } ^ { \mathrm { 3 D } } / V _ { \mathrm { b } } ^ { \mathrm { 1 D } }$ [BBSW17, KS01]. The burn volumes for high modes are significantly reduced due to the growth of cold bubbles. 

# 3 Impact of Hot-Spot Flow Anisotropy

# 3.1 Motivation

In inertial confinement fusion (ICF) implosions, 3.5-MeV alpha particles and 14.1- MeV neutrons [AtV04] are produced by deuterium and tritium (DT) nuclear fusion reactions within the high-temperature, low-density hot spot. Neutrons escape the hot spot and carry essential details about the hot-spot thermal and stagnation flow properties such as inferring hot-spot apparent ion temperatures from the width of neutron energy spectra [Bry73, Mur14, AC11, AC14, Mun16, MFH $^ +$ 17] and inferring hot-spot flow velocities from the shift of mean neutron energies. [GJCF+13, MGF $^ +$ 18] However, implosion performance are degraded due to the growth of Rayleigh–Taylor (RT) instabilities [SB05a, SGC+05], preventing the full conversion of the shell’s kinetic energy into the hot-spot internal energy. [SEH+14, KTB $^ +$ 14] Unconverted kinetic energies remain in form of residual kinetic energies of RT spikes and bubbles at stagnation, leading to the degradation of hot-spot pressures and neutron yields.[WBS $^ +$ 18a, BBSW17] In particular, magnitudes of low-mode hot-spot residual kinetic energies were showed [WCC+15, SEH+14, WBS+18a] large enough to cause significant variations in ion-temperature measurement along different lines of sight (LOS’s). 

In this chapter, we present an analytical model [WBS+18b] to quantify the behavior of ion-temperature measurement asymmetry through a systematic decomposition of the velocity variance term in Eq. (3.5). This new technique provides 

a useful analytic framework to generalize effects of hot-spot flow asymmetry on ion-temperature measurements, and successfully explains the mechanism of how hot-spot flow isotropy influence the ratio of DD to DT apparent ion temperatures, to be described in Chapter 4. 

The organization of this chapter is as follows. Section 3.2 presents a comprehensive qualitative analysis of Doppler velocity broadening of synthetic neutron energy spectra. Section 3.4.1 derives and explains the physical origin of the velocity variance term in Brysk ion temperatures [Bry73, Mur14, AC11, AC14, Mun16, MFH $^ +$ 17] from the first principle of relativistic neutron kinematics. Section 3.4.2 presents our analytic method to generalize the single-mode velocity variance, and is followed by an extension to the multi-mode velocity variance in Section 3.4.3. Section 3.5 derives the analytic relations between hot-spot residual kinetic energy and ion-temperature measurement asymmetry. Section 3.6 is our conclusion for this chapter. 

# 3.2 Ion temperature measurement

# 3.2.1 Properties of non-stagnating hot-spot flows

Figure (3.1) compares hot spot configurations with low and high modes. The core parts for the high mode $Y _ { \ell = 1 2 } ^ { m = 6 }$ and low mode $\ell = 1$ at high temperatures ${ { T } _ { \mathrm { e } } } \ge 4$ keV are shown approximately spherical in shape, so the nonradial components of the shell velocities do not participate in the core compression. In 1-D the shell decelerates as the hot-spot pressure pushes against the shell 

$$
\frac {1}{2} \rho \frac {D v ^ {2}}{D t} = - \vec {v} \cdot \vec {\nabla} P, \tag {3.1}
$$

where $D / D t$ is the material derivative. The term $- | \vec { v } | | \vec { \nabla } P | \cos \theta$ is the driving factor that decelerates the shell and converts the shell’s kinetic energy into the hot-

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/5ef0d18e7eecdad6655961578fc87433c8bc330df613ec29e9df6b91e0398c40.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/4ad12d9124f663977d61d8efc5fed03ef09919d1ae8d0b4257ceaa8c6d7fa558.jpg)



Figure 3.1: The mass density profile for (a) low mode $\ell = 1$ and (b) high mode $Y _ { \ell = 1 2 } ^ { m = 6 }$ at stagnation with electron temperature contours at 1, 2, 3 and 4 keV. Black arrows indicate the fluid velocity field. Since the cores at temperatures $X _ { \mathrm { e } } \geq 4$ keV are approximately spherical in shape, the core compression is contributed mainly by the radial velocity component. For high mode $Y _ { \ell = 1 2 } ^ { m = 6 }$ , the bubbles are cold, characterized by low neutron production rates.


spot internal energy, where $\theta$ is the inclination angle between the pressure gradient vector and the shell’s velocity vector. Velocity field components $\vec { v } _ { \bot }$ perpendicular to $\vec { \nabla } P$ do not participate in the core compression. When the shell’s velocity field is less convergent with respect to 1-D, the hot spot gains less internal energy through the $P d V$ work relation 

$$
\rho D \varepsilon / D t = - P \vec {\nabla} \cdot \vec {v} + Q _ {\text {h e a t}}, \tag {3.2}
$$

where $\varepsilon$ is the hot spot’s internal energy density and $Q _ { \mathrm { h e a t } }$ is the density of total heat sources or sinks. A less convergent perturbed velocity field in the deceleration phase is the key to explain the 3-D effect of decreasing hydrodynamic efficiency to compress the hot core among different single modes. The velocity fields at stagnation for a low mode $\ell = 1$ and a high mode $Y _ { \ell = 1 2 } ^ { m = 6 }$ are compared in Fig. (3.1). 

For high modes, nonstagnating RT spikes contribute to a large fraction of the residual kinetic energy. Figure (3.2) shows the spatial distribution of kinetic energy density for RT spikes for a high mode $Y _ { \ell = 1 2 } ^ { m = 6 }$ . Long-wavelength perturbations 

induce a large spatial scale of high-velocity flows inside the hot spot and more nonradial velocity components in the compressed shell than high modes. Figure (3.3) shows a dominant jet of the low mode $\ell = 1$ at the speed of ∼ 500 km/s flows toward the negative $z$ direction. As shown in Fig. (3.3), the hot core at $T _ { \mathrm { e } } \geq 2 . 3$ keV is advected by the jet with its maximum flow velocity. For low modes, both nonstagnating RT spikes and the high-velocity jets account for the residual kinetic energy at stagnation, as shown by Fig. (3.2). 

In ICF implosions, target offset, stalk, laser beam power inbalance and DT ice roughness are the main sources of mode $\ell = 1$ non-uniformity [SEH $^ +$ 14, IMS $^ + 1 7$ , PMTT14]. The large hot-spot residual kinetic energies associated with low modes lead to significant variations in measuring ion temperatures along different lines of sight (LOS). The signature of implosions with large hot-spot RKE’s are characterized by the shift of the mean neutron energy $\langle E _ { \mathrm { n } } \rangle$ , the broadening of the width of the neutron energy spectrum $\triangle E$ , and the deviation of neutron energy spectra shape from Gaussian distribution. Neutrons born from each fluid element are characterized by a Gaussian energy spectrum centered at the neutron birth energy $E _ { 0 } = 1 4 . 1$ MeV. For a stationary fusion plasma, the width of the spectrum is a function of thermal ion temperature $T _ { \mathrm { i } }$ given by $\bigtriangleup E = \sqrt { 2 m _ { \mathrm { n } } T _ { \mathrm { i } } E _ { 0 } / ( m _ { \mathrm { n } } + m _ { \alpha } ) }$ , where ${ m } _ { \mathrm { n } }$ and $m _ { \alpha }$ are neutron and alpha particle masses, respectively. The effect of collective bulk motion of nonstagnating fusing ions within the hot spot causes a shift $\triangle E _ { \mathrm { f l o w } }$ in the neutron mean energy [Bry73]. 

$$
\triangle E _ {\text {f l o w}} = v _ {\text {f l o w}} \cos \theta \sqrt {2 m _ {\mathrm {n}} E _ {0}}, \tag {3.3}
$$

where $\theta$ is the angle between the spectrometer LOS and the direction of the bulk velocity $v _ { \mathrm { H o w } }$ . For neutrons born from non-stationary fluid elements, the mean neutron energy is $\langle E _ { \mathrm { n } } \rangle = E _ { 0 } + \triangle E _ { \mathrm { f l o w } }$ . The sign of $\triangle E _ { \mathrm { f l o w } }$ depends on the angle $\theta$ between the direction of velocity field vectors and the LOS. From Eq. (3.3), the 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/f0cdcffa817ca7d40bd7db1b063288a709a81db7d9889f06e0046c8e991e181c.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/a487be337d493be6d7a39810b21de501bd1165c253624ec18901e396f83355b7.jpg)



Figure 3.2: The kinetic energy density profiles for (a) low mode $\ell = 1$ and (b) high mode $Y _ { \ell = 1 2 } ^ { m = 6 }$ at stagnation. For low mode $\ell = 1$ , both the RT spike and the jet contribute to residual kinetic energy. For high mode $Y _ { \ell = 1 2 } ^ { m = 6 }$ , the residual kinetic energy is dominated mainly by the nonstagnating RT spikes. The blue spherical outline is the residual kinetic energy of the unshocked part of the shell.


bulk velocity 500 km/s in the mode $\ell = 1$ simulation as shown in Fig. (3.3) causes a maximum mean energy shift of about 271 keV, when the detector is aligned in the same direction as the bulk velocity, i.e., $\theta = 0$ or 2% mean energy shift with respect to $E _ { 0 }$ . 

The mass density, velocity, and ion thermal temperature profiles for mode $\ell = 1$ in Fig. (3.1) are post-processed by IRIS3D, a Monte Carlo-based neutron transport code. The hot-spot ion temperature is inferred from the width $\triangle E _ { \mathrm { D } }$ of Doppler-broadened primary neutron energy spectrum [Bry73] 

$$
T _ {\mathrm {i}} ^ {\text {i n f e r r e d}} = \frac {m _ {\alpha} + m _ {\mathrm {n}}}{1 6 m _ {\mathrm {n}} E _ {0} \ln 2} \triangle E _ {\mathrm {D}} ^ {2}. \tag {3.4}
$$

The effect of Doppler velocity broadening of the width of neutron energy spectrum is [WRF18] 

$$
\sigma_ {\mathrm {n}} ^ {2} = \frac {2 m _ {\mathrm {n}} T _ {\mathrm {i}} E _ {0}}{m _ {\mathrm {n}} + m _ {\alpha}} + 2 m _ {\mathrm {n}} E _ {0} \sigma_ {v} ^ {2}, \tag {3.5}
$$

where $\sigma _ { v } ^ { 2 } = \mathrm { v a r } \left[ \vec { v } \cdot \hat { d } \right]$ is the variance in the component of the fluid velocity along the direction of the detector, denoted by the unit vector $\hat { d }$ , and $\sigma _ { \mathrm { n } } ^ { 2 }$ is the variance of energy for primary neutrons. The nonvanishing contribution to $\vec { v } \cdot \hat { d }$ results 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/765da01809059ac04b745f56b387084148aeb90d626b7e54be4c7e2dabd34f6c.jpg)



Figure 3.3: Velocity magnitude for low mode $\ell = 1$ at stagnation. A jet at 500 km/s is flowing along the negative $z$ direction. The yield-over-clean is 0.74 and the burn-averaged hot-spot velocity magnitude is 166 km/s in this simulation. The red curve indicates the hot core at a 2.3-keV electron temperature, showing that the hot core is advected by the jet with the maximum flow velocity.


from velocity components parallel to $\hat { d }$ . The presence of hot-spot RKE leads to larger inferred ion temperatures than thermal ion temperatures. The effect of Doppler broadening on the width of the neutron energy spectrum is the same for two detectors located at opposite directions. 

The flow effect on inferred ion temperature variations not only depends on the unique flow pattern within the hot spot of each single mode, but also on the spatial distribution of neutron productions. Because the velocity variance is weighted by the burn distribution over space, the relative importance of vortex structure depends on the difference in neutron production rates between the hot core and the bubbles. As shown in Fig. (3.2), high modes do not have explicit jets flowing inside the hot spot and have small spatial scale of vortices located inside the cold bubbles, where neutron production rates are small. Therefore high modes in general have smaller inferred ion temperature variations than low modes. 

DT ion temperatures are inferred by IRIS3D for mode $\ell = 1$ using six detectors along different LOS’s: $- X , + X , - Y , + Y , - Z , + Z$ , as shown in Fig. (3.2.1). Since 

the vortex structure for mode $\ell = 1$ has a rotational symmetry along the $z$ axis, the inferred ion temperatures at the LOS’s at $- X , + X , - Y , + Y$ are about the same. These four detectors are located perpendicular to the jet, resulting in negligible inferred ion temperature variations or zero variance. The LOS’s at $- Z , + Z$ have the largest parallel velocity components and cause about a 1.25-keV inferred ion temperature variation. The blue curve in Fig. (3.2.1) shows a small inferred ion temperature variation for the high mode $Y _ { \ell = 1 2 } ^ { m = 6 }$ . The averaged inferred ion temperature over six LOS’s on the blue curve is 3.67-keV, which is close to the neutron-averaged ion temperature (3.64 keV) measured at that time for high mode m=6 $Y _ { \ell = 1 2 } ^ { m = 6 }$ Y`=12 . 

Vortices of high modes within the cold bubbles do not contribute to significant inferred ion temperature variations because of low neutron production weights in the variance. The high-velocity jet in mode $\ell = 1$ corresponds to the largest inferred ion temperature variation in the mode spectrum. The ion temperature asymmetries in the mode spectrum $\ell = 1$ to 12 are summarized in Fig. (3.2.1), showing a decreasing inferred ion temperature variation with Legendre mode number. It implies that large ion temperature variations observed in experiments indicate the presence of low modes [GJKC+16]. Large ion temperature variations along different LOS’s are caused by large-scale high-velocity flows induced by longwavelength perturbations. 

# 3.2.2 Neutron energy spectrum model

Neutron spectrometry not only infer areal densities through measuring the downscattered neutron yield due to n-D/n-T elastic scatterings occur in the cold shell, but also capable to infer ion temperatures from the shape of primary neutron energy spectra. The thermal velocity broadening centered at DD/DT neutron birth energies is caused by the center-of-mass motion of DD/DT ion-pairs. The 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/ffeb99fdf5906ae9aa96fc2bd482fba7ba1810c867a2c54cc560fe7fd5943529.jpg)



(a)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/a96add29aa91693399d9cdc63889d911c89ace25be19d6c53617e1a68f36f80d.jpg)



(a)



Figure 3.4: (a) The inferred DT ion temperatures at stagnation by IRIS3D using six detectors at different LOS’s at $- X , + X , - Y , + Y , - Z , + Z$ for low mode $\ell = 1$ (red curve) and high mode $Y _ { \ell = 1 2 } ^ { m = 6 }$ (blue curve). (b) The inferred DT ion temperatures at stagnation by $I R I S 3 D$ using 16 detectors at different LOS’s including $- X , + X , - Y , + Y , - Z , + Z$ and other 10 typical nTOF (neutron time of flight) diagnostics on OMEGA.


shape of neutron production spectrum produced by a single fluid element within a stationary fusion plasma is approximately close to a Gaussian distribution with a variance $\sigma _ { \mathrm { B } } ^ { 2 }$ that is proportional to the thermal ion temperature $\mathrm { { T } _ { i } ^ { t h e r m a l } }$ . The subscript $_ \mathrm { B }$ denotes for Brysk, who first derived the non-relativistic Gaussian neutron energy spectrum for a stationary fusion plasma [Bry73]. Exact shape of relativistic primary neutron production spectra for non-stationary fusion plasma without neutron scatterings were first derived by Appelbe [AC11, AC14], and followed by Munro [Mun16, MFH $^ +$ 17]. 

$$
\sigma_ {\mathrm {B}} ^ {2} = \frac {2 m _ {\mathrm {n}} T _ {\mathrm {i}} ^ {\mathrm {t h e r m a l}} E _ {0}}{m _ {\mathrm {n}} + m _ {X}}, \tag {3.6}
$$

where ${ m } _ { \mathrm { n } }$ and $m _ { X }$ are rest masses of neutron and Helium isotope mass which equals to $m s _ { \mathrm { H e } }$ for DD reactions and $m _ { \mathrm { ^ { 4 H e } } }$ or $m _ { \alpha }$ for DT reactions. For nonstationary fusion plasma, center-of-mass DD/DT ion-pair velocities are boosted by the fluid velocity $\vec { v }$ , resulting in a kinematic Doppler shift in the mean neutron 

energy $\mu$ by an amount of ${ \vec { v } } \cdot { \hat { d } } { \sqrt { 2 m _ { \mathrm { n } } E _ { 0 } } }$ such that, 

$$
\mu = E _ {0} + \vec {v} \cdot \hat {d} \sqrt {2 m _ {\mathrm {n}} E _ {0}}. (3. 7)
$$

The overall effect of Doppler velocity broadening of neutron energy spectra is obtained by superposition of Doppler-shifted and burn-weighted neutron energy spectra produced by all non-stagnating fluid elements. Properties (3.6)-(3.7) are applied to generate the synthetic neutron energy spectrum [GSB $^ +$ 14] $f _ { \mathrm { L O S } } ( E _ { \mathrm { n } } )$ along with a given LOS by post-processing DEC3D hydrodynamic data as follows. 

$$
f _ {\mathrm {L O S}} (E _ {\mathrm {n}}) = \sum_ {\mathrm {c e l l}} \underbrace {\frac {Y _ {\mathrm {c e l l}} (\vec {x} , t)}{Y _ {\mathrm {t o t a l}} (t)}} _ {R _ {Y}} \mathrm {E x p} \left[ - \frac {(E _ {\mathrm {n}} - \mu) ^ {2}}{2 \sigma_ {\mathrm {B}} ^ {2}} \right], \tag {3.8}
$$

where $\begin{array} { r } { Y _ { \mathrm { t o t a l } } ( t ) = \sum _ { \mathrm { c e l l } } Y _ { \mathrm { c e l l } } ( \vec { x } , t ) } \end{array}$ is the total neutron yield produced at time $t$ obtained by summing over neutron yield production $Y _ { \mathrm { c e l l } } ( \vec { x } , t )$ of each finite-volume cell. The burn weight $R _ { Y } ( \vec { x } , t ) = Y _ { \mathrm { c e l l } } ( \vec { x } , t ) / Y _ { \mathrm { t o t a l } } ( t )$ associated with a given neutron energy spectrum produced by a fluid element located at $\vec { x }$ is defined as the ratio of neutron yield produced by the fluid element to the total neutron yield at time $t$ . Ion temperatures $T _ { \mathrm { i } } ^ { \mathrm { i n f e r r e d } }$ are inferred from the full width at the half maximum $\triangle E _ { \mathrm { F W H M } }$ for a Gaussian-fitted neutron energy spectrum, 

$$
\triangle E _ {\mathrm {F W H M}} = \sigma_ {\mathrm {B}} \sqrt {8 \ln 2}. \tag {3.9}
$$

By substituting Eq. (3.9) into Eq. (3.6), the neutron-inferred ion temperature obtained by Gaussian-fitted of neutron energy spectra in experiments is, 

$$
T _ {\mathrm {i}} ^ {\text {i n f e r r e d}} = \frac {\triangle E _ {\mathrm {F W H M}} ^ {2} \left(m _ {\mathrm {n}} + m _ {X}\right)}{E _ {0} m _ {\mathrm {n}} 1 6 \ln 2}. \tag {3.10}
$$

To benchmark the neutron energy spectrum model provided by Eq. (3.8), 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/552151dc33d0a780342168102f3fbbc37c921d53d8221064bc5f7e77440d3552.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/193c52bb56358679f31e466559d2c2f35ecafbfc96d564caa840d99cc3f87a59.jpg)



Figure 3.5: (a) DEC3D mass density profile for the single-mode $Y _ { \ell = 1 } ^ { m = 0 }$ at stagnation, simulated by 7% initial velocity perturbation. The electron temperature of 2.4 to 2.45 keV is shown by the red hemisphere located at the origin. The arrows indicate the fluid velocity vectors. A jet with neutron-averaged velocity $\left. v _ { z } \right. = - 1 7 6 ~ \mathrm { k m / s }$ is shown flowing through the central part of the hot spot toward the negative $z$ direction. $\theta$ is the angle measured in the clockwise direction between the positive $z$ axis and LOS. (b) Comparison of the inferred ion temperature measurements from the neutron energy spectrum model in Eq. (3.8) with IRIS3D[WRF18] by post-processing DEC3D mode $\ell = 1$ hydrodynamic data in (a).


neutron-inferred ion temperatures extrapolated from Eq. (3.10) are compared with the result of IRIS3D Monte-Carlo neutron transport simulations [WRF18] at varying LOS angles $\theta$ from north to south poles for the single-mode spectrum $\ell = 1 - 1 2$ . Figure (3.5)-(a) shows $D E C 3 D$ profiles for the mass density and velocity filed vectors for low mode $\ell = 1$ . The neutron-averaged hot-spot flow velocity $\left. v _ { z } \right. = - 1 7 6 ~ \mathrm { k m / s }$ s is large enough to introduce a Doppler shift by a magnitude of $\langle \vec { v } _ { z } \cdot \hat { d } \rangle \sqrt { 2 m _ { \mathrm { n } } E _ { 0 } } = - 9 7 \hat { z } \cdot \hat { d } \mathrm { ~ k e V }$ . The sign depends on the position of the detector in the $\hat { z } \cdot \hat { d }$ term in Eq. (3.7). In OMEGA experiments, typical neutron-inferred hot-spot flow velocities [MGF+18] are about $\sim$ 40 km/s and reach ∼ 100 km/s in shots with large variations in ion temperature measurements. Figure (3.5)-(b) compares neutron-inferred ion temperatures for mode $\ell = 1$ simulation by the neutron energy spectrum model with that by IRIS3D. The Doppler shift effect is shown to vanish when LOS is located at the equator, resulting in $T _ { \mathrm { i , m i n } } ^ { \mathrm { i n f e r r e d } }$ . 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/1c3f39a03c1e2c781571b659c40ca70b926b0ee13a11bfd98a238497236f9581.jpg)



Figure 3.6: Comparison of neutron-inferred ion temperatures between the synthetic neutron energy spectrum model by Eq. (3.8) and IRIS3D by Monte-Carlo simulations of neutron transport.


# 3.3 Neutron energy spectra

# 3.3.1 Single-mode ion-temperature characteristics

The effect of hot-spot flow anisotropy on ion-temperature measurement asymmetry is characterized by the ratio of maximum to minimum neutron-inferred ion temperatures, 

$$
R _ {T} = T _ {\mathrm {i , m a x}} ^ {\text {i n f e r r e d}} / T _ {\mathrm {i , m i n}} ^ {\text {i n f e r r e d}}. \tag {3.11}
$$

Figure (3.6) shows a good agreement for the result of $R _ { T }$ , which is obtained from ion temperatures inferred from north to south poles uniformly distributed over 16 LOS detectors, between synthetic neutron energy spectra by post-processing DEC3D deceleration-phase single-mode database [WBS+18a] through Eqs. (3.6) – (3.10) and IRIS3D, a Monte-Carlo based neutron transport code [WRF18]. The single modes $\ell = 1 - 1 2$ are generated by seeding 1% − 14% initial radial velocity perturbations on the inner shell surface for 1-D profiles of a cryogenic shot 77068 at the beginning of the deceleration phase, which is defined by the time with the maximum 1-D shell implosion velocity. The database contains 2-D single-modes defined by $m = 0$ for all $\ell$ ’s, and 3-D single-modes are chosen as $m = \ell _ { \mathrm { e v e n } } / 2$ only 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/8b05fed4539812074104a321998317d8fd48c491b5d842a104aac6afaa163856.jpg)



Figure 3.7: DEC3D kinetic energy density profile in (a) and $D E C 3 D$ mass density profile in (b) at stagnation for high mode $\ell = 4 0 , m = 2 0$ . The green contour surface is $T _ { \mathrm { e } } = 1$ keV. The vortices of high modes are localized within the cold bubbles, $T _ { \mathrm { e } } ~ < ~ 1$ keV, that do not produce significant amount of neutrons. The Doppler shift term resulting from vortices of large- $\ell$ single-mode perturbations is negligible because of low contribution to the burn distribution.


for even Legendre numbers $\ell _ { \mathrm { e v e n } } = 4$ , $6$ , 8, 10 and 12. 

The advantage of the neutron energy spectrum model in Eq. (3.8) provides a clear physical picture to understand qualitative effects of RT instabilities on neutron-inferred ion temperatures. For mid and high modes, neutron velocities are significantly boosted by high-velocity vorticity circulating inside the cold bubbles, in particular after the onset of nonlinear RT growths [OAK+01, SBRR04]. However, ion-temperature measurement asymmetry is observed to remain weak $R _ { T } \leq 1 . 2$ in Fig. (3.6), because of low neutron production rates within the cold bubbles. As a result, the burn weights $R _ { Y }$ for Doppler-shifted neutron energy spectra produced by fluid elements with high-velocity vorticity are negligible, and suppress their contribution in the final shape of neutron energy spectrum. 

Figure (3.7)-(a) shows an example of high-mode $\ell = 4 0$ and $m = 2 0$ in DEC3D simulation. The residual kinetic energy density of high-velocity vorticity at stagnation is shown highly localized inside the cold bubbles, and vanishes inside the clean volume enclosed by ${ { T } _ { \mathrm { e } } } = 1$ -keV contour surface. Figure (3.7)-(b) is the 3- 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/55dceb8ffaffc22ab81addd604effc25ea68b7508a8d9f816ca686a03a7c986e.jpg)



Figure 3.8: Result of neutron energy spectrum model by post-processing DEC3D hydrodynamic data to compare the Doppler shift of the mean neutron energy and the bulk velocity broadening in (a) for low mode $\ell = 1$ and (b) for high mode $\ell = 4 0 , m = 2 0$ . The unnormalized Doppler-shifted neutron energy spectra sampling inside the cold bubble $f ( E _ { \mathrm { n } } ) _ { \mathrm { b u b b l e } }$ and the high-temperature hot core $f ( E _ { \mathrm { n } } ) _ { \mathrm { c o r e } }$ are compared with the normalized neutron energy spectrum $f _ { \mathrm { L O S } } ( E _ { \mathrm { n } } ) e$ observed at $+ z$ .


D mass density profile, showing the trap of high-velocity vorticity by symmetric converging RT spikes for large- $\ell$ single-mode perturbations. 

Figures (3.8) compares the Doppler shift term $\langle \vec { v } \cdot \hat { d } \rangle \sqrt { 2 m _ { \mathrm { n } } E _ { 0 } }$ in the neutron energy spectrum model between a low mode $\ell = 1$ and a high mode $\ell = 4 0 , m = 2 0$ . For mode $\ell = 1$ , the jet flows along the negative $z$ direction. The sign of $\langle v _ { z } \rangle \hat { z } \cdot \hat { d }$ leads to a downward-shifted red curve in Fig. (3.8)-(a) when a detector is located at the north pole $\hat { d } = \hat { z }$ , and a upward-shifted blue curve when the detector is located at the south pole $\hat { d } = - \hat { z }$ . The black curve with the peak neutron mean energy at $E _ { 0 }$ in Fig. (3.8)-(a) is the neutron energy spectrum without modeling the Doppler shift term in Eq. (3.7), used to compare neutron energy spectra with and without flow effects. 

Although the high-velocity vorticity is observed to produce a large Doppler shift effect in neutron energy spectra for neutrons produced inside the cold bubbles, as shown by the blue curve in Fig. (3.8)-(b), the overall effect of Doppler velocity broadening remains dominated by the large amount of neutrons produced inside the high-temperature clean volume with negligible neutron-averaged hot-spot fluid 

velocities, as shown by the red curve in Fig. (3.8)-(b). As a result, a vanishing Doppler shift effect on the mean neutron energy is observed for large- $\ell$ single modes, which explains the unshifted overall neutron energy spectrum shown by the black curve for mode $\ell = 4 0$ in Fig. (3.8)-(b). 

The property of symmetric converging RT spikes is lost in multimode simulations. The high-velocity vorticity near the base of bubble experience a higher temperature and higher burn weight than that located at the tip of bubble. The complex phenomenon of ion-temperature measurement asymmetries in multi-mode perturbations cannot be qualititatively explained by the neutron energy spectrum model. 

Low modes behave differently, because of their unique feature with approximately uniform physical properties over space such as mass density, pressure and thermal ion temperatures within the 3-D hot spot defined by $T _ { \mathrm { e } } ~ \geq ~ 1$ keV [WBS+18a]. Consequently, the high temperature of the hot core “diffuses”into the interior region of the large bubbles. The warmer bubbles than that for mid and high modes, lead to higher burn weights for long wavelength fluid velocity disturbance within the hot spot, and is the key leading to large ion-temperature measurement asymmetries observed in experiments. 

The signature of single-mode ion-temperature measurement variations is summarized in Fig. (3.9) by a direct comparison of $R _ { T }$ against the single-mode spectrum $\ell = 1 - 1 2$ , simulated by $D E C 3 D$ with the same 7% initial velocity perturbation. T inferred $T _ { \mathrm { i , m a x } } ^ { \mathrm { i n f e r r e d } }$ and $T _ { \mathrm { i , m i n } } ^ { \mathrm { i n f e r r e d } }$ are inferred from LOS uniformly distributed from north to south poles in the neutron energy spectrum model and 16 detectors in IRIS3D. In OMEGA experiments, the averaged ion-temperature ratio is $\sim 1 . 1 8$ and the standard deviation of $T _ { \mathrm { m a x } } ^ { \mathrm { e x p } } / T _ { \mathrm { m i n } } ^ { \mathrm { e x p } }$ is $\sim 0 . 1 4$ . The study of single mode simulations at the 7% initial velocity perturbation level, corresponding to $R _ { T } \sim 1 . 4$ for mode $\ell = 1$ in Fig. (3.9)-(a), is close to observed values in shots with large ion temperature asymmetries exhibiting one standard deviation higher than the 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/66f6561c50cbf4eb3daa8731d81b07461969020b8f49fb8d2f5d3f9611de499c.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/4506565b4a3a1fcf5f44be3c309b4586d1691081cacc1c8cb4602c83d3585603.jpg)



Figure 3.9: Measurement of (a) $R _ { T }$ by the neutron energy spectrum model and IRIS3D. (b) ratios of maximum and minimum neutron-inferred ion temperatures to thermal ion temperatures by the neutron energy spectrum model. 3-D modes $Y _ { \ell = 4 } ^ { m = 2 }$ , $Y _ { \ell = 6 } ^ { m = 3 }$ , $Y _ { \ell = 8 } ^ { m = 4 }$ , $Y _ { \ell = 1 0 } ^ { m = 5 }$ and $Y _ { \ell = 1 2 } ^ { m = 6 }$ are denoted by $\ell = 4 . 5$ , 6.5, 8.5, 10.5 and 12.5 on the $x$ -axis.


averaged experimental ion-temperature ratio. Mode $\ell = 1$ is shown to exhibit the largest ion-temperature variation $T _ { \mathrm { i , m a x } } ^ { \mathrm { i n f e r r e d } } - T _ { \mathrm { i , m i n } } ^ { \mathrm { i n f e r r e d } }$ Ti,min in the single-mode spectrum, and ion-temperature variations are shown to decrease with Legendre mode number. 

Figure (3.9)-(b) studies the unique characteristics of flow structure for different single modes by comparing the maximum with minimum neutron-inferred ion temperatures. Because $T _ { \mathrm { i , m i n } } ^ { \mathrm { i n f e r r e d } }$ reflects the magnitude of isotropic velocity variance, whereas the ion-temperature variation T inferredi,max − T inferredi,min $T _ { \mathrm { i , m a x } } ^ { \mathrm { i n f e r r e d } } - T _ { \mathrm { i , m i n } } ^ { \mathrm { i n f e r r e d } }$ quantifies the magnitude of anisotropic velocity variance. For mode $\ell = 1$ , the highly directional jet contributes to a large Doppler shift in the mean neutron energy, resulting in large anisotropic velocity variance. The decreasing ion-temperature variation with Legendre mode number indicates a trend of growing isotropic velocity variance for mid and high modes. Mode $\ell = 2$ is shown to exhibit a large $T _ { \mathrm { i , m i n } } ^ { \mathrm { i n f e r r e d } }$ caused by the presence of radial flow structure within the large donut-shape warm bubble. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/33d170890e75b1de829b0c34202c2637d72f7c04e9528fadbd7cc6ef2d779669.jpg)



Figure 3.10: DEC3D low mode $\ell = 2$ simulation by 7% initial velocity perturbation. The velocity fields are plotted on the top of mass density profiles on the $x - y$ and $z - x$ planes. The middle is the 3-D electron temperature contour surface to visualize the 3-D configuration of the high-temperature hot core. Non-stagnating hot-spot fluid velocity disturbance driven by the pair of RT spikes along the $z$ -axis and the expanding radial flow structure inside the large donut-shape warm bubble contribute a significant amount of non-translational residual kinetic energies.


# 3.3.2 Skewness and kurtosis for non-Gaussian spectra

Deviations from Gaussian distribution of neutron energy spectra are intensified by the presence of large low-mode hot-spot flow asymmetry. The onset of non-Gaussian neutron energy spectra implies the failure to infer ion temperatures from a Gaussian-fit to the neutron energy spectrum in Eq. (3.10). To capture the non-Gaussian nature, the third order $M _ { 3 }$ (skewness) and the fourth order $M _ { 4 }$ (kurtosis) standardized central moments of neutron energy spectra are measured. Because a normal Gaussian distribution has zero skewness and zero excess kurtosis $( M _ { 4 } - 3 )$ . The $m$ th order standardized central moment is, 

$$
M _ {m} = \frac {\int \left(\frac {E _ {\mathrm {n}} - \langle E _ {\mathrm {n}} \rangle}{\sigma_ {E}}\right) ^ {m} f _ {\mathrm {L O S}} \left(E _ {\mathrm {n}}\right) d E _ {\mathrm {n}}}{\int f _ {\mathrm {L O S}} \left(E _ {\mathrm {n}}\right) d E _ {\mathrm {n}}}, \tag {3.12}
$$

where $\begin{array} { r } { { \langle E _ { \mathrm { n } } \rangle } = { \frac { 1 } { N } } \int E _ { \mathrm { n } } f _ { \mathrm { L O S } } ( E _ { \mathrm { n } } ) d E _ { \mathrm { n } } } \end{array}$ is the mean neutron energy and $\begin{array} { r } { \sigma _ { E } ^ { 2 } = \frac { 1 } { N } \int ( E _ { \mathrm { n } } - } \end{array}$ $\langle E _ { \mathrm { n } } \rangle ) ^ { 2 } f _ { \mathrm { L O S } } ( E _ { \mathrm { n } } ) d E _ { \mathrm { n } }$ is the variance of neutron energy, and $\begin{array} { r } { N = \int f _ { \mathrm { L O S } } ( E _ { \mathrm { n } } ) d E _ { \mathrm { n } } } \end{array}$ is the normalization factor. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/a20585ea2c6ab6956789c140fa6adfffe2d7b37cc0162461358c69ce9ad97c47.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/8e152c9c85186a229690d0d510d6b4195039ac26a0684d443f9d8da33511174d.jpg)



Figure 3.11: (a)The synthetic neutron energy spectrum for mode $\ell = 2$ simulated by 14% initial velocity perturbation for a LOS located on the equator along the positive $x$ -axis. Positive excess kurtosis are pbserved around the tails of the neutron energy spectrum. (b)The synthetic neutron energy spectrum for mode $\ell = 1$ simulated by 7% initial velocity perturbation for a LOS located on the north pole. The downward-flowing jet leads to the formation of negative skewness and negative kurtosis relative to an observer located at the north pole.


Figure (3.10) shows the 3-D electron temperature contour surface, defined by $T _ { \mathrm { e } } = 2 - 2 . 3$ keV for a high burn-weight distorted hot core, for a low mode $\ell = 2$ DEC3D simulation with a seed of 7% initial velocity perturbation. Hot-spot flow velocity fields on $x - y$ and $z - x$ planes are plotted on 2-D mass density profiles to visualize the flow structure driven by RT spikes and bubbles. For a LOS located on the equator along the $x$ -axis, a symmetric radially-outward flow structure exists within the large donut-shape warm bubble, resulting in a pair of parallel and antiparallel Doppler-shifted neutron energy spectra, $R _ { Y } \mathrm { E x p } \left[ - ( E _ { \mathrm { n } } - \mu _ { + } ) ^ { 2 } / \left( 2 \sigma _ { \mathrm { B } } ^ { 2 } \right) \right]$ and $R _ { Y } \mathrm { E x p } \left[ - ( E _ { \mathrm { n } } - \mu _ { - } ) ^ { 2 } / \left( 2 \sigma _ { \mathrm { B } } ^ { 2 } \right) \right]$ respectively, where the Doppler-shifted mean neutron energy is defined as $\mu _ { \pm } = E _ { 0 } \pm v _ { x } \sqrt { 2 E _ { 0 } / m _ { \mathrm { n } } }$ . The burn weight factor $R _ { Y }$ is the same for any fluid element located at rotational symmetric position. After ensemble-averaging of neutron energy spectra contributed from all fluid elements, the final shape of neutron energy spectrum defined by Eq. (3.8) exhibits a positive excess kurtosis for any LOS located on the equator. 

Figure (3.11)-(a) shows the profile of neutron energy spectrum observed by a LOS located on the positive $x -$ axis for mode $\ell = 2$ . Because neutron energy 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/d361b9f2d1f5072874033e9650b515de4dd14f75b22215dddbb59eb91f9d387c.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/29e7ef40f0bc7ecbc45c3914accd375986bde7df39558823fbae0262807250ea.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/a9b3911b4524acbdf2005ff593b24ec38c276bf3592e9c634600b81e483f560f.jpg)



Figure 3.12: Measurement of kurtosis $M _ { 4 }$ for the single-mode spectrum $\ell = 1 - 1 2$ simulated by $1 \% - 1 4 \%$ initial velocity perturbations. For a normal Gaussian distribution, the kurtosis is 3, while the excess kurtosis is zero which is defined by subtracting 3 from kurtosis.


spectra of fluid elements within the donut-shape warm bubble are Doppler-shifted to the left and to the right symmetrically, the blue curve with flow effects is shown higher than the red curve without flow effects near the tails of the spectrum, leading to a positive excess kurtosis. As a result, the effect of symmetric Doppler shift caused by parallel and anti-parallel flows within the large donut-shape warm bubble, leads to zero skewness but positive excess kurtosis of neutron energy spectra for any LOS located on the equatorial plane for mode $\ell = 2$ at all levels of initial velocity perturbations. 

Figure (3.11)-(b) shows the profile of neutron energy spectrum observed by a LOS located on the north pole for mode $\ell = 1$ . The neutron energy spectrum is roughly proportional to the superposition of a Doppler-shifted neutron energy spectrum due to the jet, flowing along the negative $z$ direction, and an un-shifted neutron energy spectrum due to the background, $\sim \mathrm { E x p } \left[ - ( E _ { \mathrm { n } } - \mu _ { + } ) ^ { 2 } / \left( 2 \sigma _ { \mathrm { B } } ^ { 2 } \right) \right] +$ $\mathrm { E x p } \left[ - E _ { \mathrm { n } } ^ { 2 } / \left( 2 \sigma _ { \mathrm { B } } ^ { 2 } \right) \right]$ , where $\mu _ { + } = E _ { 0 } + v _ { z } \sqrt { 2 E _ { 0 } / m _ { \mathrm { n } } }$ is the Doppler-shifted mean neutron energy relative to the LOS located at the positive $z$ -axis. This configuration of neutron energy spectrum has a negative excess kurtosis. 

Figure (3.12) shows the measurements of skewness and kurtosis for the neutron energy spectrum model by Eq. (3.8) in the single-mode spectrum $\ell = 1 - 1 2$ simulated by DEC3D with 1% − 14% initial velocity perturbations for a LOS 

located on the north pole along the positive $z$ -axis and a LOS located on the equatorial plane along the positive $x$ -axis. Figure (3.12)-(a) shows that for a LOS located at the north pole, neutron energy spectra for mode $\ell = 1$ are observed to exhibit the largest negative skewness due to the jet flowing along the negative $z$ direction, and simultaneously exhibit large negative excess kurtosis with $M _ { 4 } =$ $2 . 4 - 3$ shown in Fig. (3.12)-(b). Without flow effects in Fig. (3.12)-(a, b, c), skewness are zero and kurtosis of neutron energy spectra for all modes are about $\sim 3 . 1 - 3 . 2$ , which agrees with the property of Gaussian distribution derived for a stationary fusion plasma [Bry73]. The small deviation from 3 is caused by nonsingle thermal ion temperature distribution over space within the distorted hot spot, and is observed to increase slightly with perturbations. Figure (3.12)-(c) shows that mode $\ell = 2$ exhibits large positive excess kurtosis with $M _ { 4 } = 3 . 1 - 3 . 4$ for the LOS located on the equatorial plane, because of symmetric parallel and anti-parallel flows within the large donut-shape warm bubble. Characteristic of kurtosis are observed not clearly distinguishable for modes $\ell > 3$ in Fig. (3.12)- (b, c). Measurements of directional variation of skewness and excess kurtosis for neutron energy spectra can be used to identify the signature of mode $\ell = 1$ and 2 in experiments. 

# 3.4 A velocity variance decomposition technique

# 3.4.1 Physical origin of velocity variance

The methodology to infer ion temperatures relies on the width of neutron energy spectra, which are governed by DT ion-pair center-of-mass motions and Doppler shift effects caused by non-stationary fusion plasma. To understand the physical origin of the velocity variance term in Eq. (3.5) from the first principle of Newtonian mechanics, the main result of Munro [Mun16, MFH $^ +$ 17] and Appelbe 

[AC11, AC14] relativistic neutron kinematics treatment are reviewed. Only the case for DT is considered, but the same analysis is valid to DD. 

Consider a Lorentz boost of the neutron momentum $p _ { \mathrm { n } }$ in the center-of-mass (CM) frame of a DT ion-pair by its CM frame velocity $v _ { \mathrm { c m } } ^ { \mathrm { D T } }$ relative to the neutron momentum $p _ { \mathrm { n } } ^ { \prime }$ observed in the fluid rest frame. $v _ { \mathrm { c m } } ^ { \mathrm { D T } }$ is the velocity component parallel to the LOS unit vector $\hat { d }$ , where neutrons are detected. 

$$
p _ {\mathrm {n}} ^ {\prime} = \gamma_ {\mathrm {c m}} (p _ {\mathrm {n}} + v _ {\mathrm {c m}} ^ {\mathrm {D T}} E _ {\mathrm {n}} / c ^ {2}), \tag {3.13}
$$

where $E _ { \mathrm { n } } = m _ { \mathrm { n } } c ^ { 2 } + K _ { \mathrm { n } }$ is the total mass-energy of neutron in the CM frame and $c$ is the light speed. ${ m } _ { \mathrm { n } }$ and $K _ { \mathrm { n } }$ are the rest mass and the relativistic kinetic energy of neutron in the CM frame. The Lorentz factor $\gamma _ { \mathrm { c m } } = \left( 1 - \beta _ { \mathrm { c m } } ^ { 2 } \right) ^ { - 1 / 2 }$ is a function of the DT ion-pair CM frame velocity $\beta _ { \mathrm { c m } } = v _ { \mathrm { c m } } / c$ . The kinetic energy of neutron in the CM frame is obtained from the DT nuclear fusion energy release $Q = 1 7 . 6$ MeV [AtV04] and the relative kinetic energy $K$ or the total kinetic energy of DT fusion reactants in their CM frame [AC11]. 

$$
\frac {p _ {\mathrm {n}} ^ {2}}{2 m _ {\mathrm {n}}} = K _ {0} + \frac {\mu}{m _ {\mathrm {n}}} K, \tag {3.14}
$$

where $K _ { 0 } = m _ { \alpha } Q / ( m _ { \mathrm { n } } + m _ { \alpha } ) = m _ { \mathrm { n } } v _ { 0 } ^ { 2 } / 2$ is the neutron birth energy and $\mu =$ $m _ { \mathrm { n } } m _ { \alpha } / ( m _ { \mathrm { n } } + m _ { \alpha } )$ is the reduced mass of DT fusion products. Let $p _ { 0 } ~ = ~ { } ^ { \prime } I l _ { \mathrm { n } } v _ { 0 }$ be the neutron momentum at the zero DT relative kinetic energy limit $K = 0$ and substitute $p _ { 0 } = \sqrt { 2 \mu Q }$ into Eq. (3.14) to expand the neutron momentum $p _ { \mathrm { { n } } } = p _ { 0 } ( 1 + K / Q ) ^ { 1 / 2 }$ with positive relative kinetic energy $K \ > \ 0$ at the first order, 

$$
p _ {\mathrm {n}} \simeq \left(1 + \frac {K}{2 Q}\right) p _ {0}. \tag {3.15}
$$

Substitute Eq. (3.15) into Eq. (3.13) and expand the Lorentz factor $\gamma _ { \mathrm { c m } } \simeq 1 +$ 

$\beta _ { \mathrm { c m } } ^ { 2 } / 2$ at the first order to obtain the neutron momentum in the fluid rest frame, 

$$
p _ {\mathrm {n}} ^ {\prime} = p _ {0} + \frac {K}{2 Q} p _ {0} + m _ {\mathrm {n}} v _ {\mathrm {c m}} ^ {\mathrm {D T}} + \triangle_ {1} + \triangle_ {2}, \tag {3.16}
$$

where $\triangle _ { 1 } = \beta _ { \mathrm { c m } } K _ { \mathrm { n } } / c$ and $\triangle _ { 2 } = \beta _ { \mathrm { c m } } ^ { 2 } ( p _ { \mathrm { n } } + v _ { \mathrm { c m } } ^ { \mathrm { D T } } E _ { \mathrm { n } } / c ^ { 2 } ) / 2$ are two relativistic correction terms. The neutron velocity $p _ { \mathrm { n } } ^ { \prime \prime }$ observed in the laboratory frame along the direction of LOS unit vector $\hat { d }$ is obtained by the second Lorentz boost by the fluid velocity $\vec { v }$ observed in the laboratory frame, 

$$
\vec {p} _ {\mathrm {n}} ^ {\prime \prime} \cdot \hat {v} = \gamma_ {v} \left(\vec {p} _ {\mathrm {n}} ^ {\prime} \cdot \hat {v} + v E _ {\mathrm {n}} ^ {\prime} / c ^ {2}\right), \tag {3.17}
$$

$$
\vec {p} _ {\mathrm {n}} ^ {\prime \prime} \cdot \hat {v} _ {\perp} = \vec {p} _ {\mathrm {n}} ^ {\prime} \cdot \hat {v} _ {\perp}, \tag {3.18}
$$

where $E _ { \mathrm { n } } ^ { \prime } = \gamma _ { v } m _ { \mathrm { n } } c ^ { 2 }$ is the total mass-energy of neutron in the fluid rest frame, $\gamma _ { v } = \big ( 1 - \beta _ { v } ^ { 2 } \big ) ^ { - 1 / 2 }$ is the Lorentz factor as a function of fluid velocity $\beta _ { v } = v / c$ , and $\hat { v } _ { \perp }$ is a unit vector perpendicular to the direction of fluid velocity defined by $\hat { v } = \vec { v } / v$ . Only the component of neutron momentum $\vec { p _ { \mathrm { n } } ^ { \prime \prime } } \cdot \hat { v }$ parallel to the fluid velocity vector is Lorentz boosted. Expand the Lorentz factor $\gamma _ { v } \simeq 1 + \beta _ { v } ^ { 2 } / 2$ at the first order and add two equations above to obtain the neutron momentum vector $\vec { p } _ { \mathrm { n } } ^ { \prime \prime } = ( \vec { p } _ { \mathrm { n } } ^ { \prime \prime } \cdot \hat { v } ) \hat { v } + ( \vec { p } _ { \mathrm { n } } ^ { \prime \prime } \cdot \hat { v } _ { \perp } ) \hat { v } _ { \perp }$ in the laboratory frame, and use $\begin{array} { r } { E _ { \mathrm { n } } ^ { \prime } \simeq ( 1 + \frac { 1 } { 2 } \beta _ { v } ^ { 2 } ) m _ { \mathrm { n } } c ^ { 2 } } \end{array}$ . 

$$
\vec {p} _ {\mathrm {n}} ^ {\prime \prime} = \vec {p} _ {\mathrm {n}} ^ {\prime} + m _ {\mathrm {n}} \vec {v} + \frac {\beta_ {v} ^ {2}}{2} \left[ m _ {\mathrm {n}} v + \vec {p} _ {\mathrm {n}} ^ {\prime} \cdot \hat {v} + v E _ {\mathrm {n}} ^ {\prime} / c ^ {2} \right] \hat {v}. \tag {3.19}
$$

Magnitudes of neutron momenta $p _ { \mathrm { n } } ^ { \prime \prime } = \vec { p } _ { \mathrm { n } } ^ { \prime \prime } \cdot \hat { d }$ and $p _ { \mathrm { n } } ^ { \prime } = \vec { p } _ { \mathrm { n } } ^ { \prime } \cdot \hat { d }$ parallel to the LOS unit vector $\hat { d }$ are obtained by taking a dot product on both sides of Eq. (3.19), 

$$
p _ {\mathrm {n}} ^ {\prime \prime} = p _ {\mathrm {n}} ^ {\prime} + m _ {\mathrm {n}} \vec {v} \cdot \hat {d} + \frac {\beta_ {v} ^ {2}}{2} \left[ m _ {\mathrm {n}} v + \vec {p} _ {\mathrm {n}} ^ {\prime} \cdot \hat {v} + v E _ {\mathrm {n}} ^ {\prime} / c ^ {2} \right] \hat {v} \cdot \hat {d}. \tag {3.20}
$$

In the non-relativistic limit $\beta _ { v }  0$ , Eq. (3.20) is reduced to a simple momentum addition $p _ { \mathrm { n } } ^ { \prime \prime } \ = \ p _ { \mathrm { n } } ^ { \prime } + m _ { \mathrm { n } } { \vec { v } } \cdot { \hat { d } }$ , which implies that the neutron velocity 

$p _ { \mathrm { n } } ^ { \prime \prime } = m _ { \mathrm { n } } v _ { \mathrm { n } } ^ { \prime \prime }$ detected along with a given LOS unit vector $\hat { d }$ is the sum of the neutron birth velocity $v _ { 0 } = \sqrt { 2 K _ { 0 } / m _ { 0 } }$ , the DT center-of-mass velocity $v _ { \mathrm { c m } } ^ { \mathrm { D T } }$ in the fluid rest frame that contributes to the thermal ion temperature, a small velocity shift $\kappa = v _ { 0 } K / ( 2 Q ) > 0$ due to DT relative kinetic energy, and a Doppler velocity shift term $\vec { v } \cdot \hat { d }$ due to a non-stationary fusion plasma. The final expression of non-relativistic neutron momentum observed in the laboratory frame is obtained by substituting Eq. (3.16) with $\triangle _ { 1 } = \triangle _ { 2 }  0$ into Eq. (3.20) with $\beta _ { v }  0$ , 

$$
v _ {\mathrm {n}} ^ {\prime \prime} = v _ {0} + v _ {\mathrm {c m}} ^ {\mathrm {D T}} + \kappa + \vec {v} \cdot \hat {d}. \tag {3.21}
$$

The burn-averaged mean neutron velocity is, 

$$
\left\langle v _ {\mathrm {n}} ^ {\prime \prime} \right\rangle = v _ {0} + \left\langle \kappa \right\rangle + \left\langle \vec {v} \cdot \hat {d} \right\rangle , \tag {3.22}
$$

where the mean fluctuation of DT center-of-mass velocity is zero $\langle v _ { \mathrm { c m } } ^ { \mathrm { D T } } \rangle = 0$ , because $v _ { \mathrm { c m } } ^ { \mathrm { D T } }$ is isotropic in space. The burn-averaged bracket is defined as $\langle f \rangle =$ $\int f d N / \int d N$ , where the burn distribution $d N ( T _ { \mathrm { i } } ^ { \mathrm { t h e r m a l } } , \vec { v } )$ is a four-dimensional array that contains thermal ion temperatures and fluid velocities $v _ { x } , v _ { y }$ and $v _ { z }$ in the laboratory frame. The total number of neutrons $N _ { \mathrm { t o t a l } } = \int d N$ is obtained by integrating over the burn distribution. The mean velocity shift due to DT relative kinetic energy, 

$$
\langle \kappa \rangle = v _ {0} \frac {\langle K \rangle}{2 Q}, \tag {3.23}
$$

can be rewritten into Munro’s expression [Mun16, MFH+17] by substituting $Q =$ $p _ { 0 } ^ { 2 } / ( 2 \mu )$ and approximating the reduced mass $\mu = { m _ { \mathrm { n } } m _ { \alpha } } / { ( m _ { \mathrm { n } } + m _ { \alpha } ) }$ for DT fusion products, 

$$
\frac {m _ {\alpha}}{m _ {\mathrm {n}} + m _ {\alpha}} \simeq 1 - \frac {E _ {0}}{M _ {\mathrm {D T}} c ^ {2}}, \tag {3.24}
$$

where $M _ { \mathrm { D T } } = m _ { \mathrm { D } } + m _ { \mathrm { T } }$ is the total DT reactant mass and $E _ { 0 } = m _ { \mathrm { n } } c ^ { 2 } + K _ { 0 }$ 

is the total mass-energy of neutron at the zero DT relative kinetic energy. Two underlying assumptions are $M _ { \mathrm { D T } } \simeq m _ { \mathrm { n } } + m _ { \alpha }$ and $m _ { \mathrm { n } } c ^ { 2 } \gg K _ { 0 }$ . Substituting Eq. (3.24) into Eq. (3.25), 

$$
\langle \kappa \rangle_ {\mathrm {M u n r o}} = \left(1 - \frac {E _ {0}}{M _ {\mathrm {D T}} c ^ {2}}\right) \frac {\langle K \rangle}{p _ {0}}, \tag {3.25}
$$

where the mean DT relative kinetic energy $\langle K \rangle = 5 T _ { \mathrm { i } } ^ { \mathrm { t } }$ hermal [Bry73] is a function of the thermal ion temperature. The value of the bracket is about $[ 1 - E _ { 0 } / ( M _ { \mathrm { D T } } c ^ { 2 } ) ] \sim$ 0.8, and the mean neutron velocity shift due to DT relative kinetic energy is $\langle \kappa \rangle _ { \mathrm { M u n r o } } = 1 . 4 7 ~ \mathrm { k m / s } \cdot \mathrm { k e V ^ { - 1 } } \times 5 T _ { \mathrm { i } } ^ { \mathrm { t h } }$ hermal [Mun16, MFH $^ +$ 17], which is small for typical thermal ion temperatures in cryogenic implosions $T _ { \mathrm { i } } ^ { \mathrm { t h e r m a l } } < 5 ~ \mathrm { k } \epsilon$ eV. 

Introduce a velocity variable $w \equiv v _ { \mathrm { n } } ^ { \prime \prime } - v _ { 0 }$ relative to the neutron birth velocity, which is parallel to LOS, and take the variance of $w$ in Eq. (3.21). 

$$
\begin{array}{l} \operatorname {v a r} [ w ] = \operatorname {v a r} \left[ v _ {\mathrm {c m}} ^ {\mathrm {D T}} \right] + \operatorname {v a r} [ \kappa ] + \operatorname {v a r} \left[ \vec {v} \cdot \hat {d} \right] (3.26) \\ + 2 \operatorname {c o v} \left[ v _ {\mathrm {c m}} ^ {\mathrm {D T}} \kappa \right] + 2 \operatorname {c o v} \left[ \kappa (\vec {v} \cdot \hat {d}) \right] + 2 \operatorname {c o v} \left[ (\vec {v} \cdot \hat {d}) v _ {\mathrm {c m}} ^ {\mathrm {D T}} \right], (3.27) \\ \end{array}
$$

where the covariance for any two scalar functions $f$ and $g$ is defined as $\operatorname { c o v } [ f g ] =$ $\langle f g \rangle - \langle f \rangle \langle g \rangle$ , and the variance for any scalar function $f$ is defined as $\operatorname { v a r } \left[ f \right] =$ $\langle f ^ { 2 } \rangle - \langle f \rangle ^ { 2 }$ , using the burn-average bracket defined previously. Equation (3.26) is a general expression for the neutron velocity variance in the non-relativistic limit. Three covariance terms approach to zero in the limit when DT center-of-mass motions $v _ { \mathrm { c m } } ^ { \mathrm { D T } }$ , small thermal shifts $\kappa$ and non-stagnation hot-spot fluid motions $\vec { v } \cdot \hat { d }$ are independent variables to give zero correlation with each other. Multiply both sides of Eq. (3.26) with the DT total product mass $( m _ { \mathrm { n } } + m _ { \alpha } )$ in the limit of zero covariance terms to obtain the neutron-inferred ion temperature $T _ { \mathrm { i } } ^ { \mathrm { i n f e r r e d } }$ and 

the burn-averaged thermal ion temperature $\mathrm { { T } _ { i } ^ { t h e r m a l } }$ , 

$$
T _ {\mathrm {i}} ^ {\text {i n f e r r e d}} = \left(m _ {\mathrm {n}} + m _ {\alpha}\right) \operatorname {v a r} [ w ], \tag {3.28}
$$

$$
T _ {\mathrm {i}} ^ {\mathrm {t h e r m a l}} = \left(m _ {\mathrm {n}} + m _ {\alpha}\right) \operatorname {v a r} \left[ v _ {\mathrm {c m}} ^ {\mathrm {D T}} \right] + \left(m _ {\mathrm {n}} + m _ {\alpha}\right) \operatorname {v a r} \left[ \kappa \right], \tag {3.29}
$$

such that 

$$
T _ {\mathrm {i}} ^ {\text {i n f e r r e d}} = T _ {\mathrm {i}} ^ {\text {t h e r m a l}} + \left(m _ {\mathrm {n}} + m _ {\alpha}\right) \operatorname {v a r} \left[ \vec {v} \cdot \hat {d} \right]. \tag {3.30}
$$

The second line of Eq. (3.28) means that the thermal velocity broadening of the neutron energy spectrum is a sum of isotropic center-of-mass motions of DT ionpairs and DT ions relative motions in their CM frame. The second term in Eq. (3.30) is the velocity variance that measures the Doppler velocity broadening due to flow effects. The transformation from the energy variance in Eq. (3.6) into the velocity variance is obtained by substituting $\sigma _ { \mathrm { B } } ^ { 2 } = \left. \delta E _ { \mathrm { n } } ^ { 2 } \right. - \langle \delta E _ { \mathrm { n } } \rangle ^ { 2 } = \mathrm { v a r } \left[ E _ { \mathrm { n } } \right]$ using the linearized relation $\delta E _ { \mathrm { n } } = \sqrt { 2 m _ { \mathrm { n } } E _ { 0 } } \delta v _ { \mathrm { n } }$ . 

$$
\operatorname {v a r} \left[ E _ {\mathrm {n}} \right] = 2 m _ {\mathrm {n}} E _ {0} \cdot \operatorname {v a r} [ w ] \tag {3.31}
$$

Brysk thermal ion temperatures T thermali can be recovered from above equation by substituting $\sigma _ { \mathrm { B } } ^ { 2 } = \mathrm { v a r } \left[ E _ { \mathrm { n } } \right]$ using Eq. (3.6) and taking a zero flow effect limit var $[ \vec { v } \cdot \hat { d } ]  0$ in Eq. (3.30). 

The traditional method to infer ion temperatures from Gaussian-fitted neutron energy spectra is only valid when flow effects are small. In the presence of large hot-spot flow asymmetry, the superposition of Gaussian neutron energy spectra from all non-stagnating fluid elements by Eq. (3.8) is shown to result in non-Gaussian neutron energy spectra in Section 3.3.2. Consequently, the methodology of Gaussian-fitted neutron-inferred ion temperatures can result in poorer correlation with experimental yields when neutron energy spectra exhibit large 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/63b43683825da0518391eba5c4de62004d544f20ffac9af95c7fb5d35a24bea3.jpg)



Figure 3.13: Comparison of inferred ion temperature ratio $R _ { T }$ between the Brysk temperatures by Eq. (3.30) through measuring the velocity variance and the synthetic neutron energy spectrum through (a) inferring ion temperatures from the FWHM of the neutron energy spectrum $f _ { \mathrm { L O S } } ( E _ { \mathrm { n } } )$ and (b) inferring ion temperatures from the velocity variance of the neutron velocity spectrum $f _ { L O S } ( v _ { \mathrm { n } } )$ .


non-zero skewness and excess kurtosis. By inferring ion temperatures from the energy variance in Eq. (3.31) or the velocity variance of neutron velocity spectra $f _ { \mathrm { L O S } } \left( v _ { \mathrm { n } } = \sqrt { 2 E _ { \mathrm { n } } / m _ { \mathrm { n } } } \right)$ , a robust agreement of neutron-inferred ion temperatures was observed in Fig. (3.13) between the direction computation of the hot-spot fluid velocity variance by post-processing DEC3D hydrodynamic data by Eq. (3.30) and the measurement of velocity variance from neutron velocity spectra by Eq. (3.32), 

$$
T _ {\mathrm {i}} ^ {\mathrm {i n f e r r e d}} = (m _ {\mathrm {n}} + m _ {\alpha}) \frac {\int \left(v _ {\mathrm {n}} - \langle v _ {\mathrm {n}} \rangle\right) ^ {2} f _ {\mathrm {L O S}} (v _ {\mathrm {n}}) d v _ {\mathrm {n}}}{\int f _ {\mathrm {L O S}} (v _ {\mathrm {n}}) d v _ {\mathrm {n}}}. \tag {3.32}
$$

# 3.4.2 Single-mode velocity variance

The neutron energy spectrum model in Eq. (3.8) only provides qualitative explanation for ion-temperature measurement asymmetry. The model is too simple to explain the complicated multi-mode ion-temperature measurement variations. To quantify effects of hot-spot flow asymmetry on single-mode and multi-mode ion-temperature measurement asymmetry, physical properties of velocity variance in Eq. (3.30) are comprehensively analyzed in this section through a systematic vector decomposition of the $\vec { v } \cdot \hat { d }$ term. 

It is convenient to define the variable of normalized ion temperature to absorb the fusion product mass term for compact representation purpose. 

$$
\hat {T} _ {\mathrm {i}} = T _ {\mathrm {i}} / \left(m _ {\mathrm {n}} + m _ {X}\right), \tag {3.33}
$$

where $X = \mathrm { _ { 2 } ^ { 3 } H e }$ for DD, and $X = \mathrm { { _ { 2 } ^ { 4 } H e } }$ or $\alpha$ for DT fusion reactions. Equation (3.30) in terms of normalized temperatures is 

$$
\hat {T} _ {\mathrm {i}} ^ {\text {i n f e r r e d}} = \hat {T} _ {\mathrm {i}} ^ {\text {t h e r m a l}} + \operatorname {v a r} \left[ \vec {v} \cdot \hat {d} \right], \tag {3.34}
$$

where the velocity variance is, 

$$
\operatorname {v a r} \left[ \vec {v} \cdot \hat {d} \right] = \langle (\vec {v} \cdot \hat {d}) ^ {2} \rangle - \langle \vec {v} \cdot \hat {d} \rangle^ {2}. \tag {3.35}
$$

The hot-spot fluid velocity vector $\vec { v }$ parallel to the LOS unit vector $\hat { d }$ is decomposed into components along three orthogonal Cartesian axes. Expand the inner product, 

$$
\vec {v} \cdot \hat {d} = \left(v _ {x}, v _ {y}, v _ {z}\right) \cdot \left(g _ {x}, g _ {y}, g _ {z}\right) = v _ {x} g _ {x} + v _ {y} g _ {y} + v _ {z} g _ {z}, \tag {3.36}
$$

and the square of the inner product, 

$$
(\vec {v} \cdot \hat {d}) ^ {2} = v _ {x} ^ {2} g _ {x} ^ {2} + v _ {y} ^ {2} g _ {y} ^ {2} + v _ {z} ^ {2} g _ {z} ^ {2} + v _ {x} v _ {y} 2 g _ {x} g _ {y} + v _ {y} v _ {z} 2 g _ {y} g _ {z} + v _ {z} v _ {x} 2 g _ {z} g _ {x}. \tag {3.37}
$$

The three geometrical factors that define the LOS unit vector $\hat { d }$ is, 

$$
g _ {x} = \sin \theta \cos \phi , \tag {3.38}
$$

$$
g _ {y} = \sin \theta \sin \phi , \tag {3.39}
$$

$$
g _ {z} = \cos \theta . \tag {3.40}
$$

By substituting vector decomposition in Eqs. (3.36) and (3.37) into neutron-

averaged brackets in Eq. (3.35), 

$$
\operatorname {v a r} \left[ \vec {v} \cdot \hat {d} \right] = \sum_ {i = 1} ^ {3} \sigma_ {i i} g _ {i} g _ {i} + \sum_ {i \neq j} \sigma_ {i j} g _ {i} g _ {j}. \tag {3.41}
$$

The factor of 2 is absorbed by the symmetric summation in the second term of Eq. (3.41). The summation indices that define Cartesian coordinates $1 = x$ , $2 = y$ , and $3 = z$ will be used interchangeably. Shorthand notations for variance $\sigma _ { i i }$ (with $i = j$ ) and covariance $\sigma _ { i j }$ (with $i \neq j$ ) are defined as, 

$$
\sigma_ {i i} = \left\langle v _ {i} ^ {2} \right\rangle - \left\langle v _ {i} \right\rangle^ {2}, \tag {3.42}
$$

$$
\sigma_ {i j} = \langle v _ {i} v _ {j} \rangle - \langle v _ {i} \rangle \langle v _ {j} \rangle . (3. 4 3)
$$

Therefore, Brysk neutron-inferred ion temperatures in Eq. (3.34) are well defined by a complete set of six hot-spot flow parameters in terms of three directional variances ( $\sigma _ { x x }$ , $\sigma _ { y y }$ , and $\sigma _ { z z }$ ) and three covariances ( $\sigma _ { x y }$ , $\sigma _ { y z }$ , and $\sigma _ { z x }$ ). LOS effects are manifested through three geometrical factors $g _ { x }$ , $g _ { y }$ , and $g _ { z }$ . The terminology of “directional variance”, for variables $\sigma _ { x x }$ , $\sigma _ { y y }$ and $\sigma _ { z z }$ , are used to distinguish themselves from the velocity variance var $\left[ { \vec { v } } \cdot { \hat { d } } \right]$ defined in Eq. (3.35) before the vector decomposition. The compact form of non-relativistic Brysk ion temperature in Eq. (3.34) is, 

$$
\hat {T} _ {\mathrm {i}} ^ {\text {i n f e r r e d}} = \hat {T} _ {\mathrm {i}} ^ {\text {t h e r m a l}} + \sum_ {i, j} \sigma_ {i j} g _ {i} g _ {j}. \tag {3.44}
$$

Figure (3.14) compares the ion temperature ratio $R _ { T } \ = \ T _ { \mathrm { i , m a x } } ^ { \mathrm { i n f e r r e d } } / T _ { \mathrm { i , m i n } } ^ { \mathrm { i n f e r r e d } }$ simulated by IRIS3D with Brysk ion temperatures through a direct computation of directional-variance and covariance defined by Eq. (3.44). $R _ { T }$ for midmodes $\ell = 5 - 1 2$ is shown to be smaller than that for low modes $\ell = 1 - 4$ due to smaller magnitudes of neutron-averaged hot-spot fluid velocities. When hot spots are significantly distorted, neutron energy spectra are no longer in Gaussian shape. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/abe0ec93273adacaf1e04fc6e874d7a86920db34dfd242ade31b80a9476d13fd.jpg)



Figure 3.14: Comparison of maximum to minimum neutron-inferred ion temperatures between IRIS3D and the Brysk ion temperatures to validate Eq. (3.44) using DEC3D single-mode database $\ell = 1 - 1 2$ with different initial velocity perturbations $\delta v / v _ { 0 } = 0 . 0 1 - 0 . 1 4$ .


An increasing discrepancy is observed between Brysk ion temperatures defined in Eq. (3.44) and ion temperatures inferred from the Gaussian-fit method in IRIS3D. At $R _ { T } = 1 . 4$ , which is comparable to OMEGA experiments with large ion temperature variations, discrepancy between ion temperatures inferred from the velocity variance measurement and the Gaussian-fit method for mode $\ell = 4$ is about 5% higher than the $Y = X$ curve in Fig. 3.14. When neutron energy spectra are non-Gaussian, an accurate ion temperature is obtained from a direct measurement of energy variance or velocity variance from the neutron production spectrum outlined by Eq. (3.32). The discrepancy between Brysk temperatures and ion temperatures inferred from the variance measurement method is observed to vanish. 

The neutron-inferred ion-temperature measurement asymmetry is a collective effect of six hot-spot flow parameters through a linear superposition relation, 

$$
\triangle \hat {T} _ {\mathrm {i}} = \hat {T} _ {\mathrm {i}} ^ {\text {i n f e r r e d}} - \hat {T} _ {\mathrm {i}} ^ {\text {t h e r m a l}} = \sum_ {i, j} \sigma_ {i j} g _ {i} g _ {j}. \tag {3.45}
$$

Equation (3.45) can be generalized into a system of linear equations that connects the state vector containing the complete set of six hot-spot flow parameters $\vec { \sigma } = ( \sigma _ { x x } , \sigma _ { y y } , \sigma _ { z z } , 2 \sigma _ { x y } , 2 \sigma _ { y z } , 2 \sigma _ { z x } )$ with the ion-temperature measurement variation vector $\triangle \vec { T } _ { \mathrm { L O S } }$ observed at six LOS’s through a $6 \times 6$ invertible line-of-sight matrix $\hat { M } _ { \mathrm { L O S } } ^ { \mathrm { 3 D } }$ , which is determined by six pairs of LOS angles $\theta _ { \mathrm { L O S } }$ and $\phi _ { \mathrm { L O S } }$ , 

$$
\left[ \begin{array}{l} \triangle \hat {T} _ {\mathrm {i}, 1} \\ \triangle \hat {T} _ {\mathrm {i}, 2} \\ \triangle \hat {T} _ {\mathrm {i}, 3} \\ \triangle \hat {T} _ {\mathrm {i}, 4} \\ \triangle \hat {T} _ {\mathrm {i}, 5} \\ \triangle \hat {T} _ {\mathrm {i}, 6} \end{array} \right] = \left[ \begin{array}{c c c c c c} g _ {x, 1} g _ {x, 1} & . & . & . & . & g _ {z, 1} g _ {x, 1} \\ g _ {x, 2} g _ {x, 2} & . & . & . & . & g _ {z, 2} g _ {x, 2} \\ g _ {x, 3} g _ {x, 3} & . & . & . & . & g _ {z, 3} g _ {x, 3} \\ g _ {x, 4} g _ {x, 4} & . & . & . & . & g _ {z, 4} g _ {x, 4} \\ g _ {x, 5} g _ {x, 5} & . & . & . & . & g _ {z, 5} g _ {x, 5} \\ g _ {x, 6} g _ {x, 6} & . & . & . & . & g _ {z, 6} g _ {x, 6} \end{array} \right] _ {\triangle \vec {T} _ {\text {L O S}}} \left[ \begin{array}{l} \sigma_ {x x} \\ \sigma_ {y y} \\ \sigma_ {z z} \\ 2 \sigma_ {x y} \\ 2 \sigma_ {y z} \\ 2 \sigma_ {z x} \end{array} \right] \cdot \tag {3.46}
$$

For two-dimensional (2-D) planar flows with a translational symmetry in the $z$ direction, $\sigma _ { z z } = \sigma _ { z x } = \sigma _ { y z } = 0$ , the linear system is reduced to an invertible $3 \times 3$ matrix, 

$$
\left[ \begin{array}{l} \triangle \hat {T} _ {\mathrm {i}, 1} \\ \triangle \hat {T} _ {\mathrm {i}, 2} \\ \triangle \hat {T} _ {\mathrm {i}, 3} \end{array} \right] = \left[ \begin{array}{c c c} g _ {x, 1} g _ {x, 1} & g _ {y, 1} g _ {y, 1} & g _ {x, 1} g _ {y, 1} \\ g _ {x, 2} g _ {x, 2} & g _ {y, 2} g _ {y, 2} & g _ {x, 2} g _ {y, 2} \\ g _ {x, 3} g _ {x, 3} & g _ {y, 3} g _ {y, 3} & g _ {x, 3} g _ {y, 3} \end{array} \right] _ {\hat {M} _ {\mathrm {L O S}} ^ {\mathrm {2 D}}} \left[ \begin{array}{l} \sigma_ {x x} \\ \sigma_ {y y} \\ 2 \sigma_ {x y} \end{array} \right], \tag {3.47}
$$

with the 2-D determinant (det) of $\hat { M } _ { \mathrm { L O S } } ^ { \mathrm { 2 D } }$ 

$$
\det  _ {2 D} = f \times \sin \triangle \phi_ {1 2} \sin \triangle \phi_ {2 3} \sin \triangle \phi_ {3 1}, \tag {3.48}
$$

where $f = - \sin ^ { 2 } \theta _ { 1 } \sin ^ { 2 } \theta _ { 2 } \sin ^ { 2 } \theta _ { 3 }$ , $\triangle \phi _ { 1 2 } = \phi _ { 1 } - \phi _ { 2 }$ , $\triangle \phi _ { 2 3 } = \phi _ { 2 } - \phi _ { 3 }$ , and $\triangle \phi _ { 3 1 } =$ $\phi _ { 3 } - \phi _ { 1 }$ . The determinant $\mathrm { d e t _ { 2 D } }$ is nonzero for any two pairs of nonparallel LOS’s, i.e., $\phi _ { i } - \phi _ { j } \neq 0 , \pi$ . Similarly, in 3-D six LOS’s cannot be chosen to be parallel with each other to avoid forming a singular LOS matrix and the LOS matrix $\hat { M } _ { \mathrm { L O S } } ^ { \mathrm { 3 D } }$ 

should have a small conditional number to minimize error propagations. 

The matrix representation can be utilized to approximate the true minimum of neutron-inferred ion temperature through six ion-temperature measurements. The terminology of Brysk and matrix model will be used interchangeably. Define a column vector $\vec { T _ { \mathrm { i } } }$ to represent six neutron-inferred ion-temperature measurements, a column vector $\vec { T } _ { \mathrm { t h } }$ to represent the thermal ion temperature, and a $6 \times 6$ LOS matrix $\hat { M } _ { 0 } = \hat { M } _ { \mathrm { L O S } } ^ { \mathrm { 3 D } }$ to store six LOS geometrical factors, 

$$
\vec {T} _ {\mathrm {i}} = \sum_ {j = 1} ^ {6} \hat {T} _ {\mathrm {i}, \mathrm {j}} ^ {\text {i n f e r r e d}} \hat {e} _ {j}, \tag {3.49}
$$

$$
{\vec {T} _ {\mathrm {t h}}} = {\hat {T} _ {\mathrm {i}} ^ {\mathrm {t h e r m a l}} \sum_ {j = 1} ^ {6} \hat {e} _ {j},}
$$

where $\{ \hat { e } _ { j = 1 , \ldots , 6 } \}$ are six orthonormal base vectors. Let $\hat { T } _ { \mathrm { p } }$ be a prediction of neutron-inferred ion temperature at an arbitrary LOS $( \theta _ { \mathrm { p } } , \phi _ { \mathrm { p } } )$ , which is uniquely specified by the new LOS matrix $\hat { M } _ { \mathrm { p } }$ and its matrix elements $g _ { i } g _ { j } ( \theta _ { \mathrm { p } } , \phi _ { \mathrm { p } } )$ , 

$$
\begin{array}{r c l}\vec {T} _ {\mathrm {i}}&=&\vec {T} _ {\mathrm {t h}} + \hat {M} _ {0} \cdot \vec {\sigma},\\\vec {\rightarrow}&&\hat {\rightarrow}\end{array}\tag {3.50}
$$

$$
{\vec {T} _ {\mathrm {p}}} = {\vec {T} _ {\mathrm {t h}} + \hat {M} _ {\mathrm {p}} \cdot \vec {\sigma}.}
$$

The solution for the state vector $\vec { \sigma }$ is first obtained by multiplying the inverse of the LOS matrix with the vector of six ion-temperature measurements, 

$$
\vec {\sigma} = \hat {M} _ {0} ^ {- 1} \cdot \left(\vec {T} _ {\mathrm {i}} - \vec {T} _ {\mathrm {t h}}\right). \tag {3.51}
$$

The prediction of neutron-inferred ion temperature at the new LOS is, 

$$
\vec {T} _ {\mathrm {p}} = (\hat {\mathrm {I}} - \hat {M} _ {\mathrm {p}} \cdot \hat {M} _ {0} ^ {- 1}) \cdot \vec {T} _ {\mathrm {t h}} + \hat {M} _ {\mathrm {p}} \cdot \hat {M} _ {0} ^ {- 1} \cdot \vec {T} _ {\mathrm {i}}. \tag {3.52}
$$

The first term on the right-hand side of Eq. (3.52) is denoted by the departure 

matrix $\hat { \delta } _ { \mathrm { p , 0 } }$ , which characterizes the departure of the new LOS from the old one. 

$$
\hat {\delta} _ {\mathrm {p}, 0} = \hat {\mathrm {I}} - \hat {M} _ {\mathrm {p}} \cdot \hat {M} _ {0} ^ {- 1}. \tag {3.53}
$$

The exact relation predicting neutron-inferred ion temperatures at new LOS is, 

$$
\vec {T} _ {\mathrm {p}} = \hat {\delta} _ {\mathrm {p}, 0} \cdot \vec {T} _ {\mathrm {t h}} + \hat {M} _ {\mathrm {p}} \cdot \hat {M} _ {0} ^ {- 1} \cdot \vec {T} _ {\mathrm {i}}. \tag {3.54}
$$

When predictions of neutron-inferred ion temperatures are near the vicinity of old LOS, the term $\hat { \delta } _ { \mathrm { p } , 0 } \cdot \vec { T } _ { \mathrm { t h } }$ approaches zero. Therefore, an approximated solution of neutron-inferred ion temperatures in the neighborhood of old LOS is, 

$$
\vec {T} _ {\mathrm {p}} \simeq \hat {M} _ {\mathrm {p}} \cdot \hat {M} _ {0} ^ {- 1} \cdot \vec {T} _ {\mathrm {i}}, \tag {3.55}
$$

and the minimum neutron-inferred ion temperature in the full map at all angles $\theta _ { \mathrm { p } }$ and $\phi _ { \mathrm { p } }$ is, 

$$
\hat {T} _ {\mathrm {i , m i n}} ^ {\mathrm {a p p r o , i n f}} = \operatorname {M i n} [ \hat {M} _ {\mathrm {p}} (\theta_ {\mathrm {p}}, \phi_ {\mathrm {p}}) \cdot \hat {M} _ {0} ^ {- 1} \cdot \vec {T} _ {\mathrm {i}} ] _ {4 \pi}. \tag {3.56}
$$

To validate the vanishing contribution of the term $\hat { \delta } _ { \mathrm { p } , 0 } \cdot \vec { T } _ { \mathrm { t h } }$ , six LOS’s at the same neutron time-of-flight (nTOF) locations in OMEGA are chosen, indicated by six red dots on the sky map in Fig. (3.15). Positions of six LOS angles complete all matrix elements for $\hat { M } _ { \mathrm { L O S } } ^ { \mathrm { 3 D } }$ or its inverse $\hat { M } _ { 0 } ^ { - 1 }$ . Using OMEGA six LOS configuration, the following vector is computed, 

$$
\hat {M} _ {\mathrm {p}} (\theta_ {\mathrm {p}}, \phi_ {\mathrm {p}}) \cdot \hat {M} _ {0} ^ {- 1} \cdot \vec {T} _ {\mathrm {t h}} = \sum_ {j = 1} ^ {6} a _ {j} \hat {T} _ {\mathrm {i}} ^ {\mathrm {t h e r m a l}} \hat {e} _ {j}, \tag {3.57}
$$

$$
\hat {\delta} _ {p, 0} ^ {3 - \mathrm {D}} \equiv \left(\hat {I} - \hat {M} _ {p} \hat {M} _ {0} ^ {- 1}\right) \vec {T} _ {\mathrm {i}, 3 \mathrm {k e V}} ^ {\mathrm {t h}} \left(\times 1 0 ^ {- 1 5} \mathrm {k e V}\right)
$$

$$
- 7. 5 - 5. 0 - 2. 5 0. 0 2. 5 5. 0
$$

$$
\max  \left[ \hat {\delta} _ {p, 0} ^ {3 - D} \right] _ {4 \pi} = 5. 3 \times 1 0 ^ {- 1 5} \mathrm {k e V}, \min  \left[ \hat {\delta} _ {p, 0} ^ {3 - D} \right] 4 \pi = - 9. 3 \times 1 0 ^ {- 1 5} \mathrm {k e V}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/33587e2ad41c0d0fcb0de7bf77851f3e01e6f1f128dd5c2fffb6598e8c2db930.jpg)



Figure 3.15: A numerical test to show the term $\hat { \delta } _ { \mathrm { p } , 0 } \cdot \vec { T } _ { \mathrm { t h } }$ has vanishing contribution in 4π solid angles by taking the thermal ion temperature as 3 keV. Six LOS’s are indicated by red dots. The coordinates $( \theta , \phi )$ in unit of radians are (1.0700, 0.8314)1, (1.0847, 3.5888) $^ 2$ , (2.0345, 2.8274)3, (0.6705, 4.3563)4, (1.5334, 2.8142)5, and (1.4831, 5.4412)6


with vector elements, 

$$
{ a _ { j } } { = } { \cos ^ { 2 } \theta _ { \mathrm { p } } + \alpha _ { j } \sin 2 \theta _ { \mathrm { p } } \sin \phi _ { \mathrm { p } } + \sin ^ { 2 } \theta _ { \mathrm { p } } ( \cos ^ { 2 } \phi _ { \mathrm { p } } + \beta _ { j } \sin 2 \phi _ { \mathrm { p } } + \sin ^ { 2 } \phi _ { \mathrm { p } } ) . ~ ( 3 . 5 8 ) }
$$

By neglecting the two small constants $\alpha _ { j } = 8 . 9 \times 1 0 ^ { - 1 6 }$ and $\beta _ { j } = - 2 . 2 \times 1 0 ^ { - 1 6 }$ in Eq. (3.58), we have $a _ { j } \simeq 1$ in Eq. (3.57), resulting in a well-approximated zero vector $\hat { \delta } _ { \mathrm { p } , 0 } \cdot \vec { T } _ { \mathrm { t h } } = ( \hat { \mathrm { I } } - \hat { M } _ { \mathrm { p } } \cdot \hat { M } _ { 0 } ^ { - 1 } ) \cdot \vec { T } _ { \mathrm { t h } } \simeq \vec { 0 }$ . Figure (3.15) shows a numerical test to estimate $\widehat { \delta } _ { \mathrm { p } , 0 } ( \theta _ { \mathrm { p } } , \phi _ { \mathrm { p } } ) \cdot \vec { T } _ { \mathrm { t h } }$ at all angles $\theta _ { \mathrm { p } }$ and $\phi _ { \mathrm { p } }$ by direct computing the departure matrix $\hat { \delta } _ { \mathrm { p , 0 } }$ using Eq. (3.53). Magnitudes of the vector components $\hat { \delta } _ { \mathrm { p } , 0 } ( \theta _ { \mathrm { p } } , \phi _ { \mathrm { p } } ) \cdot \vec { T } _ { \mathrm { t h } } \sim 1 0 ^ { - 1 5 } ~ \mathrm { k e V }$ are observed negligibly small. 

Figure (3.16)-(a) shows the reconstruction of $4 \pi$ approximated prediction of neutron-inferred ion temperatures using Eq. (3.55) based on six neutron-inferred ion temperatures from six LOS’s for mode $\ell = 1$ perturbation obtained from an IRIS3D simulation. The neutron-averaged thermal ion temperature $T _ { \mathrm { i } } ^ { \mathrm { t h e r m a l } } =$ 3.551 keV for mode $\ell = 1$ is shown to be in agreement with the prediction of 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/7873414adb8a00216aa3ce796dec887e0d8c341dbf604e8515c918eb38c9ba8b.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/d25ae5cfb5f39ad7b14cc4e3cbdb0436701c9c9314ae66a1c0e670bbae778e2e.jpg)



TC14362J1



Figure 3.16: (a) Prediction of neutron-inferred ion temperatures $T _ { \mathrm { p } }$ in the full map by Eq. (3.55) using six neutron-inferred ion temperatures, shown by the red dots, at six nTOF locations in OMEGA. The neutron-averaged thermal ion ion temperature temperature $\langle T _ { \mathrm { i } } \rangle _ { \mathrm { b } } \equiv T _ { \mathrm { i } } ^ { \mathrm { t h e r m a l } }$ $T _ { \mathrm { i , p r e d } } ^ { \mathrm { t h e r m a l } } \equiv T _ { \mathrm { i , m i n } } ^ { \mathrm { a p p r o . m f } }$ is 3.551 keV, and the minimum of neutron-inferred is 3.525 keV. (b) Comparison of neutroninferred ion temperatures simulated by IRIS3D and the neutron energy spectrum model using 16 LOS’s with that of simulated by the Brysk ion temperatures by Eq. (3.55) using six LOS’s.


minimum neutron-inferred ion temperature $T _ { \mathrm { i , m i n } } ^ { \mathrm { a p p r o . i n f } } ~ = ~ 3 . 5 2 5 ~ \mathrm { k e V }$ obtained by Eq. (3.56). Figure (3.16)-(b) shows that the prediction of neutron-inferred ion measurements using Eq. (3.54) from north to south poles for mode $\ell = 1$ at the fixed angle $\phi = 0$ , agrees with the result of IRIS3D and the neutron energy spectrum model, both using 16 LOS’s uniformly distributed over azimuthal angles $\theta \in [ 0 , \pi ]$ . 

Equation (3.45) can be used to explain the symmetric variation of neutroninferred ion temperatures for mode $\ell = 1$ from north to south poles. DEC3D single-mode database shows that magnitudes of single-mode covariances are negligibly small. 3-D single-mode neutron-inferred ion-temperature variations are therefore governed by three directional variance in Eq. (3.45), 

$$
\triangle \hat {T} _ {\mathrm {i}} ^ {\mathrm {3 D}} (\theta , \phi) = \sin^ {2} \theta \left(\sigma_ {x x} \cos^ {2} \phi + \sigma_ {y y} \sin^ {2} \phi\right) + \sigma_ {z z} \cos^ {2} \theta . \tag {3.59}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/4f884a943eab589b4420ff9b0977d59009d2a3df34678d59f8d770e42c1fdf88.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/31c9948e0ebaeb825c85865a3cd5b6ad93a1b080d01e8720b58168eaf1257f13.jpg)



Figure 3.17: DEC3D single-mode simulations with $\delta v / v _ { 0 } = 0 . 0 7$ initial velocity perturbations. (a) Comparison of the directional and total velocity variances. Mode $\ell = 2$ exhibits the largest total velocity variance. (b) Comparison of the minimum inferred ion-temperature formula in Eq. (3.56) shown by open and solid red circles with the improved approximation formula in Eq. (3.66) shown by open and solid blue squares with isotropic variance separation. Significant improved prediction of thermal ion temperatures is observed.


LOS dependence enters through the square of three geometrical factors , while the magnitudes of three directional varianceof the hot-spot fluid velocity distribution determine the ion-temperature measurement variation at different angles. For 2-D rotational symmetric modes defined by $m = 0$ , the property of directional variance $\sigma _ { x x } = \sigma _ { y y }$ leads to a simplified expression for 2-D single-mode neutron-inferred ion-temperature measurement variations, 

$$
\triangle \hat {T} _ {\mathrm {i}} ^ {\mathrm {2 D}} (\theta , \phi) = (\sigma_ {x x} + \sigma_ {y y}) \sin^ {2} \theta + \sigma_ {z z} \cos^ {2} \theta . \tag {3.60}
$$

Figure (3.17)-(a) investigates the magnitude of directional variance in the single-mode spectrum. For mode $\ell = 1$ , the property of $\sigma _ { x x } + \sigma _ { y y } \ll \sigma _ { z z }$ indicates that the last term $\sigma _ { z z } \cos ^ { 2 } \theta$ in Eq. (3.60) dominates. Therefore, the neutron-inferred ion-temperature measurement variation for mode $\ell = 1$ is, 

$$
\triangle \hat {T} _ {\mathrm {i}} ^ {\ell = 1} \simeq \sigma_ {z z} \cos^ {2} \theta , \tag {3.61}
$$

which explains the sinusoidal curve shown in Fig. (3.16)-(b), with the maximum 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/241decc8c5a75755c1c11ec4ecb4f8a8d614a3204f8a12f4ad095904a33b4f9e.jpg)



Figure 3.18: The full map of ion-temperature measurement variation for low mode $\ell = 4 , m = 2$ by Eq. (3.45) through direct computation of variance and covariance from $D E C 3 D$ hydrodynamic data. Two dim regions are resulted from unequal directional variance $\sigma _ { x x } \neq \sigma _ { y y }$ in Eq. (3.59).


ion temperatures observed at north and south poles and the minimum ion temperature observed on the equator. Although magnitudes of $\sigma _ { z z }$ for modes $\ell = 1$ , 3, and 4 are close to each other, the flow structure for mode $\ell = 1$ has the least non-translational hot-spot fluid motion on the $x - y$ plane, resulting in the largest anisotropic velocity variance $\sigma _ { z z } - \sigma _ { x x }$ or $\sigma _ { z z } - \sigma _ { y y }$ in Fig. (3.17)-(a). Mode $\ell = 2$ is shown to exhibit the largest total directional variance $\sigma _ { x } ^ { 2 } + \sigma _ { y } ^ { 2 }$ in Fig. (3.17)-(a), because of the large radial flow structure within the donut-shape warm bubble. 2-D $m = 0$ modes are shown to have larger directional variance in the $z$ direction than in $x$ or $y$ directions. Because the central spike along the $z$ axis for 2-D $m = 0$ modes behaves as a 3-D spike due to the azimuthal rotational symmetry, the 3-D central spike grows faster than 2-D ring-structure RT spikes, resulting in more intense hot-spot fluid velocity disturbance in the $z$ direction. Figure (3.17)-(b) compares the minimum neutron-inferred ion temperature given by Eq. (3.56) for different single modes, in which six neutron-inferred ion temperatures are simulated by IRIS3D. The property of large isotropic velocity variance for mode $\ell = 2$ is shown to cuase more than 10% higher minimum neutron-inferred ion temperature than the neutron-averaged thermal ion temperature. 

Figure (3.18) shows that mode $\ell = 4 , m = 2$ has sinusoidal ion-temperature measurement variations in the azimuthal direction $\phi$ , which is caused by unequal directional variance $\sigma _ { x x } \neq \sigma _ { y y }$ . The period of 2 variation is due to the presence of $\sin ^ { 2 } \phi$ and $\cos ^ { 2 } \phi$ terms in Eq. (3.59). As a result, 3-D effects of azimuthal iontemperature measurement variations are caused by different values of directional variance in $x$ , $y$ and $z$ directions. Equation (3.34) implies that the variation among different directional variance is the source of neutron-inferred ion temperature measurement asymmetry. Therefore, the directional variance in Eq. (3.59) can be written into a sum of a constant part defined as a directional-independent isotropic velocity variance $\sigma _ { \mathrm { i s o } }$ and a fluctuation part of directional variance $\triangle \sigma _ { i i }$ , 

$$
\sigma_ {i i} = \sigma_ {\mathrm {i s o}} + \triangle \sigma_ {i i}. \tag {3.62}
$$

Effects of isotropic velocity variance on minimum inferred ion temperatures can be explained by substituting Eq. (3.62) into Eq. (3.45). By using the vector property $\begin{array} { r } { \sum _ { i = 1 } ^ { 3 } g _ { i } g _ { i } = \hat { d } \cdot \hat { d } = 1 } \end{array}$ , the isotropic velocity variance $\sigma _ { \mathrm { i s o } }$ is shown to exist uniformly in three different orthogonal directions, 

$$
\sum_ {i = 1} ^ {3} \left(\sigma_ {\text {i s o}} + \triangle \sigma_ {i i}\right) g _ {i} g _ {i} = \sigma_ {\text {i s o}} + \sum_ {i = 1} ^ {3} \triangle \sigma_ {i i} g _ {i} g _ {i}. \tag {3.63}
$$

The fluctuation part $\triangle \sigma _ { i i }$ of directional variance is shown to be the origin of neutron-inferred ion-temperature measurement variations, whereas the isotropic velocity variance $\sigma _ { \mathrm { i s o } }$ is shown to cause minimum neutron-inferred ion temperatures above thermal ion temperatures. Equation (3.45) can be written as, 

$$
\Delta \hat {T} _ {\mathrm {i}} (\theta , \phi) = \sigma_ {\mathrm {i s o}} + \sum_ {i = 1} ^ {3} \Delta \sigma_ {i i} g _ {i} g _ {i} + \sum_ {i \neq j} \sigma_ {i j} g _ {i} g _ {j}. \tag {3.64}
$$

The isotropic velocity variance $\sigma _ { \mathrm { i s o } }$ should be subtracted from all directional variance $\sigma _ { i i }$ in order to shift minimum neutron-inferred ion temperatures in Eq. (3.55) 

down to the level of thermal ion temperatures. By neglecting covariance terms in Eq. (3.64), the functional form of single-mode isotropic velocity variance can be taken as the minimum value among three directional variance $\sigma _ { x x }$ , $\sigma _ { y y }$ and $\sigma _ { z z }$ in order to remove the uniform background of directional variance, 

$$
\sigma_ {\mathrm {i s o}} = \operatorname {M i n} [ \sigma_ {x x}, \sigma_ {y y}, \sigma_ {z z} ]. \tag {3.65}
$$

An improved approximation for thermal ion temperatures is obtained by removing the isotropic velocity variance using Eq. (3.65) from the minimum inferred ion temperature in Eq. (3.55): 

$$
\hat {T} _ {\mathrm {i , a p p r o}} ^ {\mathrm {t h e r m a l}} \simeq \hat {T} _ {\mathrm {i , m i n}} ^ {\mathrm {a p p r o . i n f}} - \sigma_ {\mathrm {i s o}}. \tag {3.66}
$$

The improved approximation for thermal ion temperatures using Eq. (3.66) is studied by blue squares in Fig. (3.17)-(b). After removing the isotropic velocity variance, minimum neutron-inferred ion temperatures indicated by red circles are all down to the level of neutron-averaged thermal ion temperatures for all different single modes, within 3% uncertainty. Low mode $\ell = 1$ is shown to exhibit the least isotropic velocity variance, and its minimum neutron-inferred ion temperature is shown close to the neutron-averaged thermal ion temperature within 1% uncertainty, even without removing the isotropic velocity variance. 

However, solutions of isotropic velocity variance given by Eq. (3.65) are obtained from calculating directional variance using simulation data. $\sigma _ { \mathrm { i s o } }$ cannot be determined in experiments, caused by the fact that both isotropic velocity variance and thermal ion temperatures are directionally independent. As a result, the minimum inferred ion temperature includes a nonseparable contribution from isotropic hot-spot fluid motions. 

The method of vector decomposition can be applied to expand each high order 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/59d655fa7388c2974bf131f413c48e4b97b747de655bff78ff5a503585af4876.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/ce0a6fa0914d7b408fe44312cc50a811d47bfb2626b6f3576286442a04ce0900.jpg)



Figure 3.19: Investigation of 3-D effects of covariance terms for the multimode perturbation $\ell = 1 0 , m = 5$ , and $\ell = 1$ by post-processing $D E C 3 D$ hydrodynamic data to obtain directional variance and covariance: (a) use Eq. (3.44) to obtain $\begin{array} { r } { \hat { T } _ { i } ^ { \mathrm { { i n f e r r e d } } } = \hat { T } _ { i } ^ { \mathrm { { t h e r m a l } } } + \sum _ { i = 1 } ^ { 3 } \sigma _ { i i } g _ { i } g _ { i } + \sum _ { i \neq j } \sigma _ { i j } g _ { i } g _ { j } } \end{array}$ ; (b) use the same formula but neglects directional variance terms to obtain $\begin{array} { r } { \hat { T } _ { i } ^ { \mathrm { i n f e r r e d } } = \hat { T } _ { i } ^ { \mathrm { t h e r m a l } } + \sum _ { i \neq j } \sigma _ { i j } g _ { i } g _ { j } } \end{array}$ .


correction term in Munro’s relativistic treatment [Mun16], at the cost to require more LOS ion-temperature measurements. For 14 MeV neutrons, the relativistic correction is approximately to be < 5% [Mun16, MFH $^ +$ 17]. 

In applications, the true minimum inferred ion temperature provides more realistic evaluations of the true thermal temperature. The reconstruction of neutroninferred ion temperatures through six LOS ion-temperature measurements can also be used to infer orientations of neutron-averaged hot-spot flow velocities. Since the latter is well correlated with neutron-inferred areal density variations in the presence of low mode $\ell = 1 - 2$ , relations between variations in ion temperature and areal density measurements can be used to understand signature of low mode asymmetries in experiments. 

# 3.4.3 Multi-mode velocity variance

Figure (3.19)-(a) investigates the effect of rotational asymmetry on neutron-inferred ion-temperature measurement asymmetry for a two-mode simulation that has a 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/9a01dcabe825c4e6e4cc0747cb4fa5e5b8b0c2ba90753a3b07db9cc70c42e6a9.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/c153d214d1279aa2277470965130923473262003e1a6b0d26591204aa212dd83.jpg)



Figure 3.20: Plot (a) compares $R _ { T }$ in $D E C 3 D$ two-mode simulations (1) with a dominant low mode $\ell = 1$ in the blue curve and (2) with a dominant high mode $\ell = 1 0 , m = 5$ in the red curve. Plot (b) compares the corresponding yield-overclean.


dominant low mode $\ell = 1$ and a small perturbation of high mode $\ell = 1 0 , m = 5$ . The azimuthal rotational symmetry is broke slightly due to the presence of odd- $m$ mode, leading to a significant 3-D neutron-inferred ion-temperature measurement variation on the sky map observed in Fig. (3.19)-(b). 3-D effects of non-zero covariance terms denoted by $\Delta \hat { T } _ { \mathrm { i } } ^ { \mathrm { c o v } }$ in Eq. (3.41) are considered separately, 

$$
\triangle \hat {T} _ {\mathrm {i}} ^ {\mathrm {c o v}} (\theta , \phi) = \sigma_ {x y} \sin^ {2} \theta \sin 2 \phi + \sigma_ {y z} \sin 2 \theta \sin \phi + \sigma_ {z x} \sin 2 \theta \cos \phi . (3. 6 7)
$$

Two bright spots observed in Fig. (3.19)-(b) are caused by the period of 2 in polar and azimuthal angles $\theta$ and $\phi$ . Since covariance terms can be negative, the $4 \pi$ minimum of three covariance terms $\triangle \hat { T } _ { \mathrm { i } } ^ { \mathrm { c o v } }$ should be included inside the minimum bracket in Eq. (3.65) to obtain the value of isotropic velocity variance for multimode perturbations. 

Figure (3.20) compares the ion temperature ratio $R _ { T }$ and yield-over-clean (YOC = 2-D or 3-D fusion yield/ 1-D clean yield) in $D E C 3 D$ two-mode simulations: (1) a dominant low mode $\ell = 1$ with an increasing initial perturbation of a high mode $\ell = 1 0 , m = 5$ ; (2) a dominant high mode $\ell = 1 0 , m = 5$ with an increasing initial perturbation of a low mode $\ell = 1$ . The steeper rise of neutroninferred ion-temperature variation is observed for the dominant high mode $\ell = 1 0$ 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/deeae238b95e2f1643ca70fdf57858c0d86436d621649afafe90182945ade705.jpg)



Figure 3.21: (a) The neutron-inferred ion-temperature sky map simulated by IRIS3D at high resolution by post-processing a strongly distorted DEC3D multimode simulation with an initial perturbation spectrum $\begin{array} { r l } { \Delta v / v _ { 0 } \sum _ { \ell = 1 } ^ { 1 2 } Y _ { \ell } ^ { m = \ell / 2 } } \end{array}$ and YOC = 0.36. (b) $D E C 3 D$ mass density and velocity field profiles on the $x - y$ plane. A developed jet structure is observed in the $x$ direction. The shape of the distorted hot spot is indicated by $T _ { \mathrm { e } } = 0 . 5$ -keV contour surface.


simulations, because of a rapid increasing fluid velocity disturbance in the $z$ direction from $\sigma _ { z z } ^ { \mathrm { k e V } } = 0 . 4 4$ to $\sigma _ { z z } ^ { \mathrm { k e V } } = 1 . 8 2$ , whereas the fluid velocity disturbance on the $x - y$ plane was observed to increase slightly from $\sigma _ { x x } ^ { \mathrm { k e V } } = 0 . 2 9 4$ to $\sigma _ { x x } ^ { \mathrm { k e V } } = 0 . 3 2 5$ . The overall effect of superposition of a mode $\ell = 1$ with a dominant mode $\ell = 1 0$ significantly increases the flow asymmetry in the $z$ direction. 

Figure (3.21)-(a) shows an neutron-inferred ion-temperature sky map simulated by IRIS3D for a strongly perturbed 3-D multi-mode simulation with the a YOC $= 0 . 3 6$ . The initial perturbation spectrum $\begin{array} { r } { \bigtriangleup v / v _ { 0 } \sum _ { \ell = 1 } ^ { 1 2 } Y _ { \ell } ^ { m = \ell / 2 } } \end{array}$ Y m=`/2 consists of modes $\ell = 1 - 1 2$ , $m = 0$ for odd- $\ell$ modes and $m = \ell / 2$ for even- $\ell$ modes. Each mode has the same initial velocity perturbation $\triangle v / v _ { 0 } = 0 . 1 4$ . Directional variance are $\sigma _ { x x } ^ { \mathrm { k e V } } = 1 . 2$ , $\sigma _ { y y } ^ { \mathrm { k e V } } = 0 . 7 4$ , and $\sigma _ { z z } ^ { \mathrm { k e V } } = 0 . 7 6$ , indicating the existence of a dominant jet structure as shown by the velocity field profile in Fig. (3.21)-(b). The motion of the jet on the $x - y$ plane can be read from the neutron-averaged hot-spot fluid velocities $\left. v _ { x } \right. = - 7 . 6 ~ \mathrm { k m / s }$ and $\langle v _ { y } \rangle = - 5 . 4 ~ \mathrm { k m / s }$ , roughly indicating the orientation of the jet from $\theta = \pi / 2 , \phi _ { 1 } = 0 . 9 5 )$ ) to $\theta = \pi / 2 , \phi _ { 2 } = 4 . 1$ ), which accounts for two asymmetric bright spots observed in Fig. (3.21)-(a). 

Figure (3.22)-(a) shows the reconstruction of a 3-D neutron-inferred ion-temperature sky map by Eq. (3.44), and is shown in a good agreement with IRIS3D simulation 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/8c1526961a301856afd8632b170d04fcc0d00eb50beca97c1b2d434b140010d9.jpg)



(rad)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/3300bb325155ed6caadc065e077da0a0858934455555d25cd82af459da8bd02a.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/e25733d746b9fccea9d72c9ac70893de0791695d0735d919d1e0b002f2c2f52f.jpg)



$\phi$ (rad)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/d34ef1b557f2bee9efd8b02ee9d31e5499401642e44c2fde841ec42138d99a33.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/f893704af53027a5c395062a079df05a1140b2e2e91b8a3255b4493010e809e7.jpg)



$\phi$ (rad)



Figure 3.22: Reconstruction of 3-D neutron-inferred ion-temperature profiles by Brysk ion temperature model for the strongly perturbed multi-mode simulation shown in Fig. (3.21) (a) using Eq. (3.44) by including three directional variance terms denoted by $\mathrm { v a r } = \sigma _ { x x } g _ { x } g _ { x } + \sigma _ { y y } g _ { y } g _ { y } + \sigma _ { z z } g _ { z } g _ { z }$ and three covariance terms denoted by $\mathrm { c o v } = 2 \sigma _ { x y } g _ { x } g _ { y } + 2 \sigma _ { y z } g _ { y } g _ { z } + 2 \sigma _ { z x } g _ { z } g _ { x }$ , (b) including only the three directional variance terms, and (c) including only the three covariance terms. The neutron-averaged thermal ion temperature is $\langle T _ { \mathrm { i } } \rangle _ { \mathrm { b } } ^ { I R I S 3 D } = 2 . 7 \mathrm { k e V }$ .


in Fig. (3.21)-(a). Magnitudes of covariance terms σkeV $\sigma _ { x y } ^ { \mathrm { k e V } } = 0 . 2 6$ xy , σkeV $\sigma _ { y z } ^ { \mathrm { k e V } } = - 0 . 0 4 8$ , and $\sigma _ { z x } ^ { \mathrm { k e V } } ~ = ~ - 0 . 0 8 9$ in the strongly perturbed multimode simulation are large enough to cause significant 3-D ion-temperature measurement variations. Effects of directional variance and covariance on neutron-inferred ion-temperature measurement variations are investigated separately in Figs. (3.22)-(b) and (3.22)-(c), showing that the formation of two bright asymmetric spot located at $\theta = \pi / 2$ and $\phi = 0$ in Fig. (3.22)-(a) are caused by large covariance terms. 

# 3.5 Hot-spot residual kinetic energy

Effects of Doppler velocity broadening on neutron energy spectra can be used to infer properties of hot-spot nontranslational residual kinetic energies. The sum of three measurements of ion temperature at orthogonal LOS’s along with $x _ { 1 } = x$ , $x _ { 2 } = y$ , and $x _ { 3 } = z$ directions is related to the thermal ion temperature and the 

total velocity variance, in which three covariance terms vanish exactly, 

$$
\frac {1}{3} \sum_ {j = 1} ^ {3} \hat {T} _ {\mathrm {i}, j} ^ {\text {i n f e r r e d}} = \hat {T} _ {\mathrm {i}} ^ {\text {t h e r m a l}} + \frac {1}{3} \left(\sigma_ {\mathrm {i s o}} + \sum_ {i = 1} ^ {3} \triangle \sigma_ {i i}\right). \tag {3.68}
$$

The appearance of isotropic velocity variance causes larger averaged ion temperatures than thermal ion temperatures by an amount of $\sigma _ { \mathrm { i s o } }$ . Since geometrical factors or LOS effects do not enter the sum of ion-temperature measurements at three orthogonal axes, the apparent ion temperature is further increased by an amount of $\triangle \sigma _ { x x } + \triangle \sigma _ { y y } + \triangle \sigma _ { z z }$ due to anisotropic velocity variance. The total velocity variance and the hot spot nontranslational fluid velocity disturbance are closely related. The hot-spot fluid velocities in three orthogonal directions can be decomposed into a translational component of neutron-averaged linear velocity $\left. v _ { i } \right.$ and a nontranslational component of velocity fluctuation denoted by $\triangle v _ { i } ^ { \mathrm { n o n t r a n s } }$ , 

$$
v _ {i} = \langle v _ {i} \rangle + \triangle v _ {i} ^ {\text {n o n t r a n s}}. \tag {3.69}
$$

Substituting Eq. (3.69) into the definition of directional variance $\sigma _ { i i } = \langle ( v _ { i } -$ $\langle v _ { i } \rangle ) ^ { 2 } \rangle$ , the expression of total velocity variance is equivalent to the inner product of total nontranslational hot-spot fluid velocities, 

$$
\sigma_ {\text {t o t a l}} ^ {2} \equiv \sum_ {i = 1} ^ {3} \sigma_ {i i} = \langle \triangle \vec {v} \cdot \triangle \vec {v} \rangle^ {\text {n o n t r a n s}}. \tag {3.70}
$$

By introducing the total nontranslational hot-spot residual kinetic energy (KE), 

$$
\mathrm {K E} _ {\mathrm {H S}} ^ {\text {n o n t r a n s}} = \frac {1}{2} M _ {\mathrm {H S}} ^ {\mathrm {b}} \left\langle \triangle v ^ {2} \right\rangle^ {\text {n o n t r a n s}}, \tag {3.71}
$$

where $M _ { \mathrm { H S } } ^ { \mathrm { b } } = \textstyle { \frac { 1 } { 2 } } n _ { \mathrm { i } } ^ { \mathrm { b } } m _ { \mathrm { D T } } V _ { \mathrm { H S } } ^ { \mathrm { b } }$ is the neutron-averaged hot-spot mass, $n _ { \mathrm { i } } ^ { \mathrm { b } }$ is the neutronaveraged hot-spot total ion number density, $V _ { \mathrm { H S } } ^ { \mathrm { b } }$ is the neutron-averaged hot-spot volume, $m _ { \mathrm { D T } }$ is the total DT ion-pair mass, and the superscript $^ \mathrm { b }$ denotes the 

spatial neutron-averaging. By expressing the total velocity variance in terms of the hot-spot nontranslational residual kinetic energy using Eq. (3.71), Eq. (3.68) can be written into, 

$$
\frac {1}{3} \sum_ {j = 1} ^ {3} \hat {T} _ {\mathrm {i}, j} ^ {\text {i n f e r r e d}} = \hat {T} _ {\mathrm {i}} ^ {\text {t h e r m a l}} \left(1 + \frac {2 \mathrm {K E} _ {\mathrm {H S}} ^ {\text {n o n t r a n s}}}{3 M _ {\mathrm {H S}} ^ {\mathrm {b}} \hat {T} _ {\mathrm {i}} ^ {\text {t h e r m a l}}}\right). \tag {3.72}
$$

The product of neutron-averaged hot-spot mass and the normalized thermal ion temperature can be simplified into $M _ { \mathrm { H S } } ^ { \mathrm { b } } \hat { T } _ { \mathrm { i } } ^ { \mathrm { t h e r m a l } } \simeq \mathrm { I E _ { H S } ^ { b } } / 6$ . The thermodynamic equilibrium between electron and ion temperatures is assumed to give the neutronaveraged hot-spot total pressure $P _ { \mathrm { H S } } ^ { \mathrm { b } } = 2 n _ { \mathrm { i } } ^ { \mathrm { b } } T _ { \mathrm { i } } ^ { \mathrm { t h e r m a l } }$ , which is a common approximation adopted in deceleration phase hydrodynamic models [SB05a, BCBW16] to express the ion-temperature-dependent fusion reactivity as a function of the total gas pressure. $m _ { \mathrm { D T } } \simeq m _ { \mathrm { n } } + m _ { \alpha }$ is assumed by neglecting the small mass deficit in DT fusion reactions, and $\mathrm { I E _ { H S } ^ { b } = \frac { 3 } { 2 } } P _ { H S } ^ { b } V _ { H S } ^ { b }$ is the neutron-averaged hotspot internal energy. Therefore, the average of three ion temperatures measured at orthogonal directions in Eq. (3.72) is related to the ratio of neutron-averaged nontranslational hot-spot residual kinetic energy to the neutron-averaged hot-spot internal energy: 

$$
\frac {1}{3} \sum_ {j = 1} ^ {3} \hat {T} _ {\mathrm {i}, j} ^ {\text {i n f e r r e d}} = \hat {T} _ {\mathrm {i}} ^ {\text {t h e r m a l}} \left(1 + 4 \frac {\mathrm {K E} _ {\mathrm {H S}} ^ {\text {n o n t r a n s}}}{\mathrm {I E} _ {\mathrm {H S}} ^ {\mathrm {b}}}\right). \tag {3.73}
$$

Equation (3.73) is validated in Fig. (3.23)-(a), in which the nontranslational hotspot kinetic energy is replaced by the total hot-spot kinetic energy KEtotalHS because their magnitudes are close: $\mathrm { K E _ { H S } ^ { t o t a l } \simeq K E _ { H S } ^ { n o n t r a n s } }$ . This approximation is only valid for modes $\ell \ \geq \ 2$ with small amount of translational hot-spot residual kinetic energies. Figure (3.23)-(b) investigates the effect of isotropic velocity variance by comparing minimum neutron-inferred ion temperatures obtained by Eq. (3.56) with thermal ion temperatures in the mode spectrum. Mode $\ell = 1$ is shown to 

exhibit the least isotropic flow so that its minimum inferred ion temperature is close to the thermal ion temperature. Mode $\ell = 2$ is shown to exhibit large deviations between minimum inferred temperatures and thermal temperatures because of large isotropic velocity variance. The intimate relation between ion-temperature measurement variation and residual kinetic energy can be shown by considering the conservation of total energy at the time of stagnation between 3-D and 1-D implosions [WBS+18a], 

$$
\mathrm {I E} _ {\mathrm {H S}, 3 \mathrm {D}} ^ {\mathrm {b}} / \mathrm {I E} _ {\mathrm {H S}, 1 \mathrm {D}} ^ {\mathrm {b}} = 1 - \mathrm {R K E} _ {\text {t o t a l}} - \mathrm {R I E} _ {\mathrm {S H}}, \tag {3.74}
$$

where notations for the normalized total residual kinetic energy RKE $\mathrm { \ t o t a l { \Omega } }$ , the normalized neutron-averaged hot-spot residual kinetic energy RKEbHS, and the normalized shell internal energy RIESH are defined by 

$$
\mathrm {R K E} _ {\mathrm {t o t a l}} = \left(\mathrm {K E} _ {\mathrm {t o t a l}} ^ {\mathrm {3 D}} - \mathrm {K E} _ {\mathrm {t o t a l}} ^ {\mathrm {1 D}}\right) / \mathrm {I E} _ {\mathrm {H S , 1 D}} ^ {\mathrm {b}},
$$

$$
\mathrm {R K E} _ {\mathrm {H S}} ^ {\mathrm {b}} = \left(\mathrm {K E} _ {\mathrm {H S}} ^ {\mathrm {b}, 3 \mathrm {D}} - \mathrm {K E} _ {\mathrm {H S}} ^ {\mathrm {b}, 1 \mathrm {D}}\right) / \mathrm {I E} _ {\mathrm {H S}, 1 \mathrm {D}} ^ {\mathrm {b}}, \tag {3.75}
$$

$$
\mathrm {R I E} _ {\mathrm {S H}} = (\mathrm {I E} _ {\mathrm {S H}} ^ {\mathrm {3 D}} - \mathrm {I E} _ {\mathrm {S H}} ^ {\mathrm {1 D}}) / \mathrm {I E} _ {\mathrm {H S}, 1 \mathrm {D}} ^ {\mathrm {b}}.
$$

The ratio of the neutron-averaged hot-spot residual kinetic energy to the neutronaveraged hot-spot internal energy can be rewritten as $\mathrm { K E _ { H S , 3 D } ^ { b } / I E _ { H S , 3 D } ^ { b } = R K E _ { H S } ^ { b } / ( 1 - }$ RKEtotal − RIESH). Except for low mode $\ell = 1$ , the neutron-averaged nontranslational hot-spot residual kinetic energy in Eq. (3.73) can be well approximated by the neutron-averaged total hot-spot residual kinetic energy since modes $\ell \geq 2$ contain less translational hot-spot fluid motions. Therefore, the summation of neutron-inferred ion temperatures at three orthogonal directions in Eq. (3.73) is a function of residual kinetic energies: 

$$
\frac {1}{4} \left(\frac {\sum_ {j = 1} ^ {3} \hat {T} _ {\mathrm {i} , j} ^ {\mathrm {i n f e r r e d}}}{3 \hat {T} _ {\mathrm {i}} ^ {\mathrm {t h e r m a l}}} - 1\right) \simeq \frac {\mathrm {R K E} _ {\mathrm {H S}} ^ {\mathrm {b}}}{1 - \mathrm {R K E} _ {\mathrm {t o t a l}} - \mathrm {R I E} _ {\mathrm {S H}}}. \tag {3.76}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/50e9ee4ec007fb86f2f05b7d1139ad6582d3dd05fa19a944cb1d7ce67e69a6f3.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/6d6f22aee8ccb078fbaf20e77d18d8e9828efab39e0d368433bfcef75e3c7528.jpg)



Figure 3.23: DEC3D 77068 single-mode database. Plot (a) validates Eq. (3.73) by comparing the sum of three neutron-inferred ion temperatures at $\hat { x }$ , $\hat { y }$ , and $\hat { z }$ to the thermal ion temperature and total variance. The red dots are low mode $\ell = 1$ . Plot (b) compares the minimum neutron-inferred ion temperatures defined by Eq. (3.56) with the neutron-averaged thermal ion temperatures for all single modes. The red dots are low mode $\ell = 1$ and the blue dots are low mode $\ell = 2$ . Plot (c) compares the average neutron-inferred ion temperatures over three orthogonal directions with hot-spot and shell residual kinetic energies to validate Eq. (3.76).


Equation (3.76) shows that increasing neutron-inferred ion temperatures are caused by nontranslational hot-spot residual kinetic energies. Only when $\mathrm { K E _ { H S } ^ { n o n t r a n s } = 0 }$ in Eq. (3.72), the average of neutron-inferred ion temperature measured at three orthogonal axes is equal to the thermal ion temperature. Figure (3.23)(c) shows that Eq. (3.76) is valid for all modes $\ell \geq 2$ . Mode $\ell = 1$ is the outlier because it has large translational hot-spot fluid motions that violate the approximation of KEnontransHS ' KEbHS,3D in Eq. (3.76). 

# 3.6 Conclusion

A comprehensive analysis of velocity variance in Brysk ion temperatures [Bry73, Mur14, Mun16] is presented for single-mode and multimode RT instabilities in the deceleration phase of ICF implosions. 3-D effects of hot-spot flow asymmetry on variations in ion-temperature measurement is shown to be uniquely determined by a complete set of six hot-spot flow parameters including three directional variance and three covariance, which are calculated from hot spot fluid velocity distribu-

tions. An approximated solution to minimum neutron-inferred ion temperatures is derived, and the approximated neutron-inferred ion-temperature profiles over full $4 \pi$ solid angles are shown to be well reconstructed by six ion-temperature measurements. The reconstructed neutron-inferred ion-temperature profiles are in agreement with the result of the neutron energy spectrum model [GSB $^ +$ 14] and Monte Carlo neutron transport code IRIS3D. [WRF18] Predictions of minimum neutron-inferred ion temperatures for low mode $\ell = 1$ are shown to close to thermal ion temperatures because mode $\ell = 1$ exhibits the least isotropic velocity variance. Low mode $\ell = 2$ is shown to exhibit a large isotropic velocity variance due to its large radial flow structure within the donut-shape warm bubble on the equatorial plane. A large isotropic velocity variance is shown to cause the minimum neutron-inferred ion temperature well above the thermal ion temperature. An improved approximation of the thermal ion temperature is derived by separating the isotropic velocity variance. The presence of nontranslational hot-spot residual kinetic energy is shown to cause larger averaged neutron-inferred ion temperatures measured at three orthogonal directions than the thermal ion temperature. 

# 4 Impact of Isotropic Flows within the Hot Spot

# 4.1 Motivation

In National Ignition Facility (NIF) experiments [GJKC $^ +$ 16], DD neutron-inferred ion temperatures less than that of DT were inferred. Ratios of DD to DT neutroninferred ion temperatures were observed in between $T _ { \mathrm { D D } } / T _ { \mathrm { D T } } ~ = ~ 0 . 8 - 1$ , and simultaneously these shots exhibited a small ion-temperature measurement variation among different LOS. The “isotropic source ”was hypothesized [GJKC $^ +$ 16] to explain DD and DT neutron-inferred ion temperatures being well above thermal ion temperatures. Although the physical interpretation for the velocity variance was described in Murphy [Mur14] for the limit of fully developed turbulence and by Munro [Mun16, MFH $^ +$ 17] with relativistic corrections, an accurate description for the isotropic source from the structure of velocity variance is missing. The presence of isotropic sources leads to overestimation of the inferred hot-spot pressure [SGE $^ + 1 2$ , KTB+14], which is an important metric to assess implosion performance. 

In Chapter 3, the technique of decomposing the velocity variance into six hotspot flow parameters [WBS $^ +$ 18b], provides a unified framework to understand the complex 3-D flow structure analytically. In this chapter, a comprehensive analysis for the structure of velocity variance is presented, through applying the technique of Chapter 3, with the objective to derive a general expression for the isotropic 

velocity variance for single-mode and multi-mode perturbations, and investigate the impact of isotropic sources on causing the DD/DT ion-temperature ratio to fall below unity. 

# 4.2 Physical meaning of variance and covariance

In Chapter 3, the technique [WBS+18b] that decomposes the velocity variance in the non-relativistic Brysk ion temperatures in Eq. (3.30) is extended to a full relativistic treatment by retaining $\triangle _ { 1 }$ and $\triangle _ { 2 }$ terms in Eq. (3.16). However, magnitudes of relativistic corrections are too small compared with large ion-temperature measurement variations caused by the velocity variance term. The extension of Chapter 3 technique to a relativistic version is proposed as for future work, because the main component of isotropic source is hidden in the non-relativistic Brysk ion temperatures. 

By substituting the fluid velocity $\vec { v } = v _ { i } \hat { e } _ { i }$ measured in the laboratory frame and the LOS unit vector $\hat { d } = g _ { i } \hat { e } _ { i }$ into the velocity variance var $\Big [ \vec { v } \cdot \hat { d } \Big ]$ , the compact representation for non-relativistic Brysk ion temperatures in Eq. (3.30) is, 

$$
\hat {T} _ {\mathrm {i}} ^ {\text {i n f e r r e d}} = \hat {T} _ {\mathrm {i}} ^ {\text {t h e r m a l}} + g _ {i} g _ {j} \sigma_ {i j}, \tag {4.1}
$$

where Einstein summation convention for any two repeating indices is used. The summation indices running over $i \ = \ 1$ , 2 and 3 correspond to Cartesian coordinates $x$ , $y$ and $z$ respectively. $\{ \hat { e } _ { i } \}$ is a set of orthonormal unit vectors that parallel with the $x$ -axis: $\hat { e } _ { 1 } = ( 1 , 0 , 0 )$ , the $y$ -axis: $\hat { e } _ { 1 } = ( 0 , 1 , 0 )$ and the $z$ -axis: $\hat { e } _ { 3 } = ( 0 , 0 , 1 )$ respectively. The Kronecker delta $\delta _ { i j }$ is unity for $i = j$ and zero for $i \neq j$ , which naturally splits the summation term in Eq. (3.44) into three directional-variance “var” and three covariance “cov”. The notation of “var” represents the directional-variance and it should not be confused with that for the 

velocity variance var $\Big [ \vec { v } \cdot \hat { d } \Big ]$ . Define the projector operator $\mathcal { P }$ to select $i = j$ and $i \neq j$ components in the summation. 

$$
\operatorname {v a r} = \mathcal {P} _ {i = j} \left[ g _ {i} g _ {j} \sigma_ {i j} \right] = g _ {1} g _ {1} \sigma_ {1 1} + g _ {2} g _ {2} \sigma_ {2 2} + g _ {3} g _ {3} \sigma_ {3 3}, \tag {4.2}
$$

$$
\operatorname {c o v} = \mathcal {P} _ {i \neq j} \left[ g _ {i} g _ {j} \sigma_ {i j} \right] = 2 g _ {1} g _ {2} \sigma_ {1 2} + 2 g _ {2} g _ {3} \sigma_ {2 3} + 2 g _ {3} g _ {1} \sigma_ {3 1}. \tag {4.3}
$$

In Chapter 3, negligible small single-mode covariance terms were observed in DEC3D single-mode database. This property is not valid in multi-mode perturbations, because the onset of rotational asymmetry triggers a transition from zero single-mode covariance to non-zero multi-mode covariance. 

$$
\sigma_ {1 2} = \left\langle v _ {1} v _ {2} \right\rangle - \left\langle v _ {1} \right\rangle \left\langle v _ {2} \right\rangle , \tag {4.4}
$$

$$
\sigma_ {2 3} = \left\langle v _ {2} v _ {3} \right\rangle - \left\langle v _ {2} \right\rangle \left\langle v _ {3} \right\rangle , \tag {4.5}
$$

$$
\sigma_ {3 1} = \left\langle v _ {3} v _ {1} \right\rangle - \left\langle v _ {3} \right\rangle \left\langle v _ {1} \right\rangle . \tag {4.6}
$$

An accurate physical interpretation for directional-variance and covariance terms is derived as follows in order to build up the physical picture for an isotropic source systematically. The burn-averaged linear $\left. v _ { i } \right.$ and bilinear $\left. v _ { i } v _ { j } \right.$ hot-spot fluid velocities are non-zero only for odd- $m$ single modes. 

$$
\left\langle v _ {1} \right\rangle = \delta_ {m, m _ {\mathrm {o d d}}} \left\langle v _ {1} \right\rangle_ {m _ {\mathrm {o d d}}}, \tag {4.7}
$$

$$
\langle v _ {2} \rangle = \delta_ {m, m _ {\mathrm {o d d}}} \langle v _ {2} \rangle_ {m _ {\mathrm {o d d}}}, \tag {4.8}
$$

$$
\langle v _ {3} \rangle = \delta_ {\ell , \ell_ {\mathrm {o d d}}} \langle v _ {3} \rangle_ {\ell_ {\mathrm {o d d}}}, \tag {4.9}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/81129766e2c262c371da628219f4865c83e8d1993981729705f5aa97a123de78.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/700a0ab9a76c6a0a087867f0ed9e2591c2a8a7ae5ae56246857ede8e7d91ce0b.jpg)



Figure 4.1: Comparison of single-mode burn-averaged linear velocities in (a) and bilinear velocities in (b) at stagnations simulated by $D E C 3 D$ at 7% initial velocity perturbation for OMEGA shot 77068. Mode $\ell = 1$ has large burn-averaged $z$ - velocity due to the jet. Burn-averaged bilinear velocities in (b) are negligible. The combined results in (a) and (b) lead to vanishing single-mode covariance.


and 

$$
\langle v _ {1} v _ {2} \rangle = \delta_ {m, m _ {\mathrm {o d d}}} \langle v _ {1} v _ {2} \rangle_ {m _ {\mathrm {o d d}}}, \tag {4.10}
$$

$$
\langle v _ {2} v _ {3} \rangle = \delta_ {m, m _ {\mathrm {o d d}}} \langle v _ {2} v _ {3} \rangle_ {m _ {\mathrm {o d d}}}, \tag {4.11}
$$

$$
\left\langle v _ {3} v _ {1} \right\rangle = \delta_ {m, m _ {\mathrm {o d d}}} \left\langle v _ {3} v _ {1} \right\rangle_ {m _ {\mathrm {o d d}}}. \tag {4.12}
$$

Figure (4.1) shows the survey of burn-averaged linear and bilinear velocities in DEC3D deceleration-phase single-mode database, with the mode spectrum $\ell =$ $1 - 1 2$ including 2-D $m = 0$ modes and 3-D $m = \ell _ { \mathrm { e v e n } } / 2$ modes. Periodic jumps of burn-averaged linear and bilinear velocities are observed only to occur 3-D $m \neq 0$ modes, which lead to the result of zero single-mode covariance, 

$$
\sigma_ {1 2} ^ {(m _ {\mathrm {e v e n}})} = \sigma_ {2 3} ^ {(m _ {\mathrm {e v e n}})} = \sigma_ {3 1} ^ {(m _ {\mathrm {e v e n}})} = 0. \tag {4.13}
$$

To understand the physical origin of zero single-mode covariance caused by azimuthal symmetry, the fluid velocity $\vec { v }$ is decomposed into a burn-averaged component $\langle \vec { v } \rangle$ that represents the mean fluid velocity within the distorted 3-D hot spot and a fluctuation component $\triangle \vec { v }$ , in a similar manner of Reynolds decomposition 

in fluid mechanics, 

$$
\vec {v} (\vec {x}, t) = \langle \vec {v} (t) \rangle + \triangle \vec {v} (\vec {x}, t). \tag {4.14}
$$

In the spatial integration to calculate the burn-averaged or neutron-averaged bracket, the covariance term $\sigma _ { i j } = \langle ( v _ { i } - \langle v _ { i } \rangle ) ( v _ { j } - \langle v _ { j } \rangle ) \rangle$ is separated into two equal halves, defined by the left integration domain $\mathcal { D } _ { L } : 0 \le \phi < \pi$ and the right integration domain $\mathcal { D } _ { R } : \pi \le \phi < 2 \pi$ , 

$$
\sigma_ {i j} = \langle \triangle v _ {i} (\phi_ {\mathcal {D} _ {L}}) \triangle v _ {j} (\phi_ {\mathcal {D} _ {L}}) \rangle + \langle \triangle v _ {i} (\phi_ {\mathcal {D} _ {R}}) \triangle v _ {j} (\phi_ {\mathcal {D} _ {R}}) \rangle . (4. 1 5)
$$

The physical meaning of azimuthal symmetry in angle $\phi$ implies the conservation of a zero total translational momentum of the whole capsule in the direction orthogonal to the rotation axis. Equivalently, the whole capsule experiences a zero net torque, which determines the rotation of the whole capsule off the rotational axis. The azimuthal symmetry is manifested by requiring the existence of a pair of parallel and anti-parallel flow velocities rotational symmetrically on the $x - y$ plane in the single-mode perturbations, which immediately leads to the zero covariance for even- $m$ modes in Eqs. (4.13) and (4.15). 

$$
\triangle v _ {i} (\theta , \phi_ {\mathcal {D} _ {L}}) = - \triangle v _ {i} (\theta , \phi_ {\mathcal {D} _ {R}}). \tag {4.16}
$$

In the presence of odd- $m$ modes, the whole capsule rotates and translates on the plane orthogonal to the rotation axis, meaning that the loss of azimuthal symmetric configuration to maintain the perfect cancellation of the pair of parallel and anti-parallel flows in the covariance term in Eq. (4.15). However, magnitudes of non-zero covariance for odd- $m$ modes decrease with the azimuthal variation number $m$ . Because for large odd- $m$ modes, flow velocities are highly localized within the cold bubbles, which are even colder with lower burn weights as RT spikes converge to form the secondary piston. As a result, the overall contribution to the 

covariance terms due to the azimuthal symmetric high-velocity vorticity within the cold bubbles decreases as $m$ increases. Therefore, single-mode covariance terms are analytically zero for even- $m$ modes and approach zero in the limit of large odd- $m$ modes, 

$$
\lim  _ {m _ {\mathrm {o d d}} \rightarrow \infty} \sigma_ {i j} ^ {(m _ {\mathrm {o d d}})} = 0. \tag {4.17}
$$

Eq. (4.7) shows that burn-averaged linear velocities $\left. v _ { i } \right.$ are non-zero only for odd- $\ell$ and odd- $m$ modes. Since non-zero $\left. v _ { i } \right.$ can result in translation of the whole capsule in $x$ , $y$ or $z$ directions due to unbalanced momenta, the terminology of non-translational velocity is adopted to describe the term $\triangle \vec { v } = \vec { v } - \langle \vec { v } \rangle$ . The physical meaning of the directional-variance, which is defined as the square of fluctuation velocity with respect to the mean, 

$$
\sigma_ {i i} = \left\langle \left(v _ {i} - \left\langle v _ {i} \right\rangle\right) \cdot \left(v _ {i} - \left\langle v _ {i} \right\rangle\right) \right\rangle = \left\langle \triangle v _ {i} \triangle v _ {i} \right\rangle , \tag {4.18}
$$

is related to the residual kinetic energy of the non-translational hot-spot fluid velocity in the direction $\hat { e } _ { i }$ defined by, 

$$
\mathrm {K E} _ {\mathrm {h s}, i} ^ {\mathrm {n o n t r a n s}} = \frac {1}{2} M _ {\mathrm {h s}} \langle \triangle v _ {i} \triangle v _ {i} \rangle , \tag {4.19}
$$

or equivalently 

$$
\sigma_ {i i} = 2 \mathrm {K E} _ {\mathrm {h s}, i} ^ {\mathrm {n o n t r a n s}} / M _ {\mathrm {h s}}. \tag {4.20}
$$

The physical meaning of the covariance is related to the measurement of the degree of azimuthal asymmetry, 

$$
\sigma_ {i j} = \left\langle \left(v _ {i} - \left\langle v _ {i} \right\rangle\right) \cdot \left(v _ {j} - \left\langle v _ {j} \right\rangle\right) \right\rangle = \left\langle \triangle v _ {i} \triangle v _ {j} \right\rangle . \tag {4.21}
$$

For fully-developed turbulence, the fluctuation velocity $\triangle v _ { i }$ can be treated as a 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/2fbe783b9e429d280f02847c5c575318ab3c2078c66faf101a3ff78e25a6e7bf.jpg)



Figure 4.2: A level diagram to represent the transition of levels of the absolute magnitudes of covariance from single modes to multi-modes, and finally to a fullydeveloped turbulence by random mixing of multi-modes.


random variable with a zero mean $\langle v _ { i } \rangle = 0$ . Using the property of zero correlation between any two pair of $\triangle v _ { i }$ and $\triangle v _ { j }$ , the burn-averaged bracket in Eq. (4.21) can be decomposed into a product of two independent burn-averaged fluctuation velocities $\langle \triangle v _ { i } \triangle v _ { j } \rangle = \langle \triangle v _ { i } \rangle \langle \triangle v _ { j } \rangle = 0$ , resulting in zero covariance in the limit of fully-developed turbulence. Figure (4.2) shows a level diagram to summarize the transition of magnitudes of covariance from highly-structure azimuthal-symmetric single modes to multi-modes, and finally to fully-developed turbulence. Magnitudes of covariance for odd- $m$ modes are only slightly above zero, as explained by the small contribution to burn-averaged velocities due to the fact of low burn weights for azimuthal asymmetric high-velocity vorticity within the cold bubbles. Both single-mode perturbations and fully-turbulent hot spots satisfy the same Brysk ion temperatures in Eq. (4.1) with zero covariance terms. This observation agrees with Murphy work [Mur14], in which no covariance terms are present in analyzing effects of turbulent kinetic energy on inferred ion temperatures. 

# 4.3 Properties of isotropic velocity variance

In this section, the analysis of isotropic source will proceed from single modes to multi modes, following the order from the left to the right in the level diagram in Fig. (4.2). For single-mode neutron-inferred ion-temperatures, three covariance terms are analytically zero in Eq. (4.1) to conserve the rotational symmetry. 

$$
\hat {T} _ {\mathrm {i}, \text {s i n g l e - m o d e}} ^ {\text {i n f e r r e d}} = \hat {T} _ {\mathrm {i}} ^ {\text {t h e r m a l}} + \mathcal {P} _ {i = j} [ g _ {i} g _ {j} \sigma_ {i j} ]. \tag {4.22}
$$

Following the same procedure described in Chapter 3, the directional-variance $\sigma _ { i i } = \sigma _ { \mathrm { i s o } } + \triangle \sigma _ { i i }$ is decomposed into a constant part $\sigma _ { \mathrm { i s o } }$ that is uniform in space and a fluctuation part $\triangle \sigma _ { i i }$ . The exact form for single-mode neutron-inferred ion-temperature measurement is, 

$$
\hat {T} _ {\mathrm {i , s i n g l e - m o d e}} ^ {\mathrm {i n f e r r e d}} = \hat {T} _ {\mathrm {i}} ^ {\mathrm {t h e r m a l}} + \sigma_ {\mathrm {i s o}} + \mathcal {P} _ {i = j} [ g _ {i} g _ {j} \triangle \sigma_ {i j} ]. \tag {4.23}
$$

The vector property $g _ { i } g _ { j } \delta _ { i j } = 1$ is used to derive above equation. Equation (4.23) states that any directional-variance below the minimum of $\sigma _ { x x }$ , $\sigma _ { y y }$ and $\sigma _ { z z }$ does not contribute to ion-temperature measurement variations but the fluctuation parts $\triangle \sigma _ { i j }$ contribute, because $\triangle \sigma _ { i j }$ are multiplied with directional-dependent geometrical factors. This observation leads to the following general expression for single-mode isotropic velocity variance, 

$$
\sigma_ {\mathrm {i s o}} (t) = \operatorname {M i n} \left[ \sigma_ {x x} (t), \sigma_ {y y} (t), \sigma_ {z z} (t) \right]. \tag {4.24}
$$

For $m = 0$ modes, $\sigma _ { x x } ( t ) = \sigma _ { y y } ( t )$ due to the rotational symmetry, the isotropic velocity variance is simplified into, 

$$
\sigma_ {\mathrm {i s o}} ^ {m = 0} (t) = \operatorname {M i n} [ \sigma_ {x x} (t), \sigma_ {z z} (t) ], \tag {4.25}
$$

Using Eq. (4.20), the isotropic velocity variance is related to hot-spot residual kinetic energies. Therefore, the minimum neutron-inferred ion temperature in Eq. (4.23), being well above the thermal ion temperature is caused the presence of non-translational hot-spot residual kinetic energy. 

The time-dependency of isotropic velocity variance in Eq. (4.24) is used to capture the dynamical behavior of directional residual kinetic energy during the neutron production temporal history. For modest even $\ell$ -mode perturbations, RT spikes do not reach the origin at radius $r = 0$ during the whole pre-stagnation and post-stagnation phases. However, for large even $\ell$ -mode perturbations, the fastgrowing RT spikes for low modes can reach the origin during the peak compression, leading to a significant reduction of burn volume in the direction parallel to the pair of center-reaching RT spikes and reduce the local burn weight for any flow nearby the origin. An large mode $\ell = 2$ simulation is observed to have a timevarying $z$ -directional variance. Before the pair of RT spikes reaching the center, the magnitude of $\sigma _ { z z }$ is dominated by the burn-averaged fluctuation velocity $\triangle v _ { z } ^ { 2 }$ along the $z$ -axis. After the pair of RT spikes reaching the center, the burn volume along the $z$ -axis is reduced significantly, resulting in a sharp decrease of $z$ -directional variance $\sigma _ { z z } ^ { \ell = 2 } \ll \sigma _ { x x } ^ { \ell = 2 }$ . The magnitude of expanding radial symmetric flow $\sigma _ { x x }$ or $\sigma _ { y y }$ remains growing within the large donut-shape warm bubble irregardless large or modest perturbations. The general expression of mode $\ell = 2$ isotropic velocity variance is, 

$$
\sigma_ {\mathrm {i s o}} ^ {\ell = 2} (t) = \operatorname {M i n} \left[ \sigma_ {x x} ^ {\text {b u b b l e}} (t), \sigma_ {z z} ^ {\text {s p i k e}} (t) \right]. \tag {4.26}
$$

For a phase-reversed mode $\ell = 2$ , the superscript for bubble and $\mathrm { s p i k e }$ interchanges. The relative importance of residual kinetic energies contributing to isotropic source between converging RT spikes and expanding RT bubbles depends on the time-

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/7eef68f440b2ffc69ffee4f93355f0c958aab21925ff266168866ba4d0f11180.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/be60b3a53b8c039d095f9374659437c7a4d014fcef12b0ddb69daf48fbf4d266.jpg)



Figure 4.3: Comparison of time-integrated burn-averaged directional-variance $\langle \sigma _ { x x } \rangle$ and $\langle \sigma _ { z z } \rangle$ against different levels of initial velocity perturbations. The transition of $\left. \sigma _ { z z } \right. \mathrm { ~ < ~ } \left. \sigma _ { x x } \right.$ is observed in large mode $\ell = 2$ perturbations with $\triangle v / v _ { 0 } > 1 0 \%$ in $D E C 3 D$ simulations.


integrated burn-averaged isotropic velocity variance, 

$$
\left\langle \sigma_ {\mathrm {i s o}} \right\rangle = \int \dot {Y} (t) \sigma_ {\mathrm {i s o}} (t) d t / \int \dot {Y} (t) d t, \tag {4.27}
$$

where $\dot { Y } ( t )$ is the neutron production rate at time $t$ . In an large mode $\ell = 2$ perturbation, the minimum neutron-inferred ion temperature was observed perpendicular to the equatorial plane of the large donut-shape warm bubble, where the maximum neutron-inferred ion temperature is observed. In an modest mode $\ell = 2$ perturbation [WBS+18b], the minimum neutron-inferred ion temperature was observed parallel to the equatorial plane of the large donut-shape warm bubble, and the maximum neutron-inferred ion temperature was observed along the poles parallel to RT spikes. 

Figure (4.3) compares the time-integrated burn-averaged directional-variance $\langle \sigma _ { x x } \rangle$ and $\langle \sigma _ { z z } \rangle$ between modes $\ell = 1$ and $\ell = 2$ against different levels of initial velocity perturbations. The transition of $\sigma _ { z z } ^ { \ell = 2 } < \sigma _ { x x } ^ { \ell = 2 }$ for mode $\ell = 2$ occurs at $1 0 \%$ initial velocity perturbation. The non-translational residual kinetic energy $\mathrm { R K E } _ { \mathrm { e q u a t o r i a l } } ^ { \ell = 2 } = \textstyle { \frac { 1 } { 2 } } M _ { \mathrm { h s } } \sigma _ { x x , \ell = 2 } ^ { \mathrm { b u b b l e } }$ in Eq. (4.20) driven by the expanding radial flow structure, parallel to the equatorial plane, within the large donut-shape warm bubble for mode $\ell = 2$ is shown to grow monotonically with initial velocity perturbations, whereas mode $\ell = 1$ shows negligible small RKE`=1equatorial. At large ini-

tial velocity perturbations $\geq 1 0 \%$ , however, the growing non-translational residual kinetic energy RKE`=2equatorial eventually takes over the non-translational residual kinetic energy $\begin{array} { r } { \mathrm { R K E } _ { \mathrm { p o l e } } ^ { \ell = 2 } = \frac { 1 } { 2 } M _ { \mathrm { h s } } \sigma _ { z z , \ell = 2 } ^ { \mathrm { s p i k e } } } \end{array}$ driven by the pair of RT spikes along the pole. Figure (4.3) shows a sharp transition of $z$ -directional variance $\left. \sigma _ { z z } \right. \mathrm { ~ < ~ } \left. \sigma _ { x x } \right.$ for mode $\ell = 2$ at large initial perturbations. After the pair of RT spikes symmetrically reaching the origin at radius $r = 0$ , the configuration of a pair of parallel and an anti-parallel $z$ -directional flows along the pole is lost and the main contribution to the burn-averaged $z$ -directional variance is caused by $v _ { z }$ from vorticicty further away from the center, where burn weights are much lower. As a result, the burn-averaged $z$ -directional variance is significantly reduced after the transition to large mode $\ell = 2$ perturbations. Define $\begin{array} { r } { \triangle T = T _ { \mathrm { i , m a x } } ^ { \mathrm { i n f e r r e d } } - T _ { \mathrm { i , m i n } } ^ { \mathrm { i n f e r r e d } } } \end{array}$ as the iontemperature measurement asymmetry, the summary of isotropic velocity variance for mode $\ell = 2$ in Eq. (4.26) is, 

$$
\sigma_ {\text {i s o , n o n a g g}} ^ {\ell = 2} = \sigma_ {x x} ^ {\text {b u b b l e}} \rightarrow \triangle T _ {\text {s m a l l}} \text {a n d} \sigma_ {\text {i s o}} ^ {\text {l a r g e}}, \tag {4.28}
$$

$$
\sigma_ {\text {i s o , a g g}} ^ {\ell = 2} = \sigma_ {z z} ^ {\text {s p i k e}} \rightarrow \triangle T _ {\text {l a r g e}} \text {a n d} \sigma_ {\text {i s o}} ^ {\text {s m a l l}}, \tag {4.29}
$$

In Chapter 3, modest mode $\ell = 2$ simulations were shown to contribute the most isotropic source in the single-mode spectrum $\ell = 1 - 1 2$ , with small iontemperature ratios $R _ { T } ~ = ~ T _ { \mathrm { i , m a x } } ^ { \mathrm { i n t e r r e d } } / T _ { \mathrm { i , m i n } } ^ { \mathrm { i n t e r r e d } }$ and small inferred ion temperature variations among different LOS in Fig. (3.9). modest perturbations occur more frequent than large perturbations in high-performance ICF implosion experiments. Unless implosions are ideally 1-D, signature of small ion-temperature measurement variations among different LOS indicate a large non-separable isotropic source in the minimum neutron-inferred ion temperatures. 

Figure (4.4) shows that no transition is observed for the rest of mode spectrum $\ell = 3 - 1 2$ by comparing the ratio $\langle \sigma _ { z z } ^ { \mathrm { D T } } \rangle / \langle \sigma _ { x x } ^ { \mathrm { D T } } \rangle$ with different levels of initial velocity perturbations. The $x$ -directional variance is shown not overtaking the 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/34d8e6476f31153df933acad895116707af5f9e641c09a073338200ffb1b3e6f.jpg)



Figure 4.4: Comparison of the ratio $\sigma _ { z z } ^ { \mathrm { D T } } / \sigma _ { x x } ^ { \mathrm { D T } }$ for other modes $\ell = 3 - 1 2$ to investigate the transition phenomenon.


$z$ -directional variance for 1%-14% initial velocity perturbations. 2-D modes with $m = 0$ are shown to be more anisotropic, with $\langle \sigma _ { z z } ^ { \mathrm { D  T } } \rangle / \langle \sigma _ { x x } ^ { \mathrm { D T } } \rangle > 2$ , than 3-D modes with $m \neq 0$ , because the expanding radial flow structure for 3-D modes appears to be spherical shape whereas the expanding radial flow structure for 2-D modes appears to be ring shape to conserve the rotational symmetry about the $z$ -axis. The branch of 3-D modes lies inside the isotropic source regime defined by $\mathcal { R } _ { \mathrm { i s o } }$ : $1 \leq \langle \sigma _ { z z } ^ { \mathrm { D 1 } } \rangle / \langle \sigma _ { x x } ^ { \mathrm { D 1 } } \rangle \leq 2$ , which results in 3-D $m \neq 0$ modes being more isotropic source than 2-D $m = 0$ modes. 2-D modes $\ell = 5 - 6$ are shown to approach $\mathcal { R } _ { \mathrm { i s o } }$ at the cost of large perturbations. 

The conditions to maximize the isotropic velocity variance in Eq. (4.24) occurs when all directional-variance has equal and large magnitudes, which restricts the anisotropic velocity variance to be zero simultaneously, 

$$
\sigma_ {\text {i s o}} ^ {\max } = \operatorname {M i n} \left[ \left\{\sigma_ {i i} \right\} \gg 0 \right] \text {a n d} \left\{\triangle \sigma_ {i i} \right\} = 0, \tag {4.30}
$$

and the minimum ion-temperature among different LOS by Eq. (4.23) is, 

$$
\hat {T} _ {\mathrm {i}, \text {s i n g l e - m o d e}} ^ {\text {i n f e r r e d}} = \hat {T} _ {\mathrm {i}} ^ {\text {t h e r m a l}} + \sigma_ {\text {i s o}}. \tag {4.31}
$$

Mode $\ell = 2$ due to the transition property of $\sigma _ { x x } = \sigma _ { z z }$ , and 3-D $m \neq 0$ modes due to lying inside the isotropic source regime $\mathcal { R } _ { \mathrm { i s o } }$ satisfy the maximization condi-

tions, and contribute to large isotropic source in neutron-inferred ion-temperature measurements. The summary of single-mode isotropic velocity variance is, 

$$
\sigma_ {\mathrm {i s o}} ^ {\text {l a r g e}}: \ell = 2 \rightarrow \text {m o d e s t}, \tag {4.32}
$$

$$
\sigma_ {\text {i s o}} ^ {\text {m o d e r a t e}}: m \neq 0 \rightarrow 3 \mathrm {D} \text {m o d e s}, \tag {4.33}
$$

$$
\sigma_ {\text {i s o}} ^ {\text {s m a l l}}: m = 0. \tag {4.34}
$$

The third class specified by $m = 0$ represents rotational symmetric modes about the $z$ -axis including 2-D even- $\ell$ and odd- $\ell$ modes. Growth factors for RT spikes along the $z$ -axis is proportional to Legendre functions $P _ { \ell } ( \cos \theta )$ , which are peaked at the north $\theta = 0$ and south $\theta = \pi$ poles. As a result, fast-growing RT spikes along the $z$ -axis drive more non-translational residual kinetic energies in the $z$ -direction than that in $x$ or $y$ directions, which are driven by rotational-symmetric ring-structure RT spikes, so that 2-D $m = 0$ modes behave more anisotropic than 3- D $m \neq 0$ modes due to $\sigma _ { z z } \gg \sigma _ { x x }$ and $\sigma _ { z z } \gg \sigma _ { y y }$ . The mechanism to produce large isotropic source must require breaking the rotational symmetry along the $z$ -axis by introducing $m \neq 0$ modes to form more three-dimensionally expanding warm bubbles to raise the non-translational residual kinetic energy $\sigma _ { x x }$ and $\sigma _ { y y }$ so that the class of 3-D $m \neq 0$ modes are moderate isotropic. Mode $\ell = 2$ is the special case among rotational symmetric modes because the rapid expanding donut-shape warm bubble has fast-growing directional variance in $\sigma _ { x x }$ and $\sigma _ { y y }$ , which eventually overtake the non-translational residual kinetic $\sigma _ { z z }$ driven by the pair of RT spikes along the $z$ -axis in large perturbations. Around the transition from modest to large perturbations, mode $\ell = 2$ satisfies the maximization condition to develop a large isotropic velocity variance according Eq. (4.24). 

# 4.4 Effects on DD/DT ion-temperature ratio

The exact expression of isotropic velocity variance in Eq. (4.24) is applied to investigate the effect of growing isotropic non-stagnating hot-spot fluid motion on inferring minimum of DD and DT ion temperatures. 

$$
T _ {\mathrm {i , m i n}} ^ {\mathrm {D D}, \text {i n f e r r e d}} (t) = T _ {\mathrm {i , D D}} ^ {\text {t h e r m a l}} (t) + \left(m _ {\mathrm {n}} + m _ {\mathrm {3 H e}}\right) \sigma_ {\mathrm {i s o}} ^ {\mathrm {D D}} (t), \tag {4.35}
$$

$$
T _ {\mathrm {i , m i n}} ^ {\mathrm {D T , i n f e r r e d}} (t) = T _ {\mathrm {i , D T}} ^ {\mathrm {t h e r m a l}} (t) + \left(m _ {\mathrm {n}} + m _ {\alpha}\right) \sigma_ {\mathrm {i s o}} ^ {\mathrm {D T}} (t). \tag {4.36}
$$

$\sigma _ { \mathrm { i s o } } ^ { \mathrm { D D } }$ and $\sigma _ { \mathrm { i s o } } ^ { \mathrm { D T } }$ are DD and DT isotropic velocity variance defined by Eq. (4.24) respectively. Bosch-Hale fusion reactivities are used to calculate DD and DT burnaveraged brackets h...iDD/DT to define the directional-variance in Eq. (4.18). The time-dependency is used to capture the transition from modest to large singlemode perturbations. Any quantity without the notation of “ $( t ) ^ { , }$ is assumed obtained by time-integrated burn-averaging. The ratio of DD to DT time-integrated burn-averaged minimum neutron-inferred ion temperatures in Eq. (4.35) is 

$$
\frac {T _ {\mathrm {i} , \min } ^ {\mathrm {D D} , \text {i n f e r r e d}}}{T _ {\mathrm {i} , \min } ^ {\mathrm {D T} , \text {i n f e r r e d}}} = \left(1 + \frac {\sigma_ {\mathrm {i s o}} ^ {\mathrm {D D}}}{\hat {T} _ {\mathrm {i} , \mathrm {D D}} ^ {\text {t h e r m a l}}}\right) \left(1 + \frac {\sigma_ {\mathrm {i s o}} ^ {\mathrm {D T}}}{\hat {T} _ {\mathrm {i} , \mathrm {D T}} ^ {\text {t h e r m a l}}}\right) ^ {- 1}, \tag {4.37}
$$

where DD and DT normalized thermal ion temperatures are defined as $\hat { T } _ { \mathrm { i , D D } } ^ { \mathrm { t h e r m a l } } =$ $T _ { \mathrm { i , D D } } ^ { \mathrm { t h e r m a l } } / ( m _ { \mathrm { n } } + m _ { ^ { 3 } \mathrm { H e } } )$ and $\hat { T } _ { \mathrm { i , D T } } ^ { \mathrm { t h e r m a l } } = T _ { \mathrm { i , D T } } ^ { \mathrm { t h e r m a l } } / ( m _ { \mathrm { n } } + m _ { \alpha } )$ respectively. The ratio of DD to DT minimum neutron-inferred ion temperatures is a function of $f _ { \mathrm { r k e } }$ , which is defined as the ratio of DT isotropic velocity variance to DT thermal ion temperature, 

$$
f _ {\mathrm {r k e}} = (m _ {\mathrm {n}} + m _ {\alpha}) \sigma_ {\mathrm {i s o}} ^ {\mathrm {D T}} / T _ {\mathrm {i , D T}} ^ {\mathrm {t h e r m a l}}. \tag {4.38}
$$

Omit the notation $\mathrm { ( . . . ) _ { i } ^ { i n f e r r e d } }$ , an approximate expression for DD/DT minimum ion-temperature ratio in Eq. (4.37) is, 

$$
T _ {\mathrm {m i n}} ^ {\mathrm {D D}} / T _ {\mathrm {m i n}} ^ {\mathrm {D T}} = \left(1 + f _ {\mathrm {r k e}} \cdot \frac {m _ {\mathrm {n}} + m _ {\mathrm {H e} ^ {3}}}{m _ {\mathrm {n}} + m _ {\alpha}} \cdot \frac {\sigma_ {\mathrm {i s o}} ^ {\mathrm {D D}}}{\sigma_ {\mathrm {i s o}} ^ {\mathrm {D T}}}\right) (1 + f _ {\mathrm {r k e}}) ^ {- 1}, \tag {4.39}
$$

assuming a local thermodynamic equilibrium (LTE) between DD and DT ions $T _ { \mathrm { i , D D } } ^ { \mathrm { t h e r m a l } } = T _ { \mathrm { i , D T } } ^ { \mathrm { t h e r m a l } }$ . The term $f _ { \mathrm { r k e } }$ is related to the DT non-translational isotropic hot-spot residual kinetic energy $\mathrm { K E _ { h s , D T } ^ { n o n t r a n s } } = M _ { \mathrm { h s } } \sigma _ { \mathrm { i s o } } ^ { \mathrm { D T } } / 2$ by taking Eq. (4.20). 

$$
f _ {\mathrm {r k e}} = \left(\frac {6}{1 - P _ {\mathrm {e}} / P _ {\mathrm {h s}}} \cdot \frac {m _ {\mathrm {n}} + m _ {\alpha}}{m _ {\mathrm {D}} + m _ {\mathrm {T}}}\right) \cdot \frac {\mathrm {K E} _ {\mathrm {h s , D T} , i} ^ {\mathrm {n o n t r a n s}}}{\mathrm {I E} _ {\mathrm {h s}}}, \tag {4.40}
$$

where $P _ { \mathrm { h s } } = P _ { \mathrm { i } } + P _ { \mathrm { e } }$ is the total hot-spot pressure, $P _ { \mathrm { { e } } }$ is the electron pressure, $\mathrm { { I E } _ { \mathrm { { h s } } } = }$ $\frac { 3 } { 2 } P _ { \mathrm { h s } } V _ { \mathrm { h s } }$ is the total hot-spot internal energy, $V _ { \mathrm { h s } }$ is the hot-spot volume. The value within the bracket equals to 12 in Eq. (4.40), assuming LTE between electrons and ions $P _ { \mathrm { e } } = P _ { \mathrm { i } }$ , and a small mass deficit in DT fusion reactions $m _ { \mathrm { n } } + m _ { \alpha } = m _ { \mathrm { D } } + m _ { \mathrm { T } }$ . The conversion of $f _ { \mathrm { r k e } }$ into Murphy’s definition [Mur14] $f _ { \mathrm { r k e } } ^ { \mathrm { M } } = E _ { \mathrm { k } } ^ { \mathrm { M } } / E _ { \mathrm { t h } } ^ { \mathrm { M } }$ is given by, 

$$
f _ {\mathrm {r k e}} ^ {\mathrm {M}} = f _ {\mathrm {r k e}} / 4, \tag {4.41}
$$

which is obtained by substituting ${ \cal E } _ { \mathrm { k } } ^ { \mathrm { M } } \ = \ \ \frac 3 2 M _ { \mathrm { h s } } \sigma _ { \mathrm { i s o } } ^ { \mathrm { D ^ { \prime } I } }$ and $E _ { \mathrm { t h } } ^ { \mathrm { M } } \ = \ \ { \frac { 3 } { 2 } } ( n _ { \mathrm { D } } + n _ { \mathrm { T } } \ +$ $n _ { \mathrm { e } } ) T _ { \mathrm { i } } ^ { \mathrm { t h e r m a l } } V _ { \mathrm { h s } }$ . The hot-spot mass-volume relation is $M _ { \mathrm { h s } } = ( m _ { \mathrm { D } } + m _ { \mathrm { T } } ) n _ { \mathrm { i } } V _ { \mathrm { h s } } / 2$ , and the total ion number density $n _ { \mathrm { i } } = n _ { \mathrm { D } } + n _ { \mathrm { T } }$ is assumed equal to the electron number density $n _ { \mathrm { e } }$ for a fully ionized plasma. A factor of 3 is required to define Murphy’s total hot-spot kinetic energy $E _ { \mathrm { k } } ^ { \mathrm { M } } = 3 \mathrm { K } \mathrm { E } _ { \mathrm { h s } , \mathrm { D T } , i } ^ { \mathrm { n o n t r a n s } }$ , because the isotropic velocity variance $\sigma _ { \mathrm { i s o } }$ in Eq. (4.24) is defined for the non-translational residual kinetic energy in one direction only. The ratio in Eq. (4.39) exhibits a lower bound 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/82ae90f11b9c0a98c36102b950866d1a8eb14c72c9d2af45faa8d038daf40167.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/2a9a2e3bc8ad9781058386a57964ced82fd72383a1443bfa4425bd05f8f09e60.jpg)



Figure 4.5: (a) Comparison of 1-D profiles of normalized total ion number density and normalized thermal ion temperature for OMEGA shot 77068 at stagnation. The blue solid line indicates the range of influence caused short-wavelength burn-weight factor mode perturbations from within the region $Y _ { \mathrm { c e l l } } ^ { \mathrm { D D } } / Y _ { \mathrm { t o t a l } } ^ { \mathrm { D D } } > Y _ { \mathrm { c e l l } } ^ { \mathrm { D  T } } / Y _ { \mathrm { t o t a l } } ^ { \mathrm { D ^ { \prime } I } }$ $\mathcal { R } _ { \mathrm { e d g e } } : r _ { c } < r < R _ { \mathrm { h s } }$ . The blue dashed line indicates the due to a higher range of influence caused by long-wavelength mode perturbations within the region $\mathcal { R } _ { \mathrm { c o r e } } : r < r _ { c }$ due to a higher burn-weight factor $Y _ { \mathrm { c e l l } } ^ { \mathrm { D 1 T } } / Y _ { \mathrm { t o t a l } } ^ { \mathrm { D 1 } } > Y _ { \mathrm { c e l l } } ^ { \mathrm { D 1 D } } / Y _ { \mathrm { t o t a l } } ^ { \mathrm { D 1 D } }$ > Y D / D DD . (b) Comparison of DD and DT burn-weight factor defined by $Y _ { \mathrm { c e l l } } ( r ) / Y _ { \mathrm { t o t a l } }$ as function of radius. Arbitrary unit is used in the $y$ -axis to show the difference of burn-weight factors essentially.


in the limit of a large fraction of residual kinetic energy, 

$$
\lim  _ {f _ {\mathrm {r k e}} \rightarrow \infty} \left(\frac {T _ {\operatorname* {m i n}} ^ {\mathrm {D D}}}{T _ {\operatorname* {m i n}} ^ {\mathrm {D T}}}\right) _ {\text {f l o o r}} = \frac {m _ {\mathrm {n}} + m _ {3 \mathrm {H e}}}{m _ {\mathrm {n}} + m _ {\alpha}} \cdot \frac {\sigma_ {\mathrm {i s o}} ^ {\mathrm {D D}}}{\sigma_ {\mathrm {i s o}} ^ {\mathrm {D T}}}. \tag {4.42}
$$

By taking the ratio of DD to DT fusion product mass $R _ { m } ^ { \mathrm { D D / D T } } \equiv ( m _ { \mathrm { n } } + m _ { ^ 3 \mathrm { H e } } ) / ( m _ { \mathrm { n } } +$ mα) = 0.8 and assume the ratio of DD to DT isotropic velocity variance RDD/DTσ $m _ { \alpha } ) = 0 . 8$ $R _ { \sigma } ^ { \mathrm { { D D / D T } } } \equiv$ $\sigma _ { \mathrm { i s o } } ^ { \mathrm { D D } } / \sigma _ { \mathrm { i s o } } ^ { \mathrm { D T } } = 1$ , the lower bound equals to $( T _ { \mathrm { m i n } } ^ { \mathrm { D D } } T _ { \mathrm { m i n } } ^ { \mathrm { D T } } ) _ { \mathrm { f l o o r } } = 0 . 8$ . 

Figure (4.5)-(a) shows the 1-D stagnation ion number density and ion thermal temperature for OMEGA shot 77068. DD and DT fusion reactivities are obtained by fitting Bosch-Hale table over ion temperatures 1-5 keV to approximate $\langle \sigma v \rangle _ { \mathrm { D D } } ~ = ~ S _ { \mathrm { D D } } T _ { \mathrm { k e V } } ^ { \beta _ { \mathrm { D D } } } s ^ { - 1 } \mathrm { c m } ^ { 3 }$ and $\langle \sigma v \rangle _ { \mathrm { D T } } ~ = ~ { \cal S } _ { \mathrm { D T } } T _ { \mathrm { k e V } } ^ { \beta _ { \mathrm { D T } } } s ^ { - 1 } \mathrm { c m } ^ { 3 }$ respectively, where constants are $S _ { \mathrm { D D } } ~ = ~ 3 . 9 \times 1 0 ^ { - 2 2 }$ , $S _ { \mathrm { D T } } ~ = ~ 2 . 8 \times 1 0 ^ { - 2 0 }$ , $\beta _ { \mathrm { D D } } ~ = ~ 3 . 4$ and $\beta _ { \mathrm { D T } } = 3 . 9$ . Figure (4.5)-(b) shows that a transition of the DD burn-weight factor $Y _ { \mathrm { c e l l } } ^ { \mathrm { D D } } / Y _ { \mathrm { t o t a l } } ^ { \mathrm { D D } } \ = \ n _ { \mathrm { D } } n _ { \mathrm { D } } \langle \sigma v \rangle _ { \mathrm { D D } } 4 \pi r ^ { 2 } \triangle r / Y _ { \mathrm { t o t a l } } ^ { \mathrm { D D } }$ DD/ from being below to above DT burnweight factor $Y _ { \mathrm { c e l l } } ^ { \mathrm { D T } } / Y _ { \mathrm { t o t a l } } ^ { \mathrm { D T } } = n _ { \mathrm { D } } n _ { \mathrm { T } } \langle \sigma v \rangle _ { \mathrm { D T } } 4 \pi r ^ { 2 } \triangle r / Y _ { \mathrm { t o t a l } } ^ { \mathrm { D T } }$ at a critical radius $r _ { \mathrm { c } }$ . The spatial profiles for burn-weight factors in Fig. (4.5)-(b) are characterized by the 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/1d22cec5a7bd75ecca201da6cf3d1ef7d389a388c0c96da4b108f019c6366545.jpg)



Z(um) L> x(um)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/04f34430be73870d3ce62d966c80b1c48fb2382c76a96c4d5a0931eb4ad3044e.jpg)



Ti=1keV


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/6ec79de8421a34beaf9e70a548f360377a0bf4fbcf28c056c69b7852f66df1a9.jpg)



Figure 4.6: Comparison of fluid velocity field vectors on $x - z$ planes for low modes $\ell = 1 - 2$ and a high mode $\ell = 1 2 , m = 6$ . Contours of electron temperatures at are shown to outline the region of high burn-weights. All plots are in the same spatial scale shown by the rulers of $1 0 \mu \mathrm { m }$ in the corner of mode $\ell = 1$ plot. The color contour shows the spatial profile of $n _ { \mathrm { i } } ^ { 2 } T _ { \mathrm { i } } ^ { \beta } r ^ { 2 }$ , where $\beta = 3 . 9$ is taken to study the spatial distribution of DT burn-weight factor. The increasing burn weights are shown by bright regions closer to hot spot interface within warm bubbles only.


shape function $F _ { \mathrm { s } } ( r ) = n _ { \mathrm { i } } ^ { 2 } T _ { \mathrm { i , t h } } ^ { \beta } r ^ { 2 }$ by omitting the proportionality constant $S _ { \mathrm { D D / D T } }$ , the grid size in the radial direction $\triangle r$ and the total yield $Y _ { \mathrm { t o t a l } } ^ { \mathrm { D D / D T } }$ . The critical radius is obtained by solving $F _ { \mathrm { s } } ^ { \prime } ( r ) = 0$ , where the superscript of prime $\begin{array} { r } { ( \ldots ) ^ { \prime } = \frac { d } { d r } ( \ldots ) } \end{array}$ denotes for the spatial derivative with respect to the radius. 

$$
r _ {\mathrm {c}} = 1 / \left[ 0. 5 \beta L _ {T} ^ {- 1} \left(r _ {\mathrm {c}}\right) - L _ {n} ^ {- 1} \left(r _ {\mathrm {c}}\right) \right], \tag {4.43}
$$

where $L _ { T } ^ { - 1 } ( r ) = | T _ { \mathrm { i } } ^ { \prime } / T _ { \mathrm { i } } |$ and $L _ { n } ^ { - 1 } ( r ) = | n _ { \mathrm { i } } ^ { \prime } / n _ { \mathrm { i } } |$ are thermal ion temperature and ion number density gradient scale lengths in the radial direction respectively. The peaks in Fig. (4.5)-(b) through solving Eq. (4.43) are located at $r _ { \mathrm { c } } ^ { \mathrm { D D } } = 1 5 . 7$ µm and $r _ { \mathrm { c } } ^ { \mathrm { D T } } = 1 4 . 7 ~ \mu m$ . The small difference in exponents $\beta$ of temperature dependence in fusion reactivities results in a small shift in position of the maximum burn-weight, leading to slightly different DD and DT burn-averaged brackets. For example, the ratio of burn-averaged DD to DT thermal ion temperatures based on 1-D ion number density and ion temperature profiles in Fig. (4.5)-(a) is 0.971. 

Effects of different DD and DT fusion reactivities on burn-averaged isotropic velocity variance are, 

$$
\sigma_ {i i} ^ {\mathrm {D D}} / \sigma_ {i i} ^ {\mathrm {D T}} > 1, \quad \text {w h e n} \quad \operatorname {M a x} [ \triangle v _ {i} ^ {2} (r _ {c} \leq r \leq R _ {\mathrm {h s}}) ], \tag {4.44}
$$

$$
\sigma_ {i i} ^ {\mathrm {D D}} / \sigma_ {i i} ^ {\mathrm {D T}} <   1, \quad \text {w h e n} \quad \operatorname {M a x} [ \triangle v _ {i} ^ {2} (r \leq r _ {c}) ], \tag {4.45}
$$

where $R _ { \mathrm { h s } }$ is the hot-spot radius at the 1-keV electron temperature contour surface. Figure (4.6) plots the fluid velocity field pattern on $x - z$ plane between low modes $\ell = 1 - 2$ and a high mode $\ell = 1 2 , m = 6$ on the top of the DT burn weight shape function to visualize that the same hot-spot fluid velocity distribution is weighted differently in space. For high modes, hot-spot residual kinetic energies are highly localized within the cold bubbles in form of high-velocity vorticity, implying that non-translation fluctuation velocities $\triangle v _ { i } ^ { 2 }$ are maximized in the region closed to the hot-spot edge $\mathcal { R } _ { \mathrm { e d g e } } : r _ { c } \leq r \leq R _ { \mathrm { h s } }$ . As a result, the burn-averaged value for $\triangle v _ { i } ^ { 2 }$ is higher in DD than DT in the vicinity of 3-D hot-spot boundaries so that high modes favor the condition of $\sigma _ { i i } ^ { \mathrm { D D } } / \sigma _ { i i } ^ { \mathrm { D T } } > 1$ . For low modes such as modest mode $\ell = 2$ perturbations, the vortex structure is localized within the region of hot core $\mathcal { R } _ { \mathrm { c o r e } } : r < r _ { c }$ so that the burn-averaged value for $\triangle v _ { i } ^ { 2 }$ is higher for DT than DD and favors the condition of $\sigma _ { i i } ^ { \mathrm { D D } } / \sigma _ { i i } ^ { \mathrm { D T } } < 1$ . By taking the first order approximation for the ratio $T _ { \mathrm { m i n } } ^ { \mathrm { L D D } } / T _ { \mathrm { m i n } } ^ { \mathrm { D } ^ { \prime } \mathrm { I } ^ { \prime } }$ in Eq. (4.39), 

$$
T _ {\mathrm {m i n}} ^ {\mathrm {D D}} / T _ {\mathrm {m i n}} ^ {\mathrm {D T}} \simeq 1 + f _ {\mathrm {r k e}} \left(\frac {m _ {\mathrm {n}} + m _ {\mathrm {3 H e}}}{m _ {\mathrm {n}} + m _ {\alpha}} \cdot \frac {\sigma_ {\mathrm {i s o}} ^ {\mathrm {D D}}}{\sigma_ {\mathrm {i s o}} ^ {\mathrm {D T}}} - 1\right). \tag {4.46}
$$

The properties of $\sigma _ { i i } ^ { \mathrm { D D } } / \sigma _ { i i } ^ { \mathrm { D T } } > 1$ and $< ~ 1$ are in effect to perturb the bracket $( R _ { m } ^ { \mathrm { D D / D T } } R _ { \sigma } ^ { \mathrm { D D / D T } } - 1 )$ slightly above or below the value of $- 0 . 2$ . 

Figure (4.7) compares the ratio of DD to DT minimum neutron-inferred ion temperatures, $T _ { \mathrm { m i n } } ^ { \mathrm { D D } } / T _ { \mathrm { m i n } } ^ { \mathrm { D T } } = \ : T _ { \mathrm { i , m i n } } ^ { \mathrm { D D , i n t e r r e d } } / T _ { \mathrm { i , m i n } } ^ { \mathrm { D T , i n t e r r e d } }$ , against the fraction of nontranslational isotropic hot-spot residual kinetic energy $f _ { \mathrm { r k e } }$ in the single-mode 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/302f64238c36badc9b0ec1d64daaa2553032bb445b34d7776045da8a299abb9d.jpg)



Figure 4.7: Comparison of the ratio T DDmin/T DTmin against the fraction of residual kinetic energy $f _ { \mathrm { r k e } }$ for the mode spectrum $\ell = 1 - 1 2$ over all levels of 1%-14% initial velocity perturbations. The back solid line is Eq. (4.39) by taking the RDD/DTσ reaction product mass ratio $R _ { \sigma } ^ { \mathrm { { D D / D T } } } = 1$ . The blue arrows show the trend of decreasing DD/DT minimum $R _ { m } ^ { \mathrm { { D D / D T } } } = 0 . 8$ and isotropic velocity variance ratio inferred ion-temperature ratio as a result of increasing isotropic velocity variance.


spectrum. A trend of decreasing DD/DT ion-temperature ratio with isotropic velocity variance is observed. Points deviate from the black solid curve indicate effects of different DD and DT isotropic velocity variance $\sigma _ { \mathrm { i s o } } ^ { \mathrm { D D } } \neq \sigma _ { \mathrm { i s o } } ^ { \mathrm { D } ^ { \prime } \mathrm { I } ^ { \prime } }$ . The ratio of $T _ { \mathrm { m i n } } ^ { \mathrm { L D D } } / T _ { \mathrm { m i n } } ^ { \mathrm { D T } }$ for modest mode $\ell = 2$ perturbations are shown to below the black curve, because the large radial flow structure is localized within the region of hot core ${ \mathcal { R } } _ { \mathrm { c o r e } }$ , where DT neutron rates are higher than DD’s. For large mode $\ell = 2$ and other mid- and high-mode perturbations, the ratio of $T _ { \mathrm { m i n } } ^ { \mathrm { L D D } } / T _ { \mathrm { m i n } } ^ { \mathrm { L D } }$ are above the black analytic curve, because short-wavelength velocity perturbations are localized closed to the hot-spot edge $\mathcal { R } _ { \mathrm { e d g e } }$ , where DD neutron rates are higher than DT’s. Effects of isotropic velocity variance are shown to cause DD minimum ion temperatures to exhibit a maximum of 4% below than DT minimum ion temperatures. In OMEGA experiments, the averaged ratio $T _ { \mathrm { m i n } } ^ { \mathrm { L D D } } / T _ { \mathrm { m i n } } ^ { \mathrm { D 1 } }$ was about 0.967, which is closed to the result of Fig. (4.7). 

# 4.5 Multi-mode DD/DT ion-temperature ratio

Figure (4.7), suggests that single-mode perturbations can only result in a small range of $f _ { \mathrm { r k e } } = 0 - 0 . 3$ . Since $f _ { \mathrm { r k e } }$ is the ratio of velocity-square to thermal ion temperature, result of single-mode simulations remains unchanged by hydrodynamic-

equivalent scaling to NIF deceleration phase simulations. Because ion temperature is an intensive parameter whereas perturbed velocities have approximately the same RT growth factors between OMEGA and NIF. 

To approach the limit of $T _ { \mathrm { m i n } } ^ { \mathrm { D D } } / T _ { \mathrm { m i n } } ^ { \mathrm { D T } }  0 . 8$ observed in NIF experiments, more isotropic flows and lower thermal ion temperatures are required. Multi-mode perturbations can satisfy this requirement. Because one can superposition different isotropic velocity variance from non-aggress mode $\ell \ = \ 2$ and 3-D $m \ne 0$ modes as much as possible to produce large multi-mode isotropic velocity variance $\sigma$ multi−modeiso [WBS+18b]. At the same, by tuning up initial velocity perturbations of multi-mode spectra, lower thermal ion temperatures can be achieved by reducing conversion efficiencies of shell kinetic energies into hot-spot internal energies [WBS+18a] due to RT instabilities. 

To explain multi-mode ion-temperature measurement asymmetry, effects of azimuthal asymmetries caused by non-zero covariance terms must be considered to modify the expression for single-mode isotropic velocity variance in Eq. (4.24) by taking the global minimum of Eq. (4.1) with respect to angles $\theta$ and $\phi$ , 

$$
\hat {T} _ {\mathrm {i}} ^ {\text {i n f e r r e d}} = \hat {T} _ {\mathrm {i}} ^ {\text {t h e r m a l}} + \sigma_ {\text {i s o}} + \underbrace {\mathcal {P} _ {i = j} [ g _ {i} g _ {j} \triangle \sigma_ {i j} ]} _ {\geq 0} + \underbrace {\mathcal {P} _ {i \neq j} [ g _ {i} g _ {j} \sigma_ {i j} ]} _ {\leq 0 \text {o r} \geq 0}. \tag {4.47}
$$

Since the third term on the right-hand-side of Eq. (4.47) is positive-definite, the overall minimum isotropic velocity variance for multi-mode perturbations that includes three non-zero covariance terms is, 

$$
\sigma_ {\mathrm {i s o}} ^ {\mathrm {m u l t i - m o d e}} = \sigma_ {\mathrm {i s o}} + \operatorname {M i n} \left[ \mathcal {P} _ {i \neq j} [ g _ {i} g _ {j} \sigma_ {i j} ] \right]. \tag {4.48}
$$

The explicit form of covariance terms in terms of LOS angles $\mathrm { c o v } _ { 1 2 } = 2 g _ { 1 } g _ { 2 } \sigma _ { 1 2 }$ 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/0874848ed401d2ec337a822f3bbcd0ca57d108fa5738b464c5e79e47101a7f0c.jpg)



(a)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/fb7ee50b4e83e3ca3b917265ba3b28b506715a9104afb9b0faa1b3dfc7d9f93a.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/5ea551cb8d875857a4b7a879ebad21007d9bc7f9b6dc19ad9d054ab152e4710c.jpg)



(c)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/010391942bc554a8507f00024369645f8412537a3baf0399e6aac81a85bc659b.jpg)



(d)



Figure 4.8: The full map of three different covariance terms to understand the distributions of maximum and minimum values at all LOS angles $\theta$ and $\phi$ . (a) shows $x$ - $y$ covariance term $\mathrm { c o v } _ { 1 2 } / \sigma _ { 1 2 } = \sin ^ { 2 } \theta \sin 2 \phi$ , (b) shows $y$ - $z$ covariance term $\mathrm { c o v } _ { 2 3 } / \sigma _ { 2 3 } = \sin 2 \theta \sin \phi$ , (c) shows $z$ - $x$ covariance term $\cos _ { 3 1 } / \sigma _ { 3 1 } = \sin 2 \theta \cos \phi$ . (d) overlaps all contour lines with values $\pm 0 . 7 5$ to examine the influence of superposition effect.


$\mathrm { c o v _ { 2 3 } = 2 } g _ { 2 } g _ { 3 } \sigma _ { 2 3 }$ and $\mathrm { c o v } _ { 3 1 } = 2 g _ { 3 } g _ { 1 } \sigma _ { 3 1 }$ inside the square bracket of Eq. (4.48) are, 

$$
\operatorname {c o v} _ {1 2} / \sigma_ {1 2} = \sin^ {2} \theta \sin 2 \phi , \tag {4.49}
$$

$$
\operatorname {c o v} _ {2 3} / \sigma_ {2 3} = \sin 2 \theta \sin \phi , \tag {4.50}
$$

$$
\operatorname {c o v} _ {3 1} / \sigma_ {3 1} = \sin 2 \theta \cos \phi . \tag {4.51}
$$

No analytic form is available to simplify the term $\mathrm { M i n } \left[ g _ { i } g _ { j } \sigma _ { i j } ( 1 - \delta _ { i j } ) \right]$ , because the global minimum depends on values of $\sigma _ { 1 2 }$ , $\sigma _ { 2 3 }$ and $\sigma _ { 3 1 }$ after the superposition $\mathrm { c o v _ { 1 2 } + c o v _ { 2 3 } + c o v _ { 3 1 } }$ . Figures (4.8)-(a,b,c) show the variation of each covariance term in the sky map, in which locations of global minimum of each covariance term are outlined within the regions of contour values equal to $- 0 . 7 5$ . Figure (4.8)-(d) shows that locations of global minimum of each covariance term do not overlap with each other. 

As discussed in Section 4.2, the covariance terms defined by Eq. (4.21) approach zero in the fully turbulent limit. This property reduces the complexity of multi-mode isotropic velocity variance in Eq. (4.48) to seek for large fraction of residual kinetic energy $f _ { \mathrm { r k e } }$ to lower DD/DT minimum neutron-inferred iontemperature ratios in the fully turbulence limit. 

$$
\sigma_ {\mathrm {i s o}} ^ {\text {t u r b u l e n c e}} \rightarrow \sigma_ {\mathrm {i s o}}, \tag {4.52}
$$

meaning that a fully turbulent hot spot is filled with non-translational residual kinetic energies $\sim \langle \triangle v _ { i } ^ { 2 } \rangle$ driven by converging RT spikes and expanding warm bubbles that are homogeneous in space, and has the same functional form of single-mode isotropic velocity variance. Only non-fully turbulent hot spots are described by non-zero covariance terms. 

To understand the transition of increasing isotropic velocity variance from single-mode to multi-mode perturbations, and simultaneously maintain a low level of complexity of flow structures caused by mixing different modes, a series of simulations by mixing the same single mode with random phases are investigated using the initial superposition spectrum: $\begin{array} { r } { A _ { \ell } ^ { m } = \sum _ { i = 1 } ^ { N } ( \triangle v / v _ { 0 } ) Y _ { \ell } ^ { m } ( \theta + \theta _ { i } , \phi + \phi _ { i } ) } \end{array}$ , where $\theta _ { i }$ and $\phi _ { i }$ are random phases for the $i$ -th single mode, $N$ is the total number of modes in superposition, and $\triangle v / v _ { 0 }$ is the initial velocity perturbation. Figure (4.9) shows that DD/DT minimum neutron-inferred ion-temperature ratios exhibit a significant reduction to the level of $\sim 0 . 9$ , as a result of rapid increasing fraction of residual kinetic energy $f _ { \mathrm { r k e } }$ in random-phase mixing simulations. The result of single-mode random-phase simulation suggests that the flow structure of 3-D $m \neq 0$ modes including $Y _ { \ell = 4 } ^ { m = 2 }$ 2 and Y m=3 $Y _ { \ell = 6 } ^ { m = 3 }$ are candidates of large isotropic source to account for the trend of $T _ { \mathrm { m i n } } ^ { \mathrm { D D } } / T _ { \mathrm { m i n } } ^ { \mathrm { D T } } \to 0 . 8$ . 

Figure (4.10) compares the isotropic velocity variance formula in Eq. (4.24) between DD and DT using the same dataset in single-mode random-phase simu-

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/204828cc700c05bc905dfcc47ce3428db998cb8a69ead81a6edf18c4e665b980.jpg)



Figure 4.9: Comparison of DD/DT minimum inferred ion-temperature ratio with the fraction of residual kinetic energy in multi-mode perturbations. The initial spectrum $\begin{array} { r } { A _ { \ell } ^ { m } = \sum _ { i = 1 } ^ { N } ( \triangle v / v _ { 0 } ) Y _ { \ell } ^ { m } ( \theta + \theta _ { i } , \phi + \phi _ { i } ) } \end{array}$ is obtained by superposition of a given single mode with $N$ -set of random phases, where $N = 2 0$ was used. The black solid-dashed line is the analytic curve by Eq. (4.39) with product mass ratio $R _ { m } ^ { \mathrm { { D D / D T } } } = 0 . 8$ and DD/DT isotropic velocity variance ratio $R _ { \sigma } ^ { \mathrm { D D / D T } } = 1$ . The blue solid-dashed line is the same analytic curve by substituting Murphy’s definition of fraction of residual kiprovide a small range of $f _ { \mathrm { r e k } } = 4 f _ { \mathrm { r e k } } ^ { \mathrm { M } }$ Single mode perturbations onlyesponding to a weak degradation $f _ { \mathrm { r e k } } ^ { \mathrm { s i n g l e - m o d e } } : 0 - 0 . 3$ produce large isotropic source require multi-mode perturbations to fill in the hot spot with numerous isotropic flows while significantly degrade the thermal ion temperature to push $f _ { \mathrm { r e k } } = ( m _ { \mathrm { n } } + m _ { \alpha } ) \sigma _ { \mathrm { i s o } } ^ { \mathrm { D T } } / T _ { \mathrm { i } } ^ { \mathrm { t h e r m a l } }$ to transit from single-mode regime f single−mode into multi-mode regime $f _ { \mathrm { r e k } } ^ { \mathrm { s i n g l e - m o d e } }$ fmulti−mode : 0.3 − 1. The legend of $f _ { \mathrm { r e k } } ^ { \mathrm { m u l t i - m o d e } } : 0 . 3 - 1$ result in the least ratio T DDmin/T DTmin different color points is the same as in Fig. (4.7). The purple and green points that $T _ { \mathrm { m i n } } ^ { \mathrm { D D } } / T _ { \mathrm { m i n } } ^ { \mathrm { D T } } \sim 0 . 9$ are single-mode random-phase simulations for mode $Y _ { \ell = 4 } ^ { m = 2 }$ and $Y _ { \ell = 6 } ^ { m = 3 }$ respectively, meaning that their flow structure are highly isotropic under random phase mixings. Data lies on the black curve implying the accuracy of analytic formula by Eq. (4.39).


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/427a2ec90d334eea55efa423fa8fa6e945831639e067f803f700e632dc5b85ad.jpg)



Figure 4.10: Comparison of isotropic velocity variance formula in Eq. (4.24) in the single-mode random-phase simulations between DD and DT. $\sigma _ { \mathrm { i s o } } ^ { \mathrm { D D } }$ and σDT $\sigma _ { \mathrm { i s o } } ^ { \mathrm { D T } }$ isotropic velocity variance are approximately equal at the leading order. Small deviations are caused by the shifts in positions of peak burn-weight between DD and DT.


lations. Although some small deviations are observed caused by the small shift in the peaks of burn weight distribution between DD and DT, DD and DT isotropic velocity variance are shown to be approximately the same at the leading order. 

# 4.6 Diagnosing for hot spot flow isotropy

# 4.6.1 The first approximate closure

The presence of isotropic source in neutron-inferred ion temperatures leads to the challenge of inferring thermal ion temperatures from ion temperature measurements. In this section, three solution strategies are derived to infer thermal ion temperatures and hot-spot isotropies simultaneously from ion-temperature measurements by utilizing the unique expression of DD/DT minimum ion-temperature ratio in terms of the fraction of residual kinetic energy $f _ { \mathrm { r k e } }$ in Eq. (4.39). 

In the first strategy, the technique of six line-of-sight (LOS) ion-temperature measurements [WBS+18b] is applied to solve for the minimum DD and DT neutroninferred ion temperatures in Eq. (4.39), followed by a ratio method to extrapolate $\mathrm { { T } _ { i } ^ { t h e r m a l } }$ and $\sigma _ { \mathrm { i s o } } ^ { \mathrm { D T } }$ using the solution of $f _ { \mathrm { r k e } }$ . This strategy is derived from the observation that the minimum Brysk DT ion temperature is also a function of $f _ { \mathrm { r k e } }$ 

according to Eq. (4.35), 

$$
T _ {\min } ^ {\mathrm {D T}} = T _ {\mathrm {i}, \mathrm {D T}} ^ {\mathrm {t h}} \left(1 + \frac {m _ {\mathrm {n}} + m _ {\alpha}}{T _ {\mathrm {i} , \mathrm {D T}} ^ {\mathrm {t h}}} \cdot \sigma_ {\mathrm {i s o}} ^ {\mathrm {D T}}\right) = T _ {\mathrm {i}, \mathrm {D T}} ^ {\mathrm {t h}} \left(1 + f _ {\mathrm {r k e}}\right), \tag {4.53}
$$

where the superscript $\mathrm { ( . . . ) ^ { t h e r m a l } }$ is replaced with $( \ldots ) ^ { \mathrm { t h } }$ for convenience. The thermal ion temperature $T _ { \mathrm { i , D T } } ^ { \mathrm { t h } } = T _ { \mathrm { m i n } } ^ { \mathrm { D T } } / ( 1 + f _ { \mathrm { r k e } } )$ can be extrapolated from Eq. (4.53), followed by extrapolation of DT isotropic velocity variance $\sigma _ { \mathrm { i s o } } ^ { \mathrm { D T } } = ( T _ { \mathrm { m i n } } ^ { \mathrm { D T } } -$ $T _ { \mathrm { i , D T } } ^ { \mathrm { t h } } ) / ( m _ { \mathrm { n } } + m _ { \alpha } )$ . The solution for $f _ { \mathrm { r k e } }$ is obtained from the DD/DT minimum ion-temperature ratio according to Eq.(4.39), 

$$
f _ {\mathrm {r k e}} = \frac {1 - T _ {\min} ^ {\mathrm {D D}} / T _ {\min} ^ {\mathrm {D T}}}{T _ {\min} ^ {\mathrm {D D}} / T _ {\min} ^ {\mathrm {D T}} - R _ {m} ^ {\mathrm {D D / D T}} R _ {\sigma} ^ {\mathrm {D D / D T}}}. \tag {4.54}
$$

The solution for minimum DD and DT neutron-inferred ion temperatures are obtained by applying the six line-of-sight technique [WBS+18b] for both DD and DT ion-temperature measurements. By utilizing the full properties of velocity variance in the non-relativistic Brysk ion temperatures in Eq. (4.22), the hot-spot flow asymmetry is uniquely characterized by a state vector $\vec { \sigma } ^ { X } = ( \sigma _ { 1 1 } ^ { X } , \sigma _ { 2 2 } ^ { X } , \sigma _ { 3 3 } ^ { X } , 2 \sigma _ { 1 2 } ^ { X } , 2 \sigma _ { 2 3 } ^ { X } , 2 \sigma _ { 3 1 } ^ { X } )$ that contains the six hot-spot flow parameters including three directional-variance $( \sigma _ { 1 1 } ^ { X } , \sigma _ { 2 2 } ^ { X } , \sigma _ { 3 3 } ^ { X } )$ and three covariance $( \sigma _ { 1 2 } ^ { X } , \sigma _ { 2 3 } ^ { X } , \sigma _ { 3 1 } ^ { X } )$ . Six neutron-inferred normalized ion-temperature measurements $\vec { T } _ { 6 } ^ { X } = ( T _ { 1 } ^ { X } , T _ { 2 } ^ { X } , T _ { 3 } ^ { X } , T _ { 4 } ^ { X } , T _ { 5 } ^ { X } , T _ { 6 } ^ { X } ) / ( m _ { \mathrm { n } } + m _ { X } )$ forms an invertible linear system of equations, 

$$
\vec {T} _ {6} ^ {X} = \vec {T} _ {\mathrm {i}, X} ^ {\mathrm {t h}} + \hat {M} _ {\mathrm {L O S}} \cdot \vec {\sigma} ^ {X}, \tag {4.55}
$$

where $X = ^ { 3 }$ He for DD reactions and $X = { ^ { 4 } \mathrm { H e } }$ or $\alpha$ for DT reactions. The matrix elements for the $k$ -th row of the line-of-sight matrix $\hat { M } _ { \mathrm { L O S } }$ is given by $M _ { \mathrm { L O S } } ^ { ( k ) } =$ $( g _ { 1 } g _ { 1 } , g _ { 2 } g _ { 2 } , g _ { 3 } g _ { 3 } , g _ { 1 } g _ { 2 } , g _ { 2 } g _ { 3 } , g _ { 3 } g _ { 1 } ) ^ { ( k ) }$ . The index $k$ running over 1 to 6 corresponds to the $k$ -th LOS so that all geometrical factors $g _ { i } ( \theta _ { k } , \phi _ { k } ) g _ { j } ( \theta _ { k } , \phi _ { k } )$ for matrix elements of $k$ -th row $M _ { \mathrm { L O S } } ^ { ( k ) }$ requires the $k$ -th LOS angles $\theta _ { k }$ and $\phi _ { k }$ . The normalized thermal 

ion temperature vector is defined by $\vec { T } _ { \mathrm { i } , X } ^ { \mathrm { t h } } = T _ { \mathrm { i } , X } ^ { \mathrm { t h } } ( 1 , 1 , 1 , 1 , 1 , 1 ) / ( m _ { \mathrm { n } } + m _ { X } )$ . By measuring DD and DT ion temperatures along the same LOS at six different locations, the following two matrix equations for DD and DT are obtained. 

$$
\vec {T} _ {6} ^ {\mathrm {D T}} = \vec {T} _ {\mathrm {i , D T}} ^ {\mathrm {t h}} + \hat {M} _ {\mathrm {L O S}} \cdot \vec {\sigma} ^ {\mathrm {D T}}, \tag {4.56}
$$

$$
\vec {T} _ {6} ^ {\mathrm {D D}} = \vec {T} _ {\mathrm {i , D D}} ^ {\mathrm {t h}} + \hat {M} _ {\mathrm {L O S}} \cdot \vec {\sigma} ^ {\mathrm {D D}}. \tag {4.57}
$$

From Eq. (4.55), the state vector is, 

$$
\vec {\sigma} ^ {X} = \hat {M} _ {\mathrm {L O S}} ^ {- 1} \cdot \left(\vec {T} _ {6} ^ {X} - \vec {T} _ {\mathrm {i}, X} ^ {\mathrm {t h}}\right). \tag {4.58}
$$

Ion temperatures at six new LOS are $\vec { T } _ { \mathrm { n e w } } ^ { X } = \vec { T } _ { \mathrm { i } , X } ^ { \mathrm { t h } } + \hat { M } _ { \mathrm { n e w } } \cdot \vec { \sigma } ^ { X }$ T i,X , where the state vector is given by Eq. (4.58). 

$$
\vec {T} _ {\text {n e w}} ^ {X} = \left(\hat {I} - \hat {M} _ {\text {n e w}} \cdot \hat {M} _ {\text {L O S}} ^ {- 1}\right) \cdot \vec {T} _ {\mathrm {i}, X} ^ {\text {t h}} + \hat {M} _ {\text {n e w}} \cdot \hat {M} _ {\text {L O S}} ^ {- 1} \cdot \vec {T} _ {6} ^ {X}. \tag {4.59}
$$

The first term $\hat { \delta } = \hat { I } - \hat { M } _ { \mathrm { n e w } } \cdot \hat { M } _ { \mathrm { L O S } } ^ { - 1 }$ is the departure matrix described in Chapter 3, which measures how far the new six LOS’s are away from the original six LOS’s. When new and original LOS’s overlap, the departure matrix is zero. Magnitudes of peak-to-valley differences over $4 \pi$ solid angles for the departure matrix are observed [WBS+18b] to be negligible, leading to an approximate $4 \pi$ reconstruction of neutron-inferred ion temperatures using six line-of-sight ion-temperature measurements. Figure (4.11) compares the peak-to-valley difference of the departure matrix $\hat { \delta } _ { \mathrm { m a x } } - \hat { \delta } _ { \mathrm { m i n } } = \mathrm { M a x } [ \hat { \delta } ( \theta , \phi ) ] _ { 4 \pi } - \mathrm { M i n } [ \hat { \delta } ( \theta , \phi ) ] _ { 4 \pi }$ against the conditional number for various configuration of the 6-th LOS in NIF. The minimum $\hat { \delta } _ { \mathrm { m a x } } - \hat { \delta } _ { \mathrm { m i n } }$ is shown to occur at ( $\theta _ { 6 } = 6 2 , \phi _ { 6 } = 5 1$ ) with conditional number of 21 for $\hat { M } _ { \mathrm { L O S } }$ in NIF. The matrix elements $\delta ( \theta , \phi )$ vary from positive to negative values over $4 \pi$ solid angles for a given LOS $( \theta _ { 6 } , \phi _ { 6 } )$ . Figure (4.11) shows that magnitudes of $\hat { \delta } _ { \mathrm { m a x } } - \hat { \delta } _ { \mathrm { m i n } }$ are small enough to be neglected in all configurations of the 6th LOS 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/25c4e363a7dd16c189609efa99010f0e86fb00af7c5b62ce2df4cd0877faa8da.jpg)



Figure 4.11: The search for the $6 t h$ LOS angles $\theta _ { 6 }$ and $\phi _ { 6 }$ in NIF is obtained by discretizing $\theta$ and $\phi$ angles into 16 and 32 uniform mesh. Current NIF five LOS are NITOF at ( $\theta _ { 1 } = 9 0 , \phi _ { 1 } = 3 1 5$ ) for DT, Spec-A at ( $\theta _ { 2 } = 1 1 6 , \phi _ { 2 } = 3 1 6$ ) for DT and DD, Spec-SP at ( $\theta _ { 3 } = 1 6 1 , \phi _ { 3 } = 5 6$ ) for DT and DD, MRS at ( $\theta _ { 4 } = 7 3 , \phi _ { 4 } = 3 2 4$ ) for DT, and Spec-E at $\theta _ { 5 } = 9 0 , \phi _ { 5 } = 1 7 4 ) ,$ for DT and DD. The peak-to-valley are defined by $\hat { \delta } _ { \operatorname* { m a x } } = \mathrm { M a x } [ \hat { \delta } ( \theta , \phi ) ] _ { 4 \pi }$ and $\hat { \delta } _ { \mathrm { m i n } } = \mathrm { M i n } [ \hat { \delta } ( \theta , \phi ) ] _ { 4 \pi }$ , where the departure matrix is given by $\hat { \delta } = \hat { I } - \hat { M } _ { \mathrm { n e w } } \cdot \hat { M } _ { \mathrm { L O S } } ^ { - 1 }$ .


in NIF. 

$$
\vec {T} _ {\text {n e w}} ^ {X} = \hat {M} _ {\text {n e w}} \cdot \hat {M} _ {\text {L O S}} ^ {- 1} \cdot \vec {T} _ {6} ^ {X}. \tag {4.60}
$$

Using six line-of-sight DD and DT neutron-inferred ion-temperature measurements, an approximate solution with the least number of LOS to the minimum of DD and DT neutron-inferred ion temperatures are, 

$$
T _ {\mathrm {m i n}} ^ {\mathrm {D D / D T}} = \operatorname {M i n} \left[ \hat {M} _ {\mathrm {n e w}} (\theta , \phi) \cdot \hat {M} _ {\mathrm {L O S}} ^ {- 1} \cdot \vec {T} _ {6} ^ {\mathrm {D D / D T}} \right] _ {4 \pi}. \tag {4.61}
$$

Equations (4.53), (4.54) and (4.61) form a closure to the linear system of DD and DT non-relativistic Brysk ion temperatures to extrapolate the thermal ion temperature at the cost to assume $\sigma _ { \mathrm { i s o } } ^ { \mathrm { D D } } = \sigma _ { \mathrm { i s o } } ^ { \mathrm { D ^ { \prime } I } }$ in Eq. (4.54) and neglecting the departure matrix $\hat { \delta }$ in Eq. (4.60). 

$$
T _ {\mathrm {i , D T}} ^ {\text {t h e r m a l}} = T _ {\min } ^ {\mathrm {D T}} / \left(1 + f _ {\mathrm {r k e}}\right). \tag {4.62}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/b6d834507b55afe8177ed177a37472edfe0ae8b28b83bae0a1a19b4191fd6cde.jpg)



(a)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/0902c0827b688623e89ec934a11026f4f6d3a9a3b96d9643b774231525647d14.jpg)



Figure 4.12: (a) The map of standard deviation va $\cdot [ \vec { T } _ { \mathrm { n e w } } ( \theta , \phi ) ] ^ { 1 / 2 }$ for OMEGA six LOS. (b) The positions of re-allocated six LOS to minimize error propagations.


# 4.6.2 Error propagation analysis for the first closure

To analyze the propagation of ion-temperature measurement errors, an error vector for six ion temperature measurements $\delta \vec { T } _ { 6 } = ( \epsilon _ { 1 } , \epsilon _ { 2 } , \epsilon _ { 3 } , \epsilon _ { 4 } , \epsilon _ { 5 } , \epsilon _ { 6 } )$ is introduced in Eq. (4.61). The distribution of measurement error k=1,...,6 is a normal Gaussian with standard deviations equal to the magnitude of error bar in experiments. Let $\vec { T } _ { 6 } ^ { X }$ be the true measurement and $\vec { T } _ { \mathrm { n e w } } ^ { X }$ be ion temperature at new LOS subjected to a vector of random error $\delta \vec { T } _ { 6 } ^ { X }$ . 

$$
\vec {T} _ {\mathrm {n e w}} ^ {X} = \hat {M} _ {\mathrm {n e w}} \cdot \hat {M} _ {\mathrm {L O S}} ^ {- 1} \cdot (\vec {T} _ {6} ^ {X} + \delta \vec {T} _ {6} ^ {X}). \tag {4.63}
$$

By averaging different ensemble of normal Gaussian error distributions, the mean temperature $\langle \vec { T } _ { \mathrm { n e w } } ^ { X } \rangle$ is unaffected because the mean of random error is zero $\langle \delta \vec { T } _ { 6 } ^ { X } \rangle =$ 0. The variance of prediction, which is defined as var $\begin{array} { r } { \Big [ \vec { T } _ { \mathrm { n e w } } ^ { X } \Big ] = \langle \Big ( \vec { T } _ { \mathrm { n e w } } ^ { X } - \langle \vec { T } _ { \mathrm { n e w } } ^ { X } \rangle \Big ) ^ { 2 } \rangle } \end{array}$ is given by, 

$$
\mathrm {v a r} \left[ \vec {T} _ {\mathrm {n e w}} ^ {X} \right] = \langle \left(\hat {M} _ {\mathrm {n e w}} \cdot \hat {M} _ {\mathrm {L O S}} ^ {- 1} \cdot \delta \vec {T} _ {6} ^ {X}\right) ^ {2} \rangle . (4. 6 4)
$$

Figure (4.12-a) shows the sky map of standard deviation $\mathrm { v a r } [ \vec { T } _ { \mathrm { n e w } } ( \theta , \phi ) ] ^ { 1 / 2 }$ for current six LOS in OMEGA: ( $\theta _ { 1 } = 6 1 , \phi _ { 1 } = 4 8$ ) for 15.8m nTOF, ( $\theta _ { 2 } = 6 2 , \phi _ { 2 } =$ 

206) for 15.9m nTOF, ( $\theta _ { 3 } = 1 1 7 , \phi _ { 3 } = 1 6 2$ ) for 13.0m nTOF, $\theta _ { 4 } = 3 8 , \phi _ { 4 } = 2 5 0$ ) for 10.4m nTOF, ( $\theta _ { 5 } = 8 8 , \phi _ { 5 } = 1 6 1$ ) for 12.0m nTOF, ( $\theta _ { 6 } \ = \ 8 5 , \phi _ { 6 } \ = \ 3 1 2$ ) for 5.3m nTOF. $\theta$ and $\phi$ angles are discretized into 16 and 32 uniform mesh. On each grid $( \theta _ { i } , \phi _ { j } )$ , 100 random samples specified by a normal Gaussian distribution $\epsilon ( \mu , \sigma )$ with zero mean $\mu = 0$ and standard deviation $\sigma = 0 . 2$ keV are generated to represent the six components of the error vector $\delta \vec { T } _ { 6 } ( \theta _ { i } , \phi _ { j } ) = ( \epsilon _ { 1 } , \epsilon _ { 2 } , \epsilon _ { 3 } , \epsilon _ { 4 } , \epsilon _ { 5 } , \epsilon _ { 6 } )$ . The new temperature on each grid $T _ { \mathrm { n e w } } ( \theta _ { i } , \phi _ { j } )$ given by Eq. (4.60) contains 100 repeated measurements with Gaussian error distribution. On each grid, the final error propagated onto $T _ { \mathrm { n e w } } ( \theta _ { i } , \phi _ { j } )$ is obtained by ensemble-averaging to obtain the standard deviation $\mathrm { v a r } [ \vec { T } _ { \mathrm { n e w } } ( \theta _ { i } , \phi _ { j } ) ] ^ { 1 / 2 }$ . 

Figure (4.12-a) shows the distribution of error propagation for current OMEGA six LOS obtained by ensemble-averaging of 100 random Gaussian errors to represent the error vector $\delta \vec { T } _ { 6 } ^ { X }$ in Eq. (4.63). The red region shows an intensive error amplification above 1 keV due to nonlinear terms ${ { y } _ { i } } { { y } _ { j } } { { y } _ { n } } { { y } _ { m } }$ rooted from the product of $\hat { M } _ { \mathrm { n e w } } \cdot \hat { M } _ { \mathrm { L O S } } ^ { - 1 }$ in Eq. (4.60). The blue region shows a moderate error propagation between 0.2 keV to 0.5 keV. The purple region shows a small error propagation less than 0.2 keV. 

Figure (4.12-b) shows the positions of re-allocated six LOS to minimize error propagations. The first three LOS are fixed along three orthogonal directions: ( $\theta _ { 1 } = 0 , \phi _ { 1 } = 0$ ) along the $z$ -axis, ( $\theta _ { 2 } = 9 0 , \phi _ { 2 } = 0$ ) along the $x$ -axis and ( $\theta _ { 3 } =$ $9 0 , \phi _ { 3 } = 9 0 $ ) along the $y$ -axis to maximize their range of influence. The rest of three LOS are produced by random selection, and are chosen when the combination of these three LOS yields the least error propagation in the sky map defined by $\mathrm { M a x } [ \mathrm { v a r } [ \vec { T } _ { \mathrm { n e w } } ( \theta , \phi ) ] ^ { 1 / 2 } ]$ . The configuration in Fig. (4.12-b) has a maximum error of 0.47 keV within the red region and a minimum error 0.15 keV within the blue region. The upper bound for error propagator $\hat { M } _ { \mathrm { n e w } } \cdot \hat { M } _ { \mathrm { L O S } } ^ { - 1 }$ in Eq. (4.60) is shown to be controllable by re-allocating six LOS to reduce the nonlinear term ${ { y } _ { i } } { { y } _ { j } } { { y } _ { n } } { { y } _ { m } }$ . 

Figure (4.13)-(a) shows the performance of the 6-LOS method by Eq. (4.61) 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/612f888c545a8e4970072c66fdea5bbe49d15a1a665d402773c3505c625a5a04.jpg)



(a)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/2ebf730e339ec1ea28a3c33003cbadcaee6d4a190c063022bbc8370360b5d209.jpg)



(b)



Figure 4.13: (a) Implementation of 6-LOS method in Eq. (4.61) to extrapolate the true minimum DT inferred ion temperatures in OMEGA experiments according to method 1. The red data are minimum DT ion temperatures over the purple regions in Fig. (4.12-a) that has error propagation less than 0.2 keV. The gray circles are the minimum of 6 nTOF experimental DT ion temperatures. The red data has less spread and shows increasing temperature dependence of experimental yield. (b) Comparison of fitting exponent $b _ { \mathrm { { f i t } } }$ against different tolerance of error propagation. The first data at error of 0.2 keV gives the most robust fitting exponent $b _ { \mathrm { f i t } } = 3 . 4 8$ .


to extrapolate the true minimum DT neutron-inferred ion temperatures over the purple colored region in Fig. (4.12-a) that has the minimum error propagation in the first method. The fitting exponent of minimum DT neutron-inferred iontemperatures with experimental yields is observed to increase from 3.26 to 3.48. Figure (4.13)-(b) shows a decreasing fitting exponent by searching the minimum DT ion temperatures over regions with increasing error propagation. To extrapolate the global minimum of neutron-inferred ion temperatures accurately, the configuration of 6 nTOF can be re-allocated according to Fig. (4.12)-(b) to minimize error propagations. 

# 4.6.3 The second approximate closure

In the second strategy, thermal ion temperatures are extrapolated without the assumption to drop the departure matrix in Eq. (4.60). By taking the scalar products on both sides of Eq. (4.56) with a unit vector defined as ${ \hat { e } } _ { 6 } = ( 1 , 1 , 1 , 1 , 1 , 1 )$ , 

two systems of DD and DT matrix equations are reduced to two scalar equations, 

$$
\vec {T} _ {6} ^ {\mathrm {D T}} \cdot \hat {e} _ {6} = \vec {T} _ {\mathrm {i , D T}} ^ {\mathrm {t h}} \cdot \hat {e} _ {6} + (m _ {\mathrm {n}} + m _ {\alpha}) \hat {M} _ {\mathrm {L O S}} \cdot \vec {\sigma} _ {6} ^ {\mathrm {D T}} \cdot \hat {e} _ {6}, \tag {4.65}
$$

$$
\vec {T} _ {6} ^ {\mathrm {D D}} \cdot \hat {e} _ {6} = \vec {T} _ {\mathrm {i , D D}} ^ {\mathrm {t h}} \cdot \hat {e} _ {6} + (m _ {\mathrm {n}} + m _ {3 \mathrm {H e}}) \hat {M} _ {\mathrm {L O S}} \cdot \vec {\sigma} _ {6} ^ {\mathrm {D D}} \cdot \hat {e} _ {6}. \tag {4.66}
$$

where $\vec { T _ { 6 } ^ { X } } = ( T _ { 1 } , T _ { 2 } , T _ { 3 } , T _ { 4 } , T _ { 5 } , T _ { 6 } ) _ { \mathrm { i } , X } ^ { \mathrm { i n f e r r e d } }$ and $\vec { T } _ { \mathrm { i } , X } ^ { \mathrm { t h } } = T _ { \mathrm { i } , X } ^ { \mathrm { t h } } \cdot \hat { e } _ { 6 }$ are vectors containing six neutron-inferred ion temperatures and thermal ion temperatures without normalized by fusion product masses respectively. The scalar product represents the summation over six components in a compact form. For example, the scalar product $\vec { \sigma } _ { 6 } ^ { X } \cdot \hat { e } _ { 6 } = ( \sigma _ { 1 1 } ^ { X } , \sigma _ { 2 2 } ^ { X } , \sigma _ { 3 3 } ^ { X } , 2 \sigma _ { 1 2 } ^ { X } , 2 \sigma _ { 2 3 } ^ { X } , 2 \sigma _ { 3 1 } ^ { X } ) \cdot \hat { e } _ { 6 }$ . Assume DD and DT have the same state vector $\vec { \sigma } _ { 6 } ^ { \mathrm { D D } } = \vec { \sigma } _ { 6 } ^ { \mathrm { D T } } = \vec { \sigma } _ { 6 }$ and the same thermal ion temperature ~ th $\vec { T } _ { \mathrm { i , D D } } ^ { \mathrm { t h } } = \vec { T } _ { \mathrm { i , D T } } ^ { \mathrm { t h } } = \vec { T } _ { \mathrm { i } } ^ { \mathrm { t h } }$ T~ th , an invertible $2 \times 2$ matrix system $\hat { M } _ { 2 }$ is formed, 

$$
\left[ \begin{array}{c} \vec {T} _ {6} ^ {\mathrm {D T}} \cdot \hat {e} _ {6} \\ \vec {T} _ {6} ^ {\mathrm {D D}} \cdot \hat {e} _ {6} \end{array} \right] = \left[ \begin{array}{c c} 1 & R _ {m} ^ {\mathrm {D T / D D}} \\ 1 & 1 \end{array} \right] _ {\hat {M} _ {2}} \cdot \left[ \begin{array}{c} \vec {T} _ {\mathrm {i}} ^ {\mathrm {t h}} \cdot \hat {e} _ {6} \\ m _ {\mathrm {D D}} \hat {M} _ {\mathrm {L O S}} \cdot \vec {\sigma} _ {6} \cdot \hat {e} _ {6} \end{array} \right], \tag {4.67}
$$

where $m _ { \mathrm { D D } } = m _ { \mathrm { n } } + m _ { ^ { 3 } \mathrm { H e } }$ is the total DD product mass and $R _ { m } ^ { \mathrm { D T / D D } } = ( m _ { \mathrm { n } } +$ $m _ { \alpha } ) / ( m _ { \mathrm { n } } + m _ { ^ 3 \mathrm { H e } } ) \simeq 1 . 2 5$ is the ratio of DT to DD product mass. 

# 4.6.4 The third approximate closure

In the third strategy, using the multi-mode isotropic velocity variance formula given by Eq. (4.48), the neutron-inferred ion temperature at a given LOS, 

$$
\hat {T} _ {\mathrm {i}} ^ {\text {i n f e r r e d}} = \hat {T} _ {\min } + \sigma_ {\text {a n i s o}} (\theta , \phi), \tag {4.68}
$$

is a sum of the minimum ion temperature $\hat { T } _ { \mathrm { m i n } } = \hat { T } _ { \mathrm { i } } ^ { \mathrm { t h e r m a l } } + \sigma _ { \mathrm { i s o } } ^ { \mathrm { m u l t i - m o d e } }$ and a directional-dependent anisotropic velocity variance defined as, 

$$
\sigma_ {\text {a n i s o}} = \mathrm {P} _ {i = j} [ g _ {i} g _ {j} \triangle \sigma_ {i j} ] + \mathrm {P} _ {i \neq j} [ g _ {i} g _ {j} \sigma_ {i j} ] - \operatorname {M i n} [ \mathrm {P} _ {i \neq j} [ g _ {i} g _ {j} \sigma_ {i j} ] ]. \tag {4.69}
$$

From Eq. (4.68), with known DD and DT ion-temperature measurements along the same LOS at one location, 

$$
T _ {\mathrm {D T}} ^ {\mathrm {L O S}} = T _ {\min } ^ {\mathrm {D T}} + \left(m _ {\mathrm {n}} + m _ {\alpha}\right) \sigma_ {\text {a n i s o}} ^ {\mathrm {D T}}, \tag {4.70}
$$

$$
T _ {\mathrm {D D}} ^ {\mathrm {L O S}} = T _ {\min } ^ {\mathrm {D D}} + \left(m _ {\mathrm {n}} + m _ {3 \mathrm {H e}}\right) \sigma_ {\mathrm {a n i s o}} ^ {\mathrm {D D}}. \tag {4.71}
$$

Assume DD and DT have the same anisotropic velocity variance $\sigma _ { \mathrm { a n i s o } } ^ { \mathrm { D D } } = \sigma _ { \mathrm { a n i s o } } ^ { \mathrm { D ^ { \prime } I } }$ σaniso, DT the DD minimum ion temperature can be extrapolated by removing the common anisotropic term $\sigma _ { \mathrm { a n i s o } }$ . 

$$
T _ {\mathrm {m i n}} ^ {\mathrm {D D}} = T _ {\mathrm {D D}} ^ {\mathrm {L O S}} - \underbrace {(T _ {\mathrm {D T}} ^ {\mathrm {L O S}} - T _ {\mathrm {m i n}} ^ {\mathrm {D T}}) R _ {m} ^ {\mathrm {D D / D T}}} _ {\sigma_ {\mathrm {a n i s o}}}, \tag {4.72}
$$

where the ratio of DD to DT fusion product mass is given by $R _ { m } ^ { \mathrm { D D / D T } } = ( m _ { \mathrm { n } } +$ $m _ { \mathrm { { ^ 3 H e } } } ) / ( m _ { \mathrm { { n } } } + m _ { \alpha } )$ , and $T _ { \mathrm { m i n } } ^ { \mathrm { { D T } } }$ is the minimum DT ion temperature extrapolated from the 6-LOS method. From Eq. (4.39), ratios of DD/DT minimum inferred ion temperatures are below unity and approaches to the limit $\sim 0 . 8$ given by Eq. (4.42). As a result, experimental DD minimum neutron-inferred ion temperatures are closer to thermal ion temperatures and exhibit a stronger correlation with experimental yields. 

Figure (4.14-a) shows a significant improvement in the correlation between extrapolated minimum DD ion temperatures and experimental yields in OMEGA experiments by implementing the third strategy using Eq. (4.72). The minimum DT ion temperatures are taken as the minimum among all available nTOF DT ion-temperature measurements $T _ { \mathrm { m i n } } ^ { \mathrm { L D T } } ~ { = } ~ \mathrm { M i n } [ T _ { \mathrm { n T O F } } ^ { \mathrm { D T } } ]$ . The exponent of DD iontemperature dependence with experimental yields is increased significantly from 1.76 for DD ion temperatures measured by 13.4 m nTOF to 3.52 for the extrapolated minimum DD ion temperatures. A less-spread and more tight correlation with experimental yields is observed simultaneously. Figure (4.14-b) shows the 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/a7e6d13812740f32d6a700ca49acb5c6a90cc090b583d8d67a8cee51799afc38.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/f1ecd7298427bcedbdff3df586b4b8989f9ed5fd82bcc5d9fe35ff8ec01a9ff3.jpg)



Figure 4.14: Investigation of the performance of Eq. (4.72) in extrapolating the minimum DD ion temperature in OMEGA experiments. Gray circles are DD ion temperatures measured by 13.4 m nTOF while red circles are extrapolated DD minimum ion temperatures by Eq. (4.72). The minimum DT ion temperatures are taken as: (a) the minimum among all available nTOF’s DT measurements i.e. $T _ { \mathrm { m i n } } ^ { \mathrm { L D T } } = \mathrm { M i n } [ T ^ { \mathrm { D T } } ]$ , and (b) the extrapolated minimum DT ion temperature from the 6-LOS method by searching the minimum over the purple region in Fig. (4.12-a). The sample size is the same as in Fig. (4.13).


correlation between minimum DT ion temperatures extrapolated from the 6-LOS method and experimental yields. The exponent of DD minimum ion temperature dependence is even increased raised from 1.76 to 3.96, but a fewer sample size was available in this analysis. 

# 4.6.5 Integrated performance in multi-mode simulation

A strongly perturbed multi-mode simulation, as described in Figs. (4.15) and (4.16), is examined to provide a comprehensive demonstration about the procedures for three strategies to extrapolate thermal ion temperature, hot-spot isotropic velocity variance and hot-spot anisotropic velocity variance from DD and DT ion-temperature measurements along the same LOS at six different locations. The initial perturbation spectrum contains a uniform 42% initial velocity perturbation for $\ell = 1 - 1 2$ including $m \neq 0$ modes. In this strongly perturbed multi-mode simulation, a large hot-spot flow isotropy with $T _ { \mathrm { m i n } } ^ { \mathrm { L D D } } / T _ { \mathrm { m i n } } ^ { \mathrm { L D T } } = 0 . 9 3 4$ is observed, and simultaneously exhibits a small DT ion-temperature measurement variation with T DTmax $T _ { \mathrm { m a x } } ^ { \mathrm { D T } } / T _ { \mathrm { m i n } } ^ { \mathrm { D T } } = 1 . 1 3$ . 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/caf93a9f94d8ddb5688c8cb4c35db2af0de9aa0695299005a3c95d7101228d8c.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/267e406b438246a85c7743c05e38319ea17d0323c91f84f0a9dbaf9f416d57a4.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/002a092861a44253f0f7679ae0c6528a319c9768494518e523eac3364ef417d6.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/b2901ca6424e1c43750b68f37ac5ac77066fc0e3c7f0af247fcc57e0e013b6a1.jpg)



Figure 4.15: DEC3D multi-mode simulation for a strongly distorted hot-spot. Left is the 3-D electron temperature contour surface at 0.8 keV. DD and DT ion temperatures are inferred along the same LOS at six different locations. Right is the 2-D $x$ - $z$ plane for the hot spot electron temperature at stagnation. Black arrows are hot spot fluid velocity vectors. The size of arrow heads increase with the magnitude of fluid velocities. The red contour line is the electron temperature at 0.55 keV.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/42fe9bd07e1753b91ff41455b4f814148d60fbd1a8f15e96a51efde51348b1a2.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/37e623304b0495f8307b5c31274c581928436dd28cb37ca561fae64283e5f259.jpg)



Figure 4.16: The same $\it { D E C S D }$ multi-mode simulation as described in Fig. (4.15) for Brysk ion temperature in (a) and the neutron-inferred hot-spot flow velocities in (b). A strong correlation is observed between the flow velocity vector of the jet and positions of maximum neutron-inferred ion temperatures.


Figures (4.16)-(a) and (4.16)-(b) show a strong correlation between DT neutroninferred ion temperatures obtained by Eq. (3.30) and neutron-inferred hot-spot flow velocities defined by $\langle \vec { v } \cdot \hat { d } \rangle = g _ { i } \langle v _ { i } \rangle$ . Two bright spots in Fig. (4.16)-(a) are observed align with the same direction with the tail and the head of the jet as shown in Fig. (4.16)-(b). 


Table 4.1: DEC3D multi-mode perturbation (unit for $\sqrt { \sigma _ { i j } }$ is km/s)


<table><tr><td>Variables</td><td>\(T_{\text{i,keV}}^{\text{thermal}}\)</td><td>\(T_{\text{min,keV}}^{\text{inferred}}\)</td><td>\(T_{\text{max,keV}}^{\text{inferred}}\)</td><td>\(\sqrt{\sigma_{11}}\)</td><td>\(\sqrt{\sigma_{22}}\)</td><td>\(\sqrt{\sigma_{33}}\)</td><td>\(\sqrt{\sigma_{12}}\)</td><td>\(\sqrt{\sigma_{23}}\)</td><td>\(\sqrt{\sigma_{31}}\)</td></tr><tr><td>DD</td><td>2.86</td><td>3.47</td><td>3.83</td><td>134</td><td>152</td><td>122</td><td>34.9i</td><td>32.6i</td><td>17.4i</td></tr><tr><td>DT</td><td>2.95</td><td>3.72</td><td>4.20</td><td>135</td><td>154</td><td>122</td><td>36.1i</td><td>34.5i</td><td>8.86i</td></tr></table>

Table (4.1) summarizes thermal and neutron-inferred ion temperatures, and the six flow parameters between DD and DT corresponding to this multi-mode simulation. The imaginary number “i” is used to represent negative covariance. The square-root of directional-variance is above $\sim 1 0 0 \mathrm { k m / s }$ , which is about the half of DT ion-pair center-of-mass thermal velocity defined by $\sqrt { T _ { \mathrm { i , D T } } ^ { \mathrm { t h e r m a l } } / ( m _ { \mathrm { n } } + m _ { \alpha } ) } =$ 238 km/s, indicating the strong influence of the velocity variance term due to nonstagnating fluid motion. Directional-variance between DD and DT show small deviations within 2-4%. Although relative large differences are observed in DD and DT covariance, their absolute magnitudes are ∼ 4 $\times$ smaller than directionalvariance. 

Table (4.2) summarizes the performance to extrapolate thermal ion temperatures from DD and DT ion-temperature measurements along six different LOS in OMEGA as discussed in Fig. (4.12). The extrapolated thermal ion temperature by the 6-LOS method is shown to be closer to the true thermal ion temperature than that by $\hat { M } _ { 2 }$ matrix inversion in Eq. (4.67). Impact of different fusion reactivities is observed to manifest in two different ways in this multi-mode simulation: different neutron-averaged thermal ion temperatures $T _ { \mathrm { i , D D } } ^ { \mathrm { t h e r m a l } } < T _ { \mathrm { i , D T } } ^ { \mathrm { t h e r m a l } }$ in Table (4.2) and different velocity-variance state vectors $\vec { \sigma } ^ { \mathrm { D D } } \neq \vec { \sigma } ^ { \mathrm { D T } }$ in Table (4.1). 

Table 4.2: Performance of extrapolation thermal ion temperatures by inferring DD and DT along the same LOS at six different locations. Each LOS output DD and DT ion temperatures according to Brysk ion temperature formula in Eq. (4.1): $T _ { \mathrm { i } , X } ^ { \mathrm { i n f e r r e d } } = T _ { \mathrm { i } } ^ { \mathrm { t h e r m a l } } + ( m _ { \mathrm { n } } + m _ { X } ) g _ { i } g _ { j } \sigma _ { i j } ^ { X }$ , where “X”denotes $\mathrm { H e ^ { 3 } }$ for DD reactions and $\alpha$ for DT reactions. Values of six hot-spot flow parameters $\sigma _ { i j } ^ { X }$ are shown in Table (4.1). In the method 3, thermal ion temperatures can be extrapolated location given by directly from DD and DT ion temperature measured along the same LOS at one $\hat { M } _ { 2 } ^ { - 1 } \cdot ( T _ { \mathrm { L O S } } ^ { \mathrm { D T } } , T _ { \mathrm { L O S } } ^ { \mathrm { D D } } ) ^ { T } \cdot ( 1 , 0 ) ^ { T }$ , where $T$ denotes for the transpose of a row vector into a column vector. In the method 2, the extrapolated thermal ion temperature is the averaged of six extrapolated thermal ion temperature from each LOS given by $\begin{array} { r } { \frac { 1 } { 6 } \hat { M } _ { 2 } ^ { - 1 } \cdot ( \vec { T } _ { 6 } ^ { \mathrm { D T } } \cdot \hat { e } _ { 6 } , \vec { T } _ { 6 } ^ { \mathrm { D D } } \cdot \hat { e } _ { 6 } ) ^ { T } \cdot ( 1 , 0 ) ^ { T } } \end{array}$ . In method 3, the extrapolated thermal ion temperature is shown closer to the burn or neutron-averaged thermal ion temperature. 

<table><tr><td rowspan="2">T (keV)</td><td rowspan="2">T1</td><td rowspan="2">T2</td><td rowspan="2">T3</td><td rowspan="2">T4</td><td rowspan="2">T5</td><td rowspan="2">T6</td><td colspan="2">Thermal, no-flow</td></tr><tr><td>i,neutron-avg</td><td>Ti,extraoplate</td></tr><tr><td>DD</td><td>3.59</td><td>3.60</td><td>3.62</td><td>3.63</td><td>3.65</td><td>3.77</td><td>2.86</td><td>—</td></tr><tr><td>DT</td><td>3.88</td><td>3.88</td><td>3.92</td><td>3.93</td><td>3.96</td><td>4.12</td><td>2.95</td><td>—</td></tr><tr><td>Tthermali,(method 3)</td><td>2.39</td><td>2.46</td><td>2.38</td><td>2.41</td><td>2.41</td><td>2.33</td><td>—</td><td>—</td></tr><tr><td>Tthermali,(method 2)</td><td>—</td><td>—</td><td>—</td><td>—</td><td>—</td><td>—</td><td>—</td><td>2.40</td></tr><tr><td>Tthermali,(method 1)</td><td>—</td><td>—</td><td>—</td><td>—</td><td>—</td><td>—</td><td>—</td><td>2.49</td></tr></table>

# 4.6.6 Error propagation analysis in the integration test

To understand the physical origin of the small deviation of ∼ 0.5 keV from the true thermal ion temperature observed in Table (4.2), the DT velocity-variance state vector $\vec { \sigma } ^ { \mathrm { D T } } = \vec { \sigma } ^ { \mathrm { D D } } + \triangle \vec { \sigma } ^ { \mathrm { D T / D D } }$ in the exact DT 6-LOS matrix in Eq. (4.56) is expanded in terms of DD velocity-variance state vector $\vec { \sigma } ^ { \mathrm { D D } }$ and a fluctuation vector defined by $\bigtriangleup \vec { \sigma } ^ { \mathrm { D T / D D } } = \vec { \sigma } ^ { \mathrm { D T } } - \vec { \sigma } ^ { \mathrm { D D } }$ . Similarly, the DT thermal ion temperature $T _ { \mathrm { i , D T } } ^ { \mathrm { t h e r m a l } } = T _ { \mathrm { i , D D } } ^ { \mathrm { t h e r m a l } } + \Delta T _ { \mathrm { i , D T / D D } } ^ { \mathrm { t h e r m a l } }$ T thermali,DT/DD is expanded into the DD thermal ion temperature and a fluctuation component defined by $\Delta T _ { \mathrm { i , D T / D D } } ^ { \mathrm { t h e r m a l } } = T _ { \mathrm { i , D T } } ^ { \mathrm { t h e r m a l } } -$ T thermal − $T _ { \mathrm { i , D D } } ^ { \mathrm { t h e r m a l } }$ , so that the normalized DT thermal ion temperature vector is rewritten as $\vec { T } _ { \mathrm { i , D T } } ^ { \mathrm { t h } } = T _ { \mathrm { i , D D } } ^ { \mathrm { t h e r m a l } } \hat { e } _ { 6 } / m _ { \mathrm { D T } } + \triangle T _ { \mathrm { i , D T / D D } } ^ { \mathrm { t h e r m a l } } \hat { e } _ { 6 } / m _ { \mathrm { D T } }$ , where the DT fusion product mass is $m _ { \mathrm { D T } } = m _ { \mathrm { n } } + m _ { \alpha }$ . Two residual error vectors $\vec { \epsilon } _ { T } = ( \triangle T _ { \mathrm { i , D T / D D } } ^ { \mathrm { t h e r m a l } } \hat { e } _ { 6 } \cdot \hat { e } _ { 6 } , 0 )$ and 

$\vec { \epsilon } _ { \sigma } = ( m _ { \mathrm { D T } } \hat { M } _ { \mathrm { L O S } } \cdot \triangle \vec { \sigma } ^ { \mathrm { D T / D D } } \cdot \hat { e } _ { 6 } , 0 )$ are introduced in Eq. (4.67), 

$$
\vec {\mathbf {y}} = \hat {M} _ {2} \cdot \vec {\mathbf {x}} + \vec {\epsilon} _ {T} + \vec {\epsilon} _ {\sigma}, \tag {4.73}
$$

where $\vec { \bf y } = ( \vec { T } _ { 6 } ^ { \mathrm { D T } } \cdot \hat { e } _ { 6 } , \vec { T } _ { 6 } ^ { \mathrm { D D } } \cdot \hat { e } _ { 6 } )$ and $\vec { \bf x } = ( T _ { \mathrm { i , D D } } ^ { \mathrm { t h e r m a l } } \hat { e } _ { 6 } \cdot \hat { e } _ { 6 } , m _ { \mathrm { D D } } \hat { M } _ { \mathrm { L O S } } \cdot \vec { \sigma } ^ { \mathrm { D D } } \cdot \hat { e } _ { 6 } )$ are left and right column vectors in Eq. (4.67). The error propagation is obtained by multiplying both sides of Eq. (eqn:linear system residual error ) with the inverse matrix of $\hat { M } _ { 2 }$ , 

$$
\begin{array}{l} - \hat {M} _ {2} ^ {- 1} \cdot \vec {\epsilon} _ {T} \cdot \hat {e} _ {1} = - \triangle T _ {\mathrm {i , D T / D D}} ^ {\mathrm {t h e r m a l}} \hat {e} _ {6} \cdot \hat {e} _ {6} / \det  [ \hat {M} _ {2} ], (4.74) \\ - \hat {M} _ {2} ^ {- 1} \cdot \vec {\epsilon} _ {\sigma} \cdot \hat {e} _ {1} = - m _ {\mathrm {D T}} \hat {M} _ {\mathrm {L O S}} \cdot \triangle \vec {\sigma} ^ {\mathrm {D T / D D}} \cdot \hat {e} _ {6} / \det [ \hat {M} _ {2} ], (4.75) \\ \end{array}
$$

where the determinant is $\operatorname * { d e t } [ \hat { M } _ { \mathrm { 2 } } ] ~ = ~ - 1 / 4$ by taking the product mass ratio $R _ { m } ^ { \mathrm { { D T / D D } } } = 1 . 2 5$ and ${ \boldsymbol { \hat { e } } _ { 1 } } ~ = ~ ( 1 , 0 )$ is a unit vector to extract the first row element. By direct computation of right hand sides of above equations to obtain the exact values of averaged error propagation $- \hat { M } _ { 2 } ^ { - 1 } \cdot \vec { \epsilon } _ { T } \cdot \hat { e } _ { 1 } / 6 = 0 . 3 8 2 \ \mathrm { k e }$ V and $- \hat { M } _ { 2 } ^ { - 1 } \cdot \vec { \epsilon } _ { \sigma } \cdot \hat { e } _ { 1 } / 6 = 0 . 0 7 8 8 \mathrm { k e V }$ , the DD neutron-averaged thermal ion temperature 2.86 keV can be recovered by adding these two residual errors to Tthermali,(method 2) in Table (4.2). 

# 4.6.7 TDD/TDT analysis

To understand the influence of different DD and DT burn-averaged thermal ion temperatures and velocity variance, the ratio of DD to DT neutron-inferred ion temperatures T LOSDD / $^ { \prime } I$ LOSDT at one given LOS is expanded with respect to DT’s using Eq. (3.30). 

$$
\frac {T _ {\mathrm {D D}} ^ {\mathrm {L O S}}}{T _ {\mathrm {D T}} ^ {\mathrm {L O S}}} = 1 + \frac {\triangle T _ {\mathrm {i} , \mathrm {D D} / \mathrm {D T}} ^ {\text {t h e r m a l}}}{T _ {\mathrm {D T}} ^ {\mathrm {L O S}}} + \frac {g _ {i} g _ {j} \left(m _ {\mathrm {D D}} \sigma_ {i j} ^ {\mathrm {D D}} - m _ {\mathrm {D T}} \sigma_ {i j} ^ {\mathrm {D T}}\right)}{T _ {\mathrm {D T}} ^ {\mathrm {L O S}}}, \tag {4.76}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/fbb34619b38b1b4210f46dfde1ce3d98c18755690a343aa9b512cebc984419ff.jpg)



(a)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/0b5e9e61ec9453474ec9fde3d85210481857f35de930d59386b5cc1427e2a723.jpg)



Figure 4.17: (a) Comparison of relative changes between the second and the third term in Eq. (4.76). (b) Comparison of relative changes between the first and the second term in Eq. (4.77).


where the relative change of thermal ion temperatures is defined as $\Delta T _ { \mathrm { i , D D / D T } } ^ { \mathrm { t h e r m a l } } =$ $T _ { \mathrm { i , D D } } ^ { \mathrm { t h e r m a l } } - T _ { \mathrm { i , D T } } ^ { \mathrm { t h e r m a l } }$ − T thermali,DT . Substitute the expansion of DD velocity variance σi $\sigma _ { i j } ^ { \mathrm { D D } } = \sigma _ { i j } ^ { \mathrm { D T } } +$ σ Dij $\bigtriangleup \sigma _ { i j } ^ { \mathrm { D D / D T } }$ 4σij DD/DT , where the fluctuation is defined as 4σDDij $\bigtriangleup \sigma _ { i j } ^ { \mathrm { D D / D T } } = \sigma _ { i j } ^ { \mathrm { D D } } - \sigma _ { i j } ^ { \mathrm { D T } }$ /DT σij σij , the third term on the right-hand-side of Eq. (4.76) becomes, 

$$
g _ {i} g _ {j} \left[ \left(\frac {m _ {\mathrm {D D}}}{m _ {\mathrm {D T}}} - 1\right) \cdot \frac {m _ {\mathrm {D T}} \sigma_ {i j} ^ {\mathrm {D T}}}{T _ {\mathrm {D T}} ^ {\mathrm {L O S}}} + \frac {m _ {\mathrm {D D}} \triangle \sigma_ {i j} ^ {\mathrm {D D / D T}}}{T _ {\mathrm {D T}} ^ {\mathrm {L O S}}} \right]. \tag {4.77}
$$

As a result, the exact expression for the ratio of DD to DT inferred ion temperatures at arbitrary LOS in Eq. (4.76) contains a leading term caused by the product mass difference $( m _ { \mathrm { D D } } / m _ { \mathrm { D T } } - 1 ) \simeq - 0 . 2$ , which agrees with the result of Eq. (4.46). Two secondary effects $\Delta T _ { \mathrm { i , D D / D T } } ^ { \mathrm { t h e r m a l } }$ and $\bigtriangleup \sigma _ { i j } ^ { \mathrm { D D / D T } }$ are caused by different DD and DT fusion reactivities. 

Figure (4.17)-(a) compares the relative changes caused by the second and third terms in Eq. (4.76) at different LOS angles. The second term is shown to be negative and exhibits a small range of variation $\sim 0 . 4 \%$ , because the burn-averaged DD thermal ion temperature is slightly smaller than DT’s. The third term is shown to be negative, meaning that the velocity variance of DD is smaller than DT at various LOS angles in this multi-mode simulation and produces a more dominant negative slope of $\sim 3 \%$ variation. In this multi-mode simulation, the isotropic velocity variance of DD is smaller than DT’s by $\sim 1 \%$ . Figure (4.17)-(b) 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/7f0adaf6cf9f02c03d0945b0a2e051d00f102f63a60ec73a9850f5c949303780.jpg)



Figure 4.18: Comparison of DD to DT ion-temperature ratios $T _ { \mathrm { D D } } ^ { \mathrm { 1 3 . 4 m } } / T _ { \mathrm { D T } } ^ { \mathrm { P e t a l } }$ measured along the same LOS in OMEGA experiments, indicated by black circles, with the same multi-mode simulation described in Figs. (4.17-a, b and c), indicated by red circles.


shows that the first term in Eq. (4.77) causes $\sim 2 \%$ variation, the second term causes $\sim 1 \%$ variation, and the sum of these two terms explains the total of $\sim 3 \%$ variation of the slope observed in Fig. (4.17)-(a). 

Figure (4.18) compares the ratio of DD ion temperatures inferred by 13.4 m nTOF to DT ion temperatures inferred by Petal nTOF along the same LOS in OMEGA, shown by black circles. The general features are similar to NIF experiments [GJCF $^ +$ 13], showing a DD/DT ratio below unity and a well-define negative slope with DT ion temperatures. The red circles represent the same multi-mode simulation studied in Figs. (4.17), and is interpolated by the blue curve to estimate its range of influence. A similar trend of negative slope, caused by LOS variations driven by two terms in Eq. (4.77), with OMEGA experiments is reproduced. 

Figure (4.19) shows that the pattern of full-map variations of DD/DT ratios given by Eq. (4.76) exhibit a good correlation with neutron-inferred iontemperature and hot-spot flow velocity asymmetries as shown in figures (4.16)-(a) and (4.16)-(b), and are expected to exhibit strong correlations with areal density variations for low modes. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/97a3a7026d4b1265191fe62d03c43de26e1ec6a53364ee1e4a5bce9498e63b85.jpg)



Figure 4.19: (c) The full map variation of DD/DT neutron-inferred ion temperature ratios given by Eq. (4.76). General features are strongly correlated with DT neutron-inferred ion temperature and neutron-inferred hot-spot flow velocity asymmetries.


# 4.7 Effects of hot-spot flow anisotropy on yield degradation

In this section, the residual kinetic energy model [WBS+18a] developed in Chapter 2 is applied to describe the yield degradation in terms of the ratio of the maximum to minimum DT neutron-inferred ion-temperatures for low mode $\ell = 1$ 

$$
\mathrm {Y O C} \simeq (1 - \overline {{\mathrm {R K E}}} _ {\mathrm {t o t}}) ^ {\mu}, \tag {4.78}
$$

where $\mu = 4 . 4 - 5 . 5$ . Equation (4.78) is derived for low modes $\ell = 1 - 6$ in the deceleration phase of ICF implosions assuming a time-invariant hot-spot adiabatic parameter $ { P _ { \mathrm { h s } } }  { V _ { \mathrm { h s } } } ^ { 5 / 3 }$ , which is robust for OMEGA implosions because of weak alpha heatings and radiation losses. For mode $\ell = 1$ , the ratio of maximum to minimum neutron-inferred ion temperatures is dominated by the non-translational hot-spot residual kinetic energy along the direction of the jet parallel to the $z$ -axis, and exhibits the least non-translational residual kinetic energies $\triangle \sigma _ { 1 1 } , \triangle \sigma _ { 2 2 }  0$ in $x$ and $y$ directions parallel to the equatorial plane [WBS+18b, WBS+18a]. The exact 

form of Brysk ion temperature in Eq. (4.23) with zero covariance is approximated as Ti,DT $\hat { T } _ { \mathrm { i , D T } } ^ { \mathrm { i n f e r r e d } } = \hat { T } _ { \mathrm { i , D T } } ^ { \mathrm { t h e r m a l } } + \sigma _ { \mathrm { i s o } } ^ { \mathrm { D T } } + g _ { 3 } g _ { 3 } \triangle \sigma _ { 3 3 } ^ { \mathrm { D T } }$ Tˆthermali,DT + σDTiso TiDT for mode $\ell = 1$ , where ion temperatures are normalized with respect to DT reaction product masses $\hat { T } _ { \mathrm { D T } } = T / ( m _ { \mathrm { n } } + m _ { \alpha } )$ and the geometrical factor is $g _ { 3 } = \cos \theta$ . The ratio of maximum $\hat { T } _ { \mathrm { i , D T } } ^ { \mathrm { m a x } } ( \theta = 0 , \pi )$ and minimum $\hat { T } _ { \mathrm { i , D T } } ^ { \mathrm { m i n } } ( \theta = \pi / 2 )$ neutron-inferred ion-temperature of DT is, 

$$
\left(\frac {T _ {\mathrm {i , D T}} ^ {\operatorname* {m a x}}}{T _ {\mathrm {i , D T}} ^ {\operatorname* {m i n}}}\right) _ {\ell = 1} = \frac {\hat {T} _ {\mathrm {i , D T}} ^ {\operatorname* {t h e r m a l}} + \sigma_ {\mathrm {i s o}} ^ {\operatorname* {D T}} + \triangle \sigma_ {3 3} ^ {\operatorname* {D T}}}{\hat {T} _ {\mathrm {i , D T}} ^ {\operatorname* {t h e r m a l}} + \sigma_ {\mathrm {i s o}} ^ {\operatorname* {D T}}}. \tag {4.79}
$$

By substituting the property of vanishing non-translational residual kinetic energies $\bigtriangleup \sigma _ { 1 1 } = \bigtriangleup \sigma _ { 2 2 } = 0$ in $x$ and $y$ directions into the isotropic velocity variance formula in Eq. (4.24), mode $\ell = 1$ isotropic velocity variance is zero $\sigma _ { \mathrm { i s o } } ^ { \mathrm { D T } } = 0$ and the fluctuation part by definition $\bigtriangleup \sigma _ { 3 3 } ^ { \mathrm { D ^ { \prime } I ^ { \prime } } } = \sigma _ { 3 3 } ^ { \mathrm { D ^ { \prime } I ^ { \prime } } } - \sigma _ { \mathrm { i s o } } ^ { \mathrm { D ^ { \prime } I ^ { \prime } } } = \sigma _ { 3 3 } ^ { \mathrm { D ^ { \prime } I } }$ σ33 DT σiso DT σ33 is equal to the directional-variance in $z$ -direction. 

$$
\left(T _ {\mathrm {i , D T}} ^ {\mathrm {m a x}} / T _ {\mathrm {i , D T}} ^ {\mathrm {m i n}}\right) _ {\ell = 1} \simeq 1 + \sigma_ {3 3} ^ {\mathrm {D T}} / \hat {T} _ {\mathrm {i , D T}} ^ {\mathrm {t h e r m a l}}. (4. 8 0)
$$

We define the total non-translational residual kinetic energy along one direction $f _ { \mathrm { r k e } } ^ { \mathrm { t o t a l } } = f _ { \mathrm { r k e } } + f _ { \mathrm { r k e } } ^ { \mathrm { a n i s o } }$ , where the isotropic non-translational residual kinetic energy $f _ { \mathrm { r k e } }$ is defined in Eq. (4.38) and the anisotropic non-translational residual kinetic energy is defined as, 

$$
f _ {\mathrm {r k e}} ^ {\text {a n i s o}} = \triangle \sigma_ {i i} / \hat {T} _ {\mathrm {i}} ^ {\text {t h e r m a l}}. \tag {4.81}
$$

We substitute the expression for σ33 / $\sigma _ { 3 3 } ^ { \mathrm { D T } } / \hat { T } _ { \mathrm { i , D T } } ^ { \mathrm { t h e r m a l } } = f _ { \mathrm { r k e } } ^ { \mathrm { t o t a l } }$ DT ˆthermal , the large ion-temperature measurement ratio for mode $\ell = 1$ is 

$$
\left(T _ {\mathrm {i}, \mathrm {D T}} ^ {\max } / T _ {\mathrm {i}, \mathrm {D T}} ^ {\min }\right) _ {\ell = 1} \simeq 1 + f _ {\mathrm {r k e}} ^ {\text {t o t a l}}, \tag {4.82}
$$

which is driven by the total non-translation residual kinetic energy along the direction of the jet or in $z$ -directions. Two hot-spot fluid properties for mode $\ell = 1$ 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/717b6b9938dbadad82510201285257ed2a60d6f8c18462efb826e551ac13d4c5.jpg)



(a)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/e53780144bfaee7c32322cbce863ff4743032734d2eac75828255290a2530a56.jpg)



Figure 4.20: Investigation the fluid properties of mode $\ell = 1$ : (a) examine the small isotropic velocity variance in Eq. (4.85) (b) examine the fraction of the nontranslational hot-spot residual kinetic energy with respect to the total hot-spot residual kinetic energy in Eq. (4.89)


are observed from DEC3D deceleration-phase single-mode simulations, 

$$
M _ {\mathrm {H S}} ^ {\mathrm {3 D}} \sigma_ {3 3} / 2 \simeq \mathrm {K E} _ {\mathrm {H S}} ^ {\mathrm {3 D}} / 3, \tag {4.83}
$$

$$
\mathrm {R K E} _ {\mathrm {H S}} \simeq \mathrm {R K E} _ {\mathrm {S H}}. \tag {4.84}
$$

where the normalized hot-spot residual kinetic energy RKEHS and the normalized shell, including both shocked and un-shocked parts, residual kinetic energy at stagnation are defined as $\mathrm { R K E _ { H S } = ( K E _ { H S } ^ { 3 D } - K E _ { H S } ^ { 1 D } ) _ { s t a g } / K E _ { m a x } ^ { 1 D } }$ and RKESH = $( \mathrm { K E _ { S H } ^ { 3 D } - K E _ { S H } ^ { 1 D } ) _ { s t a g } / K E _ { m a x } ^ { 1 D } }$ respectively. 

The first property is caused by the unique flow structure of mode $\ell = 1$ that satisfies, 

$$
\langle \triangle v _ {3} ^ {2} \rangle \simeq \langle v _ {3} ^ {2} \rangle / 3. \tag {4.85}
$$

Equation (4.85) provides a strong trend for both modest and large perturbations as shown in Fig. (4.20)-(a). The second property, however, is a weaker correlation that is observed to be valid only for modest perturbations, and implies that 

hydrodynamics of modest mode $\ell = 1$ perturbations is, 

$$
\mathrm {K E} _ {\mathrm {H S}} ^ {\mathrm {3 D}} / \mathrm {K E} _ {\max } ^ {\mathrm {1 D}} \simeq \mathrm {R K E} _ {\mathrm {t o t}} / 2, \tag {4.86}
$$

where the normalized total residual kinetic energy at stagnation is defined as $\mathrm { R K E } _ { \mathrm { t o t } } = \mathrm { R K E } _ { \mathrm { H S } } + \mathrm { R K E } _ { \mathrm { S H } }$ . The simplification for $f _ { \mathrm { r k e } } ^ { \mathrm { t o t a l } }$ has the same form of $f _ { \mathrm { r k e } }$ given by Eq. (4.40) by replacing the isotropic residual kinetic energy KEnontranshs,DT,i in the $z$ -direction with the total $z$ -directional residual kinetic energy $M _ { \mathrm { H S } } ^ { \mathrm { 3 D } } \sigma _ { 3 3 } / 2$ . The latter is expressed in terms of the one-third of the total hot-spot residual kinetic energy according to Eq. (4.83), 

$$
f _ {\mathrm {r k e}} ^ {\text {t o t a l}} \simeq 4 \mathrm {K E} _ {\mathrm {H S}} ^ {\mathrm {3 D}} / \mathrm {I E} _ {\mathrm {H S}} ^ {\mathrm {3 D}}. \tag {4.87}
$$

The ratio of 3-D hot-spot kinetic energies to 3-D hot-spot internal energies in Eq. (4.87) is rewritten as, 

$$
\mathrm {K E} _ {\mathrm {H S}} ^ {\mathrm {3 D}} / \mathrm {I E} _ {\mathrm {H S}} ^ {\mathrm {3 D}} = \left(\frac {\mathrm {K E} _ {\mathrm {H S}} ^ {\mathrm {3 D}}}{\mathrm {K E} _ {\mathrm {m a x}} ^ {\mathrm {1 D}}}\right) \left(\frac {\mathrm {K E} _ {\mathrm {m a x}} ^ {\mathrm {1 D}}}{\mathrm {I E} _ {\mathrm {H S}} ^ {\mathrm {1 D}}}\right) \left(\frac {\mathrm {I E} _ {\mathrm {H S}} ^ {\mathrm {1 D}}}{\mathrm {I E} _ {\mathrm {H S}} ^ {\mathrm {3 D}}}\right). \tag {4.88}
$$

The three different terms in Eq. (4.88) are given by: (1) IE3DHS/IE1DHS = 1 − RKE $\mathrm { \ t o t { \Omega } }$ which is an exact result from conservation of total energy at stagnation [WBS+18a], (2) IE1DHS/KE1Dmax ' 1/2 which is observed from the 1-D implosion database meaning that about a half of the maximum 1-D shell kinetic energy is converted into the 1-D hot-spot internal energy at stagnation, and (3) the last term of KE3DHS/KE1Dmax is replaced by Eq. (4.86). Equation (4.88) is rewritten as 

$$
\mathrm {K E} _ {\mathrm {H S}} ^ {\mathrm {3 D}} / \mathrm {I E} _ {\mathrm {H S}} ^ {\mathrm {3 D}} \simeq \mathrm {R K E} _ {\mathrm {t o t}} / (1 - \mathrm {R K E} _ {\mathrm {t o t}}). \tag {4.89}
$$

Equation (4.89) is studied in Fig. (4.20)-(b), which is shown to be valid only for modest mode $\ell = 1$ perturbations defined by condition $R _ { \mathrm { m i n } } ^ { \mathrm { 3 D } } / R _ { \mathrm { 1 D } } ^ { \mathrm { s t a g } } > 1 / 2$ sta . In 

high-performance ICF experiments, mode $\ell = 1$ ion-temperature measurement asymmetries generally belong to the class of modest perturbation regime. Therefore, the ratio of maximum to minimum DT neutron-inferred ion temperatures for mode $\ell = 1$ in Eq. (4.80) is a unique function of the total residual kinetic energy, 

$$
\left(T _ {\mathrm {i , D T}} ^ {\max } / T _ {\mathrm {i , D T}} ^ {\min }\right) _ {\ell = 1} \simeq 1 + 4 \mathrm {R K E} _ {\mathrm {t o t}} / (1 - \mathrm {R K E} _ {\mathrm {t o t}}), \tag {4.90}
$$

which is inverted to give, 

$$
\mathrm {R K E} _ {\mathrm {t o t}} = \xi / (1 + \xi), \tag {4.91}
$$

where the ion-temperature measurement asymmetry parameter $\xi = ( 1 / 4 ) ( R _ { T } -$ $1 ) \geq 0$ is a function of ion-temperature ratio $R _ { T } = T _ { \mathrm { i , D T } } ^ { \mathrm { m a x } } / T _ { \mathrm { i , D T } } ^ { \mathrm { m u } }$ . Therefore, the yield degradation through Eq. (4.78) is a function of neutron-inferred ion-temperature measurement asymmetry parameter for mode $\ell = 1$ , 

$$
\mathrm {Y O C} = [ 1 - \xi / (1 + \xi) ] ^ {\mu}. \tag {4.92}
$$

When the ion temperature asymmetry parameter is small $0 \leq \xi \ll 1$ for modest mode $\ell = 1$ perturbations, the total residual kinetic energy is proportional to the ion temperature ratio in Eq. (4.91) $\mathrm { R K E } _ { \mathrm { t o t } } \simeq \xi = ( 1 / 4 ) ( R _ { T } - 1 )$ and the yield degradation through Eq. (4.78), 

$$
\mathrm {Y O C} = \left[ 1 - \frac {1}{4} \left(\frac {T _ {\mathrm {i , D T}} ^ {\mathrm {m a x}}}{T _ {\mathrm {i , D T}} ^ {\mathrm {m i n}}} - 1\right) \right] ^ {\mu}. (4. 9 3)
$$

Equations (4.90) and (4.93) are only valid for mode $\ell = 1$ in the limit of vanishing isotropic velocity variance in Eq. (4.80) and applying two approximate fluid properties in Eqs. (4.85) and (4.86) in the regime of modest perturbations. For large- $\ell$ modes or fully turbulent hot spot, the anisotropic velocity variance decrease significantly and transit into the isotropy limit: $T _ { \mathrm { m a x } } ^ { \mathrm { B r y s k } } / T _ { \mathrm { m i n } } ^ { \mathrm { B r y s k } }  1$ , which 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/128d3235c7250631adfe9da1f41f7a0e8a3ad427479c15e41e4a8972be22cc63.jpg)



Figure 4.21: Examination of Eq. (4.92) where $\xi \ = \ \textstyle { \frac { 1 } { 4 } } ( R _ { T } - 1 )$ and $\boldsymbol { R _ { T } } ~ =$ poles while T IRIS3Dmax /T IRIS3Dmin . For mode $T _ { \mathrm { m i n } } ^ { \mathrm { l R l S 3 D } }$ max  is observed at the equator. $\ell = 1$ , $T _ { \mathrm { m a x } } ^ { \mathrm { l R l S 3 D } }$ 4   are obsvered at the north and south


have characteristics of small ion-temperature measurement variations among LOS and large minimum inferred ion temperatures. Apart from mode $\ell \ = \ 1$ that has well-behave fluid properties governed by Eqs. (4.85) and (4.86), neutroninferred ion-temperature measurement asymmetries for mid and high modes do not exhibit strong correlations with total residual kinetic energies, resulting in weak correlations between YOC and ion-temperature ratios. 

To validate Eq. (4.93), DEC3D hydrodynamic data at stagnation are postprocessed by IRIS3D [WRF18], a Monte-Carlo based neutron transport code.16 detectors are set up from north to south poles uniformly at a fixed azimuthal angle $\phi = 0$ to infer ion temperatures from the width of the Gaussian-fitted neutron energy spectra. Figure (4.21) compares the yield degradation with the ratio of neutron-inferred maximum to minimum ion temperatures obtained from IRIS3D. The data of mode $\ell = 1$ curve is shown accurately bounded in between YOC = $[ 1 - \xi / ( 1 + \xi ) ] ^ { 4 . 4 }$ and Y $\mathrm { \Delta } \mathrm { T } = [ 1 - \xi / ( 1 + \xi ) ] ^ { 5 . 5 }$ curves as predicted by Eq. (4.92). IRIS3D shows that mode $\ell = 2$ exhibits a small ion-temperature ratio due to large isotropic velocity variances or large $T _ { \mathrm { m i n } }$ . Ion temperature ratios for 2-D mode $\ell = 4 , m = 0$ are shown significantly larger than 3-D mode $\ell = 4 , m = 2$ because of 3-D spherical appearance of spike-to-bubble flow structure. Mid and high modes $\ell = 5 - 1 2$ do not exhibit large ion temperature ratios because their clean volumes have negligible small hot-spot residual kinetic energy. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/d23f84da7af178da06e90a2fce4ddb3d4fde1e9ffed02500f5b8a8d1fd1b1941.jpg)



Figure 4.22: A correlation study between OMEGA experimental yields and DD/DT neutron-inferred ion-temperature measurement asymmetries. The impact of hot-spot flow isotropy represented by the term $T _ { \mathrm { m i n } } ^ { \mathrm { L D D } } / T _ { \mathrm { m i n } } ^ { \mathrm { L D } }$ is a function of non-translational hot-spot residual kinetic energies given by Eq. (4.39). The impact of hot-spot flow anisotropy represented by the term $1 - \xi$ is a function of the shell residual kinetic energies for low modes. The impativities between DD and DT is represented by the term $G _ { 1 0 } ^ { \mathrm { L I L A C } }$ is the mode 10 growth factor obtained from 1-D $1 - T _ { \mathrm { D D } } ^ { \mathrm { L O S ( 1 3 . 4 m ) } } / T _ { \mathrm { D T } } ^ { \mathrm { L O S ( p e t a l ) } }$ $L I L A C$ DT simulations. -.


Figure (4.22) shows a correlation study to summarize impacts of DD and DT ion-temperature measurement asymmetries on OMEGA experimental yields. When hot-spot flow anisotropies are strong, DT neutron-inferred ion-temperature exhibits strong directional-dependent measurements, captured by the large exponent of the hot-spot flow anisotropy term $1 - \xi$ , where $\xi = ( T _ { \mathrm { m i n } } ^ { \mathrm { L D D } } / T _ { \mathrm { m i n } } ^ { \mathrm { L D T } } - 1 ) / 4$ 1 . Impact of hot-spot flow isotropy is characterized by the ratio of DD to DT minimum inferred ion-temperatures T DDmin/T DTmin, $T _ { \mathrm { m i n } } ^ { \mathrm { L D D } } / T _ { \mathrm { m i n } } ^ { \mathrm { L D } }$ min/ which is a function of non-translational hotspot residual kinetic energies by Eq. (4.39). Although hot-spot flow isotropies are invisible in DD or DT ion-temperature measurements, they quantify the content of non-translational hot-spot residual kinetic energies. The minimum DD ion temperature is extrapolated from Eq. (4.72). Secondary effects of different fusion reactivities between DD and DT is represented by the term $1 - T _ { \mathrm { D D } } ^ { \mathrm { L O S ( 1 3 . 4 m ) } } / T _ { \mathrm { D T } } ^ { \mathrm { L O S ( p e t a l ) } }$ , and $G _ { 1 0 } ^ { \mathrm { L I L A C } }$ is the mode 10 growth factor obtained from 1-D $L I L A C$ simulations to present the dominant laser mode driven by 60-beam geometry in OMEGA. 

# 4.8 Conclusion

A comprehensive analysis of the impact of isotropic velocity variance on the DD/DT neutron-inferred ion-temperature ratio is presented. A general expression of isotropic velocity variance for single modes, multi modes and fully-developed turbulence is derived. The secondary effects of different fusion reactivities leading to differences in DD and DT neutron-averaged velocity variance are studied. The regime of large multi-mode perturbations are shown to produce large enough isotropic hot-spot flows driving the ratio $T _ { \mathrm { m i n } } ^ { \mathrm { L D D } } / T _ { \mathrm { m i n } } ^ { \mathrm { L D } }$ to the theoretical limit of 0.8. Three approximate solution strategies are proposed to diagnose fusion ion thermal temperature, hot-spot flow isotropy and flow anisotropy through utilizing the measurements of the DD/DT ion-temperature ratio. The six-LOS technique shows a promising capability to extrapolate the true minimum thermal ion temperatures, and its prediction accuracy is shown limited only by error propagation. An expression for the minimum DD ion temperature is derived by removing the anisotropic velocity variance, and is shown strongly correlated with experimental yields in OMEGA data. An analytical expression to quantify the effect of mode $\ell = 1$ ion-temperature measurement asymmetry on yield degradation in the limit of strong hot-spot flow anisotropies is derived. In multi-mode simulation, the slope for DD/DT ion-temperature ratio caused by LOS variations is shown to reproduce the trend of OMEGA experimental data. 

# 5 DEC3D Computer Code

# 5.1 Motivation

• DEC3D is a sopyhciated extension of the original code DEC2D, written by K. Anderson, R. Betti and T. A. Gardiner.[ABG01]. 

• DEC3D is a code development project with the objective to model the full three-dimensional physics in the deceleration phase of inertial confinement fusion implosions using advanced modern numerical techniques. 

In the 2-D code development, three important new features were implemented in this thesis: (1) the upgrade of hydrodynamic solver from MacCormack finitedifference scheme [Sod78] to Godunov’s scheme with MUSCL [vL79] and PPM [CW84] approximate Riemann solvers to attain a strong shock capturing capability in simulations for highly nonlinear RT instabilities, (2) the implementation of the multi-group radiation and alpha particle transport packages for accurate accounting radiation coolings and radiative ablative RT instabilities, as well as strong alpha heatings in the burning plasma regime, and (3) the upgrade of Cartesian mesh into spherical mesh, with the application of the macro-zoning technique, to attain a noise-free simulation environment due to the grid effects in spherical implosions. 

In the 3-D code development, because of the parallel architecture, most of numerical methods such as the direct solve by the Gaussian elimination method for implicit diffusions in the original DEC2D can not be applied in the 3-D parallel 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/27a48a04fd11168d27b9619927c8896918329ecadd67a323164877d4b0a621d0.jpg)



TC14047J1



Figure 5.1: The 3-D electron temperature contour at 1 keV at stagnation for the single-mode $Y _ { \ell = 1 2 } ^ { m = 6 }$ . The resolution is $1 2 8 \times 1 2 8 \times 2 5 6$ for $r , \theta , \phi$ zones, respectively.


code after the domain decomposition. The following is a summary of the most up-to-date code development status for $D E C \mathcal { Z } D$ and DEC3D accomplished in this thesis. 


Table 5.1: The summary for DEC2D and DEC3D.


<table><tr><td>Numerical methods</td><td>Original DEC2D</td><td>DEC2D</td><td>DEC3D</td></tr><tr><td rowspan="2">Mesh</td><td>Cartesian</td><td>Cartesian</td><td>Cartesian</td></tr><tr><td>-</td><td>Spherical</td><td>Spherical</td></tr><tr><td>Moving-mesh</td><td>Finite-difference</td><td>Conservative</td><td>Conservative</td></tr><tr><td>Hydro solver</td><td>MacCormak</td><td>MUSCL</td><td>PPM</td></tr><tr><td>Diffusion solver</td><td>Direct solve</td><td>Direct solve</td><td>Iterative solve bySOR &amp; HYPRE</td></tr><tr><td>Multi-group radiation</td><td>No</td><td>Yes</td><td>Yes</td></tr><tr><td>Multi-group alpha</td><td>No</td><td>Yes</td><td>Yes</td></tr></table>

In this chapter, a comprehensive description for literature review on theories, advanced modern numerical techniques, benchmark tests, and code description are presented. 

The integrated performance for $D E C 3 D$ is summarized as follows. Figure (5.1) shows the perturbed simulation of a single-mode $Y _ { \ell = 1 2 } ^ { m = 6 }$ . For perturbed simulations, 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/cfb75da57e4f698c1509c535bf83b1604630b8d6e1b1c5835425c8bd517559f3.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/779f138f69e1e050e6ed0b3857f5c18aaa531cd7d80b15f880096c78b8df9f31.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/7db2aded4d3e7cacc5b82d4959b3820fb159aebeace9aa7e8a6529cbb61c7302.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/081eded526979e1fd4624f7edba93f0857a28bcae2fa8142aa3dbb1ebdbedb82.jpg)



Figure 5.2: The benchmark test for a full 3-D, deceleration-phase, clean simulation for the shot 77068. Comparison of solutions between $L I L A C$ (black circles) and DEC3D (solid red curves) for (a) density, (b) total pressure, (c) radial velocity, and (d) electron temperature profiles at stagnation using resolution $1 2 8 \times 6 4 \times 1 2 8$ in $r , \theta , \phi$ directions. The 1-D burn radius at stagnation for shot 77068 is ∼ 19 µm, and the return shock is located at ∼ 28 µm.


neutron yields converge when the angular resolution is about $1 2 8 \times 2 5 6$ in $\theta$ and $\phi$ directions for mode $\ell = 1 2$ . The resolution of the DEC3D single-mode database is $1 2 8 \times 1 2 8 \times 2 5 6$ in $r , \theta , \phi$ directions, about 25 zones per wavelength for mode 10, which is sufficient for Legendre modes $\ell = 1$ to 12 studies. 

The benchmark test for a full 3-D, deceleration-phase, clean 1-D unperturbed simulation for the shot 77068 is shown in Figs. (5.2)–(5.3) to summarize the integrated numerical performance of $D E C 3 D$ . In the radial velocity profile, the steep spatial gradient across the return shock at $r \sim 2 8$ µm is well resolved by the third-order PPM method, while $L I L A C$ solution is shown to be more diffusive due to the use of numerical viscosity. The performance of HYPRE for the thermal diffusion are validated in the electron temperature profile. A slight larger total pressure at $r ~ = ~ 0$ in $D E C 3 D$ than $L I L A C$ is observed, which results in slight increases in the ion temperature and neutron production. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/3c9f89cea3d3ad381ccea688a8e0f071a3bf7b591b2912f005ef879e49bccb0e.jpg)



Figure 5.3: The benchmark test for neutron productions without alpha heating. Comparison of solutions between $L I L A C$ (black circles) and $D E C 3 D$ (solid red curves) for the temporal history of neutron rate. The fusion reactivity $< \sigma v >$ scales with ion temperatures $T _ { \mathrm { i } }$ in a power law $\sim T _ { \mathrm { i } } ^ { 3 . 8 5 }$ for the temperature range $0 . 2 < T _ { \mathrm { i } } ^ { \prime } < 5$ keV in the BUCKY [HMS05] fusion reactivity model.


# 5.2 Physical models

# 5.2.1 Governing equations

DEC3D is a deceleration phase ICF code that solves an inviscid single-fluid twotemperature plasma model. The single-fluid pressure $P = P _ { \mathrm { e } } + P _ { \mathrm { i } }$ is a sum of the electron $P _ { \mathrm { { e } } }$ and ion $P _ { \mathrm { i } }$ pressures. Ideal gas equation of state and fully ionized plasma are assumed i.e., with the ratio of specific heats $\gamma = 5 / 3$ , the averaged Deuterium-Tritium (DT) ion charge $Z = 1$ and the ideal gas laws for electrons and ions $P _ { \mathrm { e / i } } = n _ { \mathrm { e / i } } T _ { \mathrm { e / i } }$ , where $n _ { \mathrm { e / i } }$ and $T _ { \mathrm { e / i } }$ are the number density and thermal temperature for electrons and ions respectively. During the deceleration and disassembly phases, the high temperature hot spot is sufficient to maintain the fully ionized plasma state, and the high density shell remains weakly degenerate $P / P _ { \mathrm { F e r m i } } \sim 4 - 5$ so that the ideal gas approximation is adequate. Here $P _ { \mathrm { F e r m i } }$ is the Fermi pressure [AtV04]. 

The governing equations for the inviscid single-fluid plasma are as follows. 

$$
\partial_ {t} \rho + \vec {\nabla} \cdot (\rho \vec {v}) = 0, \tag {5.1}
$$

$$
\partial_ {t} (\rho \vec {v}) + \vec {\nabla} \cdot (\rho \vec {v} \otimes \vec {v} + \hat {\mathbb {I}} P) = 0, \tag {5.2}
$$

$$
\partial_ {t} \left(\frac {P}{\gamma - 1} + \frac {1}{2} \rho v ^ {2}\right) + \vec {\nabla} \cdot \left[ \vec {v} \left(\frac {\gamma P}{\gamma - 1} + \frac {1}{2} \rho v ^ {2}\right) + \vec {Q} \right] = S. \tag {5.3}
$$

Here $\hat { \mathbb { I } }$ is the identity matrix, $\otimes$ is the outer product, $\vec { Q }$ and $S$ represent the sum of heat flux and heat source terms caused by thermal diffusion, alpha heating and radiation transport, $\vec { v }$ is the single-fluid velocity and $\rho = n _ { \mathrm { i } } m _ { \mathrm { i } }$ is the single-fluid mass density, where $\bar { m } _ { \mathrm { i } } = ( m _ { \mathrm { D } } + m _ { \mathrm { T } } ) / 2$ is the average of D and T ion masses. The internal energy of electrons change dramatically in space and time as a result of the strong electron heat conduction that diffuses the heat of hot spot into the cold shell, the strong electron drag force that slows down alpha particles in the straight-line motion, the intense radiation-material interaction that cools down the hot-spot electrons by self-emissions of high-energy group photons or radiative heat diffusion by re-absorbing low-energy group photons, and the strong heat transfer into ions by electron-ion collisions. These rigorous plasma energy transfer processes leads to different electron and ion thermodynamic temperatures, resulting in the single-fluid two-temperature plasma model. 

Figure (5.4) shows the flow chart for executing different physical modules in one time-step by applying the time-splitting technique [Har11]. The coarsening mapping $\hat { \mathcal { M } }$ is computed first to map the fine states of fluid variables $\vec { Q }$ onto the coarser states $\vec { Q } _ { \mathrm { c } } = \hat { \mathcal { M } } \vec { Q }$ , followed by the hydrodynamic update: $\partial _ { t } \vec { Q } _ { \mathrm { c } } = \hat { \mathcal { H } } \vec { Q } _ { \mathrm { c } }$ , which operates on coarser states to attain a relaxed time-step size $\triangle t$ . Next various plasma transport phenomena are followed, which operate on the fine states $\vec { Q }$ and are solved implicitly using the relaxed time-step size $\triangle t$ obtained from the hydrostep, including the heat transfer due to the electron and ion thermal diffusions: $\partial _ { t } \vec { Q } = \hat { T } \vec { Q }$ , the electron and ion equilibration: $\partial _ { t } \vec { Q } = \hat { \varepsilon } \vec { Q }$ , the radiation-material 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/353a0a228c89fa424705e733cb641d069d1b3c988bd76b25b528e4f22c41d165.jpg)



Figure 5.4: The flow chart for executing different physical module in one time-step in DEC3D.


interaction: $\partial _ { t } \vec { Q } = \hat { \mathcal { R } } \vec { Q }$ , and the straight-line slowing-down of alpha particles due to electron drags: $\partial _ { t } \vec { Q } = \hat { A } \vec { Q }$ . 

In the operator representation, the system of Eqs. (5.1)–(5.3) 

$$
\partial_ {t} \vec {Q} = \left[ \hat {\mathcal {H}} + \hat {\mathcal {T}} + \hat {\mathcal {E}} + \hat {\mathcal {R}} + \hat {\mathcal {A}} \right] \vec {Q}, \tag {5.4}
$$

has an exact solution within the time-step size $\triangle t$ , 

$$
\begin{array}{l} \vec {Q} (t ^ {n + 1}) = e ^ {\left[ \hat {\mathcal {H}} + \hat {\mathcal {T}} + \hat {\mathcal {E}} + \hat {\mathcal {R}} + \hat {\mathcal {A}} \right] \triangle t} \vec {Q} (t ^ {n}), \\ = \underbrace {e ^ {\hat {\mathcal {H}} \triangle t} e ^ {\hat {\mathcal {T}} \triangle t} e ^ {\hat {\mathcal {E}} \triangle t} e ^ {\hat {\mathcal {R}} \triangle t} e ^ {\hat {\mathcal {A}} \triangle t} \vec {Q} \left(t ^ {n}\right)} _ {\text {L i e s p l i t t i n g} \vec {Q} _ {\text {L i e}} ^ {n + 1}} + \mathcal {O} (\triangle t ^ {2}) \tag {5.5} \\ \end{array}
$$

The solution form of $\vec { Q } _ { \mathrm { L i e } } ^ { n + 1 }$ in the second line in Eq. (5.5) is called Lie splitting, which is obtained by an operator expansion using Baker-Campbell- Hausdorff (BCH) formula [Sha94] and has a leading error proportional to $\triangle t ^ { 2 }$ . For example, an advection-diffusion process $\partial _ { t } \vec { Q } \ = \ \left[ \hat { \mathcal { H } } + \hat { \mathcal { T } } \right] \vec { Q }$ , which is only governed by a hydro $\hat { \mathcal { H } }$ and a thermal diffusion $\hat { \tau }$ steps, has the following approximated solution 

using the BCH formula, 

$$
\begin{array}{l} \vec {Q} (t ^ {n + 1}) = e ^ {[ \hat {\mathcal {H}} + \hat {\mathcal {T}} ] \triangle t} \vec {Q} (t ^ {n}), \\ = e ^ {\hat {\mathcal {H}} \triangle t} e ^ {\hat {\mathcal {T}} \triangle t} \vec {Q} (t ^ {n}) - \frac {\triangle t ^ {2}}{2} \left[ \hat {\mathcal {H}}, \hat {\mathcal {T}} \right] \vec {Q} (t ^ {n}). \tag {5.6} \\ \end{array}
$$

Unless the operators of $\hat { \mathcal { H } }$ , $\hat { \tau }$ , $\hat { \mathcal { E } }$ , $\hat { \mathcal { R } }$ and $\hat { A }$ commute with each other 

$$
\left[ \hat {\mathcal {H}}, \hat {\mathcal {T}} \right] = \left[ \hat {\mathcal {T}}, \hat {\mathcal {E}} \right] = \left[ \hat {\mathcal {E}}, \hat {\mathcal {R}} \right] = \left[ \hat {\mathcal {R}}, \hat {\mathcal {A}} \right] = 0, \tag {5.7}
$$

the Lie splitting is first-order accurate in time. However, the commutation relations in Eq. (5.7), requiring all PDEs for hydrodynamics and plasma heat transfer processes being linear operators in space, are not easily fulfilled in reality. 

In this chapter, five individual sections are presented in the following order to cover the governing equations in $D E C \mathcal { Z } D$ and $D E C 3 D$ for 

(1) the hydrodynamics module ( $\partial _ { t } \vec { Q } = \hat { \mathcal { H } } \vec { Q }$ ), 

$$
\partial_ {t} \rho + \vec {\nabla} \cdot (\rho \vec {v}) = 0, \tag {5.8}
$$

$$
\partial_ {t} (\rho \vec {v}) + \vec {\nabla} \cdot (\rho \vec {v} \otimes \vec {v} + \hat {\mathbb {I}} P) = 0, \tag {5.9}
$$

$$
\partial_ {t} \left(\frac {P}{\gamma - 1} + \frac {1}{2} \rho v ^ {2}\right) + \vec {\nabla} \cdot \left[ \vec {v} \left(\frac {\gamma P}{\gamma - 1} + \frac {1}{2} \rho v ^ {2}\right) \right] = 0, \tag {5.10}
$$

(2) the thermal diffusion module ( $\partial _ { t } \vec { Q } \ = \ \hat { \tau } \vec { Q }$ ) with electron and ion thermal conductivities $\kappa _ { \mathrm { e } }$ and $\kappa _ { \mathrm { i } }$ , 

$$
\partial_ {t} \left(\frac {P _ {\mathrm {e}}}{\gamma - 1}\right) = \vec {\nabla} \cdot \kappa_ {\mathrm {e}} \vec {\nabla} T _ {\mathrm {e}}, \tag {5.11}
$$

$$
\partial_ {t} \left(\frac {P _ {\mathrm {i}}}{\gamma - 1}\right) = \vec {\nabla} \cdot \kappa_ {\mathrm {i}} \vec {\nabla} T _ {\mathrm {i}}, \tag {5.12}
$$

(3) the electron-and-ion equilibration module ( $\partial _ { t } \vec { Q } = \hat { \varepsilon } \vec { Q }$ ) with the electron and 

ion relaxation time $\tau _ { \mathrm { e i } }$ , 

$$
\partial_ {t} T _ {\mathrm {e}} = - \frac {1}{\tau_ {\mathrm {e i}}} \left(T _ {\mathrm {e}} - T _ {\mathrm {i}}\right), \tag {5.13}
$$

$$
\partial_ {t} T _ {\mathrm {i}} = - \frac {1}{\tau_ {\mathrm {e i}}} \left(T _ {\mathrm {i}} - T _ {\mathrm {e}}\right), \tag {5.14}
$$

(4) the multi-group radiation transport module ( $\partial _ { t } \vec { Q } \ = \ \hat { \mathcal { R } } \vec { Q }$ ) with Planck and Rosseland opacities $\kappa _ { g } ^ { \mathrm { P } }$ and $\kappa _ { g } ^ { \mathrm { R } }$ , and the flux-limited diffusion coefficient $D _ { g }$ , the group-weighted self-emission factor, radiation energy density and radiation pressure $B _ { g }$ , $\langle U _ { g } \rangle$ and $\langle P _ { g } \rangle$ , 

$$
\partial_ {t} \langle U _ {g} \rangle + \vec {\nabla} \cdot \vec {v} \langle U _ {g} \rangle + \langle P _ {g} \rangle \vec {\nabla} \cdot \vec {v} = \vec {\nabla} \cdot \bar {D} _ {g} \left(\kappa_ {g} ^ {\mathrm {R}}\right) \vec {\nabla} \langle U _ {g} \rangle + c \kappa_ {g} ^ {\mathrm {P}} \left(B _ {g} - \langle U _ {g} \rangle\right), \tag {5.15}
$$

(5) the one-group or multi-group alpha particle transport module ( $\partial _ { t } \vec { Q } = \hat { \mathcal { A } } \vec { Q }$ ) with the alpha particle energy density, diffusion coefficient and birth energy $\varepsilon _ { \alpha }$ , $D _ { \alpha }$ and $E _ { \alpha 0 } = 3 . 5$ -MeV, D and T ion number densities $n _ { \mathrm { D } }$ , $n _ { \mathrm { T } }$ and fusion reactivity $< \sigma v > _ { \mathrm { D T } }$ , and the alpha-electron relaxation time $\tau _ { \alpha \mathrm { e } }$ , 

$$
\partial_ {t} \varepsilon_ {\alpha} = \vec {\nabla} \cdot D _ {\alpha} \vec {\nabla} \varepsilon_ {\alpha} + n _ {\mathrm {D}} n _ {T} <   \sigma v > _ {\mathrm {D T}} E _ {\alpha 0} - \frac {\varepsilon_ {\alpha}}{\tau_ {\alpha \mathrm {e}}}. \tag {5.16}
$$

Since the fully implicit discretization, which is unconditionally numerical stable for large time-step sizes, is applied to solve for all diffusion type equations, the relaxed time-step size $\triangle t$ obtained from the coarser states in the hydro-step is used in modules (2)–(5). 

# 5.2.2 Cartesian & spherical mesh discretization

# Cartesian mesh DEC2D and DEC3D

In Cartesian mesh versions, the cell-center coordinates $( x _ { i } , y _ { j } , z _ { k } )$ are defined as $x _ { i } = ( i - \alpha ) \triangle x _ { x }$ , $y _ { j } = ( j - \alpha ) \triangle { y _ { y } }$ and $z _ { k } = ( k - \alpha ) \triangle z _ { z }$ . The value of $\alpha$ can be 

taken as 1 or 0.5 , but it is restricted to be 0.5 in spherical coordinates to avoid singularities. A uniform discretization $\triangle x _ { x } = \triangle y _ { y } = \triangle z _ { z } = L ( t ) / N$ in $x$ , $y$ and $z$ directions is adopted to provide a uniform resolution for RT instabilities, where $L ( t )$ is the length of the simulation domain in one direction at the time $t$ and $N$ is the number of cells. At each time step $\triangle t$ , the length of simulation domain $L ( t )$ is updated explicitly according to $\begin{array} { r } { \frac { d L } { d t } = \beta v _ { \mathrm { C M } } ^ { \mathrm { s h e l l } } ( t ) } \end{array}$ , 

$$
\frac {L \left(t ^ {n + 1}\right) - L \left(t ^ {n}\right)}{\triangle t} = \beta v _ {\mathrm {C M}} ^ {\text {s h e l l}} \left(t ^ {n}\right), \tag {5.17}
$$

where $v _ { \mathrm { C M } } ^ { \mathrm { s h e l l } } ( t )$ is the center-of-mass velocity of the imploding shell in the deceleration phase or the exploding shell in the disassembly phase, and $\beta$ is a constant in between $1 - 2$ to control the moving mesh velocity. 

For 3-D simulations, the problem of increasing computational times and CPU memories due to operating with a large data size to store 3-D variables, leads to the implementation of massively parallel simulations. The domain decomposition for the Cartesian-mesh version DEC3D is implemented as follows. Let $m _ { x }$ , $m _ { y }$ and $m _ { z }$ be the number of sub-domains in each orthogonal direction respectively, the total number of sub-domains is therefore given by $m _ { x } \times m _ { y } \times m _ { z }$ , whereas the total number of cells in each sub-domain is given by $N _ { x } ^ { \mathrm { s u b } } \times N _ { y } ^ { \mathrm { s u b } } \times N _ { z } ^ { \mathrm { s u b } }$ , where $N _ { x } ^ { \mathrm { s u b } } ~ = ~ N / m _ { x }$ , $N _ { y } ^ { \mathrm { s u b } } ~ = ~ N / m _ { y }$ and $N _ { z } ^ { \mathrm { s u b } } ~ = ~ N / m _ { z }$ are respectively, the number of cells in each sub-domain along $x$ , $y$ and $z$ directions. Under the domain decomposition, the storage for a 3-D variable $N ^ { 3 }$ is split into $N _ { x } ^ { \mathrm { s u b } } \times N _ { y } ^ { \mathrm { s u b } } \times N _ { z } ^ { \mathrm { s u b } }$ stored by each core, which is labelled as “rank”. The reduced size to store a 3- D array within a sub-domain not only relaxes the memory requirement but also shorten the computational times to complete a $N _ { x } ^ { \mathrm { s u b } } \times N _ { y } ^ { \mathrm { s u b } } \times N _ { z } ^ { \mathrm { s u b } }$ do-loop. 

Figure (5.5) demonstrates the flexibility of Cartesian topology domain decomposition implemented in the Cartesian-mesh version DEC3D through the message-passing-interface (MPI). Two figures on the top show a sector of 3-D 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/4909c45c0311df6a5ad2bfcd1df24aa63488ea417a998316445e5918c7a5a339.jpg)



Figure 5.5: The Cartesian topology domain decomposition [CLS94] for the Cartesian mesh version DEC3D [WBB+15] through the message-passing-interface (MPI) for a 3-D mass density profile.


mass density profile with domain decompositions in $x$ , $y$ and $z$ directions defined by $m _ { x } \times m _ { y } \times m _ { z }$ , in which $2 \times 2 \times 2$ is shown on the left and $2 \times 2 \times 1$ on the right, Whereas two figures at the bottom show the full 3-D mass density profile with domain decompositions of $2 \times 2 \times 2$ on the left and $4 \times 4 \times 4$ on the right. 

Within each time-step $\triangle t$ , cell information at the boundaries are exchanged with neighboring sub-domains through MPI send and receive functions, whereas ghost cells are obtained from physical boundary conditions such as the azimuthal rotational symmetry along the poles at inner boundaries or the zero inflow of fluid, heat and radiation from the vacuum at outer boundaries without MPI exchanges. The exchange of boundary only operates on variables on the fine mesh. 

Figure (5.6) shows the MPI send and receive scheme implemented in the Cartesian-mesh version DEC3D. The left graph shows the exchange of cell information at boundaries, whereas the right graph show a $2 \times 2$ domain decomposi-

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/00dc94479b1077c11c972ed3c33bbe6d1fb191018d724ffb1499c95363010e43.jpg)



(a)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/193fc5b657e5f366542c98cdeba60499139ba185bea37d4067a6c83b0d9e88d2.jpg)



(b)



Figure 5.6: Exchanging boundary information between neighboring sub-domains through MPI send and receive functions in one time-step.


tion in order to visualize MPI send and receive operations to complete an explicit update for a 2-D mass density profile in a hydro-step.. 

Code descriptions for the Cartesian-mesh DEC3D were presented in [WBB+15]. The most up-to-date numerical methods and available solvers are summarized in Table 5.2. However, the multi-group radiation transport is solved by the classical red-black SOR iteration, which was observed to convergence extremely slowly for X-ray diffusions. The multi-grid solver by HYPRE was implemented only in the spherical-mesh DEC3D as descried in Table 5.2. Apart from the third-order PPM hydro solver and HYPRE diffusion solver, all up-to-date solvers such as the multigroup alpha transport were available in the Cartesian-mesh DEC2D in this thesis, because the second-order MUSCL hydro solver and the direct solve of implicit diffusions by Gaussian elimination were observed working efficiently. 

In future work, PPM and HYPRE will be implemented into this Cartesianmesh version DEC3D for the interest of high-resolution and strong shock-capturing capabilities for turbulence simulations. 


Table 5.2: The summary for Cartesian-mesh DEC2D and DEC3D.[WBB+15]


<table><tr><td>Numerical methods</td><td>DEC2D</td><td>DEC3D</td></tr><tr><td>Mesh</td><td>Cartesian</td><td>Cartesian</td></tr><tr><td>Moving-mesh</td><td>Conservative</td><td>Conservative</td></tr><tr><td>Hydro solver</td><td>MUSCL</td><td>MUSCL</td></tr><tr><td>Diffusion solver</td><td>Gaussian elimination</td><td>Iterative solve by SOR</td></tr><tr><td>Multi-group radiation</td><td>Yes</td><td>Yes</td></tr><tr><td>Multi-group alpha</td><td>Yes</td><td>No</td></tr><tr><td>One-group alpha</td><td>Yes</td><td>Yes</td></tr></table>

# Spherical mesh DEC2D and DEC3D

In spherical coordinates, the cell-center coordinates $( r _ { i } , \theta _ { j } , \phi _ { k } )$ are defined as $r _ { i } =$ $( i - 1 / 2 ) \triangle r$ , $\theta _ { j } = ( j - 1 / 2 ) \triangle \theta$ , and $\phi _ { k } = ( k - 1 / 2 ) \triangle \phi$ . The discretization in radius and angles are defined by $\triangle r = R ( t ) / N _ { r }$ , $\triangle \theta = \pi / N _ { \theta }$ , and $\bigtriangleup \phi = 2 \pi / N _ { \phi }$ , with indices $i = 1 , . . . , N _ { r }$ , $j = 1 , . . . , N _ { \theta }$ , and $k = 1 , . . . , N _ { \phi }$ . The total number of cells is $N _ { r } \times N _ { \theta } \times N _ { \phi }$ , and the length of the radial domain $R ( t )$ is updated explicitly according to $\begin{array} { r } { \frac { d } { d t } R ( t ) = \beta v _ { \mathrm { C M } } ^ { \mathrm { s h e l l } } ( t ^ { n } ) } \end{array}$ in the same manner as Eq. (5.17 ). 

At the outer radii $r _ { N }$ and $r _ { N + 1 }$ , the zero inflow boundary condition is applied to the mass density $\rho _ { N } ~ = ~ \rho _ { N + 1 }$ , the total pressure $P _ { N } = P _ { N + 1 }$ and the radial velocity $v _ { r }$ , $\mathbf { \nabla } _ { N } = v _ { r }$ , $N { + } 1$ . At the origin and along the poles, no boundary condition is imposed for 3-D perturbing flows. Because the periodicity in the polar angle $\theta$ and the azimuthal angle $\phi$ set up a 3-D relation to map the variables from the ghost cells with that from the interior cells. Equation (5.18) shows the set of 3-D relations to map values of ghost cell $Q _ { g }$ at coordinates $r _ { g }$ , $\theta _ { g }$ , $\theta _ { N _ { \theta } + g }$ , $\phi _ { g }$ and $\phi _ { N _ { \phi } + g }$ , where $g = 1 - 3$ because three ghost cells are needed for PPM third-order hydro solver. Figure (5.7) shows an example of spherical mesh for one radial sub-domain. The boundary information for three ghost cells shown in red circles in each end of 


Radial domain decomposition


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/40dfec91f9c9eecb26baf0c6165eb7e4c5b15e0730ba3e54b957c8a2549fd5c5.jpg)



Figure 5.7: Three ghost cells are required at each boundary of the radial subdomain for the third-order explicit PPM hydro update.


the radial domain are obtained by MPI boundary exchanges. 

$$
Q \left(r _ {1 - g}, \theta_ {j}, \phi_ {k}\right) = Q \left(r _ {g}, \theta_ {N _ {\theta} - j + 1}, \phi_ {k + \operatorname {S i g n} \left(N _ {\phi} / 2 - k\right) \times N _ {\phi} / 2}\right),
$$

$$
Q \left(r _ {i}, \theta_ {1 - g}, \phi_ {k}\right) = Q \left(r _ {i}, \theta_ {g}, \phi_ {k + \operatorname {S i g n} \left(N _ {\phi} / 2 - k\right) \times N _ {\phi} / 2}\right),
$$

$$
Q \left(r _ {i}, \theta_ {N _ {\theta} + g}, \phi_ {k}\right) = Q \left(r _ {i}, \theta_ {N _ {\theta} + 1 - g}, \phi_ {k + \operatorname {S i g n} \left(N _ {\phi} / 2 - k\right) \times N _ {\phi} / 2}\right), \tag {5.18}
$$

$$
Q \left(r _ {i}, \theta_ {j}, \phi_ {1 - g}\right) = Q \left(r _ {i}, \theta_ {j}, \phi_ {N _ {\phi} + 1 - g}\right),
$$

$$
Q \left(r _ {i}, \theta_ {j}, \phi_ {N _ {\phi} + g}\right) = Q \left(r _ {i}, \theta_ {j}, \phi_ {g}\right).
$$

The domain decomposition is applied in the radial direction through messagepassage-interface (MPI). [CLS94] Let $m$ be the total number of sub-domains in the radial direction, the total number of cells in each sub-domain is given by $N _ { r } ^ { \mathrm { s u b } } \times N _ { \theta } \times N _ { \phi }$ , where $N _ { r } ^ { \mathrm { s u b } } = N _ { r } / m$ . Figure (5.8) shows an example of 3-D single-mode $\ell = 6$ simulation. The radial domain is decomposed into six subdomains, showing different layer of mass density profiles for the 3-D RT spikes. The core of rank 0 stores data for the inner most cells, whereas the core of rank 5 store data for the outer most cells. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/b9eff5161b8b7a5a59e9a703eaf2ba9e74c551439114de0ed1827830f4090632.jpg)



Figure 5.8: The radial domain decomposition of $D E C 3 D$ into six sub-domains for a 3-D mass density profile for a single-mode $\ell = 6$ simulation.


# 5.2.3 Macro-zoning

Near the origin $r _ { i } \ \to \ 0$ and along the north pole $\theta _ { j } \ \to \ 0$ and the south pole $\theta _ { j } \to \pi$ , the polar arc length $\triangle S _ { \theta } = r _ { i } \triangle \theta$ and the azimuthal arc length $\triangle S _ { \phi } =$ $r _ { i } \sin \theta _ { j } \triangle \phi$ for finite-volume cells are too small, resulting in an extremely small time-step size $\triangle t$ as required by the Courant condition [CF76]. The problem of small time-step size only affects the explicit schemes for hyperbolic equations such as PPM hydrodynamics, because the distance travelled by a characteristic wave in an explicit update for the wave propagation cannot be lager than the cell size. The macro-zoning technique is applied in DEC3D to map fluid variables from the fine mesh onto a coarser mesh to relaxed the time-step size. 

$$
\Delta t = \operatorname {M i n} \left[ \triangle r / v _ {\max }, \triangle S _ {\theta} / v _ {\max }, \triangle S _ {\phi} / v _ {\max } \right] _ {i, j, k}, \tag {5.19}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/39d57450f373c04fae162112707bc45a6bcaed449ea6c84e493b5c192b63a15c.jpg)



Figure 5.9: Example of coarse mesh generation of azimuthal mesh along the poles in DEC3D.


where $v _ { \mathrm { m a x } } = \mathrm { M a x } [ v _ { r } \pm c _ { s } , v _ { \theta } \pm c _ { s } , v _ { \phi } \pm c _ { s } ]$ is the local maximum speed of hydrodynamic signal propagation and $c _ { s } = \sqrt { \gamma P / \rho }$ is the local sound speed. In each time-step, the coarsening algorithm is computed within each radial layer of cells: (I) first the fine mesh detection, (II) second the coarse mesh generation (III) and finally the prolongation and restriction [TOS01] of primitive variables between fine and coarse mesh. Figure (5.9) shows an example of coarse mesh structure near the $z$ -axis obtained by the three-step coarsening algorithm. The minimum arc length in the polar and azimuthal angles on the coarser mesh are restricted to be about $\triangle r / 2$ , so that the final time-step size is determined by the radial discretization $\triangle r$ in Eq. (5.19) 

In step (I), the fine meshes are detected whenever either one of the following fine mesh definitions are satisfied for an individual finite-volume cell. 

$$
\triangle S _ {\theta} ^ {\text {f i n e}} = r _ {i} \triangle \theta <   \triangle r / 2, \tag {5.20}
$$

$$
\triangle S _ {\phi} ^ {\mathrm {f i n e}} = r _ {i} \sin \theta_ {j} \triangle \phi <   \triangle r / 2.
$$

In step (II), a new discretization in polar and azimuthal angles is computed independently to produce a coarser polar mesh $\triangle \theta ^ { \mathrm { c o a r s e } } = \triangle \theta \times M _ { \theta }$ and a coarser 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/8da223d0fa732be573227ff4d985ecba41ae91e8a15da760252549d7a73c3056.jpg)



Figure 5.10: Illustration of 3-D coarse mesh generation near the origin and along the poles in different radial subdomains.


azimuthal mesh $\bigtriangleup \phi ^ { \mathrm { c o a r s e } } = \bigtriangleup \phi \times M _ { \phi }$ through the binary re-combination. This choice of coarse mesh generation algorithm restricts the original resolution in the polar angle $N _ { \theta } = 2 ^ { a }$ and azimuthal angle $N _ { \phi } = 2 ^ { b }$ to require positive non-zero integers $a$ and $b$ respectively. The even integer multipliers $M _ { \theta } = 2 ^ { L _ { \theta } }$ and $M _ { \phi } = 2 ^ { L _ { \phi } }$ are selected by the following coarse mesh conditions, where $L _ { \theta }$ and $L _ { \phi }$ are restricted to be positive integers. 

$$
\triangle S _ {\theta} ^ {\text {c o a r s e}} = r _ {i} \triangle \theta \times M _ {\theta} > \triangle r / 2, \tag {5.21}
$$

$$
\triangle S _ {\phi} ^ {\mathrm {c o a r s e}} = r _ {i} \sin \theta_ {j} \triangle \phi \times M _ {\phi} > \triangle r / 2.
$$

A general coarsening factor $f _ { \mathrm { c o a r s e } }$ is defined to specify the definition of fine mesh. For instance, Eqs. (5.20-5.21) refer to a coarsening factor equals to $1 / 2$ . The solutions for $L _ { \theta }$ and $L _ { \phi }$ in Eq. (5.21) with an arbitrary coarsening factor are 

$$
L _ {\theta} > \log \left[ f _ {\text {c o a r s e}} \triangle r / \left(r _ {i} \triangle \theta\right) \right] / \operatorname {L o g} 2, \tag {5.22}
$$

$$
{L _ {\phi}} > {\mathrm {L o g} [ f _ {\mathrm {c o a r s e}} \triangle r / (r _ {i} \sin \theta_ {j} \triangle \theta) ] / \mathrm {L o g} 2.}
$$

For solution $L _ { \theta } = L _ { \phi } = 0$ , no coarsening is required. Equation (5.22) is used to determine the new discretization in the coarser polar mesh $\triangle \theta ^ { \mathrm { c o a r s e } } = \triangle \theta \times 2 ^ { L _ { \theta } }$ and the coarser azimuthal mesh $\bigtriangleup \phi ^ { \mathrm { c o a r s e } } = \bigtriangleup \phi \times 2 ^ { L _ { \phi } }$ for each radial layer of cells 

from $r _ { 1 }$ to $r _ { N _ { r } }$ including the cells at the origin and along the poles. Figure (5.10) shows the coarse mesh generation in the first and second radial subdomains to illustrate the appearance of the coarse mesh structure near the origin and along the poles respectively. In step (III), physical quantities on the fine mesh are mapped onto the coarser mesh through the volume-averaging. After the explicit Riemann solver update for hydrodynamics, new physical quantities on the coarser mesh are mapped onto the fine mesh. Coarsening treatment is only applied to the explicit hydrodynamics. Other modeling equations including thermal, alpha and radiation diffusion are solved implicitly on the fine mesh because implicit schemes are numerically stable with no restriction on the time-step size. 

# 5.3 Hydrodynamics

# 5.3.1 Conservative moving mesh

# 1-D finite volume moving-mesh

In finite-volume methods, all numerical fluxes passing through cell interfaces are treated conservatively, with the advantage to conserve the total mass, momenta and energy within and leaving the simulation domain over times. DEC2D and DEC3D implements the moving-mesh method derived for one-dimensional hyperbolic conservation laws by Fazio and LeVeque [FL03]. 

In the non-relativistic limit, the transformation is Galilean 

$$
t = \tau , \tag {5.23}
$$

$$
x = \xi + v \tau , \tag {5.24}
$$

and satisfies the following chain rules that relate the time and spatial derivatives 

between the stationary and moving frames, 

$$
\left(\frac {\partial}{\partial t}\right) _ {x} = \left(\frac {\partial}{\partial \tau}\right) _ {\xi} \left(\frac {\partial \tau}{\partial t}\right) _ {x} + \left(\frac {\partial}{\partial \xi}\right) _ {\tau} \left(\frac {\partial \xi}{\partial t}\right) _ {x}, \tag {5.25}
$$

$$
\left(\frac {\partial}{\partial x}\right) _ {t} = \left(\frac {\partial}{\partial \tau}\right) _ {\xi} \left(\frac {\partial \tau}{\partial x}\right) _ {t} + \left(\frac {\partial}{\partial \xi}\right) _ {\tau} \left(\frac {\partial \xi}{\partial x}\right) _ {t}. \tag {5.26}
$$

Substitute the Galilean transformations from Eqs. (5.23)–(5.24) into above chain rules, and apply the triple product rule to rewrite the partial derivative as $\left( { \frac { \partial \xi } { \partial t } } \right) _ { x } =$ $- { \frac { \left( { \frac { \partial x } { \partial t } } \right) _ { \xi } } { \left( { \frac { \partial x } { \partial \xi } } \right) _ { t } } }$ , the transformation for partial derivatives are 

$$
\partial_ {t} = \partial_ {\tau} - \frac {x _ {t}}{x _ {\xi}} \partial_ {\xi}, \tag {5.27}
$$

$$
\partial_ {x} = \xi_ {x} \partial_ {\xi}. \tag {5.28}
$$

Therefore, a 1-D hyperbolic conservation law 

$$
\partial_ {t} q + \partial_ {x} f = 0 \tag {5.29}
$$

with the state variable $\boldsymbol { q } ( \boldsymbol { x } , t )$ and the flux variable $f ( x , t )$ described by the stationary frame coordinates $( x , t )$ is transformed into the moving frame with the state variable ${ \hat { q } } ( \xi , \tau )$ and the flux variable $\hat { f } ( \xi , \tau )$ , 

$$
\partial_ {\tau} \hat {q} - \frac {x _ {t}}{x _ {\xi}} \partial_ {\xi} \hat {q} + \xi_ {x} \partial_ {\xi} \hat {f} = 0. \tag {5.30}
$$

Equation (5.30) can be rewritten into an exact conservative form by taking the following relations $x _ { \xi } \partial _ { \tau } \hat { q } = \partial _ { \tau } \left( x _ { \xi } \hat { q } \right) - \hat { q } \partial _ { \tau } \left( x _ { \xi } \right) , x _ { \tau } \partial _ { \xi } \hat { q } = \partial _ { \xi } \left( x _ { \tau } \hat { q } \right) - \hat { q } \partial _ { \xi } \left( x _ { \tau } \right) , \partial _ { \tau } \left( x _ { \xi } \right) =$ $x _ { \xi } \partial _ { \tau } \hat { q } = \partial _ { \tau } \left( x _ { \xi } \hat { q } \right) - \hat { q } \partial _ { \tau } \left( x _ { \xi } \right)$ $\partial _ { \xi } \left( x _ { \tau } \right)$ , $x _ { \xi } \xi _ { x } = 1$ , and $x _ { \tau } = x _ { t }$ . 

$$
\partial_ {\tau} \left(x _ {\xi} \hat {q}\right) + \partial_ {\xi} \left(\hat {f} - x _ {t} \hat {q}\right) = 0, \tag {5.31}
$$

In the second step, the exact result in Eq. (5.31) is applied to a new transformation, defined by $x _ { i } = \xi _ { i } \triangle x$ , from the physical domain $( x _ { i } , t ^ { n } )$ to a computational domain $( \xi _ { i } , t ^ { n } )$ . The time-dependent spatial discretization is given by $\triangle x ( t ^ { n } ) = L ( t ^ { n } ) / N$ , where $L ( t )$ is the length of the physical domain at the time $t$ and $N$ is the number of cells. Since the variable $\xi _ { i }$ is a dimensionless number in the computational domain, an integration of Eq. (5.31) over the computational domain of a cell 

$$
\int_ {\xi_ {i - 1 / 2}} ^ {\xi_ {i + 1 / 2}} d \xi \left[ \partial_ {\tau} \left(x _ {\xi} \hat {q}\right) + \partial_ {\xi} \left(\hat {f} - x _ {t} \hat {q}\right) \right] = 0. \tag {5.32}
$$

Since the integration domain $\left\lfloor \xi _ { i - 1 / 2 } , \xi _ { i + 1 / 2 } \right\rfloor$ is time-independent, the time derivative can be pulled out from the integral giving, 

$$
\frac {d}{d \tau} \int_ {\xi_ {i - 1 / 2}} ^ {\xi_ {i + 1 / 2}} d \xi (x _ {\xi} \hat {q}) + [ \hat {f} - x _ {t} \hat {q} ] _ {\xi_ {i - 1 / 2}} ^ {\xi_ {i + 1 / 2}} = 0. \tag {5.33}
$$

Substitute the Jacobian relation $x _ { \xi } d \xi = d x$ in Eq. (5.33), and define the cellaveraged physical quantity $Q ( x _ { i } , t )$ in the physical domain. The term $x _ { \xi }$ is called the capacity function in the wave propagation algorithm. [Lev02] When the mesh is static, one can integrate Eq. (5.31) with $x _ { t } = 0$ and pulls out the time derivative in the same manner. However, once the mesh moves, one can only pull out the time derivative in the computational domain, which explains the significance of the second step. 

$$
Q \left(x _ {i}, t\right) = \frac {1}{\triangle x (t)} \int_ {x _ {i - 1 / 2}} ^ {x _ {i + 1 / 2}} q (x, t) d x, \tag {5.34}
$$

Finally, the conservative update for Eq. (5.33) is, 

$$
\frac {d}{d t} \left[ \triangle x (t) Q (x _ {i}, t) \right] + [ f - x _ {t} q ] _ {x _ {i - 1 / 2}} ^ {x _ {i + 1 / 2}} = 0. \tag {5.35}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/689fec856d3cb86fabb68388f0d32281a01578311110f759fcf9790a30e1a1e6.jpg)


$$
Q _ {i} ^ {n + 1} = \underbrace {\frac {\triangle_ {i - 1 / 2} ^ {n}}{\triangle x _ {i} ^ {n + 1}} Q _ {i - 1} ^ {n} + \frac {\triangle_ {i + 1 / 2} ^ {n}}{\triangle x _ {i} ^ {n + 1}} Q _ {i} ^ {n}} _ {\text {l i n e a r i n t e r p o l a t i o n}} = \underbrace {\frac {\triangle x _ {i} ^ {n}}{\triangle x _ {i} ^ {n + 1}} Q _ {i} ^ {n}} _ {\text {m o v i n g m e s h u p d a t e}} - \underbrace {\frac {\triangle t}{\triangle x _ {i} ^ {n + 1}} \left[ \dot {x} _ {i - 1 / 2} ^ {n} Q _ {i - 1} ^ {n} - \dot {x} _ {i + 1 / 2} ^ {n} Q _ {i} ^ {n} \right]}
$$


Figure 5.11: The principle of 1-D finite-volume moving-mesh method. The cell interfaces are moving to the left in the next time step level $t ^ { n + 1 }$ so that the cellaveraged quantity $Q _ { i } ^ { n + 1 }$ is defined by the linear interpolation to account for the contribution of $Q _ { i - 1 } ^ { n }$ and $Q _ { i } ^ { n }$ in the previous time-step level $t ^ { n }$ .


Equation (5.35) is the exact result for an explicit update for the 1-D hyperbolic conservation law in Eq. (5.29) on a moving mesh. From a geometrical point of view, the physical origin of moving mesh numerical fluxes $x _ { t } q$ appears naturally from the definition of cell-averaging. Figure (5.11) shows an example of a mesh moving to the left at the next time level $t ^ { n + 1 }$ with negative cell interface velocities ${ \dot { x } } _ { i \pm 1 / 2 } ^ { n } < 0$ . Three cell-averaged quantities $Q _ { i - 1 } ^ { n }$ , $Q _ { i } ^ { n }$ and $Q _ { i + 1 } ^ { n }$ at the time level $t ^ { n }$ are located at positions $x _ { i - 1 } ^ { n }$ , $\boldsymbol { x } _ { i } ^ { n }$ and $x _ { i + 1 } ^ { n }$ respectively. Since the cell interfaces at $x _ { i \pm 1 / 2 } ^ { \prime \imath }$ move to the left to new coordinates at $x _ { i \pm 1 / 2 } ^ { n + 1 }$ . In the finite volume method approach, the cell-averaging of the quantity $Q _ { i } ^ { n + 1 }$ bounded by space $[ x _ { i - 1 / 2 } ^ { n + 1 } , x _ { i + 1 / 2 } ^ { n + 1 } ]$ at time $t ^ { n + 1 }$ is defined as, 

$$
Q _ {i} ^ {n + 1} = \frac {\triangle_ {i - 1 / 2} ^ {n}}{\triangle x _ {i} ^ {n + 1}} Q _ {i - 1} ^ {n} + \frac {\triangle_ {i + 1 / 2} ^ {n}}{\triangle x _ {i} ^ {n + 1}} Q _ {i} ^ {n}, \tag {5.36}
$$

where the ratios $\frac { \triangle _ { i - 1 / 2 } ^ { n } } { \triangle x _ { i } ^ { n + 1 } }$ and 4 x n +1i $\frac { \triangle _ { i + 1 / 2 } ^ { n } } { \triangle x _ { i } ^ { n + 1 } }$ measure the weights of $Q _ { i - 1 } ^ { n }$ and $Q _ { i } ^ { n }$ respec-

tively, and $\dot { x } _ { i \pm 1 / 2 } ^ { n } \triangle t$ is the displacement of the cell interface over one time-step 4t. Substitute the distances traveled by the left and right cell interfaces at xni−1/2 $\triangle t$ $x _ { i - 1 / 2 } ^ { \prime \imath }$ and $x _ { i + 1 / 2 } ^ { \prime \imath }$ respectively, 

$$
\triangle_ {i - 1 / 2} ^ {n} = - \dot {x} _ {i - 1 / 2} ^ {n} \bigtriangleup t, \tag {5.37}
$$

$$
\triangle_ {i + 1 / 2} ^ {n} = \triangle x _ {i} ^ {n} + \dot {x} _ {i + 1 / 2} ^ {n} \triangle t. \tag {5.38}
$$

into the linear interpolation in Eq. (5.36) 

$$
Q _ {i} ^ {n + 1} = \frac {\triangle x _ {i} ^ {n}}{\triangle x _ {i} ^ {n + 1}} Q _ {i} ^ {n} - \frac {\triangle t}{\triangle x _ {i} ^ {n + 1}} \left(\dot {x} _ {i - 1 / 2} ^ {n} Q _ {i - 1} ^ {n} - \dot {x} _ {i + 1 / 2} ^ {n} Q _ {i} ^ {n}\right). \tag {5.39}
$$

Equation (5.39) is the same explicit update as Eq. (5.35), except with the zero physical flux $f = 0$ , and is also known as the upwind update because the linear interpolation in Eq. (5.36) uses the solutions of $Q _ { i - 1 } ^ { n }$ and $Q _ { i } ^ { n }$ from the previous time level $t ^ { n }$ . In the moving-mesh update by Eq. (5.39), motions of cell interfaces are known in prior to the left, and is first-order accurate in time and space. The explicit update of moving-mesh numerical fluxes can either be operator-split from hydrodynamics, or simultaneously being updated with the hydrodynamics. The latter results in less numerically diffusive solutions and are adopted in DEC2D and DEC3D. 

# 1-D wave-propagation algorithm

In order to attain a smooth transition to Riemann solvers, an upwind method following LeVeque wave-propagation algorithm [Lev02] was derived at second-order accurate in space and time for the operator-split moving-mesh update with zero physical fluxes. 

$$
\partial_ {\tau} \left(x _ {\xi} \hat {q}\right) - \partial_ {\xi} \left(x _ {t} \hat {q}\right) = 0. \tag {5.40}
$$

This moving mesh equation has the same feature of a 1-D scalar advection problem for a state variable $\bar { q }$ , which is described by a linear system of hyperbolic conservation laws with a constant-coefficient matrix $\hat { A }$ . 

$$
\partial_ {t} \vec {q} + \hat {A} \cdot \partial_ {x} \vec {q} = 0, \tag {5.41}
$$

where the commutator $[ \partial _ { x } ( { \hat { A } } { \cdot } { \vec { q } } ) , { \hat { A } } { \cdot } \partial _ { x } { \vec { q } } ]$ commutes. The 1-D scalar advection equation can be rewritten into equations described by a diffusion and wave dispersion terms by taking time derivatives on both sides of Eq. (5.41). 

$$
\partial_ {t t} \vec {q} = \hat {A} ^ {2} \cdot \partial_ {x x} \vec {q}, \tag {5.42}
$$

$$
\partial_ {t t t} \vec {q} = - \hat {A} ^ {3} \cdot \partial_ {x x x} \vec {q}. \tag {5.43}
$$

The eigenvalues $\lambda$ of the matrix $\hat { A }$ through the diagonalization of $\hat { A } = \hat { R } \cdot \Lambda \cdot \hat { R } ^ { - 1 }$ by the right-eignvector matrix $\hat { R }$ in Eq. (5.41) corresponds to the cell interface velocity $x _ { t }$ in Eq. (5.40). Let $\Lambda ^ { + }$ be a matrix containing the positive eigenvalues on the diagonal with negative ones replaced by zeros, and vice versa to define the $\Lambda ^ { - }$ matrix, so as to define the matrices for $\hat { A } ^ { + } = \hat { R } \cdot \Lambda ^ { + } \cdot \hat { R } ^ { - 1 }$ and $\hat { A } ^ { - } = \hat { R } \cdot \Lambda ^ { - } \cdot \hat { R } ^ { - 1 }$ , which obey 

$$
\hat {A} = \hat {A} ^ {+} + \hat {A} ^ {-}, \tag {5.44}
$$

$$
\left| \hat {A} \right| = \hat {A} ^ {+} - \hat {A} ^ {-}, \tag {5.45}
$$

$$
\hat {A} ^ {\pm} = \frac {1}{2} \left[ \hat {A} \pm | \hat {A} | \right], \tag {5.46}
$$

$$
\hat {A} ^ {2} = | \hat {A} | ^ {2}. \tag {5.47}
$$

The last equality for matrices $\hat { A } ^ { 2 } = \hat { R } \cdot \hat { \Lambda } ^ { 2 } \cdot \hat { R } ^ { - 1 }$ and $| \hat { A } | ^ { 2 }$ can be validated by comparing the elements of the diagonalized matrices $\{ \hat { \Lambda } ^ { 2 } \} _ { i j } = \lambda _ { i } ^ { 2 } \delta _ { i j }$ and $\{ ( \hat { \Lambda } ^ { + } -$ $\hat { \Lambda } ^ { - } ) { \cdot } ( \hat { \Lambda } ^ { + } { - } \hat { \Lambda } ^ { - } ) \} _ { i j } = | \lambda _ { i } | ^ { 2 } \delta _ { i j }$ , where $\delta _ { i j }$ is Kronecker delta which equals to one for $i =$ 

$j$ ; otherwise zero, and $\lambda _ { i }$ is the eigenvalue at the $i$ -row. Equations (5.44)–(5.47) will be used to derive LeVeque wave-propagation algorithm to be discussed latter. For compressible fluid dynamics, classical higher order $> 1$ finite-differencing schemes are numerically unstable[Sod78]. To understand this basic concept, consider a third-order Taylor expansion in time for a fluid variable $\boldsymbol { q } ( \boldsymbol { x } , t )$ from time levels $t ^ { n }$ to $t ^ { n + 1 }$ at a fixed space position $x$ , 

$$
\vec {q} (x, t ^ {n + 1}) = \vec {q} (x, t ^ {n}) + \triangle t \partial_ {t} q (x, t ^ {n}) + \frac {\triangle t ^ {2}}{2} \partial_ {t t} \vec {q} (x, t ^ {n}) + \frac {\triangle t ^ {3}}{6} \partial_ {t t t} \vec {q} (x, t ^ {n}) + \mathcal {O} _ {3}, \tag {5.48}
$$

which can be shown equal to 

$$
\vec {q} ^ {n + 1} = \vec {q} ^ {n} - \triangle t \underbrace {\hat {A} \cdot \partial_ {x} \vec {q} ^ {n}} _ {\text {a d v e c t i o n}} + \frac {\triangle t ^ {2}}{2} \underbrace {\hat {A} ^ {2} \cdot \partial_ {x x} \vec {q} ^ {n}} _ {\text {d i f f u s i o n}} - \frac {\triangle t ^ {3}}{6} \underbrace {\hat {A} ^ {3} \cdot \partial_ {x x x} \vec {q} ^ {n}} _ {\text {d i s p e r s i o n}} + \mathcal {O} _ {3}, \tag {5.49}
$$

by substituting the wave transport relations from Eq. (5.41) to Eq. (5.43). Equation (5.49) states that a first-order finite-difference scheme is numerically stable because the truncation error term $\partial _ { x x } \vec { q }$ is diffusive, but a second-order scheme is numerically unstable because the truncation error term $\partial _ { x x x } \vec { q }$ is dispersive. 

At the first-order finite-difference discretization, 

$$
\frac {\vec {q} ^ {n + 1} - \vec {q} ^ {n}}{\triangle t} = - \hat {A} \cdot \partial_ {x} \vec {q} ^ {n} + \frac {\triangle t}{2} \underbrace {\hat {A} ^ {2} \cdot \partial_ {x x} \vec {q} ^ {n}} _ {\text {d i f f u s i o n}}, \tag {5.50}
$$

the update equation contains not only an advection term but also includes a numerical diffusion coefficient $\begin{array} { r } { D _ { \mathcal { O } _ { 1 } ( \triangle t ) } = \frac { \triangle t } { 2 } \hat { A } ^ { 2 } \sim \frac { \triangle t } { 2 } \hat { \lambda _ { i } } ^ { 2 } } \end{array}$ due to the first-order explicit time discretization, [Lev02] where $D _ { \mathcal { O } _ { 1 } ( \triangle t ) }$ has the dimension of length2/time. The strength of the intrinsic numerical diffusion depends on the characteristic wave velocity $\lambda _ { i }$ . 

This is the starting point for the classical method of artificial numerical viscosity [VR50] by VonNeumann and Richtmyer in the early 1950s, which introduces 

one or more numerical diffusive terms into fluid equations to suppress the numerical noises triggered by second-order finite-difference discretization such as MacCormack scheme [Mac] in 1970s. See the survey of classical finite difference methods in 1980s by Sod. [Sod78] However, the main drawback is that the magnitudes of artificial numerical viscosities require users’ input, and the correct magnitudes that are large enough to suppress numerical noises are not known in prior to simulations. 

At the second-order finite-difference discretization, 

$$
\frac {\vec {q} ^ {n + 1} - \vec {q} ^ {n}}{\triangle t} = - \hat {A} \cdot \partial_ {x} \vec {q} ^ {n} + \frac {\triangle t}{2} \hat {A} ^ {2} \cdot \partial_ {x x} \vec {q} ^ {n} - \frac {\triangle t ^ {2}}{6} \underbrace {\hat {A} ^ {3} \cdot \partial_ {x x x} \vec {q} ^ {n}} _ {\text {d i s p e r s i o n}}, \tag {5.51}
$$

as the fluid transits from incompressible to compressible regimes, the steepening of nonlinear waves [Whi74] occurs so that the wave propagation from behind continuously catch up the wave propagation ahead, because of the gradually increasing characteristic wave velocity $\lambda$ from waves at behind to waves ahead, and eventually forms a shock wave characterized by a sharp discontinuity. 

The finite-difference approximation in Eq. (5.49) is known as the weak form[CF76], meaning that spatial partial derivatives such as $\partial _ { x x x } \vec { q }$ are not differentiable at any discontinuity. However, discontinuous solutions such as shock waves and contact discontinuous are admitted [CF76, Whi74] in Euler equations, so that numerical methods for computational fluid dynamics (CFD) must be numerically stable to transport a step function in time and space in a simple sense. 

Godunov [God59] in 1960s proposed the strong form to solve the integral representation of Euler equations by treating the fluid data between two adjacent finite-volume cells as a Riemann problem. However, Godunov schemes are only first-order accurate and is too diffusive in practices as explained by Eq. (5.50). 

The turning point was due to the introduction of flux limiters in 1980s that suppress the noise amplifications for high-order $> 1$ schemes in a non-linear manner, 

first observed by Roe [Roe84, Roe97, Roe86] in his flux-vector splitting scheme, as well as by Boris and Book in their flux-corrected transport (FCT) scheme. [BB73] The functional forms of flux limiters are required to satisfy a set of algebraic rules, which are called as the total variation diminishing (TVD) property first introduced by Harten [Har83] in 1980s, derived from explicit discretization of hyperbolic conservation laws, so as to impose a strong condition of monotonically varying in space for all profiles of fluid variables in the high-order schemes. Because the monotonic solution is a fundamental property for Euler equations. [CF76, Tor09a] 

The first genuine second-order Godunov scheme, called MUSCL, was formulated by van Leer in 1980s [vL79], who first successfully translated the physics of monotone solutions in terms of a clear geometric picture to design a monotone flux limiter. The third-order Godunov scheme called PPM was derived by Collela soon later in 1980s [CW84]. Variants of applying the TVD properties in finite-difference schemes emerged since 1980s, namely TVD-MacCormack scheme by Davis [Dav84, Dav87], the explicit and implicit compact schemes by Yee [Yee85, Yee87, Yee94, Yee97], and the essentially non-oscillatory (ENO) schemes by Harten [Har89]. The construction of TVD finite-difference schemes can be boosted to high order > 3 more readily than solving for geometric rules to construct high order finite-volume approaches. 

The Riemann problem between two adjacent fluid data can be solved either exactly by computational expensive non-linear methods or approximation methods. A family of approximated Riemann solvers were developed since 1980s, namely the two-wave model approximate Harten-Lax-Van Lee (HLL) Riemann solver [HLL83] and the three-wave model approximate Harten-Lax-Van Lee-Contact HLLC (C for including for the contact discontinuity wave in HLL) Riemann solver by Toro [Tor09a], and the five-wave model approximate HLLD Riemann solver by Mignone [MUB08] for magnetohydrodynamics. HLLC is implemented in DEC2D and DEC3D because it contains all three types of characteristic waves in Euler equa-

tions including the non-linear rarefaction wave, the shock wave and the contact discontinuity wave. 

In the following, the main concepts of LeVeque wave-propagation algorithm [Lev02] is recalled, which is served to provide a smooth transition to understand the origin of flux limiters. Substitute the second-order center-in-space discretization for partial derivatives $\partial _ { x }$ and $\partial _ { x x }$ in Eq. (5.49), 

$$
\partial_ {x} \vec {q} _ {i} ^ {n} = \frac {\vec {q} _ {i + 1} ^ {n} - \vec {q} _ {i - 1} ^ {n}}{2 \triangle x}, \tag {5.52}
$$

$$
\partial_ {x x} \vec {q} _ {i} ^ {n} = \frac {\vec {q} _ {i + 1} ^ {n} - 2 \vec {q} _ {i} ^ {n} + \vec {q} _ {i - 1} ^ {n}}{\triangle x}, \tag {5.53}
$$

such that the second-order finite-difference scheme 

$$
\vec {q} _ {i} ^ {n + 1} = \vec {q} _ {i} ^ {n} - \triangle t \hat {A} \cdot \partial_ {x} \vec {q} _ {i} ^ {n} + \frac {\triangle t ^ {2}}{2} \hat {A} ^ {2} \cdot \partial_ {x x} \vec {q} _ {i} ^ {n} + \mathcal {O} _ {2}, \tag {5.54}
$$

is reduced to Lax-Wendroff discretization, 

$$
\vec {q} _ {i} ^ {n + 1} = \vec {q} _ {i} ^ {n} - \frac {\triangle t}{2 \triangle x} \hat {A} \cdot \left(\vec {q} _ {i + 1} ^ {n} - \vec {q} _ {i - 1} ^ {n}\right) + \frac {\triangle t ^ {2}}{2 \triangle x} \hat {A} ^ {2} \cdot \left(\vec {q} _ {i + 1} ^ {n} - 2 \vec {q} _ {i} ^ {n} + \vec {q} _ {i - 1} ^ {n}\right) + \mathcal {O} _ {2}, \quad (5. 5 5)
$$

where $\begin{array} { r } { \mathcal { O } _ { 2 } = - \frac { \triangle t ^ { 3 } } { 6 } \hat { A } ^ { 3 } { \cdot } \partial _ { x x x } \vec { q } _ { i } ^ { n } } \end{array}$ is the dispersive truncation term. Define Lax-Wendroff numerical fluxe s F~ LW $\vec { F } _ { i \pm 1 / 2 } ^ { \mathrm { L W } }$ at the cell interfaces $x _ { i \pm 1 / 2 }$ , 

$$
\vec {F} _ {i + 1 / 2} ^ {\mathrm {L W}} = \frac {1}{2} \hat {A} \cdot \left(\vec {q} _ {i} ^ {n} + \vec {q} _ {i + 1} ^ {n}\right) - \frac {\triangle t}{2} \hat {A} ^ {2} \cdot \left(\vec {q} _ {i + 1} ^ {n} - \vec {q} _ {i} ^ {n}\right), \tag {5.56}
$$

$$
\vec {F} _ {i - 1 / 2} ^ {\mathrm {L W}} = \frac {1}{2} \hat {A} \cdot \left(\vec {q} _ {i - 1} ^ {n} + \vec {q} _ {i} ^ {n}\right) - \frac {\triangle t}{2} \hat {A} ^ {2} \cdot \left(\vec {q} _ {i} ^ {n} - \vec {q} _ {i - 1} ^ {n}\right), \tag {5.57}
$$

and expand 

$$
\frac {1}{2} \hat {A} \cdot \left(\vec {q} _ {i} ^ {n} + \vec {q} _ {i + 1} ^ {n}\right) = \hat {A} ^ {+} \vec {q} _ {i} ^ {n} + \hat {A} ^ {-} \vec {q} _ {i + 1} ^ {n} + \frac {1}{2} | \hat {A} | \left(\vec {q} _ {i + 1} ^ {n} - \vec {q} _ {i} ^ {n}\right), \tag {5.58}
$$

$$
\frac {1}{2} \hat {A} \cdot \left(\vec {q} _ {i - 1} ^ {n} + \vec {q} _ {i} ^ {n}\right) = \hat {A} ^ {+} \vec {q} _ {i - 1} ^ {n} + \hat {A} ^ {-} \vec {q} _ {i} ^ {n} + \frac {1}{2} | \hat {A} | \left(\vec {q} _ {i} ^ {n} - \vec {q} _ {i - 1} ^ {n}\right). \tag {5.59}
$$


van Leer


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/2d03a0f68decd33ae01376a8e04a1b26993754d55633592cd0a89918442cabbd.jpg)



Figure 5.12: van Leer’s geometrical interpolation for flux limiters. By defining the slopes $\triangle Q _ { i + 1 / 2 } = Q _ { i + 1 } - Q _ { i }$ and $\triangle Q _ { i - 1 / 2 } = Q _ { i } - Q _ { i - 1 }$ , a cell is considered as a numerical noise if the product of slopes is negative $\triangle Q _ { i + 1 / 2 } \triangle Q _ { i - 1 / 2 } < 0$ , and no cell reconstruction is needed $\triangle Q _ { i } = 0$ . Otherwise, the cell profile is monotonic and is reconstructed by addling a linear slope $\triangle Q _ { i } = \left( \triangle Q _ { i + 1 / 2 } + \triangle Q _ { i - 1 / 2 } \right) / 2$ .


Equations (5.58)–(5.59) imply that the first-order upwind solution contains a numerical diffusion coefficient $D _ { { \mathcal O } _ { 1 } ( \triangle x ) } = \frac { 1 } { 2 } | \hat { A } | \triangle x$ [Roe84] due to the spatial discretization in Lax-Wendroff scheme. Similarly, the first-order but not upwind, Lax-Friedrichs scheme [Tor09a] has a numerical diffusion coefficient $D _ { \mathcal { O } _ { 1 } ( \triangle x ) } =$ $\scriptstyle { \frac { 1 } { 2 } } \triangle x ^ { 2 } / \triangle t$ , [Lev02] which only depend on the mesh and time-step sizes but not characteristic wave velocities. Define the differences of state variables $\triangle { \vec { q } } _ { i + 1 / 2 } =$ $\left( \vec { q } _ { i + 1 } ^ { n } - \vec { q } _ { i } ^ { n } \right)$ and $\triangle \vec { q } _ { i - 1 / 2 } = \left( \vec { q } _ { i } ^ { n } - \vec { q } _ { i - 1 } ^ { n } \right)$ across the cell interfaces at $x _ { i \pm 1 / 2 }$ , Lax-Wendroff cell-interface numerical fluxes in Eqs. (5.56)–(5.57) are simplified into 

$$
\vec {F} _ {i + 1 / 2} ^ {\mathrm {L W}} = \left(\hat {A} ^ {+} \vec {q} _ {i} ^ {n} + \hat {A} ^ {-} \vec {q} _ {i + 1}\right) + \frac {1}{2} | \hat {A} | \left(\hat {I} - \frac {\Delta t}{\Delta x} | \hat {A} |\right) \Delta \vec {q} _ {i + 1 / 2}, \tag {5.60}
$$

$$
\vec {F} _ {i - 1 / 2} ^ {\mathrm {L W}} = \left(\hat {A} ^ {+} \vec {q} _ {i - 1} ^ {n} + \hat {A} ^ {-} \vec {q} _ {i}\right) + \frac {1}{2} | \hat {A} | \left(\hat {I} - \frac {\triangle t}{\triangle x} | \hat {A} |\right) \triangle \vec {q} _ {i - 1 / 2}. \tag {5.61}
$$

By replacing the state difference $\triangle \vec { q } _ { i \pm 1 / 2 }$ with a total-variation-diminishing (TVD) flux-limited state difference $\triangle \vec { q } _ { i \pm 1 / 2 } ^ { \mathrm { r } \mathrm { 1 } \mathrm { \upsilon } \mathrm { 1 } }$ in Eqs. (5.60)–(5.61), the resulting fluxlimited second-order Lax-Wendroff scheme in Eq. (5.55) is called LeVeque wave propagation algorithm. Comprehensive analysis for functional forms of TVD flux limiters were reported by Sweby [Swe84]. Figure (5.12) shows the exam-

ple of van Leer TVD flux limiter to explain flux limiters from a simple geometrical point of view. The basic idea is to detect whether a cell belongs to numerical noises or solutions of Euler equations that exhibiting monotone profiles in space for all fluid variables. The detection is accomplished by checking the product of slopes $\triangle q _ { i - 1 / 2 } \triangle q _ { i + 1 / 2 }$ between two adjacent cells, where $\triangle q _ { i - 1 / 2 } =$ $q _ { i } \mathrm { ~ - ~ } q _ { i - 1 }$ and $\triangle q _ { i + 1 / 2 } = q _ { i + 1 } - q _ { i }$ . A cell is noise if the product is negative $\triangle q _ { i - 1 / 2 } \triangle q _ { i + 1 / 2 } < 0$ , otherwise the cell is reconstructed by adding a linear slope $\triangle q _ { i } = \left( \triangle q _ { i + 1 / 2 } + \triangle q _ { i - 1 / 2 } \right) / 2$ . Therefore, one of the TVD flux limiters derived by van Leer is, 

$$
\begin{array}{l} \triangle \bar {q} _ {i} ^ {\mathrm {T V D}, (\text {v a n L e e r})} = \frac {1}{2} \left(\triangle q _ {i + 1 / 2} + \triangle q _ {i - 1 / 2}\right), \quad \text {i f} \triangle q _ {i - 1 / 2} \triangle q _ {i + 1 / 2} > 0, (5.62) \\ = 0, \quad \text {i f} \triangle q _ {i - 1 / 2} \triangle q _ {i + 1 / 2} <   0. (5.63) \\ \end{array}
$$

This is called MUSCL scheme. [vL79] For those cell being detected as numerical noises are updated with zero slopes $\triangle \vec { q } _ { i \pm 1 / 2 }$ in the upwind numerical fluxes in Eqs. (5.60)–(5.61), implying that all cells with a noisy background are updated at the first-order numerically diffusive method whereas all other cells with monotone profiles are updated at second-order. Figure (5.13) compares the numerical stability in resolving the mass density profile across a sharp inner interface in a NIF implosion between a classical finite-difference MacCormack scheme with artificial numerical viscosities and MUSCL scheme. The principle of flux limiter is to switch second-order solution into first-order numerically diffusive solution to damp the growth of numerical noises. 

$$
\vec {q} _ {i} ^ {n + 1} = \vec {q} _ {i} ^ {n} - \frac {\triangle t}{\triangle x} \left(\vec {F} _ {i + 1 / 2} ^ {\mathrm {L W}} - \vec {F} _ {i - 1 / 2} ^ {\mathrm {L W}}\right). \tag {5.64}
$$

The last piece to attain a second-order TVD scheme to update the operator-


5% initial perturbation


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/f028f44ed01fab41215457f7343f8e50da5c9d5d8982c7ece7298dbd6c4eaf78.jpg)



Figure 5.13: Significant numerical noises are observed in MacCormack scheme across the sharp shell interface, whereas the solution of MUSCL scheme HLLC approximate Riemann solver is numerically stable.


split moving-mesh equation in Eq. (5.40) is to introduce the second-order Lax-Wendroff TVD numerical fluxes F~ LW(O2)i±1/2 in Eq. (5.39) by replacing the matrix $\vec { F } _ { i \pm 1 / 2 } ^ { \mathrm { L W } ( \mathcal { O } _ { 2 } ) }$ of characteristic wave velocities $| \hat { A } |$ with the cell interface mesh velocity $| { \dot { x } } _ { i \pm 1 / 2 } ^ { n } |$ , because the cell interface velocity is the characteristic wave velocity in the movingmesh equation in Eq. (5.40). 

$$
\vec {F} _ {i \pm 1 / 2} ^ {\mathrm {L W} (\mathcal {O} _ {2})} = \frac {1}{2} | \hat {A} | \left(\hat {I} - \frac {\triangle t}{\triangle x} | \hat {A} |\right) \triangle \bar {q} _ {i \pm 1 / 2} ^ {\mathrm {T V D}}. \tag {5.65}
$$

Because the upwind numerical fluxes in Eq. (5.39) in the example of left-moving cell interfaces can be generalized as the first-order Lax-Wendroff numerical fluxes F~ LW(O1)i±1/2 in Eqs. (5.60)–(5.61), $\vec { F } _ { i \pm 1 / 2 } ^ { \mathrm { L W } ( \mathcal { O } _ { 1 } ) }$ 

$$
\vec {F} _ {i + 1 / 2} ^ {\mathrm {L W} \left(\mathcal {O} _ {1}\right)} = \hat {A} ^ {+} \vec {q} _ {i} ^ {n} + \hat {A} ^ {-} \vec {q} _ {i + 1}, \tag {5.66}
$$

$$
\vec {F} _ {i - 1 / 2} ^ {\mathrm {L W} \left(\mathcal {O} _ {1}\right)} = \hat {A} ^ {+} \vec {q} _ {i - 1} ^ {n} + \hat {A} ^ {-} \vec {q} _ {i}. \tag {5.67}
$$

For a moving mesh with all cell interfaces move to the left systematically, the solution for the first-order Lax-Wendroff numerical fluxes is given by substituting 

$$
\hat {A} ^ {+} = - \dot {x} _ {i \pm 1 / 2} ^ {n} > 0 \mathrm {a n d} \hat {A} ^ {-} = 0,
$$

$$
\vec {F} _ {i + 1 / 2} ^ {\mathrm {L W} \left(\mathcal {O} _ {1}\right), \mathrm {L}} = - \dot {x} _ {i + 1 / 2} ^ {n} \vec {q} _ {i} ^ {n}, \tag {5.68}
$$

$$
\vec {F} _ {i - 1 / 2} ^ {\mathrm {L W} \left(\mathcal {O} _ {1}\right), \mathrm {L}} = - \dot {x} _ {i - 1 / 2} ^ {n} \vec {q} _ {i - 1} ^ {n}, \tag {5.69}
$$

and $\hat { A } ^ { + } = 0$ and $\hat { A } ^ { - } = - \dot { x } _ { i \pm 1 / 2 } ^ { n } < 0$ for a right-moving mesh, 

$$
\vec {F} _ {i + 1 / 2} ^ {\mathrm {L W} \left(\mathcal {O} _ {1}\right), \mathrm {R}} = - \dot {x} _ {i + 1 / 2} ^ {n} \vec {q} _ {i + 1}, \tag {5.70}
$$

$$
\vec {F} _ {i - 1 / 2} ^ {\mathrm {L W} \left(\mathcal {O} _ {1}\right), \mathrm {R}} = - \dot {x} _ {i - 1 / 2} ^ {n} \vec {q} _ {i}. \tag {5.71}
$$

The first-order upwind scheme for Eq. (5.39) in the example of a left moving-mesh can be written in first-order Lax-Wendroff form, 

$$
Q _ {i} ^ {n + 1} = \frac {\triangle x _ {i} ^ {n}}{\triangle x _ {i} ^ {n + 1}} Q _ {i} ^ {n} - \frac {\triangle t}{\triangle x _ {i} ^ {n + 1}} \left(F _ {i + 1 / 2} ^ {\mathrm {L W} (\mathcal {O} _ {1}), \mathrm {L}} - F _ {i - 1 / 2} ^ {\mathrm {L W} (\mathcal {O} _ {1}), \mathrm {L}}\right), \tag {5.72}
$$

and the second-order form of LeVeque wave-propagation algorithm is, 

$$
Q _ {i} ^ {n + 1} = \underbrace {\frac {\triangle x _ {i} ^ {n}}{\triangle x _ {i} ^ {n + 1}}} _ {\text {g e o m e t r y} 1} Q _ {i} ^ {n} - \underbrace {\frac {\triangle t}{\triangle x _ {i} ^ {n + 1}}} _ {\text {g e o m e t r y} 2} \left(F _ {i + 1 / 2} ^ {\mathrm {L W} (\mathcal {O} _ {1}), \mathrm {L}} + \vec {F} _ {i + 1 / 2} ^ {\mathrm {L W} (\mathcal {O} _ {2})} - F _ {i - 1 / 2} ^ {\mathrm {L W} (\mathcal {O} _ {1}), \mathrm {L}} - \vec {F} _ {i - 1 / 2} ^ {\mathrm {L W} (\mathcal {O} _ {2})}\right). \tag {5.73}
$$

Figure (5.14) compares the performance of 1-D LeVeque wave-propagation algorithm by Eq. (5.73) between a slow moving-mesh at ${ \dot { x } } = - 1$ and a fast movingmesh at $\dot { x } = - 2$ . Good agreements with the exact solution on static mesh are observed for the fluid variables of the mass density $\rho$ , velocity $u$ and pressure $P$ in a 1-D planner shock tube problem. Although DEC2D and DEC3D do not implement the operator-split moving-mesh update because of relatively more diffusive solutions than solving the moving-mesh update in Riemann solver simultaneously, the value of Eq. (5.73) explains how two places of the cell geometry enter the 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/a6bcf9c8c33344d87e1621608bce2975aed9ebe50ff02655b88702913f237394.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/c535c79b48955d7d2296d968399a8771f448dc14420fab17c541fddb935e624f.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/9a5c6eb4b5b7f2dfa71ec33261055ce927d54245efbb47405a0b66f7a715e03a.jpg)



Figure 5.14: Comparison of 1-D LeVeque wave-propagation algorithm between a slow moving-mesh at ${ \dot { x } } = - 1$ and a fast moving-mesh at $\dot { x } = - 2$ that updates a operator-split moving mesh using Eq. (5.73), followed by a hydro update using the second-order MUSCL-HLLC Riemann solver in a shock tube problem.


explicit update of finite-volume moving mesh. These geometrical factors result from the capacity function $x _ { \xi } = \triangle x$ in Eq. (5.33). Without loss of generality, this result is applied to the 3-D finite-volume moving-mesh update. 

To understand the physical meaning of constructing an upwind solution, the total change for a cell content $Q _ { i } ^ { n + 1 } - Q _ { i } ^ { n }$ over one time-step $\triangle t$ for a first-order Lax-Wendroff update in Eq. (5.72) on a static mesh $\triangle x _ { i } ^ { n } = \triangle x _ { i } ^ { n + 1 } = \triangle x$ can be shown as a linear interpolation in four terms, 

$$
Q _ {i} ^ {n + 1} - Q _ {i} ^ {n} = - \underbrace {\frac {\triangle t \hat {A} _ {i + 1 / 2} ^ {+}}{\triangle x}} _ {\text {l e a v e}} Q _ {i} ^ {n} - \underbrace {\frac {\triangle t \hat {A} _ {i + 1 / 2} ^ {-}}{\triangle x}} _ {\text {e n t e r}} Q _ {i + 1} ^ {n} + \underbrace {\frac {\triangle t \hat {A} _ {i - 1 / 2} ^ {+}}{\triangle x}} _ {\text {e n t e r}} Q _ {i - 1} ^ {n} + \underbrace {\frac {\triangle t \hat {A} _ {i - 1 / 2} ^ {-}}{\triangle x}} _ {\text {l e a v e}} Q _ {i} ^ {n}. \tag {5.74}
$$

On the right hand side of Eq. (5.74), the first term measures the leaving cell content $Q _ { i } ^ { n }$ with a weight of $\frac { \triangle t { \hat { A } } _ { i + 1 / 2 } ^ { + } } { \triangle x } > 0$ and the total sign is negative meaning 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/1ba16d8ef5c6727c66cb43e90a7bc6b89829b71e6f69494a461d2d6cf944939d.jpg)



Figure 5.15: The upwind solution in the first-order Lax-Wendroff update is equivalent to a linear interpolation to weigh all possible cell contents leaving and entering into a cell.


for a loss at the cell interface $x _ { i + 1 / 2 }$ , the second term measures the entering cell content $Q _ { i + 1 } ^ { n }$ with a weight of 4tAˆ−i+1/2 $\frac { \triangle t { \hat { A } } _ { i + 1 / 2 } ^ { - } } { \triangle x } < 0$ 4x < 0 and the total sign is positive meaning for a gain at the cell interface $x _ { i + 1 / 2 }$ , the third term measures the entering cell content $Q _ { i - 1 } ^ { n }$ with a weight of $\frac { \triangle t { \hat { A } } _ { i - 1 / 2 } ^ { + } } { \triangle x } > 0$ > 0 and the total sign is positive meaning for a gain at the cell interface $x _ { i - 1 / 2 }$ , the forth term measures the leaving cell content $Q _ { i } ^ { n }$ with a weight of 4tAˆ−i−1/2 $\frac { \triangle t { \hat { A } } _ { i - 1 / 2 } ^ { - } } { \triangle x } < 0$ 4x < 0 and the total sign is negative meaning for a loss at the cell interface $x _ { i - 1 / 2 }$ . 

Figure (5.15) illustrate the wave propagations leaving and entering a cell $Q _ { i } ^ { n }$ over the time levels from $t ^ { n }$ to $t ^ { n + 1 }$ with the position centered at $x _ { i }$ on a static mesh. The sign of characteristic wave velocity is indicated by the superscript of $\pm$ . The contribution of adjacent cell contents $Q _ { i \pm 1 } ^ { n }$ is manifested by waves propagating into the cell at $x _ { i }$ . In hydrodynamics, these waves are categorized into three types including the shock wave, the rarefaction wave and the contact discontinuity wave. 

# 3-D finite-volume moving-mesh in Cartesian geometry

In this thesis, the finite-volume method approach was applied to discretize the inviscid compressible Euler equations on both Cartesian and spherical movingmesh for DEC2D and DEC3D. The most up-to-date code development status is 

reported in Table 5.1. To discretize Euler equations in a systematic manner, the following notations for time and space variables are introduced. 

$$
x _ {0} = t, \quad x _ {1} = x, \quad x _ {2} = y, \quad x _ {3} = z, \tag {5.75}
$$

$$
\xi_ {0} = \tau , \quad \xi_ {1} = \xi_ {x}, \quad \xi_ {2} = \xi_ {y}, \quad \xi_ {3} = \xi_ {z}. \tag {5.76}
$$

The set of variables $\{ x _ { i = 0 , 1 , 2 , 3 } \}$ are coordinates in the physical domain, whereas $\{ \xi _ { i = 0 , 1 , 2 , 3 } \}$ are coordinates in the computational domain. Along a given $i$ -direction, let $L _ { i } ( t )$ be the length of the physical domain at time $t$ and let $N _ { i }$ be the number of cells for the discretization $\triangle x _ { i } ( t ) = L _ { i } ( t ) / N _ { i }$ . The general transformation of coordinates between physical and computational domains is 

$$
x _ {0} = \xi_ {0}, \tag {5.77}
$$

$$
x _ {i} = \xi_ {i} \triangle x _ {i} (t), \text {w h e r e} \triangle x _ {i} (t) = L _ {i} (t) / N _ {i}, \forall i = 1, 2, 3. \tag {5.78}
$$

Using the general chain rule relation 

$$
\frac {\partial}{\partial x _ {i}} = \sum_ {i = 0} ^ {3} \frac {\partial \xi_ {i}}{\partial x _ {i}} \frac {\partial}{\partial \xi_ {i}}, \tag {5.79}
$$

the time and spatial derivatives in the physical domain are transformed into the computational domain, 

$$
\frac {\partial}{\partial t} = \frac {\partial}{\partial \tau} + \sum_ {i = 1} ^ {3} \frac {\partial \xi_ {i}}{\partial t} \frac {\partial}{\partial \xi_ {i}}, \tag {5.80}
$$

$$
\frac {\partial}{\partial x _ {i}} = \frac {\partial \xi_ {i}}{\partial x _ {i}} \frac {\partial}{\partial \xi_ {i}}, \quad \forall i = 1, 2, 3. \tag {5.81}
$$

Define the total energy density of the fluid $\begin{array} { r } { \varepsilon = \frac { \mathrm { ~ \varsigma ~ } } { \gamma - 1 } + \frac { 1 } { 2 } \rho v ^ { 2 } } \end{array}$ , where $P$ is the fluid pressure, $\gamma = 5 / 3$ is the ratio of specific heats for an ideal gas, $\rho$ is the fluid mass density and $\boldsymbol { v }$ is the fluid velocity. Euler equations in the physical domain 

described by a set stationary coordinates $\{ x _ { i = 0 , 1 , 2 , 3 } \}$ , 

$$
\frac {\partial}{\partial t} \left[ \begin{array}{l} \rho \\ \rho \vec {v} \\ \varepsilon \end{array} \right] + \vec {\nabla} \cdot \left[ \begin{array}{c} \rho \vec {v} \\ \rho \vec {v} \otimes \vec {v} + \hat {I} P \\ \vec {v} (\varepsilon + P) \end{array} \right] = \vec {0}, \tag {5.82}
$$

can be represented by a system of vector equations that contains a vector of state variables $\vec { q } ( \vec { x } , t )$ and three vectors of physical fluxes $\vec { f _ { x } } ( \vec { x } , t )$ , $\vec { f } _ { y } ( \vec { x } , t )$ and $\vec { f _ { z } } ( \vec { x } , t )$ in $x$ , $y$ and $z$ directions respectively. 

$$
\frac {\partial \vec {q}}{\partial t} + \frac {\partial \vec {f _ {x}}}{\partial x} + \frac {\partial \vec {f _ {y}}}{\partial y} + \frac {\partial \vec {f _ {z}}}{\partial z} = 0, \tag {5.83}
$$

with the explicit form 

$$
\frac {\partial}{\partial t} \left[ \begin{array}{l} \rho \\ \rho v _ {x} \hat {x} \\ \rho v _ {y} \hat {y} \\ \rho v _ {z} \hat {z} \\ \varepsilon \end{array} \right] + \frac {\partial}{\partial x} \left[ \begin{array}{l} \rho v _ {x} \\ \rho v _ {x} v _ {x} + P \\ \rho v _ {x} v _ {y} \\ \rho v _ {x} v _ {z} \\ v _ {x} (\varepsilon + P) \end{array} \right] + \frac {\partial}{\partial y} \left[ \begin{array}{l} \rho v _ {y} \\ \rho v _ {y} v _ {x} \\ \rho v _ {y} v _ {y} + P \\ \rho v _ {y} v _ {z} \\ v _ {y} (\varepsilon + P) \end{array} \right] + \frac {\partial}{\partial z} \left[ \begin{array}{l} \rho v _ {z} \\ \rho v _ {z} v _ {x} \\ \rho v _ {z} v _ {y} \\ \rho v _ {z} v _ {z} + P \\ v _ {z} (\varepsilon + P) \end{array} \right] = \vec {0} (5. 8 4)
$$

In Cartesian coordinates, the time derivatives of stationary unit vectors $\hat { x }$ , $\hat { y }$ and $\hat { z }$ are zero, so that $\partial _ { t } \hat { x } _ { i } = 0$ for $i = 1 , 2 , 3$ which, however, is not valid for curvilinear coordinates. 

$$
\partial_ {t} \left(\rho v _ {i} \hat {x} _ {i}\right) = \hat {x} _ {i} \partial_ {t} \left(\rho v _ {i}\right) + \left(\rho v _ {i}\right) \underbrace {\partial_ {t} \hat {x} _ {i}} _ {\text {z e r o}}. \tag {5.85}
$$

Substitute the chain rules of Eqs. (5.80)–(5.81), and the property of zero time derives for Cartesian unit vectors of Eq. (5.85), into Eq. (5.83), 

$$
\frac {\partial \vec {q}}{\partial \tau} + \sum_ {i = 1} ^ {3} \frac {\partial \xi_ {i}}{\partial t} \frac {\partial \vec {q}}{\partial \xi_ {i}} + \sum_ {i = 1} ^ {3} \frac {\partial \xi_ {i}}{\partial x _ {i}} \frac {\partial \vec {f} _ {i}}{\partial \xi_ {i}} = \vec {0}. \tag {5.86}
$$

In above equation, substitute the triple product rule of $\begin{array} { r } { \frac { \partial \xi _ { i } } { \partial t } = - \frac { \partial x _ { i } } { \partial t } / \frac { \partial x _ { i } } { \partial \xi _ { i } } } \end{array}$ ∂xi to replace the time derivatives, and multiply the both sides with a term Q j =1 $\textstyle \prod _ { j = 1 } ^ { 3 } { \frac { \partial x _ { j } } { \partial \xi _ { j } } }$ ∂xj , 

$$
\left(\prod_ {j = 1} ^ {3} \frac {\partial x _ {j}}{\partial \xi_ {j}}\right) \frac {\partial \vec {q}}{\partial \tau} - \sum_ {i = 1} ^ {3} \left(\prod_ {j = 1, \neq i} ^ {3} \frac {\partial x _ {j}}{\partial \xi_ {j}}\right) \frac {\partial x _ {i}}{\partial t} \frac {\partial \vec {q}}{\partial \xi_ {i}} + \sum_ {i = 1} ^ {3} \left(\prod_ {j = 1, \neq i} ^ {3} \frac {\partial x _ {j}}{\partial \xi_ {j}}\right) \frac {\partial \vec {f _ {i}}}{\partial \xi_ {i}} = \vec {0}, (5. 8 7)
$$

The sum of first two terms in Eq. (5.87), 

$$
\left(\prod_ {j = 1} ^ {3} \frac {\partial x _ {j}}{\partial \xi_ {j}}\right) \frac {\partial \vec {q}}{\partial \tau} = \frac {\partial}{\partial \tau} \left(\prod_ {j = 1} ^ {3} \frac {\partial x _ {j}}{\partial \xi_ {j}} \vec {q}\right) - \vec {q} \frac {\partial}{\partial \tau} \left(\prod_ {j = 1} ^ {3} \frac {\partial x _ {j}}{\partial \xi_ {j}}\right),
$$

$$
\sum_ {i = 1} ^ {3} \left(\prod_ {j = 1, \neq i} ^ {3} \frac {\partial x _ {j}}{\partial \xi_ {j}}\right) \frac {\partial x _ {i}}{\partial t} \frac {\partial \vec {q}}{\partial \xi_ {i}} = \sum_ {i = 1} ^ {3} \left[ \frac {\partial}{\partial \xi_ {i}} \left(\prod_ {j = 1, \neq i} ^ {3} \frac {\partial x _ {j}}{\partial \xi_ {j}} \frac {\partial x _ {i}}{\partial t} \vec {q}\right) - \vec {q} \frac {\partial}{\partial \xi_ {i}} \left(\prod_ {j = 1, \neq i} ^ {3} \frac {\partial x _ {j}}{\partial \xi_ {j}} \frac {\partial x _ {i}}{\partial t}\right) \right],
$$

is equals to 

$$
\frac {\partial}{\partial \tau} \left(\prod_ {j = 1} ^ {3} \frac {\partial x _ {j}}{\partial \xi_ {j}} \vec {q}\right) - \sum_ {i = 1} ^ {3} \frac {\partial}{\partial \xi_ {i}} \left(\prod_ {j = 1, \neq i} ^ {3} \frac {\partial x _ {j}}{\partial \xi_ {j}} \frac {\partial x _ {i}}{\partial t} \vec {q}\right), \tag {5.88}
$$

because of the following chain rule, 

$$
\sum_ {i = 1} ^ {3} \frac {\partial}{\partial \xi_ {i}} \left(\prod_ {j = 1, \neq i} ^ {3} \frac {\partial x _ {j}}{\partial \xi_ {j}} \frac {\partial x _ {i}}{\partial t}\right) = \sum_ {i = 1} ^ {3} \left(\prod_ {j = 1, \neq i} ^ {3} \frac {\partial x _ {j}}{\partial \xi_ {j}}\right) \frac {\partial}{\partial t} \left(\frac {\partial x _ {i}}{\partial \xi_ {i}}\right) = \frac {\partial}{\partial \tau} \left(\prod_ {j = 1} ^ {3} \frac {\partial x _ {j}}{\partial \xi_ {j}}\right).
$$

The third term in Eq. (5.87) is, 

$$
\sum_ {i = 1} ^ {3} \left(\prod_ {j = 1, \neq i} ^ {3} \frac {\partial x _ {j}}{\partial \xi_ {j}}\right) \frac {\partial \vec {f} _ {i}}{\partial \xi_ {i}} = \sum_ {i = 1} ^ {3} \frac {\partial}{\partial \xi_ {i}} \left(\prod_ {j = 1, \neq i} ^ {3} \frac {\partial x _ {j}}{\partial \xi_ {j}} \vec {f} _ {i}\right) \tag {5.89}
$$

Using results of Eqs. (5.88) and (5.89), Euler equations in the computational domain in Eq. (5.87) can be written into a conservative form as follows, 

$$
\frac {\partial}{\partial \tau} \left(\prod_ {j = 1} ^ {3} \frac {\partial x _ {j}}{\partial \xi_ {j}} \vec {q}\right) + \sum_ {i = 1} ^ {3} \frac {\partial}{\partial \xi_ {i}} \left(\prod_ {j = 1, \neq i} ^ {3} \frac {\partial x _ {j}}{\partial \xi_ {j}}\right) \left(\vec {f _ {i}} - \frac {\partial x _ {i}}{\partial t} \vec {q}\right) = \vec {0}. \tag {5.90}
$$

To obtain the finite-volume discretization, integrate Eq. (5.90) over the volume of a cell in the computational domain $\begin{array} { r } {  { \mathcal { D } } _ { \xi } ~ = ~ [ \xi _ { i - 1 / 2 } , \xi _ { i + 1 / 2 } ] \times [ \xi _ { j - 1 / 2 } , \xi _ { j + 1 / 2 } ] \times } \end{array}$ $[ \xi _ { k - 1 / 2 } , \xi _ { k + 1 / 2 } ]$ , 

$$
\int_ {\mathcal {D} _ {\xi}} d \xi_ {1} d \xi_ {2} d \xi_ {3} \left[ \frac {\partial}{\partial \tau} \left(\prod_ {j = 1} ^ {3} \frac {\partial x _ {j}}{\partial \xi_ {j}} \vec {q}\right) + \sum_ {i = 1} ^ {3} \frac {\partial}{\partial \xi_ {i}} \left(\prod_ {j = 1, \neq i} ^ {3} \frac {\partial x _ {j}}{\partial \xi_ {j}}\right) \left(\vec {f _ {i}} - \frac {\partial x _ {i}}{\partial t} \vec {q}\right) \right] = \vec {0}. \tag {5.91}
$$

The integration of the first term in Eq. (5.94) is obtained by Leibniz integral rule such that the time derivative is pulled out directly, without introducing an extra surface integral term due to the property of time-independent boundaries in the computational domain. 

$$
\int_ {\mathcal {D} _ {\xi}} d \xi_ {1} d \xi_ {2} d \xi_ {3} \left[ \frac {\partial}{\partial \tau} \left(\prod_ {j = 1} ^ {3} \frac {\partial x _ {j}}{\partial \xi_ {j}} \vec {q}\right) \right] = \frac {d}{d t} \left(V _ {i j k} \vec {Q} _ {i j k}\right), \tag {5.92}
$$

where $V _ { i j k }$ is the cell volume in the physical domain $\mathcal { D } _ { x } ~ = ~ \left[ x _ { i - 1 / 2 } , x _ { i + 1 / 2 } \right] \ \times$ $[ y _ { j - 1 / 2 } , y _ { j + 1 / 2 } ] \times [ z _ { k - 1 / 2 } , z _ { k + 1 / 2 } ]$ , and $\vec { Q } _ { i j k }$ is the cell-volume-averaged state variable, 

$$
\vec {Q} _ {i j k} = \frac {1}{V _ {i j k}} \int_ {\mathcal {D} _ {x}} d x _ {1} d x _ {2} d x _ {3} \vec {q}. \tag {5.93}
$$

The integration of the second term in Eq. (5.94) for the $i = 1$ component is, 

$$
\int_ {\mathcal {D} _ {\xi}} d \xi_ {1} d \xi_ {2} d \xi_ {3} \left[ \frac {\partial}{\partial \xi_ {1}} \left(\frac {\partial x _ {2}}{\partial \xi_ {2}} \frac {\partial x _ {3}}{\partial \xi_ {3}}\right) \left(\vec {f} _ {1} - \frac {\partial x _ {1}}{\partial t} \vec {q}\right) \right] = \left[ A _ {1} \left(\vec {F} _ {1} ^ {*} - \frac {\partial x _ {1}}{\partial t} \vec {Q} ^ {*}\right) \right] _ {i - 1 / 2} ^ {i + 1 / 2}, \tag {5.94}
$$

where $A _ { 1 } = \triangle x _ { 2 } \triangle x _ { 3 }$ is the cell surface area normal to the $i = 1$ direction, where $\vec { F } _ { 1 } ^ { * }$ and $\vec { Q } ^ { * }$ are the cell-surface-averaged cell-interface numerical flux and the cellsurface-averaged cell-interface state variable respectively, 

$$
\vec {F} _ {1} ^ {*} = \frac {1}{A _ {1}} \int_ {\partial \mathcal {D} _ {x}} d x _ {2} d x _ {3} \vec {f} _ {1}, \tag {5.95}
$$

$$
\vec {Q} ^ {*} = \frac {1}{A _ {1}} \int_ {\partial \mathcal {D} _ {x}} d x _ {2} d x _ {3} \vec {q}, \tag {5.96}
$$

and the same manner to define other two cell-interface numerical fluxes $\vec { F } _ { 2 } ^ { * }$ and $\vec { F } _ { 3 } ^ { * }$ and the cell-interface state variables $\vec { Q } ^ { * }$ in $y$ and $z$ directions. The notation of $\partial { \mathcal { D } } _ { x }$ referes to the surface of the physical domain $\mathcal { D } _ { x }$ . The final form of finite-volume moving-mesh update for Euler equations in 3-D Cartesian geometry is, 

$$
\frac {\left(V _ {i j k} \vec {Q} _ {i j k}\right) ^ {n + 1} - \left(V _ {i j k} \vec {Q} _ {i j k}\right) ^ {n}}{\triangle t} = \sum_ {i = 1} ^ {3} \left[ A _ {i} ^ {n} \left(\vec {F} _ {i} ^ {*} - \dot {x} _ {i} \vec {Q} ^ {*}\right) \right] _ {i - 1 / 2} ^ {i + 1 / 2}. \tag {5.97}
$$

Figure (5.16) shows the performance of 3-D parallel simulations for mode $\ell = 1 0$ and $\ell = 2 0$ using Cartesian moving-mesh algorithm described by Eq. (5.97). 

# 3-D finite-volume moving-mesh in spherical geometry

In spherical coordinates, Euler equations are, 

$$
\frac {\partial \vec {Q}}{\partial t} + \frac {1}{r ^ {2}} \frac {\partial \left(r ^ {2} \vec {F} _ {r}\right)}{\partial r} + \frac {1}{r \sin \theta} \frac {\partial (\sin \theta \vec {F} _ {\theta})}{\partial \theta} + \frac {1}{r \sin \theta} \frac {\partial \vec {F} _ {\phi}}{\partial \phi} = \vec {S}. \tag {5.98}
$$

# Configuration of electron temperature at stagnation

$$
m = 1 0 \quad \ell = 1 0
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/88ec3c9e2da436e6591c82c6fec2a2877a3e81e8a01620e22920793551a4af30.jpg)


$$
m = 2 0 \quad \ell = 2 0
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/e07f1be82caaf9c706a1eaf07cc6ddce476967781604ea78c395f09a880fa1fc.jpg)



Figure 5.16: Performance of 3-D parallel simulations for modes $\ell = 1 0$ and $\ell = 2 0$ using the Cartesian-mesh version DEC3D implemented with the Cartesian moving-mesh algorithm described by Eq. (5.97).


Here $\vec { Q }$ is the state vector and $\vec { S }$ is the geometry source vector representing the pseudo forces due to the time derivative of momentum vectors in the spherical coordinates. 

$$
\vec {Q} = \left[ \begin{array}{c} \rho \\ \rho v _ {r} \\ \rho v _ {\theta} \\ \rho v _ {\phi} \\ \varepsilon \end{array} \right], \quad \vec {S} = \left[ \begin{array}{c} 0 \\ 2 P / r + \rho (v _ {\theta} ^ {2} + v _ {\phi} ^ {2}) / r \\ - \rho v _ {\theta} v _ {r} / r + \frac {\cos \theta}{\sin \theta} (\rho v _ {\phi} ^ {2} + P) / r \\ - \rho v _ {r} v _ {\phi} / r - \frac {\cos \theta}{\sin \theta} \rho v _ {\theta} v _ {\phi} / r \\ 0 \end{array} \right]. \tag {5.99}
$$

Here $\rho$ is the fluid mass density, $v _ { r } , v _ { \theta } , v _ { \phi }$ are fluid velocity components in $r , \theta , \phi$ directions respectively. $\varepsilon = P / ( \gamma - 1 ) + \rho ( v _ { r } ^ { 2 } + v _ { \theta } ^ { 2 } + v _ { \phi } ^ { 2 } ) / 2$ is the total energy density. 

The flux vectors $\vec { F } _ { r }$ , $\vec { F } _ { \theta }$ and $\vec { F } _ { \phi }$ in Eq. (5.82) are 

$$
\vec {F} _ {r} = \left[ \begin{array}{c} \rho v _ {r} \\ \rho v _ {r} ^ {2} + P \\ \rho v _ {r} v _ {\theta} \\ \rho v _ {r} v _ {\phi} \\ h v _ {r} \end{array} \right], \qquad \vec {F} _ {\theta} = \left[ \begin{array}{c} \rho v _ {\theta} \\ \rho v _ {r} v _ {\theta} \\ \rho v _ {\theta} ^ {2} + P \\ \rho v _ {\theta} v _ {\phi} \\ h v _ {\theta} \end{array} \right], \qquad \vec {F} _ {\phi} = \left[ \begin{array}{c} \rho v _ {\phi} \\ \rho v _ {r} v _ {\phi} \\ \rho v _ {\theta} v _ {\phi} \\ \rho v _ {\phi} ^ {2} + P \\ h v _ {\phi} \end{array} \right], \qquad (5. 1 0 0)
$$

where $h = \varepsilon + P$ is the enthalpy density. The pseudo force terms in the fluid momentum equations are result of non-zero time derivative of rotating unit vectors. Substitute the following time derivative relations of rotating unit vectors in spherical coordinates into Eq. (5.85), 

$$
\dot {\hat {r}} = \dot {\theta} \hat {\theta} + \dot {\phi} \sin \theta \hat {\phi}, \tag {5.101}
$$

$$
\dot {\hat {\theta}} = - \dot {\theta} \hat {r} + \dot {\phi} \cos \theta \hat {\phi}, \tag {5.102}
$$

$$
\dot {\hat {\phi}} = - \dot {\phi} \sin \theta \hat {r} - \dot {\phi} \cos \theta \hat {\theta}, \tag {5.103}
$$

such that the time derivatives of momentum vectors in spherical coordinates are 

$$
\partial_ {t} \left(\rho v _ {r} \hat {r}\right) = \hat {r} \partial_ {t} \left(\rho v _ {r}\right) + \frac {\rho v _ {r} v _ {\theta}}{r} \hat {\theta} + \frac {\rho v _ {r} v _ {\phi}}{r} \hat {\phi}, \tag {5.104}
$$

$$
\partial_ {t} \left(\rho v _ {\theta} \hat {\theta}\right) = \hat {\theta} \partial_ {t} \left(\rho v _ {\theta}\right) - \frac {\rho v _ {\theta} ^ {2}}{r} \hat {r} + \frac {\rho v _ {\theta} v _ {\phi} \cos \theta}{r \sin \theta} \hat {\phi}, \tag {5.105}
$$

$$
\partial_ {t} \left(\rho v _ {\phi} \hat {\phi}\right) = \hat {\phi} \partial_ {t} \left(\rho v _ {\phi}\right) - \frac {\rho v _ {\phi} ^ {2}}{r} \hat {r} - \frac {\rho v _ {\phi} ^ {2} \cos \theta}{r \sin \theta} \hat {\theta}, \tag {5.106}
$$

where the angular velocities are $v _ { \theta } = r \dot { \theta }$ and $v _ { \phi } = r \sin \theta \dot { \phi }$ . Extra terms in Eqs. (5.104)–(5.106) are caused by non-zero $\rho v _ { r } { \boldsymbol { \dot { r } } }$ , $\rho v _ { \theta } \dot { \hat { \theta } }$ and $\rho v _ { \phi } \dot { \hat { \phi } }$ . When the pressure tensor is strongly anisotropic in space, the matrix elements for $\{ \hat { P } \} _ { i j }$ are arranged in the same manner as that for $\{ \vec { v } \otimes \vec { v } \} _ { i j }$ in the flux vectors $\vec { F } _ { r }$ , $\vec { F _ { \theta } }$ and $\vec { F } _ { \phi }$ in Eq. (5.100). When the pressure tensor is isotropic in space $\{ \hat { P } \} _ { i j } = P \delta _ { i j }$ , the 

component-wise of the gradient of pressure can be rewritten as, 

$$
\vec {\nabla} P = \left(\frac {1}{r ^ {2}} \frac {\partial (r ^ {2} P)}{\partial r} - \frac {2 P}{r}\right) \hat {r} + \left(\frac {1}{r \sin \theta} \frac {\partial (\sin \theta P)}{\partial \theta} - \frac {P \cos \theta}{r \sin \theta}\right) \hat {\theta} + \left(\frac {1}{r \sin \theta} \frac {\partial P}{\partial \phi}\right) \hat {\phi},
$$

which introduces two components in the radial and polar directions in the source vector in Eq. (5.99). 

The technique of mapping to computational domains to derive the movingmesh numerical fluxes in curvilinear coordinates is non-trivial because of geometric factors in Eq. (5.98). In a more straightforward approach as described below, Leibniz integral rule is generalized to solve for the cell-interface numerical fluxes for any arbitrary unstructured moving mesh. Integrate the Euler equations in Eq. (5.98) over a cell volume $d V = r ^ { 2 } \sin \theta d \theta d \phi d r$ in the physical domain $\mathcal { D } =$ $[ r _ { i - 1 / 2 } , r _ { i + 1 / 2 } ] \times [ \theta _ { j - 1 / 2 } , \theta _ { j + 1 / 2 } ] \times [ \phi _ { k - 1 / 2 } , \phi _ { k + 1 / 2 } ]$ . The volume integral over the divergence term ${ \vec { \nabla } } \cdot { \vec { F } }$ for any flux vector $\vec { F }$ is represented by an equivalent surface integral by Gauss law for all types of geometries such that 

$$
\int_ {\mathcal {D}} d V \left(\frac {1}{r ^ {2}} \frac {\partial \left(r ^ {2} \vec {F} _ {r}\right)}{\partial r}\right) = \left[ A _ {r} \vec {F} _ {r} \right] _ {i - 1 / 2} ^ {i + 1 / 2}, \tag {5.107}
$$

$$
\int_ {\mathcal {D}} d V \left(\frac {1}{r \sin \theta} \frac {\partial (\sin \theta \vec {F} _ {\theta})}{\partial \theta}\right) = \left[ A _ {\theta} \vec {F} _ {\theta} \right] _ {j - 1 / 2} ^ {j + 1 / 2}, \tag {5.108}
$$

$$
\int_ {\mathcal {D}} d V \left(\frac {1}{r \sin \theta} \frac {\partial (\vec {F} _ {\phi})}{\partial \phi}\right) = \left[ A _ {\phi} \vec {F} _ {\phi} \right] _ {k - 1 / 2} ^ {k + 1 / 2}, \tag {5.109}
$$

where $A _ { r } = r ^ { 2 } \sin \theta d \theta d \phi$ , $A _ { \theta } = r \sin \theta d \phi d r$ and $A _ { \phi } = r d \theta d r$ are cell surface areas normal to $r$ , $\theta$ and $\phi$ directions respectively. The volume integral over the time derivative is given by Leibniz integral rule, 

$$
\int_ {\mathcal {D}} d V \left(\frac {\partial \vec {Q}}{\partial t}\right) = \frac {d}{d t} \int_ {\mathcal {D}} d V \vec {Q} - \int_ {\partial \mathcal {D}} \vec {Q} \vec {v} _ {\partial \mathcal {D}} \cdot d \vec {A}, \tag {5.110}
$$

where $\vec { v } _ { \partial D }$ is the cell interface velocity for any non-static mesh. The surface integral 


DEC2D and DEC3D



Cartesian moving mesh



DEC2D and DEC3D



spherical moving mesh


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/3730f67624a1ea5c1d6aff81e140f69968843fde26d5f6d2f4dcb55813ea24eb.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/9603acee6eb057fffd3d0322fd881a22e4d62c572bec88764f2f726836ff10aa.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/a25767e1d87dc224328a2a933bfc3eb71020dcf38d5ddf89f740994c503c5877.jpg)



Figure 5.17: The implementation of moving-mesh in Cartesian and spherical geometries in $D E C \mathcal { Z } D$ and $D E C 3 D$ . Leibniz integral rule states that the total rate of change of a cell content $\begin{array} { r } { \frac { d } { d t } \int _ { \mathcal { D } ( t ) } Q ( \vec { x } , t ) d V } \end{array}$ is the sum of the rate of change of the cell content $\begin{array} { r } { \int _ { D ( t ) } \frac { \partial } { \partial t } Q ( \vec { x } , t ) d V } \end{array}$ within the cell volume $\mathcal { D } ( t )$ and the rate of change of the cell content $\begin{array} { r } { \int _ { \partial D ( t ) } Q ( \vec { x } , t ) \vec { v } _ { \partial D ( t ) } \cdot \hat { n } d S } \end{array}$ due to the moving cell boundary $\partial \mathcal { D } ( t )$ , where $\hat { n }$ is the unit vector normal to the cell surface $d S$ .


is decomposed into three components in $r$ , $\theta$ and $\phi$ directions to measure the total cell-interface numerical fluxes. In $D E C \mathcal { Z } D$ and $D E C 3 D$ , the mesh moves along the radial direction only. 

$$
\int_ {\partial \mathcal {D}} \vec {Q} \vec {v} _ {\partial \mathcal {D}} \cdot d \vec {A} = \left[ A _ {r} \dot {r} \vec {Q} \right] _ {i - 1 / 2} ^ {i + 1 / 2} \tag {5.111}
$$

The final form of finite-volume moving-mesh update for Euler equations in 3-D spherical geometry is, 

$$
\begin{array}{l} \frac {\left(V _ {i j k} \vec {Q} _ {i j k}\right) ^ {n + 1} - \left(V _ {i j k} \vec {Q} _ {i j k}\right) ^ {n}}{\triangle t} = \left[ A _ {r} ^ {n} \left(\vec {F} _ {r} ^ {*} - \dot {r} \vec {Q} ^ {*}\right) \right] _ {i - 1 / 2} ^ {i + 1 / 2} + \left[ A _ {\theta} ^ {n} \vec {F} _ {\theta} ^ {*} \right] _ {j - 1 / 2} ^ {j + 1 / 2} + \left[ A _ {\phi} ^ {n} \vec {F} _ {\phi} ^ {*} \right] _ {k - 1 / 2} ^ {k + 1 / 2} \\ + \left(V _ {i j k} \vec {S} _ {i j k}\right) ^ {n}, \tag {5.112} \\ \end{array}
$$

where $\vec { F } ^ { * }$ and $\vec { Q } ^ { * }$ are upwind solution from Riemann solvers at the cell interfaces. 

The generalization to treat 3-D moving cell interfaces for an unstructured-mesh in any curvilinear coordinates is summarized into two simple steps. In the first 

step, collect all time derivatives of unit vectors in rotating coordinates $( \hat { e } _ { 1 } , \hat { e } _ { 2 } , \hat { e } _ { 3 } )$ , such as $( \hat { e } _ { r } , \hat { e } _ { \theta } , \hat { e } _ { \phi } )$ in spherical coordinates and $( \hat { e } _ { r } , \hat { e } _ { \theta } , \hat { e } _ { z } )$ in cylindrical coordinates, as source terms in the fluid momentum equations, 

$$
\frac {\partial}{\partial t} \left[ \begin{array}{l} \rho \\ \rho v _ {1} \\ \rho v _ {2} \\ \rho v _ {3} \\ \varepsilon \end{array} \right] + \vec {\nabla} \cdot \left[ \begin{array}{l} \rho \vec {v} \\ \rho \vec {v} \otimes \vec {v} + \hat {P} \\ \vec {v} (\varepsilon + P) \end{array} \right] = \left[ \begin{array}{c} 0 \\ - \left(\rho v _ {2} \dot {\hat {e}} _ {2} + \rho v _ {3} \dot {\hat {e}} _ {3}\right) \cdot \hat {e} _ {1} \\ - \left(\rho v _ {3} \dot {\hat {e}} _ {3} + \rho v _ {1} \dot {\hat {e}} _ {1}\right) \cdot \hat {e} _ {2} \\ - \left(\rho v _ {1} \dot {\hat {e}} _ {1} + \rho v _ {2} \dot {\hat {e}} _ {2}\right) \cdot \hat {e} _ {3} \\ 0 \end{array} \right]. \tag {5.113}
$$

The divergence of the pressure tensor $\hat { P }$ can be treated in the same manner as the fluid velocity outer product $\vec { v } \otimes \vec { v }$ or written as a pressure gradient in the same form of the divergence with the same geometrical factors. The latter introduces extra source terms as functions of pressures on the right hand side of Eq. (5.113). The dot product in the source terms in Eq. (5.113) is a result of producing other rotating unit vectors after taking the time derivative in the original base vector. 

In the second step, integrate the conservation laws in Eq. (5.113) over a cell volume $\mathcal { D } ( t )$ by (1) applying Leibniz integral rule to pull out the time derivative to account for the total moving-mesh cell-interface numerical flux across all surfaces $\partial \mathcal { D } ( t )$ with cell interface velocities $\vec { v } _ { \partial D ( t ) }$ for any structured or unstructured mesh and (2) applying Gauss law to convert the volume integral of divergence into a surface integral of physical fluxes. 

$$
\frac {d}{d t} \int_ {\mathcal {D} (t)} Q (\vec {x}, t) d V = \int_ {\mathcal {D} (t)} \frac {\partial}{\partial t} Q (\vec {x}, t) d V + \int_ {\partial \mathcal {D} (t)} Q (\vec {x}, t) \vec {v} _ {\partial \mathcal {D} (t)} \cdot \hat {n} d S. \qquad (5. 1 1 4)
$$

Figure (5.17) show a sketch to illustrate the implementation of moving-mesh in DEC2D and DEC3D for the step two. The moving-mesh cell-interface numerical fluxes obtained by Leibniz integral rule, in above equation is absorbed into the surface integral of physical fluxes obtained by Gauss law. The final form of finite-

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/ca177bd4e1226d43bb94436c723b790f7dbe7408ffbd911b390f090a934058ea.jpg)



Figure 5.18: The wave propagation diagram on the $x - t$ plane to explain the two-wave model of HLL approximate Riemann solver.


volume moving-mesh update for Euler equations in 3-D curvilinear coordinates for an unstructured mesh with $N$ cell surfaces is, 

$$
\frac {\left(V \vec {Q}\right) ^ {n + 1} - \left(V \vec {Q}\right) ^ {n}}{\triangle t} = \sum_ {i = 1} ^ {N} \left[ A _ {i} ^ {n} \left(\vec {F} _ {i} ^ {*} - v _ {\partial \mathcal {D} _ {i}} \vec {Q} ^ {*}\right) \right] _ {i - 1 / 2} ^ {i + 1 / 2}, \tag {5.115}
$$

where $\left( V \vec { Q } \right)$ is evaluated at coordinates of the cell center. 

# 5.3.2 Approximate Riemann solvers

# Approximate HLL Riemann solver

In this section, the approximate solutions of HLL and HLLC for a Riemann problem in 1-D Euler equations are derived following Toro’s work. [Tor09a] DEC2D and DEC3D do not use LeVque wave-propagation algorithm to solve Riemann problem, although it is feasible to track down upwind solutions by performing characteristic wave tracings. The reason is that LeVque’s method assumes the linearized wave approximation, which does not satisfy Rankine-Hugoniot jump conditions for entropies across the shock front when waves become highly nonlinear. An extra step called entropy fixed [Lev02] is required to remedy the solution, and the same problem applies to all linearized wave approximated Riemann solvers such as Roe’s. The solutions of HLL and HLLC, however, are valid for nonlinear 

waves. Consider a 1-D Riemann problem that contains two initial states $\vec { q _ { L } }$ and $\vec { q } _ { R }$ , which are constant within the cells, separated by a cell interface at the position of $x _ { 0 } = 0$ at an initial time $t _ { 0 } = 0$ . As shown in Fig. (5.18), a set of linear and nonlinear waves are initiated at the cell interface and propagate to a maximum range of influence, after time $t _ { T } = T$ , on the left at position $T S _ { L }$ and the right at $T S _ { R }$ , where $S _ { L }$ and $S _ { R }$ are the slowest and fastest wave signals estimated by Davis. [Dav88, Tor09a] 

$$
S _ {L} = \min \left[ v _ {i, L} - c _ {L}, v _ {i, R} - c _ {R} \right], \tag {5.116}
$$

$$
S _ {R} = \min \left[ v _ {i, L} + c _ {L}, v _ {i, R} + c _ {R} \right], \tag {5.117}
$$

where $c \sqrt { \gamma P / \rho }$ is the ideal gas sound speed. The estimation of signal velocities is the only place where the general equation of states matters, because the following derivation for HLL approximate Riemann solver is independent of material properties. 

$$
\frac {\partial \vec {q}}{\partial t} + \frac {\partial \vec {f _ {x}}}{\partial x} = 0, \quad \vec {q} (\vec {x} \leq x _ {0}, t _ {0}) = \vec {q} _ {L}, \quad \vec {q} (\vec {x} > x _ {0}, t _ {0}) = \vec {q} _ {R}. \tag {5.118}
$$

The explicit form for the state vector $\vec { q }$ and the physical flux ${ \vec { f } } _ { x }$ in the $x$ -direction are the same as in Eq. (5.84) in Cartesian geometry. Integrate Eq. (5.118) over the control volume defined by $[ x _ { L } , x _ { R } ] \times [ t _ { 0 } , t _ { T } ]$ , 

$$
\int_ {x _ {L}} ^ {x _ {R}} \vec {q} (x, T) d x = \int_ {x _ {L}} ^ {x _ {R}} \vec {q} (x, 0) d x + \int_ {0} ^ {T} \vec {f} _ {x} (x _ {L}, t) d t - \int_ {0} ^ {T} \vec {f} _ {x} (x _ {R}, t) d t. \tag {5.119}
$$

The first two integrals are expanded into different parts that cover the domain of influence and the domain of un-influence, 

$$
\begin{array}{l} \int_ {x _ {L}} ^ {x _ {R}} \vec {q} (x, T) d x = \int_ {x _ {L}} ^ {T S _ {L}} \vec {q} (x, T) d x + \int_ {T S _ {L}} ^ {T S _ {R}} \vec {q} (x, T) d x + \int_ {T S _ {R}} ^ {x _ {R}} \vec {q} (x, T) d x, \\ = \left(T S _ {L} - x _ {L}\right) \vec {q} _ {L} + \int_ {T S _ {L}} ^ {T S _ {R}} \vec {q} (x, T) d x + \left(x _ {R} - T S _ {L}\right) \vec {q} _ {R}, \tag {5.120} \\ \end{array}
$$

and 

$$
\int_ {x _ {L}} ^ {x _ {R}} \vec {q} (x, 0) d x = \int_ {x _ {L}} ^ {0} \vec {q} (x, 0) d x + \int_ {0} ^ {x _ {R}} \vec {q} (x, 0) d x = - x _ {L} \vec {q} _ {L} + x _ {R} \vec {q} _ {R}, \tag {5.121}
$$

whereas the third and forth integrals are 

$$
\int_ {0} ^ {T} \vec {f} _ {x} \left(x _ {L}, t\right) d t = T \vec {f} _ {x} (\vec {q} _ {L}), \tag {5.122}
$$

$$
\int_ {0} ^ {T} \vec {f} _ {x} (x _ {R}, t) d t = T \vec {f} _ {x} (\vec {q} _ {R}). \tag {5.123}
$$

Substitute the results in Eqs. (5.120)–(5.123) into Eq. (5.119), and solve for the integral $\begin{array} { r } { \int _ { T S _ { L } } ^ { T ^ { \prime } S _ { R } } \vec { q } ( x , T ) d x } \end{array}$ to obtain the approximated state vector $\vec { q } ^ { \mathrm { H L L } }$ within the range of influence defined as follows, 

$$
\vec {q} ^ {\mathrm {H L L}} = \frac {1}{T (S _ {R} - S _ {L})} \int_ {T S _ {L}} ^ {T S _ {R}} \vec {q} (x, T) d x = \frac {\vec {q} _ {R} S _ {R} - \vec {q} _ {L} S _ {L} + \vec {f} _ {L} - \vec {f} _ {R}}{S _ {R} - S _ {L}}, \tag {5.124}
$$

where $\vec { f } _ { L } = \vec { f } _ { x } ( \vec { q } _ { L } )$ and $\vec { f } _ { R } = \vec { f } _ { x } ( \vec { q } _ { R } )$ are fluxes vectors evaluated using the initial states of $\vec { q _ { L } }$ and $\vec { q } _ { R }$ respectively. Rankine-Hugoniot conditions are applied to relate the jump of HLL flux $\vec { f } ^ { \mathrm { H L L } } = \vec { f } _ { x } ( \vec { q } ^ { \mathrm { H L L } } )$ with respect to $\vec { f } _ { L / R }$ across the sharp discontinuities in the frame of moving at the slowest/fastest wave signal velocity $S _ { L / R }$ . Integrate the 1-D conservation law over space $[ \xi _ { 0 } - d \xi , \xi _ { 0 } + d \xi ]$ in Eq. (5.31) for a frame $\xi _ { 0 }$ moving at the wave signal velocity $S _ { K }$ , and take the limit of 

vanishing spatial integration, where $K = L / R$ denotes for left and right states. 

$$
\lim  _ {d \xi \rightarrow 0} \int_ {\xi_ {0} - d \xi} ^ {\xi_ {0} + d \xi} \left[ \partial_ {\tau} (x _ {\xi} \hat {q}) + \partial_ {\xi} (\hat {f} - S _ {K} \hat {q}) \right] = 0, \tag {5.125}
$$

to obtain the general jump conditions between two arbitrary states separated by a moving discontinuity $\xi _ { 0 }$ , 

$$
\vec {f} _ {K} - S _ {K} \vec {q} _ {K} = \vec {f} ^ {\mathrm {H L L}} - S _ {K} \vec {q} ^ {\mathrm {H L L}}. (5. 1 2 6)
$$

Therefore, the upwind solution for a time-averaged HLL numerical flux at a cell interface moving at a velocity $\dot { x } _ { i \pm 1 / 2 }$ is, 

$$
\begin{array}{l} \vec {f} _ {i \pm 1 / 2} ^ {*} = \vec {f} _ {L}, \quad \dot {x} _ {i \pm 1 / 2} \leq S _ {L} (5.127) \\ = \vec {f} ^ {\mathrm {H L L}}, \quad S _ {L} <   \dot {x} _ {i \pm 1 / 2} \leq S _ {R} (5.128) \\ = \vec {f} _ {R}. \quad \dot {x} _ {i \pm 1 / 2} <   S _ {R} (5.129) \\ \end{array}
$$

# Approximate HLLC Riemann solver

The two-wave model of HLL is only complete for 1-D hydrodynamics. In 3-D, however, the presence of non-zero tangential fluid velocities parallel to cell interfaces, are responsible for shear flows. The complete description for 3-D hydrodynamics is derived by Toro’s three-wave model of HLLC approximate Riemann solver [Tor09a] by adding the contact discontinuity wave. Figure (5.19) shows the main features of a wave diagram for 3-D Euler equations on the $x - t$ plane. The slowest wave signal $S _ { L }$ and the fastest wave signal $S _ { R }$ , corresponding to the propagation of nonlinear waves including the rarefaction and shock waves, are the same as HLL’s determined by Davis’ wave speed estimations. The contact discontinuity wave signal $S ^ { * }$ is traveling at the same fluid velocity normal to the cell interface. In the 1-D example along the $x$ -direction, the wave signal $S ^ { * } = v _ { x } ^ { * }$ separates two 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/1275438b2ac05ea4ccdfed35cbf6d00c0d2531e8aed51123ad4ce6d46b5b39c7.jpg)



Figure 5.19: The wave propagation diagram on the $x - t$ plane to illustrate the main features for the three-wave model of HLLC approximate Riemann solver.


unknown constant states $\vec { q } _ { L } ^ { * }$ and $\vec { q } _ { R } ^ { * }$ , with unequal mass densities $\rho _ { L } ^ { * } \neq \rho _ { R } ^ { * }$ and unequal tangential fluid velocities $v _ { y , L } ^ { * } \neq v _ { y , R } ^ { * }$ , $v _ { z , L } ^ { * } \neq v _ { z , R } ^ { * }$ in $y$ and $z$ directions respectively. The property gives rise to the shear flow along the transport of contact discontinuity wave in the $x$ -direction. However, the fluid pressure and the fluid velocity normal to the cell interface are continuous across the contact discontinuity wave $S ^ { * }$ to maintain the contact between two fluid elements or material interfaces. 

$$
P _ {L} ^ {*} = P _ {R} ^ {*} = P ^ {*}, \tag {5.130}
$$

$$
v _ {x, L} ^ {*} = v _ {x, R} ^ {*} = v ^ {*}. \tag {5.131}
$$

Since the intermediate states $\vec { q } _ { L } ^ { * }$ and $\vec { q } _ { R } ^ { * }$ have equal tangential fluid velocities with the left state $K = L$ and the right state $K = R$ , 

$$
v _ {y, K} ^ {*} = v _ {y, K}, \tag {5.132}
$$

$$
v _ {z, K} ^ {*} = v _ {z, K}, \tag {5.133}
$$

a total of four unknown variables $P ^ { * }$ , $v ^ { * }$ , $\rho _ { L } ^ { * }$ and $\rho _ { R } ^ { * }$ are introduced. Any two consecutive states are connected by Rankine-Hugoniot conditions across the jumps at wave signals $S _ { L }$ and $S _ { R }$ , 

$$
\left[ \vec {f} - S _ {K} \vec {q} \right] = 0, \tag {5.134}
$$

for the fluid mass density, 

$$
\rho_ {K} \left(v _ {x, K} - S _ {K}\right) = \rho_ {K} ^ {*} \left(v ^ {*} - S _ {K}\right), \tag {5.135}
$$

for the fluid $x$ -momentum density, 

$$
\rho_ {K} v _ {x, K} \left(v _ {x, K} - S _ {K}\right) + P _ {K} = \rho_ {K} ^ {*} v ^ {*} \left(v ^ {*} - S _ {K}\right) + P ^ {*}. \tag {5.136}
$$

Combining Eqs. (5.135)–(5.136), the pressure for both intermediate states $\vec { q } _ { L } ^ { * }$ and $\vec { q } _ { R } ^ { * }$ satisfies 

$$
P ^ {*} = P _ {K} + \rho_ {K} \left(v _ {x, K} - S _ {K}\right) \left(v _ {x, K} - v ^ {*}\right). \tag {5.137}
$$

By equating pressures $P _ { L } ^ { * } = P _ { R } ^ { * }$ using above equation, the solution for the wave signal of the contact discontinuity wave $S ^ { * } = v ^ { * }$ is 

$$
S ^ {*} = \frac {P _ {R} - P _ {L} + \rho_ {R} v _ {x , R} \left(v _ {x , R} - S _ {R}\right) - \rho_ {L} v _ {x , L} \left(v _ {x , L} - S _ {L}\right)}{\rho_ {R} \left(v _ {x , R} - S _ {R}\right) - \rho_ {L} \left(v _ {x , L} - S _ {L}\right)}. \tag {5.138}
$$

Substitute the solution of $S ^ { * } = v ^ { * }$ from Eq. (5.138) into the jump condition for the mass density and the $x$ -momentum density in Eqs. (5.135)–(5.137), solutions for the left and right intermediate mass densities $\rho _ { L } ^ { * }$ , $\rho _ { P } ^ { * }$ , $P ^ { * }$ and the intermediate state pressure $P ^ { * }$ are obtained, which also define the left and right intermediate 

state fluid internal energy densities 

$$
\varepsilon_ {K} ^ {*} = \frac {P ^ {*}}{\gamma - 1} + \frac {1}{2} \rho_ {K} ^ {*} \left(v ^ {*} v ^ {*} + v _ {y, K} ^ {2} + v _ {z, K} ^ {2}\right). \tag {5.139}
$$

The solution for HLLC intermediate states are completed, 

$$
\bar {q} _ {K} ^ {*}, ^ {\mathrm {H L L C}} = \left[ \rho , \rho v _ {x, K}, \rho v _ {y, K}, \rho v _ {z, K}, \varepsilon_ {K} \right] _ {K} ^ {*} \tag {5.140}
$$

and can be used to determine the left and right HLLC cell-interface numerical fluxes using Rankine-Hugoniot conditions in Eq. (5.134), 

$$
\vec {f} _ {K} ^ {*}, ^ {\mathrm {H L L C}} = \vec {f} _ {K} ^ {*} - S _ {K} \vec {q} _ {K} ^ {*} + S _ {K} \vec {q} _ {K} ^ {*}, ^ {\mathrm {H L L C}} \tag {5.141}
$$

Therefore, the upwind solution for a time-averaged HLLC numerical flux at a cell interface moving at a velocity $\dot { x } _ { i \pm 1 / 2 }$ is, 

$$
\begin{array}{l} \vec {f} _ {i \pm 1 / 2} ^ {*} = \vec {f} _ {L}, \quad \dot {x} _ {i \pm 1 / 2} \leq S _ {L} (5.142) \\ = \bar {f} _ {L} ^ {\mathrm {H L L C}}, \quad S _ {L} <   \dot {x} _ {i \pm 1 / 2} \leq S ^ {*} (5.143) \\ = \vec {f} _ {R} ^ {\mathrm {H L L C}}, \quad S ^ {*} <   \dot {x} _ {i \pm 1 / 2} \leq S _ {R} (5.144) \\ = \vec {f} _ {R}. \quad \dot {x} _ {i \pm 1 / 2} <   S _ {R} (5.145) \\ \end{array}
$$

# Multiple species’ internal energy advection

In DEC2D and DEC3D, the advection of multiple species’ internal energies e.g., electron pressures, alpha particle pressures and radiation pressures, is treated as a scalar advection that has a characteristic wave travels at the same fluid velocity. 

The Lagrangian forms of Euler equations are [ZR02] 

$$
D _ {t} \rho + \rho \vec {\nabla} \cdot \vec {v} = 0, \tag {5.146}
$$

$$
\rho D _ {t} \vec {v} + \vec {\nabla} P = 0, \tag {5.147}
$$

$$
D _ {t} \varepsilon + (\varepsilon + P) \vec {\nabla} \cdot \vec {v} = Q, \tag {5.148}
$$

where $D _ { t } = \partial _ { t } + \vec { v } \cdot \vec { \nabla }$ is the material time derivative and $Q$ is the rate of heat energy density for sink or source, $\varepsilon = P / ( \gamma - 1 )$ for an ideal gas. The mass density in Eq. (5.146) satisfies, 

$$
D _ {t} \ln \rho + \vec {\nabla} \cdot \vec {v} = 0, \tag {5.149}
$$

Since fluid particles $\varepsilon _ { \mathrm { f l u i d } } = \frac { 3 } { 2 } P _ { \mathrm { f l u i d } }$ and photons $\varepsilon _ { \mathrm { p h o t o n } } = 3 P _ { \mathrm { p h o t o n } }$ have different equation of states, a general form for the internal energy density $\varepsilon = 1 ^ { \prime } P$ is used to rewrite Eq. (5.148), where $\Gamma = 3 / 2$ for fluids and $\Gamma = 3$ for radiation. 

$$
D _ {t} \ln P ^ {\Gamma / (\Gamma + 1)} + \vec {\nabla} \cdot \vec {v} = \frac {Q}{P (\Gamma + 1)}. \tag {5.150}
$$

The advection of multiple species’ internal energies is operator-split from the heat source or sink terms, 

$$
D _ {t} \ln P ^ {\Gamma / (\Gamma + 1)} + \vec {\nabla} \cdot \vec {v} = 0, \tag {5.151}
$$

which is equivalent to the scalar advection for the mass density in Eq. (5.149) so that the time evolution for species internal energies is equivalent to a scalar advection of species pressure in a conservative form. 

$$
\partial_ {t} P ^ {\Gamma / (\Gamma + 1)} + \vec {\nabla} \cdot \vec {v} P ^ {\Gamma / (\Gamma + 1)} = 0, \tag {5.152}
$$

The new components in HLLC intermediate state with multiple species pressures $P _ { \mathrm { { s } } }$ , where s = e for electrons and $\mathrm { s } = \alpha$ for alpha particles, 

$$
\bar {q} _ {K} ^ {\ast , \mathrm {H L L C}} = \left[ \rho , \rho v _ {x, K}, \rho v _ {y, K}, \rho v _ {z, K}, \varepsilon_ {K}, P _ {\mathrm {e}} ^ {3 / 5}, P _ {\alpha} ^ {3 / 5}, P _ {\mathrm {p h o t o n}} ^ {3 / 4} \right] _ {K} ^ {\ast}, \tag {5.153}
$$

are determined from Rankine-Hugoniot conditions in Eq. (5.134), 

$$
P _ {\mathrm {s}, * , K} ^ {\Gamma / (\Gamma + 1)} = \left(\frac {v _ {x , K} - S _ {K}}{v ^ {*} - S _ {K}}\right) P _ {\mathrm {s}, K} ^ {\Gamma / (\Gamma + 1)}. \tag {5.154}
$$

# MUSCL second-order high resolution method

The upwind solutions of HLL and HLLC are only first-order accurate in space and time, and is too numerically diffusive to be used in practice. Although LeVeque introduced TVD flux limiters to limit the second-order Lax-Wendroff correction fluxes in Eq. (5.60)–(5.61), the assumption of linearized wave propagation in LeVeque’s approach does not satisfy Rankine-Hugoniot conditions, resulting in non-physical entropy jumps across strong shocks for highly nonlinear hydrodynamic simulations. The genuine second-order high-resolution method that boosts all types of first-order upwind solutions, was derived by van Leer in his MUSCL scheme, [vL79] which stands for Monotonic Upwind Scheme for Conservation Laws. Consider the primitive forms of Euler equations in Eqs. (5.146)–(5.148), 

$$
\partial_ {t} \rho + \vec {\nabla} \cdot (\rho \vec {v}) = 0, \tag {5.155}
$$

$$
\hat {e} _ {i} \partial_ {t} v _ {i} + (\vec {v} \cdot \vec {\nabla}) \vec {v} + \vec {\nabla} P / \rho = - v _ {i} \partial_ {t} \hat {e} _ {i}, \tag {5.156}
$$

$$
\partial_ {t} P + (\vec {v} \cdot \vec {\nabla}) P + \gamma P \vec {\nabla} \cdot \vec {v} = 0, \tag {5.157}
$$

where Einstein summation for any two repeating indices is assumed in Eq. (5.156) to denote for non-zero time derivatives of unit vectors in rotating coordinates. Introduce a state vector for primitive variables $\vec { W } = ( \rho , v _ { 1 } , v _ { 2 } , v _ { 3 } , P )$ and a state 

vector for source terms $\vec { S } = ( 0 , S _ { 1 } , S _ { 2 } , S _ { 3 } , 0 )$ such as due to pseudo forces or gravity. In spherical geometry, the primitive forms of Euler equations in Eqs. (5.155)– (5.157) can be written into 

$$
\frac {\partial \vec {W}}{\partial t} + \sum_ {i = 1} ^ {3} \hat {A} _ {i} \cdot \frac {\partial \vec {W}}{g _ {i} \partial x _ {i}} = \vec {S}, \tag {5.158}
$$

where $g _ { 1 } = 1$ , $y _ { 2 } = r$ and $g _ { 3 } = r \sin \theta$ are metric factors. The short-hand notations for coordinates are $x _ { 1 } = r$ , $x _ { 2 } = \theta$ and $x _ { 3 } = \phi$ . Consider a directional-splitting in the $i$ -direction with one-third of the original source vector, 

$$
\frac {\partial \vec {W}}{\partial t} + \hat {A} _ {i} \cdot \frac {\partial \vec {W}}{g _ {i} \partial x _ {i}} = \frac {1}{3} \vec {S}, \tag {5.159}
$$

First, an explicit predictor step with a half time-step size $\triangle t / 2$ is applied to produce an intermediate state $\tilde { W } _ { 0 } ^ { n + 1 / 2 }$ with the cell center at $x _ { 0 }$ and at the midtime level $t ^ { n + 1 / 2 }$ , 

$$
\tilde {W} _ {0} ^ {n + 1 / 2} = \vec {W} _ {0} ^ {n} - \frac {\triangle t}{2 g _ {i} \triangle x _ {i}} \hat {A} _ {i} ^ {n} \cdot \delta \vec {W} _ {0} ^ {n} + \frac {\triangle t}{6} \vec {S} _ {0} ^ {n}. \tag {5.160}
$$

Second, a corrector step with a half grid-size is applied to produce a linear slope for a cell centered at $x _ { 0 }$ with edges $\boldsymbol { x } _ { 0 \pm 1 / 2 }$ , 

$$
\bar {W} _ {0 + 1 / 2, L} ^ {n + 1 / 2} = \tilde {W} _ {0} ^ {n + 1 / 2} + \frac {1}{2} \delta \vec {W} _ {0} ^ {n}, \tag {5.161}
$$

$$
\bar {W} _ {0 - 1 / 2, R} ^ {n + 1 / 2} = \tilde {W} _ {0} ^ {n + 1 / 2} - \frac {1}{2} \delta \vec {W} _ {0} ^ {n}. \tag {5.162}
$$

The classical predictor-corrector steps in Eqs. (5.160)–(5.162) boost upwind solutions to second-order in time and space but also trigger the propagation of numerical noises whenever discontinuous solutions occur such as shocks and contact discontinuities. The last step of MUSCL scheme is to replace the slope difference $\delta \vec { W }$ with the TVD limited-slope δW~ TVD,(vanLeer) in Eqs. (5.62)–(5.63), to update 


First order HLLC


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/0caf24bd1e253019925bed9e7e1b05301f9d1137e8951ad7000cc88815fe16dc.jpg)



Second order HLLC


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/6221e86829775ba4f998d2f6c2e0afa58ef594f1283bda5cbc30923abfd719fd.jpg)



Figure 5.20: Good agreements were obtained in the benchmark tests between first and second order HLLC approximate Riemann solvers and the exact solutions in a shock tube problem. The first-order upwind HLLC solution is numerically diffusive, and is applied to update any cell with non-monotonic profile to damp the numerical noises across .


any cell with a noisy background at first-order diffusive upwind schemes. MUSCL is implemented without directional splitting in DEC2D and DEC3D. 

Figure (5.20) shows the good agreement of benchmark tests between first-order upwind, second-order MUSCL for HLLC approximate Riemann solvers and the exact solution in a shock tube problem. The diffusive first-order upwind solution is applied to update any cell with non-monotone profiles to damp the growth of numerical noises. 

Figure (5.21) shows the performance of a MUSCL scheme for HLLC approximate Riemann solver for a 2-D mode $\ell = 2 0$ . The sharp mass density profiles across RT spikes introduce significant numerical noises in MacCormack scheme even with the application of artificial numerical viscosities. However, numerical noises are damped effectively by switching to the first-order diffusive upwind solution update using the slope limiter described by Eq. (5.62)–(5.63). 

Figure (5.22) shows the performance of second-order MUSCL scheme for HLLC approximate Riemann solver for a highly nonlinear 2-D mode $\ell = 2 0$ simulation. 


Mode 20 and 5 % initial velocity perturbation


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/3c7ad121897e1016e13cf804eaffbc9deec9c6eca8ec3a74987db3b73e380547.jpg)



Figure 5.21: Comparison of numerical noise damping capabilities on a mass density profile between MUSCL scheme for HLLC and MacCormack scheme with artificial numerical viscosities for a mode $\ell = 2 0$ for a NIF implosion simulation.



Velocity at stagnation


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/5bf5c84f68c947a65ae672a9a8c84f3c56a4470e475d6db43bb8badc48eddb31.jpg)



Density at stagnation


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/bd4a3c6a5295a8940b8668ea9f4c6302cb4d93449b86583031ed2f61dc7619a6.jpg)



Figure 5.22: Performance of second-order MUSCL scheme for HLLC approximate Riemann solver for a 2-D mode $\ell = 2 0$ on resolving the fluid velocity on the left and the mass density on the right.


Highly well resolved RT mushroom structures and numerically stable solutions are obtained. 

# PPM third-order high resolution method

PPM solves the component-wise characteristic wave equations by using a parabolic cell profile to define the upwind solution. Following the same notation as LeVeque in Eqs. (5.44)–(5.47) to diagonalize the matrix $\hat { A } _ { i } = \hat { R } _ { i } \cdot \hat { \Lambda } _ { i } \cdot \hat { R } _ { i } ^ { - 1 }$ in Eq. (5.158) using the right-matrix $\hat { R } _ { i }$ in different directions $i = 1 , 2 , 3$ , whereas the left-matrix is defined as $\hat { L } _ { i } = \hat { R } _ { i } ^ { - 1 }$ . The eigenvalue of the matrix ${ \hat { A } } _ { i }$ is a set of characteristic wave velocities $\{ v _ { i } \pm c , v _ { i } \}$ for 3-D Euler equations, where $c = \sqrt { \gamma P / \rho }$ is the fluid sound speed. $v _ { i } \pm c$ represents the nonlinear wave such as rarefaction and shock wave, whereas $v _ { i }$ with the multiplicity of $D + m$ represents the linear wave of contact discontinuity. Here $D$ is the dimension of the single-fluid Euler equations and $m = 3$ is the total number of additional species pressure advection equations for electrons, alpha particles and photons added to the single-fluid Euler equations. Multiply the both sides of Eq. (5.158) with the left-matrix $\hat { L } _ { i }$ to map the primitive state vector $\vec { W } _ { i }$ in the $i$ -direction onto the characteristic state vector $\vec { C } _ { i } = \hat { L } _ { i } \cdot \vec { W } _ { i }$ , 

$$
\frac {\partial \vec {C}}{\partial t} + \sum_ {i = 1} ^ {3} \hat {\Lambda} _ {i} \cdot \frac {\partial \vec {C}}{g _ {i} \partial x _ {i}} = \hat {L} _ {i} \cdot \vec {S}. \tag {5.163}
$$

By performing directional-splitting, with each direction shares one-third of the original source term $\hat { L } _ { i } \cdot \vec { S }$ , and consider the $j$ -component of the 1-D characteristic wave equation along the $i$ -direction. 

$$
\frac {\partial C _ {j}}{\partial t} + \lambda_ {j i} \frac {\partial C _ {j}}{g _ {i} \partial x _ {i}} = \frac {1}{3} \hat {L} _ {i} \cdot S _ {j}, \tag {5.164}
$$

where $\lambda _ { j i }$ is the characteristic wave velocity and $g _ { 1 } = 1$ , $g _ { 2 } = r$ and $g _ { 3 } = r \sin \theta$ are metric factors. The upwind solution for Eq. (5.164) is first solved by LeVeque 

characteristic wave tracing method, 

$$
\langle C _ {j} (x _ {k}) \rangle^ {n + 1} = \langle C _ {j} (x _ {k}) \rangle^ {n} - \frac {\triangle t}{g _ {i} \triangle x _ {i}} \left[ F _ {j i} ^ {\downarrow} (x _ {k + 1 / 2}) - F _ {j i} ^ {\downarrow} (x _ {k - 1 / 2}) - \frac {1}{3} \hat {L} _ {i} \cdot S _ {j} ^ {n} (x _ {k}) \right],
$$

where the upwind cell interface numerical fluxes depend on the sign of characteristic wave velocities $\lambda _ { j i } ^ { \pm }$ computed within the cell centered at $x _ { k }$ at the time level $t ^ { n }$ , in which $^ +$ stands for right-traveling waves and - for left-traveling waves. 

$$
F _ {j i} ^ {\downarrow} \left(x _ {k + 1 / 2}\right) = \lambda_ {j i} ^ {+} \left\langle C _ {j} \left(x _ {k}\right) \right\rangle^ {n} + \lambda_ {j i} ^ {-} \left\langle C _ {j} \left(x _ {k + 1}\right) \right\rangle^ {n}, \tag {5.165}
$$

$$
F _ {j i} ^ {\downarrow} \left(x _ {k - 1 / 2}\right) = \lambda_ {j i} ^ {+} \left\langle C _ {j} \left(x _ {k - 1}\right) \right\rangle^ {n} + \lambda_ {j i} ^ {-} \left\langle C _ {j} \left(x _ {k}\right) \right\rangle^ {n}. \tag {5.166}
$$

In the first step of PPM, all cell-averaged contents $\langle C _ { j } ( x _ { k } ) \rangle ^ { n }$ on the right-hand-side of Eqs. (5.165)–(5.166) are reconstructed, 

$$
\langle C _ {j} \left(x _ {k}\right) \rangle^ {n} = \frac {1}{\triangle x _ {i}} \int_ {x _ {i, k - 1 / 2}} ^ {x _ {i, k + 1 / 2}} C _ {j} \left(x _ {i}\right) d x _ {i} \tag {5.167}
$$

using a parabolic profile, 

$$
C \left(x _ {i}\right) = C _ {L, k} + \sigma \left(\triangle C _ {k} + C _ {6, k} (1 - \sigma)\right), \tag {5.168}
$$

$$
\sigma = \frac {x _ {i} - x _ {i , k - 1 / 2}}{\triangle x _ {i}}, \tag {5.169}
$$

where the variables of $\triangle C _ { k }$ and $C _ { 6 , k }$ are expressed in terms of two unknowns of cell edge values ${ C } _ { L / R , k }$ , 

$$
\triangle C _ {k} = C _ {R, k} - C _ {L, k}, \tag {5.170}
$$

$$
C _ {6, k} = 6 \left(\left\langle C _ {j} \left(x _ {k}\right) \right\rangle^ {n} - \left(C _ {R, k} + C _ {L, k}\right) / 2\right). \tag {5.171}
$$

The formula of parabola in Eq. (5.168) has a property that the total cell content is conserved after the high order parabola cell reconstruction, and is con-

tinuous across any two adjacent cells $C ( x _ { i , k - 1 / 2 } ) = C _ { L , k }$ and $C ( x _ { i , k + 1 / 2 } ) = C _ { R , k }$ . In the second-step of PPM, the two unknowns of left and right cell edge values are assigned with high order finite-difference representations $C _ { L , k } = C _ { k - 1 / 2 }$ and $C _ { R , k } = C _ { k + 1 / 2 }$ given by, 

$$
C _ {k + 1 / 2} = \frac {7}{1 2} \left(\langle C _ {j} (x _ {k}) \rangle^ {n} + \langle C _ {j} (x _ {k + 1}) \rangle^ {n}\right) - \frac {1}{1 2} \left(\langle C _ {j} (x _ {k + 2}) \rangle^ {n} + \langle C _ {j} (x _ {k - 1}) \rangle^ {n}\right),
$$

The following monotone limited-slope [SGT $^ +$ 08] is applied for more numerically stable results, 

$$
C _ {k + 1 / 2} = \frac {\langle C _ {j} (x _ {k}) \rangle^ {n} + \langle C _ {j} (x _ {k + 1}) \rangle^ {n}}{2} - \frac {\delta \langle C _ {j} (x _ {k}) \rangle_ {\mathrm {T V D}} ^ {n} + \delta \langle C _ {j} (x _ {k + 1}) \rangle_ {\mathrm {T V D}} ^ {n}}{6}. \tag {5.172}
$$

In the last step of PPM, the left and right cell edge values are reset using the following two conditions to ensure the interpolated parabola profile is monotone. The first one is the same noise detection as in MUSCL, meaning any cell exhibits a wiggle profile in the neighborhood of adjacent cells requires no cell reconstruction and is solved at the first-order diffusive upwind scheme, 

(5.173) 

The second one is to avoid overshooting the parabola above or below the left and right cell edge values; otherwise, a local extrema within the cell is produced, which happens when the assignment of two edge values are too closed to each other $C _ { L , k } \sim C _ { R , k }$ but the curvature of the parabola is large. The condition of 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/545a3ec1122e3f709fbc6de53662ae26988281eba9fbdb87f31f4f1d7ed3f15b.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/c12f837fb92012e0d6e3551b7adbf341ce6250812c843774274fc4591878abdf.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/b9f28157d110eb7fe5489728b0f33ee0277ad4f9b4d85ea7ff21c4301a91526c.jpg)



Figure 5.23: Benchmark test of PPM high resolution method to boost HLLC upwind solution to the third order in the spherical-mesh version DEC3D with LILAC’s result in a shock tube problem.


overshooting occurs when $\lvert \triangle C _ { k } \rvert < \lvert C _ { 6 , k } \rvert$ , 

$$
(2) \text {I f}: | C _ {R / L, k} - \langle C _ {j} (x _ {k}) \rangle^ {n} | \geq 2 | C _ {L / R, k} - \langle C _ {j} (x _ {k}) \rangle^ {n} |, \tag {5.174}
$$

$$
\mathrm {T h e n :} C _ {R / L, k} = \langle C _ {j} (x _ {k}) \rangle^ {n} - 2 \left(C _ {L / R, k} - \langle C _ {j} (x _ {k}) \rangle^ {n}\right). \tag {5.175}
$$

Figure (5.23) shows a good agreement in the benchmark test of PPM high resolution method to boost HLLC solution to the third order in the sphericalmesh DEC3D with LILAC in a shock tube problem. 

# 5.4 Radiation transport

In this section, a comprehensive literature review on the radiation transport theory is presented. The objective is to derive the essential governing equation for the multi-group flux-limited radiation diffusion approximation implemented in DEC2D [WED+14] and DEC3D [WBB+15]. 

# 5.4.1 Multi-group flux-limited radiation diffusion

Boltzmann transport equation [Bru02] is a hyperbolic partial differential equation (PDE) that 

$$
\frac {1}{v} \frac {\partial f}{\partial t} + \vec {\Omega} \cdot \vec {\nabla} f = C (f). \tag {5.176}
$$

describes the transport of any type of particles through a material at a finite particle speed $\boldsymbol { v }$ . The phase space density $f \left( \vec { x } , \nu , \vec { \Omega } , t \right)$ is a seven dimensional function: three space, two angles, one energy and one time. $C$ is the collision operator that represents all interactions with material. The loss of hyperbolic PDE nature in the radiation diffusion approximation [Bru02] incurs the flux limiter technique [Pom73] to avoid the numerical flux in the diffusion equation being larger than the actual particle flux $f v$ . 

Boltzmann transport equation is applied to describe the time evolution of a spectral radiation intensity $I _ { \nu } \left( \vec { x } , \nu , \vec { \Omega } , t \right)$ , which represents the radiant energy in a spectral interval $d \nu$ , passing per unit time through a unit area, with the direction of radiant energy propagation contained within an element of solid angle $d \Omega$ along the unit vector $\vec { \Omega }$ . [ ZR02] 

$$
\frac {1}{c} \frac {\partial I _ {\nu}}{\partial t} + \vec {\Omega} \cdot \vec {\nabla} I _ {\nu} = \varepsilon_ {\nu} - \kappa_ {\nu} I _ {\nu}. (5. 1 7 7)
$$

The left hand side of Eq. (5.177) describes the transport of the spectral radiation intensity $I _ { \nu }$ at the light speed $c$ , while the right hand side terms describe the interaction between radiation and matter in terms of the emissivity $\varepsilon _ { \nu } \left( \vec { x } , \nu , t \right)$ , which is the sum of spontaneous and stimulated emissions, and the opacity of matter $\kappa _ { \nu } \left( \vec { x } , \nu , t \right)$ , which is the sum of scattering and true absorptions. In general, emissivity and opacity are functions of angles due to relativistic Doppler effects. However, the angular dependence of $\varepsilon _ { \nu }$ and $\kappa _ { \nu }$ is lost for static medium, since 

emitting and absorbing photons in the atomic level take place without preferential direction. Because the scattering mean-free-path for a photon is much longer than the absorption mean-free-path, it is reasonable to ignore the scattering effects between free electrons and photons such as Thompson, Compton and Raman scatterings in Eq. (5.177) in the deceleration phase of ICF implosions. 

The dynamics of radiation transport in Eq. (5.177) is the following. Emissivity $\varepsilon _ { \nu }$ is the source term that corresponds to emit photons. The absorption term $- \kappa _ { \nu } I _ { \nu }$ with a minus sign indicates the attenuation of spectra radiation intensity as the light wave or photons travel through the medium. The left hand side of Eq. (5.177) is a differential operator that transports a photon over an infinitesimal path length ds at the light speed $c$ in the direction of unit vector $\vec { \Omega }$ , 

$$
\frac {d I _ {\nu}}{d s} = \left[ \frac {1}{c} \frac {\partial}{\partial t} + \vec {\Omega} \cdot \vec {\nabla} \right] I _ {\nu}. (5. 1 7 8)
$$

In Cartesian coordinates, the streaming term [Bru02] 

$$
\vec {\Omega} \cdot \vec {\nabla} I _ {\nu} = \left[ \sin \theta \cos \phi \frac {\partial}{\partial x} + \sin \theta \sin \phi \frac {\partial}{\partial y} + \cos \theta \frac {\partial}{\partial z} \right] I _ {\nu} \tag {5.179}
$$

measures the incident spectral radiation flux along the direction of unit vector $\vec { \Omega }$ . Equation (5.177) is integrable along a straight line trajectory 

$$
\frac {d I _ {\nu}}{d s} = \varepsilon_ {\nu} - \kappa_ {\nu} I _ {\nu} \tag {5.180}
$$

by introducing variables of the infinitesimal optical path $d \tau _ { \nu } = \kappa _ { \nu } d s$ and the source function $S _ { \nu } = \varepsilon _ { \nu } / \kappa _ { \nu }$ such that the straight line radiation transport 

$$
\frac {d I _ {\nu}}{d \tau_ {\nu}} = S _ {\nu} - I _ {\nu} \tag {5.181}
$$

has a simple formal solution 

$$
I _ {\nu} (\tau) = I _ {\nu} (0) e ^ {- \tau_ {\nu}} + \int_ {0} ^ {\tau_ {\nu}} d \tau_ {\nu} ^ {\prime} S _ {\nu} \left(\tau_ {\nu} ^ {\prime}\right) e ^ {- \tau_ {\nu} + \tau_ {\nu} ^ {\prime}}, \tag {5.182}
$$

where $\begin{array} { r } { \tau _ { \nu } = \int _ { s _ { 0 } } ^ { s } \kappa _ { \nu } d s } \end{array}$ is the optical path, which is dimensionless because the opacity $\kappa _ { \nu }$ has the dimension of length $^ - { 1 }$ in Eq. (5.177). The optical thickness is determined by the fraction of spectral radiation intensity $I _ { \nu } ( \tau ) / I _ { \nu } ( 0 )$ passing through a material from initial position $s _ { 0 }$ to the final point $s$ . [Sch65] A optically thick plasma is defined by $\tau _ { \nu } > 1$ , and a optically thin plasma by $\tau _ { \nu } ~ < ~ 1$ . The geometrical interpretation for $\tau _ { \nu } = 1$ means that a photon travels on average of one diffusion mean-free-path $1 / \kappa _ { \nu }$ before absorption. A realistic radiation transport must seek for quantum mechanical treatments [Dir01, MS10] to quantize the electromagnetic field interacting with charged particles through Coulomb potential. During the radiation-matter interaction, a photon is born and is emitted isotropically in space, followed by an annihilation shortly. The averaged probability $p ( \tau _ { \nu } , \tau _ { \nu } + d \tau _ { \nu } )$ for a photon to travel from the point $s$ in space to another point $s + d s$ before absorption is determined by the attenuation of the radiation intensity $d I _ { \nu } / I _ { \nu } = d \tau _ { \nu } = \kappa _ { \nu } ^ { \mathrm { a } } d s$ . Otherwise, the photon is scattered isotropically in angles with the probability $1 - p ( \tau _ { \nu } , \tau _ { \nu } + d \tau _ { \nu } )$ at the point $s + d s$ . 

When the real part of index of refraction $\mu _ { \nu }$ varies strongly in space, the ray tracing for optical paths is needed. From Snell’s law, the quantity $I _ { \nu } / \mu _ { \nu } ^ { 2 }$ is constant along every ray path, provided there is only negligible energy losses or gains due to refraction, absorption, or emission. [CG68] Otherwise, 

$$
\frac {d}{d s} \left(\frac {I _ {\nu}}{\mu_ {\nu} ^ {2}}\right) = \frac {\varepsilon_ {\nu} - \kappa_ {\nu} I _ {\nu}}{\mu_ {\nu} ^ {2}}, \tag {5.183}
$$

and the material derivative is replaced with 

$$
\frac {d}{d s} = \frac {1}{c} \frac {\partial}{\partial t} + \vec {\Omega} \cdot \vec {\nabla} + \frac {d \vec {\Omega}}{d s} \cdot \vec {\nabla} _ {\vec {\Omega}}. (5. 1 8 4)
$$

The last term $\frac { d \vec { \Omega } } { d s } \cdot \vec { \nabla } _ { \vec { \Omega } }$ in Eq. (5.184) measures the additional net change of the spectral radiation intensity due to refractions of neighboring light rays. Born and Wolf showed that the rate of change [BW80] 

$$
\frac {d \vec {\Omega}}{d s} = \frac {1}{\mu_ {\nu}} \left[ \vec {\nabla} \mu_ {\nu} - \vec {\Omega} \left(\vec {\Omega} \cdot \vec {\nabla} \mu_ {\nu}\right) \right] (5. 1 8 5)
$$

is zero whenever $\vec { \Omega }$ is parallel with $\pm \vec { \nabla } \mu _ { \nu }$ or $\vec { \nabla } \mu _ { \nu } = 0$ . The latter means for the constant index of refraction in space. Equations (5.183)-(5.185) can be applied to model the laser energy deposition onto the plasma, near the corona region at the quarter of critical density, by using Kramers free-free opacity $\kappa _ { \nu } ^ { \mathrm { f f } }$ [ZR02] to account for the inverse bremsstrahlung radiation absorption. 

$$
\kappa_ {\nu} ^ {\mathrm {f f}} = \frac {4}{3} \sqrt {\frac {2 \pi}{3 m _ {\mathrm {e}} k _ {\mathrm {B}} T _ {\mathrm {e}}}} \frac {n _ {\mathrm {i}} n _ {\mathrm {e}} Z ^ {2} e ^ {6}}{h c m _ {\mathrm {e}} \nu^ {3}} \left[ \frac {1}{\mathrm {c m}} \right], \qquad (5. 1 8 6)
$$

where $n _ { \mathrm { i } }$ and $n _ { \mathrm { e } }$ are ion and electron number densities. Quantities in Eq. (5.186) are in C.G.S. units. 

The interaction of radiation with matter is an exchange of energy between the radiation field and the energy levels of molecules and atoms, which are defined by the Boltzmann temperatures. For a plasma in local thermal equilibrium (LTE) [CG68], the Boltzmann temperature is in equilibrium with the kinetic temperature Ti=1,2,3,... for all energy levels, but not in equilibrium with the radiation field. In a complete thermodynamic equilibrium (TE) [CG68], the Boltzmann temperature is in equilibrium with the kinetic temperature for all energy levels, as well as with the radiation field. The Boltzmann temperature is the electron temperature $T _ { \mathrm { e } }$ because the total number of states for a system containing ionized atoms or 

molecules embedded within a huge free electron sea is dominated by the statistics of electrons. [Coo66] 

A TE plasma is characterized by the mean free path for absorption of radiation to be less than the dimension of the plasma and a black-body radiation temperature $T _ { \mathrm { r } }$ . A full detailed balance [Coo66] is required for radiative and collisional processes to establish a TE plasma $T _ { i = 1 , 2 , 3 , \dots } = T _ { \mathrm { e } } = T _ { \mathrm { r } }$ . However, the escape of radiation from plasma disturbs the full detailed balance, resulting in deviations from TE. 

For a collision-dominated plasma, LTE is established, which is defined by Griem’s criterion [Gri62, Coo66] requiring collisional rates being much larger than radiative transition rates, meaning that a LTE plasma satisfies $T _ { i = 1 , 2 , 3 , \dots } = T _ { \mathrm { e } }$ and the black-body radiation temperature $T _ { \mathrm { r } }$ does not exist. The assumption of LTE in multi-group radiation transport allows the existence of multiple radiation temperatures, defined by $T _ { \mathrm { r } } ^ { ( \nu ) } = h \nu / k _ { \mathrm { B } }$ , varying from low to high energy group photons. 

For a non-LTE plasma $T _ { i = 1 , 2 , 3 , \dots } \neq T _ { \mathrm { e } }$ and the black-body temperature $T _ { \mathrm { r } }$ does not exist, kinetic temperatures do not in thermal equilibrium with the electron temperature such as due to imbalanced bound-bound and bound-free transition rates. The modelings of non-LTE radiation transport require solving for different kinetic temperatures $T _ { i }$ =1,2,3,... from a set of rate equations for all possible in-line transitions. 

In LTE, emissivity and opacity are related by Kirchoff’s law 

$$
\varepsilon_ {\nu} = \kappa_ {\nu} B _ {\nu} \tag {5.187}
$$

through a universal function of frequency and material temperature $T _ { \mathrm { e } }$ , the Planck 

spectrum $B _ { \nu } ( \nu , T _ { \mathrm { e } } )$ 

$$
B _ {\nu} (\nu , T _ {\mathrm {e}}) = \frac {2 h \nu^ {3}}{c ^ {2}} \frac {1}{\exp \left(\frac {h \nu}{K _ {\mathrm {B}} T _ {\mathrm {e}}}\right) - 1}, \tag {5.188}
$$

where $h$ and $K _ { \mathrm { B } }$ are Planck and Boltzmann constants respectively. The $- 1$ in the denominator of Eq. (5.188) is due to the summation of non-interacting partition functions for a boson gas or photons, which is a remarkable feature that causes Bose-Einstein condensation for boson particles in low-temperature condensed matter. Kirchoff’s law simplifies the radiation transport equation, being not to depend on the input of tabular emissivity data, 

$$
\frac {1}{c} \frac {\partial I _ {\nu}}{\partial t} + \vec {\Omega} \cdot \vec {\nabla} I _ {\nu} = \kappa_ {\nu} (B _ {\nu} - I _ {\nu}) \tag {5.189}
$$

Similar to fluid descriptions of plasma, macroscopic variables are obtained from taking various moment of Boltzmann equation. The spectral radiant energy density $U _ { \nu } ( \vec { x } , \nu , t )$ , the spectral radiant energy flux vector $\vec { F } _ { \nu } ( \vec { x } , \nu , t )$ and the spectral radiant pressure tensor $\hat { P } _ { \nu } ( \vec { x } , \nu , t )$ are defined as the zeroth, the first and the second angular moment with respect to the spectral radiation intensity $I _ { \nu } ( \vec { x } , \nu , \vec { \Omega } , t )$ respectively. This is a key step to remove the angular dependence. Because multigroup radiation diffusions solve for the fifth-dimensional variable $U _ { \nu }$ , but not the seven-dimensional variable $I _ { \nu }$ , which is too computational expensive. 

$$
U _ {\nu} (\vec {x}, \nu , t) = \frac {1}{c} \int_ {4 \pi} I _ {\nu} (\vec {x}, \nu , \vec {\Omega}, t) d \Omega , \tag {5.190}
$$

$$
\vec {F} _ {\nu} (\vec {x}, \nu , t) = \int_ {4 \pi} I _ {\nu} (\vec {x}, \nu , \vec {\Omega}, t) \vec {\Omega} d \Omega , \tag {5.191}
$$

$$
\hat {P} _ {\nu} (\vec {x}, \nu , t) = \frac {1}{c} \int_ {4 \pi} I _ {\nu} (\vec {x}, \nu , \vec {\Omega}, t) \vec {\Omega} \otimes \vec {\Omega} d \Omega , \tag {5.192}
$$

where $\vec { \Omega } \otimes \vec { \Omega }$ denotes for the outer product between two vectors. In Cartesian coordinates, the unit vector of radiation propagation $\vec { \Omega } = g _ { i } \hat { e } _ { i }$ is specified by three 

geometrical factors $g _ { 1 } = \sin \theta \cos \phi$ , $g _ { 2 } = \sin \theta \sin \phi$ and $g _ { 3 } = \cos \theta$ in the $x$ -direction $\hat { e } _ { 1 }$ , the $y$ -direction ${ \hat { e } } _ { 2 }$ and the $z$ -direction ${ \hat { e } } _ { 3 }$ respectively. The matrix element of the outer product is $\{ \vec { \Omega } \otimes \vec { \Omega } \} _ { i j } = g _ { i } g _ { j }$ . Integrate the equation of radiation transport in Eq. (5.189) over all solid angles gives the zeroth moment equation 

$$
\frac {\partial U _ {\nu}}{\partial t} + \vec {\nabla} \cdot \vec {F} _ {\nu} = c \kappa_ {\nu} \left(4 \pi B _ {\nu} / c - U _ {\nu}\right). \tag {5.193}
$$

The factor of $4 \pi$ is because the emission of photons is isotropic in angles stated by Planck spectrum in Eq. (5.188). Multiply the equation of radiation transport in Eq. (5.189) with the unit vector $\vec { \Omega }$ and integrate again over all solid angles. The first moment of Planck spectrum vanishes because it is independent of angles. 

$$
\frac {1}{c} \frac {\partial \vec {F} _ {\nu}}{\partial t} + c \vec {\nabla} \cdot \hat {P} _ {\nu} = - \kappa_ {\nu} \vec {F} _ {\nu} \tag {5.194}
$$

Equations (5.193) and (5.194) are respectively, the energy and momentum equations of radiation transport. However, the spectral pressure tensor $\hat { P } _ { \nu }$ is an unknown variable, the closure problem in the high order moment analysis. 

The closure is to represent the spectral radiant pressure tensor $\hat { P } _ { \nu }$ in terms of spectral radiant energy density $U _ { \nu }$ by introducing the spherical harmonic expansion to the spectral radiation intensity $I _ { \nu }$ up to the first order, so-called Eddington or P1 approximation. [Bru02, Pom73] 

$$
I _ {\nu} (\vec {x}, \nu , \vec {\Omega}, t) = \frac {1}{4 \pi} I _ {\nu} ^ {(0)} (\vec {x}, \nu , t) + \frac {3}{4 \pi} \vec {\Omega} \cdot \vec {I} _ {\nu} ^ {(1)} (\vec {x}, \nu , t), \tag {5.195}
$$

where the scalar $I _ { \nu } ^ { ( 0 ) }$ and the vector $\vec { I } _ { \nu } ^ { ( 1 ) }$ are the zero and the first moment of the 

spectral radiation intensity $I _ { \nu }$ with respect to angles respectively. 

$$
I _ {\nu} ^ {(0)} (\vec {x}, \nu , t) = \int_ {4 \pi} I _ {\nu} (\vec {x}, \nu , \vec {\Omega}, t) d \Omega , \tag {5.196}
$$

$$
\vec {I} _ {\nu} ^ {(1)} (\vec {x}, \nu , t) = \int_ {4 \pi} I _ {\nu} (\vec {x}, \nu , \vec {\Omega}, t) \vec {\Omega} d \Omega . \tag {5.197}
$$

The closure assumes the spectral radiation intensity $I _ { \nu }$ varies weakly with respect to angles by neglecting high order anisotropic terms. Define the angular-averaged spectral radiant energy density $\langle U _ { \nu } ( \vec { x } , \nu , t ) \rangle$ and the angular-averaged spectral radiant energy flux $\langle \vec { F } _ { \nu } ( \vec { x } , \nu , t ) \rangle$ 

$$
\langle U _ {\nu} (\vec {x}, \nu , t) \rangle = \frac {1}{c} \int_ {4 \pi} I _ {\nu} (\vec {x}, \nu , \vec {\Omega}, t) d \Omega , \tag {5.198}
$$

$$
\langle \vec {F} _ {\nu} (\vec {x}, \nu , t) \rangle = \int_ {4 \pi} I _ {\nu} (\vec {x}, \nu , \vec {\Omega}, t) \vec {\Omega} d \Omega . \tag {5.199}
$$

Substitute the relations of $\begin{array} { r } { \langle U _ { \nu } \rangle = \frac { 1 } { c } I _ { \nu } ^ { ( 0 ) } } \end{array}$ and $\langle \vec { F } _ { \nu } \rangle = \vec { I } _ { \nu } ^ { ( 1 ) }$ into Eq. (5.195) 

$$
I _ {\nu} (\vec {x}, \nu , \vec {\Omega}, t) = \frac {c}{4 \pi} \langle U _ {\nu} (\vec {x}, \nu , t) \rangle + \frac {3}{4 \pi} \vec {\Omega} \cdot \langle \vec {F} _ {\nu} (\vec {x}, \nu , t) \rangle , (5. 2 0 0)
$$

to obtain P1 approximation for $U _ { \nu }$ , $\vec { F } _ { \nu }$ and $\hat { P } _ { \nu }$ in Eq. (5.190) 

$$
U _ {\nu} (\vec {x}, \nu , t) = \left\langle U _ {\nu} \right\rangle + \frac {3}{4 \pi c} \left\langle \vec {F} _ {\nu} \right\rangle \cdot \int_ {4 \pi} \vec {\Omega} d \Omega , \tag {5.201}
$$

$$
\vec {F} _ {\nu} (\vec {x}, \nu , t) = \frac {c}{4 \pi} \left\langle U _ {\nu} \right\rangle \int_ {4 \pi} \vec {\Omega} d \Omega + \frac {3}{4 \pi} \left\langle \vec {F} _ {\nu} \right\rangle \cdot \int_ {4 \pi} \vec {\Omega} \otimes \vec {\Omega} d \Omega , \tag {5.202}
$$

$$
\hat {P} _ {\nu} (\vec {x}, \nu , t) = \frac {1}{4 \pi} \langle U _ {\nu} \rangle \int_ {4 \pi} \vec {\Omega} \otimes \vec {\Omega} d \Omega + \frac {3}{4 \pi c} \langle \vec {F} _ {\nu} \rangle \cdot \int_ {4 \pi} \vec {\Omega} \otimes \vec {\Omega} \otimes \vec {\Omega} d \Omega . \tag {5.203}
$$

Using the properties of odd integrands $\begin{array} { r } { \int _ { 4 \pi } \vec { \Omega } d \Omega = \int _ { 4 \pi } \vec { \Omega } \otimes \vec { \Omega } \otimes \vec { \Omega } d \Omega = 0 } \end{array}$ and the even integrand $\begin{array} { r } { \int _ { 4 \pi } \vec { \Omega } \otimes \vec { \Omega } d \Omega \ = \ \frac { 4 \pi } { 3 } } \end{array}$ [ZR02], P1 approximation for $U _ { \nu }$ , $\vec { F } _ { \nu }$ and $\hat { P } _ { \nu }$ are $U _ { \nu } ( \vec { x } , \nu , t ) = \langle U _ { \nu } ( \vec { x } , \nu , t ) \rangle$ , $\vec { F _ { \nu } } ( \vec { x } , \nu , t ) = \langle \vec { F _ { \nu } } ( \vec { x } , \nu , t ) \rangle$ , and $\hat { P } _ { \nu } ( \vec { x } , \nu , t ) =$ $\langle U _ { \nu } ( \vec { x } , \nu , t ) \rangle / 3$ . The spectral radiant pressure tensor $\hat { P } _ { \nu } = \langle U _ { \nu } \rangle / 3$ is observed losing its dependence on angles completely. Substitute P1 approximations into the 

energy equation in Eq. (5.193) and the momentum equation in Eq. (5.194), 

$$
\frac {\partial \langle U _ {\nu} \rangle}{\partial t} + \vec {\nabla} \cdot \langle \vec {F} _ {\nu} \rangle = c \kappa_ {\nu} \left(4 \pi B _ {\nu} / c - \langle U _ {\nu} \rangle\right), \tag {5.204}
$$

$$
\frac {1}{c} \frac {\partial \langle \vec {F} _ {\nu} \rangle}{\partial t} + \frac {c}{3} \vec {\nabla} \langle U _ {\nu} \rangle = - \kappa_ {\nu} \langle \vec {F} _ {\nu} \rangle . \tag {5.205}
$$

Above closure equations are hyperbolic that constitute the fundamental transport of radiation. However, the cost to assume a weakly anisotropic spectral radiation intensity $I _ { \nu }$ leads to a reduction of characteristic speed of transport from the light speed $c$ to $c / \sqrt { 3 }$ . The shortcoming of P1 approximations is illustrated by considering an electromagnetic wave propagation in vacuum without attenuation in the limit of zero opacities $\kappa _ { \nu }  0$ in the closure equations, i.e., $\begin{array} { r } { \frac { \partial ^ { 2 } \langle U _ { \nu } \rangle } { \partial t ^ { 2 } } = \frac { c ^ { 2 } } { 3 } \nabla ^ { 2 } \langle U _ { \nu } \rangle } \end{array}$ . The second approximation is the Fick’s law by assuming the spectral radiant flux varies slowly in time compared with the spatial gradient for the spectral radiant energy density, i.e., $\begin{array} { r } { \frac { \partial \langle \vec { F } _ { \nu } \rangle } { \partial t } \ll \frac { c ^ { 2 } } { 3 } \vec { \nabla } \langle U _ { \nu } \rangle } \end{array}$ in the momentum equation of radiation transport. By dropping the time derivative in Eq. (5.205) 

$$
\langle \vec {F} _ {\nu} \rangle = - \frac {c}{3 \kappa_ {\nu}} \vec {\nabla} \langle U _ {\nu} \rangle \tag {5.206}
$$

Substitute Eq. (5.206) into Eq. (5.204), the two of radiation transport equations in the closure is reduced into a single parabolic diffusion equation 

$$
\frac {\partial \left\langle U _ {\nu} \right\rangle}{\partial t} = \vec {\nabla} \cdot D _ {\nu} \vec {\nabla} \left\langle U _ {\nu} \right\rangle + c \kappa_ {\nu} \left(4 \pi B _ {\nu} / c - \left\langle U _ {\nu} \right\rangle\right) \tag {5.207}
$$

with a local diffusion coefficient $D _ { \nu }$ with the dimension of length $^ 2$ /time 

$$
D _ {\nu} = \frac {c}{3 \kappa_ {\nu}}. \tag {5.208}
$$

In free-streaming limit for a optically thin plasma, however, the transport of photons described by a diffusion PDE can violate the causality relation, because a 

parabolic PDE by definition has an infinite characteristic speed of transport. Consequently, the shortcoming incurred by Fick’s law leads to no upper bound for the spectral radiant flux $\langle \vec { F _ { \nu } } \rangle = - D _ { \nu } \vec { \nabla } \langle U _ { \nu } \rangle$ , so that the propagation of radiant energy carried by the numerical flux $\left. F _ { \nu } \right. > > \left. F _ { \nu } ^ { \mathrm { m a x } } \right. = c \langle U _ { \nu } \rangle$ could be faster than the actual transport by the group velocity $c$ , when the local spatial gradient of $\vec { \nabla } \langle U _ { \nu } \rangle$ is too large. The diffusion coefficient $\mathcal { D } _ { \nu }$ in Eq. (5.207) is replaced by a flux-limited diffusion coefficient $D _ { \nu }$ by introducing a nonlinear switch between the diffusion limit, $\kappa _ { \nu } \ll | \vec { \nabla } \langle U _ { \nu } \rangle / \langle U _ { \nu } \rangle |$ , and the free-streaming limit, $\kappa _ { \nu } \gg | \vec { \nabla } \langle U _ { \nu } \rangle / \langle U _ { \nu } \rangle |$ , to maintain the causality 

$$
\bar {D} _ {\nu} = \frac {c}{3 \kappa_ {\nu} + | \vec {\nabla} \langle U _ {\nu} \rangle / \langle U _ {\nu} \rangle |}. (5. 2 0 9)
$$

Equation (5.209) is called the harmonic flux limiter. Other ad hoc designs include Larsen flux limiter [OAH00] with the practical choice of $n = 2$ 

$$
\bar {D} _ {\nu} ^ {\text {L a r s e n}} = \frac {c}{\left[ \left(3 \kappa_ {\nu}\right) ^ {n} + | \vec {\nabla} \langle U _ {\nu} \rangle / \langle U _ {\nu} \rangle | ^ {n} \right] ^ {1 / n}}, \tag {5.210}
$$

and the min/max flux limiter 

$$
\bar {D} _ {\nu} ^ {\min / \max } = \min  \left[ D _ {\nu}, \langle F _ {\nu} ^ {\max } \rangle / | \vec {\nabla} \langle U _ {\nu} \rangle | \right]. \tag {5.211}
$$

A nature flux limiting technique for numerical solutions of transport equations was reported by Kershaw [Ker76]. The problem of overflow numerical flux is encountered whenever a particle transport equation is approximated by a parabolic type PDE such as the Spitzer-H¨arm electron thermal heat flux [SH53, SDMV81] in steep temperature gradients in laser-induced plasma, in which the upper bound of electron heat flux in the free-streaming limit is taken with $\alpha _ { \mathrm { e } } \sim 0 . 6 5$ in the 

min/max flux limiter, by replacing $\vec { \nabla } \langle U _ { \nu } \rangle$ with $\vec { \nabla } T _ { \mathrm { e } }$ . 

$$
\langle F _ {\mathrm {e}} ^ {\max} \rangle = \alpha_ {\mathrm {e}} n _ {\mathrm {e}} k _ {\mathrm {B}} T _ {\mathrm {e}} \sqrt {k _ {\mathrm {B}} T _ {\mathrm {e}} / m _ {\mathrm {e}}} \tag {5.212}
$$

In LTE, radiation temperatures of photons at different energy levels are not in equilibrium, leading to deviations from Planck spectrum. In the deceleration phase of ICF implosions, the mean free path for low energy group photons is less than the shell thickness, the strong absorption of low energy group photons leads to a rapid radiative thermal mass ablation at the inner shell surface that increases the density scale length and stabilize RT growth, whereas the escape of high energy group photons cools down the hot spot temperature, lower the neutron yield and degrades the alpha heating. To account for the transport of low and high energy group photons accurately in the deceleration phase, the multi-group radiation transport modeling is necessary. 

The multi-group flux-limited radiation diffusion approximation is obtained by integrating Eq. (5.207) over a frequency interval $( \nu _ { g } , \nu _ { g + 1 } )$ , where $g$ is the index of the $g$ th-group to eliminate the dependence on the frequency variable $\nu$ . 

$$
\frac {\partial \left\langle U _ {g} \right\rangle}{\partial t} = \vec {\nabla} \cdot \bar {D} _ {g} \left(\kappa_ {g} ^ {\mathrm {R}}\right) \vec {\nabla} \left\langle U _ {g} \right\rangle + c \kappa_ {g} ^ {\mathrm {P}} \left(B _ {g} - \left\langle U _ {g} \right\rangle\right), \tag {5.213}
$$

where the group-integrated radiation energy density $\langle U _ { g } \rangle$ and the group-integrated Planck spectrum $B _ { g }$ are 

$$
\langle U _ {g} (\vec {x}, t) \rangle = \int_ {\nu_ {g}} ^ {\nu_ {g + 1}} \left\langle U _ {\nu} (\vec {x}, \nu , t) \right\rangle d \nu , \tag {5.214}
$$

$$
B _ {g} \left(T _ {\mathrm {e}}\right) = \int_ {\nu_ {g}} ^ {\nu_ {g + 1}} \frac {4 \pi}{c} B _ {\nu} (\nu , T _ {\mathrm {e}}) d \nu . \tag {5.215}
$$

The group-weighted self-emission $B _ { g } ( T _ { \mathrm { e } } ; \vec { x } , t )$ is a function of the local electron temperature $T _ { \mathrm { e } } ( \vec { x } , t )$ , and is reduced to the black-body emission by integrating the 

frequency domain over $( 0 , \infty )$ , 

$$
B _ {\mathrm {b l a c k - b o d y}} (T _ {\mathrm {e}}) = \frac {8 k _ {\mathrm {B}} ^ {4} \pi^ {5}}{1 5 h ^ {3} c ^ {3}} T _ {\mathrm {e}} ^ {4} = \frac {4}{c} \sigma_ {\mathrm {S B}} T _ {\mathrm {e}} ^ {4}, \tag {5.216}
$$

where σSB = $\begin{array} { r } { \sigma _ { \mathrm { S B } } = \frac { 2 \pi ^ { 5 } k _ { \mathrm { B } } ^ { 4 } } { 1 5 h ^ { 3 } c ^ { 2 } } } \end{array}$ is the Stefen-Boltzmann constant. 

When interaction between radiation and hydrodynamics is strong, an advection term $\vec { \nabla } \cdot \vec { v } \langle U _ { g } \rangle$ and a work done term $\langle P _ { g } \rangle \vec { \nabla } \cdot \vec { v }$ must be added on the left hand side of Eq. (5.213), respectively, to transport radiation energy densities by the fluid motion, and to allow gaining radiation internal energies from plasma for converging flows $\vec { \nabla } \cdot \vec { v } < 0$ or releasing radiation internal energies into plasma for diverging flows $\vec { \nabla } \cdot \vec { v } > 0$ . Exact derivation for relativistic and non-relativistic radiation-hydrodynamics equations were reported by Pomraning [Pom73] and Castor [Cas04]. 

$\kappa _ { g }$ is group-averaged opacity. LTE tabular opacities are taken from Astrophysical Opacity Library [HMMA77] in DEC2D and DEC3D. Rosseland group-averaged opacity $\kappa _ { g } ^ { \mathrm { R } }$ weighs energetic photons more than Planck group-averaged opacity $\kappa _ { g } ^ { \mathrm { P } }$ , so that Rosseland opacities are accurate for optically thick plasma whereas Planck opacities are accurate for optically thin plasma. 

Rosseland group-averaged opacity: 

$$
\kappa_ {g} ^ {\mathrm {R}} (\rho , T _ {\mathrm {e}}) = \int_ {\nu_ {g}} ^ {\nu_ {g + 1}} \frac {\partial B _ {\nu} (\nu , T _ {\mathrm {e}})}{\partial T _ {\mathrm {e}}} d \nu / \int_ {\nu_ {g}} ^ {\nu_ {g + 1}} \frac {1}{\kappa_ {\nu} (\rho , T _ {\mathrm {e}})} \frac {\partial B _ {\nu} (\nu , T _ {\mathrm {e}})}{\partial T _ {\mathrm {e}}} d \nu , \tag {5.217}
$$

Planck group-averaged opacity: 

$$
\kappa_ {g} ^ {\mathrm {P}} (\rho , T _ {\mathrm {e}}) = \int_ {\nu_ {g}} ^ {\nu_ {g + 1}} \kappa_ {\nu} (\rho , T _ {\mathrm {e}}) B _ {\nu} (\nu , T _ {\mathrm {e}}) d \nu / \int_ {\nu_ {g}} ^ {\nu_ {g + 1}} B _ {\nu} (\nu , T _ {\mathrm {e}}) d \nu . \tag {5.218}
$$

# 5.4.2 Implementation in DEC2D and DEC3D

The following radiation transport package for the multi-group flux-limited radiation diffusion approximation model was implemented in DEC2D [WED+14] and DEC3D [WBB+15] by the author. The first application of the radiation transport package in DEC2D to compare effects of radiative thermal mass ablation on modifying the density gradient scale length on the inner shell surface between OMEGA and NIF deceleration phase simulations was published in Ref. [BWNB15]. The first code description of the package was presented in [WED+14]. 

$$
\frac {\partial \langle U _ {g} \rangle}{\partial t} + \vec {\nabla} \cdot \vec {v} \langle U _ {g} \rangle + \langle P _ {g} \rangle \vec {\nabla} \cdot \vec {v} = \vec {\nabla} \cdot \bar {D} _ {g} (\kappa_ {g} ^ {\mathrm {R}}) \vec {\nabla} \langle U _ {g} \rangle + c \kappa_ {g} ^ {\mathrm {P}} (B _ {g} - \langle U _ {g} \rangle). (5. 2 1 9)
$$

The group-weighted self-emission energy density $B _ { g }$ is a product of the black-body emission energy density $a \left( k _ { \mathrm { B } } T _ { \mathrm { e } } \right) ^ { 4 }$ , where $a = 8 \pi ^ { 5 } / ( 1 5 h ^ { 3 } c ^ { 3 } )$ is a constant, and the weight factor $0 \leq b ( u _ { g } , u _ { g + 1 } ) \leq 1$ to measure the fraction of emission of each energy group relative to the full frequency spectrum. 

$$
B _ {g} = a \left(k _ {\mathrm {B}} T _ {\mathrm {e}}\right) ^ {4} b (u _ {g}, u _ {g + 1}), \tag {5.220}
$$

where $u _ { g } = h \nu _ { g } / \left( k _ { \mathrm { B } } T _ { \mathrm { e } } \right)$ is the ratio of the photon energy $h \nu _ { g }$ to the local electron thermal energy $k _ { \mathrm { B } } T _ { \mathrm { e } } ( \vec { x } , t )$ . The weight factor is an integral of the normalized Planck distribution function, 

$$
b \left(u _ {g}, u _ {g + 1}\right) = \frac {1 5}{\pi^ {4}} \int_ {u _ {g}} ^ {u _ {g + 1}} \frac {u ^ {3}}{e ^ {u} - 1} d u. \tag {5.221}
$$

In the one-group radiation diffusion model, the normalized Planck integral in Eq. (5.221) is unity $b ( 0 , \infty ) = 1$ . The integrand gives the peak self-emission at the photon energy $h \nu _ { g } = 2 . 8 2 k _ { \mathrm { B } } T _ { \mathrm { e } }$ , and decays to negligible small values at $u > 1 5$ . The definite Planck integral is evaluated approximately using the Polylog functions 

$\begin{array} { r } { L i _ { n } z _ { g } = \sum _ { k = 1 } ^ { \infty } z _ { g } ^ { k } / k ^ { n } } \end{array}$ [Cla87] with $z _ { g } = e ^ { - u _ { g } }$ for large $u \geq 2 . 6$ , 

$$
\begin{array}{l} b \left(u _ {g} \geq 2. 6, u _ {g + 1}\right) \simeq \left[ \ln (z) ^ {3} L i _ {1} z - 3 \ln (z) ^ {2} L i _ {2} z \right. \tag {5.222} \\ + 6 \mathrm {L n} (z) L i _ {3} z - 6 L i _ {4} z ] _ {z _ {g}} ^ {z _ {g + 1}}, \\ \end{array}
$$

and Taylor series expansion for small $u < 2 . 6$ , 

$$
b \left(u _ {g}, u _ {g + 1} <   2. 6\right) \simeq \left[ \frac {u ^ {3}}{3} - \frac {u ^ {4}}{8} + \frac {u ^ {5}}{6 0} + \dots \right] _ {u _ {g}} ^ {u _ {g + 1}}. \tag {5.223}
$$

The radiation diffusion equation in Eq. (5.219) is solved fully implicitly in time to provide numerically stable solution. In every time-step, high energy photons are born within the high-temperature hot spot by the strong self-emission term $B _ { g } \sim T _ { \mathrm { e } } ^ { 4 }$ , and rapidly diffuse away from the hot spot into the cold shell, and finally escape into vacuum due to the long diffusion mean free path, and simultaneously free-free opacities in Eq. (5.186) vary significantly in space from the hot spot to the cold shell and vacuum. Figure (5.24) shows the spatial variation of Karmers free-free diffusion mean free path $\ell _ { \nu } = 1 / \kappa _ { \nu } ^ { \mathrm { H } }$ , defined in Eq. (5.186), at stagnation for a NIF implosion at photon energy level $h \nu = 2 . 1$ -keV, implies that the diffusion coefficient is strongly anisotropic in space. Since the hot spot is optically thin with large $\ell _ { \nu }$ , as shown in Fig. ((5.24)), the implementation of flux limiter is necessary to account for a proper radiation cooling within the hot spot due to the escape of high-energy-group radiation energy flux in the X-ray regime. 

To model the stiff source, the fast diffusion and the strong spatial anisotropy of opacities in radiation transport, HYPRE algebraic multi-grid [TOS01] preconditioned GMRES [SS86] solver is used to invert the large non-symmetric sparse matrix, resulted from the fully implicit discretization. Classical red-black successiveover-relaxation iteration solver [Eva84] was observed to converge extremely slowly in our numerical tests. The tridiagonal matrix for the implicit discretization in Eq. (5.219), without the advection and $P d V$ work terms, is solved by Gaussian 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/a4029f66c401662d52ec4504348c19a68c424330dfcb062ddedd811b96500235.jpg)



Figure 5.24: The plot of Kramers free-free diffusion mean free path $\ell _ { \nu } = 1 / \kappa _ { \nu }$ in Eq. (5.186) using the mass density and electron profiles for a 1-D NIF implosion at stagnation. The diffusion coefficient is changed significantly in space from the hot spot, to the cold shell and finally to the outer vacuum.


elimination in DEC2D and by HYPRE using the multi-grid iteration in DEC3D. 

Equation (5.219) is solved by a hydro-step and a diffusion-step by the operator splitting technique: [Har11] (1) update the radiation energy densities for each group due to the advection by flow and the work done on fluid 

$$
\frac {\partial \langle U _ {g} \rangle}{\partial t} + \vec {\nabla} \cdot \vec {v} \langle U _ {g} \rangle + \langle P _ {g} \rangle \vec {\nabla} \cdot \vec {v} = 0, \tag {5.224}
$$

by HLLC Riemann solver in hydrodynamics, and (2) update the resulting multigroup radiation diffusion equation with stiff source terms by a fully implicit in time and second-order in space discretization. 

$$
\frac {\langle U _ {g} \rangle^ {n + 1} - \langle U _ {g} \rangle^ {n}}{\triangle t} = \vec {\nabla} \cdot \bar {D} _ {g} ^ {n} \vec {\nabla} \langle U _ {g} \rangle^ {n + 1} + c \kappa_ {g} ^ {n} \left(B _ {g} ^ {n} - \langle U _ {g} \rangle^ {n + 1}\right). (5. 2 2 5)
$$

Here $\kappa _ { g } ^ { n }$ is the tabular Planck group-averaged LTE opacities, [HMMA77] which is a function of the local material density and the local electron temperature at the time level $t ^ { n }$ . The radiation-material interaction is modeled by coupling the net change of the radiation emission source term and the absorption sink term for a 

total number of $N$ groups to the electron temperature equation explicitly. 

$$
\frac {3}{2} n _ {\mathrm {e}} ^ {n} k _ {\mathrm {B}} \left(\frac {T _ {\mathrm {e}} ^ {n + 1} - T _ {\mathrm {e}} ^ {n}}{\triangle t}\right) = - \sum_ {g = 1} ^ {N} c \kappa_ {g} ^ {n} \left(B _ {g} ^ {n} - \langle U _ {g} \rangle^ {n + 1}\right) \tag {5.226}
$$

Currently $N = 4$ , 12 and 48 groups are available, and can be extended to more energy groups directly. The 12-group radiation transport in OMEGA and NIF deceleration phase simulations was observed works well to resolve the absorption of low energy photons on the inner shell surface. 

# 5.4.3 Boundary condition

The boundary condition of zero incoming radiation flux is imposed at the boundary surface $A ( \vec { x } , t )$ of the simulation domain at time $t$ . 

$$
F _ {\nu} ^ {\text {i n}} = - \hat {n} \cdot \vec {F} _ {\nu} (\hat {n} \cdot \Omega <   0; A (\vec {x}, t)) = 0. \tag {5.227}
$$

$\hat { n } = \vec { \nabla } A ( \vec { x } , t ) / A ( \vec { x } , t )$ is a unit vector normal to the boundary surface. Using P1 approximation in Eq. (5.202) to solve for the incoming radiation flux $F _ { \nu } ^ { \mathrm { i n } }$ . Without loss of generality, $\hat { n } = \hat { z }$ can be assumed to perform the integration and write $\langle \vec { F _ { \nu } } \rangle = \hat { n } \cdot \langle \vec { F _ { \nu } } \rangle \hat { n } + \hat { n } _ { \perp } \cdot \langle \vec { F _ { \nu } } \rangle \hat { n } _ { \perp }$ , where $\hat { n } _ { \perp }$ is a unit vector orthogonal to $\hat { n }$ . [AMJKM+05] 

$$
- \hat {n} \cdot \vec {F} _ {\nu} (\hat {n} \cdot \Omega <   0, A (t)) = - \frac {1}{4 \pi} \int_ {\hat {n} \cdot \Omega <   0} \hat {n} \cdot \vec {\Omega} \left(c \langle U _ {\nu} \rangle + 3 \vec {\Omega} \cdot \langle \vec {F} _ {\nu} \rangle\right) d \Omega , \tag {5.228}
$$

which is equal to $\begin{array} { r l } & { F _ { \nu } ^ { \mathrm { i n } } = - \frac { 1 } { 4 \pi } \int _ { - 1 } ^ { 0 } \int _ { 0 } ^ { 2 \pi } \mu \left( c \langle U _ { \nu } \rangle + 3 \mu \langle \vec { F } _ { \nu } \rangle \cdot \hat { n } \right) d \mu d \phi } \end{array}$ . Here the symmetric integration for the component of $\vec { \Omega } \cdot \hat { n } _ { \perp }$ over $2 \pi$ in angle $\phi$ vanishes, so that the solution for the incoming radiation flux is 

$$
F _ {\nu} ^ {\text {i n}} = \frac {c}{4} \left\langle U _ {\nu} \right\rangle - \frac {1}{2} \left\langle \vec {F} _ {\nu} \right\rangle \cdot \hat {n}. \tag {5.229}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/cdc7f109a2e8aa765f71789e2e32c0d01c57ad4a716f16968ab629f5545281e8.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/691621702168d4c7ae02fe61004da31e3036a7cc1e0f0235ddb9407914cfdfec.jpg)



Figure: No flux limiter


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/03ac1c9d2f95585d6c09e3fb119cc11b03288b2bdf0a8d91bbf05778ba8ad71a.jpg)



Figure:Harmonic flux limiter


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/6c0166f2989bb8d2b474e40a9b2b865ab7799665c9e0767b9bdd1c95f0d65c3b.jpg)



Figure: Larsen flux limiter



Figure: Min/Max flux limiter



Figure 5.25: Benchmark tests of the 12-group radiation transport in DEC2D with $L I L A C$ in a 1-D spherical slab. The initial electron temperature and mass density profiles are $T _ { \mathrm { e } } ( r < r _ { 0 } ) = 5 \ ]$ keV and $T _ { \mathrm { e } } ( r \geq r _ { 0 } ) = 0 . 5 ~ \mathrm { k e V }$ ; $\rho ( r < r _ { 0 } ) = 5 0 ~ \mathrm { g / c m ^ { 3 } }$ and $\rho ( r \geq r _ { 0 } ) = 1 0 0 ~ \mathrm { g / c m ^ { 3 } }$ .


The vacuum boundary condition of $F _ { \nu } ^ { \mathrm { i n } } = 0$ with $\hat { n } = \hat { r }$ in $D E C \mathcal { Z } D$ and DEC3D is reduced to the Marshak boundary condition, [Bru02] meaning that a half of the maximum radiation flux $\langle F _ { \nu } ^ { \mathrm { m a x } } \rangle = c \langle U _ { \nu } \rangle$ leaves the simulation domain. 

$$
\langle \vec {F} _ {\nu} \rangle \cdot \hat {r} = \frac {c \langle U _ {\nu} \rangle}{2}. \tag {5.230}
$$

Figure (5.25) shows the result of four benchmark tests of the 12-group radiation transport in DEC2D with $L I L A C$ between without and with harmonic, Larsen and min/max flux limiters, defined in Eqs. (5.209), (5.210) and (5.211) respectively. The geometry is a 1-D spherical slab with a radius of $R = 1 0 0 \ \mu \mathrm { m }$ , separated by two initial states at the radius $r _ { 0 } = 5 0 ~ \mu \mathrm { m }$ . The initial electron temperature profile: $T _ { \mathrm { e } } ( r < r _ { 0 } ) = 5 ~ \mathrm { k e V }$ and $T _ { \mathrm { e } } ( r \geq r _ { 0 } ) = 0 . 5 \ \mathrm { k e V }$ , and the initial material mass density profile: $\rho ( r < r _ { 0 } ) = 5 0 ~ \mathrm { g / c m ^ { 3 } }$ and $\rho ( r \geq r _ { 0 } ) = 1 0 0 ~ \mathrm { g / c m ^ { 3 } }$ , are chosen as typical values for the hot spot and the cold shell at stagnation for OMEGA and NIF implosions. The black curve is the $L I L A C$ solution with the same flux-limited radiation transport in four plots in Fig. (5.25). A good agreement is observed in 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/479c6f6ded3dc1e2bcb51d5369ca7eae1ee86dc9cf38c7e3d11f937635d4cb8f.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/b78564c34a7a12f4796f27cf108a0a3500d91b1ca79a673b89964e1ff13ed98f.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/06a72b533e244062fb249d53b525b21ca767c72565cc8ef0c2ab3dd22fe08a7b.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/43f9e5be9c4b71d45f919e7c9a9c8d8e75119588a5a37fea506193abba2ba1b1.jpg)



Figure 5.26: The benchmark test of multi-group radiation transport implemented in DEC2D. NIF* standards for National Ignition Facility.


simulating the radiative cooling within the hot spot, with $T _ { \mathrm { e } } ( r < r _ { 0 } )$ dropping from 5-keV to about 2-keV through the Marshak boundary condition at the edge $R$ . All forms of flux-limiters work well to transport the radiation flux properly from a optically thin hot spot ( $r < r _ { 0 }$ ) to a optically thick shell ( $r \leq r _ { 0 }$ ). Without any flux limiter, the simulated hot spot temperature $T _ { \mathrm { e } } ( r \textless r _ { 0 } )$ by $D E C \mathcal { Z } D$ is shown to be lower than that for $L I L A C$ , implying a overflow of high energy group radiation flux into the vacuum. 

Figure (5.26) shows the result of benchmark tests of the 12-group radiation transport implemented in DEC2D with LILAC. Effects of radiative mass ablation on the inner shell surface are compared between OMEGA and NIF implosions without and with the 12-group radiation transport. For OMEGA implosions, the effect of radiation transport is observed only leading to radiation cooling of the hot spot by decreasing the electron temperature, but negligible radiative mass ablation. For NIF implosions, in addition to radiation cooling, the diffusion mean free path of low energy photons is shorter than the dimensions of the hot spot and the cold shell, resulting in rapid radiation absorptions on the inner shell surface and significant radiative mass ablation. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/a41874e89d944081d66cda6c8cd841ecc69e673b49ee6399b41c001b86b32c22.jpg)


Figure 5.27: The flow chart of parallel multi-group radiation transport implemented in DEC2D. The diffusion equation for each group is solved independently by multiple core in parallel for each time step, followed by summing over all groups to obtained the net energy exchange with the plasma through an explicit update on the electron temperature. 

<table><tr><td>Single time-step run time (s)</td><td>Grid number (N × N)</td><td>N = 200</td><td>N = 400</td><td>N = 600</td><td>N = 800</td></tr><tr><td rowspan="3">Parallel code</td><td>Four-group in parallel</td><td>0.075</td><td>0.19</td><td>0.39</td><td>0.67</td></tr><tr><td>Other (hydro, thermal, and alpha)</td><td>0.100</td><td>0.46</td><td>1.00</td><td>1.80</td></tr><tr><td>Other + four-group in parallel</td><td>0.175</td><td>0.65</td><td>1.39</td><td>2.47</td></tr><tr><td rowspan="3">Serial code</td><td>Four-group in serial</td><td>0.130</td><td>0.63</td><td>1.4</td><td>2.6</td></tr><tr><td>Other (hydro, thermal, and alpha)</td><td>0.092</td><td>0.40</td><td>0.94</td><td>1.7</td></tr><tr><td>Other + four-group in serial</td><td>0.222</td><td>1.03</td><td>2.34</td><td>4.3</td></tr></table>


Figure 5.28: The performance of parallel multi-group radiation transport implemented in DEC2D. The parallel code for 4-group radiation transport is about $\sim 2 \times$ $\times$ faster than the serial code. The total amount of real computational time saved by the parallel code is increased with the resolution from $N = 2 0 0$ to $N = 8 0 0$ in the test.


# 5.4.4 Parallel multi-group radiation transport in DEC2D

The parallel simulation strategy to solve multi-group diffusions such as for radiation and alpha particle transport is essential for high performance 3-D simulations. Figure 5.27 shows the parallel architecture implemented in DEC2D to solve the radiation diffusion equation in Eq. (5.225) of each energy group in multiple cores in every time step through message-passing-interface (MPI) [CLS94]. After the direct solve of the tridiagonal matrix of the implicit discretization by Gaussian elimination, the net energy exchange with the plasma due to emission and absorption is summing over all groups in the core 0. Figure 5.28 shows that the performance for a parallel code is about $\sim 2 \times$ $\times$ faster than a serial code in solving 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/1c72ec921d4d8b4c66a5c48a99a6f034fb045511a94e06cbdc79679843d33122.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/a091afa19a4a5c1150d910baff62cbf00855381f7cecf108a1e9a503aa5b29af.jpg)



Figure 5.29: Comparison of effects of radiation ablative stabilization on NIF implosions between with and without radiation transport.


a 4-group radiation transport. 

# 5.4.5 Radiative ablative stabilization of RT instabilities

In decelerating ICF shells, RT instabilities are mitigated by the mass ablation on the inner shell surface caused by the heat flow from the central hot spot. The ablative stabilization of deceleration phase RT instabilities reduces the RT growth rate $\gamma _ { \ell }$ as a unique function of the shell deceleration $\ddot { R }$ , the minimum density gradient scale length $L _ { \mathrm { m i n } }$ , the mass ablation velocity on the inner shell velocity $v _ { \mathrm { a } }$ , and a numerical coefficient $\beta _ { \mathrm { m i n } } \sim 1 . 5$ . [AtV04] 

$$
\gamma_ {\ell} = \sqrt {\frac {\ddot {R} k}{1 + L _ {\mathrm {m i n}} k}} - \beta_ {\mathrm {m i n}} k v _ {\mathrm {a}}, \tag {5.231}
$$

where $k = \ell / R$ is the wavenumber for a given mode $\ell$ . The importance of an increasing density gradient scale length in Eq. (5.231) is shown reducing the RT growth rate for high mode, but vanishing effect for low modes in the limit of $k  0$ . Sources of heat flow in the deceleration phase of ICF implosions include the radiative transfer of low energy photons, electron and ion heat transfer and alpha particle energy deposition. Effect of radiation ablative stabilization for a NIF 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/67b03977c512a0cb6f90c39a709d2e0a2c1547057817fbcc290ddd9e493a9833.jpg)



$\rho _ { \mathrm { N I F , s t a g } } ^ { \mathrm { r a d } } ( \ell = 6 , m = 0 )$


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/dc1c99667766f11fe9dc1ab3361e3dfb0737856e772f88d4c210e97ca147195b.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/99058d20dc0f249c1e0cae9691af581d960fe4006e465aafd5d86b7772e197b7.jpg)



$\rho _ { \mathrm { N I F , s t a g } } ^ { \mathrm { r a d } } ( \ell = 2 0 , m = 0 )$


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/7eb5adc64cb6a623b65ce8d923ac8aa73ff8e60e08cf7a5a3ffa06dc7b99dde6.jpg)



Figure 5.30: Effect of radiation ablative stabilization for a 2-D low mode $\ell = 6$ (above), and a high mode $\ell = 2 0$ (below).


implosion is studied in Fig. (5.29), showing a sharp decrease of peak mass densities for RT spikes in the presence of radiation transport. Figure (5.30) compares effects of radiation ablative stabilization between a low mode $\ell = 6$ and a high mode $\ell = 2 0$ . The rapid absorption of low energy photons near the inner shell surface ablates the mass off from the cold shell at velocity $v _ { \mathrm { a } }$ , leading to the appearance of lower peak mass densities shown in Fig. (5.29) and thicker RT spikes shown in Figs. (5.30). The mushroom structure in the non-linear stage for the high mode $\ell = 2 0$ is shown significantly stabilized, but less effect on the low mode $\ell = 6$ . 

# 5.5 Electron and ion heat conductions

# 5.5.1 Spitzer H¨arm thermal diffusion

The basic theory for electron heat flux $\vec { F } _ { \mathrm { S p i t z e r } } = - \kappa _ { \mathrm { e } } \vec { \nabla } T _ { \mathrm { e } }$ within a fully ionized plasma are derived from the first-order perturbation in Fokker-Planck equation. [SH53] When the gradient of electron temperature $T _ { \mathrm { e } }$ is too large such as at the thermal heat front in a laser induced plasma, the electron heat flow must be fluxlimited to avoid the overflow of numerical fluxes $F _ { \mathrm { S p i t z e r } } > \langle F _ { \mathrm { e } } ^ { \mathrm { m a x } } \rangle$ , see Eq. (5.212). A proper modeling of electron and ion heat conductions plays an important role to determine accurate ablative RT instabilities, as well as accurate neutron yield productions in strong alpha heating studies. The transport coefficients for electron thermal conductivity $\kappa _ { \mathrm { e } }$ and ion thermal conductivity $\kappa _ { \mathrm { i } } = \kappa _ { \mathrm { e } } \sqrt { m _ { \mathrm { e } } / m _ { \mathrm { i } } }$ by Spitzer [SH53, Spi90], where $m _ { \mathrm { e } }$ and $m _ { \mathrm { i } }$ are electron and ion masses respectively, are adopted in DEC2D and DEC3D. 

$$
\kappa_ {\mathrm {e}} = \underbrace {\frac {2 0 (2 / \pi) ^ {3 / 2} (k _ {\mathrm {B}} T _ {\mathrm {e}}) ^ {5 / 2} k _ {\mathrm {B}}}{m _ {\mathrm {e}} ^ {1 / 2} \bar {Z} e ^ {4} \ln \Lambda_ {\mathrm {e i}} ^ {\mathrm {L M}}}} _ {\kappa_ {L}} \underbrace {\frac {0 . 0 9 5 (\bar {Z} + 0 . 2 4)}{1 + 0 . 2 4 \bar {Z}}} _ {\delta} f _ {\mathrm {L M}}, \qquad (5. 2 3 2)
$$

where $k _ { \mathrm { B } }$ is the Boltzmann constant, $e$ is the electron charge and $Z$ is the effective charge of ions. $\kappa _ { L }$ is the thermal conductivity for a Lorentz gas, [Spi90] which is an ideal fully ionized gas in which electrons do not interact with each other, and all the positive ions are at rest. $\delta \leq 0 . 3 9 6$ is the Spitzer prefactor to reduce $\kappa _ { L }$ to the effective coefficient of thermal conductivity $\delta \kappa _ { L }$ due to the thermoelectric effect. Because the secondary electric field $\vec { E }$ induces an opposite current to cancel the flow of current caused by the temperature gradient $\vec { \nabla } T$ in order to maintain a finite total current $\vec { j } = \sigma \vec { E } + \alpha \vec { \nabla } T$ in the steady state. The original value [Spi90] for the upper bound of $\delta$ with $\bar { Z } \to \infty$ reported by Spitzer is given by the produce of $\delta _ { T } \times \epsilon = 0 . 4$ , where $\delta _ { T }  1$ and $\epsilon  0 . 4$ . Derivations for the electrical 

conductivity $\sigma$ and the coefficient $\alpha$ were reported by Spitzer in Ref. [SH53]. fLM is the Lee-More degeneracy correction factor [HCB+14, LM84], 

$$
f _ {\mathrm {L M}} (n _ {\mathrm {e}}, T _ {\mathrm {e}}) = 1 + \frac {3 \pi^ {5}}{5 1 2 0 0} \left(\frac {T _ {F}}{T _ {\mathrm {e}}}\right) ^ {3} \times \delta^ {- 2}, \tag {5.233}
$$

where TF = 2mekB $\begin{array} { r } { T _ { F } \ = \ \frac { \hbar ^ { 2 } } { 2 m _ { \mathrm { e } } k _ { \mathrm { B } } } \left( 3 \pi ^ { 2 } n _ { \mathrm { e } } \right) ^ { 2 / 3 } } \end{array}$ ~2 is the Fermi temperature for electrons in a fully ionized plasma, and $n _ { \mathrm { e } }$ is the electron number density. Consider a cold shell at stagnation with a mass density $\rho = 1 0 0 ~ \mathrm { g / c m ^ { 3 } }$ and an electron temperature $T _ { \mathrm { e } } = 5 0 0 ~ \mathrm { e V }$ , Lee-More degeneracy factor is $f _ { \mathrm { L M } } \sim 1 . 4$ . $\mathrm { l n } \Lambda _ { \mathrm { e i } } ^ { \mathrm { L M } }$ is the Lee-More electron-ion Coulomb logarithm with a minimum cutoff of 2, [LM84] 

$$
\ln \Lambda_ {\mathrm {e i}} ^ {\mathrm {L M}} = \operatorname {M a x} \left[ \frac {1}{2} \ln \left[ 1 + (b _ {\max} / b _ {\min}) ^ {2} \right], 2 \right]. \tag {5.234}
$$

In the classical theory of a scattering event between an incoming electron and a stationary ion, the electron-ion Coulomb logarithm $\begin{array} { r } { \mathrm { l n } \Lambda _ { \mathrm { e i } } = \int _ { 0 } ^ { b _ { \mathrm { m a x } } / b _ { \mathrm { m i n } } } x ^ { 3 } ( 1 + x ^ { 2 } ) ^ { - 2 } d x } \end{array}$ is defined by integrating the square of the transverse velocity of the incoming electron from the minimum $b _ { \mathrm { m i n } }$ to the maximum $b _ { \mathrm { m a x } }$ impact parameters. At the large ratio of $b _ { \operatorname* { m a x } } / b _ { \operatorname* { m i n } }$ , the solution for the electron-ion Coulomb logarithm 

$$
\ln \Lambda_ {\mathrm {e i}} = \frac {1}{2} \left[ \frac {1}{1 + (b _ {\mathrm {m a x}} / b _ {\mathrm {m i n}}) ^ {2}} - 1 + \ln [ 1 + (b _ {\mathrm {m a x}} / b _ {\mathrm {m i n}}) ^ {2} ] \right] \tag {5.235}
$$

is approximated by Lee-More in Eq. (5.234), and by Spitzer [Spi90] as 

$$
\ln \Lambda_ {\mathrm {e i}} ^ {\mathrm {S p i t z e r}} = \ln [ b _ {\mathrm {m a x}} / b _ {\mathrm {m i n}} ]. \tag {5.236}
$$

The maximum impact parameter $b _ { \mathrm { m a x } }$ is set to be the Debye-H¨uckel screening length $\lambda _ { \mathrm { D H } }$ as a function of electron and ion plasma properties [LM84] 

$$
\frac {1}{\lambda_ {\mathrm {D H}} ^ {2}} = \frac {4 \pi n _ {\mathrm {e}} e ^ {2}}{k _ {\mathrm {B}} (T _ {\mathrm {e}} ^ {2} + T _ {F} ^ {2}) ^ {1 / 2}} + \frac {4 \pi n _ {\mathrm {i}} (\bar {Z} e) ^ {2}}{k _ {\mathrm {B}} T _ {\mathrm {i}}}. \tag {5.237}
$$

However, typical Fermi temperatures are much smaller than electron temperatures in the hot spot and the cold shell in the entire deceleration phase, and can be neglected in Eq. (5.237). The minimum impact parameter $b _ { \mathrm { m i n } }$ is set to be the classical distance of closest approach by Spitzer 

$$
b _ {\mathrm {m i n}} ^ {\mathrm {S p i t z e r}} = \bar {Z} e ^ {2} / (m _ {\mathrm {e}} v _ {\mathrm {e}} ^ {2}), \tag {5.238}
$$

where the electron thermal velocity is $v _ { \mathrm { e } } = \sqrt { 3 k _ { \mathrm { B } } T _ { \mathrm { e } } / m _ { \mathrm { e } } }$ [Spi90], and by Lee-More $b _ { \operatorname* { m i n } } ^ { \mathrm { L M } }$ by considering the minimum distance between $b _ { \mathrm { m i n } } ^ { \mathrm { S p i t z e r } }$ and the half of the de-Broglie wavelength $\lambda = h / ( m _ { \mathrm { e } } v _ { \mathrm { e } } )$ of an electron [LM84] 

$$
b _ {\mathrm {m i n}} ^ {\mathrm {L M}} = \min \left[ b _ {\mathrm {m i n}} ^ {\mathrm {S p i t z e r}}, \lambda / 2 \right]. (5. 2 3 9)
$$

The local theory of heat transport for electrons and ions in a plasma is a collective effect of successive small and large angle scatterings described by a diffusion approximation. 

$$
\frac {3}{2} n _ {\mathrm {s}} k _ {\mathrm {B}} \frac {\partial T _ {\mathrm {s}}}{\partial t} = \vec {\nabla} \cdot \kappa_ {\mathrm {s}} \vec {\nabla} T _ {\mathrm {s}}, \quad \mathrm {w h e r e s = e , i} \tag {5.240}
$$

# 5.5.2 Electron and ion equilibration

The energy exchange between electrons and ions is determined by the Spitzer electron-ion collision time scale, [Spi90] 

$$
\tau_ {\mathrm {e i}} = \frac {3 m _ {\mathrm {e}} m _ {\mathrm {i}} k _ {\mathrm {B}} ^ {3 / 2} (T _ {\mathrm {e}} / m _ {\mathrm {e}} + T _ {\mathrm {i}} / m _ {\mathrm {i}}) ^ {3 / 2}}{8 \sqrt {2 \pi} n _ {\mathrm {i}} \bar {Z} ^ {2} e ^ {4} \ln \Lambda_ {\mathrm {e i}} ^ {\mathrm {S p i t z e r}}}, \tag {5.241}
$$

where $\ln \Lambda _ { \mathrm { e i } } ^ { \mathrm { S p i t z e r } } = 3 / ( 2 \bar { Z } e ^ { 3 } ) \sqrt { ( k _ { \mathrm { B } } T _ { \mathrm { e } } ) ^ { 3 } / ( \pi n _ { \mathrm { e } } ) }$ is the Spitzer electron-ion Coulomb logarithm defined in Eq. (5.236). 

$$
\frac {\partial T _ {\mathrm {e}}}{\partial t} = - \frac {T _ {\mathrm {e}} - T _ {\mathrm {i}}}{\tau_ {\mathrm {e i}}}, \quad \frac {\partial T _ {\mathrm {i}}}{\partial t} = - \frac {T _ {\mathrm {i}} - T _ {\mathrm {e}}}{\tau_ {\mathrm {e i}}}. \tag {5.242}
$$

By subtracting the rate of ion heat transfer from the electron’s, the resulting coupled equation $\begin{array} { r } { \frac { \dot { \sigma } } { \partial t } ( T _ { \mathrm { e } } - T _ { \mathrm { i } } ) = - 2 ( T _ { \mathrm { e } } - T _ { \mathrm { i } } ) / { \tau _ { \mathrm { e i } } } } \end{array}$ has a simple analytical solution in each time step, 

$$
(T _ {\mathrm {e}} - T _ {\mathrm {i}}) ^ {n + 1} = \left(T _ {\mathrm {e}} - T _ {\mathrm {i}}\right) ^ {n} \operatorname {E x p} \left(- \frac {2 \triangle t}{\tau_ {\mathrm {e i}}}\right). \tag {5.243}
$$

The minus sign in the exponent in Eq. (5.243) means that the unbalanced heat flow between electrons and ions eventually vanishes when two species are put in contact for long enough time. 

# 5.5.3 Fully implicit second-order in space discretization

# 1-D discretization

A fast 1-D simulation code, written in Mathematica, was developed to test various numerical algorithms in code development for hydrodynamics, electron and ion thermal diffusions and equilibration in 1-D spherical geometry. In every timestep, the operator-splitting technique is applied to update (1) hydrodynamics: $\partial _ { t } \vec { Q } + \partial _ { r } \vec { F } = 0$ , (2) geometrical sources: $\partial _ { t } \vec { Q } = \vec { S }$ , (3) electron-and-ion thermal diffusions in Eq. (5.240), and (4) the electron-and-ion equilibration in Eq. (5.242), where partial derivatives are denoted by $\begin{array} { r } { \partial _ { t } = \frac { \partial } { \partial t } } \end{array}$ and $\begin{array} { r } { \partial _ { r } = \frac { \partial } { \partial r } } \end{array}$ . Radiation and alpha particle transport packages are not installed. 

The first approximation is the time-lagging to treat the material density $n _ { \mathrm { s } }$ as a constant in time 

$$
\partial_ {t} \left(\frac {3}{2} n _ {\mathrm {s}} k _ {\mathrm {B}} T _ {\mathrm {s}}\right) = \left(\frac {3}{2} n _ {\mathrm {s}} k _ {\mathrm {B}}\right) \partial_ {t} T _ {\mathrm {s}} + \mathcal {O} _ {2}, \tag {5.244}
$$

where the second-order term $\mathcal { O } _ { 2 } = \left( \textstyle \frac { 3 } { 2 } k _ { \mathrm { B } } T _ { \mathrm { s } } \right) \partial _ { t } n _ { \mathrm { s } }$ is ignored. Omit the species labels and absorb Boltzmann constant into the temperature variable $T = k _ { \mathrm { { B } } } T _ { \mathrm { { s } } }$ and the thermal conductivity variable $\kappa = \kappa _ { \mathrm { s } } / k _ { \mathrm { B } }$ . The fully implicit second-order in space 

discretization for a thermal diffusion problem in 1-D spherical geometry, 

$$
\frac {3}{2} n _ {i} ^ {n} \left(\frac {T _ {i} ^ {n + 1} - T _ {i} ^ {n}}{\triangle t}\right) = r _ {i} ^ {- 2} \partial_ {r} \left(r ^ {2} \kappa^ {n} \partial_ {r} T ^ {n + 1}\right), \tag {5.245}
$$

forms a tri-diagonal matrix 

$$
- T _ {i} ^ {n} = Z _ {+} ^ {n} T _ {i + 1} ^ {n + 1} - \left(1 + Z _ {+} ^ {n} + Z _ {-} ^ {n}\right) T _ {i} ^ {n + 1} + Z _ {-} ^ {n} T _ {i - 1} ^ {n + 1}, \tag {5.246}
$$

$$
R _ {i} ^ {n} = C _ {i} ^ {n} T _ {i + 1} ^ {n + 1} + B _ {i} ^ {n} T _ {i} ^ {n + 1} + A _ {i} ^ {n} T _ {i - 1} ^ {n + 1}, \tag {5.247}
$$

where $R _ { i } ^ { n } = - I _ { i } ^ { n }$ , $C _ { i } ^ { n } = Z _ { + } ^ { n }$ , $B _ { i } ^ { n } = - ( 1 + Z _ { + } ^ { n } + Z _ { - } ^ { n } )$ and $A _ { i } ^ { n } = Z _ { - } ^ { n }$ . The variable $Z _ { \pm } ^ { n }$ computed at the time level $t ^ { n }$ is defined as 

$$
Z _ {\pm} ^ {n} = \frac {(r ^ {2} \kappa) _ {i \pm 1 / 2} ^ {n}}{r _ {i} ^ {2} \triangle r ^ {2} (\frac {3}{2} n _ {i} ^ {n} / \triangle t)}. (5. 2 4 8)
$$

The diffusion coefficient at the cell interface must be treated accurately using the simple mean between two adjacent cells at radius $r _ { i }$ and $r _ { i \pm 1 }$ : $\left( r ^ { 2 } \kappa \right) _ { i \pm 1 / 2 } ^ { n } =$ $\left\lfloor \left( r ^ { 2 } \kappa \right) _ { i } ^ { n } + \left( r ^ { 2 } \kappa \right) _ { i \pm 1 } ^ { n } \right\rfloor / 2$ . Otherwise, significant errors are resulted when the thermal conductivity is treated as a constant in space using this form $\left( r ^ { 2 } \kappa \right) _ { i \pm 1 / 2 } ^ { \pi } =$ $\kappa _ { i } ^ { n } \left\lfloor \left( r ^ { 2 } \right) _ { i } ^ { n } + \left( r ^ { 2 } \right) _ { i \pm 1 } ^ { n } \right\rfloor / 2$ , which was observed leading to a large overflow of heat flux into the cold shell. Equations (5.245)–(5.248) are derived using the second-order center-in-space discretization for the spatial derives. 

$$
\partial_ {r} \left(r ^ {2} \kappa^ {n} \partial_ {r} T ^ {n + 1}\right) = \left[ \left(r ^ {2} \kappa^ {n} \partial_ {r} T ^ {n + 1}\right) _ {i + 1 / 2} - \left(r ^ {2} \kappa^ {n} \partial_ {r} T ^ {n + 1}\right) _ {i - 1 / 2} \right] / \triangle r,
$$

$$
\partial_ {r} T _ {i + 1 / 2} ^ {n + 1} = \left(T _ {i + 1} ^ {n + 1} - T _ {i} ^ {n + 1}\right) / \triangle r,
$$

$$
\partial_ {r} T _ {i - 1 / 2} ^ {n + 1} = \left(T _ {i} ^ {n + 1} - T _ {i - 1} ^ {n + 1}\right) / \triangle r,
$$

The boundary condition of zero heat flux $\vec { \nabla } T = 0$ is applied at the origin only for DEC2D because of the azimuthal rotational symmetry; and is applied at the edges of the simulation domain in both $D E C \mathcal { Z } D$ and $D E C 3 D$ , meaning for no heat 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/107587151c9d4071706cb87c4f6cffebf9cb29d9ea435c99d8047456381a954a.jpg)



Figure 5.31: The boundary condition of zero temperature gradient is applied at the origin by replacing $B _ { 1 }$ with $B _ { 1 } ^ { * }$ and at the edge of the simulation domain by replacing $B _ { N }$ with $B _ { N } ^ { * }$ for a 1-D thermal diffusion in spherical geometry. Temperatures are discretized at the center of a cell i.e., $r _ { i } = ( i - 1 / 2 ) \triangle r$ .


exchange with the outer vacuum. Figure (5.31) shows the implementation of the zero heat flux boundary conditions in 1-D thermal diffusion in spherical geometry. 

# 2-D discretization

The thermal diffusion equation with the azimuthal rotational symmetry in DEC2D 

$$
\frac {3}{2} n _ {i, j} ^ {n} \left(\frac {T _ {i , j} ^ {n + 1} - T _ {i , j} ^ {n}}{\triangle t}\right) = r _ {i} ^ {- 1} \partial_ {r} \left(r \kappa^ {n} \partial_ {r} T ^ {n + 1}\right) + \partial_ {z} \left(\kappa^ {n} \partial_ {z} T ^ {n + 1}\right), \tag {5.249}
$$

is directional-splitting into two steps: (1) a thermal diffusion update in the $r$ - direction, (2) followed by another thermal diffusion update in the $z$ -direction, 

$$
\frac {3}{2} n _ {i, j} ^ {n} \left(\frac {T _ {i , j} ^ {n + 1} - T _ {i , j} ^ {n}}{\Delta t}\right) = r _ {i} ^ {- 1} \partial_ {r} \left(r \kappa^ {n} \partial_ {r} T ^ {n + 1}\right), \tag {5.250}
$$

$$
\frac {3}{2} n _ {i, j} ^ {n} \left(\frac {T _ {i , j} ^ {n + 1} - T _ {i , j} ^ {n}}{\triangle t}\right) = \partial_ {z} \left(\kappa^ {n} \partial_ {z} T ^ {n + 1}\right). \tag {5.251}
$$

During the thermal diffusion in each direction, results of 1-D discretization can be applied directly by replacing the power of 2 on the right hand side of Eq. (5.245) with 1 for the $r$ -diffusion and with 0 for the $z$ -diffusion. The same boundary condition as in 1-D with the zero heat flux is imposed at the origin and the edges of the simulation domain. 

# 3-D discretization

A general fully implicit second-order in space diffusion solver, with the option to specify source terms, is implemented in DEC3D, in order to support all types of diffusion approximations for thermal, radiation and alpha particle transport. The boundary condition of zero heat flux is applied only at the edge of the simulation domain, whereas a full 3-D condition is applied to allow heat fluxes flowing through the origin. 

The fully implicit discretization for a general 3-D diffusion equation with sink and source terms is 

$$
A ^ {n} \left(\frac {T ^ {n + 1} - T ^ {n}}{\triangle t}\right) = \vec {\nabla} \cdot D ^ {n} \vec {\nabla} T ^ {n + 1} + C ^ {n} T ^ {n + 1} + B ^ {n}. \tag {5.252}
$$

All coefficients of $A ^ { n }$ , $B ^ { n }$ , $C ^ { n }$ and $D ^ { n }$ are evaluated at the time level $t ^ { n }$ , and are read off directly from the multi-group radiation transport in Eq. (5.225) and the electron and ion thermal diffusions in Eq. (5.240). Physical quantities on the spherical mesh in $D E C 3 D$ is defined at the cell-centered $( r _ { i } , \theta _ { j } , \phi _ { k } )$ with the radius $r _ { i } = ( i - 1 / 2 ) \triangle r$ , the polar angle $\theta _ { j } = ( j - 1 / 2 ) \triangle \theta$ and theazimuthal angle $\phi _ { k } = ( k = 1 / 2 ) \triangle \phi$ . The mesh is discretized uniformly with $\triangle r = R / N _ { r }$ , $\triangle \theta = \pi / N _ { \theta }$ and $\bigtriangleup \phi = 2 \pi / N _ { \phi }$ . The resolution of a simulation $N _ { r } \times N _ { \theta } \times N _ { \phi }$ is defined by the number of cells in the radial direction $N _ { r }$ , the polar direction $N _ { \theta }$ and the azimuthal direction $N _ { \phi }$ . The cell center discretization is adopted to avoid the singularities at the origin and along the poles. The second-order center-in-space discretization for the Laplacian operator, 

$$
\begin{array}{l} \vec {\nabla} \cdot D ^ {n} \vec {\nabla} T ^ {n + 1} = \frac {1}{r _ {i} ^ {2}} \frac {\partial}{\partial r} \left(r ^ {2} D ^ {n} \frac {\partial T ^ {n + 1}}{\partial r}\right) + \frac {1}{r _ {i} \sin \theta_ {j}} \frac {\partial}{\partial \theta} \left(\frac {\sin \theta}{r} D ^ {n} \frac {\partial T ^ {n + 1}}{\partial \theta}\right) \\ + \frac {1}{r _ {i} \sin \theta_ {j}} \frac {\partial}{\partial \phi} \left(\frac {D ^ {n}}{r \sin \theta} \frac {\partial T ^ {n + 1}}{\partial \phi}\right) (5. 2 5 3) \\ \end{array}
$$

is given by 

$$
\vec {\nabla} \cdot D ^ {n} \vec {\nabla} T ^ {n + 1} = \bar {Z} _ {+} ^ {r} \left(T _ {i + 1} ^ {n + 1} - T _ {i} ^ {n + 1}\right) - \bar {Z} _ {-} ^ {r} \left(T _ {i} ^ {n + 1} - T _ {i - 1} ^ {n + 1}\right) + \tag {5.254}
$$

$$
\bar {Z} _ {+} ^ {\theta} \left(T _ {j + 1} ^ {n + 1} - T _ {j} ^ {n + 1}\right) - \bar {Z} _ {-} ^ {\theta} \left(T _ {j} ^ {n + 1} - T _ {j - 1} ^ {n + 1}\right) + \tag {5.255}
$$

$$
\bar {Z} _ {+} ^ {\phi} \left(T _ {k + 1} ^ {n + 1} - T _ {k} ^ {n + 1}\right) - \bar {Z} _ {-} ^ {\phi} \left(T _ {k} ^ {n + 1} - T _ {k - 1} ^ {n + 1}\right). \tag {5.256}
$$

The short-hand notations for the variables $\bar { Z } _ { \pm } ^ { r }$ , $Z _ { \pm } ^ { \theta }$ and $\bar { Z } _ { \pm } ^ { \phi }$ are defined as 

$$
\bar {Z} _ {\pm} ^ {r} = \left(r ^ {2} D\right) _ {i \pm 1 / 2} ^ {n} / \left(r _ {i} ^ {2} \triangle r ^ {2}\right), \tag {5.257}
$$

$$
\bar {Z} _ {\pm} ^ {\theta} = \left(\frac {\sin \theta}{r} D\right) _ {j \pm 1 / 2} ^ {n} / \left(r _ {i} \sin \theta_ {j} \triangle \theta^ {2}\right), \tag {5.258}
$$

$$
\bar {Z} _ {\pm} ^ {\phi} = \left(\frac {D}{r \sin \theta}\right) _ {k \pm 1 / 2} ^ {n} / \left(r _ {i} \sin \theta_ {j} \triangle \phi^ {2}\right). \tag {5.259}
$$

The same variables multiplied with $\frac { \triangle t } { A ^ { n } }$ are defined without the bar, 

$$
Z _ {\pm} ^ {r} = \frac {\triangle t}{A ^ {n}} \bar {Z} _ {\pm} ^ {r}, \quad Z _ {\pm} ^ {\theta} = \frac {\triangle t}{A ^ {n}} \bar {Z} _ {\pm} ^ {\theta}, \quad Z _ {\pm} ^ {\phi} = \frac {\triangle t}{A ^ {n}} \bar {Z} _ {\pm} ^ {\phi}. \tag {5.260}
$$

Multiply the both sides of Eq. (5.252) with $\frac { \triangle t } { A ^ { n } }$ , and arrange known variables at the previous time level $t ^ { n }$ to the left hand side, whereas unknown variables at time level $t ^ { n + 1 }$ are put at the right hand side. The following large sparse matrix resulting from the fully implicit second-order in space discretization for Eq. (5.252) is solved in HYPRE without the directional-splitting in DEC3D, and can be reduced to the simple 1-D tridiagonal matrix in Eq. (5.246) by taking $C ^ { n } = B ^ { n } = Z _ { \pm } ^ { \theta } = Z _ { \pm } ^ { \phi } = 0$ . 

$$
\begin{array}{l} - T ^ {n} - \frac {\triangle t}{A ^ {n}} B ^ {n} = Z _ {+} ^ {r} T _ {i + 1} ^ {n + 1} + Z _ {-} ^ {r} T _ {i - 1} ^ {n + 1} + Z _ {+} ^ {\theta} T _ {j + 1} ^ {n + 1} + Z _ {-} ^ {\theta} T _ {j - 1} ^ {n + 1} + Z _ {+} ^ {\phi} T _ {k + 1} ^ {n + 1} + Z _ {-} ^ {\phi} T _ {k - 1} ^ {n + 1} \\ - \left(1 + Z _ {+} ^ {r} + Z _ {-} ^ {r} + Z _ {+} ^ {\theta} + Z _ {-} ^ {\theta} + Z _ {+} ^ {\phi} + Z _ {-} ^ {\phi} - \frac {\triangle t}{A ^ {n}} C ^ {n}\right) T ^ {n + 1} \tag {5.261} \\ \end{array}
$$

All quantities $Z _ { \pm }$ evaluated at the cell interface are approximated by the simple 

# Spitzer Coulomb logrithm and without degeneracy

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/4654b51181f2e6aeb7dc28c0b3375ed58124a1a47d24d91705a7226795f18799.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/37791d4e12096c4e90f913df5e3695068bd9d1116b567a7fdc12b5866ed57182.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/5cf12c07054265c4e85095cd699b32bb8e4c0d1003b5c68969d2d0732135b72c.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/3a544319c9c42a95120b5e011736242e3828022d561660497120f6430f5167e3.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/fba592553008a506f8d56552820fc1da103af95df126635f26f16f82351514a6.jpg)



Figure 5.32: A 1-D benchmark test with $L I L A C$ for electron and ion thermal diffusions within a spherical slab using Spitzer Coulomb logarithm $\ln \Lambda _ { \mathrm { e i } } ^ { \mathrm { s p i t z e r } }$ in Eq. (5.236), and without Lee-More degeneracy correction factor $f _ { \mathrm { L M } }$ in Eq. (5.233).


mean. Throughout Eqs. (5.252)–(5.261), the spatial indices for a given variable $Q _ { i , j , k } ^ { n } \to Q ^ { n }$ are omitted and is denoted as $Q ^ { n }$ . Unless the variable is evaluated at the adjacent cell, and is denoted by $Q _ { i , j , k + 1 } ^ { n } \to Q _ { k + 1 } ^ { n }$ . 

# 5.5.4 Benchmark tests

Three tests, in the same geometry of 1-D spherical slab in radiation transport benchmark tests, were studied to benchmark the thermal diffusions using Spitzer and Lee-More coefficients of thermal conductivity, as well as the electron-and-ion equilibration. The radius of the slab is $R = 1 0 0 \ \mu \mathrm { m }$ , separated by two initial states at the radius $r _ { 0 } = 5 0 ~ \mu \mathrm { m }$ , with the initial electron temperature profile: $T _ { \mathrm { e } } ( r < r _ { 0 } ) = 5 \ \mathrm { k e V }$ and $T _ { \mathrm { e } } ( r \geq r _ { 0 } ) = 0 . 5 \ \mathrm { k e V }$ , and the initial material mass density profile: $\rho ( r < r _ { 0 } ) = 5 0 ~ \mathrm { g / c m ^ { 3 } }$ and $\rho ( r \geq r _ { 0 } ) = 1 0 0 ~ \mathrm { g / c m ^ { 3 } }$ . 

Figs. (5.32) and (5.33) are pure thermal diffusion tests, without solving the electron and ion equilibration in Eq. (5.242) and without update the fluid momentum equation in Euler equations, irregardless the building up of pressure gradient, 

# Lee-More Coulomb logrithm and with degeneracy

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/935a1b89e5a20fde4277c6d158eecb964bcd98dd16cf422fe4912cc3d32b1042.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/3462243ee7e6676e2898e145638addbe821ce3b7f7eb45a212ec389972c8a559.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/7a9759af28f8f88d5a2a7980afbbc4f3e42f947aa103eedfc06556fc49f4026f.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/ca1ea2b7f340ff8f7f62e91fd185def8838ee4ecaba1131e45aa8f560d205d0d.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/01e7c9b156fe62393b5f4d6e0a875fc395ac4d7311e0b8422418ebc0ed6ecff2.jpg)



Figure 5.33: A 1-D benchmark test with $L I L A C$ for electron and ion thermal diffusions within a spherical slab using Lee-More Coulomb logarithm $\mathrm { l n } \Lambda _ { \mathrm { e i } } ^ { \mathrm { L M } }$ in Eq. (5.234), and with Lee-More degeneracy correction factor $f _ { \mathrm { L M } }$ in Eq. (5.233).


# Lee-More Coulomb logrithm and with degeneracy and Te-Ti coupling

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/a3e8f0160b9fe5a2a2c6ec8a13f3546460bd67f877d7864852db15ed487e3597.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/9dfc6c6d466f11746f809392104d46564f56d2b5190c430cd22649698cd1bad4.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/bf5b1fe44cc6df1981e30829c2f50b21c1b516ba1d2d8fbc7c85f0d39c6a6f73.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/9ea8f2c2c6cf1f2297277ba78b3d6f378dd4ce118dda3d5df23cce700c205072.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/89431a58abf92c0f37089dcb2f1c95d6c9bef5a41bad438d02884c5775d478b0.jpg)



Figure 5.34: A 1-D benchmark test with $L I L A C$ for electron and ion thermal diffusions and equilibration using Lee-More Coulomb logarithm and degeneracy correction factor. Ions are quasi static in the plasma due to heavier mass than electrons. However, there are rapid heat exchanges between electrons and ions due to fast electron and ion collisions.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/94ff2124e624abee3a10c062f224d477977b9178810ea84869d6a9d29266ba42.jpg)



(a)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/d0c03a7a51c80a532f40ddb3d44bceefcd5b7b371d1b55e40ae431cf9b7af503.jpg)



(b)



Figure 5.35: (a) Effects of space-lagging approximations in computing the cell interface diffusion coefficient $( r ^ { 2 } \kappa ) _ { i \pm 1 / 2 } ^ { n }$ on the mass density profile at stagnation for a NIF implosion. The resolution is 200 cells in the 1-D code. (b) A benchmark test to validate the simple mean approximation to compute the cell interface coefficient at high resolution with 1000 cells in the 1-D code.


in order to maintain the fluid artificially at rest. Good agreement between LILAC solutions and 1-D simulations with Lee-More Coulomb logarithm and the degeneracy correction factor is observed in Fig. (5.33), whereas the electron temperature is higher than that for $L I L A C$ when Spitzer Coulomb logarithm is used in Fig. (5.32). Result of Fig. (5.33) confirms the accuracy of the fully implicit secondorder in space discretization of thermal diffusions. 

Figure (5.34) shows a good agreement with $L I L A C$ in a benchmark test for electron and ion thermal diffusions and equilibration, without updating the fluid momentum velocities. Because of large masses, ions are quasi static in the plasma due to much shorter thermal diffusion mean free paths than that for electrons. However, the collision rates $1 / \tau _ { \mathrm { e i } }$ defined in Eq. (5.241) between electrons and ions are high enough to exchange internal energies between ions and electrons to approach thermodynamic equilibrium. In ICF implosions, electrons and ions are not in thermodynamic equilibrium at the beginning of the deceleration phase, but they approach to about thermodynamic equilibrium at stagnation. 

Figure (5.35) compare effects of computing the cell interface diffusion coefficient term $( r ^ { 2 } \kappa ) _ { i \pm 1 / 2 } ^ { n }$ by pulling out: (1) a small part defined by $\bar { \kappa } = \kappa / T ^ { 5 / 2 }$ , 

# Lee-More Coulomb logrithm and with degeneracy and Te-Ti coupling

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/dbfe1cdebd9ecc09ebf0295add4adada6b5c063e25c3bffb62df2954711d4ac6.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/8e970ede42b55c0afdefbb8cac7712b58df74197cfcf77bff29a57a7229432b9.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/88f7a6c15fc453a28f95cf654d7ff931329792a5df4b9b79cff381633f92867d.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/1fd2ab43e82bab65c8376c691dfeafb2b15af05f0f658d08d308c343018a0290.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-09/649ac5f8-f817-45d6-84f4-bb14ff2f09a0/abf459a182e8df8d23aa1fb209de6c5db671022288f360d07b60e2140545f564.jpg)



Figure 5.36: The $D E C \mathcal { Z } D$ benchmark test with $L I L A C$ for a 1-D NIF implosion to validate hydrodynamics, electron and ion thermal diffusions and equilibration.


which is a function of Coulomb logarithm, and (2) a bigger part that contains the whole thermal conductivity $\kappa$ , from the diffusion coefficient calculation at the cell interface $r _ { i \pm 1 / 2 }$ . The dotted line is the exact solution from $L I L A C$ , the red curve is the simple mean approximation to compute the cell interface diffusion coefficient $( r ^ { 2 } \kappa ) _ { i \pm 1 / 2 } ^ { n } = [ ( r ^ { 2 } \kappa ) _ { i } ^ { n } + ( r ^ { 2 } \kappa ) _ { i \pm 1 } ^ { n } ] / 2$ , the blue curve assume a space-lagging to pull out the small part $\kappa$ , and the orange curve assume a space-lagging to pull out the whole $\kappa$ from the cell interface diffusion coefficient. The orange curve in Fig. (5.35)-(a) is clearly a wrong solution. However, the approximation of space-lagging the Coulomb logarithm in Fig. (5.35)-(a) is observed to cause $\sim 1 0$ µm error in modeling the return shock front, because of the overflow of numerical heat flux from the high-temperature hot spot into the cold shell. Figure (5.35)-(b) is a high resolution simulation for the 1-D code to validate the simple mean approximation to compute the cell interface diffusion coefficient. 

Figure (5.36) is the most up-to-date summary for the integrated performance in solving for hydrodynamics, electron and ion thermal diffusions and equilibration in DEC2D and DEC3D. Good agreement with $L I L A C$ is observed for hydro variables 

including the mass density $\rho$ , the total plasma pressure $P$ , the fluid velocity $u$ , as well as thermo variables for electron and ion temperatures $T _ { \mathrm { e } }$ , $\mathrm { \Delta } T _ { \mathrm { i } }$ . 

# 5.6 Alpha Particle Transport

The alpha particle transport is modeled by Atzeni one-group diffusion model [Atz81] to describe the straight-line motion of the slowing-down alpha particles by the electron drag force. The DT fusion reactivity $< \sigma v > _ { D T }$ is calculated from Bosch and Hale model [BH92]. Alpha particles are born within the hot spot with the initial birth energy $E _ { \alpha 0 } = 3 . 5$ MeV and the initial velocity ${ { \upsilon } _ { \alpha 0 } }$ . In the early stage, a significant fraction of alpha particle energies are transferred to electrons by the electron drag force for the temperature range $0 . 5 \leq T _ { \mathrm { e } } \leq 3 0 \mathrm { k e V }$ , characterizing by the Spitzer alpha-electron relaxation time scale $\tau _ { \alpha \mathrm { e } }$ [LS90]. 

$$
\tau_ {\alpha \mathrm {e}} = \frac {3 m _ {\alpha} (K _ {\mathrm {B}} T _ {\mathrm {e}}) ^ {3 / 2}}{8 \sqrt {2 \pi m _ {\mathrm {e}}} n _ {\mathrm {e}} \bar {Z} _ {\alpha} ^ {2} e ^ {4} \ln \Lambda_ {\alpha \mathrm {e}} ^ {\mathrm {S p i t z e r}}}, \tag {5.262}
$$

where the Spitzer alpha-electron Coulomb logarithm has the same form as $\ln \Lambda _ { \mathrm { e i } } ^ { \mathrm { s _ { f } } }$ pitzer by replacing $Z$ with $Z _ { \alpha } = 2$ , and $m _ { \alpha }$ is the mass of alpha particle. The diffusion mean-free-path because of the electron drag is $\lambda _ { \mathrm { d r a g } } = v _ { \alpha 0 } \tau _ { \alpha \mathrm { e } } / 9$ [Atz81] and the diffusion coefficient is $D _ { \alpha } = v _ { \alpha 0 } \lambda _ { \mathrm { d r a g } }$ . In the late-time with lower alpha particle kinetic energy $E _ { \alpha } < 0 . 5$ MeV, the motion of alpha particles is dominated by the transverse-scattering with ions. 

$$
\frac {\partial \varepsilon_ {\alpha}}{\partial t} = \vec {\nabla} \cdot D _ {\alpha} \vec {\nabla} \varepsilon_ {\alpha} + n _ {\mathrm {D}} n _ {T} <   \sigma v > _ {\mathrm {D T}} E _ {\alpha 0} - \frac {\varepsilon_ {\alpha}}{\tau_ {\alpha \mathrm {e}}}, \tag {5.263}
$$

where $\varepsilon _ { \alpha }$ is alpha particle energy density and $n _ { \mathrm { D } } , n _ { \mathrm { T } }$ are ion number densities of Deuterium and Tritium respectively. 

# 5.7 Conclusion

In summary, a single-fluid two-temperature 3-D parallel hydrocode DEC3D was developed. The electron pressure equation $\partial _ { t } P _ { \mathrm { e } } ^ { 3 / 5 } + \vec { \nabla } \cdot \vec { v } P _ { \mathrm { e } } ^ { 3 / 5 } = 0$ is solved as a scalar advection equation [HD06] in the single-fluid HLLC (Harten-Lax-van Leercontact) [Tor09b] approximate Riemann solver, which provides the first-order hydrodynamic solution and is followed by a boost to the third-order high resolution using the piecewise-parabolic method (PPM) [CW84], which is robust for strong shock-capturing hydrodynamic simulations. 

DEC3D is a full spherical, 3-D parallel Eulerian code. The spherical mesh provides a low level of numerical noise for simulating RT instabilities in spherical geometry. The macro-zoning technique is applied in DEC3D to map the fine mesh onto a coarser mesh, defined by the size of arc lengths of finite-volume cells comparable to the size of discretized radius, to avoid the issue of small-time-step size restricted by Courant condition [CF76]. During the macro-zoning calculation in every time step, two cells are recombined to form a coarser cell in the $\theta$ direction, followed by another recombination in the $\phi$ direction. The grid spacing in the radial direction is uniform, but it shrinks during the deceleration phase and expands in the disassembly phase. The cell-interface numerical fluxes due to the moving mesh motion in the radial direction are integrated with HLLC approximate Riemann solver by implementing the finite-volume moving-mesh algorithm [Lev02]. 

DEC3D is integrated with HYPRE [FY02], a library of high-performance preconditioners, to solve large, sparse nonsymmetric systems for implicit thermal, radiation, and alpha diffusions. The nature of strong spatial variation of opacities results in highly nonsymmetric linear systems in the implicit radiation diffusion solver; the hybrid generalized minimal residual solver [SS86] preconditioned with algebraic multigrid [HY02] is used to solve the 3-D multigroup radiation diffusion, thermal and one-group alpha diffusions without directional splitting. The 

multigrid method [TOS01] provides a fast convergent rate in diffusing high-energy photons with long mean free paths. 

# 6 Conclusion

In conclusion, an advanced modern computer code DEC3D was developed to model the hydrodynamics in the deceleration phase of inertial confinement fusion (ICF) implosions. A massive single-mode deceleration-phase simulation database for various shots was generated by DEC3D code to study 3-D hydrodynamic relations systematically. A adiabatic 3-D hot-spot model was developed to explain the impact of residual kinetic energies on yield degradations. The conservation of the hot-spot adiabatic parameter $P V ^ { \gamma }$ between 1-D and 3-D deceleration-phase implosions build up the connection to explain the role of hydrodynamics instabilities in degrading the hot-spot pressures and neutron yield productions through increasing the content of total residual kinetic energies at stagnation. An analytical model was developed to explain the impact of three-dimensional hot-spot flow asymmetry on neutron-inferred ion-temperature measurements. The velocity variance in the non-relativistic Brysk ion temperature model is decomposed into six hot-spot flow parameters that uniquely determine effects of three-dimensional hot-spot flow asymmetry on neutron-inferred ion temperatures. This analysis shows that ion-temperature measurement variations among different line of sights (LOS) are connected by a geometrical LOS matrix, so that the a full map of neutron-inferred ion temperature can be reconstructed approximately using iontemperature measurements at six different LOS, and to infer the global minimum of neutron-inferred ion temperature over the $4 \pi$ angles. The analytical method of velocity variance decomposition explains the behavior of large ion-temperature 

measurement variation for mode $\ell = 1$ , in which the flow structure is shown to exhibit highly anisotropic due to the large highly directional jet flowing through the hot spot. For modest mode $\ell = 2$ perturbations, it exhibits a large content of isotropic flow structure within the hot spot, which leads to small ion-temperature measurement variations among LOS’s and higher apparent neutron-inferred ion temperatures than thermal ion temperatures. For large mode $\ell = 2$ perturbations with the pair of RT spikes colliding with each other, it exhibits a small content of isotropic flow structure within the hot spot, which leads to the opposite trend with small ion-temperature measurement variations among LOS’s. For high performance ICF implosion experiments, the occurrence of modest mode $\ell = 2$ perturbation is more frequent than large perturbations. The impact of isotropic flow on causing DD neutron-inferred ion temperatures below than DT’s was studied using the analytical method of velocity variance decomposition. The ratio of DD to DT neutron-inferred ion temperatures was shown approaching to the limit of 0.8 in multi-mode simulations, which provides a good evidence to explain experimental observations for DD/DT ratio being below the unity. An approximate closure for non-relativistic DD and DT neutron-inferred ion temperatures was developed, which provides a systematic procedure to diagnostic the hot-spot flow isotropy and anisotropy through dedicated DD and DT ion-temperature measurements at specified LOS. 

# Bibliography



[AB04] K. Anderson and R. Betti. Laser-induced adiabat shaping by relaxation in inertial fusion implosions. Physics of Plasmas, 11(1):5–8, 2004. xviii, 5 





[ABG01] K. Anderson, R. Betti, and T. A. Gardiner. Two-dimensional computational model of energy gain in nif capsules. In Bull. Am. Phys. Soc., BAPS.2001, 2001. 131 





[AC11] B. Appelbe and J. Chittenden. The production spectrum in fusion plasmas. Plasma Physics and Controlled Fusion, 53(4):045002, 2011. 38, 39, 45, 57 





[AC14] B. Appelbe and J. Chittenden. Relativistically correct dd and dt neutron spectra. High Energy Density Physics, 11:30 – 35, 2014. 38, 39, 45, 57 





[AMJKM+05] Thomas Alan Mehlhorn, Christopher J. Kurecka, Ryan McClarren, Thomas A. Brunner, and James Holloway. Advances in radiation modeling in alegra :a final report for ldrd-67120, efficient implicit mulitgroup radiation calculations. 01 2005. 204 





[AtV04] S. Atzeni and J. Meyer ter Vehn. The Physics of Inertial Fusion: Beam Plasma Interaction, Hydrodynamics, Hot Dense Matter. International Series of Monographs on Physics. Oxford, 2004. 1, 3, 10, 16, 38, 57, 134, 208 





[Atz81] S. Atzeni. A diffusive model for alpha particle energy transport in a laser plasma. NUOVO CIMENTO, 64(2), 1981. 222 





[BB73] Jay P Boris and David L Book. Flux-corrected transport. i. shasta, a fluid transport algorithm that works. Journal of Computational Physics, 11(1):38 – 69, 1973. 155 





[BBC+97] T.R Boehly, D.L Brown, R.S Craxton, R.L Keck, J.P Knauer, J.H Kelly, T.J Kessler, S.A Kumpan, S.J Loucks, S.A Letzring, F.J Marshall, R.L McCrory, S.F.B Morse, W Seka, J.M Soures, and C.P Verdon. Initial performance results of the omega laser system. Optics Communications, 133(1):495 – 506, 1997. 2 





[BBSW17] A. Bose, R. Betti, D. Shvarts, and K. M. Woo. The physics of long- and intermediate-wavelength asymmetries of the hot spot: Compression hydrodynamics and energetics. Phys. Plasmas, 24(10):102704, 2017. 22, 23, 35, 37, 38 





[BCBW16] R. Betti, A.R. Christopherson, A. Bose, and K.M. Woo. Alpha heating and burning plasmas in inertial confinement fusion. Journal of Physics: Conference Series, 717(1):012007, 2016. 31, 80 





[BH92] H.-S. Bosch and G.M. Hale. Improved formulas for fusion crosssections and thermal reactivities. Nuclear Fusion, 32(4), 1992. 222 





[Bru02] T.A. Brunner. Forms of approximate radiation transport, 2002. 189, 190, 195, 205 





[Bry73] H. Brysk. Fusion neutron energies and spectra. Plasma Phys., 15:611–617, 1973. 11, 38, 39, 41, 42, 45, 56, 60, 82 





[BUL+01] R. Betti, M. Umansky, V. Lobatchev, V. N. Goncharov, and R. L. McCrory. Hot-spot dynamics and deceleration-phase rayleigh– 





taylor instability of imploding inertial confinement fusion capsules. Phys. Plasmas, 8(12):5257–5267, 2001. 27, 31 





[BW80] B. Born and E. Wolf. Principles of Optics. Pergamon, 1980. 192 





[BWB+16] A. Bose, K. M. Woo, R. Betti, E. M. Campbell, D. Mangino, A. R. Christopherson, R. L. McCrory, R. Nora, S. P. Regan, V. N. Goncharov, T. C. Sangster, C. J. Forrest, J. Frenje, M. Gatu Johnson, V. Yu Glebov, J. P. Knauer, F. J. Marshall, C. Stoeckl, and W. Theobald. Core conditions for alpha heating attained in directdrive inertial confinement fusion. Phys. Rev. E, 94:011201, 2016. 13 





[BWNB15] A. Bose, K. M. Woo, R. Nora, and R. Betti. Hydrodynamic scaling of the deceleration-phase rayleigh taylor instability. Physics of Plasmas, 22(7):072702, 2015. 201 





[CAB+15] R. S. Craxton, K. S. Anderson, T. R. Boehly, V. N. Goncharov, D. R. Harding, J. P. Knauer, R. L. McCrory, P. W. McKenty, D. D. Meyerhofer, J. F. Myatt, A. J. Schmitt, J. D. Sethian, R. W. Short, S. Skupsky, W. Theobald, W. L. Kruer, K. Tanaka, R. Betti, T. J. B. Collins, J. A. Delettrez, S. X. Hu, J. A. Marozas, A. V. Maximov, D. T. Michel, P. B. Radha, S. P. Regan, T. C. Sangster, W. Seka, A. A. Solodov, J. M. Soures, C. Stoeckl, and J. D. Zuegel. Direct-drive inertial confinement fusion: A review. Physics of Plasmas, 22(11):110501, 2015. 2 





[Cas04] J. Castor. Radiation Hydrodynamics. Cambridge, 2004. 200 





[CBS+10] P.-Y. Chang, R. Betti, B. K. Spears, K. S. Anderson, J. Edwards, M. Fatenejad, J. D. Lindl, R. L. McCrory, R. Nora, and D. Shvarts. 





Generalized measurable ignition criterion for inertial confinement fusion. Phys. Rev. Lett., 104:135002, 2010. 34 





[CF76] R. Courant and K. O. Friedrichs. Supersonic Flow and Shock Waves. Springer, 1976. 144, 154, 155, 223 





[CG68] J. P. Cox and R. T. Giuli. Principles of Stellar Structure. Gordon and Breach, Science Publishers, 1968. 191, 192 





[Cla87] B. A. Clark. Computing multigroup radiation integrals using polylogarithm-based methods. Journal of Computational Physics, 70:311–329, 1987. 202 





[CLS94] W. Cropp, E. Lusk, and A. Skjellum. Using MPI Portable Parallel Prpgramming with the Message-Passing-Interface. The MIT Press, 1994. xxxiv, 140, 143, 207 





[CMD15] Duc Cao, Gregory Moses, and Jacques Delettrez. Improved nonlocal electron thermal transport model for two-dimensional radiation hydrodynamics simulations. Physics of Plasmas, 22(8):082308, 2015. 14 





[Coo66] J. Cooper. Plasma spectroscopy. Reports on Progress in Physics, 29(1):35, 1966. 193 





[CW84] P. Colella and P. R. Woodward. The piecewise parabolic method (ppm) for gas-dynamical simulations. Journal of Computational Physics, 54(1), 1984. 131, 155, 223 





[Dav84] S. F. Davis. Tvd finite difference schemes and artificial viscosity. (NASA-CR-172373, ICASE-84-20, NAS 1.26:172373), 1984. 155 





[Dav87] S. Davis. A simplified tvd finite difference scheme via artificial viscosity. SIAM Journal on Scientific and Statistical Computing, 8(1):1–18, 1987. 155 





[Dav88] S. F. Davis. Simplified second-order godunov-type methods. SIAM J. Sci. Stat. Comput., 9(3):445–473, May 1988. 174 





[Dav15] P. A. Davidson. Turbulence An Introduction For Scientists And Engineers. Oxford University Press, 2015. 18 





[DER+87] J. Delettrez, R. Epstein, M. C. Richardson, P. A. Jaanimagi, and B. L. Henke. Effect of laser illumination nonuniformity on the analysis of time-resolved x-ray measurements in uv spherical transport experiments. Phys. Rev. A, 36:3926, 1987. 14 





[DG90] Jill P. Dahlburg and John H. Gardner. Ablative rayleigh-taylor instability in three dimensions. Phys. Rev. A, 41:5695–5698, 1990. 20 





[DGDH93] J. P. Dahlburg, J. H. Gardner, G. D. Doolen, and S. W. Haan. The effect of shape in the three-dimensional ablative rayleigh-taylor instability. i: Single-mode perturbations. Physics of Fluids B, 5(2):571–584, 1993. 20 





[Dir01] Paul A. M. Dirac. Lectures on Quantum Mechanics. DOVER PUB-LICATIONS, INC., 2001. 191 





[Eva84] D.J. Evans. Parallel s.o.r. iterative methods. Parallel Computing, 1(1):3 – 18, 1984. 202 





[FL03] R. Fazio and R. J. LeVeque. Moving-mesh methods for onedimensional hyperbolic problems using {CLAWPACK}. Computers & Mathematics with Applications, 45(1-3):273 – 298, 2003. 147 





[FY02] R. D. Falgout and U. M. Yang. hypre: A Library of High Performance Preconditioners, volume 2331, page 632. Springer, Berlin Heidelberg, 2002. 223 





[GDF+14] Jianfa Gu, Zhensheng Dai, Zhengfeng Fan, Shiyang Zou, Wenhua Ye, Wenbing Pei, and Shaoping Zhu. A new metric of the lowmode asymmetry for ignition target designs. Physics of Plasmas, 21(1):012704, 2014. 12 





[GJCF+13] M. Gatu Johnson, D. T. Casey, J. A. Frenje, C.-K. Li, F. H. S´eguin, R. D. Petrasso, R. Ashabranner, R. Bionta, S. LePape, M. McKernan, A. Mackinnon, J. D. Kilkenny, J. Knauer, and T. C. Sangster. Measurements of collective fuel velocities in deuterium-tritium exploding pusher and cryogenically layered deuterium-tritium implosions on the nif. Physics of Plasmas, 20(4):042707, 2013. 38, 122 





[GJKC+16] M. Gatu Johnson, J. P. Knauer, C. J. Cerjan, M. J. Eckart, G. P. Grim, E. P. Hartouni, R. Hatarik, J. D. Kilkenny, D. H. Munro, D. B. Sayre, B. K. Spears, R. M. Bionta, E. J. Bond, J. A. Caggiano, D. Callahan, D. T. Casey, T. D¨oppner, J. A. Frenje, V. Yu. Glebov, O. Hurricane, A. Kritcher, S. LePape, T. Ma, A. Mackinnon, N. Meezan, P. Patel, R. D. Petrasso, J. E. Ralph, P. T. Springer, and C. B. Yeamans. Indications of flow near maximum compression in layered deuterium-tritium implosions at the national ignition facility. Phys. Rev. E, 94:021202, 2016. 44, 84 





[God59] S. K. Godunov. A difference method for numerical calculation of discontinuous solutions of the equations of hydrodynamics. Mat. Sb. (N.S.), pages 271–306, 1959. 154 





[Gri62] Hans R. Griem. High-density corrections in plasma spectroscopy. Phys. Rev., 128:997–1003, Nov 1962. 193 





[GSB+14] V. N. Goncharov, T. C. Sangster, R. Betti, T. R. Boehly, M. J. Bonino, T. J. B. Collins, R. S. Craxton, J. A. Delettrez, D. H. Edgell, R. Epstein, R. K. Follett, C. J. Forrest, D. H. Froula, V. Yu. Glebov, D. R. Harding, R. J. Henchen, S. X. Hu, I. V. Igumenshchev, R. Janezic, J. H. Kelly, T. J. Kessler, T. Z. Kosc, S. J. Loucks, J. A. Marozas, F. J. Marshall, A. V. Maximov, R. L. Mc-Crory, P. W. McKenty, D. D. Meyerhofer, D. T. Michel, J. F. Myatt, R. Nora, P. B. Radha, S. P. Regan, W. Seka, W. T. Shmayda, R. W. Short, A. Shvydky, S. Skupsky, C. Stoeckl, B. Yaakobi, J. A. Frenje, M. Gatu-Johnson, R. D. Petrasso, and D. T. Casey. Improving the hot-spot pressure and demonstrating ignition hydrodynamic equivalence in cryogenic deuterium-tritium implosions on omega. Phys. Plasmas, 21(5):056315, 2014. 46, 83 





[Har83] Ami Harten. High resolution schemes for hyperbolic conservation laws. Journal of Computational Physics, 49(3):357 – 393, 1983. 155 





[Har89] Ami Harten. Eno schemes with subcell resolution. Journal of Computational Physics, 83(1):148 – 184, 1989. 155 





[Har11] R. C. Harwood. Operator Splitting Method and Application for Semilinear Parabolic Partial Differential Equations. PhD dissertation, Washington State University, Department of Mathematics, May 2011. 135, 203 





[HCB+14] S. X. Hu, L. A. Collins, T. R. Boehly, J. D. Kress, V. N. Goncharov, and S. Skupsky. First-principles thermal conductivity of 





warm-dense deuterium plasmas for inertial confinement fusion applications. Phys. Rev. E, 89:043105, 2014. 211 





[HD06] J. A. F. Hittinger and M. R. Dorr. Improving the capabilities of a continuum laser plasma interaction code. J. Phys.: Conf. Ser., 46(1):422, 2006. 223 





[HLL83] A. Harten, P. Lax, and B. Leer. On upstream differencing and godunov-type schemes for hyperbolic conservation laws. SIAM Review, 25(1):35–61, 1983. 155 





[HMGS11] S. X. Hu, B. Militzer, V. N. Goncharov, and S. Skupsky. Firstprinciples equation-of-state table of deuterium for inertial confinement fusion applications. Phys. Rev. B, 84:224109, 2011. 14 





[HMMA77] W. F. Huebner, A. L. Merts, N. H. Magee, Jr., and M. F. Argo. Astrophysical Opacity Library, 1977. 200, 203 





[HMS05] T. A. Heltemes, G. A. Moses, and J. F. Santarius. Analysis of an improved fusion reaction rate model for use in fusion plasma simulations. In Fusion Technology Institute Report UWFDM-1268. 2005. xxxiv, 134 





[HY02] V. E. Henson and U. M. Yang. Boomeramg: A parallel algebraic multigrid solver and preconditioner. Applied Numerical Mathematics, 41(1):155 – 177, 2002. 223 





[IMS+17] I. V. Igumenshchev, D. T. Michel, R. C. Shah, E. M. Campbell, R. Epstein, C. J. Forrest, V. Yu. Glebov, V. N. Goncharov, J. P. Knauer, F. J. Marshall, R. L. McCrory, S. P. Regan, T. C. Sangster, C. Stoeckl, A. J. Schmitt, and S. Obenschain. Three-dimensional 





hydrodynamic simulations of omega implosions. Physics of Plasmas, 24(5):056307, 2017. 41 





[JFC+12] M. Gatu Johnson, J. A. Frenje, D. T. Casey, C. K. Li, F. H. S´eguin, R. Petrasso, R. Ashabranner, R. M. Bionta, D. L. Bleuel, E. J. Bond, J. A. Caggiano, A. Carpenter, C. J. Cerjan, T. J. Clancy, T. Doeppner, M. J. Eckart, M. J. Edwards, S. Friedrich, S. H. Glenzer, S. W. Haan, E. P. Hartouni, R. Hatarik, S. P. Hatchett, O. S. Jones, G. Kyrala, S. Le Pape, R. A. Lerche, O. L. Landen, T. Ma, A. J. MacKinnon, M. A. McKernan, M. J. Moran, E. Moses, D. H. Munro, J. McNaney, H. S. Park, J. Ralph, B. Remington, J. R. Rygg, S. M. Sepke, V. Smalyuk, B. Spears, P. T. Springer, C. B. Yeamans, M. Farrell, D. Jasion, J. D. Kilkenny, A. Nikroo, R. Paguio, J. P. Knauer, V. Yu Glebov, T. C. Sangster, R. Betti, C. Stoeckl, J. Magoon, M. J. Shoup, G. P. Grim, J. Kline, G. L. Morgan, T. J. Murphy, R. J. Leeper, C. L. Ruiz, G. W. Cooper, and A. J. Nelson. Neutron spectrometry-an essential tool for diagnosing implosions at the national ignition facility (invited). Review of Scientific Instruments, 83(10):10D308, 2012. 7 





[Kai00] Thomas B. Kaiser. Laser ray tracing and power deposition on an unstructured three-dimensional grid. Phys. Rev. E, 61:895–905, Jan 2000. 3 





[Ker76] D. S. Kershaw. Flux limiting nature’s own way. UCRL-78378, 1976. 198 





[Kru88] W. L. Kruer. The physics of laser plasma interactions. Addison-Wesley, 1988. 3 





[KS01] R. Kishony and D. Shvarts. Ignition condition and gain prediction for perturbed inertial confinement fusion targets. Physics of Plasmas, 8(11):4925–4936, 2001. 22, 34, 35, 37 





[KTB+14] A. L. Kritcher, R. Town, D. Bradley, D. Clark, B. Spears, O. Jones, S. Haan, P. T. Springer, J. Lindl, R. H. H. Scott, D. Callahan, M. J. Edwards, and O. L. Landen. Metrics for long wavelength asymmetries in inertial confinement fusion implosions on the national ignition facility. Physics of Plasmas, 21(4):042708, 2014. 12, 38, 84 





[Lev02] R. J. Leveque. Finite Volume Methods for Hyperbolic Problems. Cambridge Texts in Applied Mathematics. Cambridge University Press, 2002. 149, 151, 153, 156, 157, 173, 223 





[Lin95] John Lindl. Development of the indirect-drive approach to inertial confinement fusion and the target physics basis for ignition and gain. Physics of Plasmas, 2(11):3933–4024, 1995. 1 





[LM84] Y. T. Lee and R. M. More. An electron conductivity model for dense plasmas. The Physics of Fluids, 27(5):1273–1286, 1984. 211, 212 





[LS90] Jr. L. Spitzer. Physics of Fully Ionized Gases. Dover Publications Inc., Mineola, NY, 2nd rev. ed. edition, 1990. 222 





[Mac] Robert W. MacCormack. The Effect of Viscosity in Hypervelocity Impact Cratering, pages 27–43. 154 





[MFH+17] D. H. Munro, J. E. Field, R. Hatarik, J. L. Peterson, E. P. Hartouni, B. K. Spears, and J. D. Kilkenny. Impact of temperature-velocity 





distribution on fusion neutron peak shape. Physics of Plasmas, 24(5):056301, 2017. 38, 39, 45, 56, 59, 60, 75, 84 





[MFS+17] J. F. Myatt, R. K. Follett, J. G. Shaw, D. H. Edgell, D. H. Froula, I. V. Igumenshchev, and V. N. Goncharov. A wave-based model for cross-beam energy transfer in direct-drive inertial confinement fusion. Physics of Plasmas, 24(5):056308, 2017. 3 





[MGF+18] O. M. Mannion, V. Yu. Glebov, C. J. Forrest, J. P. Knauer, V. N. Goncharov, S. P. Regan, T. C. Sangster, C. Stoeckl, and M. Gatu Johnson. Calibration of a neutron time-of-flight detector with a rapid instrument response function for measurements of bulk fluid motion on omega. Review of Scientific Instruments, 89(10):10I131, 2018. 38, 47 





[MS10] F. Mandl and G. Shaw. Quantum Field Theory. WILEY, 2010. 191 





[MUB08] A. Mignone, M. Ugliano, and G. Bodo. A five-wave hll riemann solver for relativistic mhd. arXiv:0811.1483, 2008. 155 





[Mun16] David H. Munro. Interpreting inertial fusion neutron spectra. Nuclear Fusion, 56(3):036001, 2016. 11, 38, 39, 45, 56, 59, 60, 75, 82, 84 





[Mur14] T. J. Murphy. The effect of turbulent kinetic energy on inferred ion temperature from neutron spectra. Physics of Plasmas, 21(7):072701, 2014. 8, 11, 38, 39, 82, 84, 90, 98 





[OAH00] Gordon L. Olson, Lawrence H. Auer, and Michael L. Hall. Diffusion, p1, and other approximate forms of radiation transport. Jour-





nal of Quantitative Spectroscopy and Radiative Transfer, 64(6):619 – 634, 2000. 198 





[OAK+01] D. Oron, L. Arazi, D. Kartoon, A. Rikanati, U. Alon, and D. Shvarts. Dimensionality dependence of the rayleigh-taylor and richtmyer-meshkov instability late-time scaling laws. Physics of Plasmas, 8(6):2883–2889, 2001. 49 





[PMTT14] J. L. Peterson, P. Michel, C. A. Thomas, and R. P. J. Town. The impact of laser plasma interactions on three-dimensional drive symmetry in inertial confinement fusion implosions. Physics of Plasmas, 21(7):072712, 2014. 41 





[Pom73] G. C. Pomraning. The Equations of Radiation Hydrodynamics. Pergamon Press, 1973. 189, 195, 200 





[Ray83] Lord Rayleigh. The form of standing waves on the surface of running water. Proc. London Math. Soc., 14:170, 1883. 5, 6 





[RGI+16] S. P. Regan, V. N. Goncharov, I. V. Igumenshchev, T. C. Sangster, R. Betti, A. Bose, T. R. Boehly, M. J. Bonino, E. M. Campbell, D. Cao, T. J. B. Collins, R. S. Craxton, A. K. Davis, J. A. Delettrez, D. H. Edgell, R. Epstein, C. J. Forrest, J. A. Frenje, D. H. Froula, M. Gatu Johnson, V. Yu. Glebov, D. R. Harding, M. Hohenberger, S. X. Hu, D. Jacobs-Perkins, R. Janezic, M. Karasik, R. L. Keck, J. H. Kelly, T. J. Kessler, J. P. Knauer, T. Z. Kosc, S. J. Loucks, J. A. Marozas, F. J. Marshall, R. L. McCrory, P. W. McKenty, D. D. Meyerhofer, D. T. Michel, J. F. Myatt, S. P. Obenschain, R. D. Petrasso, P. B. Radha, B. Rice, M. J. Rosenberg, A. J. Schmitt, M. J. Schmitt, W. Seka, W. T. Shmayda, M. J. Shoup, A. Shvydky, S. Skupsky, A. A. Solodov, C. Stoeckl, W. Theobald, 





J. Ulreich, M. D. Wittman, K. M. Woo, B. Yaakobi, and J. D. Zuegel. Demonstration of fuel hot-spot pressure in excess of 50 gbar for direct-drive, layered deuterium-tritium implosions on omega. Phys. Rev. Lett., 117:025001, 2016. 13 





[Roe84] P. L. Roe. Generalized formulation of tvd lax-wendroff schemes. (NASA-CR-172478, ICASE-84-53, NAS 1.26:172478), 1984. 155, 157 





[Roe86] P. L. Roe. Characteristic-based schemes for the euler equations. Annual Review of Fluid Mechanics, 18(1):337–365, 1986. 155 





[Roe97] P. L. Roe. Approximate riemann solvers, parameter vectors, and difference schemes. Journal of Computational Physics, 135(2):250 – 258, 1997. 155 





[Saf92] P. G. Saffman. Vortex Dynamics. Cambridge University Press, 1992. 18 





[SB05a] J. Sanz and R. Betti. Analytical model of the ablative rayleigh– taylor instability in the deceleration phase. Phys. Plasmas, 12(4):042704, 2005. 12, 38, 80 





[SB05b] J. Sanz and R. Betti. Analytical model of the ablative rayleightaylor instability in the deceleration phase. Physics of Plasmas, 12(4):042704, 2005. 17 





[SBRR04] J. Sanz, R. Betti, R. Ramis, and J. Ram´ırez. Nonlinear theory of the ablative rayleigh-taylor instability. Plasma Physics and Controlled Fusion, 46(12B):B367, 2004. 20, 49 





[SCB+13] R. H. H. Scott, D. S. Clark, D. K. Bradley, D. A. Callahan, M. J. Edwards, S. W. Haan, O. S. Jones, B. K. Spears, M. M. Marinak, 





R. P. J. Town, P. A. Norreys, and L. J. Suter. Numerical modeling of the sensitivity of x-ray driven implosions to low-mode flux asymmetries. Phys. Rev. Lett., 110:075001, 2013. 12 





[Sch65] Martin Schwarzschild. Structure and evolution of the stars. Dover, 1965. 191 





[SDMV81] D. Shvarts, J. Delettrez, R. L. McCrory, and C. P. Verdon. Selfconsistent reduction of the spitzer-harm electron thermal heat flux in steep temperature gradients in laser-produced plasmas. Phys. Rev. Lett., 47:247–250, 1981. 4, 198 





[SEH+14] Brian K. Spears, M. J. Edwards, S. Hatchett, J. Kilkenny, J. Knauer, A. Kritcher, J. Lindl, D. Munro, P. Patel, H. F. Robey, and R. P. J. Town. Mode 1 drive asymmetry in inertial confinement fusion implosions on the national ignition facility. Physics of Plasmas, 21(4):042702, 2014. 12, 38, 41 





[SGC+05] J. Sanz, J. Garnier, C. Cherfils, B. Canaud, L. Masse, and M. Temporal. Self-consistent analysis of the hot spot dynamics for inertial confinement fusion capsules. Phys. Plasmas, 12(11):112702, 2005. 12, 38 





[SGE+12] Brian K. Spears, S. Glenzer, M. J. Edwards, S. Brandon, D. Clark, R. Town, C. Cerjan, R. Dylla-Spears, E. Mapoles, D. Munro, J. Salmonson, S. Sepke, S. Weber, S. Hatchett, S. Haan, P. Springer, E. Moses, J. Kline, G. Kyrala, and D. Wilson. Performance metrics for inertial confinement fusion implosions: Aspects of the technical framework for measuring progress in the national ignition campaign. Physics of Plasmas, 19(5):056316, 2012. 84 





[SGT+08] J. M. Stone, T. A. Gardiner, P. Teuben, J. F. Hawley, and J. B. Simon. Athena: A new code for astrophysical mhd. Astrophysical Journal, 178, 2008. 187 





[SH53] Lyman Spitzer and Richard H¨arm. Transport phenomena in a completely ionized gas. Phys. Rev., 89:977–981, 1953. 198, 210, 211 





[Sha94] R. Shankar. Principles of Quantum Mechanics. Springer, 1994. 136 





[SNB00] G. P. Schurtz, Ph. D. Nicola, and M. Busquet. A nonlocal electron conduction model for multidimensional radiation hydrodynamics codes. Physics of Plasmas, 7(10):4238–4249, 2000. 3 





[Sod78] Gary A. Sod. A survey of several finite difference methods for systems of nonlinear hyperbolic conservation laws. Journal of Computational Physics, 27(1):1–31, 1978. 131, 153, 154 





[Spi90] Jr. Lyman Spitzer. Physics of Fully Ionized Gases. Dover, 1990. 210, 211, 212 





[SS86] Y. Saad and M. H. Schultz. Gmres: A generalized minimal residual algorithm for solving nonsymmetric linear systems. SIAM Journal on Scientific and Statistical Computing, 7:856 – 869, 1986. 202, 223 





[Swe84] P. Sweby. High resolution schemes using flux limiters for hyperbolic conservation laws. SIAM Journal on Numerical Analysis, 21(5):995–1011, 1984. 157 





[Tor09a] E. F. Toro. Riemann Solvers and Numerical Methods for Fluid Dynamics. Springer, 2009. 155, 157, 173, 174, 176 





[Tor09b] E. F. Toro. Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction. Springer, Berlin, 2009. 223 





[TOS01] U. Trottenberg, C. W. Oosterlee, and A. Schuller. Multigrid. Academic Press, Inc., San Diego, 2001. 145, 202, 224 





[vL79] B. v. Leer. Towards the ultimate conservative difference scheme. v. a second-order sequel to godunov’s method. Journal of Computational Physics, 32(1), 1979. 131, 155, 158, 181 





[VR50] J. VonNeumann and R. D. Richtmyer. A method for the numerical calculation of hydrodynamic shocks. Journal of Applied Physics, 21(3):232–237, 1950. 153 





[VRMV92] G. Velarde, Y. Ronen, and J. M. Martinez-Val. Nuclear Fusion by Inertial Confinement A Comprehensive Treatise. CRC Press, 1992. 1 





[WBB+15] K. M. Woo, R. Betti, A. Bose, R. Epstein, J.A. Delettrez, K.S. Anderson, R. Yan, P.-Y. Chang, D. Jonathan, and M. Charissis. Three-dimensional simulations of the deceleration phase of inertial fusion implosions using dec3d. In Bull. Am. Phys. Soc., BAPS.2015.DPP.GO5.3, volume 60, 2015. xvii, xxxiv, 140, 141, 142, 188, 201 





[WBS+18a] K. M. Woo, R. Betti, D. Shvarts, A. Bose, D. Patel, R. Yan, P.-Y. Chang, O. M. Mannion, R. Epstein, J. A. Delettrez, M. Charissis, K. S. Anderson, P. B. Radha, A. Shvydky, I. V. Igumenshchev, V. Gopalaswamy, A. R. Christopherson, J. Sanz, and H. Aluie. Effects of residual kinetic energy on yield degradation and ion temperature asymmetries in inertial confinement fusion implosions. Physics of Plasmas, 25(5):052704, 2018. 10, 12, 38, 48, 51, 81, 103, 123, 126 





[WBS+18b] K. M. Woo, R. Betti, D. Shvarts, O. M. Mannion, D. Patel, V. N. Goncharov, K. S. Anderson, P. B. Radha, J. P. Knauer, A. Bose, V. Gopalaswamy, A. R. Christopherson, E. M. Campbell, J. Sanz, and H. Aluie. Impact of three-dimensional hot-spot flow asymmetry on ion-temperature measurements in inertial confinement fusion experiments. Physics of Plasmas, 25(10):102710, 2018. 10, 11, 38, 84, 85, 93, 103, 107, 108, 109, 123 





[WBY+16] K. M. Woo, R. Betti, R. Yan, H. Aluie, A. Bose, D. X. Zhao, and V. Gopalaswamy. Three-dimensional study of yield degradation for direct-drive inertial confinement. In Bull. Am. Phys. Soc., BAPS.2016.DPP.TO5.15, volume 61, 2016. 12, 27, 28 





[WCC+15] C. R. Weber, D. S. Clark, A. W. Cook, D. C. Eder, S. W. Haan, B. A. Hammel, D. E. Hinkel, O. S. Jones, M. M. Marinak, J. L. Milovich, P. K. Patel, H. F. Robey, J. D. Salmonson, S. M. Sepke, and C. A. Thomas. Three-dimensional hydrodynamics of the deceleration stage in inertial confinement fusion. Physics of Plasmas, 22(3):032702, 2015. 38 





[WED+14] K. M. Woo, R. Epstein, J. A. Delettrez, A. Bose, R. Betti, and K. S. Anderson. The three-dimensional hydrocode dec3d with multigroup radiation transport. In Bull. Am. Phys. Soc., BAPS.2014.DPP.UP8.87, volume 59, 2014. 188, 201 





[Whi74] G. B. Whitham. Linear and Nonlinear Waves. Wiley, 1974. 154 





[WRF18] F. Weilacher, P. B. Radha, and C. Forrest. Three-dimensional modeling of the neutron spectrum to infer plasma conditions in cryogenic inertial confinement fusion implosions. Physics of Plasmas, 25(4):042704, 2018. xxii, 42, 47, 48, 83, 128 





[Yee85] H.C. Yee. Generalized formulation of a class of explicit and implicit TVD schemes. 08 1985. 155 





[Yee87] H.C. Yee. Construction of explicit and implicit symmetric tvd schemes and their applications. Journal of Computational Physics, 68(1):151 – 179, 1987. 155 





[Yee94] H.C. Yee. A class of high-resolution explicit and implicit shockcapturing methods. 02 1994. 155 





[Yee97] H.C. Yee. Explicit and implicit multidimensional compact highresolution shock-capturing methods:formulation. Journal of Computational Physics, 131(1):216 – 232, 1997. 155 





[ZB07] C. D. Zhou and R. Betti. Hydrodynamic relations for direct-drive fast-ignition and conventional inertial confinement fusion implosions. Phys. Plasmas, 14(7):072703, 2007. 3, 32 





[ZR02] Ya. B. Zeldovich and Yu. P. Raizer. Physics of Shock Waves and High-Temperature Hydrodynamic Phenomena. DOVER PUBLI-CATIONS, 2002. 3, 180, 189, 192, 196 

