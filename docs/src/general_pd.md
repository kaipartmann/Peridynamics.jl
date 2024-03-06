# Basics

Peridynamics is a nonlocal formulation of continuum mechanics whereby in the discrete version a body is represented by a group of interconnected material points.
The basic idea is that each point interacts with any other point within the distance of the horizon $\delta$. 
These interactions called *bonds* provoke forces between points which can be calculated and result in internal forces.

```@raw html
<img src="/assets/PDBody.png?raw=true" width="400"/>
```

Each point $i$ of the body $\mathcal{B}_0$ in initial configuration is characterized by the vector $\boldsymbol{X}^{i}$. The neighborhood $\mathcal{H}^i$ of $i$ contains all points $j$ that are within the distance of horizon $\delta$ to $i$:

```math
    \mathcal{H}^i = \{ \boldsymbol{X}^j \in \mathcal{B}_0 \; | \; 0 < \left| \boldsymbol{X}^j - \boldsymbol{X}^i \right| \leq \delta \} \quad \forall \; \boldsymbol{X}^i \in \mathcal{B}_0 \; .
```

All of those combinations of points $i$ and $j$ form bonds, characterized by 
$\boldsymbol{\Delta X}^{ij}= \boldsymbol{X}^j - \boldsymbol{X}^i$.

In current configuration $\mathcal{B}_t$, lowercase letters are used to describe points and bonds: $\boldsymbol{x} ^ {i} $ and $\boldsymbol{\Delta x}^{ij}$.

The peridynamic equation of motion is established for each material point $i$ and equals [Silling2000](@cite)
```math
\rho \ddot{\boldsymbol{u}} (\boldsymbol{X}^{i} , t) = \boldsymbol{b}^{\mathrm{int}} (\boldsymbol{X}^{i} , t) + \boldsymbol{b}^{\mathrm{ext}} (\boldsymbol{X}^{i} , t) 
    \qquad \forall \; \boldsymbol{X} ^ {i} \in \mathcal{B}_0 \; , \; t \geq 0 \; .
```

Assuming the external forces $\boldsymbol{b}^{\mathrm{ext}}$ are known, there are different peridynamic approaches to calculate the internal forces $\boldsymbol{b}^{\mathrm{int}}$, for example the [bond-based formulation](@ref bondbased) or the [continuum-kinematics-inspired formulation](@ref continuumbased).