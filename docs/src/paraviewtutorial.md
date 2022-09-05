# Visualization with ParaView

The following section explains a few visualization basics with [ParaView](https://www.paraview.org) (ParaView 5.10.1 on macOS Monterey).

## Basics
### 1. Load the results
Load the resulting `.vtu` files into ParaView with **File** $\rightarrow$ **Open**.
```@raw html
<img src="https://github.com/kfrb/Peridynamics.jl/blob/main/docs/src/assets/ParaViewtutorial01.png?raw=true" width="600" />
```

### 2. Select the time array
Select the time array and all the parameters you want to analyze and then **Apply**.
```@raw html
<img src="https://github.com/kfrb/Peridynamics.jl/blob/main/docs/src/assets/ParaViewtutorial02.png?raw=true" width="1000" />
```

### 3. Representation and coloring
Change the representation to **Points** and then choose the coloring parameter.
```@raw html
<img src="https://github.com/kfrb/Peridynamics.jl/blob/main/docs/src/assets/ParaViewtutorial03.png?raw=true" width="1000" />
```

### 4. Point styling
Activate the setting to render points as spheres and set an appropriate point size.
```@raw html
<img src="https://github.com/kfrb/Peridynamics.jl/blob/main/docs/src/assets/ParaViewtutorial04.png?raw=true" width="1000" />
```

### 5. Legend styling
By default, the legend limits are set for the current range of the coloring parameter.
In this example, for the initial time step all damage values are zero so strange legend limits appear.
```@raw html
<img src="https://github.com/kfrb/Peridynamics.jl/blob/main/docs/src/assets/ParaViewtutorial05.png?raw=true" width="1000" />
```

### 6. Save animation
To generate a animation, use **File** $\rightarrow$ **Save Animation...** and follow the instructions.

## Additional Learning Resources
ParaView has a great documentation and a lot of resources for learning.
For example, see:

- [ParaView Documentation](https://docs.paraview.org/en/latest/)
- [ParaView User's Guide](https://docs.paraview.org/en/latest/UsersGuide/index.html)