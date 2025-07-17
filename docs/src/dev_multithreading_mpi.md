# Multithreading & MPI  

!!! warning "Draft & Work in Progress"
    This documentation is a draft and work in progress. It will be extended and improved in the future.

Multithreading and MPI are implemented based on a similar approach. The body is split into several chunks, that each have their own process.
```@raw html
<img src="https://github.com/user-attachments/assets/86caad93-8a56-495b-822a-ffccc5c43c9e" width="650"/>
```
There are two types of points. Each chunk consists of the local points, which lay inside the specific body chunk, as well as halo points, which lay on the surface of neighboring chunks and interact with local points of the considered chunk, because they lay inside their horizons. 

This means that every point of a body exists locally in the body chunk that it originally lies in. Additionally one or more copies of the point can exist as halo points in other chunks.
```@raw html
<img src="https://github.com/user-attachments/assets/79340f51-cb62-42b2-8f1b-fd474f7fe8fb" width="350"/>
```
Since there are points that exist in multiple chunks (in one chunk as local points and possibly in one or more chunks as halo points) information has to be exchanged between these copies of specific points.  
Two types of exchanges are of importance for both methods:
- local-to-halo exchange: Data from calculations in the chunk where the point is local is transferred to the halo versions of this point which exist in other chunks.
```@raw html
<img src="https://github.com/user-attachments/assets/34434920-a758-438a-a463-a259397f1bba" width="350"/>
```
- halo-to-local exchange: Data from calculations of the one or multiple halo versions of a point is transferred to the local version of the point.
```@raw html
<img src="https://github.com/user-attachments/assets/6265e5b8-bff5-4aa3-83b7-62a385578f37" width="350"/>
```

A body chunk contains all information that is necessary for the simulation:

```@raw html
<img src="https://github.com/user-attachments/assets/c2cd75de-065e-4568-abde-e753216f3e85" width="650"/>
```

## Multithreading

To manage the body chunks that compose a body in multithreading a `ThreadsDataHandler` is employed. It combines information about all the body chunks of the system and all halo exchanges between them.
```@raw html
<img src="https://github.com/user-attachments/assets/48e61d30-cfd5-424c-bc36-c52485a00dc2" width="250"/>
```

## MPI

For MPI simulations an `MPIDataHandler` is used for each body chunk, processing its halo exchanges and further MPI-related information.
```@raw html
<img src="https://github.com/user-attachments/assets/ab389e50-05f7-4d6e-8b0d-992fb6c79a75" width="250"/>
```
