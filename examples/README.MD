These are all going to be scenes in an upcoming 5 minute video for SPScon.

Need to concisely explain all research done in 5 minutes for a mid level undergraduate audience... will be a challenge.

Roadmap:
1: Start from bulk gold system, we have a chunk of gold and we want to know how it responds to mechanical stresses, maybe show it stretching until it breaks, with the s/s curve?

2: Zoom in on the gold to show the crystal lattice
     We can simulate this from the perspective (maybe change this wording?) of the atoms inside the material. 

     Show the lattice stretching with the same s/s curve and then breaking.

3: molecular dynamics explanation, we now do the scene of the example atoms in a representative box showing just the pair potential forces.

    Molecular dynamics treats the atoms as particles, and solves their equations of motion to get insights into the larger system. Use the temperature example here, average kinetic energy

4: Interatomic Potentials explanation, scene moves from the KE of the particles, (maybe add in the math animations here?) to the potential energy between them.

    We have the basic example of a pair potential like the lennard jones potential (show graph + particle animation)

    Then we have the embedded atom method, this is the industry standard (show animation) it adds a correction term to account for the work done to put an electron density into a host of other electron densities.

    The last method fails to capture certain effects, one of the large contributors to this is that the electron density isn't necessarily symmetrical. The ML potential can improve on this by including an angular dependence (show electron cloud animation)

5: Using these potentials and concepts we can simulate the larger system. (if I have time ergodicity here, maybe even leave a link to other vid?)

    Show some of the animations from ovito here

    We use a radial distribution function to quantify if we're looking at an ordered crystal structure or a liquid. (show the animations) 
    
    An ordered FCC lattice like gold will lead to a graph with defined peaks such as this, a liquid will have a more varied assortment of distances creating this continuous looking graph.

6: Conclusions / links?



** This probably won't fit into the 5 minute time frame :(

using the velocity example can lead into the phase space of a single particle... ergodicity conga line
 
    In modeling a system like this we need this concept of ergodicty, which would mean that in an ideal case a single atom will explore all points of this phase space given infinite time. 

    What if we discretize all possible velocities of the particle, we get these 64 possible microstates. There are two loops here, something like a potential energy barrier  
    
    In our case we can get close to the larger system with a smaller representative system in a reasonable amount of time.