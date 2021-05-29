# What do I need to write in order to have a functional molecular dynamics simulation?
#TMP Chem's github repo may be a good resource for this, as well as my mol sim textbook

#keep it simple for now

#input xyz file (atomic coordinates)
#calculate bonds
#input timestep
#input simulation length
#input forcefield? probably use a default forcefield
#calculate atomic radii and check for position clashes
#equations for full energy hamiltonian (bonded + nonbond + electrostatic)
#equations for force (gradient of energy)
#calculate initial energy, forces using velocity verlet
#advance simulation by one timestep
#repeat calculations until 