text doc to accompany 1d pyrolysis solver by Diba Behnoudfar and 
Michelle Gee 7/13/23

this project contains a few files under the umbrella of the live fuel consumption model. 

the 1d pyrolysis solver iteratively solves the conservation equations wrt destruction/production rate of the reactants (omega dot triple prime) during pyrolysis. (DB)

the TGA simulation has been modified for hard and softwood and iteratively solves the continuity equation with detailed chemical kinetics from species pyrolysis (no energy equation). (DB)

the tga_sugars file has been edited for a acetic acid TGA simulation as a initial step in adding the  kinetic scheme to the pyrolysis solver. (MG)

to run this file, a file containing stored variables ('ranzi_pyro_kinetics_gentile2017.mat') is needed.
