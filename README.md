text doc to accompany 1d pyrolysis solver by Diba Behnoudfar and 
Michelle Gee 7/13/23

this project contains a few files under the umbrella of the live fuel consumption model. 

the 1d pyrolysis solver iteratively solves the conservation equations wrt destruction/production rate of the reactants (omega dot triple prime) during pyrolysis. (DB)

the TGA(2) simulation has been modified and iteratively solves the continuity equation with detailed chemical kinetics from species pyrolysis (no energy equation). (DB)

the sugar_tga file has been simplified as a initial step in adding the surrogate kinetic scheme to the pyrolysis solver. (MG)

to run these files, files containing stored variables are needed. currently, pyrolysis_ponpo.m requires 'solid_kinetics_data_v1.mat,' while the tga models (.m) require 'ranzi_pyro_kinetics_gentile2017.mat.' 'ranzi_pyro_kinetics_gentile2017.mat' is outdated and may be phased out eventually.
