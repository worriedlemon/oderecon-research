function setup_research
    addpath real_data
    addpath orthpoly-research
    addpath iosindy-research
    addpath fractional-research
    
    olddir = pwd;
    cd oderecon
    setup_oderecon path
    cd(olddir)
end