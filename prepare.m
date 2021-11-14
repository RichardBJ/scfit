addpath("models","Qmatrix","etc","exp_mix_dist","SCAN","scripts","src");
%if compiled those functions... then also
addpath(genpath("codegen"),genpath("build"));
%otherwise drops a harmless error, but compiling/mexifying really recommended!!