function [DendOut] = NormalDend(Input, BU, BL, AU, AL, Distance, branchtosoma, ISI)
%Compute using the same biophysical parameters as Singh&Zald 2015
%   Input has one column per synapse and distance has one entry per synapse
%   with the first entry spacing to the branchpoint and other entries
%   spacing between synapses
%   Except for double-pulse stimulation protocols ISI should be zero

Cm=3.14*10^-14;
Rm=(1/3.14)*10^11;
gNMDA=3.9*10^-9;
NMDAmid=.0463;
NMDAks=2.5*10^-3;
LengthConstant=77*10^-6;
TnmdaOnVec=[4.86 28.9 7472]*10^-3;
AOnVec=[.17 .08 .13]/(.17+.08+.13);
Enmda=.07;

DendOut=DendTranBranch(Input, BU, BL, AU, AL, Distance, Rm, gNMDA,  Cm, Enmda, NMDAks, NMDAmid, LengthConstant,  branchtosoma, TnmdaOnVec,AOnVec, ISI);
end

