function [DendOut] = DendTranBranch(Input, BU, BL, AU, AL, Distance, Rm, gNMDA,  Cm, Enmda, NMDAks, NMDAmid, LengthConstant,  branchtosoma, TnmdaOnVec,AOnVec, ISI)
%%
%Generates an approximate transformation of a dendritic input combination to
%the peak induced somatic EPSP using the simple transfer function in Singh&Zald, 2015
%Please note that ALL POTENTIALS ARE TRANSLATED FOR A RESTING POTENTIAL OF ZERO
%Synapses are always ordered proximal to distal.
%For double-pulse stimulation, set ISI equal the inter-stimulus interval otherwise set ISI=0

%Input is a txn vector of inputs, however the function does not include
%temporal interactions between stimuli within the dendrites

% Bu,Bl,Au,Al are the constants in Singh&Zald,2015 giving boundaries
% (Bu=upper, Bl=lower)and curvatures (Au, Al) of upper and lower limits in the linear boundary function
% (Equation 1)
%Distance is a row vector of separations with Distance(1) equal the
%separation from the dendritic branch point and Distance(i+1) equal the
%distance from synapse (i+1) to synapse(i)

%Original simulations in Sigh&Zald, 2015 were run in base form (e.g. Resistance in Ohms,
%Conductance in Seimens, Potential in Volts, length in Meters, Capacitance
%in Farads).
%% Biophysical Parameters: 
%Rm [Membrane Resistance], 
%gNMDA[NMDA max conductance], 
%Cm [Capacitance], 
%Enmda[Translated NMDAR reversal potential]

%NMDAks gives the denominator of the exponent in Mg2+ blockade (Ks in
%equation 7), while NMDAmid gives the midpoint (half opening potential).

%BranchtoSoma gives the distance from the start of the dendritic branch to
%the soma

%TNMDAOnVec is the row vector of exponential coefficients for the NMDAR
%inter-opening distribution
%AOnVec is the corresponding row vector of coefficient weights (between 0 and 1)


n=length(Distance); %Number of input sites
Nin=size(Input,1); %Number of input cases

TotalDist=zeros(1,n);
for i=1:n
    TotalDist(i)=sum(Distance(1:i))+branchtosoma; %Converts separtions into distance from the soma
end

Tm=Rm*Cm; %Time Constant


PhiArray=LengthMat(TotalDist, LengthConstant); %Calculates the spatial decay matrix of equation 5.1

nonlocal=(PhiArray*Input')'; %Gives spatially decayed inputs to other sites.

GluTrue=Input>0; %Boolean matrix used to restrict spike generation to glutamate-activated synapses (Non-zero input)

Spikes=spike(Input,nonlocal,GluTrue, Tm, DecayConst, Rm, gNMDA, Enmda, NMDAks, NMDAmid, TnmdaOnVec, AOnVec); %Generate spiking (nonlinear) components

if ISI>0
    Linear=Input*(1+exp(-ISI/Tm)); %Double Pulse Contribution to fast currents
else
    Linear=Input; %Single Pulse
end

SomaDecay=repmat(exp(-TotalDist/LengthConstant),Nin,1); %Spatial decay in reaching the soma

Prebounds=SomaDecay.*(Linear+Spikes); %Before application of boundary function for dendritic saturation

DendOut=SigZ(Prebounds*ones(n,1), BU, BL, AU, AL); %Final approximation of peak induced depolarization at the soma
end