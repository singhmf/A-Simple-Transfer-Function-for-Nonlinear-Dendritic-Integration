function [ PhiMat] = LengthMat(TotDist,LengthConstant)
%Generates the matrix of  spatial decay components in equation
%5.1

%TotDist is the row vector of distances from each synapse to the soma

n=length(Dist); %Number of input sites

%As in paper, the LengthConstant for contributiing in spike generation is
%halved to reflect the greater spatial dependence of contribution to spikes
%than subthreshold transients.

PhiMat=zeros(n,n);%Pre-allocate for speed
for finish=1:n
    for start=finish+1:n
        PhiMat(finish, start)=exp((TotDist(start)-TotDist(finish))/(-LengthConstant/2));%To ensure esponent is negative, it is crucial that synapses are ordered proximal to distal
    end
end
PhiMat=PhiMat+PhiMat'; %Adds backprop 
end

