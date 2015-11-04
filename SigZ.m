function [SigZOut] = SigZ( Input, Top, Bottom, SlopeU, SlopeL )
%"Linear Bounds" function of equation 1
SigBottom=log(1+exp(SlopeL*(Input-Bottom)))/SlopeL;
SigTop=log(1+exp(SlopeU*(Input-Top)))/SlopeU;
SigZOut=SigBottom-SigTop+Bottom;
end

