function [Spike] = Spike(Input,nonlocal,GluTrue, Tm, Rm, gNMDA, Enmda, NMDAks, NMDAmid, TnmdaOnVec, AOnVec, ISI)
%Spike Generation Function, with input variables as defined in DendTranBranch
TDecay=AOnVec*Tm*((TnmdaOnVec+Tm).^-1); %Expected local temporal decay from Equation 4

Depol1=Input*TDecay+nonlocal; %Total depolarization at channel opening

Constants1=gNMDA*Enmda/(gNMDA+Rm^-1); %Numerator of equation 13.2
Exp1=(Depol1-NMDAmid)/NMDAks+log(1+gmaxNMDA*Rm); %Exponent of equation 13.2
Spike1=Constants1*GluTrue./(1+exp(-1*Exp1)); %Spikes of first pulse


if ISI==0;
    Spike2=0; %single pulse case
elseif ISI>0
    Depol2=Depol1*(1+exp(-1*ISI/Tm)); %In paper made Depolarization 2 even greater by decreasing TDecay for the first pulse's contribution to the second spike
    g2=gmaxNMDA*max(1,1.6*(1-ISI/.1)); % Using method of Gomez-Gonzalez to increase NMDAR conductance for second pulse with short ISI (ISI parameter in seconds)
    Constants2=g2*Enmda/(g2+Rm^-1);
    Exp2=(Depol2-NMDAmid)/NMDAks+log(1+g2*Rm);
    Spike2=Constants2*GluTrue./(1+exp(-1*Exp2));
end
Spike=min(Enmda,Spike1+Spike2); %Prevents paired-pulse summation from passing reversal potential
end

