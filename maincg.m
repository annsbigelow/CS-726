clear 
global numf numg

n=8; m=n/4;
x0=[];
for i=1:m
    x0(4*i+1)=3; x0(4*i+1)=-2;
    x0(4*i+3)=0; x0(4*i+4)=1;
end

numf=0; numg=0;

[inform,x] = SteepDescent(@xpowsing,x0,sdparams);

