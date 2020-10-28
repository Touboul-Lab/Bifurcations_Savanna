% close all

sigma=0.0001;
NEW_ICS=1;
DO_PLOTS=0;

% Basic parameters
mu=0.1;
nu=0.05;
beta=0.4;
alpha=0.57;
s1=0.01;
s2=0.05;
t1=0.4;
t2=0.4;
w0=0.9;
w1=.4;
f0=0.1;
f1=0.9;
MC=20;


% Simulation Parameters:
dt=0.0001;
T=1500;
timespan=0:dt:T;
Ntimes=length(timespan);
if NEW_ICS
    G=zeros(1,Ntimes);
    S=zeros(1,Ntimes);
    TT=zeros(1,Ntimes);

    G(1)=0.1;
    S(1)=0.1;
    TT(1)=0.1;
else
    G(1)=G(end);
    S(1)=S(end);
    TT(1)=TT(end);
end

sss=0;
s=-4:0.2:-0.4;
PerSaved=zeros(1,length(s));
AmplSaved=zeros(1,length(s));
VarPerSaved=zeros(1,length(s));
VarAmplSaved=zeros(1,length(s));

% 
for sigma=10.^s
   sss=sss+1; 
   tic
   Ampls=zeros(1,MC);
   for reps=1:MC
for k=2:Ntimes
    dG=(mu*S(k-1) + nu*TT(k-1) + (f0+(f1-f0)/(1+exp(-(G(k-1)-t2)/s2)))*(1-(G(k-1)+S(k-1)+TT(k-1))) ...
                  -beta*G(k-1)*TT(k-1)-alpha*G(k-1)*(1-(G(k-1)+S(k-1)+TT(k-1))));
    dS=(beta*G(k-1)*TT(k-1)-(w0 + (w1-w0)/(1+exp(-(G(k-1)-t1)/s1)))*S(k-1)-mu*S(k-1) ...
                        -alpha*S(k-1)*(1-(G(k-1)+S(k-1)+TT(k-1))));
	dTT=((w0+(w1-w0)/(1+exp(-(G(k-1)-t1)/s1)))*S(k-1)-nu*TT(k-1)-alpha*TT(k-1)*(1-(G(k-1)+S(k-1)+TT(k-1))));
    noiseG=sqrt(dt)*sigma*randn();
    noiseS=sqrt(dt)*sigma*randn();
    noiseT=sqrt(dt)*sigma*randn();
    
    if ((G(k-1)+dG*dt+noiseG<1) && (G(k-1)+dG*dt+noiseG>0))
        G(k)=G(k-1)+dG*dt+noiseG;
        ng=1;
    else
        G(k)=G(k-1)+dG*dt-(noiseG);
        ng=-1;
    end
    
    if ((S(k-1)+dS*dt+noiseS<1) && (S(k-1)+dS*dt+noiseS>0))
        S(k)=S(k-1)+dS*dt+noiseS;
        ns=1;
    else
        S(k)=S(k-1)+dS*dt-noiseS;
        ns=-1;
    end
    
    if (TT(k-1)+dTT*dt+noiseT<1) && (TT(k-1)+dTT*dt+noiseT>0)
        TT(k)=TT(k-1)+dTT*dt+noiseT;
        ntt=1;
    else
        TT(k)=TT(k-1)+dTT*dt-noiseT;
        ntt=-1;
    end 
    
    if (G(k)+S(k)+TT(k))>1
        G(k)=min(1,max(0,G(k)-2*ng*(noiseG)));
        S(k)=min(1,max(0,S(k)-2*ns*(noiseS)));
        TT(k)=min(1,max(0,TT(k)-2*ntt*(noiseT)));
    end
end
toc
% figure()
% plot(timespan,G);
% hold on
% plot(timespan,S);
% plot(timespan,TT);
% legend('Grass','Sapplings','Trees');
% 
% F=1-(G+S+TT);
% figure()
% plot3(G(end/2:end),S(end/2:end),F(end/2:end))

L=length(timespan);
f = (1:(L/2))/(L*dt);
periods=1./f;
ft=abs(fftshift(fft(G)))*dt;
halfft=ft(end/2+1:end);
if DO_PLOTS
    figure()
    plot(periods,halfft)
end
Ampls(reps)=max(halfft);
Pers(reps)=periods(halfft==Ampls(reps));
   end
AmplSaved(sss)=mean(Ampls);
VarAmplSaved(sss)=var(Ampls);
PerSaved(sss)=mean(Pers);
VarPerSaved(sss)=var(Pers);
end

% figure()
% SINE=sin(2*pi*timespan/100);
% plot(timespan,SINE)
% FF=abs(fftshift(fft(SINE)));
% figure()
% plot(1./f,FF(end/2:end))
% plot(G)
