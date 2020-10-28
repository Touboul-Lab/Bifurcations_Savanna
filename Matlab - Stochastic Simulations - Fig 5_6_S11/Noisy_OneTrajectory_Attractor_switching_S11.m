% close all

sigma=0.1;
DO_PLOTS=1;
NEW_ICS=1;

% Basic parameters
mu=0.1;
nu=0.05;
beta=0.2;
alpha=4;
alpha_het=0.53;
s1=0.01;
s2=0.05;
t1=0.4;
t2=0.4;
w0=0.9;
w1=.4;
f0=0.1;
f1=0.9;
MC=1;


% Simulation Parameters:
dt=0.01;
T=10000;
timespan=0:dt:T;
Ntimes=length(timespan);
if NEW_ICS
    G=zeros(1,Ntimes);
    S=zeros(1,Ntimes);
    TT=zeros(1,Ntimes);

    G(1)=0.1;
    S(1)=0.1;
    TT(1)=0.1;
    
    G_det=zeros(1,Ntimes);
    S_det=zeros(1,Ntimes);
    TT_det=zeros(1,Ntimes);

    G_det(1)=0.1;
    S_det(1)=0.1;
    TT_det(1)=0.1;
    
    G_het=zeros(1,Ntimes);
    S_het=zeros(1,Ntimes);
    TT_het=zeros(1,Ntimes);

    G_het(1)=0.1;
    S_het(1)=0.1;
    TT_het(1)=0.1;
else
    G(1)=G(end);
    S(1)=S(end);
    TT(1)=TT(end);
end

time_below=[];
for reps=1:MC
    progressbar(reps,MC)
    for alpha=.9 %0.4:0.02:1.6
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
    
%     G_det(k)=min(1,max(0,G_det(k-1)+dt*(mu*S_det(k-1) + nu*TT_det(k-1) + (f0+(f1-f0)/(1+exp(-(G_det(k-1)-t2)/s2)))*(1-(G_det(k-1)+S_det(k-1)+TT_det(k-1))) ...
%                   -beta*G_det(k-1)*TT_det(k-1)-alpha*G_det(k-1)*(1-(G_det(k-1)+S_det(k-1)+TT_det(k-1))))));
%     S_det(k)=min(1,max(0,S_det(k-1)+dt*(beta*G_det(k-1)*TT_det(k-1)-(w0 + (w1-w0)/(1+exp(-(G_det(k-1)-t1)/s1)))*S_det(k-1)-mu*S_det(k-1) ...
%                         -alpha*S_det(k-1)*(1-(G_det(k-1)+S_det(k-1)+TT_det(k-1))))));
% 	TT_det(k)=min(1,max(0,TT_det(k-1)+dt*((w0+(w1-w0)/(1+exp(-(G_det(k-1)-t1)/s1)))*S_det(k-1)-nu*TT_det(k-1)-alpha*TT_det(k-1)*(1-(G_det(k-1)+S_det(k-1)+TT_det(k-1))))));
%     if((G_det(k)+S_det(k)+TT_det(k))>1)
%         G_det(k)=G_det(k)-0.01;
%         S_det(k)=S_det(k)-0.01;
%         TT_det(k)=TT_det(k)-0.01;
%     end
%     
%     G_het(k)=min(1,max(0,G_het(k-1)+dt*(mu*S_het(k-1) + nu*TT_het(k-1) + (f0+(f1-f0)/(1+exp(-(G_het(k-1)-t2)/s2)))*(1-(G_het(k-1)+S_het(k-1)+TT_het(k-1))) ...
%                   -beta*G_het(k-1)*TT_het(k-1)-alpha_het*G_het(k-1)*(1-(G_het(k-1)+S_het(k-1)+TT_het(k-1))))));
%     S_het(k)=min(1,max(0,S_het(k-1)+dt*(beta*G_het(k-1)*TT_het(k-1)-(w0 + (w1-w0)/(1+exp(-(G_het(k-1)-t1)/s1)))*S_het(k-1)-mu*S_het(k-1) ...
%                         -alpha_het*S_het(k-1)*(1-(G_het(k-1)+S_het(k-1)+TT_het(k-1))))));
% 	TT_het(k)=min(1,max(0,TT_het(k-1)+dt*((w0+(w1-w0)/(1+exp(-(G_het(k-1)-t1)/s1)))*S_het(k-1)-nu*TT_het(k-1)-alpha_het*TT_het(k-1)*(1-(G_het(k-1)+S_het(k-1)+TT_het(k-1))))));
end
F=1-G-S-TT;
time_below(end+1,:)=[alpha,sum(F-G >0)/sum(F>-1)];
    end
    
end


%
az=45;
el=20;

close all
F=1-(G+S+TT);

figure()
plot(timespan,G);
hold on
plot(timespan,S);
plot(timespan,TT);
plot(timespan,F);
%legend('Grass','Sapplings','Trees','Forest');
ylim([0 1])

% axis([0 500 0 1])

% F=1-(G+S+TT);
% figure()
% plot3(G(ceil(end/2):end),S(ceil(end/2):end),F(ceil(end/2):end))
% 
% 
% figure()
% plot(timespan,G_det);
% hold on
% plot(timespan,S_det);
% plot(timespan,TT_det);
% legend('Grass','Sapplings','Trees');
% axis([0 500 0 1])
% 
% F=1-(G+S+TT);
% F_det=1-(G_det+S_det+TT_det);
% figure()
% plot3(G(ceil(end/2):end),S(ceil(end/2):end),F(ceil(end/2):end))
% hold on
% plot3(G_det(end/2:end),S_det(end/2:end),F_det(end/2:end),'r','LineWidth',1.2)
% view([az el])
% 
% F_het=1-(G_het+S_het+TT_het);
% figure()
% plot3(G(ceil(end/2):end),S(ceil(end/2):end),F(ceil(end/2):end))
% hold on
% plot3(G_det(end/2:end),S_det(end/2:end),F_det(end/2:end),'g','LineWidth',1.2)
% plot3(G_het(end/2:end),S_het(end/2:end),F_het(end/2:end),'r','LineWidth',1.2)
% view([az el])     
%%
% L=length(timespan);
% f = (1:(L/2))/(L*dt);
% periods=1./f;
% ft=abs(fftshift(fft(G)))*dt;
% ft_det=abs(fftshift(fft(G_det)))*dt;
% halfft=ft(end/2+1:end);
% halfft_det=ft_det(end/2+1:end);
% if DO_PLOTS
%     figure()
%     plot(periods,halfft)
%     hold on
%     plot(periods,halfft_det,'r')
% end
% fprintf('Period=%f, Deterministic=%f\n',periods(halfft==max(halfft)),periods(halfft_det==max(halfft_det)))
% 
% 
% L=length(timespan);
% f = (1:(L/2))/(L*dt);
% periods=1./f;
% ft=abs(fftshift(fft(G)))*dt;
% halfft=ft(end/2+1:end);
% if DO_PLOTS
%     figure()
%     plot(periods,halfft)
% end
% Ampls(reps)=max(halfft);
% Pers(reps)=periods(halfft==Ampls(reps));
%    end
% AmplSaved(sss)=mean(Ampls);
% VarAmplSaved(sss)=var(Ampls);
% PerSaved(sss)=mean(Pers);
% VarPerSaved(sss)=var(Pers);
% end

% figure()
% SINE=sin(2*pi*timespan/100);
% plot(timespan,SINE)
% FF=abs(fftshift(fft(SINE)));
% figure()
% plot(1./f,FF(end/2:end))
% plot(G)
