% close all

sigma=0.05;
DO_PLOTS=1;
NEW_ICS=1;

% Basic parameters
mu=0.1;
nu=0.05;
beta=1.5;
alpha=0.34;
s1=0.01;
s2=0.05;
t1=0.4;
t2=0.4;
w0=0.9;
w1=.4;
f0=0.1;
f1=0.9;
MC=1;

close all

color_custom=[133,215,144; 0,114,189;253,140,0;0,00,0]/255;

% Simulation Parameters:
dt=0.00005;
T=1000;
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
% SIGMA_VECT=linspace(log(0.0001),log(0.01),10);
SIGMA_VECT=[0.0001, 0.025, 0.001, 0];
Spikes=zeros(MC,length(SIGMA_VECT));
kalpha=0;
for sigma=(SIGMA_VECT)
    kalpha=kalpha+1;
    
    for mmc=1:MC
        tic
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
% Spikes(mmc,kalpha)=size(findpeaks(G,'MinPeakHeight',0.33,'MinPeakDistance',1000),2);

    end
    figure()
    plot(timespan, G)
    title(sprintf('sigma=%f',sigma))
    ylim([0,1])
    
    figure(10)
    hold on
%     if kalpha==2
%         plot3(G(ceil(end/2):40:end),S(ceil(end/2):40:end),TT(ceil(end/2):40:end))
%     else
%         plot3(G(ceil(end/2):40:end),S(ceil(end/2):40:end),TT(ceil(end/2):40:end),'LineWidth',2)
%     end
F=1-G-S-TT;
    if kalpha<=2
        plot3(G(ceil(end/2):40:end),S(ceil(end/2):40:end),F(ceil(end/2):40:end),'Color',color_custom(kalpha,:))
    elseif k==3
        plot3(G(1:40:end),S(1:40:end),F(1:40:end),'Color', color_custom(kalpha,:),'LineWidth',2)
    else
        plot3(G(ceil(end/2):40:end),S(ceil(end/2):40:end),F(ceil(end/2):40:end),'Color', color_custom(kalpha,:),'LineWidth',2)
    end
end
figure(10)
az=45;
el=20;
view([az,el])

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
