alpha=0.4:0.02:1.6;
M=[];
v=[];

for a=alpha
    M(end+1)=mean(time_below(time_below(:,1)==a,2));
    v(end+1)=var(time_below(time_below(:,1)==a,2));
end

errorbar(alpha,M,sqrt(v),'linewidth',1.4)
hold on
errorbar(alpha,1-M,sqrt(v),'linewidth',1.4)
ylim([0,1])
