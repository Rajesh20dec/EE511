function Acc_Rej_Method()
clc
clear all
P=[0.06 0.06 0.06 0.06 0.06 0.15 0.13 0.14 0.15 0.13 0 0 0 0 0 0 0 0 0 0];
Q(1:20)=0.05;
c=max([P./Q]);
sprintf('c=%f',c)
C=[];
No_of_Samples=10000;
X=[];

for i=1:No_of_Samples 
    flag=-1;
    k=0;
    while flag~=0
        k=k+1;
        U=rand();
        Y1=floor(20*rand)+1;
        if U<(P(Y1)/(c*Q(Y1)))
            flag=0;
            X(i)=Y1;
            C(i)=k;
        end
    end
 
end

Sample_mean=mean(X);
Sample_Variance=var(X);
therotical_mean=sum(P.*([1:20]));
temp=([1:20]-therotical_mean).^2;
therotical_variance=sum(P.*temp);
estimated_efficency=10000/sum(C(1:end));
Theoretical_efficency=1/c;
sprintf('Sample mean =%f and therotical mean=%f',Sample_mean,therotical_mean)
sprintf('Sample vairance =%f and therotical vairance=%f',Sample_Variance,therotical_variance)
sprintf('estimated_efficency =%f and Theoretical_efficency=%f',estimated_efficency,Theoretical_efficency)
disp(mean(C));

No_bins = 1:20;
No_counts = histc(X,No_bins);
handle1 = bar(1:20,No_counts/sum(No_counts))
title('Plot of target distribution(Pj) and random samples generated using Accept-reject method')
hold on;
handle2 = plot(P,'r')
legend([handle2,handle1],{'Target distribution(Pj)','samples generated using accept-reject method'});
xlabel('X=xj'); ylabel('P(X=xj)');
xlim([0 20 + 1]);
end
