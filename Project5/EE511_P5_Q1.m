function Single_Server_Queque()
clc;
clear all;
mean=1/25;
t=0;
NA=0;
ND=0;
n=0;
i=1;
T=100;
Ts=0;
A=0;
D=0;

breakduration=0;
totalbreak=0;
[ta,Ts]=generate_Ts(T,Ts);
%Ts(i)=temp;
breaktime=ta;
td=inf;
while ta<=T %
        if ta<=td &&ta<=T
            t=ta;
            NA=NA+1;
            n=n+1;
            s(1,NA+ND)=t;
            s(2,NA+ND)=n;
            [ta,Ts]=generate_Ts(T,Ts);
            %Ts(i)=temp;
            if n==1
                td=breaktime+exprnd(mean);%generates service time
            end
            A(NA)=t;
        elseif ta>td && td<=T
            t=td;
            n=n-1;
            ND=ND+1;
            s(1,NA+ND)=t;
            s(2,NA+ND)=n;
            D(ND)=t;
            if n==0
                td=inf;
                breaktime=t;
                while breaktime<ta
                    breakduration=0.3*rand;%generates Break duration 
                    breaktime=breaktime+breakduration;
                    totalbreak=totalbreak+breakduration;
                    
                end
                
            else
             td=t+exprnd(mean);
            end
 
            end
            
            
        end
    

while min(ta,td)>T&& n>0
    t=td;
    n=n-1;
    ND=ND+1;
    td=t+exprnd(mean);
    D(ND)=t;
     s(1,NA+ND)=t;
     s(2,NA+ND)=n;

end
if n==0 && min(ta,td)>T
    Tp=max(T-t,0);
end
figure(1)

stairs([0 s(1,:)],[0 s(2,:)],'linewidth',1.5);
xlim([0 T]);
title('Single server Queque model') 
xlabel('t--->')
ylabel('Queque at server')

sprintf('# of arrival %f=',NA)
sprintf('# of daparture %f=',ND)
sprintf('Arrival time of the customer')
disp(A)
sprintf('Daparture time of the customer')
disp(D)
sprintf('Ts=')
disp(Ts)
sprintf('total break time=%f',totalbreak)


end
%------------------------------------------------------------------%
%Funtion name:generate_Ts(T,Ts)
%Input parameters:T,Ts
%Output Parameters:None
%Defination:This function generated the arrival time using inverse
%transform method .T is time interval of obervation,Ts is previous 
%arrival time 
%--------------------------------------------------------------------%



function[t,Ts]= generate_Ts(T,Ts)

t=Ts;
lamda=19;
t1=0;
%T=100;
while(t<T)
       u1 = rand ();
       t = t- log(u1)/lamda;
       u2 = rand();
       if mod(t,10)<=5
           t1=mod(t,10);
           if u2 <= (mod(4+3*t1,10)/lamda)
           Ts = t;
           break;
           end
         else 
           t1=mod(t,10);
          if u2 <= (mod(19-3*(t1-5),10)/lamda)
          Ts = t;
          %i=i+1;
          break;
          end
      end    
   

end
end
