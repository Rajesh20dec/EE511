function Proj_P5_Q2()
clc;
clear all;
Rij=0.5;
figure_Count=1;
%this function call is when Rij=0.5 
sprintf('below print is for the case when Rij=0.5')
figure_Count=HOL_Simulation_Algo(figure_Count,Rij);
Rij=0.75;
%this function call is when Rij=0.75 
sprintf('below print is for the case when R1=0.75 and r2=0.25')
figure_Count=HOL_Simulation_Algo(figure_Count,Rij);



end

%------------------------------------------------------------------%
%Funtion name:HOL_Simulation_Algo(figure_Count,Rij)
%Input parameters:figure_Count,Rij
%Output Parameters:figure_Count
%Defination:This function simulates all condition involves In 2x2 HOL switch and
%calculates follwoing parmaeters
%           Mean #of packets at buffer Inport1
%           Mean #of packets at buffer Inport2
%           Mean #of packets processed by switched time perslot
%--------------------------------------------------------------------%




function[figure_Count]= HOL_Simulation_Algo(figure_Count,Rij)
% No of time slots
Time_slot=20;
Prorobality_of_arrival=[0:0.05:0.95];
% buffers for input ports 1 and 2
b1=0; 
b2=0; 
buffer1=[];
buffer2=[];
processed_packets=[];
outport1 =zeros(1,Time_slot);
outport2=zeros(1,Time_slot);
for i=1:Time_slot
Processed_packet_count=0;
Probable_Inport1(i)=rand(1,1);
Probable_Inport2(i)=rand(1,1);
if(Probable_Inport1(i)<=Prorobality_of_arrival(i) && Probable_Inport2(i)<=Prorobality_of_arrival(i)) 
 Inport1(i)=Probable_Inport1(i);
 Inport2(i)=Probable_Inport2(i);

if (Inport1(i)<=Rij)
    Inport1(i)=0;%this is for outport0  
else
    Inport1(i)=1;%this is for outport1
end
if (Inport2(i)<=Rij)
    Inport2(i)=0;%this is for outport0
else
    Inport2(i)=1;%this is for outport1
end

if(i>1)
if(outport1(i-1)== 1)          %Packet on outport was previously reached
  Inport2(i)=Inport2(i-1);
elseif(outport2(i-1)==1)       %Packet on outport was previously reached
  Inport1(i)=Inport1(i-1);
end
end
%This conditions check HOL 
 if(Inport1(i)==Inport2(i))   
  m=rand(1,1);
   if(m<=0.5) %this condition will decides which packet wil be sent 
     b2=b2+1; % Packet on port2 is buffered ( in b2 )
     outport1(i)=1; %Packet will be transferred to outport1
     Processed_packet_count = Processed_packet_count +1;
   else
     b1=b1+1; %Packet on port1 is buffered ( in b1 )
     outport2(i)=1; %Packet will be transferred to outport2
     Processed_packet_count = Processed_packet_count +1;
   end
 else
     Processed_packet_count = Processed_packet_count +2;
 end
 elseif(Probable_Inport1(i)<=Prorobality_of_arrival(i) && Probable_Inport2(i)>Prorobality_of_arrival(i))
     Processed_packet_count=1;
 elseif(Probable_Inport1(i)>Prorobality_of_arrival(i) & Probable_Inport2(i)<=Prorobality_of_arrival(i))
     Processed_packet_count=1;
else
    Processed_packet_count=0;
end
 processed_packets = [processed_packets Processed_packet_count];
 buffer1 = [buffer1 b1];
 buffer2= [buffer2 b2];

end


%Below code will plot Distribution of packets at inputs Vs probability of arrival
figure(figure_Count)

plot(Prorobality_of_arrival,buffer1,'b',Prorobality_of_arrival,buffer2,'r--o');
title('Distribution of packets at inputs Vs probability of arrival')
xlabel('Probality');
ylabel('# of packets at inport1 and inport2');
legend('Inport1','Inport2')
Mean_No_of_packets_atInport1=mean(buffer1);
Mean_No_of_packets_atInport2=mean(buffer2);
sprintf('mean of the number of packets in the buffer at input 1 =%f',Mean_No_of_packets_atInport1)
sprintf('mean of the number of packets in the buffer at input 2 =%f',Mean_No_of_packets_atInport2)

%Below code will plot # number of packets processed per time slot
figure_Count =figure_Count+1;
figure(figure_Count)
histogram(processed_packets)
Mean_No_of_packets_a_processed_by_siwtch_perslot=mean(processed_packets);
sprintf('Mean_No_of_packets_a_processed_by_switch_perslot =%f',Mean_No_of_packets_a_processed_by_siwtch_perslot)
histogram(processed_packets)
title('Distribution of number of # processed packet per time slot');
xlabel('# of packets processed');
ylabel('Frequency');
figure_Count=figure_Count +1;
efficency=processed_packets/2;
Confidence_interval=bootci(Time_slot,@mean,efficency);
sprintf('Confidence interval =[%f %f]',Confidence_interval(1,1),Confidence_interval(2,1))

end


