function Proj5_Q3()
warning off;
clear all;
clc;
N=10;
No_of_experiment=1;
steady_state=[];
Initial_input=zeros(1,2*N+1);
Initial_input(101)=1;
markov_chain(No_of_experiment,Initial_input);


Num_of_1=[1,2,5];
%this function generates initial P matrix in which 1 at 5 diffrent positions are placed out of 201
%elements
for i=1:length(Num_of_1)
    Initial_input=zeros(1,2*N+1);
    for j=1:Num_of_1(i)
    Index_of_1=randi(1:2*N+1); %it will calculate random position where 1 wii be inserted
    Initial_input(Index_of_1)=1/Num_of_1(i);
    end
    markov_chain(1,Initial_input);
end
%------------------------------------------------------------------%
%Funtion name:markov_chain(No_of_experiment,Initial_input)
%Input parameters:No_of_experiment,Initial_input
%Output Parameters:None
%Defination:This function simulates steady state condtion using markov 
%using markov.The input argument Initial_input is intial P matrix reuired
%for Markov chain.No_of_experiment is how many times ,this steady states 
%conditions needs to be simulated using same intial P matrix.
%--------------------------------------------------------------------%


function markov_chain(No_of_experiment,Initial_input)
N=100;
for iteration=1:No_of_experiment

    P=zeros(2*N+1,2*N+1);
    for i=1:2*N+1
        for j = 1:2*N+1
         P(i,j) = nchoosek(2*N,j-1)*((i-1)/(2*N))^(j-1)*(1-(i-1)/(2*N))^(2*N-j+1);
     end
    end
n=5000;           % number of time steps to take
output=zeros(n+1,2*N+1); % clear out any old values
t=0:n;% time indices
output(1,:)=Initial_input;
i = 0;
    for i=1:n,
    output(i+1,:) = output(i,:)*P;
        %a tolerance check to  automatically stop the simulation when the density is close to its steady-state
     LIT = ismembertol(output(i+1,:),output(i,:));
        if all(LIT == 1)     
             break;
        end
    end
    sprintf('time steps required to reach steady state =%f',i)
    steady_state(iteration)=i;
end
%histogram(steady_state);
end