%------------------------------------------------------------------%
%Funtion name:Flip_coins()
%Input parameters:None
%Output Parameters:None
%Defination:This function simulates until nos. of heads outcome matches user
%           input's nos. of heads         
%--------------------------------------------------------------------

function trial_outcomes =Flip_coins()
    No_of_Heads_Required=input('Specified postive no of heads:');
    No_of_Tosses=0;
    outcome_of_trials= [];
    trials= [];
    No_of_Heads=0;
    while(1)  
        
        temp=rand();
        No_of_Tosses=No_of_Tosses+1;
        
        trials(No_of_Tosses)=temp;
        if trials(No_of_Tosses)>0.5000
            outcome_of_trials(No_of_Tosses)=1;
        else
          outcome_of_trials(No_of_Tosses)=0;
        end
         if outcome_of_trials(No_of_Tosses)==1 %condition to get heads
            No_of_Heads=No_of_Heads+1;
         end
        
        if No_of_Heads_Required == No_of_Heads
            break;
        end
       
    end
        
    sprintf('Benoulli outcome(n trials)=')
    disp(outcome_of_trials)
    sprintf('No.of tosses required until reaching %d no. of heads is=%d',No_of_Heads_Required,No_of_Tosses)
   
end




    
    