function Proj3_Q1()
clc

Prob_of_Rejection=lot_sampling(5);
theoretical_prob_rejection = hygepdf(1,125,6,5);
sprintf('Prob_of_Rejection(with 5 samples at a time) =%f and theoretical_prob_rejection(with 5 samples at a time)=%f',Prob_of_Rejection,theoretical_prob_rejection)
Ninety_five_Perc_rejection();
end

function [Prob_of_Rejection]=  lot_sampling(sample_size)
outcome=[];
Total_experiments=1000;
No_of_Reject=0;

    for i=1:Total_experiments
        outcome=datasample(1:125,sample_size,'Replace',false);
        if min(outcome)<=6 %it compares if in lot atleast one value is <= 6 ,then its rejcted lot
            No_of_Reject=No_of_Reject+1;
        end
    end
Prob_of_Rejection=No_of_Reject/Total_experiments;
end



function Ninety_five_Perc_rejection()
sample_size=1;
Expected_Prob_of_rejection=0.95; % this is 95% chance to reject lot
Calculated_Prob_of_rejection=0;
flag=-1;
    while (Expected_Prob_of_rejection ~= Calculated_Prob_of_rejection) & (sample_size<125)
            Calculated_Prob_of_rejection=lot_sampling(sample_size);
            sample_size=sample_size+1;
        if Calculated_Prob_of_rejection<=(0.95+0001) & Calculated_Prob_of_rejection>=(0.95-0.0001)%since floating point camparison is not allowed i have taken tolerance of+-(0.001)with0.95
            sprintf('the fewest number of microchips to reject this lot 95 Percentage of the time=%f',(sample_size-1))
            break;
        end
    end
end


