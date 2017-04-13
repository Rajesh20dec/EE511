function P2_Q2()
No_of_samples=input('enter the number of sample required:');
Xk=rand(1,No_of_samples);
shift_value_by=[1 2 3];%use to shift Xk by 1,2 and 3 units.
Xk1 =zeros(1,No_of_samples);
Xk2 =zeros(1,No_of_samples);
Xk3 =zeros(1,No_of_samples);
Xk1(shift_value_by(1)+1:end)=Xk(1:end-shift_value_by(1));%it computesX(k-1)
Xk2(shift_value_by(2)+1:end)=Xk(1:end-shift_value_by(2));%It computesX(k-2)
Xk3(shift_value_by(3)+1:end)=Xk(1:end-shift_value_by(3));%It computes X(k-3)
Yk=[Xk-2*Xk1+0.5*Xk2-Xk3];%It computes Y(k)= X(k)-2.X(k-1)+0.5*X(k-2)-X(k-3)
Covariance_Result=cov(Xk,Xk1);
%disp(Covariance_Result)
%disp(var(Xk))
sprintf('Covariance of Xk and Xk+1 is=%f',Covariance_Result(1,2))
%disp(corrcoef(Xk,Xk1));
Covariance_Result1=cov(Xk,Yk);
sprintf('Covariance of Xk and Yk is=%f',Covariance_Result1(1,2))
end
