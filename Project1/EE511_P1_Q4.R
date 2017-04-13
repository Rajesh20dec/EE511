Flip_coins<-function(){
  No_of_Heads<-0
  No_of_Tosses<-0
  No_of_Heads<-0
  trial_outcome<-vector()
  No_of_Heads_Required<-readline(prompt="Specified postive no of heads:")
  while(1){
    
    No_of_Tosses<-No_of_Tosses+1
    trial_outcome[No_of_Tosses]<-sample(c(0,1),1,replace=FALSE)
    if(trial_outcome[No_of_Tosses]==1){
        No_of_Heads=No_of_Heads+1
      }
      if (No_of_Heads_Required==No_of_Heads){
            break
          }
          
        }
    
  
  cat("Outcome of   N flips","\n", trial_outcome,sep=" ","\n")
  print(paste("No.of tosses required until reaching ",No_of_Heads_Required," no. of heads is",No_of_Tosses))
}
