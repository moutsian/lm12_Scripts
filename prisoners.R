prisoners= function(n_people,seq_in){
state='00'
counts_from_key_person=0;
appearances=matrix(nrow=n_people,ncol=2,0) #second column is the state
state_sequences=matrix(nrow=length(seq_in),ncol=1)
colnames(appearances)=c("True counts","been_or_not")
#for the second row, 0 is not been and 1 is been. The first column counts the actual (true) times each person has gone in
for(i in 1:length(seq_in)){
state_sequences[i]=state
appearances[seq_in[i],1]=appearances[seq_in[i],1]+1;
	if(seq_in[i]==1){ #key person
		if(state=='00'){
			state='01';
		}else if(state=="01"){
			state='00';
		}else if(state=="10"){
			state='00'
		}else{
			state='01'
			counts_from_key_person=counts_from_key_person+1;
			print(paste("counts: ",counts_from_key_person,sep=""))
			print(appearances);
			#print(cbind(seq_in[1:i],state_sequences[1:i]));
		}
	}
	else{
		if(appearances[seq_in[i],2]<1){ #not been in yet, plus special cases
			if(state=='00'){
				state='10';
			}else if(state=="01"){
				state='11';
				appearances[seq_in[i],2]=appearances[seq_in[i],2]+1 
			}else if(state=="10"){
				state='11';
				appearances[seq_in[i],2]=appearances[seq_in[i],2]+1
			}else{
				state='10';
				appearances[seq_in[i],2]=appearances[seq_in[i],2]-1
			}			
		}
		else{ #already in
			if(state=='00'){
				state='01';
			}else if(state=="01"){
				state='00';
			}else if(state=="10"){
				state='00'
			}else{
				state='10'
				appearances[seq_in[i],2]=appearances[seq_in[i],2]-1
			}	
		}
	}
	if(counts_from_key_person==n_people-1){
	print(paste("Ready to go out, after a total of ",appearances[1,1]," visits of the room from key person, and a total of ",sum(appearances[,1])," in the room with the switches in total",sep=""))
	if(sum(appearances[,1]>0)==dim(appearances)[1]){
	print("correct");
	}else{
	print("wrong!");
	print(seq_in[1:i])}
	break;
	}
}
print(paste("counts from key person: ",counts_from_key_person,sep=""))
return (appearances)
}

###now run
N=23
randomseq=sample(1:N,30000,replace=T)
tmp=prisoners(max(randomseq),randomseq)

