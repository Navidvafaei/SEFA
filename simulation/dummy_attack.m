key_t=100;%number of different keys that might attacker consider
sample_t=500;%number of random plaintext 
fb=2;%number of random bit-faults
k=2;% it shows one dummy round
Sel_R=10*k;% fault is injected to the last dummy round

%%%%%%%%%%%%%%%following function returns faulty and non faulty ciphertext
%%%%%%%%%%%%%%%by considereing number of dummy round which can be specify
%%%%%%%%%%%%%%%by k, k=2 have one dummy round. Sel_R specify the round
%%%%%%%%%%%%%%%which an attacker want to inject faults
[key_col,cipherc,cipherf]=dummyround(sel_R,k,sample_t,key_t,fb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%following function returns calculation of LLR and SEI based
%%%%%%%%%%%%%%%on the input and rank the corrected key after saample_t plaintext. 
[rank_eff_sei,rank_ineff_sei,rank_joint_sei,rank_eff_llr,rank_ineff_llr,rank_joint_llr]=sifa_sefa_calc(key_col,cipherc,cipherf,sample_t,key_t,fb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%following function returns the plot of average for SEI and
%%%%%%%%%%%%%%%LLR

plot(mean(rank_eff_sei(1:key_t-1,:),1),'b')
hold on
plot(mean(rank_ineff_sei(1:key_t-1,:),1),'r')
hold on
plot(mean(rank_joint_sei(1:key_t-1,:),1),'g')
