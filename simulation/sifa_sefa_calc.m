function [rank_eff_sei,rank_ineff_sei,rank_joint_sei,rank_eff_llr,rank_ineff_llr,rank_joint_llr]=sifa_sefa_calc(key_col,cipherc,cipherf,sample_t,key_t,fb)
mask=(2^fb)-1;
array_eff_ineff_n=zeros(key_t,2);
if fb==2
    q_i=[4/9,2/9,2/9,1/9];
    q_e=[0,2/7,2/7,3/7];
elseif fb==4
    q_i=[0.20,0.10,0.10,0.05,0.10,0.05,0.05,0.02,0.10,0.05,0.05,0.02,0.05,0.02,0.02,0.01];
    q_e=[0,0.045703839,0.045703839,0.068555759,0.045703839,0.068555759,0.068555759,0.079981718,0.045703839,0.068555759,...
    0.068555759,0.079981718,0.068555759,0.079981718,0.079981718,0.085694698];
end
% key_g=zeros(16);
key_base=zeros(16);
key_n=1;
data_d_coll_e=zeros(key_t,256,nt);
data_d_coll_i=zeros(key_t,256,nt);
delta_d_e=zeros(key_t,256,256);
delta_d_i=zeros(key_t,256,256);
SEI_i=zeros(key_t,nt+1,256);
SEI_e=zeros(key_t,nt+1,256);
SEI_joint=zeros(key_t,nt+1,256);

LLR_E=zeros(key_t,nt+1,256);
LLR_I=zeros(key_t,nt+1,256);
LLR_joint=zeros(key_t,nt+1,256);

rank_eff_sei=zeros(key_t,nt);
rank_ineff_sei=zeros(key_t,nt);
rank_joint_sei=zeros(key_t,nt);


rank_eff_llr=zeros(key_t,nt);
rank_ineff_llr=zeros(key_t,nt);
rank_joint_llr=zeros(key_t,nt);

while key_n<key_t
   key_hex=(cell2mat([key_col{key_n}]));
   [s_box, inv_s_box, w, poly_mat, inv_poly_mat] = aes_init(key_hex);
   key_b_0=w(41,1);
        ni=0;
        ne=0;
        for i=1:sample_t
            if isequal(cipherc(:,i,key_n),cipherf(:,i,key_n))
                ni=ni+1;
            else
                ne=ne+1;
            end
            for key_g=0:255
                if isequal(cipherc(:,i,key_n),cipherf(:,i,key_n))
                    data_d_coll_i(key_n,key_g+1,ineff)=bitand(sub_bytes(bitxor((cipherc(1,i,key_n)),key_g),inv_s_box),mask);
                    delta_d_i(key_n,data_d_coll_i(key_n,key_g+1,ineff)+1,key_g+1)=delta_d_i(key_n,data_d_coll_i(key_n,key_g+1,ineff)+1,key_g+1)+1;
                else
                    data_d_coll_e(key_n,key_g+1,eff)=bitand(sub_bytes(bitxor((cipherc(1,i,key_n)),key_g),inv_s_box),mask);
                    delta_d_e(key_n,data_d_coll_e(key_n,key_g+1,eff)+1,key_g+1)=delta_d_e(key_n,data_d_coll_e(key_n,key_g+1,eff)+1,key_g+1)+1;
                end 
            end 
            for m=1:(mask+1)
                if ni>0
                    SEI_i(key_n,i,:)=(((delta_d_i(key_n,m,:)/ni)-(1/(mask+1))).^2)+SEI_i(key_n,i,:);
                end
                if ne>0
                    SEI_e(key_n,i,:)=(((delta_d_e(key_n,m,:)/ne)-(1/(mask+1))).^2)+SEI_e(key_n,i,:);
                    if ni>0
                        for ms=1:(mask+1)
                            SEI_joint(key_n,i,:)=((((delta_d_i(key_n,ms,:)).*delta_d_e(key_n,m,:)/(ne*ni))-(1/((mask+1)^2))).^2)+SEI_joint(key_n,i,:);
                        end
                    end
                end
                for key_g=1:256
                     if (delta_d_e(key_n,m,key_g))>0 && ne>0
                        LLR_E(key_n,i,key_g)=ne*q_e(m)*log2((delta_d_e(key_n,m,key_g)/ne)/(1/(mask+1)))+ LLR_E(key_n,i,key_g);
                     end
                     if (delta_d_i(key_n,m,key_g))>0 && ni>0
                        LLR_I(key_n,i,key_g)=ni*q_i(m)*log2((delta_d_i(key_n,m,key_g)/ni)/(1/(mask+1)))+LLR_I(key_n,i,key_g);
                        if (delta_d_e(key_n,m,key_g))>0 && ne>0
                            for ms=1:(mask+1)
                                LLR_joint(key_n,i,key_g)=ni*q_i(ms)*ne*q_e(m)*log2((((delta_d_i(key_n,ms,key_g)).*delta_d_e(key_n,m,key_g)/(ne*ni))/(1/((mask+1)^2))))+LLR_joint(key_n,i,key_g);
                            end
                        end
                     end
                end
            end
            rank_eff_sei(key_n,i)=sum(SEI_e(key_n,i,:)>=SEI_e(key_n,i,key_b_0+1));
            rank_ineff_sei(key_n,i)=sum(SEI_i(key_n,i,:)>=SEI_i(key_n,i,key_b_0+1));
            rank_joint_sei(key_n,i)=sum(SEI_joint(key_n,i,:)>=SEI_joint(key_n,i,key_b_0+1));
            rank_eff_llr(key_n,i)=sum(LLR_E(key_n,i,:)>=LLR_E(key_n,i,key_b_0+1));
            rank_ineff_llr(key_n,i)=sum(LLR_I(key_n,i,:)>=LLR_I(key_n,i,key_b_0+1)); 
            rank_joint_llr(key_n,i)=sum(LLR_joint(key_n,i,:)>=LLR_joint(key_n,i,key_b_0+1));
        end
    array_eff_ineff_n(key_n,1)=ne;
    array_eff_ineff_n(key_n,2)=ni;
    key_n=key_n+1;
end
end
