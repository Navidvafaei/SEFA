function [exactmissrate,key_col,cipherc,cipherf]=missrateF(missrate,sample_t,key_t,fb)
key_col=cell(1,1000);
cipherc=zeros(16,15000,100);
cipherf=zeros(16,15000,100);
Pe=100*missrate;
random_f= round(100*rand(key_t,sample_t));
random_f(:)=(random_f(:)<Pe);
exactmissrate=sum(sum(random_f(1:key_t,:)))/sample_t*key_t;
    for key_n=1:key_t
        key(1,1)=round(rand*255)  ;
        key(2,1)=round(rand*255)  ;
        key(3,1)=round(rand*255)  ;
        key(4,1)=round(rand*255)  ;
        key(5,1)=round(rand*255)  ;
        key(6,1)=round(rand*255)  ;
        key(7,1)=round(rand*255)  ;
        key(8,1)=round(rand*255)  ;
        key(9,1)=round(rand*255)  ;
        key(10,1)=round(rand*255) ;
        key(11,1)=round(rand*255) ;
        key(12,1)=round(rand*255) ;
        key(13,1)=round(rand*255) ;
        key(14,1)=round(rand*255) ;
        key(15,1)=round(rand*255) ;
        key(16,1)=round(rand*255) ;
        key_col{key_n} = {dec2hex(key)};
        SS=reshape(key,1,16);
        key_hex=dec2hex(SS);
        [s_box, inv_s_box, w, poly_mat, inv_poly_mat] = aes_init(key_hex);
            for i=1:sample_t
                faultee=1;
                if random_f(key_n,i)==1
                    faultee=0;
                end
                plaintext(1,1)=round(rand*255);
                plaintext(2,1)=round(rand*255);
                plaintext(3,1)=round(rand*255);
                plaintext(4,1)=round(rand*255);
                plaintext(5,1)=round(rand*255);
                plaintext(6,1)=round(rand*255);
                plaintext(7,1)=round(rand*255);
                plaintext(8,1)=round(rand*255);
                plaintext(9,1)=round(rand*255);
                plaintext(10,1)=round(rand*255);
                plaintext(11,1)=round(rand*255);
                plaintext(12,1)=round(rand*255);
                plaintext(13,1)=round(rand*255);
                plaintext(14,1)=round(rand*255);
                plaintext(15,1)=round(rand*255);
                plaintext(16,1)=round(rand*255);
                faultvalue=bitxor((256-(2^fb)),round(rand*((2^fb)-1)));
                [ciphertext,] = cipher (plaintext, w, s_box, poly_mat,0,0,faultvalue);
                cipherc(:,i,key_n)=ciphertext;
                [ciphertextf,] = cipher (plaintext, w, s_box, poly_mat,0,faultee,faultvalue);
                cipherf(:,i,key_n)=ciphertextf;
            end
    end
end







