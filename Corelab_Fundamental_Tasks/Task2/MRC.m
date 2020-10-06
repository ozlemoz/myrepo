function [BERsim_MRC1, BERsim_MRC2] = MRC()

NumBit = 10^6;
index = 1;

data = rand(1,NumBit)>0.5; % generating 0,1 with equal probability

Nr = [2 4];
for jj = 1:length(Nr)
    index = 1;
    for k = 0:1:45               
        SNR = 10^(k/10);           
        
        h = (sqrt(0.5))*(randn(Nr(jj),NumBit) + 1i*randn(Nr(jj),NumBit));
        n = (sqrt(0.5))*(randn(Nr(jj),NumBit) + 1i*randn(Nr(jj),NumBit));
        
        y = h.*(2*data-1) + (sqrt(1/SNR))*n;
        
        for kk = 1:NumBit
            temp = 0;
            temp1 = 0;
            for i = 1:Nr(jj)
                temp = temp + conj(h(i,kk))*y(i,kk);
                temp1 = temp1 + (norm((h(i,kk)),'fro').^2);  
            end
            d(kk) = temp./temp1;
            data_detect(kk) = real(d(kk))>0;
        end
        
        error = xor(data,data_detect);
        BERsim_MRC(jj,index) = sum(error)/NumBit;      
        snr(index) = k;
        
        index = index+1;
    end
end
BERsim_MRC1=BERsim_MRC(1,:);
BERsim_MRC2=BERsim_MRC(2,:);
end

