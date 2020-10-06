function [BERsim_mqam, BERthe_mqam] = awgn_ber_mqam(NumBit,M,SNRdB)

M_qam = 2^M;
bit_size = log2(sqrt(M_qam));
const_size = 2^bit_size; % constellation size  
scaling_factor = (2/3)*(M_qam-1);

% defining the real and imaginary PAM constellation
% for 16-QAM
alphaRe = [-(2*sqrt(M_qam)/2-1):2:-1 1:2:2*sqrt(M_qam)/2-1];
alphaIm = [-(2*sqrt(M_qam)/2-1):2:-1 1:2:2*sqrt(M_qam)/2-1];
k_QAM = 1/sqrt(scaling_factor);

Es_N0_dB  = SNRdB + 10*log10(log2(M_qam));

% Mapping for binary <--> Gray code conversion
ref = 0:const_size-1;
map = bitxor(ref,floor(ref/2));
[~, ind] = sort(map);                                

for k = 1:length(SNRdB)
    
    SNR_mqam = 10.^(k/10);
    
    % symbol generation
    % ------------------
    ipBit = rand(1,NumBit*log2(M_qam),1)>0.5; % random 1's and 0's
    ipBitReshape = reshape(ipBit,log2(M_qam),NumBit).';
    bin2DecMatrix = ones(NumBit,1)*(2.^((log2(M_qam)/2-1):-1:0)) ; % conversion from binary to decimal
    % real
    ipBitRe =  ipBitReshape(:,(1:log2(M_qam)/2));
    ipDecRe = sum(ipBitRe.*bin2DecMatrix,2);
    ipGrayDecRe = bitxor(ipDecRe,floor(ipDecRe/2));
    % imaginary
    ipBitIm =  ipBitReshape(:,(log2(M_qam)/2+1:log2(M_qam)));
    ipDecIm = sum(ipBitIm.*bin2DecMatrix,2);
    ipGrayDecIm = bitxor(ipDecIm,floor(ipDecIm/2)); 
    % mapping the Gray coded symbols into constellation
    modRe = alphaRe(ipGrayDecRe+1);
    modIm = alphaIm(ipGrayDecIm+1);
    % complex constellation
    mod = modRe + 1i*modIm;
    s = k_QAM*mod; % normalization of transmit power to one 
    
    % noise
    % -----
    n = 1/sqrt(2)*(randn(1,NumBit) + 1i*randn(1,NumBit)); % white guassian noise, 0dB variance 
    y = s + 10^(-Es_N0_dB(k)/20)*n; % additive white gaussian noise
    
    % demodulation
    % ------------
    y_re = real(y)/k_QAM; % real part
    y_im = imag(y)/k_QAM; % imaginary part

    % rounding to the nearest alphabet
    ipHatRe = 2*floor(y_re/2)+1;
    ipHatRe(ipHatRe>max(alphaRe)) = max(alphaRe);
    ipHatRe(ipHatRe<min(alphaRe)) = min(alphaRe);
    ipHatIm = 2*floor(y_im/2)+1;
    ipHatIm(ipHatIm>max(alphaIm)) = max(alphaIm);
    ipHatIm(ipHatIm<min(alphaIm)) = min(alphaIm);

    % Constellation to Decimal conversion
    ipDecHatRe = ind(floor((ipHatRe+const_size)/2+1))-1; % LUT based
    ipDecHatIm = ind(floor((ipHatIm+const_size)/2+1))-1; % LUT based

    % converting to binary string
    ipBinHatRe = dec2bin(ipDecHatRe,log2(M_qam)/2);
    ipBinHatIm = dec2bin(ipDecHatIm,log2(M_qam)/2);

    % converting binary string to number
    ipBinHatRe = ipBinHatRe.';
    ipBinHatRe = ipBinHatRe(1:end).';
    ipBinHatRe = reshape(str2num(ipBinHatRe).',log2(M_qam)/2,NumBit).' ;
    
    ipBinHatIm = ipBinHatIm.';
    ipBinHatIm = ipBinHatIm(1:end).';
    ipBinHatIm = reshape(str2num(ipBinHatIm).',log2(M_qam)/2,NumBit).' ;

    % counting errors for real and imaginary
    nBitErr(k) = size(find((ipBitRe- ipBinHatRe)),1) + size(find((ipBitIm - ipBinHatIm)),1) ;
    
    BERsim_mqam(k) = nBitErr(k)/(NumBit*log2(M_qam));    
    
    %BERthe_mqam(k) = (4/M)*(1-1/sqrt(M_qam))*(1/2)*erfc(sqrt(3*M*SNR_mqam/(M_qam-1))/sqrt(2));
    
%     if(M_qam == 16)
%     %BERthe_mqam(k) = (1/log2(M_qam))*3/2*erfc(sqrt(log2(M_qam)*(1/scaling_factor)*(SNR_mqam)));
%     BERthe_mqam(k) = (4/M)*(1-1/sqrt(M_qam))*(1/2)*erfc(sqrt(3*M*SNR_mqam/(M_qam-1))/sqrt(2));
%     end
%     if(M_qam == 64)
%     %BERthe_mqam(k) = (1/log2(M_qam))*7/4*erfc(sqrt((log2(M_qam)*(1/scaling_factor))*(SNR_mqam)));
%     BERthe_mqam(k) = (4/M)*(1-1/sqrt(M_qam))*(1/2)*erfc(sqrt(3*M*SNR_mqam/(M_qam-1))/sqrt(2));
%     end
    
end 

BERthe_mqam=berawgn(SNRdB,'qam',M_qam,'nondiff');

end







