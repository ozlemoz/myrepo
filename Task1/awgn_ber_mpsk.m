
function [BERsim_mpsk, BERthe_mpsk] = awgn_ber_mpsk(NumBit,M,SNRdB)

M_mpsk = 2^M;   % constellation size
p = log2(M_mpsk); % bits per symbol

thetaMpsk = (0:M_mpsk-1)*2*pi/M_mpsk; % reference phase values

Es_N0_dB  = SNRdB + 10*log10(p);
ss = pskmod(0 : M_mpsk-1, M_mpsk, 0, 'Gray');
Es = ss * ss'/ M_mpsk;
Eb = Es/log2(M_mpsk);

% Mapping for binary <--> Gray code conversion
ref = (0:M_mpsk-1); 
map = bitxor(ref,floor(ref/2));
[~, ind] = sort(map);                                

ipPhaseHat = zeros(1,NumBit);
for k = 1:length(SNRdB)
    
    SNR_mpsk = 10.^(k/10);
    
    % symbol generation
    % ------------------
    ipBit = rand(1,NumBit*p,1)>0.5; % random 1's and 0's
    ipBitReshape =  reshape(ipBit,p,NumBit).'; % grouping to N symbols having k bits each
    ipGray = bi2de(ipBitReshape,'left-msb')';
    
    
    % Gray coded constellation mapping
    ipDec = ind(ipGray+1)-1; % bit group to constellation point 
    ipPhase = ipDec*2*pi/M_mpsk; % conversion to phase 
    ip = exp(1i*ipPhase);  % modulation 
    s = ip; 
    
    % noise
    % -----
    n = 1/sqrt(2)*(randn(1,NumBit) + 1i*randn(1,NumBit)); % white guassian noise, 0dB variance 
    y = s + 10^(-Es_N0_dB(k)/20)*n; % additive white gaussian noise
    
    % demodulation
    % ------------
    % finding the phase from [-pi to +pi]
    opPhase = angle(y); 
    % unwrapping the phase i.e. phase less than 0 are 
    % added 2pi
    opPhase(opPhase<0) = opPhase(opPhase<0) + 2*pi;

    % rounding the received phase to the closest constellation
    ipPhaseHat = 2*pi/M_mpsk*round(opPhase/(2*pi/M_mpsk))	;
    % as there is phase ambiguity for phase = 0 and 2*pi,
    % changing all phases reported as 2*pi to 0.
    % this is to enable comparison with the transmitted phase
    ipPhaseHat(ipPhaseHat==2*pi) = 0;
    ipDecHat = round(ipPhaseHat*M_mpsk/(2*pi));
 
    % Decimal to Gray code conversion
    ipGrayHat = map(ipDecHat+1); % converting to decimal 
    ipBinHat1 = de2bi(ipGrayHat','left-msb',p)';
    ipBinHat2 =  reshape(ipBinHat1,1,[]);
%     ipBinHat = dec2bin(ipGrayHat,k) ; % decimal to binary

%     % converting binary string to number
%     ipBinHat = ipBinHat.';
%     ipBinHat = ipBinHat(1:end).';
%     ipBinHat = str2num(ipBinHat).' ;
    
    % counting errors
    nBitErr(k) = size(find((ipBit- ipBinHat2)),2); % couting the number of errors
    BERsim_mpsk(k) = nBitErr(k)/(NumBit*p);
    

end 

 BERthe_mpsk=berawgn(SNRdB,'psk',M_mpsk,'nondiff');
 
end