function [BERsim_alamouti1, BERsim_alamouti2] = Alamouti_T2R1_T2R2()
NumBit = 10^6; % number of bits or symbols
SNRdB = 0:1:45; % multiple Eb/N0 values
nRx = 2;
for k = 1:length(SNRdB)
    
    % Transmitter
    s_bit = rand(1,NumBit)>0.5; % generating 0,1 with equal probability
    s = 2*s_bit-1; % BPSK modulation 0 -> -1; 1 -> 0
    
    % Alamouti STBC
    sCode = 1/sqrt(2)*kron(reshape(s,2,NumBit/2),ones(1,2)) ;
    % sCode
    %  [x1 x1
    %   x2 x2 ];
    % channel
    h = 1/sqrt(2)*(randn(nRx,NumBit) + 1i*randn(nRx,NumBit)); % Rayleigh channel
    
    n = 1/sqrt(2)*(randn(nRx,NumBit) + 1i*randn(nRx,NumBit)); % white gaussian noise, 0dB variance
    
    y = zeros(nRx,NumBit);
    yMod = zeros(nRx*2,NumBit);
    hMod = zeros(nRx*2,NumBit);
    %     N=6
    %     h= [1+i 2+i 3+i 4+i 5+i 6+i];
    
    % for first receiver
    hMod = kron(reshape(h(1,:),2,NumBit/2),ones(1,2)); % repeating the same channel for two symbols
    
    %  [h11 h11]
    %  [h12 h12];
    temp = hMod;
    hMod(1,[2:2:end]) = conj(temp(2,(2:2:end)));
    hMod(2,[2:2:end]) = -conj(temp(1,(2:2:end)));
    
    y(1,:) = sum(sCode.*hMod,1) + 10^(-SNRdB(k)/20)*n(1,:);

    % Receiver
    yMod([1:2],:) = kron(reshape(y(1,:),2,NumBit/2),ones(1,2));

    
    hEq((1:2),:) = hMod;
    hEq(1,(1:2:end)) = conj(hEq(1,(1:2:end)));
    hEq(2,(2:2:end)) = conj(hEq(2,(2:2:end)));
    
    hEqPower = sum(hEq.*conj(hEq),1);
    yHat = sum(hEq(1:2,:).*yMod(1:2,:),1)./hEqPower;
    yHat(2:2:end) = conj(yHat(2:2:end));
    
    % receiver - hard decision decoding
    ipHat = real(yHat)>0;
    
    % counting the errors
    nErr1(k) = size(find((s_bit- ipHat)),2);

    % for second receiver
    hMod = kron(reshape(h(2,:),2,NumBit/2),ones(1,2)); % repeating the same channel for two symbols
    %  [h21 h21]
    %  [h22 h22];
    
    temp = hMod;
    hMod(1,[2:2:end]) = conj(temp(2,[2:2:end]));
    hMod(2,[2:2:end]) = -conj(temp(1,[2:2:end]));
    
    % Channel and noise Noise addition
    y(2,:) = sum(hMod.*sCode,1) + 10^(-SNRdB(k)/20)*n(2,:);

    % Receiver
    yMod([3:4],:) = kron(reshape(y(2,:),2,NumBit/2),ones(1,2));
    
    % forming the equalization matrix
    hEq([3:4],:) = hMod;
    %  [h21  h22*]
    %  [h22 -h21*];
    
    hEq(3,[1:2:end]) = conj(hEq(3,[1:2:end]));
    hEq(4,[2:2:end]) = conj(hEq(4,[2:2:end]));
    
    % equalization
    hEqPower = sum(hEq.*conj(hEq),1);
    yHat = sum(hEq.*yMod,1)./hEqPower;
    yHat(2:2:end) = conj(yHat(2:2:end));
    
    % receiver - hard decision decoding
    ipHat = real(yHat)>0;
    
    % counting the errors
    nErr2(k) = size(find([s_bit- ipHat]),2);
    
end
BERsim_alamouti1 = nErr1/NumBit;
BERsim_alamouti2 = nErr2/NumBit; % simulated ber
end




