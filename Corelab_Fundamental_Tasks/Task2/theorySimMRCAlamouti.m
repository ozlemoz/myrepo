clear all;
close all;
clc;

Narray = 4;

tic

SNRdB = 0:1:45;

% Theory_MRC
BERthe_MRC = zeros(Narray,length(SNRdB));

for N = 1:Narray
    for i = 1:length(SNRdB)
        
        SNR(i) = 10.^(SNRdB(i)/10);
        p = 1/2 - (1/2)*(1+1./SNR(i)).^(-1/2);
        
        for k = 0:N-1
            Pe = p^N .* nchoosek(N-1+k,k).*(1-p).^k;
            BERthe_MRC(N,i) = BERthe_MRC(N,i) + Pe;
        end 
    end
end

% Alamouti
BERthe_Alamouti = zeros(Narray,length(SNRdB));

for N = 1:Narray
    for i = 1:length(SNRdB)
        
        SNR(i) = 0.5*10.^(SNRdB(i)/10);
        pA = 1/2 - (1/2)*(1+1./SNR(i)).^(-1/2); 
        
        for k = 0:N-1
            PeA = pA^N .* nchoosek(N-1+k,k).*(1-pA).^k;
            BERthe_Alamouti(N,i) = BERthe_Alamouti(N,i) + PeA;
        end
    end
end
[BERsim_MRC1, BERsim_MRC2] = MRC();
[BERsim_alamouti1, BERsim_alamouti2] = Alamouti_T2R1_T2R2();
figure
semilogy(SNRdB,BERthe_MRC(1,:),'k--');
hold on
semilogy(SNRdB,BERsim_alamouti1,'r');
semilogy(SNRdB,BERthe_Alamouti(2,:),'r*');
semilogy(SNRdB,BERsim_alamouti2,'b');
semilogy(SNRdB,BERthe_Alamouti(4,:),'bo');
semilogy(SNRdB,BERsim_MRC1,'g');
semilogy(SNRdB,BERthe_MRC(2,:),'g+');
semilogy(SNRdB,BERsim_MRC2,'m');
semilogy(SNRdB,BERthe_MRC(4,:),'mx');
axis([0 45 10^-5 1])
grid on
legend('No Diversity (1Tx - 1Rx', 'Alamouti Simulation (2Tx - 1Rx)', 'Alamouti Theoretical (2Tx - 1Rx)',...
    'Alamouti Theoretical (2Tx - 2Rx)', 'Alamouti Simulation (2Tx - 2Rx)', ...
    'MRC Simulation (1Tx - 2Rx)', 'MRC Theoretical (1Tx - 2Rx)', ...
    'MRC Simulation (1Tx - 4Rx)', 'MRC Theoretical (1Tx - 4Rx)');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('Alamouti vs Maximal Ratio Combining');
toc