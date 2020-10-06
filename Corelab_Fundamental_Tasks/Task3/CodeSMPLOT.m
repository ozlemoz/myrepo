clc;
clear all;
close all;
tic;
M=2;
Nt=4;
Nr=4;

SigConsDiag = [-sqrt(-1) ;sqrt(-1)];
diag = SigConsDiag'*SigConsDiag/M;
SigConsDiag = SigConsDiag/sqrt(diag);
eta = log2(Nt) + log2(M);
SNR_Vector = 0:2:25;
ABER_Ana = zeros(size(SNR_Vector));
for snr = 1: length(SNR_Vector)
    sigma_n = (1/(10^(SNR_Vector(snr)/10)));
    for ell_t = 1:Nt
        for ell = 1:Nt
            for i_t = 1:M
                for i = 1:M
                    if ell_t ==ell
                        Psi = abs(SigConsDiag(i_t) - SigConsDiag(i))^2;
                        BarGamma_SM = (1/(2*sigma_n))*Psi;
                    else
                        Psi = abs(SigConsDiag(i_t))^2 +  abs(SigConsDiag(i))^2;
                        BarGamma_SM = (1/(2*sigma_n))*Psi;
                    end
                    alpha_a = (1/2)*(1-sqrt((BarGamma_SM/2)/(1+(BarGamma_SM/2))));
                    PEP = 0;
                    for nr = 0:(Nr-1)
                        PEP = PEP + nchoosek(Nr-1+nr,nr)*(1-alpha_a)^nr;
                    end
                    PEP = alpha_a^Nr*PEP;
                    bitError = (biterr(ell_t-1,ell-1,log2(Nt)) + biterr(i_t-1,i-1,log2(M)))/eta;
                    ABER_Ana(snr) = ABER_Ana(snr) + bitError * PEP/8;
                    
                end
            end
        end
    end
    ABER_Ana = ABER_Ana;
end

SNRdB = 0:2:25;
[SER_opt_conv, SER_subopt_conv] = ConvSM();
[SER_opt_cons, SER_subopt_cons] = ConsSM();
semilogy(SNRdB, SER_subopt_conv,'b-d');
hold on
semilogy(SNRdB, SER_subopt_cons, 'r:d');
semilogy(SNRdB, SER_opt_conv,'y-*');
semilogy(SNRdB, SER_opt_cons,'m:*');
semilogy(SNRdB, ABER_Ana, 'g--');
semilogy(SNRdB, [1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7], '-');
semilogy(SNRdB, [1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7], ':');
grid on
title('Spatial Modulation');
xlabel('\rho (dB)');
ylabel('P_{e,bit}');
legend('SM (Mesleh), N_{t}=4','SM (Mesleh), N_{t}=4','SM (Optimal), N_{t}=4','SM (Optimal), N_{t}=4','SM (Analytical), Nt=4','conventional H','constrained H');
axis([0 25 10^-5 1])
toc;