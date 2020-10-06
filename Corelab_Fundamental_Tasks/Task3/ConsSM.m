function [SER_opt_cons, SER_subopt_cons] = ConsSM()

NumBit = 10^6;
M = 2; % Modulation order
Nt = 4;  % Tx antennas
Nr = 4; % Rx antennas
antenna_index_num = log2(Nt);
symbol_index_num = log2(M);
bit_array_num = antenna_index_num + symbol_index_num;
symbol = [-sqrt(-1) sqrt(-1)];
antenna_index = 0;
SNRdB = 0:2:25;
SER1(1,length(SNRdB)) = zeros;
SER2(1,length(SNRdB)) = zeros;
SER_opt_cons(1,length(SNRdB)) = zeros;
SER4(1,length(SNRdB)) = zeros;
SER5(1,length(SNRdB)) = zeros;
SER_subopt_cons(1,length(SNRdB)) = zeros;
symbol_errors1 = zeros(length(SNRdB),1);
symbol_errors2 = zeros(length(SNRdB),1);
symbol_errors3 = zeros(length(SNRdB),1);
symbol_errors4 = zeros(length(SNRdB),1);
symbol_errors5 = zeros(length(SNRdB),1);
symbol_errors6 = zeros(length(SNRdB),1);
y = zeros(Nr,1);

for k = 1:length(SNRdB)
    
    sigma = 10^(-SNRdB(k)/20);
    for a=1:NumBit
        bit_array =[1 0 1];
        
        for j = 1 : antenna_index_num
            antenna_bits(j) = bit_array(j);
        end
        antenna_index = bi2de(antenna_bits,'left-msb') + 1;
        for j = antenna_index_num + 1 : bit_array_num
            signal_bits(j-antenna_index_num) = bit_array(j);
        end
        signal_index = bi2de(signal_bits,'left-msb') + 1;
        x = zeros(Nt, 1);
        x(antenna_index) = symbol(signal_index);
        
        % The channel
        n = (1/sqrt(2))*(randn(Nr,1)+sqrt(-1)*randn(Nr,1));   % white Gaussian Noise 0dB variance
        Hv = (1/sqrt(2))*(randn(Nr, Nt) + sqrt(-1)*randn(Nr, Nt)); % Conventional Channel
        
        Hc = [];
        for j = 1 : Nt
            h_norm = norm(Hv(:, j), 'fro');
            normalized_column = Hv(:, j) ./ h_norm;
            Hc = [Hc normalized_column]; % Constrained Channel
        end
        y = Hc*x +sigma*n;
        
        % ML   Optimal Detection Scheme
        p = 1;
        
        for j = 1:Nt
            for i = 1:M
                g = Hc(:,j)*symbol(i);
                g_f = 0;
                g_f = g_f + norm(g,'fro')^2;
                argmin(i,j) = sqrt(p)*g_f - 2*(real(y'*g));
            end
        end
        
        MinValue = min(min(argmin));
        [row,col] = find(MinValue == argmin);
        detected_antenna_index = col;
        detected_signal_index = row ;
        
        detected_antenna_bits = de2bi(detected_antenna_index - 1, antenna_index_num, 'left-msb');
        detected_signal_bits = de2bi(detected_signal_index - 1, symbol_index_num, 'left-msb');
        
        % error
        l = isequal(antenna_bits, detected_antenna_bits);
        if l == 0
            symbol_errors1(k) = symbol_errors1(k) + 1;
        end
        
        q = isequal(signal_bits, detected_signal_bits);
        if q == 0
            symbol_errors2(k) = symbol_errors2(k) + 1;
            
        end
         symbol_errors3(k) = (symbol_errors1(k) + symbol_errors2(k));
         
        detected_antenna_bits = [];
        detected_signal_bits = [];
        detected_antenna_index = [];
        detected_signal_index = [];

        % Sub Optimal(MRC) Detection Scheme
        for j = 1:Nt
            z(j) = Hc(:,j)'*y;
        end
        for i = 1:Nt
            z_norm(i) = norm(z(:,i));
        end
        for i = 1:Nt
            H_norm(i) = norm(Hc(:,i),'fro')^2;
        end
        
        Z = z./H_norm;
        Z_norm = z_norm./H_norm;
        
        [~,detected_antenna_index] = max(Z_norm);
        Z_dist_symbol = [];
        
        for j = 1 : M
            Z_dist_symbol(j) = norm(Z(detected_antenna_index) - symbol(j))^2;
        end
        [~, detected_signal_index] = min(Z_dist_symbol);
        detected_signal_bits = de2bi(detected_signal_index - 1, symbol_index_num, 'left-msb');
        detected_antenna_bits = de2bi(detected_antenna_index-1, antenna_index_num, 'left-msb');

        % error
        l = isequal(antenna_bits, detected_antenna_bits);
        if l == 0
            symbol_errors4(k) = symbol_errors4(k) + 1;
        end
        
        q = isequal(signal_bits, detected_signal_bits);
        if q == 0
            symbol_errors5(k) = symbol_errors5(k) + 1;
            
        end
        symbol_errors6(k) = (symbol_errors4(k) + symbol_errors5(k));
    end

    SER_opt_cons(k) = symbol_errors3(k)/(3*NumBit);
    symbol_errors3 = zeros(length(SNRdB),1);
     SER_subopt_cons(k) = symbol_errors6(k)/(3*NumBit);
    symbol_errors6 = zeros(length(SNRdB),1);
end
end
