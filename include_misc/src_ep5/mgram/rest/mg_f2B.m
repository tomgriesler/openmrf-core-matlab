function [B] = mg_f2B(f)
    % w = gamma * B
    % f = gamma * B /2/pi
    % B = w / gamma
    % B = 2*pi* f / gamma
    gamma = 2.6752218744 *1e8;
    B     = 2*pi* f / gamma;
end

