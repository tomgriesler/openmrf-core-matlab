function [phi] = mg_phispace(phi1, phi2, n)
    phi = linspace(phi1, phi2, n+1);
    phi(end) = [];
end

