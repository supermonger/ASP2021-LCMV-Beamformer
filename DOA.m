function [locs] = DOA(R)
    theta_range = 90;
    decimal = 1;
    list_range = theta_range*2*10^decimal + 1;
    zero_index = theta_range*10^decimal + 1;
    P = zeros(1,list_range);
    for theta = 1:list_range
        size_R = size(R);
        a = a_of_theta((theta-zero_index)/10^(decimal),size_R(1));
        P_inv = a'/R*a;
        P(theta) = 1/P_inv;
    end
    r_P = real(P);
    [~, locs] = findpeaks(r_P,"SortStr","descend");
    locs = (locs-zero_index)/10^(decimal);
    locs = locs(1:2);
end