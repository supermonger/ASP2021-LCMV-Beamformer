function [loc] = DOA2(R, last_loc, degree_range)
    theta_range = 90;
    decimal = 2;
    zero_index = theta_range*10^(decimal) + 1;
    last_loc_index = zero_index + last_loc*10^decimal;
    P = zeros(1,2*degree_range*10^(decimal)+1);
    count = 1;
    for theta = last_loc_index-degree_range*10^(decimal):last_loc_index+degree_range*10^(decimal)
        size_R = size(R);
        a = a_of_theta(last_loc + (theta-last_loc_index)/10^(decimal),size_R(1));
        P_inv = a'/R*a;
        P(count) = 1/P_inv;
        count = count + 1;
    end
    r_P = real(P);
    [~, locs] = findpeaks(r_P,"SortStr","descend");
    if ~isempty(locs)
        locs = last_loc + (locs + last_loc_index-degree_range*10^(decimal) - 1 - last_loc_index)/10^(decimal);
        loc = locs(1);
    else
        loc = last_loc;
    end
    
end