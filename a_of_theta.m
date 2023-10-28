function a = a_of_theta(theta,N)
    a = (1:N)';
    a = exp(1i*pi*sind(theta)*(a-1));
end