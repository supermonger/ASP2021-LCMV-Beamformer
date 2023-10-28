clear;clc;
load('ASP_Final_Data.mat')
%% parameter setting
[N,L] = size(matX);
R = zeros(N,N);
lambda = 0.99;
gi = 0.001;
theta_s_hat = zeros(1,L);
theta_i_hat = zeros(1,L);
y = zeros(1,L);
s_t_hat = zeros(1,L);
temp_s = 0;
%% DOA and beamformer algorithm
% source_loc = 2;
% interference_loc = -13;
tic
for i = 1:L
    R = lambda*R + matX(:,i)*matX(:,i)';
    R_of_i = R + 0.0001*lambda^(i)*eye(N);
    if i <= 5
        [locs] = DOA(R_of_i);
        source_loc = locs(1);
        interference_loc = locs(2);
    else
        [source_loc] = DOA2(R_of_i,source_loc, 1);
        [interference_loc] = DOA2(R_of_i,interference_loc, 1);
    end
    theta_s_hat(i) = source_loc;
    theta_i_hat(i) = interference_loc;
    W = LCMV(R_of_i,source_loc,interference_loc,gi,N);
    s_t_hat(i) = temp_s*0.5+W'*matX(:,i)*0.5;
    y(i) = W'*matX(:,i);
    temp_s = s_t_hat(i);
end
toc
%% plot source and interference theta
figure(1)
plot(theta_s_hat)
hold on
plot(theta_i_hat)
legend('\theta_s','\theta_i')
title('\theta(t) Estimation')
%% plot LCMV output
figure(2)
plot(abs(y))
title('|s(t)| of LCMV beamformer')
xlabel('time (t)')
ylabel('|s(t)|')
figure(3)
plot(unwrap(angle(y)))
title('angle(s(t)) of LCMV beamformer')
xlabel('time (t)')
ylabel('angle(s(t))')
%% plot my beamformer
figure(4)
plot(abs(s_t_hat))
ylim([0 6])
title('|s(t)| of My Beamformer')
xlabel('time (t)')
ylabel('|s(t)|')
figure(5)
plot(unwrap(angle(s_t_hat)))
title('angle(s(t)) of My beamformer')
xlabel('time (t)')
ylabel('angle(s(t))')
%% output .mat file
save('theta_s_hat.mat','theta_s_hat')
save('theta_i_hat.mat','theta_i_hat')
save('s_t_hat.mat','s_t_hat')