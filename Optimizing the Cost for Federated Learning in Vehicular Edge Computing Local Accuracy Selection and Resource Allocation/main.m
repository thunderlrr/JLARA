clear;clc;
beta_t = 0.7;
beta_e = 1 - beta_t;
[vehicle_num, server_num] = deal(20, 5);
N_m = [4,4,4,4,4];
delta = 40 * 2 ^ 10 * 8 * ones(1, vehicle_num);
k_n = 5e-27;
B = 20e6;
N0 = -174;
c_n = randi([900, 1600], 1, vehicle_num);
D = randi([60, 150], 1, vehicle_num) * 2 ^ 10 * 8;
L = 100;
v = 25 + sqrt(7.5) .* randn(1, vehicle_num);
E_max = 2;
h_nm = zeros(vehicle_num, server_num);
pl_nm = zeros(vehicle_num, server_num);
g = zeros(vehicle_num, server_num);
for i = 1:vehicle_num
    for j = 1:server_num
        pl_nm(i, j) = 135 + rand() * 10;
        h_nm(i, j) = 10 ^ (-(pl_nm(i, j)-10) / 10);
        g(i, j) = h_nm(i, j) / (B * 10 ^ (N0/10) / N_m(j));
    end
end

[C, E, T] = ProposedAlgorithm(vehicle_num, server_num, N_m, B, L, E_max, delta, D, beta_t, beta_e, c_n, v, g, k_n);

