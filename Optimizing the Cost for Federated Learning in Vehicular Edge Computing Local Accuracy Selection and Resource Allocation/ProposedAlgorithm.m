function [C, E, Time] = ProposedAlgorithm(vehicle_num, server_num, N_m, B, L, E_max, delta, D, beta_t, beta_e, c_n, v, g, k_n)
epislon = 0.0001;
[t, T] = deal(1, 20);
[C, E, Time] = deal(zeros(1, T), zeros(1, T), zeros(1, T));
[temp_freq, temp_theta] = deal(zeros(T, vehicle_num), zeros(T, vehicle_num));
[temp_rate, temp_x] = deal(zeros(vehicle_num, server_num, T), zeros(vehicle_num, server_num, T));
temp_freq(1, :) = randi([2e8, 1e9], 1, vehicle_num);
temp_theta(1, :) = rand(1, vehicle_num);
[temp_x(:, :, 1), temp_rate(:, :, 1)] = Algorithm1(temp_theta(1, :), temp_freq(1, :), vehicle_num, server_num, N_m, B, L, E_max, delta, D, beta_t, beta_e, c_n, v, g);
C(1) = compute_C(temp_theta(1, :), temp_freq(1, :), temp_rate(:, :, 1), temp_x(:, :, 1), vehicle_num, server_num, B, delta, D, beta_t, beta_e, c_n, g, k_n);
while t < T
    t = t + 1;
    for i = 1:vehicle_num
        f1 = L / v(i);
        f2 = power(beta_t / (2 * beta_e * k_n), 1/3);
        for j = 1:server_num
            f1 = f1 - temp_x(i, j, t-1) * delta(i) / temp_rate(i, j, t-1);
        end
        f1 = log2(1 / temp_theta(t-1, i)) * c_n(i) * D(i) / f1;
        temp_freq(t, i) = max(f1, f2);
    end
    temp_theta(t, :) =  Algorithm2(temp_freq(t-1, :), temp_x(:, :, t-1), temp_rate(:, :, t-1), vehicle_num, server_num, B, L, delta, D, beta_t, beta_e, c_n, v, g, k_n);
    [temp_x(:, :, t), temp_rate(:, :, t)] = Algorithm1(temp_theta(t, :), temp_freq(t, :), vehicle_num, server_num, N_m, B, L, E_max, delta, D, beta_t, beta_e, c_n, v, g);
    [C(t), E(t), Time(t)] = compute_C(temp_theta(t, :), temp_freq(t, :), temp_rate(:, :, t), temp_x(:, :, t), vehicle_num, server_num, B, delta, D, beta_t, beta_e, c_n, g, k_n);
    if abs(C(t) - C(t - 1)) < epislon
        break
    end
end