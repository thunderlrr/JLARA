function [C,e,t] = compute_C(theta, freq, rate, x, vehicle_num, server_num, B, delta, D, beta_t, beta_e, c_n, g, k_n)
C = 0;
[E, T] = deal(zeros(1, vehicle_num), zeros(1, vehicle_num));
for i = 1:vehicle_num
    [C_t, C_e] = deal(log2(1/theta(i)) * c_n(i) * D(i) / freq(i), ...
        log2(1/theta(i)) * k_n * c_n(i) * D(i) * freq(i) ^ 2);
    E(i) = E(i) + log2(1 / theta(i)) * k_n * c_n(i) * D(i) * freq(i) ^ 2;
    T(i) = T(i) + log2(1 / theta(i)) * c_n(i) * D(i) / freq(i);
    for j = 1:server_num
        C_t = C_t + x(i, j) * delta(i) / rate(i, j);
        C_e = C_e + x(i, j) * delta(i) * (2 ^ (rate(i, j) / B) - 1) / (rate(i, j) * g(i, j));
        E(i) = E(i) + x(i, j) * delta(i) * (2 ^ (rate(i, j) / B) - 1) / (rate(i, j) * g(i, j));
        T(i) = T(i) + x(i, j) * delta(i) / rate(i, j);
    end
    C = C + (beta_t * C_t + beta_e * C_e) / (1 - theta(i));
    E(i) = E(i) / (1 - theta(i));
    T(i) = T(i) / (1 - theta(i));
    
end
e = sum(E);
t = sum(T);