function [theta] = Algorithm2(freq, x, rate, vehicle_num, server_num, B, L, delta, D, beta_t, beta_e, c_n, v, g, k_n)
theta = zeros(1, vehicle_num);
epislon = 0.01;
T = 100;
gama = zeros(T,vehicle_num);
for i = 1:vehicle_num
    t = 1;
    tao = L / v(i);
    A = beta_t * c_n(i) * D(i) / freq(i) + beta_e * k_n * c_n(i) * D(i) * freq(i)^2;
    [C_1, C_2] = deal(0, 0);
    for j = 1:server_num
        C_1 = x(i, j) * delta(i) / rate(i, j);
        C_2 = x(i, j) * delta(i) * (2^(rate(i,j) / B) - 1) / (rate(i, j) * g(i, j));
    end
    C = beta_t * C_1 + beta_e * C_2;
    while t < T
        temp1 = gama(t, i) * A;
        temp2 = 0;
        for j = 1:server_num
            temp2 = temp2 + x(i, j) * delta(i) / rate(i, j);
        end
        temp2 = 2 ^ ((temp2 - tao) * freq(i) / (c_n(i) * D(i)));
        theta(i) = max(temp1, temp2);
        theta(i) = min(theta(i), 1);
        Cost = log2(1 / theta(i)) * A + C;
        if 1 - theta(i) - gama(t, i) * Cost <= epislon
            for l = t+1: T
                gama(l, i) = gama(t, i);
            end
            break
        else
            t = t + 1;
            gama(t, i) = (1 - theta(i)) / Cost;
        end
    end

end
