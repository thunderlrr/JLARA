function [x, rate] = Algorithm1(theta, freq, vehicle_num, server_num, N_m, B, L, E_max, delta, D, beta_t, beta_e, c_n, v, g)
[t, T]= deal(1, 100);
precision = 0.1;
lambda = ones(T, vehicle_num) * 10;
mu = ones(T, vehicle_num) * 2;
phi = ones(T, server_num) * 1;
[dlambda, dmu, dphi] = deal(zeros(1, vehicle_num), zeros(1, vehicle_num), zeros(1, server_num));
H_nm = zeros(vehicle_num, server_num);

while t <= T
    rate = zeros(vehicle_num, server_num);
    x = zeros(vehicle_num, server_num);
    a = 0.1 / t;
    for i = 1:vehicle_num
        for j = 1:server_num
            temp = ((lambda(t, i) + beta_t / (1-theta(i))) * g(i, j) - beta_e / (1-theta(i)) - mu(t, i)) / ((beta_e / (1-theta(i)) + mu(t, i))* exp(1));
            rate(i, j) = B * (lambertw(temp) + 1) / log(2);
            H_nm(i, j) = (beta_e / (1-theta(i)) + mu(t, i)) * delta(i) * log(2) * 2^(rate(i, j) / B) / (B * g(i, j));
        end
    end
    for i = 1:vehicle_num
        [minH, n, m] = deal(max(H_nm(i, :)), i, 1);
        for j = 1:server_num
            if (H_nm(i, j) <= minH) && (sum(x(:, j)) < N_m(j))
                minH = H_nm(i, j);
                [n, m]= deal(i, j);
            end
        end
        x(n, m) = 1;
    end
    for i = 1:vehicle_num
        dlambda(i) = log2(1 / theta(i)) * c_n(i) * D(i) / freq(i) - L / v(i);
        dmu(i) = -E_max;
        for j = 1:server_num
            dlambda(i) = dlambda(i) + x(i, j) * delta(i) / rate(i, j);
            dmu(i) = dmu(i) + x(i, j) * delta(i) * (2^(rate(i, j) / B) - 1) / (rate(i, j) * g(i, j));
        end
        lambda(t + 1, i) = max(0, lambda(t, i) - a * dlambda(i));
        mu(t + 1, i) = max(0, mu(t, i) - a * dmu(i));
    end
    for j = 1:server_num
        dphi(j) = N_m(j);
        for i = 1:vehicle_num
            dphi(j) = dphi(j) - x(i, j);
        end
        phi(t + 1, j) = max(0, phi(t, j) - a * dphi(j));
    end
    if (norm(lambda(t + 1, :)-lambda(t, :)) < precision) && (norm(phi(t + 1, :)-phi(t, :)) < precision) 
        for l = t+1: T
            lambda(l, :) = lambda(t, :);
            mu(l, :) = mu(t, :);
        end
        break;
    end
    t = t + 1;
end