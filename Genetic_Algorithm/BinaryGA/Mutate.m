function y = Mutate(x, mu)

    flag = (rand(size(x)) < mu);

    y = x;
    y(flag) = 1 - y(flag);

end