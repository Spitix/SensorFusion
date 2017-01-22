function result=h(x_k,s_1, s_k, h_0)
    u=double((x_k(1)-s_1(1))^2+(x_k(2)-s_1(2))^2+h_0^2);
    v=double((x_k(1)-s_k(1))^2+(x_k(2)-s_k(2))^2+h_0^2);
    result=u/v;
end