function k=get_settle(x_state,x_t_vec,presision)
    k=NaN;
    tf=length(x_state);
    i=tf;
    while norm(x_state(:,i)-x_t_vec')< presision && i>1
        k=i;
        i=i-1;
    end
end