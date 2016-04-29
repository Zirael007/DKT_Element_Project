function f = El_Force(nnel,shape,P)
    fef = shape*P;
    for i = 1:nnel
        i1=(i-1)*3+1;
        i2=i1+1;
        i3=i2+1;
        f(i1,1) = fef(1);
        f(i2,1) = 0;
        f(i3,1) = 0;
    end
end