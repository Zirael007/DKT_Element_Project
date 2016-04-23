
nnodes
coords
x = coords(:,1);
y = coords(:,2);
if nnodes == 6
    N = Tri_shape_fxn(x,y);
elseif nnodes == 8
    N = Quad_shape_fxn(x,y);
end
b_xi(5) = 0.5*(b_xi(4) + b_xi(1));
b_xi(6) = 0.5*(b_xi(1) + b_xi(2));
b_xi(7) = 0.5*(b_xi(2) + b_xi(3));
b_xi(8) = 0.5*(b_xi(3) + b_xi(4));
b_x = N*b_xi';
b_y = N*b_yi';
Hx_x
Hy_y
Hx_y
Hy_x
J
    