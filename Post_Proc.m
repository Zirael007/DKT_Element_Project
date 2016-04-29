function [w, thetax, thetay] = Post_Proc(displacement)

    w = displacement(1:3:end);
    thetax = displacement(2:3:end);
    thetay = displacement(3:3:end);

end