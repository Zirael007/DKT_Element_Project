function B=strain_displacement(jacob, H_x_eps, H_x_eta, H_y_eta, H_y_eps)
j11 = jacob(2,2)/det(jacob);
j12 = -jacob(1,2)/det(jacob);
j21 = -jacob(2,1)/det(jacob);
j22 = jacob(2,2)/det(jacob);

B=[j11*H_x_eps + j12*H_x_eta; 
    j21*H_y_eps+j22*H_y_eta;
    j11*H_y_eps+j12*H_y_eta+j21*H_x_eps+j22*H_x_eta];

end
