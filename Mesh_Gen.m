function [coords, connectivity] = Mesh_Gen(l, b, corner, seedx, seedy)

    % corner : Assumed to be bottom left corner of rectangle

    coords = [];
    connectivity = [];
    flag = 0;
    for j = corner(2):b/seedy:corner(2)+b
        for i = corner(1):l/seedx:corner(1)+l
            coords = [coords; i, j];
            flag = flag+1;
        end
    end

    x = 1:1:seedx+1;
    y = 1:1:seedy+1;
    
    for i = 1:1:(seedx+1)*seedy-1
        if ~~mod(i,seedx+1)
            connectivity = [connectivity; i i+1 i+seedx+2 i+seedx+1];
        else
            continue;
        end
    end
    
end