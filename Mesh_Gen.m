function [coords, connectivity] = Mesh_Gen(l, b, corner, seedx, seedy, elType)

    % corner : Assumed to be bottom left corner of rectangle

    coords = [];
    connectivity = [];

    for j = corner(2):b/seedy:corner(2)+b
        for i = corner(1):l/seedx:corner(1)+l
            coords = [coords; i, j];
        end
    end

    switch lower(elType)
        case 'quad'

    for i = 1:1:(seedx+1)*seedy-1
        if ~~mod(i,seedx+1)
            connectivity = [connectivity; i i+1 i+seedx+2 i+seedx+1];
        else
            continue;
        end
    end

        case 'tri'
        
            for i = 1:1:(seedx+1)*seedy-1
                if ~~mod(i,seedx+1)
                    connectivity = [connectivity; i i+1 i+seedx+1];
                    connectivity = [connectivity; i+1 i+seedx+2 i+seedx+1];
                else
                    continue;
                end
            end
    
    end
end