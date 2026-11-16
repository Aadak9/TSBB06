%%Linear Triangulation
function xb = prepQ12(y1, y2, C1, C2)
    y1_tri = liu_crossop(y1);
    y2_tri = liu_crossop(y2);
    A = [y1_tri * C1; 
        y2_tri * C2];
    [~, ~, V]  = svd(A);
    xb = V(:, end);
    xb = xb / xb(4);
end

    