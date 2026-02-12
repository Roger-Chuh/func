function [Q, R] = QrHouseholder(A)
[Q_gt, R_gt] = qr(A);
AA = A;
Q = [];
R = A;
rows = size(R, 1);
cols = size(R, 2);
%assert(rows >= cols);
if rows < cols
    cols = rows;
end
R = R(:,1:cols);
Q = eye(rows);
offset = 0;
for col = 1 : cols
    diag_ele = R(col - offset, col);
    col_norm = norm(R(col - offset:rows,col));
    e = zeros(rows - col + 1 + offset,1);
    sign_diag = sign(diag_ele);
    if (abs(e(1)) < 1e-10)
        sign_diag = 1;
    end
    e(1) = 1 * sign_diag * col_norm;
    if (abs( R(col - offset:rows,col)) < 1e-10)
        %         offset = offset + 1;
        %        continue;
    end
    v = R(col - offset:rows,col) + e;
    v_normalized = v./norm(v);
    
    %     v_normalized * v_normalized' - (v * v') / (v' * v)
    
    H = eye(rows - col + 1 + offset) - 2 * (v * v') / (v' * v);
    R(col - offset:rows,:) = H * R(col - offset:rows,:);
    Q = Q * [eye(col-1) zeros(col-offset-1, rows -col + 1);
        zeros(rows - col + 1 + offset, col-1) H];
end

if size(AA, 2) ~= cols
    R_check = Q' * AA;
    R = [R R_check(:,end-((size(AA, 2) - cols) - 1):end)];
end

Q' * AA;
diff = R' * R - AA' * AA;
max(abs(diff(:)))
diff_gt = R_gt' * R_gt - AA' * AA;
max(abs(diff_gt(:)))
end