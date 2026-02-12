function [G_big, Q, R, A] = QrGivens(A)
AA = A;
[Q_gt, R_gt] = qr(A);
G = [];
Q = [];
R = A;
rows = size(R, 1);
cols = size(R, 2);

%assert(rows >= cols);
if rows < cols
    cols = rows;
end
R = R(:,1:cols);
Gs = {};
G_big = eye(rows);
for col = 1 : cols
    for row = rows : -1 : col+1
        diag_ele = R(col, col);
        ele = R(row, col);
        if ele == 0
            continue;
        end
        G = eye(rows);
        r = sqrt(diag_ele^2 + ele^2);
        c = diag_ele / r;
        s = ele / r;
        Q = [c -s; s c];
        R([col row],:) = Q' * R([col row],:);
        G(col, col) = c;
        G(col, row) = -s;
        G(row, col) = s;
        G(row, row) = c;
        G_big = G_big * G;
        Gs = [Gs;G];
    end
end
if size(AA, 2) ~= cols
    R_check = G_big' * AA;
    R = [R R_check(:,end-((size(AA, 2) - cols) - 1):end)];
end
% for i = 1 : length(Gs)
%     G_big = Gs{i} * G_big;
% end

% G_big = G_big';
R_check = G_big' * A;
R_diff = R_check - R;

diff = R' * R - AA' * AA;
max(abs(diff(:)))
diff_gt = R_gt' * R_gt - AA' * AA;
max(abs(diff_gt(:)))
end