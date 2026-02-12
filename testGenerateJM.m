function testGenerateJM()
thr = 1;


pose_num =10;

A = rand(pose_num * 6, pose_num * 6);
vec = diag(A) * 10 ;
diag_indices = GetDiagInd(A);
A(diag_indices) = vec;

% [Q, R] = QrHouseholder(A);
[Q,~, J] = QrGivens(A);


mat_dim = size(J, 2);
M_inv = zeros(mat_dim, mat_dim);
d = diag(J);
d(abs(d) < thr) = thr;
M_inv(diag_indices) = (1.0 ./ d);

for  j = pose_num : -1: 1
    for  i = j - 1 : -1:1
        M_temp = M_inv((i-1) * 6+1:i * 6, (i-1) * 6+1:i * 6);
        diag_indices = GetDiagInd(M_temp);
        Jii_inv = M_temp(diag_indices);
        
        tmp_vec =zeros(6,1);
        for (k = i + 1: j)
            J_temp = J((i-1) * 6+1:i*6, (k-1) * 6+1:k*6);
            diag_indices_J = GetDiagInd(J_temp);
            M_temp = M_inv((k-1) * 6+1:k*6, (j-1) * 6+1:j*6);
            diag_indices_M = GetDiagInd(M_temp);
            tmp_vec = tmp_vec+ (J_temp(diag_indices_J) .* M_temp(diag_indices_M));
        end
        
        M_temp = M_inv((i-1) * 6+1:i*6, (j-1) * 6+1:j*6);
        diag_indices_M = GetDiagInd(M_temp);
        M_temp(diag_indices_M) = -Jii_inv .* tmp_vec;
        M_inv((i-1) * 6+1:i*6, (j-1) * 6+1:j*6) = M_temp;
    end
end

JM_inv = zeros(mat_dim, mat_dim);

for ( i = 1: 1:pose_num)
    for ( j = i:pose_num)
        cur_mat = JM_inv((i-1) * 6+1:i*6, (j-1) * 6+1:j*6);
        for ( k = i:j)
            M_temp =  M_inv((k-1) * 6+1:k*6, (j-1) * 6+1:j*6);
            ind = GetDiagInd(M_temp);
            cur_mat = cur_mat + J((i-1) * 6+1:i*6, (k-1) * 6+1:k*6).* repmat(M_temp(ind)', 6, 1);
        end
        JM_inv((i-1) * 6+1:i*6, (j-1) * 6+1:j*6) = cur_mat;
    end
end

for (i = 1:mat_dim)
    max([norm(JM_inv(:,i)) thr]);
    col_norm = max([norm(JM_inv(:,i)) thr]);
    JM_inv(:,i) = JM_inv(:,i)./col_norm;
    M_inv(:,i) = M_inv(:,i)./col_norm;
end


end
function diag_indices = GetDiagInd(A)
n = size(A, 1);  % »ñÈ¡¾ØÕó½×Êý
diag_indices = [1:n+1:n^2]';
end
