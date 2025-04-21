function H = FillMatrix(H, start_row, start_col, row_size, col_size, value)

if (start_row < 1 || start_col < 1)
    return
end

H(start_row:start_row+row_size-1, start_col:start_col+col_size-1) = H(start_row:start_row+row_size-1, start_col:start_col+col_size-1) + value;
end