function at = mytranspose(a)
   [rows, columns] = size(a);
   at = zeros(columns, rows);
   
   for i = 1:rows
        for j = 1:columns
            at(j, i) = a(i,j);
        end
   end
end