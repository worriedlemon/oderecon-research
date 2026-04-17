function c = mymult(a, b)
   [r1, c1] = size(a); [~, c2] = size(b);

   c = zeros(r1, c2);
   for i = 1:r1
        for j = 1:c2
            for k = 1:c1
                c(i, j) = c(i, j) + a(i, k) * b(k, j);
            end
        end
   end
end