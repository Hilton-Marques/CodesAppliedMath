function heav = heaviside(n, qsi, eta, d1, d2, d3, d4)
  heav = zeros(n, n);
  Phi = 0.5*(d2*(ones(n, n)+qsi)+d4*(ones(n, n)+eta)-d1*(eta+qsi));
  for i = 1:n
    for j = 1:n
      if Phi (i,j) > 0
        heav (i,j) = 1;
      elseif Phi (i,j) <= 0
        heav (i,j) = -1;
      end
    end
  end
end
