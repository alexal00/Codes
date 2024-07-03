function y = blade(R,y0,n)
    % Discretization of the blade
    lb = R-y0; % Lenght of the actual blade
    % n = 100;
    deltay = lb/n;
    
    y =zeros(1,n);
    for j= 0:n-1
        y(j+1)=y0+deltay/2+j*deltay;
    end

end
