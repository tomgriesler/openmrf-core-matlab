function [Tukey_Window] = mg_get_tukey_window(yMatrix, xMatrix, tukey_R, MOD)

%% Modus: Kreis
if MOD==0
    D = min([xMatrix yMatrix]);
    Tukey_Window_D = zeros(D,D);
    tukey_fun = tukeywin(D,tukey_R);
    tukey_fun = tukey_fun(D/2+1:D);
    x = linspace(-D/2, D/2, D);
    y = linspace(-D/2, D/2, D);
    for i=1:D
    for j=1:D   
      if (round(sqrt(x(i)^2 + y(j)^2)) <= D/2)
          Tukey_Window_D(i,j) = tukey_fun(round(sqrt(x(i)^2+y(j)^2)));
      end
    end                         
    end
    
[X,Y] = meshgrid(1:D,1:D);
[Xq,Yq] = meshgrid(1:xMatrix,1:yMatrix);
X = X./max(max(X));
Y = Y./max(max(Y));
Xq = Xq./max(max(Xq));
Yq = Yq./max(max(Yq));
Tukey_Window = interp2(X,Y,Tukey_Window_D,Xq,Yq);
Tukey_Window(isnan(Tukey_Window))=0;     
    
end

%% Modus: Rechteck
if MOD==1
    tukey_X = zeros(yMatrix, xMatrix);
    tukey_Y = zeros(yMatrix, xMatrix);
    for i=1:yMatrix
        tukey_X(i,:) = tukeywin(xMatrix, tukey_R);
    end
    for i=1:xMatrix
        tukey_Y(:,i) = tukeywin(yMatrix, tukey_R);
    end
    Tukey_Window = tukey_X .* tukey_Y;    
end

end



