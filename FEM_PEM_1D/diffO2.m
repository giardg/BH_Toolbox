function [vX,aX]=diffO2(t,X)
h=(max(t)-min(t))/(length(t)-1);
n=length(X);
vX=zeros(size(X));
aX=zeros(size(X));

for i=1:n
    if i==1
        vX(i)=-3*X(i)+4*X(i+1)-1*X(i+2);
        aX(i)=2*X(i)-5*X(i+1)+4*X(i+2)-1*X(i+3);
    elseif i==n
        vX(i)=3*X(i)-4*X(i-1)+1*X(i-2);
        aX(i)=2*X(i)-5*X(i-1)+4*X(i-2)-1*X(i-3);
    else
        vX(i)=-1*X(i-1)+1*X(i+1);
        aX(i)=1*X(i-1)-2*X(i)+1*X(i+1);
    end
end

vX = vX/(2*h);
aX = aX/h^2;