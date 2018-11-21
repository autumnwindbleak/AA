function [G1,S,T,W] = mrf1(G0,f11,f00,b)
f10 = 1 - f00;
f01 = 1 - f11;
w1 = log(f11/f10);
w0 = log(f01/f00);
[height,width] = size(G0);
imax = height * width;
count = 1;
source = width * height + 1;
sink = width * height + 2;
for i = 1 : 1 : height
    for j = 1 : 1 : width
        if G0(i,j) == 0
            ci = w0;
        else
            ci = w1;
        end
        if ci > 0
            S(count) = source;
            T(count) = (j - 1) * height + i; 
            W(count) = ci;
            count = count + 1;
        else
            S(count) = (j - 1) * height + i; 
            T(count) = sink;
            W(count) = -ci;
            count = count + 1;
        end
        if j+1 <= width
            S(count) = j * height + i;
            T(count) = (j - 1) * height + i; 
            W(count) = b;
            count = count + 1;
            S(count) = (j - 1) * height + i; 
            T(count) = j * height + i;
            W(count) = b;
            count = count + 1;
        end
        if i + 1 <= height
            if j - 1 > 0
                S(count) = (j - 1 - 1) * height + i + 1;
                T(count) = (j - 1) * height + i; 
                W(count) = b;
                count = count + 1;
                S(count) = (j - 1) * height + i; 
                T(count) = (j - 1 - 1) * height + i + 1;
                W(count) = b;
                count = count + 1;
            end
            S(count) = (j - 1) * height + i + 1;
            T(count) = (j - 1) * height + i; 
            W(count) = b;
            count = count + 1;
            S(count) = (j - 1) * height + i; 
            T(count) = (j - 1) * height + i + 1;
            W(count) = b;
            count = count + 1;
            if j + 1 <= width
                S(count) = j * height + i + 1;
                T(count) = (j - 1) * height + i; 
                W(count) = b;
                count = count + 1;
                S(count) = (j - 1) * height + i; 
                T(count) = j * height + i + 1;
                W(count) = b;
                count = count + 1;
            end
        end
    end
end
G = digraph(S,T,W);
[mf,~,cs,ct] = maxflow(G,source,sink);
G1 = ones(height,width);
for i = 1:1:length(ct)
    if ct(i) <= imax
        G1(ct(i)) = 0;
    end
end

