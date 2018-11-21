function points = simplex(c,A,b,x0)
c = c';
A = [A;-1 0;0 -1];
b = [b;0;0]; 

points = x0;
pointx = [x0(1)];
pointy = [x0(2)];
[M,d] = size(A);
intsx = [];
intsy = [];

wucha = 1e-10;


xinf = 0;
yinf = 0;

for i = 1:M-1
    for j = i+1:M
        subA = [A(i,1) A(i,2);A(j,1) A(j,2)];
        subb = [b(i);b(j)];
        point = subA\subb;
        if point(1) == Inf | point(2) == Inf
           
        else
            if xinf < point(1)
                xinf = point(1);
            end
            if yinf < point(2)
                yinf = point(2);
            end
            if A * [point(1);point(2)] - b <= wucha
                intsx = [intsx point(1)];
                intsy = [intsy point(2)];
            end
        end
    end
end

xinf = xinf + 1;
yinf = yinf + 1;

notunique = [];
for i = 1:length(intsx)
    for j =i + 1:length(intsx)
        if intsx(i) == intsx(j) & intsy(i) == intsy(j)
            notunique = [notunique;j];
        end
    end
end
notunique = unique(notunique);
intsx(notunique) = [];
intsy(notunique) = [];


if length(intsx) == 0
    points = NaN;
    return;
end
    if length(intsx) == 1
        points = [points [intsx(1);intsy(1)]];
        if intsx(1) > x0(1)
            xhigh = intsx(1) + 1;
        else
            xhigh = x0(1) + 1;
        end
        if intsy(1) > x0(2)
            yhigh = intsy(1) + 1;
        else
            yhigh = x0(2) + 1;
        end
        
        last = [intsx(1);intsy(1)];
        
        if checkunbounded(c,A,b,last,xinf,yinf) == 1
            points = Inf;
            return;
        end
        
        %plot x(0) point
        plot(x0(1),x0(2),'r.','Markersize',25);
        %set text to x(0) point
        hold on;
        quiver(x0(1),x0(2)/2,0,x0(2)/2,'black');
        text(x0(1),x0(2)/2,'x^{(0)}');
        
        %plot x(*)
        plot(intsy(1),intsy(1),'g.','Markersize',25);
        quiver(intsy(1),intsy(1)/2,0,intsy(1)/2,'black');
        text(intsy(1),intsy(1)/2,0,'x^{*}');
        axis([0 xhigh 0 yhigh]);
        hold off;
        return;
    end
    if length(intsx) == 2
        xhigh = max([intsx x0(1)]) + 1;
        yhigh = max([intsy x0(1)]) + 1;
 
        index = getclose(intsx,intsy,x0);
        points = [points [intsx(index);intsy(index)]];
        ct1 = c * [intsx(index);intsy(index)];
        bintsx = intsx;
        bintsy = intsy;
        intsx(index) = [];
        intsy(index) = [];
        ct2 = c * [intsx(1);intsy(1)];
        
        
        if ct1 > ct2
            points = [points [intsx(1);intsy(1)]];
            last = [points(1,3);points(2,3)];
            if checkunbounded(c,A,b,last,xinf,yinf) == 1
                points = Inf;
                return;
            end
            plot(bintsx,bintsy,'b-');
            hold on;
            %plot x(0) point
            plot(x0(1),x0(2),'r.','Markersize',25);
            %set text to x(0) point
            quiver(x0(1),x0(2)/2,0,x0(2)/2,'black');
            text(x0(1),x0(2)/2,'x^{(0)}');
            %plot x(t)
            plot(points(1,2),points(2,2),'b.','Markersize',25);
            quiver(points(1,2),points(2,2)/2,0,points(2,2)/2,'black');
            text(points(1,2),points(2,2)/2,0,'Simplex path');
            %plot x(*)
            plot(points(1,3),points(2,3),'g.','Markersize',25);
            quiver(points(1,3),points(2,3)/2,0,points(2,3)/2,'black');
            text(points(1,3),points(2,3)/2,0,'x^{*}');
            axis([0 xhigh 0 yhigh]);
        else
            last = [points(1,2);points(2,2)];
            if checkunbounded(c,A,b,last,xinf,yinf) == 1
                points = Inf;
                return;
            end
            plot(bintsx,bintsy,'b-');
            hold on;
            %plot x(0) point
            plot(x0(1),x0(2),'r.','Markersize',25);
            %set text to x(0) point
            quiver(x0(1),x0(2)/2,0,x0(2)/2,'black');
            text(x0(1),x0(2)/2,'x^{(0)}');
            %plot x(*)
            plot(points(1,2),points(2,2),'g.','Markersize',25);
            quiver(points(1,2),points(2,2)/2,0,points(2,2)/2,'black');
            text(points(1,2),points(2,2)/2,0,'x^{*}');
            axis([0 xhigh 0 yhigh]);
        end
        hold off;
        return;
    end

    
    %if all in a line
    slope = (intsy(2)-intsy(1))/(intsx(2) - intsx(1));
    flag = 0;
    for i = 3: length(intsx)
        if (intsy(i)-intsy(1))/(intsx(i) - intsx(1)) ~= slope
            flag =1;
            break;
        end
    end
    
    
    
    if flag == 1
        h0 = convhull(intsx,intsy);
    else
        h0 =[1:length(intsx)];
    end

h = h0;
h(length(h)) = [];
x = [];
y = [];
for i = 1:length(h)
    x = [x intsx(h(i))];
    y = [y intsy(h(i))];
end
pppp = [x;y];

index = getclose(x,y,x0);

indexs = [];
indexs = [indexs;index];
ii = 1;
while 1
    index = getbest(x,y,indexs(ii),c);
    if   c * [x(indexs(ii));y(indexs(ii))] - c * [x(index);y(index)] <= wucha
        break;
    else
        indexs = [indexs;index];
        ii = ii + 1;
    end
end


for i = 1:length(indexs)
    pointx = [pointx; x(indexs(i))];
    pointy = [pointy; y(indexs(i))];
    points = [points [x(indexs(i));y(indexs(i))]];
end







%get the firgure axis
xhigh = max(x) + 1;
xlow = min(x) - 1;
yhigh = max(y) + 1;
ylow = min(y) -1;
%tx = [x0(1)/xhigh x0(1)/xhigh];
%ty = [x0(2)/yhigh x0(2)/yhigh];
tx = [.5 0];
ty = [.5 0];
%unbounded
last = [pointx(length(pointx)); pointy(length(pointy))];
if checkunbounded(c,A,b,last,xinf,yinf) == 1
    points = Inf;
    return;
end

plot(intsx(h0),intsy(h0),'b-');
hold on;

%plot x(0) point
plot(x0(1),x0(2),'r.','Markersize',25);
%set text to x(0) point
quiver(x0(1),x0(2)/2,0,x0(2)/2,'black');
text(x0(1),x0(2)/2,'x^{(0)}');

%get x(0) out of points
pointx(1) = [];
pointy(1) = [];
%plot x(1) to x(* - 1)
plot(pointx,pointy,'b.','Markersize',25);
cx = mean(x);
cy = mean(y);
%plot arrow
for i = 1: length(pointx)-1
    quiver(cx,cy,pointx(i) - cx,pointy(i) - cy,'black');
end
%put text
text(cx,cy,'Simplex path');

%plot x(*)
plot(pointx(length(pointx)),pointy(length(pointx)),'g.','Markersize',25);
quiver(pointx(length(pointx)),(xhigh + pointy(length(pointx)))/2,0,-(xhigh - pointy(length(pointx)))/2,'black');
text(pointx(length(pointx)),(xhigh + pointy(length(pointx)))/2,0,'x^{*}');
hold off;
axis([0 xhigh 0 yhigh]);


    function index = getbest(x,y,index,c)
        if index == 1
            left = length(x);
        else
            left = index - 1;
        end
        if index == length(x)
            right = 1;
        else
            right = index + 1;
        end
        if c * [x(left);y(left)] - c * [x(right);y(right)] <= wucha
            index = left;
        else
            index = right;
        end 
    end
end
    function index = getclose(x,y,x0)
        distances = [];
        for i = 1: length(x)
            distance = sqrt((x(i) - x0(1)) * (x(i) - x0(1)) + (y(i) - x0(2)) * (y(i) - x0(2)));
            distances = [distances;distance];
        end
        [val,index] = min(distances);
    end

    function flag = checkunbounded(c,A,b,last,xinf,yinf)
        wucha = 1e-10;
        flag = 0;
        index = [];
        [M,d] = size(A);
        for i = 1:M
            if abs([A(i,1),A(i,2)] * last - b(i)) <= wucha
                index = [index;i];
            end
        end
        for i = 1:length(index)
            % if slope ~= 0
            if A(index(i),1) ~=0
                subA = [A(index(i),1) A(index(i),2);0 1];
                subb = [b(index(i));yinf];
                point = subA\subb;
                if A * point - b <= wucha
                    % c*point < c*last
                    if c * point -  c * last < 0
                        flag = 1;
                        return;
                    end
                end

            end
            %if slope ~= inf
            if A(index(i),2) ~=0
                subA = [A(index(i),1) A(index(i),2);1 0];
                subb = [b(index(i));xinf];
                point = subA\subb;
                if A * point - b <= wucha
                    if c * point -  c * last <= wucha
                        flag = 1;
                        return;
                    end
                end
            end
        end
    end

