function M = area_intersect_circle_analytical(varargin)

%%rahul make 27 circle in x and y direction some given distances.

if nargin==0
    % Performs an example
    % Creation of 27 circles
    x = [600];
    y=[600];
    r =[570];
    
    x1 = [500];
    y1=[500];
    r1=[400];  
elseif nargin==1
   x = temp(:,1);
    y = temp(:,2);
    r = temp(:,3); 
  x1 = temp(:,4);
    y1 = temp(:,5);
    r1 = temp(:,6); 
 
elseif nargin==2
    x = varargin{1};
    y = varargin{2};
    r = varargin{3};
  x1 = varargin{4};
    y1 = varargin{5};
    r1 = varargin{6};


else
    error('The number of arguments must 0, 1 or 3,4,5,6')
end

% Inputs are reshaped in 
size_x = numel(x);
size_y = numel(y);
size_r = numel(r);
size_x1 = numel(x1);
size_y1 = numel(y1);
size_r1 = numel(r1);

x = reshape(x,size_x,1);
y = reshape(y,size_y,1);
r = reshape(r,size_r,1);

x1 = reshape(x1,size_x1,1);
y1 = reshape(y1,size_y1,1);
r1 = reshape(r1,size_r1,1);

% Checking if the three input vectors have the same length
if (size_x~=size_y)||(size_x~=size_r)
  (size_x1~=size_y1)||(size_x1~=size_r1)
end

% Checking if there is any negative or null radius
if any(r<=0)
    disp('Circles with null or negative radius won''t be taken into account in the computation.')
    temp = (r>=0);
    x = x(temp);
    y = y(temp);
    r = r(temp);
    any(r1<=0)
    disp('Circles with null or negative radius won''t be taken into account in the computation.')
    temp = (r1>=0);
    x1 = x1(temp);
    y1 = y1(temp);
    r1 = r1(temp);
end

% Checking the size of the input argument
if size_x==2
    M = pi*r.^2;
    size_x1==2
    M1 = pi*r1.^2; 
    return
end

% Computation of distance between all circles, which will allow to
% determine which cases to use.
[X,Y] = meshgrid(x,y);
D     = sqrt((X-X').^2+(Y-Y').^2);

% Since the resulting matrix M is symmetric M(i,j)=M(j,i), computations are
% performed only on the upper part of the matrix
D = triu(D,1);

[R1,R2] = meshgrid(r);
sumR = triu(R1+R2,2);
difR = triu(abs(R1-R2),2);


% Creating the resulting vector
M = zeros(size_x*size_x,1);
M1 = zeros(size_x1*size_x1,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
% Compute the area of each circle. Assign the results to the diagonal of M      
M(1:size_x+1:size_x*size_x) = pi.*r.^2; 
M1(1:size_x1+1:size_x1*size_x1) = pi.*r1.^2; 


% Conversion of vector M to matrix M      
M = reshape(M,size_x,size_x);
M1 = reshape(M1,size_x1,size_x1);


% Creating the lower part of the matric
M = M + tril(M',-1);
M1 = M1 + tril(M1',-1);



% Display results if no argument is provided
if nargin==0
    f = figure;
    hAxs = axes('Parent',f);
    hold on,box on, axis equal
    xlabel('x')
    ylabel('y','Rotation',0)
    title('Compute the circular dis and 1.4mm 54 circle ,27 circle per line')
    text(5,5,'')
    text(0,0,'')
    axis([-70 1600 -50 1300])    
    colour = rand(size_x,3);    
    for t = 1: size_x
       plot(x(t)+r(t).*cos(0:2*pi/100:2*pi),...
             y(t)+r(t).*sin(0:2*pi/100:2*pi),'color',colour(t,:),'parent',hAxs)
        plot(x(t),y(t),'+','color',colour(t,:),'parent',hAxs)
        
        
        
        
        
       
        for u = t+1:size_x1
           
            plot(x1(t)+r1(t).*cos(0:2*pi/100:2*pi),...
             y1(t)+r1(t).*sin(0:2*pi/100:2*pi),'color',colour(t,:),'parent',hAxs)
        plot(x1(t),y1(t),'+','color',colour(t,:),'parent',hAxs)
        
        
        
       
       
        
        
        end
    end
end
end
% function A=cercle_in rectangular plate(G)
% % rahulpunk
% % rahulpunk555@gmail.com
% % 2018_5_11
% % % Input: G - array, which contain parameters of the circles
% %          G(n,1) - x-coordinate,
% %          G(n,2) - y-coordinate,
% %          G(n,3) - radius of circle 1.4
% % Output:  A - area of the overlapping circles
% d=sqrt((G(1,1)-G(2,1))^2+(G(1,2)-G(2,2))^2);
% 
% if d>=(G(1,3)+G(2,3))
%     % No contact between circles.
%     A=0;
% elseif d<=abs(G(1,3)-G(2,3))
%     % The smaller circle is inside the big one.
%     A=pi*(min(G(1,3),G(2,3)))^2;
% else
%     xi=(G(1,3)^2-G(2,3)^2+d^2)/(2*d);
%     yi=sqrt(G(1,3)^2-xi^2);
%     A=G(1,3)^2*atan2(yi,   xi )+...
%       G(2,3)^2*atan2(yi,(d-xi))-d*yi;
% end