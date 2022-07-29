% Generating Radom XY values Within range 0 to 100

rng default;
% m = no. of partions required
m = 500;
m = m + 0.15*m;
m = round(m);
%k = 0 + (100-0).*rand(m,1);
x = rand([m 1]);
y = rand([m 1]);
z = rand([m 1]);
%z = rand([m 1]);
x = x*100;
y = y*100;
z = z*100;
X = [x y];

%C = [x,y];
%C = meshgrid(1:100);
%s = pcolor(C);
%s.FaceColor = 'interp';
%hold on;
%colouring the grid


%dt = delaunayTriangulation(x(:),y(:));
%[V1,C1] = voronoiDiagram(dt);

% Generating Voronoi Diagram using randomly Genarated XY co-ordinates and
% retriving cell vertices and cell topology using Voronoin Function same
% can be achieved using delaunay triangulation


[V,C] = voronoin(X);
%voronoi(x,y)
p = [];
p1 = [];
C1 =C ;
% loop to find Cells No.s Containing inf points as vertices 
for i = 1:1:size(C,1)

   k = C{i};

   for j =1:1:size(k,2)

     k1 = k(1,j);

     if k1 == 1
         p(end+1) = i;
        
     end
   end
end

%loop to find cell no.s that have vertices outside boundaries (-ve XY and above 100)

for i=1:1:size(V,1)
    if V(i,1) > 100 || V(i,1) < 0 || V(i,2) > 100 || V(i,2) < 0

        p1(end + 1) = i;
        
    end
    
end

for i = 1:1:size(C,1)

   k = C{i};


   for j =1:1:size(k,2)

     k1 = k(j);
     for r=1:1:size(p1,2)
         if k1 == p1(1,r)
             %disp(k1)
             p(end+1) = i; 
         end
     end
   end
end

p = unique(p,'stable');
p = sort(p);
r = 0;
%Deleting Cells with Inf as veritices and vertices out of boundaries
for i = 1:1:size(C,1) - size(p,2)
   %disp(C(i,:))
   k = C{i};
   for j=1:1:size(p,2)
       if i == p(1,j) - r
          r =r + 1;
          C(i,:) = [];
          X(i,:) = [];
       end
   end

end

% Assigning XY co-ordinates to Cell No.s

for i=1:1:size(C,1)
    vx = [];
    k = C{i};
    for j=1:1:size(k,2)
        %disp(k(1,j))
        vx(j,:) =  V(k(1,j),:);
    end
    vx = {vx};
    Vc(i,:) = vx;
    
end

%loop Defining Properties of Cells and their Assignment

properties = [];
T = [];
S = [];
for i=1:1:size(C,1)
    k=C{i};
    Ti = 25;
    Si = 0 ;
    T(end+1) = Ti;
    S(end+1)= Si;
    phi1 = 0 + (360-0).*rand(size(C,1),1);
    phi2 = 0 + (180-0).*rand(size(C,1),1);
    phi3 = 0 + (360-0).*rand(size(C,1),1); 
    
    
    
    %z = zeros(length(N),1) + 25;
    
    %property_1 = (T S phi);
    
    %properties(end + 1) = property_1;
    
end 
T = T.';
S = S.';
phi = [phi1 phi2 phi3];
properties = [T S phi];


%Generating Mesh
%mesh_no = no of mesh elements required
mesh_no = 6400*4;
max_lim = 90;
min_lim = 10;

mesh_size = (max_lim - min_lim)/sqrt(mesh_no);
mesh_co = [];
start_co = mesh_size/2;
x_c = [];
y_c = [];
for i=1:1:sqrt(mesh_no)
    for j = 1:1:sqrt(mesh_no)
        xi = min_lim + start_co + (j-1)*mesh_size;
        yi = min_lim + start_co + (i-1)*mesh_size;
        x_c(end+1) = xi;
        y_c(end+1) = yi;
    end 
    
       
end
x_c = x_c.';
y_c = y_c.';

mesh_co = [x_c y_c];
%mesh_co = co-ordinates of mesh centres

for i=1:1:size(mesh_co,1)
    xt = mesh_co(i,1);
    yt = mesh_co(i,2);
    
    for j=1:1:size(Vc,1)
        vxy = Vc{j,1};
        vx  = vxy(:,1);
        vy = vxy(:,2);
        try
            ip = inpolygon(xt,yt,vx,vy);
            if ip ==1
                mesh_p(i,:) = [i,j];
            end
        end
    end
    
end


%mesh_p = Data of positon of mesh 
% ex 1 991 means mesh no.1 lies in cell no.991
mesh_properties = [];
for i =1:1:size(mesh_p,1)
    h = mesh_p(i,1);
    %disp(h)
    n = mesh_p(i,2);
    mesh_properties(i,:) = properties(n,:);
end

%pt1=[];
%pt2 =[];
%pt3 = [];
%pt4 = [];
mesh_vertices = [];
for i = 1:1:size(mesh_co,1)
    x1 = mesh_co(i,1) - mesh_size/2;
    x2 = mesh_co(i,1) + mesh_size/2;
    y1 = mesh_co(i,2) - mesh_size/2;
    y2 = mesh_co(i,2) + mesh_size/2;
        %mesh_v = {p}
    pt1{i} = [x1,y1];
    pt2{i} = [x1,y2];
    pt3{i} = [x2,y2];
    pt4{i} = [x2,y1];
    
end
mesh_vertices = [pt1 pt2 pt3 pt4];
mesh_vertices = mesh_vertices.';
mesh_vertices = cell2mat(mesh_vertices);
mesh_vertices = unique(mesh_vertices,'rows','stable');
%p = sort(p);

mesh_v = [pt1; pt2; pt3; pt4];
%mesh_v = mesh_v.';
%mesh_t = [];
%for i=1:1:size(mesh_v,1)
 %   for j=1:1:4
  %      z = mesh_v{i,j};
   %     for k = 1:1:size(mesh_vertices,1)
    %       if z == mesh_vertices(k,:)
     %          mesh_t(i,j) = k;
      %         
       %    end
        %end
    %end
    
%end

mesh_v=mesh_v.';
mesh_v=cell2mat(mesh_v);

ar1 = mesh_v(:,1);
ar2 = mesh_v(:,3);
ar3 = mesh_v(:,5);
ar4 = mesh_v(:,7);
ar = [ar1 ar2 ar3 ar4];
ar = ar.';

br1 = mesh_v(:,2);
br2 = mesh_v(:,4);
br3 = mesh_v(:,6);
br4 = mesh_v(:,8);
br = [br1 br2 br3 br4];
br = br.';

prompt = "Pick the property needed for colouring before the process:";
prompt = prompt + newline + "1.Orientation";
prompt = prompt + newline + "2.Temperature" + newline;
x = input(prompt);
        
if x == 1
patch(ar,br,mesh_properties(:,3),'Edgecolor','black');
elseif x == 2
patch(ar,br,mesh_properties(:,1),'Edgecolor','black');
end
%patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'FaceColor','flat');
%patch('Faces',mesh_v,'Vertices',mesh_vertices,'FaceVertexCData',mesh_properties(:,3),'FaceColor','flat','Edgecolor','black')



%Code to Plot Figure
A = C;
B = cellfun('length',A);
Cp = cellfun(@(v,n)[v,nan(1,max(B)-n)],A,num2cell(B),'UniformOutput',false);
Cp = vertcat(Cp{:});
patch('Faces',Cp,'Vertices',V,'FaceColor','none','Edgecolor','black')
%grid on
[x,y,z] = meshgrid(x,y,z);

%Code to name every cell in figure
%numv = size(X,1);
%vlabels = arrayfun(@(n) {sprintf('C%d', n)}, (2:numv)');
%hold on
%Hpl = text(X(2:end,1), X(2:end,2)+.2, vlabels, ...
  %    'FontWeight', 'normal', 'HorizontalAlignment',...
 %     'center', 'BackgroundColor', 'none');
%hold off

axis equal

disp("_____________________________________________________________________________________")
disp("Code Executed SuccessFully")
disp("_____________________________________________________________________________________")


