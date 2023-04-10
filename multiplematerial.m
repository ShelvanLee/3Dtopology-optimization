function multiplematerial()

% 3D Topology Optimization with stress constraints

clear all;
close all;
clc;

addpath 'F:\mosek\9.3\toolbox\R2015a'  %Path for MOSEK
addpath 'F:\mosek\9.3\toolbox\R2015a'  %Path for MOSEK

%%                       Geometry of domain

M_size = 0.2;

ax = [0:M_size:10];
ay = [-5:M_size:0];
az = [0:0.01:0.01];
% az = [0:M_size:3];

M_size_target = 0.04;
Refine_Target={['Edge'], [0.15], [0.85]}; % type - Volume/Edge; size; density_threhold.
% Refine_Target={['Volume'], [0.02]};

filter_radius = M_size*1.5;

traction_b = 0.01*382.2/2; % 200N
% traction_b = 3*382.2/2; % 200N

[x,y,z] = meshgrid(ax,ay,az); % a cube
x = [x(:)];
y = [y(:)];
z = [z(:)];

%%           Generate mesh and store node and element information
DT = delaunayTriangulation(x,y,z);

%{
% draw surface of the body
[F,P] = freeBoundary(DT);
trisurf(F,P(:,1),P(:,2),P(:,3),'FaceColor','cyan','FaceAlpha',0.8);
%}

node_xyz = DT.Points';
element = DT.ConnectivityList';

Tnode = length(node_xyz);
Telement = length(element);

node_xyz = [[1:Tnode];node_xyz];
element  = [[1:Telement];element];

fprintf(1,'Total number of node - Tnode =  %d  \n', Tnode);
fprintf(1,'Total number of element - Telement =  %d \n', Telement);
    
%%                  Initial stress and displacement                 %
Sigm = zeros(Tnode*6,1); %[sig_x sig_y sig_z sig_xy sig_yz sig_xz];

Penalty_w  = zeros(Tnode,1);
Topo_w     = ones(Tnode,1);


%%                  Material properties

% density
pho = 0;

% Plastic properties
phy = 0.0/180*3.1415926;
cohesion = 300; 
%cohesion = cohesion/1.732;


%%               Control parameters

du=0.005; % in-effect
istep=1;
nstep=1;

iteration_num = 10;
iter = 1;
while (iter <= iteration_num)

%%             Loop over
while (istep<=nstep)
    
    if (istep==1)
        
        xtan =tan(phy);
    %    A = xtan/sqrt(9+12*xtan^2);  
    %    B = 3*cohesion/sqrt(9+12*xtan^2);
        
        A=2*sin(phy)/sqrt(3)/(3+sin(phy));

        B=6*cohesion*cos(phy)/sqrt(3)/(3+sin(phy));
     
        B = B/2;
        
        
        D_g=CalD_d_3DNSFEM(Tnode);
        
        % Determine elements adjacent to each node
        Adjacent_elements=vertexAttachments(DT);
        
        % Calculate volume of each element
        Velement=Cal_Volumeelement(element,node_xyz,Telement);
        Vcell=Cal_Volumecell(Tnode,Velement,Adjacent_elements);
        
        fprintf(1, '\n Filter_info = CalPreFilter started \n');
        
        Filter_info = CalPreFilter(Tnode, node_xyz, element, Adjacent_elements, filter_radius);    % information for Filtering (struct:)
                
        
        % Calculate BT for NSFEM
        
        fprintf(1, '\n BT_node_3DNSFEM = CalBT_node_3DNSFEM \n');
        
        BT_node_3DNSFEM = CalBT_node_3DNSFEM(Tnode,Telement,node_xyz,element,Velement,Vcell,Adjacent_elements);
        
        %% Displacement boundary conditions
        [Au, Omeg] = calAu_Omeg(Tnode,node_xyz,du);
        
    end
    
    
     %% Calculate the penalty for the density and LT_g
    Penalty_w(:) = exp(10*(1-Topo_w(:)));
    
    LT_g = Penalty_w.*Vcell;
    
    %% Calculate A_g_xi - a factor relating xi and theta (density)
    A_g_xi = sparse(1:Tnode,1:Tnode,1);
    
    for i=1:Tnode
       
        index = 6*(i-1)+1;
        
        A_g_xi(i,i) = -A*sum(Sigm(index:index+2)) + B ;
        
    end
    
    
    %% debug from here
    %   [Sigma  pho  theta(material)  xi ru]
    
    
    Aeq1 = [BT_node_3DNSFEM  spalloc(3*Tnode,6*Tnode,0),  spalloc(3*Tnode,Tnode,0), ...
        spalloc(3*Tnode,Tnode,0),  -Au'];
    
    Aeq2 = [ -D_g            sparse(1:6*Tnode,1:6*Tnode,1)   spalloc(6*Tnode,Tnode,0)     ...
        spalloc(6*Tnode,Tnode,0)   spalloc(6*Tnode,3*Tnode,0)];
    

    
    Aeq3 = [ spalloc(Tnode,6*Tnode,0)  spalloc(Tnode,6*Tnode,0)  A_g_xi ...
          -sparse(1:Tnode,1:Tnode,1)     spalloc(Tnode,3*Tnode,0)];
    
    Aeq=[Aeq1; Aeq2; Aeq3];

    
    Ftraction = zeros(Tnode*3,1);
    
    for inode=1:Tnode
        if(node_xyz(2,inode)<=1.25 && abs(node_xyz(3,inode)-0)<0.001)
            
                Ftraction(3*(inode-1)+2)= -1;
            
        end
    end
    

    Ft=Ftraction/sum(Ftraction)*traction_b;

    beq = [Ft; zeros(6*Tnode,1); zeros(Tnode,1)];
    
    

%% call Mosek for solution

prob = [];

%% [Sigm rho  xi h y ru]
c = [zeros(6*Tnode,1); zeros(6*Tnode,1); LT_g; zeros(Tnode,1); -Omeg];

prob.c = c;
prob.a = Aeq;
prob.blc = beq;
prob.buc = beq;

prob.blx = -inf*eye(6*Tnode*2+5*Tnode,1);
prob.bux = inf*eye(6*Tnode*2+5*Tnode,1);

for i=(6*Tnode*2+1):(6*Tnode*2+Tnode)
    prob.blx(i)=0;
    prob.bux(i)=1;
end
%%
for inode=1:Tnode
    
    if(node_xyz(2,inode)<=1.25 && abs(node_xyz(3,inode)-0)<0.001)
        
        tem = 6*Tnode*2+inode;
        prob.blx(tem)=1;
        prob.bux(tem)=1;
        
    end
    
end

prob.cones = cell(Tnode,1);

for i=1:Tnode
    prob.cones{i}.type = 'MSK_CT_QUAD';
    
    prob.cones{i}.sub(1) = 6*Tnode*2 + Tnode + i;
    
    prob.cones{i}.sub(2) = 6*Tnode + 6*(i-1) + 1;
    prob.cones{i}.sub(3) = 6*Tnode + 6*(i-1) + 2;
    prob.cones{i}.sub(4) = 6*Tnode + 6*(i-1) + 3;
    prob.cones{i}.sub(5) = 6*Tnode + 6*(i-1) + 4;
    prob.cones{i}.sub(6) = 6*Tnode + 6*(i-1) + 5;
    prob.cones{i}.sub(7) = 6*Tnode + 6*(i-1) + 6;
    
end


param = [];

%

param.MSK_DPAR_MIO_TOL_ABS_RELAX_INT =1e-1;
param.MSK_DPAR_MIO_TOL_ABS_GAP = 1e-1;
param.MSK_DPAR_MIO_TOL_REL_GAP = 1e-1;
param.MSK_DPAR_MIO_REL_GAP_CONST = 1e-1;
%}
param.MSK_IPAR_PRESOLVE_USE = 'MSK_PRESOLVE_MODE_OFF';

fprintf(1, '\n Mosek started \n');

tic;

[exitflag, res]=mosekopt('minimize echo(0)',prob,param);

toc;

%if (exitflag ~=0)
    fprintf(1, '\n %d\n', exitflag);
%end

sol = res.sol;

Sigm = res.sol.itr.xx(1:6*Tnode);
Un =  res.sol.itr.y(1:3*Tnode);
    
Topo_w = res.sol.itr.xx(2*6*Tnode+1:2*6*Tnode+Tnode);

Volume_ratio = sum(Topo_w'*Vcell)/sum(Vcell);
fprintf(1,'istep = %d;  Volume ratio is %f  \n',istep, Volume_ratio);
 %
Filter_flag = 1;

if(Filter_flag)
    
    Topo_w = CalFilter(Tnode, Vcell, Topo_w, Filter_info);

end
%

if(istep == nstep)
  %  Plot_pq(A,B,Tnode,Sigm);
end


%NODE_XYZ=datanode(:,2:end);
tem = node_xyz';
NODE_XYZ=tem(:,2:end);


%ELEMENT=element(:,2:end);

tem = element';
ELEMENT=tem(:,2:end);

if(istep == nstep)
 output_paraview(iter, Tnode, Telement, NODE_XYZ, ELEMENT, Sigm, Un, Topo_w);
end


istep=istep+1;


end


%% mesh adaptivity

istep=1;

[DT, Tnode, Telement, Topo_w, Sigm] = Mesh_Adaptivity (Tnode, Telement, node_xyz, element, Velement, Topo_w, Sigm, Refine_Target, DT);

Penalty_w  = zeros(Tnode,1);

node_xyz = DT.Points';
element = DT.ConnectivityList';

node_xyz = [[1:Tnode];node_xyz];
element  = [[1:Telement];element];

fprintf(1,'Total number of node - Tnode =  %d  \n', Tnode);
fprintf(1,'Total number of element - Telement =  %d \n', Telement);
               

iter = iter +1;

end



end


function Plot_pq(A,B,Tnode,Sigm)

xx = zeros(Tnode,1);
yy = zeros(Tnode,1);

for i=1:Tnode
   sig_node = Sigm(6*i-5:6*i);
   
   xx(i) = sig_node(1) + sig_node(2) + sig_node(3);
   
   yy(i) = sqrt( ((sig_node(1) - sig_node(2))^2 + (sig_node(2) - sig_node(3))^2 + (sig_node(1) - sig_node(3))^2)/6 ...
                 + (sig_node(4))^2 + (sig_node(5))^2 + (sig_node(6))^2);
    
end

x_min = min(xx);
y_min = min(yy);
x_max = max(xx);
y_max = max(yy);


plot([x_min, x_max],[B - A*x_min,  B - A*x_max]);

hold on; 
plot(xx,yy,'.');

hold on;

%clf;

end


function [D_g]=CalD_d_3DNSFEM(Tnode)

ConstitutiveModel = 'DP';

switch(ConstitutiveModel)
    
    case 'DP'
        % D_node (6x6 matrix for 3D)
        SQRT_SIX = sqrt(6);
        D = [1/SQRT_SIX  -1/SQRT_SIX             0    0   0   0;...
                      0   1/SQRT_SIX   -1/SQRT_SIX    0   0   0;...
             1/SQRT_SIX            0   -1/SQRT_SIX    0   0   0;...
                      0            0             0    1   0   0;...
                      0            0             0    0   1   0;...
                      0            0             0    0   0   1];
        % d_node = [B  0  0  0  0  0  0]';
        index_r = zeros(6*6*Tnode,1);
        index_c = zeros(6*6*Tnode,1);
        value = zeros(6*6*Tnode,1);
        
        index = 1;
        for inode=1:Tnode
            
            jrow = 6*(inode-1);
            jcol = 6*(inode-1);
            
            for k=1:6
                index_r(index:(index+5)) = jrow + k;
                index_c(index:(index+5)) = (jcol+1):(jcol+6);
                
                value(index:(index+5)) = D(k,:);
                index = index +6;
            end
            
        end
        D_g = sparse (index_r,index_c,value);
        
    otherwise
        fprintf(1, 'ConstitutiveModel %s is not available!!\n', ConstitutiveModel);
end

end

function [Velement]=Cal_Volumeelement(element,node_xyz,Telement)

%  Calculations of volume of tetrahedral elements  %

%  ----input data includes-------%

%    element(node_i node_j node_k node_l): using 4-node tetrahedral elements
%    Telement=length(element): the total number of elements in the domain
%    gcoord_node (xi yi zi)


%                      |1 xi yi zi|
% Velement_i=(1/6) det|1 xj yj zj|
%                      |1 xk yk zk|
%                      |1 xl yl zl|

Velement =zeros(Telement,1);

for iel=1:Telement
    index_node=element(2:end,iel)'; % node_i node_j node_k node_l
    V_iel1= [1, node_xyz(2:end,index_node(1))'];
    V_iel2= [1, node_xyz(2:end,index_node(2))'];
    V_iel3= [1, node_xyz(2:end,index_node(3))'];
    V_iel4= [1, node_xyz(2:end,index_node(4))'];
    V_iel = [V_iel1; V_iel2; V_iel3; V_iel4];
    Velement(iel) = det(V_iel)/6; %?
end

end


function [Vcell]=Cal_Volumecell(Tnode,Velement,Adjacent_elements)

Vcell = zeros(Tnode,1);

for i=1:Tnode
    
    Vcell(i) = 0.25*sum(Velement(Adjacent_elements{i}(:)));
    
end

end


function [Filter_info] = CalPreFilter(Tnode, node_xyz, element, Adjacent_elements, filter_radius)

%   Filter_info =
%
%       struct with fields:
%
%                   TN: [Tnode×1 double]
%                Point: {Tnode×1 cell}
%               Weight: {Tnode×1 cell}


Filter_info = struct;

Filter_info.TN = zeros(Tnode,1);

Filter_info.Point  = cell(Tnode,1);

Filter_info.Weight = cell(Tnode,1);

Max_neighbour_num = 100;


SearchCase = '02';

switch (SearchCase)
    
    case('01')
        Neighbour_num = zeros(Tnode,1);
        
        Filter_Point_index =  zeros(Tnode,Max_neighbour_num);
        
        Filter_Weight_index = zeros(Tnode,Max_neighbour_num);
        
        for i=1:Tnode
            
            Point_i_xyz = node_xyz(2:4,i);
            
            distance = zeros(Tnode,1);
            
            
            parfor j=1:Tnode
                
                Point_j_xyz = node_xyz(2:4,j);
                
                distance(j) = norm(Point_i_xyz - Point_j_xyz);
                
            end
            
            index = 0;
            for j=1:Tnode
                
                if distance(j)<=filter_radius
                    
                    index = index+1;
                    Filter_Point_index(i,index) = j;
                    Filter_Weight_index(i,index) = exp(-0.5*(distance(j)/filter_radius)^2);
                end
                
            end
            Neighbour_num(i) = index;
            
        end
        
        
        Filter_info.TN = Neighbour_num;
       
        
        for i=1:Tnode
            
            TN = Neighbour_num(i);
            
            Filter_info.Point{i}  = Filter_Point_index(i,1:TN);
            Filter_info.Weight{i} = Filter_Weight_index(i,1:TN);
            
        end
        
    case('02')
        
        for i = 1:Tnode
           
            Elem_Neighbour = cell2mat(Adjacent_elements(i));
            
            TE = length(Elem_Neighbour);
            
            PointList = zeros(4*TE,1);
            

            for j=1:TE
            
                PointList(4*j-3:4*j) = element(2:5,Elem_Neighbour(j));
            
            end

            PointList= unique(PointList);
            
            TN = length(PointList);
            
 
            Point_i_xyz = node_xyz(2:4,i);
            
            distance = zeros(TN,1);
            weight   = zeros(TN,1);
            
            for j=1:TN 
                
                Point_j_xyz = node_xyz(2:4,PointList(j));

                distance(j) = norm(Point_i_xyz - Point_j_xyz);
                
                weight(j) = exp(-0.5*(distance(j)/filter_radius)^2);

            end 
            
            Filter_info.Point{i} = PointList;
            
            Filter_info.TN(i) = TN;
            
            Filter_info.Weight{i} = weight;

        end
        
        
end


end

function [Topo_w] = CalFilter(Tnode, Vcell, Topo_w, Filter_info)

Topo_w_Filter = zeros(Tnode,1);

for i=1:Tnode
    
    TN = Filter_info.TN(i);
    
    point_list = Filter_info.Point{i};
    
    weight = Filter_info.Weight{i};
    
    Upper = 0;
    Lower = 0;
    
    
    for j = 1:TN
        
        inode = point_list(j);
        
        tem = weight(j)*Vcell(inode);
        
        Upper = Upper + tem*Topo_w(inode);
        
        Lower = Lower +tem;
        
    end
    
    Topo_w_Filter (i) =  Upper/Lower;
    
end


Topo_w = Topo_w_Filter;


end

function [BT_node_3DNSFEM]=CalBT_node_3DNSFEM(Tnode,Telement,node_xyz,element,Velement,Vcell,Adjacent_elements)
    

%BT_node_3DNSFEM

% 4-node terahedral elements to approximate displacement
% 1-node terahedral elements to approximate stress
% 4 node Gauss integration
%              | x1 y1 z1|
%              | x2 y2 z2|
% node_xyz  =  | x3 y3 z3|
%              | x4 y4 z4|
%====================================================%
BT_e=cell(Tnode,1);

BT_e = calBT3D_e_4node(node_xyz,element,Telement,Velement);

index = 1;

for iel=1:Telement
    
    BT_iel=0.25*BT_e{iel};
    
    for inode = 1:4
        
        node_index = element(inode+1,iel);
        icol = node_index*6-5;          % (stress)
        
        for i=1:4
            
            irow=element(i+1,iel)*3-2; % node 1 (displacement)
            
            index_r(index:(index+5)) = irow;
            index_c(index:(index+5)) = icol:icol+5;
            value(index:index+5) = BT_iel(3*i-2,:);
            index=index+6;
            
            index_r(index:(index+5)) = irow+1;
            index_c(index:(index+5)) = icol:icol+5;
            value(index:index+5) = BT_iel(3*i-1,:);
            index=index+6;
            
            index_r(index:(index+5)) = irow+2;
            index_c(index:(index+5)) = icol:icol+5;
            value(index:index+5) = BT_iel(3*i,:);
            index=index+6;
            
        end
        
    end
    
end

BT_node_3DNSFEM = sparse(index_r,index_c,value);

end

function [BT_e]=calBT3D_e_4node(node_xyz,element,Telement,Velement)

% 4-node terahedral elements to approximate displacement
% 1-node terahedral elements to approximate stress
% 4 node Gauss integration
%              | x1 y1 z1|
%              | x2 y2 z2|
% node_xyz  =  | x3 y3 z3|
%              | x4 y4 z4|
%====================================================%
BT_e=cell(Telement,1);


for iel=1:Telement
    
    node_element = element(2:end,iel)';
    node_xyz_e = zeros(4,3);
    node_xyz_e(1,:) = node_xyz(2:end,node_element(1))';
    node_xyz_e(2,:) = node_xyz(2:end,node_element(2))';
    node_xyz_e(3,:) = node_xyz(2:end,node_element(3))';
    node_xyz_e(4,:) = node_xyz(2:end,node_element(4))';

    
    dNdr = zeros(3,4);
    %{
            dNdr(1,:) = [1 0 0 -1];
            dNdr(2,:) = [0 1 0 -1]; % N1=xi; N2=eta; N3=zeta; N4=1-xi-eta-zeta
            dNdr(3,:) = [0 0 1 -1];
    %}
    
    dNdr(1,:) = [-1 1 0 0 ];
    dNdr(2,:) = [-1 0 1 0 ]; % N2=xi; N3=eta; N3=zeta; N1=1-xi-eta-zeta
    dNdr(3,:) = [-1 0 0 1 ];
    
    J = dNdr * node_xyz_e; % node_xy_e' [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4];
    
    dNdx = J\dNdr;
    
    
    Bu_e = zeros(6,12); % 4-node terahedral elements to approximate displacement ; each node has 3 degree of freedom
    
    for j=1:4
        Bu_e(1,3*j-2) = dNdx(1,j);
        Bu_e(2,3*j-1) = dNdx(2,j);
        Bu_e(3,3*j)   = dNdx(3,j);
        
        Bu_e(4,3*j-2) = dNdx(2,j); %
        Bu_e(4,3*j-1) = dNdx(1,j); %
        
        Bu_e(5,3*j-1) = dNdx(3,j);
        Bu_e(5,3*j)   = dNdx(2,j);
        
        Bu_e(6,3*j-2) = dNdx(3,j);
        Bu_e(6,3*j)   = dNdx(1,j);
    end
    
    % 4-node terahedral elements to approximate stress, each node has 6 stress
    % components
    
    
    BT_e{iel}=Bu_e'*Velement(iel);
    
    
end
end


function [Au, Omeg] = calAu_Omeg(Tnode,node_xyz,du)
%   ------------------   % %  zmin=-10; zmax=0; 
%  |                 |
%  |                 |
%  |                 |
%  |                 |
%   ------------------  %  xmin=ymin=0; xmax=ymax=10; 

% Bx=2;By=2;  

index = 1;
index_r = zeros(3*Tnode,1);
index_c = zeros(3*Tnode,1);
value = zeros(3*Tnode,1);

%  Au  = zeros(2*Tnode);
Omeg = zeros(3*Tnode,1);

%
max_xi =max(node_xyz(2,:)); % max_xi
min_xi =min(node_xyz(2,:)); % min_xi

max_yi =max(node_xyz(3,:)); % max_yi
min_yi =min(node_xyz(3,:)); % min_yi

max_zi =max(node_xyz(4,:)); % max_zi
min_zi =min(node_xyz(4,:)); % min_zi
%}

for i=1:Tnode
    
    index_r(3*i-2) = 3*i-2;
    index_c(3*i-2) = 3*i-2;
    index_r(3*i-1) = 3*i-1;
    index_c(3*i-1) = 3*i-1;
    index_r(3*i)   = 3*i;
    index_c(3*i)   = 3*i;
   
    xi=node_xyz(2,i); %  xmin=ymin=0; xmax=ymax=10; 
    yi=node_xyz(3,i);
    zi=node_xyz(4,i);
    % Plane-strain condition
    % uy=0 at y=ymin and y=ymax
    
    if (abs(xi-min_xi)<0.001)    % ux=uz=0 at x=xmin
        value(3*i-2)  = 1;  Omeg(3*i-2)=0;  
      %  value(3*i-1)  = 1;  Omeg(3*i-1)=0;
        value(3*i-0)  = 1;  Omeg(3*i-0)=0; 
    end
    
    if (abs(xi-max_xi)<0.001)  % ux=uy=uz=0 at x=xmax
        value(3*i-2)  = 1;  Omeg(3*i-2)=0;  
        value(3*i-1)  = 1;  Omeg(3*i-1)=0;
        value(3*i-0)  = 1;  Omeg(3*i-0)=0; 
    end

     
end

Au = sparse (index_r,index_c,value);

end

function [DT, Tnode, Telement, Topo_w, Sigm] = Mesh_Adaptivity (Tnode, Telement, node_xyz, element, Velement, Topo_w, Sigm, Refine_Target, DT)


MeshAdaptCase = Refine_Target{1}; 

V_target = Refine_Target{2};

Density_Threshold = Refine_Target{3};

switch MeshAdaptCase
    
    case 'Volume'
        
        node_xyz_add = zeros (4, Telement);
        
        Topo_w_add = zeros(Telement,1);
        
        Sigm_add = zeros(6*Telement,1);
        
        Tnode_add = 0;
        
        for ielement = 1:Telement
            
            node_ID = element(2:5, ielement);
            
            Top_sum = sum (Topo_w(node_ID));
            
            if (Top_sum > Density_Threshold && Velement(ielement) > V_target)
                
                Tnode_add = Tnode_add +1;
                
                Topo_w_add (Tnode_add) = 0.25*Top_sum ;
                
                Sigm_add (6*Tnode_add-5) = 0.25*(Sigm(6*node_ID(1)-5) + Sigm(6*node_ID(2)-5) + Sigm(6*node_ID(3)-5) + Sigm(6*node_ID(4)-5));
                Sigm_add (6*Tnode_add-4) = 0.25*(Sigm(6*node_ID(1)-4) + Sigm(6*node_ID(2)-4) + Sigm(6*node_ID(3)-4) + Sigm(6*node_ID(4)-4));
                Sigm_add (6*Tnode_add-3) = 0.25*(Sigm(6*node_ID(1)-3) + Sigm(6*node_ID(2)-3) + Sigm(6*node_ID(3)-3) + Sigm(6*node_ID(4)-3));
                Sigm_add (6*Tnode_add-2) = 0.25*(Sigm(6*node_ID(1)-2) + Sigm(6*node_ID(2)-2) + Sigm(6*node_ID(3)-2) + Sigm(6*node_ID(4)-2));
                Sigm_add (6*Tnode_add-1) = 0.25*(Sigm(6*node_ID(1)-1) + Sigm(6*node_ID(2)-1) + Sigm(6*node_ID(3)-1) + Sigm(6*node_ID(4)-1));
                Sigm_add (6*Tnode_add-0) = 0.25*(Sigm(6*node_ID(1)-0) + Sigm(6*node_ID(2)-0) + Sigm(6*node_ID(3)-0) + Sigm(6*node_ID(4)-0));
                
                node_xyz_add(2, Tnode_add) = sum(node_xyz(2,node_ID))*0.25;
                node_xyz_add(3, Tnode_add) = sum(node_xyz(3,node_ID))*0.25;
                node_xyz_add(4, Tnode_add) = sum(node_xyz(4,node_ID))*0.25;
                
            end
            
        end
        
    case 'Edge'
        
        Edge_index = edges(DT);
        
        Tedge = length (Edge_index);
        

        node_xyz_add = zeros (4, Tedge);
        
        Topo_w_add = zeros(Tedge,1);
        
        Sigm_add = zeros(6*Tedge,1);
        
        
        Tnode_add = 0;
        
        for iedge = 1:Tedge
            
            node_ID = Edge_index(iedge, 1:2);
            
            Edge_length = norm(node_xyz(2:4,node_ID(1)) - node_xyz(2:4,node_ID(2)));
            
            Top_sum = sum (Topo_w(node_ID));
            
            if (Top_sum > Density_Threshold && Edge_length > V_target)
                
                Tnode_add = Tnode_add +1;
                
                Topo_w_add (Tnode_add) = 0.5*Top_sum ;
                
                Sigm_add (6*Tnode_add-5) = 0.25*(Sigm(6*node_ID(1)-5) + Sigm(6*node_ID(2)-5) );
                Sigm_add (6*Tnode_add-4) = 0.25*(Sigm(6*node_ID(1)-4) + Sigm(6*node_ID(2)-4) );
                Sigm_add (6*Tnode_add-3) = 0.25*(Sigm(6*node_ID(1)-3) + Sigm(6*node_ID(2)-3) );
                Sigm_add (6*Tnode_add-2) = 0.25*(Sigm(6*node_ID(1)-2) + Sigm(6*node_ID(2)-2) );
                Sigm_add (6*Tnode_add-1) = 0.25*(Sigm(6*node_ID(1)-1) + Sigm(6*node_ID(2)-1) );
                Sigm_add (6*Tnode_add-0) = 0.25*(Sigm(6*node_ID(1)-0) + Sigm(6*node_ID(2)-0) );
                
                node_xyz_add(2, Tnode_add) = sum(node_xyz(2,node_ID))*0.5;
                node_xyz_add(3, Tnode_add) = sum(node_xyz(3,node_ID))*0.5;
                node_xyz_add(4, Tnode_add) = sum(node_xyz(4,node_ID))*0.5;
                
            end
            
        end
        
        

end

%% form new matrix for vectors

Topo_w_add = Topo_w_add(1:Tnode_add);
Sigm_add   = Sigm_add(1:6*Tnode_add);
node_xyz_add = node_xyz_add(:,1:Tnode_add);

node_xyz = [node_xyz, node_xyz_add];

Topo_w = [Topo_w ; Topo_w_add];

Sigm = [Sigm; Sigm_add];


%% remesh

x = [node_xyz(2,:)]';
y = [node_xyz(3,:)]';
z = [node_xyz(4,:)]';

%%           Generate mesh and store node and element information
DT = delaunayTriangulation(x,y,z);

%{
% draw surface of the body
[F,P] = freeBoundary(DT);
trisurf(F,P(:,1),P(:,2),P(:,3),'FaceColor','cyan','FaceAlpha',0.8);
%}

node_xyz = DT.Points';
element = DT.ConnectivityList';

Tnode = length(node_xyz);
Telement = length(element);


end

function output_paraview(istep, Tnode, TELEMENT, NODE_XYZ, ELEMENT, Sigm, Un, Topo_w)
% output post-processing file for Paraview

    Num_str = int2str(istep);
    filename_res = strcat ( 'ParaviewResult', Num_str, '.vtk' );
    
    %**********************************************************************
    % output mesh information
    NT998 = fopen( filename_res,'w');
    
    % output mesh information
    fprintf(NT998,'# vtk DataFile Version 4.1 \n');
    fprintf(NT998,'This file was created by matlab source code \n');
    fprintf(NT998,'ASCII\n');
    fprintf(NT998,'DATASET UNSTRUCTURED_GRID \n');
    
    fprintf(NT998,'POINTS %d  FLOAT \n', Tnode);
    
    fprintf(NT998,'%d     %f     %f\n', NODE_XYZ');
    
    fprintf(NT998,'\n');
    fprintf(NT998,'CELLS  %d  %d \n', TELEMENT, TELEMENT*5);

    ELEMENT(:,:)=ELEMENT(:,:)-1;

    for i=1:TELEMENT
    fprintf(NT998,'%d    %d     %d     %d     %d  \n',4, ELEMENT(i,1:4));
    end 
    
    
    fprintf(NT998,'CELL_TYPES  %d\n', TELEMENT);
    fprintf(NT998,'%d   \n',ones(TELEMENT,1)*10);
    

%**********************************************************************
    % output displacement (velocity) and stresses (tensor) for each node. 

    fprintf(NT998,'\n');
    fprintf(NT998,'POINT_DATA  %d\n', Tnode); 
    
    
    fprintf(NT998,'\n');
    fprintf(NT998,'VECTORS Displacement float \n');
    
    for i=1:Tnode
        fprintf(NT998,'%f  %f  %f\n', Un(3*(i-1)+1:3*(i-1)+3));
    end
    
    
    fprintf(NT998,'\n');
    fprintf(NT998,'VECTORS TopoData float \n');
    
    for i=1:Tnode
       tem = Sigm(6*(i-1)+1:6*(i-1)+6);
       Sig_VM(i) = sqrt((tem(1)-tem(2))^2 +(tem(2)-tem(3))^2 + (tem(1)-tem(3))^2 ...
                                           + tem(4)^2+ tem(5)^2+ tem(6)^2);
    end
    
    for i=1:Tnode
        fprintf(NT998,'%f  %f  %f\n', Topo_w(i), Sig_VM(i), 0);
    end

    fprintf(NT998,'\n');
    fprintf(NT998,'TENSORS Stress float \n');
    for i=1:Tnode
        index = 6*(i-1);
        Sig_tem = [Sigm(index+1)  Sigm(index+4)  Sigm(index+6) ...
                   Sigm(index+4)  Sigm(index+2)  Sigm(index+5) ...   
                   Sigm(index+6)  Sigm(index+5)  Sigm(index+3)];
        fprintf(NT998,'%f  %f  %f  \n', Sig_tem');
    end
    
    
    
    
    fprintf(NT998,'\n');
    fprintf(NT998,'CELL_DATA  %d\n', TELEMENT);
    fprintf(NT998,'SCALARS TopoW  float  1\n');
    fprintf(NT998,'LOOKUP_TABLE default\n'); 
    
    for i=1:TELEMENT
        
        node_ID = ELEMENT(i,1:4) + [1 1 1 1];
        
        
        fprintf(NT998,'%f \n', sum(Topo_w(node_ID))*0.25);
        
    end

    
    fclose (NT998);
    fprintf(1,'Output for paraview completed successfully - step  %d  \n', istep);
end
