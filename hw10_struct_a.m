clear all;
%Read in data
nodesRead = readmatrix('nodesR.txt','FileType','text'); %read in nodes
nodesSize = size(nodesRead); %number of nodes
elementsRead = readmatrix('elementsR.txt','FileType','text'); %read in elements
fid = fopen('elementsR.txt');
varEleLineOne = strsplit(fgetl(fid), ' ');
fclose(fid);
eleSize = size(elementsRead); %number of elements
forcesRead = readmatrix('forcesR.txt','FileType','text'); %read in forces
forceSize = size(forcesRead); %number of forces
displacementsRead = readmatrix('displacementsR.txt','FileType','text'); %read in displacements
dispSize = size(displacementsRead); %number of displacements

ndim = 2; %number of dimesions
nnds = 3; %number of nodes per element
nnodes = nodesSize(1); %number of nodes
nfbcs = forceSize(1); %number of forces
ndbcs = dispSize(1); %number of displacements

gcon = zeros(nodesSize(1),ndim); %initialize gcon

%populate gcon
for i = 1:nnodes
    for j = 1:ndim
        gcon(i,j) = 2*(i-1)+j;
    end
end
nodes = nodesRead(:,2:3);

elements = elementsRead(:,2:4); %elements matrix
Emat = ones(eleSize(1),1);
Emat = Emat.*str2double(varEleLineOne(3)); %Young's modulus
poisson = str2double(varEleLineOne(4));
fbcnd = forcesRead(:,1); %force number
fbcdof = forcesRead(:,2); %force dof
fbcval = forcesRead(:,3); %force value

dbcnd = displacementsRead(:,1); %displacement number
dbcdof = displacementsRead(:,2); %displacement dof
dbcval = displacementsRead(:,3); %displacement value
ndof = ndim*nodesSize(1); %number of dof

%initialize dof matrix and calculate gcon matrix
dofnum = zeros(ndbcs,nnodes);
for i = 1:ndbcs
    dofnum = gcon(dbcnd(i),dbcdof(i));
    for j = 1:nnodes
       for k = 1:ndim
          if(gcon(j,k) > dofnum)
              gcon(j,k) = gcon(j,k)-1;
          end
       end
    end
    gcon(dbcnd(i),dbcdof(i)) = nnodes*ndim;
    ndof = ndof-1;
end

%initialize global force matrix
Fglobal = zeros(nnodes*ndim, 1);
for i = 1:forceSize(1)
    nd = fbcnd(i);
    dof = fbcdof(i);
    gdof = gcon(nd,dof);
    Fglobal(gdof) = Fglobal(gdof) + fbcval(i);
end

%initialize global u matrix
uglobal = ones(nnodes, ndim);
for i = 1:dispSize(1)
    uglobal(dbcnd(i),dbcdof(i)) = dbcval(i);
end

%find reduced matrix size
usize = size(uglobal);
reducednum = nnodes*ndim-ndbcs;
%make reduced force matrix
Freduced = Fglobal(1:reducednum,:);

%initialize reduced K matrix
Kreduced = zeros(reducednum, reducednum);
Ktest = zeros(reducednum, reducednum);
count = 0;
%populate initial K matrix
 for iele = 1:eleSize(1)
     node1 = elements(iele,1);
     node2 = elements(iele,2);
     node3 = elements(iele,3);
     x1 = nodes(node1,1);
     y1 = nodes(node1,2);
     x2 = nodes(node2,1);
     y2 = nodes(node2,2);
     x3 = nodes(node3,1);
     y3 = nodes(node3,2);
     
     Adelta = 0.5*((x2*y3)-(x3*y2)+(x3*y1)-(x1*y3)+(x1*y2)-(x2*y1));
     b1 = (1/(2*Adelta))*(y2-y3);
     b2 = (1/(2*Adelta))*(y3-y1);
     b3 = (1/(2*Adelta))*(y1-y2);
     c1 = (1/(2*Adelta))*(x3-x2);
     c2 = (1/(2*Adelta))*(x1-x3);
     c3 = (1/(2*Adelta))*(x2-x1);
     B = [b1 0 b2 0 b3 0;
          0 c1 0 c2 0 c3;
          c1 b1 c2 b2 c3 b3];
     E = Emat(iele);
     v = poisson;
     Cstress = [E/(1-v^2) v*E/(1-v^2) 0;
                v*E/(1-v^2) E/(1-v^2) 0;
                0 0 E/(2*(1+v))];
     Kele = Adelta*transpose(B)*Cstress*B;
     for lnode1 = 1:nnds
         for lxy1 = 1:ndim
             ldof1 = ndim*(lnode1-1) + lxy1;
             gnode1 = elements(iele,lnode1);
             gdof1 = gcon(gnode1, lxy1);
             if(gdof1 > ndof)
                 continue
             end
             for lnode2 = 1:nnds
                for lxy2 = 1:ndim
                   ldof2 = ndim*(lnode2-1) + lxy2;
                   gnode2 = elements(iele, lnode2);
                   gdof2 = gcon(gnode2, lxy2);
                   if (gdof2 > ndof)
                       Freduced(gdof1) = Freduced(gdof1) - Kele(ldof1, ldof2)*uglobal(gnode2, lxy2);
                   else
                       Kreduced(gdof1,gdof2) = Kreduced(gdof1,gdof2) + Kele(ldof1, ldof2);
                   end
                end
             end
         end
     end
 end

%solve for displacements
ureduced = Kreduced\Freduced;
%populate uglobal
%r = size(ureduced);
for i= 1:nnodes
    for j= 1:ndim
        dof = gcon(i,j);
        if(dof<ndof + 1)
            uglobal(i,j) = ureduced(dof);
        end
    end
end

disp('Displacements = ');
disp(uglobal);

% figure

 s = 0.05; %scaling factor needs to be changed
 
strain = zeros(eleSize(1),3);
stress = zeros(eleSize(1),3);
 for i = 1:eleSize(1)
 
     node1 = nodes(elements(i, 1), :);
     x1 = node1(1);
     y1 = node1(2);
     u1 = uglobal(elements(i,1),1);
     v1 = uglobal(elements(i,1),2);
 
     node2 = nodes(elements(i, 2), :);
     x2 = node2(1);
     y2 = node2(2);
     u2 = uglobal(elements(i,2),1);
     v2 = uglobal(elements(i,2),2);
 
     node3 = nodes(elements(i, 3), :); %check the node 3 addition
     x3 = node3(1);
     y3 = node3(2);
     u3 = uglobal(elements(i,3),1);
     v3 = uglobal(elements(i,3),2);
 
     Adelta = 0.5*((x2*y3)-(x3*y2)+(x3*y1)-(x1*y3)+(x1*y2)-(x2*y1));
     b1 = (1/(2*Adelta))*(y2-y3);
     b2 = (1/(2*Adelta))*(y3-y1);
     b3 = (1/(2*Adelta))*(y1-y2);
     c1 = (1/(2*Adelta))*(x3-x2);
     c2 = (1/(2*Adelta))*(x1-x3);
     c3 = (1/(2*Adelta))*(x2-x1);
     
     B = [b1 0 b2 0 b3 0;
          0 c1 0 c2 0 c3;
          c1 b1 c2 b2 c3 b3];
     u = [u1;
          v1;
          u2;
          v2;
          u3;
          v3];
     strainEle = B*u;
     for a = 1:3
         strain(i,a) = strainEle(a);
     end
     E = Emat(iele);
     v = poisson;
     Cstress = [E/(1-v^2) v*E/(1-v^2) 0;
                v*E/(1-v^2) E/(1-v^2) 0;
                0 0 E/(2*(1+v))];
     stressEle = Cstress*strainEle;
     for b = 1:3
         stress(i,b) = stressEle(b);
     end
 end

%middle of triangle
edgesx = zeros(eleSize(1),1);
edgesy = zeros(eleSize(1),1);
numEdgesx = 0;
numEdgesy = 0;

% figure
for i = 1:eleSize(1)
     node1 = elements(i,1);
     node2 = elements(i,2);
     node3 = elements(i,3);
     x1 = nodes(node1,1);
     y1 = nodes(node1,2);
     x2 = nodes(node2,1);
     y2 = nodes(node2,2);
     x3 = nodes(node3,1);
     y3 = nodes(node3,2);
     if((x1 == 0 && x2 == 0) || (x1 ==0 && x3 == 0) || (x2 == 0 && x3 == 0))
        edgesy(i) = i; 
        numEdgesy = numEdgesy + 1;
     end
     if((y1 == 0 && y2 == 0) || (y1 ==0 && y3 == 0) || (y2 == 0 && y3 == 0))
        edgesx(i) = i; 
        numEdgesx = numEdgesx + 1;
     end
  
end
triangle = zeros(eleSize(1),1);
triangleRedx24 = zeros(numEdgesx,1);
triangleRedy24 = zeros(numEdgesy,1);
tri2y24 = zeros(numEdgesx,1);
tri2x24 = zeros(numEdgesy,1);
strainRedx = zeros(numEdgesx,1);
strainRedy = zeros(numEdgesy,1);
stressRedx24 = zeros(numEdgesx,1);
stressRedy24 = zeros(numEdgesy,1);
redValx = 1;
redValy = 1;
for i = 1:eleSize(1)
     node1 = elements(i,1);
     node2 = elements(i,2);
     node3 = elements(i,3);
     x1 = nodes(node1,1);
     y1 = nodes(node1,2);
     x2 = nodes(node2,1);
     y2 = nodes(node2,2);
     x3 = nodes(node3,1);
     y3 = nodes(node3,2);
     triangle(i,1) = (x1+x2+x3)/3;
     triangle(i,2) = (y1+y2+y3)/3;
     if(edgesx(i)>0)
         triangleRedx24(redValx) = (x1+x2+x3)/3;
         strainRedx(redValx) = strain(i,1);
         stressRedx24(redValx) = stress(i,1);
         tri2y24(redValx) = stress(i,2);
         redValx = redValx + 1;
     end
     if(edgesy(i)>0)
         triangleRedy24(redValy) = (y1+y2+y3)/3;
         strainRedy(redValy) = strain(i,2);
         stressRedy24(redValy) = stress(i,2);
         tri2x24(redValy) = stress(i,1);
         redValy = redValy + 1;
     end
end

disp('Stress = ');
disp(stress);
disp('Strain = ');
disp(strain);
