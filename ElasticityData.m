%Geometry
%%Square
ng = 90; %number of geometries
name = "square";
L = 0.3; %0.35; %0.35; %0.3; %0.5;
x1 = -L; y1 = -L;
x2 = L; y2 = -L;
x3 = L; y3 = L;
x4 = -L; y4 = L;
min_points = 10000; %Gross Points
L_out = 1.0; %1;

for i=1:ng
    theta = i*pi/180.0;
    x1r = x1*cos(theta) - y1*sin(theta);
    y1r = x1*sin(theta) + y1*cos(theta);    
    x2r = x2*cos(theta) - y2*sin(theta);
    y2r = x2*sin(theta) + y2*cos(theta);
    x3r = x3*cos(theta) - y3*sin(theta);
    y3r = x3*sin(theta) + y3*cos(theta);
    x4r = x4*cos(theta) - y4*sin(theta);
    y4r = x4*sin(theta) + y4*cos(theta);
    r1 = [3 4 -L_out L_out L_out -L_out -L_out -L_out L_out L_out];
    r2 = [3 4 x1r x2r x3r x4r y1r y2r y3r y4r];
    gdm = [r1; r2]';
    g = decsg(gdm,'R1-R2',['R1'; 'R2']');
    
    %thermal model
    thermalmodel = createpde('thermal','steadystate');
    gm = geometryFromEdges(thermalmodel,g);
     pdegplot(thermalmodel,'EdgeLabels','on');
    thermalBC(thermalmodel,'Edge',1,'Temperature',1);
    thermalBC(thermalmodel,'Edge',2,'Temperature',1);
    thermalBC(thermalmodel,'Edge',3,'Temperature',0);
    thermalBC(thermalmodel,'Edge',4,'Temperature',0);
    thermalBC(thermalmodel,'Edge',5,'Temperature',0);
    thermalBC(thermalmodel,'Edge',6,'Temperature',0);
    thermalBC(thermalmodel,'Edge',7,'Temperature',1);
    thermalBC(thermalmodel,'Edge',8,'Temperature',1);                   
    thermalProperties(thermalmodel,"ThermalConductivity",1.0);
    mesh = generateMesh(thermalmodel);
    %mesh = generateMesh(thermalmodel,'Hmax',0.105);
    thermalresults = solve(thermalmodel);
    T = thermalresults.Temperature;
    gradTx =thermalresults.XGradients;
    gradTy =thermalresults.YGradients;
    %Linear elasticity model (Plane Stress)
    structuralmodel = createpde('structural','static-planestress');
    structuralmodel.Geometry = gm;
    structuralmodel.Mesh = mesh;
    structuralProperties(structuralmodel,'YoungsModulus',1, ...
                                     'PoissonsRatio',0.3, ...'
                                     'CTE',1.0);
    
    %%%Outter cylinder
    structuralBC(structuralmodel,'Edge',1,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',2,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',7,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',8,'Constraint','fixed');
    %%%Inner cylinder
    structuralBC(structuralmodel,'Edge',3,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',4,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',5,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',6,'Constraint','fixed');
    
    structuralBodyLoad(structuralmodel,'Temperature',thermalresults);
    structuralmodel.ReferenceTemperature = 0; %No affect on the solution!
    thermalstressresults = solve(structuralmodel);
    u = thermalstressresults.Displacement.x ;
    v = thermalstressresults.Displacement.y ;
    stress = thermalstressresults.VonMisesStress ;
    coord = thermalresults.Mesh.Nodes';
    x_coord = coord(:,1);
    y_coord = coord(:,2);
    %disp(size(x_coord));
    if (size(x_coord)<min_points)
        aaa = size(x_coord);
        min_points = aaa(1);
    end
    %figure
    %pdeplot(structuralmodel,"XYData",v,"Contour","on","ColorMap","hot");
    %title 'v displacement';

    %writing the results
    id = fopen(name+'T'+string(i)+'.txt','w');
    fprintf(id,'%f\n',T);
    fclose(id);

    id = fopen(name+'dTdx'+string(i)+'.txt','w');
    fprintf(id,'%f\n',gradTx);
    fclose(id);

    id = fopen(name+'dTdy'+string(i)+'.txt','w');
    fprintf(id,'%f\n',gradTy);
    fclose(id);
    
    id = fopen(name+'u'+string(i)+'.txt','w');
    fprintf(id,'%f\n',u);
    fclose(id);
    
    id = fopen(name+'v'+string(i)+'.txt','w');
    fprintf(id,'%f\n',v);
    fclose(id);
    
    id = fopen(name+'x'+string(i)+'.txt','w');
    fprintf(id,'%f\n',x_coord);
    fclose(id);
    
    id = fopen(name+'y'+string(i)+'.txt','w');
    fprintf(id,'%f\n',y_coord);
    fclose(id);
end

display(min_points);

%clc;
%clear;

%%Triangle
ng = 120; %number of geometries
name = "triangle";
L = 0.3; %0.35; %0.35; %0.3; %0.5;
x1 = -L; y1 = -L;
x2 = L; y2 = -L;
x3 = 0; y3 = L;
min_points = 10000; %Gross Points
L_out = 1;

for i=1:ng
    theta = i*pi/180.0;
    x1r = x1*cos(theta) - y1*sin(theta);
    y1r = x1*sin(theta) + y1*cos(theta);    
    x2r = x2*cos(theta) - y2*sin(theta);
    y2r = x2*sin(theta) + y2*cos(theta);
    x3r = x3*cos(theta) - y3*sin(theta);
    y3r = x3*sin(theta) + y3*cos(theta);
    r1 = [3 4 -L_out L_out L_out -L_out -L_out -L_out L_out L_out];
    r2 = [2 3 x1r x2r x3r y1r y2r y3r 0 0];
    gdm = [r1; r2]';
    g = decsg(gdm,'R1-P1',['R1'; 'P1']');
    
    %thermal model
    thermalmodel = createpde('thermal','steadystate');
    gm = geometryFromEdges(thermalmodel,g);
    pdegplot(thermalmodel,'EdgeLabels','on');
    thermalBC(thermalmodel,'Edge',1,'Temperature',1);
    thermalBC(thermalmodel,'Edge',2,'Temperature',1);
    thermalBC(thermalmodel,'Edge',3,'Temperature',0);
    thermalBC(thermalmodel,'Edge',4,'Temperature',0);
    thermalBC(thermalmodel,'Edge',5,'Temperature',0);
    thermalBC(thermalmodel,'Edge',6,'Temperature',1);
    thermalBC(thermalmodel,'Edge',7,'Temperature',1);                   
    thermalProperties(thermalmodel,"ThermalConductivity",1.0);
    mesh = generateMesh(thermalmodel);
    %mesh = generateMesh(thermalmodel,'Hmax',0.105);
    thermalresults = solve(thermalmodel);
    T = thermalresults.Temperature;
    gradTx =thermalresults.XGradients;
    gradTy =thermalresults.YGradients;
    %Linear elasticity model (Plane Stress)
    structuralmodel = createpde('structural','static-planestress');
    structuralmodel.Geometry = gm;
    structuralmodel.Mesh = mesh;
    structuralProperties(structuralmodel,'YoungsModulus',1, ...
                                     'PoissonsRatio',0.3, ...'
                                     'CTE',1.0);
    
    %%%Outter cylinder
    structuralBC(structuralmodel,'Edge',1,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',2,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',7,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',6,'Constraint','fixed');
    %%%Inner cylinder
    structuralBC(structuralmodel,'Edge',3,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',4,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',5,'Constraint','fixed');
    
    structuralBodyLoad(structuralmodel,'Temperature',thermalresults);
    structuralmodel.ReferenceTemperature = 0; %No affect on the solution!
    thermalstressresults = solve(structuralmodel);
    u = thermalstressresults.Displacement.x ;
    v = thermalstressresults.Displacement.y ;
    stress = thermalstressresults.VonMisesStress ;
    coord = thermalresults.Mesh.Nodes';
    x_coord = coord(:,1);
    y_coord = coord(:,2);
    %disp(size(x_coord));
    if (size(x_coord)<min_points)
        aaa = size(x_coord);
        min_points = aaa(1);
    end

    %writing the results
    id = fopen(name+'T'+string(i)+'.txt','w');
    fprintf(id,'%f\n',T);
    fclose(id);

    id = fopen(name+'dTdx'+string(i)+'.txt','w');
    fprintf(id,'%f\n',gradTx);
    fclose(id);

    id = fopen(name+'dTdy'+string(i)+'.txt','w');
    fprintf(id,'%f\n',gradTy);
    fclose(id);
    
    id = fopen(name+'u'+string(i)+'.txt','w');
    fprintf(id,'%f\n',u);
    fclose(id);
    
    id = fopen(name+'v'+string(i)+'.txt','w');
    fprintf(id,'%f\n',v);
    fclose(id);
    
    id = fopen(name+'x'+string(i)+'.txt','w');
    fprintf(id,'%f\n',x_coord);
    fclose(id);
    
    id = fopen(name+'y'+string(i)+'.txt','w');
    fprintf(id,'%f\n',y_coord);
    fclose(id);
end

display(min_points);

%%%%%%%

%clc;
%clear;

%%Hexagon
ng = 60; %number of geometries
name = "hexagon";
L = 0.3; %0.35; %0.35; %0.3; %0.5;
alpha = 60*pi/180.0;    
x1 = L; y1 = 0.0;
x2 = x1*cos(alpha) - y1*sin(alpha); y2 = x1*sin(alpha) + y1*cos(alpha);
x3 = x2*cos(alpha) - y2*sin(alpha); y3 = x2*sin(alpha) + y2*cos(alpha);  
x4 = x3*cos(alpha) - y3*sin(alpha); y4 = x3*sin(alpha) + y3*cos(alpha);
x5 = x4*cos(alpha) - y4*sin(alpha); y5 = x4*sin(alpha) + y4*cos(alpha);
x6 = x5*cos(alpha) - y5*sin(alpha); y6 = x5*sin(alpha) + y5*cos(alpha);
min_points = 10000; %Gross Points
L_out = 1;

for i=1:ng
    theta = i*pi/180.0;
    x1r = x1*cos(theta) - y1*sin(theta);
    y1r = x1*sin(theta) + y1*cos(theta);    
    x2r = x2*cos(theta) - y2*sin(theta);
    y2r = x2*sin(theta) + y2*cos(theta);
    x3r = x3*cos(theta) - y3*sin(theta);
    y3r = x3*sin(theta) + y3*cos(theta);
    x4r = x4*cos(theta) - y4*sin(theta);
    y4r = x4*sin(theta) + y4*cos(theta);
    x5r = x5*cos(theta) - y5*sin(theta);
    y5r = x5*sin(theta) + y5*cos(theta);
    x6r = x6*cos(theta) - y6*sin(theta);
    y6r = x6*sin(theta) + y6*cos(theta);
    r1 = [3 4 -L_out L_out L_out -L_out -L_out -L_out L_out L_out 0 0 0 0];
    r2 = [2 6 x1r x2r x3r x4r x5r x6r y1r y2r y3r y4r y5r y6r];
    gdm = [r1; r2]';
    g = decsg(gdm,'R1-P1',['R1'; 'P1']');
    
    %thermal model
    thermalmodel = createpde('thermal','steadystate');
    gm = geometryFromEdges(thermalmodel,g);
    pdegplot(thermalmodel,'EdgeLabels','on');
    thermalBC(thermalmodel,'Edge',1,'Temperature',1);
    thermalBC(thermalmodel,'Edge',2,'Temperature',1);
    thermalBC(thermalmodel,'Edge',3,'Temperature',0);
    thermalBC(thermalmodel,'Edge',4,'Temperature',0);
    thermalBC(thermalmodel,'Edge',5,'Temperature',0);
    thermalBC(thermalmodel,'Edge',6,'Temperature',0);
    thermalBC(thermalmodel,'Edge',7,'Temperature',0);
    thermalBC(thermalmodel,'Edge',8,'Temperature',1);
    thermalBC(thermalmodel,'Edge',9,'Temperature',1); 
    thermalBC(thermalmodel,'Edge',10,'Temperature',0); 
    thermalProperties(thermalmodel,"ThermalConductivity",1.0);
    mesh = generateMesh(thermalmodel);
    %mesh = generateMesh(thermalmodel,'Hmax',0.105);
    thermalresults = solve(thermalmodel);
    T = thermalresults.Temperature;
    gradTx =thermalresults.XGradients;
    gradTy =thermalresults.YGradients;
    %Linear elasticity model (Plane Stress)
    structuralmodel = createpde('structural','static-planestress');
    structuralmodel.Geometry = gm;
    structuralmodel.Mesh = mesh;
    structuralProperties(structuralmodel,'YoungsModulus',1, ...
                                     'PoissonsRatio',0.3, ...'
                                     'CTE',1.0);
    
    %%%Outter cylinder
    structuralBC(structuralmodel,'Edge',1,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',2,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',8,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',9,'Constraint','fixed');
    %%%Inner cylinder
    structuralBC(structuralmodel,'Edge',3,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',4,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',5,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',6,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',7,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',10,'Constraint','fixed');
    
    structuralBodyLoad(structuralmodel,'Temperature',thermalresults);
    structuralmodel.ReferenceTemperature = 0; %No affect on the solution!
    thermalstressresults = solve(structuralmodel);
    u = thermalstressresults.Displacement.x ;
    v = thermalstressresults.Displacement.y ;
    stress = thermalstressresults.VonMisesStress ;
    coord = thermalresults.Mesh.Nodes';
    x_coord = coord(:,1);
    y_coord = coord(:,2);
    %disp(size(x_coord));
    if (size(x_coord)<min_points)
        aaa = size(x_coord);
        min_points = aaa(1);
    end

    %writing the results
    id = fopen(name+'T'+string(i)+'.txt','w');
    fprintf(id,'%f\n',T);
    fclose(id);

    id = fopen(name+'dTdx'+string(i)+'.txt','w');
    fprintf(id,'%f\n',gradTx);
    fclose(id);

    id = fopen(name+'dTdy'+string(i)+'.txt','w');
    fprintf(id,'%f\n',gradTy);
    fclose(id);
    
    id = fopen(name+'u'+string(i)+'.txt','w');
    fprintf(id,'%f\n',u);
    fclose(id);
    
    id = fopen(name+'v'+string(i)+'.txt','w');
    fprintf(id,'%f\n',v);
    fclose(id);
    
    id = fopen(name+'x'+string(i)+'.txt','w');
    fprintf(id,'%f\n',x_coord);
    fclose(id);
    
    id = fopen(name+'y'+string(i)+'.txt','w');
    fprintf(id,'%f\n',y_coord);
    fclose(id);
end

display(min_points);

%%%%
%clc;
%clear;


%%Pentagon
ng = 72; %number of geometries
name = "pentagon";
L = 0.3; %0.35; %0.35; %0.3; %0.5;
alpha = 72*pi/180.0;    
x1 = L; y1 = 0.0;
x2 = x1*cos(alpha) - y1*sin(alpha); y2 = x1*sin(alpha) + y1*cos(alpha);
x3 = x2*cos(alpha) - y2*sin(alpha); y3 = x2*sin(alpha) + y2*cos(alpha);  
x4 = x3*cos(alpha) - y3*sin(alpha); y4 = x3*sin(alpha) + y3*cos(alpha);
x5 = x4*cos(alpha) - y4*sin(alpha); y5 = x4*sin(alpha) + y4*cos(alpha);
min_points = 10000; %Gross Points
L_out = 1;

for i=1:ng
    theta = i*pi/180.0;
    x1r = x1*cos(theta) - y1*sin(theta);
    y1r = x1*sin(theta) + y1*cos(theta);    
    x2r = x2*cos(theta) - y2*sin(theta);
    y2r = x2*sin(theta) + y2*cos(theta);
    x3r = x3*cos(theta) - y3*sin(theta);
    y3r = x3*sin(theta) + y3*cos(theta);
    x4r = x4*cos(theta) - y4*sin(theta);
    y4r = x4*sin(theta) + y4*cos(theta);
    x5r = x5*cos(theta) - y5*sin(theta);
    y5r = x5*sin(theta) + y5*cos(theta);
    r1 = [3 4 -L_out L_out L_out -L_out -L_out -L_out L_out L_out 0 0];
    r2 = [2 5 x1r x2r x3r x4r x5r y1r y2r y3r y4r y5r];
    gdm = [r1; r2]';
    g = decsg(gdm,'R1-P1',['R1'; 'P1']');
    
    %thermal model
    thermalmodel = createpde('thermal','steadystate');
    gm = geometryFromEdges(thermalmodel,g);
    pdegplot(thermalmodel,'EdgeLabels','on');
    thermalBC(thermalmodel,'Edge',1,'Temperature',1);
    thermalBC(thermalmodel,'Edge',2,'Temperature',1);
    thermalBC(thermalmodel,'Edge',3,'Temperature',0);
    thermalBC(thermalmodel,'Edge',4,'Temperature',0);
    thermalBC(thermalmodel,'Edge',5,'Temperature',0);
    thermalBC(thermalmodel,'Edge',6,'Temperature',0);
    thermalBC(thermalmodel,'Edge',7,'Temperature',0);
    thermalBC(thermalmodel,'Edge',8,'Temperature',1);
    thermalBC(thermalmodel,'Edge',9,'Temperature',1); 
    thermalProperties(thermalmodel,"ThermalConductivity",1.0);
    mesh = generateMesh(thermalmodel);
    %mesh = generateMesh(thermalmodel,'Hmax',0.105);
    thermalresults = solve(thermalmodel);
    T = thermalresults.Temperature;
    gradTx =thermalresults.XGradients;
    gradTy =thermalresults.YGradients;
    %Linear elasticity model (Plane Stress)
    structuralmodel = createpde('structural','static-planestress');
    structuralmodel.Geometry = gm;
    structuralmodel.Mesh = mesh;
    structuralProperties(structuralmodel,'YoungsModulus',1, ...
                                     'PoissonsRatio',0.3, ...'
                                     'CTE',1.0);
    
    %%%Outter cylinder
    structuralBC(structuralmodel,'Edge',1,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',2,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',8,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',9,'Constraint','fixed');
    %%%Inner cylinder
    structuralBC(structuralmodel,'Edge',3,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',4,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',5,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',6,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',7,'Constraint','fixed');

    structuralBodyLoad(structuralmodel,'Temperature',thermalresults);
    structuralmodel.ReferenceTemperature = 0; %No affect on the solution!
    thermalstressresults = solve(structuralmodel);
    u = thermalstressresults.Displacement.x ;
    v = thermalstressresults.Displacement.y ;
    stress = thermalstressresults.VonMisesStress ;
    coord = thermalresults.Mesh.Nodes';
    x_coord = coord(:,1);
    y_coord = coord(:,2);
    %disp(size(x_coord));
    if (size(x_coord)<min_points)
        aaa = size(x_coord);
        min_points = aaa(1);
    end

    %writing the results
    id = fopen(name+'T'+string(i)+'.txt','w');
    fprintf(id,'%f\n',T);
    fclose(id);

    id = fopen(name+'dTdx'+string(i)+'.txt','w');
    fprintf(id,'%f\n',gradTx);
    fclose(id);

    id = fopen(name+'dTdy'+string(i)+'.txt','w');
    fprintf(id,'%f\n',gradTy);
    fclose(id);
    
    id = fopen(name+'u'+string(i)+'.txt','w');
    fprintf(id,'%f\n',u);
    fclose(id);
    
    id = fopen(name+'v'+string(i)+'.txt','w');
    fprintf(id,'%f\n',v);
    fclose(id);
    
    id = fopen(name+'x'+string(i)+'.txt','w');
    fprintf(id,'%f\n',x_coord);
    fclose(id);
    
    id = fopen(name+'y'+string(i)+'.txt','w');
    fprintf(id,'%f\n',y_coord);
    fclose(id);
end

display(min_points);

%%%%
%clc;
%clear;

%%Heptagon
ng = 51; %number of geometries
name = "heptagon";
L = 0.3; %0.35; %0.35; %0.3; %0.5;
alpha = (360./7)*pi/180.0;    
x1 = L; y1 = 0.0;
x2 = x1*cos(alpha) - y1*sin(alpha); y2 = x1*sin(alpha) + y1*cos(alpha);
x3 = x2*cos(alpha) - y2*sin(alpha); y3 = x2*sin(alpha) + y2*cos(alpha);  
x4 = x3*cos(alpha) - y3*sin(alpha); y4 = x3*sin(alpha) + y3*cos(alpha);
x5 = x4*cos(alpha) - y4*sin(alpha); y5 = x4*sin(alpha) + y4*cos(alpha);
x6 = x5*cos(alpha) - y5*sin(alpha); y6 = x5*sin(alpha) + y5*cos(alpha);
x7 = x6*cos(alpha) - y6*sin(alpha); y7 = x6*sin(alpha) + y6*cos(alpha);

min_points = 10000; %Gross Points
L_out = 1;

for i=1:ng
    theta = i*pi/180.0;
    x1r = x1*cos(theta) - y1*sin(theta);
    y1r = x1*sin(theta) + y1*cos(theta);    
    x2r = x2*cos(theta) - y2*sin(theta);
    y2r = x2*sin(theta) + y2*cos(theta);
    x3r = x3*cos(theta) - y3*sin(theta);
    y3r = x3*sin(theta) + y3*cos(theta);
    x4r = x4*cos(theta) - y4*sin(theta);
    y4r = x4*sin(theta) + y4*cos(theta);
    x5r = x5*cos(theta) - y5*sin(theta);
    y5r = x5*sin(theta) + y5*cos(theta);
    x6r = x6*cos(theta) - y6*sin(theta);
    y6r = x6*sin(theta) + y6*cos(theta);
    x7r = x7*cos(theta) - y7*sin(theta);
    y7r = x7*sin(theta) + y7*cos(theta);
    r1 = [3 4 -L_out L_out L_out -L_out -L_out -L_out L_out L_out 0 0 0 0 0 0];
    r2 = [2 7 x1r x2r x3r x4r x5r x6r x7r y1r y2r y3r y4r y5r y6r y7r];
    gdm = [r1; r2]';
    g = decsg(gdm,'R1-P1',['R1'; 'P1']');
    
    %thermal model
    thermalmodel = createpde('thermal','steadystate');
    gm = geometryFromEdges(thermalmodel,g);
    pdegplot(thermalmodel,'EdgeLabels','on');
    thermalBC(thermalmodel,'Edge',1,'Temperature',1);
    thermalBC(thermalmodel,'Edge',2,'Temperature',1);
    thermalBC(thermalmodel,'Edge',3,'Temperature',0);
    thermalBC(thermalmodel,'Edge',4,'Temperature',0);
    thermalBC(thermalmodel,'Edge',5,'Temperature',0);
    thermalBC(thermalmodel,'Edge',6,'Temperature',0);
    thermalBC(thermalmodel,'Edge',7,'Temperature',0);
    thermalBC(thermalmodel,'Edge',8,'Temperature',1);
    thermalBC(thermalmodel,'Edge',9,'Temperature',1); 
    thermalBC(thermalmodel,'Edge',10,'Temperature',0); 
    thermalBC(thermalmodel,'Edge',11,'Temperature',0); 
    thermalProperties(thermalmodel,"ThermalConductivity",1.0);
    mesh = generateMesh(thermalmodel);
    %mesh = generateMesh(thermalmodel,'Hmax',0.105);
    thermalresults = solve(thermalmodel);
    T = thermalresults.Temperature;
    gradTx =thermalresults.XGradients;
    gradTy =thermalresults.YGradients;
    %Linear elasticity model (Plane Stress)
    structuralmodel = createpde('structural','static-planestress');
    structuralmodel.Geometry = gm;
    structuralmodel.Mesh = mesh;
    structuralProperties(structuralmodel,'YoungsModulus',1, ...
                                     'PoissonsRatio',0.3, ...'
                                     'CTE',1.0);
    
    %%%Outter cylinder
    structuralBC(structuralmodel,'Edge',1,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',2,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',8,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',9,'Constraint','fixed');
    %%%Inner cylinder
    structuralBC(structuralmodel,'Edge',3,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',4,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',5,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',6,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',7,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',10,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',11,'Constraint','fixed');
    
    structuralBodyLoad(structuralmodel,'Temperature',thermalresults);
    structuralmodel.ReferenceTemperature = 0; %No affect on the solution!
    thermalstressresults = solve(structuralmodel);
    u = thermalstressresults.Displacement.x ;
    v = thermalstressresults.Displacement.y ;
    stress = thermalstressresults.VonMisesStress ;
    coord = thermalresults.Mesh.Nodes';
    x_coord = coord(:,1);
    y_coord = coord(:,2);
    %disp(size(x_coord));
    if (size(x_coord)<min_points)
        aaa = size(x_coord);
        min_points = aaa(1);
    end

    %writing the results
    id = fopen(name+'T'+string(i)+'.txt','w');
    fprintf(id,'%f\n',T);
    fclose(id);

    id = fopen(name+'dTdx'+string(i)+'.txt','w');
    fprintf(id,'%f\n',gradTx);
    fclose(id);

    id = fopen(name+'dTdy'+string(i)+'.txt','w');
    fprintf(id,'%f\n',gradTy);
    fclose(id);
    
    id = fopen(name+'u'+string(i)+'.txt','w');
    fprintf(id,'%f\n',u);
    fclose(id);
    
    id = fopen(name+'v'+string(i)+'.txt','w');
    fprintf(id,'%f\n',v);
    fclose(id);
    
    id = fopen(name+'x'+string(i)+'.txt','w');
    fprintf(id,'%f\n',x_coord);
    fclose(id);
    
    id = fopen(name+'y'+string(i)+'.txt','w');
    fprintf(id,'%f\n',y_coord);
    fclose(id);
end

display(min_points);

%%%%%

%clc;
%clear;

%%Octagon
ng = 45; %number of geometries
name = "octagon";
L = 0.3; %0.35; %0.35; %0.3; %0.5;
alpha = 45*pi/180.0;    
x1 = L; y1 = 0.0;
x2 = x1*cos(alpha) - y1*sin(alpha); y2 = x1*sin(alpha) + y1*cos(alpha);
x3 = x2*cos(alpha) - y2*sin(alpha); y3 = x2*sin(alpha) + y2*cos(alpha);  
x4 = x3*cos(alpha) - y3*sin(alpha); y4 = x3*sin(alpha) + y3*cos(alpha);
x5 = x4*cos(alpha) - y4*sin(alpha); y5 = x4*sin(alpha) + y4*cos(alpha);
x6 = x5*cos(alpha) - y5*sin(alpha); y6 = x5*sin(alpha) + y5*cos(alpha);
x7 = x6*cos(alpha) - y6*sin(alpha); y7 = x6*sin(alpha) + y6*cos(alpha);
x8 = x7*cos(alpha) - y7*sin(alpha); y8 = x7*sin(alpha) + y7*cos(alpha);

min_points = 10000; %Gross Points
L_out = 1;

for i=1:ng
    theta = i*pi/180.0;
    x1r = x1*cos(theta) - y1*sin(theta);
    y1r = x1*sin(theta) + y1*cos(theta);    
    x2r = x2*cos(theta) - y2*sin(theta);
    y2r = x2*sin(theta) + y2*cos(theta);
    x3r = x3*cos(theta) - y3*sin(theta);
    y3r = x3*sin(theta) + y3*cos(theta);
    x4r = x4*cos(theta) - y4*sin(theta);
    y4r = x4*sin(theta) + y4*cos(theta);
    x5r = x5*cos(theta) - y5*sin(theta);
    y5r = x5*sin(theta) + y5*cos(theta);
    x6r = x6*cos(theta) - y6*sin(theta);
    y6r = x6*sin(theta) + y6*cos(theta);
    x7r = x7*cos(theta) - y7*sin(theta);
    y7r = x7*sin(theta) + y7*cos(theta);
    x8r = x8*cos(theta) - y8*sin(theta);
    y8r = x8*sin(theta) + y8*cos(theta);
    
    r1 = [3 4 -L_out L_out L_out -L_out -L_out -L_out L_out L_out 0 0 0 0 0 0 0 0];
    r2 = [2 8 x1r x2r x3r x4r x5r x6r x7r x8r y1r y2r y3r y4r y5r y6r y7r y8r];
    gdm = [r1; r2]';
    g = decsg(gdm,'R1-P1',['R1'; 'P1']');
    
    %thermal model
    thermalmodel = createpde('thermal','steadystate');
    gm = geometryFromEdges(thermalmodel,g);
    pdegplot(thermalmodel,'EdgeLabels','on');
    thermalBC(thermalmodel,'Edge',1,'Temperature',1);
    thermalBC(thermalmodel,'Edge',2,'Temperature',1);
    thermalBC(thermalmodel,'Edge',3,'Temperature',0);
    thermalBC(thermalmodel,'Edge',4,'Temperature',0);
    thermalBC(thermalmodel,'Edge',5,'Temperature',0);
    thermalBC(thermalmodel,'Edge',6,'Temperature',0);
    thermalBC(thermalmodel,'Edge',7,'Temperature',0);
    thermalBC(thermalmodel,'Edge',8,'Temperature',0);
    thermalBC(thermalmodel,'Edge',9,'Temperature',0); 
    thermalBC(thermalmodel,'Edge',10,'Temperature',0); 
    thermalBC(thermalmodel,'Edge',11,'Temperature',1);
    thermalBC(thermalmodel,'Edge',12,'Temperature',1);
    thermalProperties(thermalmodel,"ThermalConductivity",1.0);
    mesh = generateMesh(thermalmodel);
    %mesh = generateMesh(thermalmodel,'Hmax',0.105);
    thermalresults = solve(thermalmodel);
    T = thermalresults.Temperature;
    gradTx =thermalresults.XGradients;
    gradTy =thermalresults.YGradients;
    %Linear elasticity model (Plane Stress)
    structuralmodel = createpde('structural','static-planestress');
    structuralmodel.Geometry = gm;
    structuralmodel.Mesh = mesh;
    structuralProperties(structuralmodel,'YoungsModulus',1, ...
                                     'PoissonsRatio',0.3, ...'
                                     'CTE',1.0);
    
    %%%Outter cylinder
    structuralBC(structuralmodel,'Edge',1,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',2,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',11,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',12,'Constraint','fixed');
    %%%Inner cylinder
    structuralBC(structuralmodel,'Edge',3,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',4,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',5,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',6,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',7,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',10,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',8,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',9,'Constraint','fixed');
     
    structuralBodyLoad(structuralmodel,'Temperature',thermalresults);
    structuralmodel.ReferenceTemperature = 0; %No affect on the solution!
    thermalstressresults = solve(structuralmodel);
    u = thermalstressresults.Displacement.x ;
    v = thermalstressresults.Displacement.y ;
    stress = thermalstressresults.VonMisesStress ;
    coord = thermalresults.Mesh.Nodes';
    x_coord = coord(:,1);
    y_coord = coord(:,2);
    %disp(size(x_coord));
    if (size(x_coord)<min_points)
        aaa = size(x_coord);
        min_points = aaa(1);
    end

    %writing the results
    id = fopen(name+'T'+string(i)+'.txt','w');
    fprintf(id,'%f\n',T);
    fclose(id);

    id = fopen(name+'dTdx'+string(i)+'.txt','w');
    fprintf(id,'%f\n',gradTx);
    fclose(id);

    id = fopen(name+'dTdy'+string(i)+'.txt','w');
    fprintf(id,'%f\n',gradTy);
    fclose(id);
    
    id = fopen(name+'u'+string(i)+'.txt','w');
    fprintf(id,'%f\n',u);
    fclose(id);
    
    id = fopen(name+'v'+string(i)+'.txt','w');
    fprintf(id,'%f\n',v);
    fclose(id);
    
    id = fopen(name+'x'+string(i)+'.txt','w');
    fprintf(id,'%f\n',x_coord);
    fclose(id);
    
    id = fopen(name+'y'+string(i)+'.txt','w');
    fprintf(id,'%f\n',y_coord);
    fclose(id);
end

display(min_points);

%%%%

%clc;
%clear;

%%%%


%%Nanogan
ng = 40; %number of geometries
name = "nanogan";
L = 0.3; %0.35; %0.35; %0.3; %0.5;
alpha = 40*pi/180.0;    
x1 = L; y1 = 0.0;
x2 = x1*cos(alpha) - y1*sin(alpha); y2 = x1*sin(alpha) + y1*cos(alpha);
x3 = x2*cos(alpha) - y2*sin(alpha); y3 = x2*sin(alpha) + y2*cos(alpha);  
x4 = x3*cos(alpha) - y3*sin(alpha); y4 = x3*sin(alpha) + y3*cos(alpha);
x5 = x4*cos(alpha) - y4*sin(alpha); y5 = x4*sin(alpha) + y4*cos(alpha);
x6 = x5*cos(alpha) - y5*sin(alpha); y6 = x5*sin(alpha) + y5*cos(alpha);
x7 = x6*cos(alpha) - y6*sin(alpha); y7 = x6*sin(alpha) + y6*cos(alpha);
x8 = x7*cos(alpha) - y7*sin(alpha); y8 = x7*sin(alpha) + y7*cos(alpha);
x9 = x8*cos(alpha) - y8*sin(alpha); y9 = x8*sin(alpha) + y8*cos(alpha);

min_points = 10000; %Gross Points
L_out = 1;

for i=1:ng
    theta = i*pi/180.0;
    x1r = x1*cos(theta) - y1*sin(theta);
    y1r = x1*sin(theta) + y1*cos(theta);    
    x2r = x2*cos(theta) - y2*sin(theta);
    y2r = x2*sin(theta) + y2*cos(theta);
    x3r = x3*cos(theta) - y3*sin(theta);
    y3r = x3*sin(theta) + y3*cos(theta);
    x4r = x4*cos(theta) - y4*sin(theta);
    y4r = x4*sin(theta) + y4*cos(theta);
    x5r = x5*cos(theta) - y5*sin(theta);
    y5r = x5*sin(theta) + y5*cos(theta);
    x6r = x6*cos(theta) - y6*sin(theta);
    y6r = x6*sin(theta) + y6*cos(theta);
    x7r = x7*cos(theta) - y7*sin(theta);
    y7r = x7*sin(theta) + y7*cos(theta);
    x8r = x8*cos(theta) - y8*sin(theta);
    y8r = x8*sin(theta) + y8*cos(theta);
    x9r = x9*cos(theta) - y9*sin(theta);
    y9r = x9*sin(theta) + y9*cos(theta);
    r1 = [3 4 -L_out L_out L_out -L_out -L_out -L_out L_out L_out 0 0 0 0 0 0 0 0 0 0];
    r2 = [2 9 x1r x2r x3r x4r x5r x6r x7r x8r x9r y1r y2r y3r y4r y5r y6r y7r y8r y9r];
    gdm = [r1; r2]';
    g = decsg(gdm,'R1-P1',['R1'; 'P1']');
    
    %thermal model
    thermalmodel = createpde('thermal','steadystate');
    gm = geometryFromEdges(thermalmodel,g);
    pdegplot(thermalmodel,'EdgeLabels','on');
    thermalBC(thermalmodel,'Edge',1,'Temperature',1);
    thermalBC(thermalmodel,'Edge',2,'Temperature',1);
    thermalBC(thermalmodel,'Edge',3,'Temperature',0);
    thermalBC(thermalmodel,'Edge',4,'Temperature',0);
    thermalBC(thermalmodel,'Edge',5,'Temperature',0);
    thermalBC(thermalmodel,'Edge',6,'Temperature',0);
    thermalBC(thermalmodel,'Edge',7,'Temperature',0);
    thermalBC(thermalmodel,'Edge',8,'Temperature',0);
    thermalBC(thermalmodel,'Edge',9,'Temperature',0); 
    thermalBC(thermalmodel,'Edge',10,'Temperature',0); 
    thermalBC(thermalmodel,'Edge',11,'Temperature',1);
    thermalBC(thermalmodel,'Edge',12,'Temperature',1);
    thermalBC(thermalmodel,'Edge',13,'Temperature',0);
    
    thermalProperties(thermalmodel,"ThermalConductivity",1.0);
    mesh = generateMesh(thermalmodel);
    %mesh = generateMesh(thermalmodel,'Hmax',0.105);
    thermalresults = solve(thermalmodel);
    T = thermalresults.Temperature;
    gradTx =thermalresults.XGradients;
    gradTy =thermalresults.YGradients;
    %Linear elasticity model (Plane Stress)
    structuralmodel = createpde('structural','static-planestress');
    structuralmodel.Geometry = gm;
    structuralmodel.Mesh = mesh;
    structuralProperties(structuralmodel,'YoungsModulus',1, ...
                                     'PoissonsRatio',0.3, ...'
                                     'CTE',1.0);
    
    %%%Outter cylinder
    structuralBC(structuralmodel,'Edge',1,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',2,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',11,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',12,'Constraint','fixed');
    %%%Inner cylinder
    structuralBC(structuralmodel,'Edge',3,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',4,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',5,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',6,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',7,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',10,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',8,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',9,'Constraint','fixed');
    structuralBC(structuralmodel,'Edge',13,'Constraint','fixed');
     
    structuralBodyLoad(structuralmodel,'Temperature',thermalresults);
    structuralmodel.ReferenceTemperature = 0; %No affect on the solution!
    thermalstressresults = solve(structuralmodel);
    u = thermalstressresults.Displacement.x ;
    v = thermalstressresults.Displacement.y ;
    stress = thermalstressresults.VonMisesStress ;
    coord = thermalresults.Mesh.Nodes';
    x_coord = coord(:,1);
    y_coord = coord(:,2);
    %disp(size(x_coord));
    if (size(x_coord)<min_points)
        aaa = size(x_coord);
        min_points = aaa(1);
    end

    %writing the results
    id = fopen(name+'T'+string(i)+'.txt','w');
    fprintf(id,'%f\n',T);
    fclose(id);

    id = fopen(name+'dTdx'+string(i)+'.txt','w');
    fprintf(id,'%f\n',gradTx);
    fclose(id);

    id = fopen(name+'dTdy'+string(i)+'.txt','w');
    fprintf(id,'%f\n',gradTy);
    fclose(id);
    
    id = fopen(name+'u'+string(i)+'.txt','w');
    fprintf(id,'%f\n',u);
    fclose(id);
    
    id = fopen(name+'v'+string(i)+'.txt','w');
    fprintf(id,'%f\n',v);
    fclose(id);
    
    id = fopen(name+'x'+string(i)+'.txt','w');
    fprintf(id,'%f\n',x_coord);
    fclose(id);
    
    id = fopen(name+'y'+string(i)+'.txt','w');
    fprintf(id,'%f\n',y_coord);
    fclose(id);
end

display(min_points);
