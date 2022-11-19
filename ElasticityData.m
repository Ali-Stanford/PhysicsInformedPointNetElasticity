%Geometry
L = 0.35; %0.35; %0.3; %0.5;
x1 = -L;
y1 = -L;
x2 = L;
y2 = -L;
x3 = L;
y3 = L;
x4 = -L;
y4 = L;
min_points = 10000;

%for i=1:90
for i=91:180
    theta = i*pi/180.0;
    x1r = x1*cos(theta) - y1*sin(theta);
    y1r = x1*sin(theta) + y1*cos(theta);    
    x2r = x2*cos(theta) - y2*sin(theta);
    y2r = x2*sin(theta) + y2*cos(theta);
    x3r = x3*cos(theta) - y3*sin(theta);
    y3r = x3*sin(theta) + y3*cos(theta);
    x4r = x4*cos(theta) - y4*sin(theta);
    y4r = x4*sin(theta) + y4*cos(theta);
    r1 = [3 4 -1 1 1 -1 -1 -1 1 1];
    r2 = [3 4 x1r x2r x3r x4r y1r y2r y3r y4r];
    gdm = [r1; r2]';
    g = decsg(gdm,'R1-R2',['R1'; 'R2']');
    
    %thermal model
    thermalmodel = createpde('thermal','steadystate');
    gm = geometryFromEdges(thermalmodel,g);
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
    disp(size(x_coord));
    if (size(x_coord)<min_points)
        aaa = size(x_coord);
        min_points = aaa(1);
    end
    %figure
    %pdeplot(structuralmodel,"XYData",v,"Contour","on","ColorMap","hot");
    %title 'v displacement';

    %writing the results
    id = fopen('T'+string(i)+'.txt','w');
    fprintf(id,'%f\n',T);
    fclose(id);

    id = fopen('dTdx'+string(i)+'.txt','w');
    fprintf(id,'%f\n',gradTx);
    fclose(id);

    id = fopen('dTdy'+string(i)+'.txt','w');
    fprintf(id,'%f\n',gradTy);
    fclose(id);
    
    id = fopen('u'+string(i)+'.txt','w');
    fprintf(id,'%f\n',u);
    fclose(id);
    
    id = fopen('v'+string(i)+'.txt','w');
    fprintf(id,'%f\n',v);
    fclose(id);
    
    id = fopen('x'+string(i)+'.txt','w');
    fprintf(id,'%f\n',x_coord);
    fclose(id);
    
    id = fopen('y'+string(i)+'.txt','w');
    fprintf(id,'%f\n',y_coord);
    fclose(id);
end

display(min_points);
