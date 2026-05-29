classdef EEE_PROJECT_FROM_95_AND_96 < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        ShowPlotButton                  matlab.ui.control.Button
        BoundedTemperatureInitialEditField  matlab.ui.control.NumericEditField
        BoundedTemperatureInitialLabel  matlab.ui.control.Label
        ThicknessinmmEditField          matlab.ui.control.NumericEditField
        ThicknessinmmEditFieldLabel     matlab.ui.control.Label
        NumberofBlocksEnterbetween1to4EditField  matlab.ui.control.NumericEditField
        NumberofBlocksEnterbetween1to4Label  matlab.ui.control.Label
        HeatConvectionCoefficientEditField  matlab.ui.control.NumericEditField
        HeatConvectionCoefficientLabel  matlab.ui.control.Label
        ThermalConductivityEditField    matlab.ui.control.NumericEditField
        ThermalConductivityLabel        matlab.ui.control.Label
        Image                           matlab.ui.control.Image

        % NEW: Real-World Application UI Components
        RealWorldPanel                  matlab.ui.container.Panel
        AppScenarioLabel                matlab.ui.control.Label
        AppScenarioDropDown             matlab.ui.control.DropDown
        RunRealWorldButton              matlab.ui.control.Button
        RealWorldResultsLabel           matlab.ui.control.Label
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: ShowPlotButton
        function ShowPlotButtonPushed(app, event)
         
clear global
global nnd nel nodes DOF elDOF n Block
global geom conodesc convect nf Nodal_loads
global Length Width NoElemInXdirec NoElemInYdirec ElemSizeInXdirec ElemSizeInYdirec X_origin Y_origin
format long g

%% Inputs
% DOF
    nodes = 3;            % Number of nodes each element
    DOF = 1;              % Number of degrees of freedom
    elDOF = nodes*DOF;    % Element degrees of freedom
% Material
    ThermalConductivity =app.ThermalConductivityEditField.Value;           %Coefficient of thermal conductivity, Btu/(h-ft-F)
    Kyy = ThermalConductivity;
    h = app.HeatConvectionCoefficientEditField.Value;             %Heat convection coefficient, Btu/(h-ft-F)
    thick = app.ThicknessinmmEditField.Value;          %Beam thickness in mm
% Geometries
    Block = app.NumberofBlocksEnterbetween1to4EditField.Value;          % Number of rectangular block in the model
    a = 70/7;           % Constant geometry parameter
    Length = [a 60 60 60]; % Length of blocks
    Width = [70 a a a]; % Width of blocks
    NoElemInAparam = 2;            % Number of elements in each "a" parameter
    NoElemInXdirec = [1 6 6 6]*NoElemInAparam;% Number of elements in the x direction of each block
    NoElemInYdirec = [7 1 1 1]*NoElemInAparam;% Number of elements in the y direction of each block
    ElemSizeInXdirec = Length./NoElemInXdirec;	% Element size in the x direction
    ElemSizeInYdirec = Width./NoElemInYdirec;	% Element size in the x direction
    X_origin = [0 10 10 10]; % X origin of the global coordinate system
    Y_origin = [0 10 30 50]; % Y origin of the global coordinate system
    T3_mesh_final_term_axisymmetric

% Boundaries & Loading
    % Nodes subjected to given temperature source
    InitialBoundedTemp =app.BoundedTemperatureInitialEditField.Value; % *****F
    nf = ones(nnd,DOF);
    Nodal_loads= zeros(nnd,DOF);
    for i=1:nnd
        % Find nodes at x = 0
        if geom(i,1) == 0
            Nodal_loads(i,1) = InitialBoundedTemp;
        end
    end
    % Element edges subjected to convection of free-steam temperature
    T_inf = 20; % F
    convect = zeros(nel,nodes);
    k = [2 3 1];
    for i=1:nel
        for j=1:nodes
            % Find x = 70
            if all(geom(conodesc(i,[j k(j)]),1) == 70)
                convect(i,j) = 1;
            end
            
            % Find y = long block top/bottom & x >= 10
            find_y = [0 70];
            for p = 1:length(find_y)
                isConvect1(p) = all(geom(conodesc(i,[j k(j)]),2) == find_y(p));
            end
            if any(isConvect1) && all(geom(conodesc(i,[j k(j)]),1) >= 0)
                  convect(i,j) = 1;
            end
            
            % Find y = long block top/bottom & x >= 10
            find_y = [10 20 30 40 50 60];
            for p = 1:length(find_y)
                isConvect2(p) = all(geom(conodesc(i,[j k(j)]),2) == find_y(p));
            end
            if any(isConvect2) && all(geom(conodesc(i,[j k(j)]),1) >= 10)
                convect(i,j) = 1;
            end
            
            % Find x = 10 & y is not inside block
            find_y = [0 10; 20 30; 40 50; 60 70];
            for p = 1:4
                isConvect3(p) = all([geom(conodesc(i,[j k(j)]),2) >= find_y(p,1);
                                    geom(conodesc(i,[j k(j)]),2) <= find_y(p,2)]);
            end
            if any(isConvect3) && all(geom(conodesc(i,[j k(j)]),1) == 10)
                convect(i,j) = 1;
            end 
        end
    end


%% Pre-processing
% Form the material property matrix
dee = [ThermalConductivity    0 ;
         0  Kyy];
     
% Counting of the free degrees of freedom
n=0;
for i=1:nnd
    for j=1:DOF
        if nf(i,j) ~= 0
            n=n+1;
            nf(i,j)=n;
        end
    end
end

%% Matrix assembly
% Assemble the global stiffness matrix and convection force vector
kk = zeros(n,n); %
fg = zeros(n,DOF);
for i=1:nel
    [bee,cee,fee,g,A] = elem_T3_heat(i); % Form stiffness matrix(bee), convection vector (cee)
                                         % element force function(fee), and
                                         % steering vector(g)
    
    % Stiffness matrix                                     
    kc = thick*A*bee'*dee*bee;      % Compute conduction matrix
    kh = thick*h*cee;               % Compute convection matrix
    ke = kc + kh;                   % Compute local stiffness matrix
    kk = form_kk(kk,ke,g);          % Assemble global stiffness matrix
    
    % Force vector
    fe = thick*h*T_inf*fee;
    fg = form_fg(fg,fe,g);
end

% Modify the global force vector due to boundaries
idT = find(Nodal_loads~=0);
for i=1:nnd
    if ~any(idT==i)
        fg(i) = fg(i) - sum(kk(i,idT)*Nodal_loads(idT));
    end
end
for i=1:nnd
    if Nodal_loads(i) ~= 0
        fg(i) = Nodal_loads(i); % Given temperature source
    end
end
[kk(idT,:),kk(:,idT)] = deal(0);
for i=1:length(idT)
    kk(idT(i),idT(i)) = deal(1);
end

%% Solve Problem
tg = kk\fg ; % solve for unknown displacements

figure(1)
    patch('Faces', conodesc, 'Vertices', geom, 'FaceVertexCData',tg, ...
    'Facecolor','interp','Marker','.')
    colorbar;
    set(gcf,'position',[400 500 500 400])
    title('Temperature Profile ^{o}C')
    xlabel('X (mm)')
    ylabel('Y (mm)')    



function T3_mesh_final_term_axisymmetric
%global nnd NoElemInXdirec NoElemInYdirec ElemSizeInXdirec ElemSizeInYdirec X_origin Y_origin geom conodesc nel Block
    k = 0;
    for b = 1:Block
        for i = 1:NoElemInXdirec(b)
            for j=1:NoElemInYdirec(b)
                k = k + 1;
                if b == 1
                    n1 = j + (i-1)*(2*NoElemInYdirec(b)+1);
                    n2 = j + i*(2*NoElemInYdirec(b)+1);
                    n3 = n1 + 1;
                    n4 = n2 + 1;
                    n5 = j + (2*i-1)*NoElemInYdirec(b)+i;
                    geom(n1,:) = [(i-1)*ElemSizeInXdirec(b)+X_origin(b)	(j-1)*ElemSizeInYdirec(b)+Y_origin(b)];
                    geom(n3,:) = [(i-1)*ElemSizeInXdirec(b)+X_origin(b)	j*ElemSizeInYdirec(b)+Y_origin(b)    ];
                else
                    if i == 1
                        n1 = find(and(geom(:,1)==X_origin(b), ...
                             round(geom(:,2),2)==round(Y_origin(b)+(j-1)*ElemSizeInYdirec(b),2)));
                        n2 = nnd + j + NoElemInYdirec(b);
                        n3 = n1 + 1;
                        n4 = n2 + 1;
                        n5 = nnd + j;
                    else
                        n1 = nnd + j + (i-2)*(2*NoElemInYdirec(b)+1) + NoElemInYdirec(b);
                        n2 = nnd + j + (i-1)*(2*NoElemInYdirec(b)+1) + NoElemInYdirec(b);
                        n3 = n1 + 1;
                        n4 = n2 + 1;
                        n5 = nnd + j + (2*i-3)*NoElemInYdirec(b) + i-1 + NoElemInYdirec(b);
                    end
                end
                nnd_b(b) = n4;
                geom(n2,:) = [i*ElemSizeInXdirec(b)+X_origin(b)      (j-1)*ElemSizeInYdirec(b)+Y_origin(b)];
                geom(n4,:) = [i*ElemSizeInXdirec(b)+X_origin(b)   	j*ElemSizeInYdirec(b)+Y_origin(b)    ];
                geom(n5,:) = mean([geom(n1,:); geom(n4,:)]);
                nel = 4*k;
                m = nel -3;
                conodesc(m,:) = [n1 n2 n5];
                m = nel -2;
                conodesc(m,:) = [n2 n4 n5];
                m = nel -1;
                conodesc(m,:) = [n4 n3 n5];
                conodesc(nel,:) = [n3 n1 n5];
            end
        end
        nnd = n4;
    end
end

% This function returns the coordinates of the nodes of element i
% and its steering vector
function [bee,cee,fee,g,A] = elem_T3_heat(i)
%global nodes DOF geom conodesc nf convect

    x1 = geom(conodesc(i,1),1); y1 = geom(conodesc(i,1),2);
    x2 = geom(conodesc(i,2),1); y2 = geom(conodesc(i,2),2);
    x3 = geom(conodesc(i,3),1); y3 = geom(conodesc(i,3),2);
    
    A = (0.5)*det([1 x1 y1; ...
                   1 x2 y2; ...
                   1 x3 y3]);
               
    L = [((x2-x1)^2 + (y2-y1)^2)^.5 ...
         ((x3-x2)^2 + (y3-y2)^2)^.5 ...
         ((x1-x3)^2 + (y1-y3)^2)^.5];
    
    m11 = (x2*y3 - x3*y2)/(2*A);
    m21 = (x3*y1 - x1*y3)/(2*A);
    m31 = (x1*y2 - y1*x2)/(2*A);
    m12 = (y2 - y3)/(2*A);
    m22 = (y3 - y1)/(2*A);
    m32 = (y1 - y2)/(2*A);
    m13 = (x3 - x2)/(2*A);
    m23 = (x1 - x3)/(2*A);
    m33 = (x2 -x1)/(2*A);
    
    % Conduction stiffness matrix
    bee = [m12 m22 m32;
           m13 m23 m33];
    
    % Convection stiffness matrix
    cee = 1/6*(convect(i,1)*L(1)*[2 1 0; 1 2 0; 0 0 0] + ...
               convect(i,2)*L(2)*[0 0 0; 0 2 1; 0 1 2] + ...
               convect(i,3)*L(3)*[2 0 1; 0 0 0; 1 0 2]);  
  
    % Heat force vector
    fee = 1/2*(convect(i,1)*L(1)*[1 1 0]' + ...
               convect(i,2)*L(2)*[0 1 1]' + ...
               convect(i,3)*L(3)*[1 0 1]');
    
    % Steering vector
    l=0;
    for k=1:nodes
        for j=1:DOF
            l=l+1;
            g(l)=nf(conodesc(i,k),j);
        end
    end
end

% This function assembles the global stiffness matrix
function kk = form_kk(kk,ke,g)
%global elDOF
    for i=1:elDOF
        if g(i) ~= 0
            for j=1:elDOF
                if g(j) ~= 0
                    kk(g(i),g(j))= kk(g(i),g(j)) + ke(i,j);
                end
            end
        end
    end
end

% This function assembles the global force vector
function fg = form_fg(fg,fe,g)
%global elDOF
    for i=1:elDOF
        if g(i) ~= 0
            fg(g(i))= fg(g(i)) + fe(i);
        end
    end
end
        end

        % NEW: Real-World Heat Transfer Application Button Callback
        function RunRealWorldButtonPushed(app, event)
            scenario = app.AppScenarioDropDown.Value;

            % Shared FEM solver (same numerical method as main solver)
            % Uses T3 triangular elements, conduction + convection
            % Each scenario sets its own physics parameters and geometry,
            % then runs through the same FEM assembly and solve pipeline.

            switch scenario

                % SCENARIO 1: Electronic Device Cooling
                % A CPU chip (heat source at base) cooled by forced air
                % convection on the exposed surfaces.
                % Geometry: 20mm x 10mm silicon substrate (2D cross-section)
                case 'Electronic Device Cooling'

                    % Physics Parameters 
                    k_si        = 148;      % Thermal conductivity of silicon [W/(m·K)]
                    h_forced    = 250;      % Forced-air convection coefficient [W/(m²·K)]
                    T_ambient   = 25;       % Ambient air temperature [°C]
                    T_chip_base = 85;       % Heat source (chip junction) temperature [°C]
                    thick_dev   = 0.005;    % Device thickness (5 mm) [m]

                    % Mesh Parameters 
                    Lx = 0.020;   % Domain width  [m] (20 mm)
                    Ly = 0.010;   % Domain height [m] (10 mm)
                    nx = 10;      % Elements in x
                    ny = 5;       % Elements in y

                    titleStr  = 'Electronic Device Cooling — Temperature (°C)';
                    xLabel    = 'X (m)  [Chip Width]';
                    yLabel    = 'Y (m)  [Chip Thickness]';
                    figNum    = 2;

                    % Boundary conditions:
                    %   Bottom edge (y=0)  → fixed T = T_chip_base (heat source)
                    %   Top / Left / Right → convection with h_forced, T_ambient
                    fixedEdgeFn   = @(x,y) (abs(y) < 1e-10);          % y = 0
                    convectEdgeFn = @(x,y) (abs(y-Ly)<1e-10) | ...    % y = Ly
                                           (abs(x)<1e-10)    | ...    % x = 0
                                           (abs(x-Lx)<1e-10);         % x = Lx
                    T_fixed = T_chip_base;
                    T_inf   = T_ambient;

                    resultMsg = sprintf( ...
                        ['Electronic Device Cooling\n' ...
                         '  Material : Silicon\n' ...
                         '  k = %.0f W/(m·K)\n' ...
                         '  h = %.0f W/(m²·K)\n' ...
                         '  T_base = %.0f °C  |  T_ambient = %.0f °C\n' ...
                         '  Domain: %.0f mm × %.0f mm\n' ...
                         '  → See Figure %d for temperature map.'], ...
                        k_si, h_forced, T_chip_base, T_ambient, ...
                        Lx*1e3, Ly*1e3, figNum);

                % SCENARIO 2: Building Wall Heating
                % A concrete wall with interior heating surface and
                % exterior cold-air exposure (winter insulation study).
                % Geometry: 300mm thick × 1000mm tall wall cross-section
                case 'Building Wall Heating'

                    % Physics Parameters
                    k_concrete  = 1.7;      % Thermal conductivity of concrete [W/(m·K)]
                    h_ext       = 20;       % Exterior wind-driven convection [W/(m²·K)]
                    T_interior  = 20;       % Interior surface (heated side) [°C]
                    T_exterior  = -10;      % Exterior ambient temperature [°C]
                    thick_dev   = 1.0;      % Out-of-plane depth (1 m strip) [m]

                    % Mesh Parameters 
                    Lx = 0.300;   % Wall thickness [m] (300 mm)
                    Ly = 1.000;   % Wall height    [m] (1000 mm)
                    nx = 6;       % Elements in x (through thickness)
                    ny = 10;      % Elements in y (along height)

                    titleStr  = 'Building Wall Heating — Temperature (°C)';
                    xLabel    = 'X (m)  [Through Wall Thickness]';
                    yLabel    = 'Y (m)  [Wall Height]';
                    figNum    = 3;

                    % Boundary conditions:
                    %   Left edge  (x=0)  → fixed T = T_interior (heated room side)
                    %   Right edge (x=Lx) → convection h_ext, T_exterior (outside)
                    %   Top / Bottom      → adiabatic (no flux, default)
                    fixedEdgeFn   = @(x,y) (abs(x) < 1e-10);          % x = 0 (interior)
                    convectEdgeFn = @(x,y) (abs(x-Lx) < 1e-10);       % x = Lx (exterior)
                    T_fixed = T_interior;
                    T_inf   = T_exterior;
                    h_forced    = h_ext;
                    k_si        = k_concrete;

                    resultMsg = sprintf( ...
                        ['Building Wall Heating\n' ...
                         '  Material : Concrete\n' ...
                         '  k = %.1f W/(m·K)\n' ...
                         '  h_ext = %.0f W/(m²·K)\n' ...
                         '  T_interior = %.0f °C  |  T_exterior = %.0f °C\n' ...
                         '  Wall: %.0f mm thick × %.0f mm tall\n' ...
                         '  → See Figure %d for temperature map.'], ...
                        k_concrete, h_ext, T_interior, T_exterior, ...
                        Lx*1e3, Ly*1e3, figNum);

                otherwise
                    app.RealWorldResultsLabel.Text = 'Please select a valid scenario.';
                    return
            end

            %  SHARED FEM SOLVER  (same numerical method as main code)
            %  T3 triangular elements — conduction + convection
            % Build rectangular mesh of T3 elements (2 triangles / quad)
            dx = Lx / nx;
            dy = Ly / ny;

            % Node coordinates
            numNodes = (nx+1)*(ny+1);
            nodeCoords = zeros(numNodes, 2);
            idx = 1;
            for ix = 0:nx
                for iy = 0:ny
                    nodeCoords(idx,:) = [ix*dx, iy*dy];
                    idx = idx + 1;
                end
            end

            nodeID = @(ix,iy) ix*(ny+1) + iy + 1; % 1-based

            % Element connectivity (split each quad into 2 T3 triangles)
            numElem = 2*nx*ny;
            elemConn = zeros(numElem, 3);
            eIdx = 1;
            for ix = 0:nx-1
                for iy = 0:ny-1
                    n1 = nodeID(ix,   iy);
                    n2 = nodeID(ix+1, iy);
                    n3 = nodeID(ix+1, iy+1);
                    n4 = nodeID(ix,   iy+1);
                    elemConn(eIdx,:)   = [n1 n2 n3];  % lower-right triangle
                    elemConn(eIdx+1,:) = [n1 n3 n4];  % upper-left  triangle
                    eIdx = eIdx + 2;
                end
            end

            % Material matrix
            D = [k_si 0; 0 k_si];

            % Identify boundary conditions
            nDOF    = numNodes;
            isFixed = false(numNodes,1);
            T_bc    = zeros(numNodes,1);
            for nd = 1:numNodes
                xn = nodeCoords(nd,1);
                yn = nodeCoords(nd,2);
                if fixedEdgeFn(xn, yn)
                    isFixed(nd) = true;
                    T_bc(nd)    = T_fixed;
                end
            end

            % Global stiffness and force
            KK = zeros(nDOF, nDOF);
            FF = zeros(nDOF, 1);

            for e = 1:numElem
                ni = elemConn(e,1);
                nj = elemConn(e,2);
                nk = elemConn(e,3);

                xi = nodeCoords(ni,:);
                xj = nodeCoords(nj,:);
                xk = nodeCoords(nk,:);

                % Element area
                Ae = 0.5 * abs(det([1 xi; 1 xj; 1 xk]));
                if Ae < 1e-20, continue; end

                % Shape function gradients
                b1 = (xj(2)-xk(2))/(2*Ae);
                b2 = (xk(2)-xi(2))/(2*Ae);
                b3 = (xi(2)-xj(2))/(2*Ae);
                c1 = (xk(1)-xj(1))/(2*Ae);
                c2 = (xi(1)-xk(1))/(2*Ae);
                c3 = (xj(1)-xi(1))/(2*Ae);

                B = [b1 b2 b3;
                     c1 c2 c3];

                % Conduction matrix
                Kc = thick_dev * Ae * (B' * D * B);

                % Convection on edges — check each edge of this element
                Kh = zeros(3,3);
                fh = zeros(3,1);
                edges = [1 2; 2 3; 3 1];  % local node pairs per edge
                locNodes = [ni nj nk];

                for ed = 1:3
                    ea = locNodes(edges(ed,1));
                    eb = locNodes(edges(ed,2));
                    xa = nodeCoords(ea,:);
                    xb = nodeCoords(eb,:);

                    % Midpoint test for convection edge
                    xmid = 0.5*(xa + xb);
                    if convectEdgeFn(xmid(1), xmid(2))
                        Le = norm(xb - xa);
                        % Convection stiffness (consistent)
                        Kh_e = (h_forced * Le * thick_dev / 6) * [2 1; 1 2];
                        fh_e = (h_forced * T_inf * Le * thick_dev / 2) * [1; 1];
                        locs = edges(ed,:);
                        Kh(locs, locs) = Kh(locs, locs) + Kh_e;
                        fh(locs)       = fh(locs)       + fh_e;
                    end
                end

                % Assemble into global system
                dofs = [ni nj nk];
                KK(dofs, dofs) = KK(dofs, dofs) + Kc + Kh;
                FF(dofs)       = FF(dofs)       + fh;
            end

            % Apply fixed-temperature BCs (penalty / direct substitution)
            fixedDOFs = find(isFixed);
            for nd = fixedDOFs'
                FF = FF - KK(:, nd) * T_bc(nd);
            end
            for nd = fixedDOFs'
                KK(nd, :) = 0;
                KK(:, nd) = 0;
                KK(nd, nd) = 1;
                FF(nd)     = T_bc(nd);
            end

            % Solve
            T_sol = KK \ FF;

            % --- Plot ---
            figure(figNum)
            patch('Faces', elemConn, 'Vertices', nodeCoords, ...
                  'FaceVertexCData', T_sol, ...
                  'FaceColor', 'interp', 'EdgeColor', 'none')
            colorbar
            colormap(jet)
            title(titleStr)
            xlabel(xLabel)
            ylabel(yLabel)
            set(gcf, 'Position', [450 100 560 420])

            % Add max/min annotation
            T_max = max(T_sol);
            T_min = min(T_sol);
            annotation_str = sprintf('T_{max}=%.1f°C\nT_{min}=%.1f°C', T_max, T_min);
            text(0.02, 0.95, annotation_str, 'Units','normalized', ...
                 'FontSize', 10, 'Color', 'w', 'FontWeight','bold', ...
                 'VerticalAlignment','top')

            % Update result label in UI
            app.RealWorldResultsLabel.Text = resultMsg;
        end

    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Get the file path for locating images
            pathToMLAPP = fileparts(mfilename('Project'));

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 650];
            app.UIFigure.Name = 'MATLAB App';

            % Create Image
            app.Image = uiimage(app.UIFigure);
            app.Image.BackgroundColor = [0.1882 0.0118 0.1882];
            app.Image.Position = [1 1 640 650];
            app.Image.ImageSource = fullfile(pathToMLAPP, '360_F_487461421_L0q4P5iS2nTd4JfRlYJjKlgLRPZXP4sY.png');

            % Create ThermalConductivityLabel
            app.ThermalConductivityLabel = uilabel(app.UIFigure);
            app.ThermalConductivityLabel.HorizontalAlignment = 'center';
            app.ThermalConductivityLabel.FontSize = 14;
            app.ThermalConductivityLabel.FontWeight = 'bold';
            app.ThermalConductivityLabel.FontColor = [1 0.6 0];
            app.ThermalConductivityLabel.Position = [45 428 116 51];
            app.ThermalConductivityLabel.Text = {'Thermal '; 'Conductivity:'};

            % Create ThermalConductivityEditField
            app.ThermalConductivityEditField = uieditfield(app.UIFigure, 'numeric');
            app.ThermalConductivityEditField.BackgroundColor = [1 0.6 0];
            app.ThermalConductivityEditField.FontColor = [1 1 1];
            app.ThermalConductivityEditField.Position = [29 378 149 51];

            % Create HeatConvectionCoefficientLabel
            app.HeatConvectionCoefficientLabel = uilabel(app.UIFigure);
            app.HeatConvectionCoefficientLabel.HorizontalAlignment = 'center';
            app.HeatConvectionCoefficientLabel.FontSize = 14;
            app.HeatConvectionCoefficientLabel.FontWeight = 'bold';
            app.HeatConvectionCoefficientLabel.FontColor = [1 0.6 0];
            app.HeatConvectionCoefficientLabel.Position = [464 435 162 34];
            app.HeatConvectionCoefficientLabel.Text = {'Heat '; 'Convection Coefficient:'};

            % Create HeatConvectionCoefficientEditField
            app.HeatConvectionCoefficientEditField = uieditfield(app.UIFigure, 'numeric');
            app.HeatConvectionCoefficientEditField.BackgroundColor = [1 0.6 0];
            app.HeatConvectionCoefficientEditField.FontColor = [1 1 1];
            app.HeatConvectionCoefficientEditField.FontSize = 14;
            app.HeatConvectionCoefficientEditField.FontWeight = 'bold';
            app.HeatConvectionCoefficientEditField.Position = [464 378 162 50];

            % Create NumberofBlocksEnterbetween1to4Label
            app.NumberofBlocksEnterbetween1to4Label = uilabel(app.UIFigure);
            app.NumberofBlocksEnterbetween1to4Label.HorizontalAlignment = 'right';
            app.NumberofBlocksEnterbetween1to4Label.FontSize = 14;
            app.NumberofBlocksEnterbetween1to4Label.FontColor = [1 1 1];
            app.NumberofBlocksEnterbetween1to4Label.Position = [160 242 115 52];
            app.NumberofBlocksEnterbetween1to4Label.Text = {'Number of Blocks'; '(Enter between '; '1 to 4):'};

            % Create NumberofBlocksEnterbetween1to4EditField
            app.NumberofBlocksEnterbetween1to4EditField = uieditfield(app.UIFigure, 'numeric');
            app.NumberofBlocksEnterbetween1to4EditField.BackgroundColor = [0.1 0.05 0.15];
            app.NumberofBlocksEnterbetween1to4EditField.FontColor = [1 1 1];
            app.NumberofBlocksEnterbetween1to4EditField.FontWeight = 'bold';
            app.NumberofBlocksEnterbetween1to4EditField.Position = [290 242 104 48];

            % Create ThicknessinmmEditFieldLabel
            app.ThicknessinmmEditFieldLabel = uilabel(app.UIFigure);
            app.ThicknessinmmEditFieldLabel.HorizontalAlignment = 'right';
            app.ThicknessinmmEditFieldLabel.FontSize = 14;
            app.ThicknessinmmEditFieldLabel.FontColor = [1 0.5 0.3];
            app.ThicknessinmmEditFieldLabel.Position = [26 221 119 22];
            app.ThicknessinmmEditFieldLabel.Text = 'Thickness(in mm):';

            % Create ThicknessinmmEditField
            app.ThicknessinmmEditField = uieditfield(app.UIFigure, 'numeric');
            app.ThicknessinmmEditField.BackgroundColor = [1 0.5 0.3];
            app.ThicknessinmmEditField.FontColor = [1 1 1];
            app.ThicknessinmmEditField.Position = [13 180 148 42];

            % Create BoundedTemperatureInitialLabel
            app.BoundedTemperatureInitialLabel = uilabel(app.UIFigure);
            app.BoundedTemperatureInitialLabel.HorizontalAlignment = 'center';
            app.BoundedTemperatureInitialLabel.FontSize = 14;
            app.BoundedTemperatureInitialLabel.FontColor = [1 0.5 0.3];
            app.BoundedTemperatureInitialLabel.Position = [473 221 143 34];
            app.BoundedTemperatureInitialLabel.Text = {'Bounded Temperature'; '(Initial):'};

            % Create BoundedTemperatureInitialEditField
            app.BoundedTemperatureInitialEditField = uieditfield(app.UIFigure, 'numeric');
            app.BoundedTemperatureInitialEditField.BackgroundColor = [1 0.5 0.3];
            app.BoundedTemperatureInitialEditField.FontColor = [1 1 1];
            app.BoundedTemperatureInitialEditField.Position = [473 180 143 42];

            % Create ShowPlotButton
            app.ShowPlotButton = uibutton(app.UIFigure, 'push');
            app.ShowPlotButton.ButtonPushedFcn = createCallbackFcn(app, @ShowPlotButtonPushed, true);
            app.ShowPlotButton.FontName = 'Bell MT';
            app.ShowPlotButton.FontSize = 24;
            app.ShowPlotButton.FontWeight = 'bold';
            app.ShowPlotButton.BackgroundColor = [0.9 0.4 0];  
            app.ShowPlotButton.FontColor = [1 1 1];
            app.ShowPlotButton.Position = [233 321 176 58];
            app.ShowPlotButton.Text = 'Show Plot:';

            % NEW: Real-World Application Panel
            % Positioned below the existing controls

            % Panel container
            app.RealWorldPanel = uipanel(app.UIFigure);
            app.RealWorldPanel.Title = 'Real-World Heat Transfer Applications';
            app.RealWorldPanel.FontSize = 12;
            app.RealWorldPanel.FontWeight = 'bold';
            app.RealWorldPanel.ForegroundColor = [1 1 1];
            app.RealWorldPanel.BackgroundColor = [0.12 0.06 0.18];
            app.RealWorldPanel.Position = [10 1 620 115];

            % Scenario label
            app.AppScenarioLabel = uilabel(app.RealWorldPanel);
            app.AppScenarioLabel.Text = 'Select Scenario:';
            app.AppScenarioLabel.FontSize = 11;
            app.AppScenarioLabel.FontWeight = 'bold';
            app.AppScenarioLabel.FontColor = [0 1 1];
            app.AppScenarioLabel.Position = [5 60 110 22];

            % Scenario drop-down
            app.AppScenarioDropDown = uidropdown(app.RealWorldPanel);
            app.AppScenarioDropDown.Items = {'Electronic Device Cooling', 'Building Wall Heating'};
            app.AppScenarioDropDown.Value = 'Electronic Device Cooling';
            app.AppScenarioDropDown.FontSize = 11;
            app.AppScenarioDropDown.Position = [118 57 210 26];

            % Run button
            app.RunRealWorldButton = uibutton(app.RealWorldPanel, 'push');
            app.RunRealWorldButton.Text = 'Run Application';
            app.RunRealWorldButton.ButtonPushedFcn = createCallbackFcn(app, @RunRealWorldButtonPushed, true);
            app.RunRealWorldButton.FontSize = 11;
            app.RunRealWorldButton.FontWeight = 'bold';
            app.RunRealWorldButton.FontColor = [0.6353 0.0784 0.1843];
            app.RunRealWorldButton.Position = [340 55 120 30];

            % Results label (shows scenario summary after run)
            app.RealWorldResultsLabel = uilabel(app.RealWorldPanel);
            app.RealWorldResultsLabel.Text = 'Results will appear here after running a scenario.';
            app.RealWorldResultsLabel.FontSize = 8;
            app.RealWorldResultsLabel.FontColor = [1 1 0.5];
            app.RealWorldResultsLabel.Position = [485 2 135 103];
            app.RealWorldResultsLabel.WordWrap = 'on';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = EEE_PROJECT_FROM_95_AND_96

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end