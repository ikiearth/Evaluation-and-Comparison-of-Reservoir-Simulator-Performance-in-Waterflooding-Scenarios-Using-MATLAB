% == SKENARIO A2-G3-M == %

mrstModule add ad-core ad-blackoil ad-props mrst-gui
close all;

gravity on
profile on -timer cpu

% Conversion factors
psi_to_pascal    = 6894.76;      % PSI to Pascal
pascal_to_psi    = 0.000145038;  % Pascal to PSI
feet_to_meter    = 0.3048;       % Feet to Meters
stb_to_m3        = 0.158987;     % Stb to cubic meter
m3_to_stb        = 6.28981;      % Cubic meter to stb
cp_to_pas        = 0.001;        % Centipoise to Pascal seconds

% == Reservoir Geometry == %

% Number of Grid Bloks
% Grid dimensions
Nx = 70;  % Number of divisions in the x-direction
Ny = 70;  % Number of divisions in the y-direction
Nz = 3;   % Number of divisions in the z-direction

% Length calculations (in meters)
Lx = Nx * 142.8571428571 * feet_to_meter;  % Convert from Feet to Meters
Ly = Ny * 142.8571428571 * feet_to_meter;  % Convert from Feet to Meters
Lz = 1;                          % Length in the z-direction (in meters)

% Dimension of Grid Blocks
gridDimension = [Nx, Ny, Nz];        % Grid dimensions (Nx, Ny, Nz)
pdim = [Lx, Ly, Lz];                 % Physical dimensions (Lx, Ly, Lz)

% Create Grid (Cartesian)
G = cartGrid(gridDimension, pdim);   % Create the Cartesian Grid
G = computeGeometry(G);              % Compute the geometry of the grid

% Grid Thickness (in meters)
LayerThickness = [20, 30, 50] * feet_to_meter;  % Convert thickness from feet to meters

% Grid Top (in meters)
gridTop = 8325 * feet_to_meter;      % Convert from feet to meters


% Calculate the absolute layer edges
LayerEdgesAbs = cumsum([0, LayerThickness]) + gridTop;  % Cumulative sum of layer thicknesses, adjusted by gridTop

% Get the z-coordinate values of the nodes
z = G.nodes.coords(:, 3);  

% Adjust the z-coordinate for each layer
for k = 1:Nz
    % Find nodes in the current layer
    nodesInLayer = (z >= (k-1)/Nz) & (z <= k/Nz);
    
    % Local z-coordinate within the layer
    z_local = (z(nodesInLayer) - (k-1)/Nz) * Nz;
    
    % Adjust the z-coordinates of the nodes within the layer
    G.nodes.coords(nodesInLayer, 3) = LayerEdgesAbs(k) + z_local * (LayerEdgesAbs(k+1) - LayerEdgesAbs(k));
end

% Recompute the geometry of the grid
G = computeGeometry(G);

% Calculate total volume reservoir
totalVolume = sum(G.cells.volumes);

% Tampilkan hasil
fprintf('Total volume reservoir: %g m^3\n', totalVolume);

% == Reservoir Properties == %

% Porosity per layer (dimensionless)
poroPerLayer = [0.3, 0.3, 0.3];  

% Permeability in the I, J, and K directions (in milli-Darcy)
permI = [500, 50, 200] * milli * darcy;   % Permeability in I direction
permJ = permI;                            % Permeability in J direction (same as I)
permK = [300, 30, 50] * milli * darcy;    % Permeability in K direction

% Number of cells and cells per layer
numCells = G.cells.num;                   % Total number of cells in the grid
cellsPerLayer = Nx * Ny;                  % Number of cells per layer (Nx x Ny)

% Initialize arrays for permeability and porosity
perm = zeros(numCells, 3);                % Permeability in three directions (I, J, K)
poro = zeros(numCells, 1);                % Porosity for each cell

% Assign permeability and porosity values to each layer
for k = 1:Nz
    % Calculate the start and end indices for the current layer
    idxStart = (k-1) * cellsPerLayer + 1;
    idxEnd = k * cellsPerLayer;
    
    % Assign permeability values for each direction (I, J, K)
    perm(idxStart:idxEnd, 1) = permI(k);  % Permeability in the I direction
    perm(idxStart:idxEnd, 2) = permJ(k);  % Permeability in the J direction
    perm(idxStart:idxEnd, 3) = permK(k);  % Permeability in the K direction
    
    % Assign porosity for the current layer
    poro(idxStart:idxEnd) = poroPerLayer(k);
end

% Create rock properties based on the grid, permeability, and porosity
rock = makeRock(G, perm, poro);

% == Fluid Model Setup == %

% Set up the SPE1 fluid model
[~, ~, fluid, ~, ~] = setupSPE1();

% Show the bubble point pressure at reservoir
BubblePoinPressure = fluid.pb(226.1966570852417);
fprintf('PB = %.2f Psi\n', BubblePoinPressure .* pascal_to_psi);   % Output bubble point pressure
disp(' ');
disp('---------------------------');
disp(' ');

% Create a Generic Black Oil Model with the specified fluid properties
model = GenericBlackOilModel(G, rock, fluid, ...
                             'disgas', true, ...   % Enable dissolved gas
                             'vapoil', false, ...  % Disable vaporized oil
                             'water', true, ...    % Enable water
                             'oil', true, ...      % Enable oil
                             'gas', true);         % Enable gas

% Set the minimum pressure for the model
model.minimumPressure = 0;

% == Initial Pressure Values == %

% Define the pressure values for each layer (converted to Pascals)
p0 = [4783.4995; 4789.847; 4800] * psi_to_pascal;  % Pressure in each layer

% Define the number of layers
numLayers = length(p0);

% Calculate the number of cells per layer
numCellsPerLayer = G.cells.num / numLayers;

% Initialize an array to store the pressure values
pressure = zeros(G.cells.num, 1);

% Distribute the pressure values across the grid layers
pressure(1:numCellsPerLayer) = p0(1);                     % Top layer pressure
pressure(numCellsPerLayer+1:2*numCellsPerLayer) = p0(2);  % Middle layer pressure
pressure(2*numCellsPerLayer+1:end) = p0(3);               % Bottom layer pressure

% Initial saturations (s0) and other fluid properties
s0 = repmat([0.12, 0.88, 0.0], [G.cells.num, 1]);   % Initial saturations for each phase
rs0 = repmat(226.1966570852417, [G.cells.num, 1]);  % Initial solution gas oil ration
rv0 = 0;                                            % Initial gas density (dry gas)

% Create the state structure for the fluid properties
state = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', pressure);

% Get the initial state from the simulation results
stateInit = state;

% == Compute Initial In-Place Volumes == %

% Calculate Original Oil In Place (OOIP), Original Gas In Place (OGIP), and Original Water In Place (OWIP)
[ooip, ogip, owip] = computeOriginalInPlace(G, rock, fluid, stateInit);

% Display the computed values
fprintf('OOIP = %.2f STB\n', ooip);   % Output Original Oil In Place in Stock Tank Barrels
fprintf('OGIP = %.2e SCF\n', ogip);   % Output Original Gas In Place in Standard Cubic Feet
fprintf('OWIP = %.2f STB\n', owip);   % Output Original Water In Place in Stock Tank Barrels

disp(' ');
disp('---------------------------');
disp(' ');

% Clean up intermediate variables
clear k p0 s0 rs0;

% == Define Wells == %

% Initialize an empty array to store wells
W = [];

BHPINJ = 4500 * psi_to_pascal;          % BHP from psi to pascal
BHPPROD = 2000 * psi_to_pascal;         % BHP from psi to pascal

% Production Well: Prod
W = verticalWell(W, G, rock, ...
    35, 35, 1:3, ...                      % Location: row (I=1), column (J=1), full perforation (K=[])
    'name', 'Prod 1', ...               % Well name
    'compi', [0 1 0], ...               % Production of oil only
    'val', BHPPROD, ...                 % Oil rate in meter^3/s
    'type', 'bhp', ...                  % Control type: BHP
    'sign', -1);                        % Production well (Sign = -1)

% Injection Well: Inj
W = verticalWell(W, G, rock, ...
    1, 1, 1:3, ...                      % Location: row (I=25), column (J=25), perforation at K=1
    'name', 'Inj 1', ...                % Well name
    'compi', [1 0 0], ...               % Injection of water only
    'val', BHPINJ, ...                  % Injection rate (volume per second)
    'type', 'bhp', ...                  % Control type: Rate
    'sign', 1);                         % Injection well (Sign = 1)

W = verticalWell(W, G, rock, ...
    70, 1, 1:3, ...                     % Location: row (I=25), column (J=25), perforation at K=1
    'name', 'Inj 2', ...                % Well name
    'compi', [1 0 0], ...               % Injection of water only
    'val', BHPINJ, ...                  % Injection rate (volume per second)
    'type', 'bhp', ...                  % Control type: Rate
    'sign', 1);                         % Injection well (Sign = 1)

W = verticalWell(W, G, rock, ...
    1, 70, 1:3, ...                     % Location: row (I=25), column (J=25), perforation at K=1
    'name', 'Inj 3', ...                % Well name
    'compi', [1 0 0], ...               % Injection of water only
    'val', BHPINJ, ...                  % Injection rate (volume per second)
    'type', 'bhp', ...                  % Control type: Rate
    'sign', 1);                         % Injection well (Sign = 1)

W = verticalWell(W, G, rock, ...
    70, 70, 1:3, ...                    % Location: row (I=25), column (J=25), perforation at K=1
    'name', 'Inj 4', ...                % Well name
    'compi', [1 0 0], ...               % Injection of water only
    'val', BHPINJ, ...                  % Injection rate (volume per second)
    'type', 'bhp', ...                  % Control type: Rate
    'sign', 1);                         % Injection well (Sign = 1)

data = struct();
data.Porosity = rock.poro;                   % Porosity
data.Kx = rock.perm(:, 1) / (milli * darcy); % Permeability horizontal (Kx)
data.Ky = rock.perm(:, 2) / (milli * darcy); % Permeability horizontal (Ky)
data.Kz = rock.perm(:, 3) / (milli * darcy); % Permeability vertikal (Kz)
data.SaturationWater = state.s(:, 1);        % Water saturation (sW)
data.SaturationOil = state.s(:, 2);          % Oil saturation (sO)
data.SaturationGas = state.s(:, 3);          % Gas saturation (sG)
data.Pressure = state.pressure;              % Reservoir pressure

%figure;
%plotToolbar(G, data,'EdgeColor', 'k', 'LineWidth', 0.2);
%axis tight
%xlabel('X (m)');
%ylabel('Y (m)');
%zlabel('Z (m)');
%colorbar;

%hold on;
%plotWell(G, W, 'color', 'r', 'LineWidth', 2);
%hold off;

% Simulation time parameters
simTimeProd = 11 * year;     % Total simulation time for production (10 years)
simTimeInj = 9 * year;       % Total simulation time for injection
nstep   = 60 * day;          % Number of time steps (60 days)
refine  = 5;                 % Refinement factor for time steps

% Definisikan timestep untuk produksi dan injeksi menggunakan fungsi ramp-up
dt1 = rampupTimesteps(simTimeProd, nstep, refine); % Timestep untuk periode produksi tanpa injeksi
dt2 = rampupTimesteps(simTimeInj, nstep, refine);  % Timestep untuk periode produksi dengan injeksi

% Initialize the schedule structure for the simulation
schedule = struct();                % Create an empty schedule structure
schedule.control = struct([]);      % Initialize an empty control array for well operations
schedule.step = struct();           % Initialize an empty step array to store timesteps

% Combine timesteps for production and injection into one array
schedule.step.val = [dt1; dt2];     % Concatenate timesteps from production and injection phases

% Duplicate the well list for two different control scenarios
W1 = W;  % Scenario 1: Production wells active, injection wells inactive
W2 = W;  % Scenario 2: Both production and injection wells active

% Deactivate injection wells during the production-only period
for i = 1:numel(W1)
    if startsWith(W1(i).name, 'Inj')
        W1(i).val = 0;             % Set injection rate to zero
        W1(i).type = 'rate';       % Use rate control mode
        W1(i).compi = [1, 0, 0];   % Set composition to water phase only
    end
end

% Reactivate injection wells during the injection period
for i = 1:numel(W2)
    if startsWith(W2(i).name, 'Inj') 
        W2(i).val = BHPINJ;        % Set injection rate to predefined RATEWATER
        W2(i).type = 'bhp';        % Use rate control mode
        W2(i).compi = [1, 0, 0];   % Set composition to water phase only
    end
end

% Create control structures for each phase
schedule.control = repmat(struct('W', []), 2, 1); % Initialize two control steps
schedule.control(1).W = W1;                       % Control 1: Only production wells active
schedule.control(2).W = W2;                       % Control 2: Both production and injection wells active

% Map control indices to each timestep
schedule.step.control = [ ...
    ones(numel(dt1), 1);         % Assign control 1 to all production timesteps
    2 * ones(numel(dt2), 1)      % Assign control 2 to all injection timesteps
];

% Create the schedule using the calculated time steps
%schedule = simpleSchedule(dt1, 'W', W);

% == Plotting and Simulation Setup == %

% Define the function to plot after each simulation step
fn = getPlotAfterStep(state, model, schedule, 'plotWell', false, 'plotReservoir', false);

% Set up the nonlinear solver with relaxation enabled
nls = NonLinearSolver('useRelaxation', true);

% Simulate the schedule with the given state, model, and schedule
disp('Running for simulation');
disp(' ');
[wellSols, states, report] = simulateScheduleAD(state, model, schedule, ...
    'nonlinearsolver', nls, ...           % Use nonlinear solver
    'afterStepFn', fn);                   % Plot after each step

% Extract reservoir pressure from each state in the simulation
pRes = cellfun(@(x) x.pressure, states, 'UniformOutput', false);

% Initialize array to store average pressure for each time step
avgpRes = zeros(numel(states), 1);

% Loop through all time steps to calculate average reservoir pressure
for t = 1:numel(states)
    avgpRes(t) = mean(pRes{t});  % Compute the average pressure at time step t
end

% Find index injection wells (sign = 1)
wellINJ = find([wellSols{1}.sign] == 1);

% Find index production wells (sign = -1)
wellPROD = find([wellSols{1}.sign] == -1);

% Convert cumulative sum of schedule steps to years
T = convertTo(cumsum(schedule.step.val), year);

% Convert well solutions to vectors
[qWs, qOs, qGs, bhp] = wellSolToVector(wellSols);

%figure;
%plotToolbar(G, states); axis tight
%xlabel('X (m)');
%ylabel('Y (m)');
%zlabel('Z (m)');
%view(50, 50)
%colorbar;

%plotWellSols(wellSols, report.ReservoirTime)

%figure;
%plot(T, avgpRes .* pascal_to_psi, 'Color', 'k', 'LineStyle','-.', 'LineWidth', 1.5);
%xlabel('Time (Year)');
%ylabel('Reservoir Pressure (Psi)');
%set(gca, 'YColor', 'k');

%legend({'Reservoir Pressure (Psi)'});
%grid on

profile viewer
profile off