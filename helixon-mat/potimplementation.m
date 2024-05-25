%% section for data handling

% Load pressure and DOA data for locations 2 and 4
pressure_files = {'SRIR40k/SRIR40k/ground_1_0_p.csv', 'SRIR40k/SRIR40k/ground_15_0_p.csv'};
doa_files = {'SRIR40k/SRIR40k/ground_1_0_doa.csv', 'SRIR40k/SRIR40k/ground_15_0_doa.csv'};

p_w_pos = {};
p_w_neg = {};
doa_data = {};
p_original = {};
for i = 1:length(pressure_files)
    df_pressure = csvread(pressure_files{i});
    df_pressure = df_pressure(401:20000); % Adjust indices according to your data
    disp("PRESSURE SHAPE: ");
    disp(size(df_pressure));
    disp(df_pressure(1));
    p_original{i} = df_pressure;


    
    df_doa = csvread(doa_files{i});
    df_doa = df_doa(401:20000,:); % Adjust indices according to your data
    disp("DOA SHAPE: ");
    disp(size(df_doa));
    disp(df_doa(1,:));


    % Splitting positive and negative pressure
    [p_pos, p_neg] = split_pos_neg(df_pressure);
    p_w_pos{i} = p_pos;
    p_w_neg{i} = p_neg;
    doa_data{i} = df_doa;
end

doa_data{1} = doa_data{1}(1:1000,:);
doa_data{2} = doa_data{2}(1:1000,:);
p_w_pos{1} = p_w_pos{1}(1:1000);
p_w_pos{2} = p_w_pos{2}(1:1000);
p_w_neg{1} = p_w_neg{1}(1:1000);
p_w_neg{2} = p_w_neg{2}(1:1000);

% Create point clouds for positive and negative pressures
PC1_pos = struct('pos', doa_data{1}, 'mass', p_w_pos{1}, 'n', length(p_w_pos{1}));
PC2_pos = struct('pos', doa_data{2}, 'mass', p_w_pos{2}, 'n', length(p_w_pos{2}));

PC1_neg = struct('pos', doa_data{1}, 'mass', p_w_neg{1}, 'n', length(p_w_neg{1}));
PC2_neg = struct('pos', doa_data{2}, 'mass', p_w_neg{2}, 'n', length(p_w_neg{2}));

PC1_original = struct('pos', doa_data{1}, 'mass', p_original{1}(1:1000), 'n', length(p_original{1}(1:1000)));
PC2_original = struct('pos', doa_data{2}, 'mass', p_original{2}(1:1000), 'n', length(p_original{2}(1:1000)));

disp("CHECK IF DIMENSIONS ARE RIGHT: ");
disp(size(PC2_neg.pos));
disp(size(PC2_pos.mass));
disp(PC2_pos.n);
%%  Cost Matrix
mu = 0.5;
C_euclidean = pdist2(PC1_pos.pos, PC2_pos.pos,"squaredeuclidean");
%save("potCe.mat",'C_eu', '-mat');
sign_mask = (PC1_original.mass.*PC2_original.mass)<0;
%disp(size(sign_mask));
pressure_diff = abs(PC1_original.mass-PC2_original.mass);
%disp(size(pressure_diff));
penalty = mu * pressure_diff .* sign_mask;
%disp(size(penalty));
C = C_euclidean+penalty;


%% POT
% Set parameters for VSOT
distEx = 1.20; % Example value, adjust as needed
distTol = 0.01; % Example value, adjust as needed

pot_pos = VSOT(PC1_pos, PC2_pos, C, distEx, distTol);

%%
pot_neg = VSOT(PC1_neg, PC2_neg, C, distEx, distTol);

%% 
% Interpolation with k = 0.5
k = 0.5;
PCk_pos_better = pot_pos.interpPC(k);
PCk_neg_better = pot_neg.interpPC(k);
%% 
writematrix(PCk_pos_better.mass, 'potMatMu5_8_1-15_p_pos.csv');
writematrix(PCk_pos_better.pos, 'potMatMu5_8_1-15_doa_pos.csv');
writematrix(PCk_pos_better.mass, 'potMatMu5_8_1-15_p_neg.csv');
writematrix(PCk_pos_better.pos, 'potMatMu5_8_1-15_doa_neg.csv');

%% 
all_doas = [PCk_pos_better.pos;PCk_neg_better.pos];
all_pressures = [PCk_pos_better.mass; (PCk_neg_better.mass)*-1];

%% 
c = 343;
timeshifts = vecnorm(all_doas,2,2)/c ;
time = linspace(0,1,100000);
ir = zeros(size(time));

for i = 1:length(all_pressures)
    index = round(timeshifts(i)*length(time))+1;
    ir(index) = ir(index)+all_pressures(i);
end

figure;
plot(time,ir);
title("IR");
grid on;

