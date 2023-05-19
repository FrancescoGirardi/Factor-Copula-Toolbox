

% This code is used to estimate a factor copula with GAS dynamics.

% Innovations of the latent variables can be choosen from the Normal or
% Student's t distribution.

% The distribution of the factors can be freely chosen among the Normal,
% Student's t and the Skew Student's t.

%It's also possible to choose the dependence structure in this same file.
%That is: heterogeneous, block or homogeneous dependence.

%The following code is an adaptation of the one provided by Oh and Patton
%in their seminal paper on factor copulas with GAS dynamics (2013).

%We allow for different distributions, the possibility to use a static or dynamic model,
% we improve clarity and provide an interface for better usage.  
% If your computer has problems running the app, please use the file GASFacCop_NoGUI.m

%We note that the procedure is extremely slow, even with few variables and
%few groups. For this reason, at the cost of losing some precision, we
%allow for the possibility of choosing the number of iterations used in the
%numerical procedure. For better results, please use more than 1000
%iterations.

%Total runtime for a portfolio of 44 stocks, as the one presented by
%Ziegelmann (2016), can exceed 10 hours.  We strongly advise anyone trying
%to understand the code to limit the number of variables and the total
%number of iterations used in the numeric procedure. This, of course, come
%with a cost in terms of precision.


%In order to fit the factor copula follow this simple procedure:
%1. Input your data frame and select the columns you want to use.
%2. Run the code, a window will appear where you can choose all the
%different specifications of the model.
%3. Unless you choose "Block dependence", you're all set and the algorithm
%will start.
%4. If you choose "Block dependence", you are required to change the
%"group_code" array below, before running the code.
%5. Again, depending on the choices made, it may take several seconds just
%to see the result of the first few interations in the console. Be patient.
%6. Note: The program is done running when the variable "esti_time" is
%printed to the console. If this has not happenned, wait.
%6. At the end of the fitting process, you will be asked whether to simulate the
%fitted copula, as well as the number of simulations. If you want to
%simulate at another time, use the "Simulation.m" file provided.
%7. It takes a few seconds for the simulation to start, you will see in the
%console the progress. Be patient.

% Choose the input and output folders
input_path = '/Users/francesco_girardi/Desktop/FRM Project 2/';  
output_path = '/Users/francesco_girardi/Desktop/FRM Project 2/'; 

%%% Load Data;  (T by N) matrix. This is the data used to fit the copula.
%%% Each column contains a time series.
% The marginals should be uniformilly distributed; i.e. all elements should
% be in [0,1].
pit = readtable([input_path, 'pit.csv']);
copula_data = pit(:,2:45); %select the relevant columns in your data.
copula_data = copula_data{:,:};

% Only for Block dependence, for each variable, assign its corresponding group. 
% For example, if you have 3 stocks and the first and last one belong to the same
% group, use group_code = [2 1 2]. 

group_code = zeros(1, size(copula_data, 2)); %leave this intact
group_code = [1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 5 5 5 5 6 6 6 6 6 6 6 6 6 6];  %%%UN-COMMENT THIS IN CASE OF BLOCK DEPENDENCE.

% Create the main figure
fig = figure('Name', 'Select Options', 'MenuBar', 'none', 'ToolBar', 'none', 'NumberTitle', 'off', 'Position', [0 0 400 350], 'UserData', []);
movegui(fig, 'center');

% Create the popup menus
pm1 = createPopupMenu(fig, 'Factor distribution:', [50 280 300 20], {'Normal', 't', 'skew_t'});
pm2 = createPopupMenu(fig, 'Innovation distribution:', [50 230 300 20], {'Normal', 't'});
pm3 = createPopupMenu(fig, 'Model:', [50 180 300 20], {'Static', 'Dynamic'});
pm4 = createPopupMenu(fig, 'Dependence structure:', [50 130 300 20], {'Homogeneous', 'Block', 'Heterogeneous'});

% Create an input field for max_iter
uicontrol(fig, 'Style', 'text', 'String', 'Number of iterations:', 'FontSize', 14, 'Position', [50 80 200 20]);
max_iter_input = uicontrol(fig, 'Style', 'edit', 'String', '5', 'FontSize', 12, 'Position', [260 80 100 30]);

% Add a submit button
submit_button = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Submit', 'FontSize', 14, 'Position', [150 30 100 30]);

% Callback function for the submit button
submit_button.Callback = @(~,~) submitCallback(pm1, pm2, pm3, pm4, max_iter_input, fig);

% Wait for the user to submit
uiwait(fig);

% Retrieve the output values
if isvalid(fig)  % Check if figure still exists
    out = get(fig, 'UserData');
    % Delete the figure
    delete(fig);
else
    out = {'Normal', 'Normal', 'Static', 'Homogeneous', 5}; % Default values
end
 
factor_dist = out{1};
innov_dist = out{2};
model = out{3};
dependence = out{4};
max_iter = out{5};

if strcmp(dependence, 'Homogeneous')
    group_code = ones(1, size(copula_data, 2));
    [theta_NMLE, lam] = GASFacCop_G(copula_data, factor_dist, innov_dist, model, max_iter, group_code);
end

if strcmp(dependence, 'Block')
    [theta_NMLE, lam] = GASFacCop_G(copula_data, factor_dist, innov_dist, model, max_iter, group_code);
end

if strcmp(dependence, 'Heterogeneous')
    group_code = 1:size(copula_data, 2);
    [theta_NMLE, lam] = GASFacCop_N(copula_data, factor_dist, innov_dist, model, max_iter);
end

% Create the main figure for simulation
fig = figure('Name', 'Simulation Menu', 'Position', [100 100 400 200]);
movegui(fig, 'center');

% Create a message panel with text
msgPanel = uipanel(fig, 'Position', [0.05 0.25 0.9 0.5]);
msgText = uicontrol(msgPanel, 'Style', 'text', 'String', 'Would you like to simulate the fitted copula?', 'FontSize', 14, 'Position', [20 30 360 50]);

% Add "Yes" button
yes_button = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Yes', 'FontSize', 14, 'Position', [85 40 100 30]);
yes_button.Callback = @(~,~) yesCallback(fig, group_code, theta_NMLE, lam, model);

% Add "No" button
no_button = uicontrol(fig, 'Style', 'pushbutton', 'String', 'No', 'FontSize', 14, 'Position', [215 40 100 30]);
no_button.Callback = @(~,~) noCallback(fig);

% Function for "Yes" button callback
function yesCallback(fig, group_code, theta_NMLE, lam, model)
    % Close the main figure
    close(fig);
    
    % Create a new dialog for input
    dlg = dialog('Name', 'Simulation Input', 'Position', [100 100 400 200]);
    movegui(dlg, 'center');

    % Create an input field for number of days
    uicontrol(dlg, 'Style', 'text', 'String', 'Number of days:', 'FontSize', 14, 'Position', [50 130 150 20]);
    days_input = uicontrol(dlg, 'Style', 'edit', 'String', '120', 'FontSize', 12, 'Position', [220 130 150 30]);

    % Create an input field for number of simulations
    uicontrol(dlg, 'Style', 'text', 'String', 'Number of simulations:', 'FontSize', 14, 'Position', [50 80 150 20]);
    sim_input = uicontrol(dlg, 'Style', 'edit', 'String', '300', 'FontSize', 12, 'Position', [220 80 150 30]);

    % Add a submit button
    submit_button = uicontrol(dlg, 'Style', 'pushbutton', 'String', 'Submit', 'FontSize', 14, 'Position', [160 20 100 30]);

    % Callback function for the submit button
    submit_button.Callback = @(~,~) submitCallback2(dlg, days_input, sim_input, group_code, theta_NMLE, lam, model);

    % Wait for the user to submit
    %uiwait(dlg);
end

% Function for "No" button callback
function noCallback(fig)
    % Close the main figure
    close(fig);
    
    % Display a message
    fprintf('Simulation not selected.\n');
end

% The inner submit callback function
function submitCallback2(dlg, days_input, sim_input, group_code, theta_NMLE, lam, model)
    
    % Get the input values
    T = str2double(days_input.String);
    nsim = str2double(sim_input.String);

    close(dlg),
    %set(dlg, 'Visible', 'off')

    dim = length(group_code);

    if strcmp(model, 'Static')
        theta_hat = [1/theta_NMLE(end-2) theta_NMLE(end) 1/theta_NMLE(end-1)];
        [U] = generate_copula_static(nsim, T, dim, theta_hat, lam(length(lam), :), group_code);
        
        csvwrite('U.csv',U);
    end

    if strcmp(model, 'Dynamic')
        theta_hat = [1/theta_NMLE(end-2) theta_NMLE(end) 1/theta_NMLE(end-1) theta_NMLE(end-4) theta_NMLE(end-3)];
        omega_hat = theta_NMLE(1:(max(group_code)));
        [U] = generate_copula_dynamic(nsim, T, dim, theta_hat, omega_hat, lam(length(lam), :), group_code);
        
        csvwrite('U.csv',U);
    end

    % Display the input values (or use them in your code)
    fprintf('Number of days: %d\n', T);
    fprintf('Number of simulations: %d\n', nsim);
end





% Function to create a popup menu
function pm = createPopupMenu(parent, label, position, options)
    uicontrol(parent, 'Style', 'text', 'String', label, 'FontSize', 14, 'Position', [position(1) position(2)+25 200 20]);
    pm = uicontrol(parent, 'Style', 'popupmenu', 'String', options, 'FontSize', 12, 'Position', [position(1) position(2) position(3) position(4)]);
end

% The submit callback function
function submitCallback(pm1, pm2, pm3, pm4, max_iter_input, fig)
    % Get the selected options
    factor_dist = pm1.String{pm1.Value};
    innov_dist = pm2.String{pm2.Value};
    model = pm3.String{pm3.Value};
    dependence = pm4.String{pm4.Value};
    max_iter = str2double(max_iter_input.String);

    % Store the selected options in a cell array
    out = {factor_dist, innov_dist, model, dependence, max_iter};

    % Save the output in the figure's 'UserData' property
    set(fig, 'UserData', out);

    % Display the selected options (or use them in your code)
    fprintf('factor_dist: %s\n', factor_dist);
    fprintf('innov_dist: %s\n', innov_dist);
    fprintf('model: %s\n', model);
    fprintf('dependence: %s\n', dependence);
    fprintf('max_iter: %d\n', max_iter);

    % Show a message if Block dependence is selected
    if strcmp(dependence, 'Block')
    % Create a custom dialog
    dlg = dialog('Name', 'Block Dependence', 'Position', [100 100 400 200]);
    movegui(dlg, 'center');
    
    % Create a message panel with text
    msgPanel = uipanel(dlg, 'Position', [0 0.2 1 0.8]);
    msgText = uicontrol(msgPanel, 'Style', 'text', 'String', ...
        {'You have selected Block dependence. We assume that you have changed the variable "group_code" with the group of each variable.'}, ...
        'HorizontalAlignment', 'left', 'Units', 'normalized', 'Position', [0 0 1 1]);
    set(msgText, 'FontSize', 14);
    
    % Create an OK button
    okButton = uicontrol(dlg, 'Style', 'pushbutton', 'String', 'OK', ...
        'Units', 'normalized', 'Position', [0.4 0.05 0.2 0.1], ...
        'Callback', @(~,~) delete(dlg));

    end

    % Clear variables and close the dialog
    delete(pm1);
    delete(pm2);
    delete(pm3);
    delete(pm4);
    delete(max_iter_input);

    set(fig, 'Visible', 'off');
    uiresume(fig)
end


