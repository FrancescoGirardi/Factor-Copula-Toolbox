

% This code is used to estimate a factor copula with GAS dynamics.

% Innovations of the latent variables can be choosen from the Normal or
% Student's t distribution.

% The distribution of the factors can be freely chosen among the Normal,
% Student's t and the Skew Student's t.

%The following code is an adaptation of the one provided by Oh and Patton
%in their seminal paper on factor copulas with GAS dynamics (2013).

%We allow for different distributions and improve clarity.

%We note that the procedure is extremely slow, even with few variables and
%few groups.

%Total runtime for a portfolio of 44 stocks, as the one presented by
%Ziegelmann (2016), can exceed 10 hours.  We strongly advise anyone trying
%to understand the code to limit the number of variables and the total
%number of iterations used in the numeric procedure. This, of course, come
%with a cost in terms of precision.

%It's also possible to choose the dependence structure in this same file.
%That is: heterogeneous, block or homogeneous dependence.

% Choose the distribution of the factor. 
% Possibilites: 'Normal', 't', 'skew_t'
factor_dist = 't';

% Choose the distribution of factor innovations. 
% Possibilities: 'Normal' or 't'
innov_dist = 'Normal';

%Choose between a Static or Dynamic Model
model = 'Dynamic';

% Choose the input and output folders
input_path = '/Users/victorcozer/Documents/MATLAB/FRM Project/';  
output_path = '/Users/victorcozer/Documents/MATLAB/FRM Project/'; 

%%% Load Data;  (T by N) matrix. This is the data used to fit the copula.
%%% Each column contains a time series.
% The marginals should be uniformilly distributed; i.e. all elements should
% be in [0,1].
pit = readtable([input_path, 'pit.csv']);
copula_data = pit(:,2:45); %select the relevant columns in your data.
copula_data = copula_data{:,:};

% Choose dependence structure: Homogeneous, Block, Heterogeneous 
dependence = 'Heterogeneous';

% Only for Block dependence, for each variable, assign its corresponding group. 
% For example, if you have 3 stocks and the first and last one belong to the same
% group, use group_code = [2 1 2]. 

%group_code = [2, 1, 2];  %%%UN-COMMENT THIS IN CASE OF BLOCK DEPENDENCE.

% Choose number of iterations, for precise results we recommend > 1000.
max_iter = 5;


%%% Input check for dependence
dep_check = strcmp(dependence, 'Homogeneous') | strcmp(dependence, 'Heterogeneous') | strcmp(dependence, 'Block');
if dep_check == false
    error('Invalid dependence structure')
end


if strcmp(dependence, 'Homogeneous')
    group_code = ones(size(copula_data, 2), 1);
    GASFacCop_G(copula_data, factor_dist, innov_dist, model, max_iter, group_code)
end

if strcmp(dependence, 'Block')
    GASFacCop_G(copula_data, factor_dist, innov_dist, model, max_iter, group_code)
end

if strcmp(dependence, 'Heterogeneous')
    GASFacCop_N(copula_data, factor_dist, innov_dist, model, max_iter)
end
