% Created by Yapeng Wu on 2024/10/29
% matlab code
% To calculate a series of intermediate parameters (fAPAR£¬APAR£¬GPP£¬AGPP)
% and the final AGB with two strategies in the UAV-LUE method
% clc
% clear

% To read data
% Please replace it with your input file path
up_data = readtable('F:\LUEtoAGB_core_code_test_data\ACpar_data_agronomic_test.csv');
lai_agb_h_data = readtable('F:\LUEtoAGB_core_code_test_data\LAI_AGB_heading_test.csv');

% To ensure that the plot_Num column is in string format
up_data.plot_Num = string(up_data.plot_Num);
lai_agb_h_data.plot_Num = string(lai_agb_h_data.plot_Num);

% To initialize the result array
num_rows = height(up_data);
k_values = zeros(num_rows, 1); 
fAPAR = zeros(num_rows, 1);
APAR = zeros(num_rows, 1);
GPP = zeros(num_rows, 1); 
% stage-skipping,to calculate the AGPP from the heading stage to the current stage
AGPP = NaN(num_rows, 1); 
% stage-progressive£¬to calculate the AGPP of two adjacent stages (denoted as AAGPP)
AAGPP = NaN(num_rows, 1);    
GPP_AGB = NaN(num_rows, 1);   % stage-skipping
AGPP_AGB = NaN(num_rows, 1);  % stage-progressive

for i = 1:num_rows
    RowSpacing = up_data.RowSpacing(i);
    N = up_data.N{i};
    prosail_LAI = up_data.bp_prosail_LAI(i);
    Stage = up_data.Stage{i};
    plot_Num = up_data.plot_Num(i);
    Label_exp = up_data.Label_exp{i}; 
    
    % To obtain the calculated data
    PAR = up_data.PAR(i);
    PMIN = up_data.PMIN(i);
    APMIN = up_data.APMIN(i);
    
    % To determine the value of k
    if RowSpacing > 25
        k = 0.5;
    elseif RowSpacing <= 25 && strcmp(N, 'N0')
        k = 0.5;
    else
        k = 0.65;
    end
    k_values(i) = k;
    
    % To calculate fAPAR
    fAPAR(i) = 0.95 * (1 - exp(-k * prosail_LAI));
    
    % To calculate APAR
    APAR(i) = PAR * fAPAR(i);
    
    % To calculate GPP
    GPP(i) = fAPAR(i) * 2.85 * PMIN;
  
    % To calculate AGPP
    if ~strcmp(Stage, 'heading')
        % To find the corresponding bp_prosail_LAI value under the conditions
        % where the Label_exp and plot_Num are identical and the Stage column
        % is "heading", and record it as H_prosail_LAI.
        idx_heading = strcmp(up_data.Label_exp, Label_exp) & ...
                      strcmp(up_data.plot_Num, plot_Num) & ...
                      strcmp(up_data.Stage, 'heading');
        H_prosail_LAI = up_data.bp_prosail_LAI(idx_heading);
      
        if ~isempty(H_prosail_LAI)
            AGPP(i) = GPP(i) + APMIN * 2.85 * (0.95 * (1 - exp(-k * (H_prosail_LAI * 0.5 + prosail_LAI * 0.5))));
        end
    end
    
    % stage-progressive
    if strcmp(Stage, 'heading')
        AAGPP(i) = NaN;
    elseif i > 1 
        previous_idx = find(strcmp(up_data.Label_exp, Label_exp) & ...
                            strcmp(up_data.plot_Num, plot_Num) & ...
                            up_data.Date < up_data.Date(i));  

        if ~isempty(previous_idx)
            % To obtain bp_prosail_LAI for the previous stage
            previous_prosail_LAI = up_data.bp_prosail_LAI(previous_idx(end));           
            AAGPP(i) = GPP(i) + APMIN * 2.85 * (0.95 * (1 - exp(-k * (previous_prosail_LAI * 0.5 + prosail_LAI * 0.5))));
        else           
            AAGPP(i) = NaN;
        end
    end
        
    % To calculate AGB (denoted as GPP_AGB) with the stage-skipping strategy
    idx_label_exp_match = strcmp(lai_agb_h_data.Label_exp, Label_exp) & lai_agb_h_data.plot_Num == plot_Num;
    if any(idx_label_exp_match)
        if strcmp(Stage, 'heading')
            % To extract the LAI_AGB_H values from the file LAI_AGB_heading.csv under the condition that Label_exp and plot_Num are identical.
            GPP_AGB(i) = lai_agb_h_data.LAI_AGB_H(idx_label_exp_match);
        else
            heading_idx = strcmp(up_data.Label_exp, Label_exp) & ...
                      strcmp(up_data.plot_Num, plot_Num) & ...
                      strcmp(up_data.Stage, 'heading');
            H_GPP_AGB = GPP_AGB(heading_idx);
            
            if any(~isnan(H_GPP_AGB)) 
                GPP_AGB(i) = H_GPP_AGB + AGPP(i) * 0.533 * 2.22 * (1 / (1 + 0.1)) / 100;
            end
        end
    end
    
    % Using the stage-progressive strategy to calculate AGB (denoted as AGPP_AGB)        
    % For the first stage (heading), AGPP_AGB is calculated using the same method as GPP_AGB.
    idx_label_exp_match = strcmp(lai_agb_h_data.Label_exp, Label_exp) & lai_agb_h_data.plot_Num == plot_Num;
    if any(idx_label_exp_match)
       if strcmp(Stage, 'heading')            
           AGPP_AGB(i) = lai_agb_h_data.LAI_AGB_H(idx_label_exp_match);
       elseif i > 1                      
            % For subsequent dates, look up the AGPP_AGB value from the previous date under the same Label_exp and plot_Num.
            previous_idx = find(strcmp(up_data.Label_exp, Label_exp) & ...
                                strcmp(up_data.plot_Num, plot_Num) & ...
                                up_data.Date < up_data.Date(i));  

            if ~isempty(previous_idx)
                % To obtain the AGPP_AGB values for the previous date
                previous_AGPP_AGB = AGPP_AGB(previous_idx(end));

                % To calculate the AGPP_AGB values for the current date
                AGPP_AGB(i) = previous_AGPP_AGB + AAGPP(i) * 0.533 * 2.22 * (1 / (1 + 0.1)) / 100;
            else               
                AGPP_AGB(i) = NaN;
            end
        end             
    end
end

% To add the results to the table
up_data.k = k_values; 
up_data.fAPAR = fAPAR;
up_data.APAR = APAR;
up_data.GPP = GPP;
up_data.AGPP = AGPP;
up_data.AAGPP = AAGPP;
up_data.GPP_AGB = GPP_AGB;
up_data.AGPP_AGB = AGPP_AGB;

% To save the results
% Please replace it with your output file path
writetable(up_data, 'F:\LUEtoAGB_core_code_test_data\GPP_AGB_data_agronomic_test.csv'); 
