function [y, AR_order, MA_order] = getData(settings)
switch settings.useData
    case 0
        %Generate AR
        sampleSize         =   100
        sigma_sample        =   1.5;
        AR_coefficients     =   [0; 0; -0.75];
        MA_coefficients     =   [1; -1.5; 0.5625];
        
        %         AR_coefficients     =   [0.9; -0.23; 0.3]
        %         MA_coefficients     =   [1; -0.4; 0.2]
        %         MA_coefficients     =   [1; -0.4; 0.6; -0.5]
        %         MA_coefficients     =   [1]
        
        %         disp(['Roots AR Terms: ']);
        %         abs(roots([1; -AR_coefficients]'))
        %         disp(['Roots MA Terms: ']);
        %         abs(roots([MA_coefficients]'))
        AR_order = size(AR_coefficients,1);
        MA_order = size(MA_coefficients,1)-1;
        y = zeros(sampleSize+1000,1);
        epsilon = normrnd(zeros(sampleSize+1000+ max(MA_order,AR_order)+1,1),sigma_sample);
        
        for AR_gen_step = (max(AR_order,MA_order)+1):sampleSize+1000+max(AR_order,MA_order)+1
            y(AR_gen_step) = transpose(y(AR_gen_step - AR_order:AR_gen_step - 1,1))*AR_coefficients(end:-1:1) +...
                transpose(MA_coefficients)*epsilon(AR_gen_step:-1:AR_gen_step-MA_order);
        end;
        y = y(end-sampleSize+1:end);
        disp(['Length Data Vector: ' num2str(length(y)) ]);
    case 1
        %Generate Pure MA
        sampleSize         =   200;
        sigma_sample        =   0.8;
        MA_coefficients     =   [1]
        %         MA_pacs = [
        %    -0.8049;
        %    -0.4430;
        %     0.6]
        
        %         MA_coefficients =  [1; -getARParametersFromPACs(MA_pacs, 3)]
        
        disp(['Roots MA Terms: ']);
        abs(roots([MA_coefficients]'))
        MA_order = size(MA_coefficients,1)-1;
        AR_order = 0;
        y = zeros(sampleSize,1);
        epsilon = normrnd(zeros(sampleSize+MA_order,1),sigma_sample);
        
        for gen_step = 1:sampleSize
            y(gen_step) = transpose(MA_coefficients)*epsilon(gen_step+ MA_order:-1:gen_step);
        end;
        
    case 2 %Use Log Deviation
        AR_order = -1/0;
        MA_order = -1/0;
        [ y ] = getGDPData(1);
    case 3 %Use Log growth
        AR_order = -1/0;
        MA_order = -1/0;
        [ y ] = getGDPData(2);
    case 4 %Use Log GDP minus linear trend
        [ y ] = getGDPData(3);
        AR_order = -1/0;
        MA_order = -1/0;
    case 5 %Use synthetic data from Uhlig; Two-sided HP Filter
        load('OutputUhligExampl1_AR(1) Driving Force_Two-sided HP.mat');
        [ y ] = output(end-250:end);
        AR_order = -1/0;
        MA_order = -1/0;
    case 6 %Use US Cyclical GDP around 2-sided HP Trend
        load('GDPPC_twoSided_HP.mat');
        [ y ] = ycycle;
        AR_order = -1/0;
        MA_order = -1/0;
    case 7 %Use synthetic data Uhlig ARMA(3,3)
        load('Data RBC Model AR 0.9 -0.23 0.3 MA -0.4 0.6 -0.5 std(z) 0_712 Shock 250 Obs Synthetic.mat');
        y = transpose(y);
        AR_order = -1/0;
        MA_order = -1/0;
    case 8 %Use synthetic data from Uhlig; One-sided HP Filter
        load('OutputUhligExampl1_AR(1) Driving Force_One-sided HP.mat');
        AR_order = -1/0;
        MA_order = -1/0;
    case 9 %Use synthetic data from Uhlig; Two-sided HP Filter
        load('OutputUhligExampl1_AR(1) Driving Force_Two-sided HP.mat');
        AR_order = -1/0;
        MA_order = -1/0;
    case 10 %Use synthetic data from Uhlig; Raw
        load('OutputUhligExampl1_AR(1)DrivingForce_Raw.mat');
        y = data_raw;
        AR_order = -1/0;
        MA_order = -1/0;
    case 11 %Use US Cyclical GDP around 1-sided HP Trend
        load('GDPPC_oneSided_HP.mat');
        y = ycycle_onesided;
        AR_order = -1/0;
        MA_order = -1/0;
end;

end