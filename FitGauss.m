function Gauss_fit_coeff = FitGauss(WL_start_index, WL_stop_index, t_stop_ind, Data)
    time = Data(2:end,1);
    WL = Data(1,2:end);
    time_range = time(1:t_stop_ind);
    WL_range = WL(WL_start_index:WL_stop_index);
    %% Estimate t0
    Est_t0 = zeros(length(WL_range),1); %intialize t0 estimation vector
    Gauss_fit_coeff = zeros(length(WL_range), 3); %Intialize fit coefficients matrix
    signal_moving_avg = zeros(t_stop_ind,1); % initialize moving avg signal for t0 est
    est_FWHM = input("Enter an estimate value for the CA FWHM: ");
    for i = WL_start_index:1:WL_stop_index
        signal = Data(2:t_stop_ind+1, i);
        %% moving avg for t0 estimation (avoids picking outliers)
        for k = 1:3
            signal_moving_avg(k,1) = signal(k);
        end
        for k = length(signal)-3:length(signal)
            signal_moving_avg(k,1) = signal(k);
        end
        for j = 4:length(signal)-3
            signal_moving_avg(j,1) = mean(signal(j-3:j+3)); %mean over seven point (can be changed)
        end
        %% estimate t0 from maximum signal and kick out outliers
        [max_row, ~]= find(signal_moving_avg==max(signal_moving_avg)); %find position of max value in moving average signal
        in = i-(WL_start_index)+1;
        Est_t0(in) = time_range(max_row); %use the max postion as estimated t0
        if in ~= 1 & abs(Est_t0(in,1)-Est_t0(in-1,1))> 0.4 %kick out outliers in t0 estimation
            Est_t0(in) = Est_t0(in-1);
        end
        [xData, yData] = prepareCurveData( time_range, signal );
         %% Set CA fit function for t0-estimation. gaussian fits for my data t0 better
        ft_G0 = fittype( 'A0* exp(-4*log(2)*(t-t0)^2/FWHM^2)', 'independent', 't', 'dependent', 'y', 'coefficients', {'FWHM', 'A0', 't0'} ); %Guassian fit func
    
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' ); 
        %include into opts for constraints:
        %'lower',[0 0 -inf], 'upper', [tau_lim inf inf]
        %% estimate start par
        est_t0 = Est_t0(i-(WL_start_index)+1);
        est_A0 = max_row;
        opts.StartPoint = [est_FWHM est_A0, est_t0];
        
        %% Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft_G0, opts);
        %% save fit coeff
        Gauss_fit_coeff(i-(WL_start_index)+1,:) = [fitresult.FWHM, fitresult. A0, fitresult.t0];
    end
end 