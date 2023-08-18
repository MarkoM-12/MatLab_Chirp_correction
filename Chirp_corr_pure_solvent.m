%% Enter filename
filename = 'CA_600nm_bkg_corr_cut';
Data = readmatrix(filename);
WL = Data(1,2:end);
time = Data(2:end,1);
%% Kick out duplicates in time steps if existent
for i = 2:length(time)
    diff = time(i,1)-time(i-1,1);
    if diff == 0
        Data = [Data(1:i,:); Data((i+2):end, :)];
        str = sprintf('Duplicate in line %d kicked out', i);
        disp(str)
    end
end

%% Choose time range
stop_time= input('Until which time delay was the data recorded in short time steps?  ');
[~,stop_ind] =  min(abs(time-stop_time));
time_range = time(1:stop_ind);
%plot to choose WL
f1 = figure(1);
clf(1)
h = pcolor(WL,time_range, Data(2:stop_ind+1, 2:end));
xlabel('wavelength [nm]')
ylabel('time [ps]')
set(h, 'EdgeColor', 'none');
colormap(jet(256)); 
colorbar;
caxis([-0.04, 0.04]); %limits for the colorbar
%% Choose wavelength range
WL_start = input('Enter start wavelength where artifact appears: ');
WL_stop = input('Enter end wavelength where artifact dissapears: ');
% tau_lim = input('set upper limit for FWHM: ');
[~, WL_start_index] = min(abs( WL-WL_start));
[~, WL_stop_index] = min(abs(WL-WL_stop));
WL_range = WL(WL_start_index:WL_stop_index);

%% Apply function in file FitGauss.m get t0 points for dispersion fit
Gauss_fit_coeff = FitGauss(WL_start_index, WL_stop_index, stop_ind, Data);
t0 = Gauss_fit_coeff(:,3);

%% correct outliers im fitted t0
for i =1:length(WL_range)
    if i ~= 1 & abs(t0(i)-t0(i-1))> 0.08
        t0(i) = t0(i-1);
    end
end

%% fit t0 with a reciprocal function
sqrt_ft = fittype('K.*wl.^(-1)+L.*wl.^(-1/2)+a', 'independent', 'wl', 'dependent', 't0', 'coefficients', {'K', 'L',  'a'});
opts = fitoptions( 'Method', 'NonlinearLeastSquares'); 
startK = -1;
startL = -1;

starta = -1;
[xData, yData] = prepareCurveData( WL_range, t0 );
for i = 1:10
    opts.StartPoint = [startK startL starta];
    [fitresult_t0, gof] = fit( xData, yData, sqrt_ft, opts );
    startK = fitresult_t0.K;
    startL = fitresult_t0.L;
    starta = fitresult_t0.a;
end
Fitted_t0 = extrapolateFit(fitresult_t0.K, fitresult_t0.L, fitresult_t0.a, WL);


%% plot dispersion fit and t0 values from CA fitting
f1 = figure(1);
clf(1)
h = pcolor(WL,time_range, Data(2:stop_ind+1, 2:end)); 
xlabel('wavelength [nm]')
ylabel('time [ps]')
set(h, 'EdgeColor', 'none');
colormap(jet(256)); 
colorbar;
caxis([-0.04, 0.04]); %limits for the colorbar
hold on
Z = zeros(length(WL),1);
for i =1:length(WL) %this is just to be able to see the stuff plotted on top the color plot
    Z(i) = 0.5;
end
Z2 = zeros(length(WL_range),1);
for i =1:length(WL_range)
    Z2(i) = 0.5;
end
plot3(WL_range, t0, Z2, 'r+') % Plot t0 s from Gaussian fit as red crosses
plot3(WL, Fitted_t0, Z, 'k-', 'Linewidth', 1) % Plot dispersion fit as black line


answer = input('Happy with the fit (Yes [Y] or No [N])? : ', 's');
%% Dispersion correcion

if answer == 'Y'
    [val,zero_idx]=min(abs(Fitted_t0)); %for cutting
    new_Data = zeros(size(Data)); %initialize
    new_time = time - Fitted_t0(zero_idx); % new time axis to which we have to bring the data
    new_Data(1,2:end) = WL;
    individual_time = zeros(length(time),length(WL)); % individual time initialization
    new_Data(2:end, 1) = new_time; % put the new_time part into new_DAta
    %% Interpolate data onto new time axis
    for i = 1:length(WL)
        individual_time(:, i) = time - Fitted_t0(i); % for each wavelength
        %interpolation using the Piecewise Cubic Hermite Interpolating Polynomial method:
        new_Data(2:end,i+1) = interp1(individual_time(:,i), Data(2:end,i+1), new_time, 'pchip'); 
    end
    %% Cutting the interpolated Data
    smallest_neg_time = min(Fitted_t0);
    cut_time = min(time)-smallest_neg_time;
    [val,cut_idx] =  min(abs(new_time-cut_time));
    new_time = new_time(cut_idx:end);
    Cut_Data = zeros(size(new_Data(cut_idx:end,:)));
    Cut_Data(1, 2:end) = WL;
    Cut_Data(2:end, :) = new_Data((cut_idx+1):end,:);
    %% Plotting Dispersion corrected Data
    [~,new_stop_ind] =  min(abs(new_time-stop_time));
    f2 = figure(2);
    clf(2)
    h = pcolor(WL(2:end),new_time(1:new_stop_ind), Cut_Data(2:new_stop_ind+1, 3:end));
    xlabel('wavelength [nm]')
    ylabel('time [ps]')
    set(h, 'EdgeColor', 'none');
    colormap(jet(256)); 
    colorbar;
    caxis([-0.04, 0.04]); %limits for the colorbar
    
    CA_substraction = input('Do you want to substract the coherent artifact (Yes [Y] or No[N])?: ', 's');
    
    %% CA-fitting of dispersion corrected Data
    if CA_substraction == 'Y'
        final_CA_coeff = zeros(length(WL), 6); %intialize
        fitted_CA = Cut_Data;
        fitted_CA(2:end,2:end) = zeros(size(Cut_Data(2:end, 2:end))); %intialize
        [~, t_start_index] = min(abs( new_time+0.1));
        [~, t_stop_index] = min(abs(new_time-0.1));
        for i =  1:length(WL)
            [xData_corr, yData_corr] = prepareCurveData( new_time(t_start_index:t_stop_index), Cut_Data(t_start_index+1:t_stop_index+1,i+1) ); %Data until stop time for fit
            % Function: combined Gaussian and derivative + sinusoidal model (4components)
            ft = fittype('cos((t-t0)^2)*A0*exp((-4*log(2)*(t-t0)^2)/tau^2)-A1*(8*log(2)/(tau^2))*(t-t0)*exp((-4*log(2)*(t-t0)^2)/tau^2)-A2*(8*log(2))/tau^2 *((8*(t-t0)^2 * log(2))/tau^2 - 1)*exp((-4*log(2)*(t-t0)^2)/tau^2)'...
                ,'independent', 't', 'dependent', 'y', 'coefficients', {'tau', 'A0', 'A1', 'A2', 't0'});
            opts = fitoptions( 'Method', 'NonlinearLeastSquares');
            % 'lower',  [0 0 -inf -inf], "upper", [tau_lim inf inf inf] 
            % take values from t0 Gaussian fitting for estimation, but if
            % outside of wavelength range take first or last value resp.:
            if i < WL_start_index 
                CA_est_FWHM = Gauss_fit_coeff(1,1);
                CA_est_A0 = Gauss_fit_coeff(1,3);
            elseif i > WL_stop_index
                CA_est_FWHM = Gauss_fit_coeff(end,1);
                CA_est_A0 = Gauss_fit_coeff(end,3);
            else
                CA_est_FWHM = Gauss_fit_coeff(i-WL_start_index+1,1);
                CA_est_A0 = Gauss_fit_coeff(i-WL_start_index+1,3);
            end
            CA_est_FWHM = 0.04;
            CA_est_A1 = CA_est_A0*10^(-1);
            CA_est_A2 = CA_est_A0*10^(-2);
            opts.StartPoint = [CA_est_FWHM CA_est_A0, CA_est_A1, CA_est_A2, 0]; %start points (time is zero)
            % Fit the data:
            [fitresult_final, gof] = fit( xData_corr, yData_corr, ft, opts);
            %write the fit coefficients into the initalized matrix:
            final_CA_coeff(i, :) = [WL(i), fitresult_final.tau, fitresult_final.A0, fitresult_final.A1, fitresult_final.A2, fitresult_final.t0];
            %Calculate the fitted CA by feeding the fit parameters into the
            %fit function CAfitting4(), this is for substraction:
            fitted_CA(2:end,i+1) = CAfitting4(new_time, fitresult_final.tau, fitresult_final.A0, fitresult_final.A1, fitresult_final.A2, fitresult_final.t0);
        end
        Corr_Data = Cut_Data - fitted_CA; %substract the CA
        %% Plot the CA corrected Data
        f3 = figure(3);
        clf(3)
        h = pcolor(WL(2:end),new_time(1:new_stop_ind), Corr_Data(2:new_stop_ind+1, 3:end));
        set(h, 'EdgeColor', 'none');
        colormap(jet(256)); 
        colorbar;
        caxis([-0.04, 0.04]); %limits for the colorbar
        
        f8 = figure(8);
        clf(8)
        hold on
        plot(new_time(t_start_index:t_stop_index), fitted_CA(t_start_index+1:t_stop_index+1, 51), '-r')
        plot(new_time(t_start_index:t_stop_index), fitted_CA(t_start_index+1:t_stop_index+1, 187), '-g')
        plot(new_time(t_start_index:t_stop_index), fitted_CA(t_start_index+1:t_stop_index+1, 265), '-b')
        
        plot(new_time(t_start_index:t_stop_index), Cut_Data(t_start_index+1:t_stop_index+1, 51), 'xr')
        plot(new_time(t_start_index:t_stop_index), Cut_Data(t_start_index+1:t_stop_index+1, 187), 'xg')
        plot(new_time(t_start_index:t_stop_index), Cut_Data(t_start_index+1:t_stop_index+1, 265), 'xb')
        legend('550 nm',' 620 nm', '660 nm')
        %For trace by trace analysis of fit and data use this
%         figure(5)
%         clf(5)
%         hold on
%         plot(new_time(69:end), Cut_Data(70:end, 50), '+r')
%         plot(new_time(69:end), fitted_CA(70:end, 50), '-k')
        %% Saving the CA fit to correc t other data
        saving_CA = input('Do you want to save this CA fit?: ', 's');
        if saving_CA == 'Y'
            fitted_CA(2:end,1) = new_time;
            fitted_CA(1,2:end) = WL;
            CA_filename = append(filename, '_fittedCA.dat');
            writematrix(final_CA_coeff, CA_filename, 'Delimiter', 'tab');
            disp(append('The CA fit parameters have been saved as: ', CA_filename));
        end
    end
end

saving_dispersion = input('Do you want to save the fitted dispersion to apply to other data?: ', 's');
if saving_dispersion == 'Y'
    dispersion = [WL;Fitted_t0];
    CA_t0 = [WL_range;transpose(t0)];
    disp_filename = append(filename, '_t0Fit.dat');
    writematrix(dispersion, disp_filename, 'Delimiter', 'tab');
    disp(append(disp_filename, ' has been saved.'));
end



if answer == 'N'
    disp('Change the fit parameters, or perform outlier correction.')
end


% plot(WL_range, transpose(t0))
% hold on
% plot(WL, Fitted_t0)


% Function of dispersion correction
function out = extrapolateFit(K, L, a, wl)
    out = K.*wl.^(-1)+L.*wl.^(-1/2)+a;
end

% Fuction of CA fit
function out2 = CAfitting4(t, tau, A0, A1, A2, t0)
    out2 = cos((t-t0).^2).*A0.*exp((-4.*log(2).*(t-t0).^2)./tau.^2)-A1.*(8.*log(2)./(tau.^2)).*(t-t0).*exp((-4.*log(2).*(t-t0).^2)./tau.^2)-A2.*(8.*log(2))./tau.^2 .*((8.*(t-t0).^2 .* log(2))./tau.^2 - 1).*exp((-4.*log(2).*(t-t0).^2)./tau.^2);
end