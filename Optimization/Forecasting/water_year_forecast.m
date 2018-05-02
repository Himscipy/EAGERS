function [solution, gen] = water_year_forecast(gen,buildings,cool_tower,subnet,options,date,hist_prof,prev_data,now_data,future_data)
global Plant TestData
solution = [];
hydroforecast = false;
for i = 1:1:length(gen)
    if strcmp(gen(i).Type,'Hydro Storage')
        hydroforecast = true;
    end
end
if isfield(Plant,'WYForecast') 
    solution = Plant.WYForecast{end};
    if solution.Timestamp(1)<=date && solution.Timestamp(end)>=date
        hydroforecast = false; %don't need to re-do
    end
    n_y = length(Plant.WYForecast)+1;
else 
    n_y = 1;
end
if hydroforecast
    %first create the yearly dispatch data 
    % i.e. run Dispatch loop with updated information
    %these will be used for set points in the actual dispatch
    hold_test_data = TestData;%save to avoid re-interpolating at simulation frequency
    hold_date = date;
    options.Horizon = 364*24; %Yearly Horizon
    options.Resolution = 7*24; %Week Resolution
    D = datevec(date);
    if D(2)<10
        year = D(1)-1;
    else
        year = D(1);
    end
    date = datenum([year 10 1 1 0 0]);
    forecast_time = date+[0;build_time_vector(options)/24];
    dt = forecast_time(2:end) - forecast_time(1:end-1);
    solution.Timestamp = forecast_time;
    Plant.WYForecast{n_y} = solution;
    TestData.RealTimeData = interpolate_data(TestData,options.Resolution*3600,0.00);%create test data at correct frequency
    op_mat_a = load_matrices(gen,buildings,cool_tower,subnet,options,'A',dt); %build quadratic programming matrices for FitA
    op_matb = load_matrices(gen,buildings,cool_tower,subnet,options,'B',dt);%build quadratic programming matrices for FitB
    one_step = load_matrices(gen,buildings,cool_tower,subnet,options,'B',[]);%build quadratic programming matrices for single time step
    if n_y == 1
        hydro_soc_init = zeros(1,length(subnet.Hydro.nodes));
        for n=1:1:length(subnet.Hydro.nodes)
            i = subnet.Hydro.Equipment{n};
            if strcmp(gen(i).Type,'Hydro Storage')
                gen(i).CurrentState(2) = gen(i).VariableStruct.StartWYstate*gen(i).QPform.Stor.UsableSize;
                hydro_soc_init(n) = gen(i).CurrentState(2);
            end 
        end 
        [data_t0,gen,~] = update_forecast(gen,buildings,cool_tower,subnet,options,date(1),hist_prof,prev_data,now_data,future_data);
        [gen,cool_tower] = automatic_ic(gen,buildings,cool_tower,subnet,date(1),one_step,options,data_t0);% set the initial conditions
    end
    [forecast,gen,buildings] = update_forecast(gen,buildings,cool_tower,subnet,options,forecast_time(2:end),hist_prof,prev_data,now_data,future_data);
    if n_y>1
        solution = dispatch_loop(gen,buildings,cool_tower,subnet,op_mat_a,op_matb,one_step,options,forecast_time,forecast,Plant.WYForecast{n_y-1});
    else 
        solution = dispatch_loop(gen,buildings,cool_tower,subnet,op_mat_a,op_matb,one_step,options,forecast_time,forecast,[]);
    end 
    solution.Timestamp = forecast_time;
    solution.hydroSOC = [hydro_soc_init;solution.hydroSOC];
    disp(strcat('Water Year Forecast Completed for ',num2str(year),':',num2str(year+1)))
    date = hold_date;
    Plant.WYForecast{n_y} = solution;
    TestData = hold_test_data;
    if n_y == 1
        for i = 1:1:length(gen)
            if strcmp(gen(i).Type,'Hydro Storage')
                n = gen(i).QPform.Hydro.subnetNode;%dam #
                gen(i).CurrentState(2) = interp1(solution.Timestamp,solution.hydroSOC(:,n),date(1));
            end
        end
    end
end
end%Ends function water_year_forecast