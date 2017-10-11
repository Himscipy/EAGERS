function Out = Capacitor(varargin)
% Simple capcitor model
%1 inlets: V
%6 outlets: IcapMax, IcapMin, VcapMax, VcapMin, E, EMax
%2 states: Temp, Voltage
if length(varargin) ==1 %first initialization
    block = varargin{1};
    block.Tamb = 305;%ambient air temp
    block.AmbConvC = 5;%ambient convection coefficient
    block.Epsilon = .8;%radiation epsilon value
    block.Sigma = 5.67e-8;%Radiation sigma value
    block.SpecHeat = .5;%Specific Heat of Capacitor kJ/(kg*K)

    %%
    %can be pulled out of initialization function
    block.Mass = 1;% in kg
    block.Capacitance = 1;%in farrads
    block.VMax = 1;
    block.IMax = 1;
    block.Width = 1;
    block.Length = 1;
    block.Height = 1;
    block.InitialCharge = 1;%in volts
    block.Resistance = 1;%in ohms
    %%
    block.SurfA = 2*block.Width*block.Length + 2* block.Width*block.Height + 2*block.Length*block.Height;%Assume rectangular prism
    
    block.Scale = [block.Tamb, block.VMax];
    block.IC = [1,(block.InitialCharge/block.VMax)];
    block.UpperBound = [inf,1];
    block.LowerBound = [0,0];
    
    block.InletPorts = {'V'};
    block.V.IC = block.Scale(2)*block.IC(2);
    block.V.Saturation = [-inf,inf];
    
    block.OutletPorts = {'IcapMax','IcapMin','VcapMax','VcapMin','E','EMax'};
    block.IcapMax.IC = -block.IMax;%current at max applied voltage
    block.IcapMin.IC = block.IMax;%current at min applied voltage
    block.VcapMax.IC = block.Scale(2)*block.IC(2) + (block.IcapMax.IC*block.Resistance);
    block.VcapMin.IC = block.Scale(2)*block.IC(2) + (block.IcapMin.IC*block.Resistance);
    block.E.IC = (0.5*block.Capacitance*(block.Scale(2)*block.IC(2))^2)/1000;
    block.EMax.IC = (0.5*block.Capacitance*(block.VMax)^2)/1000;
    Out = block;
elseif length(varargin) ==2 %have inlets
    block = varargin{1};
    Inlet = varargin{2};
    
    Inlet = checkSaturation(Inlet,block);
    Y = block.IC.* block.Scale;
    %V = Y(2) - Inlet.V;%negative value means charging
    %I = V/block.Resistance;%negative value means charging
    
    QConv = (Y(1) - block.Tamb)*block.AmbConvC*block.SurfA/1000;
    QRad = ((Y(1))^4 - (block.Tamb)^4)*block.Epsilon*block.Sigma*block.SurfA/1000;
    
    delE = block.Mass*block.SpecHeat*(block.TMax - Y(1))/block.tscale;%Remaining Heat Sink energy in kJ divided over block.tscale seconds;
    IMaxT = (1000*(QConv + QRad + delE)/block.Resistance)^0.5;%Positive Current(amps) that would hit TMax after block.tscale seconds;
    
    IMax = min(block.IMax, IMaxT) - block.ILeak;%apply temperature and current limits
    
    VMax = min((Y(2)+IMax*block.Resistance), block.VMax);%apply voltage limits
    VMin = max(0,(Y(2)-IMax*block.Resistance));
    
    IMaxV = (Y(2) - VMax)/block.Resistance;
    IMinV = (Y(2) - VMin)/block.Resistance;
    
    E = (0.5*block.Capacitance*Y(2)^2)/1000;%Energy stored in capacitor in kJ
    EMax = (0.5*block.Capacitance*block.VMax^2)/1000;
    
    block.IcapMax.IC = IMaxV;
    block.IcapMin.IC = IMinV;
    block.VcapMax.IC = VMax;
    block.VcapMin.IC = VMin;
    block.E.IC = E;
    block.EMax.IC = EMax;
    % need to update block.Scale
    Out = block;
else%running the model
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    
    Inlet = checkSaturation(Inlet,block);
    V = Y(2) - Inlet.V;%negative value means charging
    I = V/block.Resistance;%negative value means charging


    QConv = (Y(1) - block.Tamb)*block.AmbConvC*block.SurfA/1000;
    QRad = ((Y(1))^4 - (block.Tamb)^4)*block.Epsilon*block.Sigma*block.SurfA/1000;
    if strcmp(string1,'Outlet')
        delE = block.Mass*block.SpecHeat*(block.TMax - Y(1))/block.tscale;%Remaining Heat Sink energy in kJ divided over block.tscale seconds;
        IMaxT = (1000*(QConv + QRad + delE)/block.Resistance)^0.5;%Positive Current(amps) that would hit TMax after block.tscale seconds;

        IMax = min(block.IMax, IMaxT) - block.ILeak;%apply current and temperature limitations

        VMax = min((Y(2)+IMax*block.Resistance), block.VMax);%apply voltage limitations
        VMin = max(0,(Y(2)-IMax*block.Resistance));

        IMaxV = (Y(2) - VMax)/block.Resistance;%current at max V
        IMinV = (Y(2) - VMin)/block.Resistance;%current at min V

        E = (0.5*block.Capacitance*Y(2)^2)/1000;%Energy stored in capacitor in kJ
        EMax = (0.5*block.Capacitance*block.VMax^2)/1000;

        Out.IcapMax = IMaxV;%Current from max circuit voltage(charge)
        Out.IcapMin = IMinV;%Current from min circuit voltage(discharge)
        Out.VcapMax = VMax;%max circuit voltage
        Out.VcapMin = VMin;%min circuit voltage
        Out.E = E;
        Out.EMax = EMax;
    elseif strcmp(string1,'dY')
        dY = 0*Y;

        %include effects of current leakage
        Itotal = 0;
        if I >= 0
            if Y(2) > 0
                Itotal = I + block.ILeak;
            end
        else%I < 0
            if Y(2) > 0
                Itotal= I - block.ILeak;
            end
        end

        QElec = Itotal*V/1000;%Waste heat in kW
        if I < 0  && Y(2) > 0%leakage current decreases effectiveness of charging current
            Itotal = Itotal +2*block.ILeak;
        end

        dY(1) = (QElec - QConv - QRad)/(block.Mass*block.SpecHeat);%temp change
        dY(2) = -Itotal/block.Capacitance;%voltage change
        Out = dY;
    end
end
end%Ends function Capacitor