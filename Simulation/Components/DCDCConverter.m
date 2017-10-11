function Out = DCDCConverter(varargin)
%A very simple DC/DC Converter
%0 States
%6 Inlets: ILoad, VLoad, VMax, Ivmax, VMin, Ivmin
%3 Outlet: Imin, Imax, VSource
if length(varargin) ==1
    block = varargin{1};
    block.Efficiency = 0.90;%can be moved to model file
    block.IC = [];%no states
    
    block.InletPorts = {'Ivmax','VMax','Ivmin','VMin','ILoad','VLoad'};
    block.Ivmax.IC = -1;
    block.Ivmax.Saturation = [-inf,inf];
    block.VMax.IC = 50;
    block.Vmax.Saturation = [0,inf];
    block.Ivmin.IC = 1;
    block.Ivmin.Saturation = [0,inf];
    block.VMin.IC = 30;
    block.Vmin.Saturation = [0,inf];
    block.ILoad.IC = 0;
    block.ILoad.Saturation = [0,inf];
    block.VLoad.IC = 80;
    block.Vload.Saturation = [0,inf];
    
    block.OutletPorts = {'Imin','Imax','VSource'};
    block.Imin.IC = -1;
    block.Imax.IC = 1;
    block.VSource.IC = 40;
    Out = block;
elseif length(varargin) ==2
    block = varargin{1};
    Inlet = varargin{2};
    Inlet = checkSaturation(Inlet,block);
    Ivmax = block.Effic*Inlet.Ivmax*Inlet.VMax/Inlet.VLoad;%DC bus current when v = vmax
    Ivmin = block.Effic*Inlet.Ivmin*Inlet.VMin/Inlet.VLoad;%DC bus current when v = vmin
    
    %if current is negative, efficiency losses are applied at dc bus
    if Inlet.Ivmax < 0
        Ivmax = Inlet.Ivmax*Inlet.VMax/(Inlet.VLoad*block.Effic);
    end
    if Inlet.Ivmin < 0
        Ivmin = Inlet.Ivmin*Inlet.VMin/(Inlet.VLoad*block.Effic);
    end
    
    V = Inlet.VMin +(Inlet.VMax - Inlet.VMin)*(Inlet.ILoad - Ivmin)/(Ivmax - Ivmin);%voltage on source circuit that creates ILoad on DC bus
    
    IMin = min(Ivmax,Ivmin);
    IMax = max(Ivmax,Ivmin);
    
    block.Imin.IC = IMin;
    block.Imax.IC = IMax;
    block.VSource.IC = V;
    
    % need to update block.Scale
    Out = block;
else%running the model
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    Inlet = checkSaturation(Inlet,block);
    if strcmp(string1,'Outlet')
        Ivmax = block.Effic*Inlet.Ivmax*Inlet.VMax/Inlet.VLoad;%DC bus current when v = vmax
        Ivmin = block.Effic*Inlet.Ivmin*Inlet.VMin/Inlet.VLoad;%DC bus current when v = vmin

        %if current is negative, efficiency losses are applied at dc bus
        if Inlet.Ivmax < 0
            Ivmax = Inlet.Ivmax*Inlet.VMax/(Inlet.VLoad*block.Effic);
        end
        if Inlet.Ivmin < 0
            Ivmin = Inlet.Ivmin*Inlet.VMin/(Inlet.VLoad*block.Effic);
        end

        V = Inlet.VMin +(Inlet.VMax - Inlet.VMin)*(Inlet.ILoad - Ivmin)/(Ivmax - Ivmin);%voltage on source circuit that creates ILoad on DC bus

        IMin = min(Ivmax,Ivmin);
        IMax = max(Ivmax,Ivmin);

        Out.IMin = IMin;%min current applied to DC bus by source
        Out.IMax = IMax;%max current applied to DC bus by source

        Out.VSource = V;%voltage on source circuit
    elseif strcmp(string1,'dY')
        %no states
    end
end
end%Ends function DCDCConverter